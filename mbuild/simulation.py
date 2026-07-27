"""Simulation methods that operate on mBuild compounds and paths."""

import logging
import os

import gmso
import hoomd
import numpy as np
from gmso.parameterization import apply

from mbuild import Compound
from mbuild.utils.io import import_

logger = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────────────────────
# Default force field parameters (used when no forcefield is provided)
# ─────────────────────────────────────────────────────────────────────────────

DEFAULT_LJ_PARAMS = {"sigma": 0.34, "epsilon": 0.36}  # nm, kJ/mol
DEFAULT_BOND_PARAMS = {"k": 300000.0, "r0": 0.15}  # kJ/mol/nm^2, nm
DEFAULT_ANGLE_PARAMS = {"k": 400.0, "theta0": 1.9111}  # kJ/mol/rad^2, ~109.5 deg


# ─────────────────────────────────────────────────────────────────────────────
# Input resolution helper
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_input(system):
    """Convert a Path or Compound input into an mBuild Compound.

    Parameters
    ----------
    system : mb.Compound or mbuild.path.build.Path
        The input system.

    Returns
    -------
    compound : mb.Compound
    is_path : bool
        Whether the original input was a Path.
    original : the original input object
    """
    import mbuild.path.build

    if isinstance(system, mbuild.path.build.Path):
        return system.to_compound(), True, system
    elif isinstance(system, Compound):
        return system, False, system
    else:
        raise TypeError(
            f"Expected an mb.Compound or mbuild.path.Path, got {type(system)}."
        )


def _sync_back(compound, is_path, original):
    """If the original was a Path, copy positions back."""
    if is_path:
        original.coordinates = compound.xyz


# ─────────────────────────────────────────────────────────────────────────────
# Property Logger
# ─────────────────────────────────────────────────────────────────────────────


class PropertyLogger:
    """Collects energies and other scalar properties during a simulation.

    Parameters
    ----------
    log_every : int
        Frequency (in steps) at which properties are recorded.
    properties : list of str
        Which properties to log. Options depend on the backend:
        - HOOMD: "potential_energy", "kinetic_energy", "temperature", "pressure"
        - OpenMM: "potential_energy", "kinetic_energy", "temperature", "volume"
    """

    def __init__(self, log_every=100, properties=None):
        self.log_every = log_every
        if properties is None:
            properties = ["potential_energy"]
        self.properties = properties
        self.data = {p: [] for p in properties}

    def reset(self):
        self.data = {p: [] for p in self.properties}

    def as_dict(self):
        return {k: np.array(v) for k, v in self.data.items()}


# ─────────────────────────────────────────────────────────────────────────────
# HOOMD Simulation
# ─────────────────────────────────────────────────────────────────────────────


class HoomdSimulation(hoomd.simulation.Simulation):
    """A custom class for hoomd-based simulation methods.

    Accepts either an mb.Compound or an mbuild.path.Path as input.

    Parameters
    ----------
    system : mb.Compound or mbuild.path.build.Path
        The system to simulate. If a Path is given, it is converted
        to a Compound internally and positions are synced back after runs.
    forcefield : foyer/gmso Forcefield or None
        The forcefield to apply. If None, default LJ + harmonic bond/angle
        parameters are used.
    r_cut : float (nm)
        The cutoff distance used in the non-bonded pair neighborlist.
    fixed_compounds : list of mb.Compound, default None
    integrate_compounds : list of mb.Compound, default None
    run_on_gpu : bool or None, default None
    seed : int, default 1
    automatic_box : bool, default False
    box_buffer : float, default 10
    logger : PropertyLogger or None
    """

    def __init__(
        self,
        system,
        forcefield=None,
        r_cut=1.0,
        run_on_gpu=None,
        seed=1,
        automatic_box=False,
        box_buffer=10,
        integrate_compounds=None,
        fixed_compounds=None,
        gsd_file_name=None,
        logger=None,
    ):
        # Resolve input type
        compound, self._is_path, self._original_system = _resolve_input(system)

        if run_on_gpu is None:
            device = hoomd.device.auto_select()
        elif run_on_gpu:
            try:
                device = hoomd.device.GPU()
                print(f"GPU found, running on device {device.device}")
            except RuntimeError:
                print(
                    "Unable to find compatible GPU device. "
                    "Set `run_on_gpu = False` or see HOOMD documentation."
                )
                device = hoomd.device.CPU()
        else:
            device = hoomd.device.CPU()

        self.compound = compound
        self.forcefield = forcefield
        self.r_cut = r_cut
        self.integrate_compounds = integrate_compounds
        self.fixed_compounds = fixed_compounds
        self.automatic_box = automatic_box
        self.box_buffer = box_buffer
        self.energies = []
        self.logger = logger

        if compound._hoomd_data:
            last_snapshot, last_forces, last_forcefield = compound._get_sim_data()
            if forcefield == last_forcefield:
                snapshot = last_snapshot
                self.forces = last_forces
            else:
                snapshot, self.forces = self._to_hoomd_snap_forces()
                compound._add_sim_data(
                    state=snapshot, forces=self.forces, forcefield=forcefield
                )
        else:
            snapshot, self.forces = self._to_hoomd_snap_forces()
            compound._add_sim_data(
                state=snapshot, forces=self.forces, forcefield=forcefield
            )

        self.active_forces = []
        self.inactive_forces = []
        self._orig_force_params = {}
        super(HoomdSimulation, self).__init__(device=device, seed=seed)
        self.create_state_from_snapshot(snapshot)
        self._orig_force_params = self._snapshot_force_params()

    def _to_hoomd_snap_forces(self):
        """Convert compound to HOOMD snapshot and forces."""
        if not self.compound.box:
            self.compound.box = self.compound.get_boundingbox(pad_box=self.box_buffer)

        top = self.compound.to_gmso()
        top.identify_connections()

        if self.forcefield is not None:
            apply(
                top,
                forcefields=self.forcefield,
                ignore_params=["dihedral", "improper"],
            )
            forces, _ = gmso.external.to_hoomd_forcefield(top, r_cut=self.r_cut)
            snap, _ = gmso.external.to_gsd_snapshot(top=top)
            forces = list(set().union(*forces.values()))
        else:
            snap, _ = gmso.external.to_gsd_snapshot(top=top)
            forces = self._build_default_forces(top, snap)

        return snap, forces

    def _build_default_forces(self, top, snap):
        """Build default HOOMD forces (LJ, bonds, angles) from topology."""
        particle_types = list(set(snap.particles.types))

        nlist = hoomd.md.nlist.Cell(buffer=0.4)
        lj = hoomd.md.pair.LJ(nlist=nlist, default_r_cut=self.r_cut)
        for i, t1 in enumerate(particle_types):
            for t2 in particle_types[i:]:
                lj.params[(t1, t2)] = dict(
                    sigma=DEFAULT_LJ_PARAMS["sigma"],
                    epsilon=DEFAULT_LJ_PARAMS["epsilon"],
                )
        forces = [lj]

        if top.n_bonds > 0:
            bond_types = list(
                set(
                    tuple(
                        sorted(
                            [b.connection_members[0].name, b.connection_members[1].name]
                        )
                    )
                    for b in top.bonds
                )
            )
            bond_force = hoomd.md.bond.Harmonic()
            for bt in bond_types:
                type_str = f"{bt[0]}-{bt[1]}"
                bond_force.params[type_str] = dict(
                    k=DEFAULT_BOND_PARAMS["k"], r0=DEFAULT_BOND_PARAMS["r0"]
                )
            forces.append(bond_force)

        if top.n_angles > 0:
            angle_types = list(
                set(
                    tuple(
                        [
                            a.connection_members[0].name,
                            a.connection_members[1].name,
                            a.connection_members[2].name,
                        ]
                    )
                    for a in top.angles
                )
            )
            angle_force = hoomd.md.angle.Harmonic()
            for at in angle_types:
                type_str = f"{at[0]}-{at[1]}-{at[2]}"
                angle_force.params[type_str] = dict(
                    k=DEFAULT_ANGLE_PARAMS["k"], t0=DEFAULT_ANGLE_PARAMS["theta0"]
                )
            forces.append(angle_force)

        return forces

    def _snapshot_force_params(self):
        """Snapshot current force params for later restoration."""
        orig = {}
        for force in self.forces:
            if not hasattr(force, "params"):
                continue
            orig[id(force)] = {
                param: dict(force.params[param]) for param in force.params
            }
        return orig

    def get_integrate_group(self):
        """Return the HOOMD filter for the integration group."""
        if self.integrate_compounds and not self.fixed_compounds:
            if all([isinstance(i, Compound) for i in self.integrate_compounds]):
                integrate_indices = []
                for comp in self.integrate_compounds:
                    integrate_indices.extend(
                        list(self.compound.get_child_indices(comp))
                    )
            elif all([isinstance(i, int) for i in self.integrate_compounds]):
                integrate_indices = self.integrate_compounds
            return hoomd.filter.Tags(integrate_indices)
        elif self.fixed_compounds and not self.integrate_compounds:
            if all([isinstance(i, Compound) for i in self.fixed_compounds]):
                fix_indices = []
                for comp in self.fixed_compounds:
                    fix_indices.extend(list(self.compound.get_child_indices(comp)))
            elif all([isinstance(i, int) for i in self.fixed_compounds]):
                fix_indices = self.fixed_compounds
            return hoomd.filter.SetDifference(
                hoomd.filter.All(), hoomd.filter.Tags(fix_indices)
            )
        elif self.integrate_compounds and self.fixed_compounds:
            raise RuntimeError(
                "You can specify only one of integrate_compounds and fixed_compounds."
            )
        else:
            return hoomd.filter.All()

    def get_force(self, instance, active_only=False):
        if active_only:
            iterForcesSet = set(self.active_forces)
        else:
            iterForcesSet = set(self.forces)
        for force in iterForcesSet:
            if isinstance(force, instance):
                return force
        raise ValueError(f"No force of {instance} was found.")

    def get_dpd_from_lj(self, A):
        """Make a best-guess DPD force from types and parameters of an LJ force."""
        lj_force = self.get_force(hoomd.md.pair.LJ)
        dpd = hoomd.md.pair.DPDConservative(nlist=lj_force.nlist)
        for param in lj_force.params:
            dpd.params[param] = dict(A=A)
            dpd.r_cut[param] = lj_force.params[param]["sigma"]
        return dpd

    def set_fire_integrator(
        self,
        dt,
        force_tol,
        angmom_tol,
        energy_tol,
        methods,
        finc_dt,
        fdec_dt,
        alpha_start,
        fdec_alpha,
        min_steps_adapt,
        min_steps_conv,
    ):
        fire = hoomd.md.minimize.FIRE(
            dt=dt,
            force_tol=force_tol,
            angmom_tol=angmom_tol,
            energy_tol=energy_tol,
            finc_dt=finc_dt,
            fdec_dt=fdec_dt,
            alpha_start=alpha_start,
            fdec_alpha=fdec_alpha,
            min_steps_adapt=min_steps_adapt,
            min_steps_conv=min_steps_conv,
            methods=methods,
        )
        fire.forces = self.active_forces
        self.operations.integrator = fire

    def set_integrator(self, dt, method):
        integrator = hoomd.md.Integrator(dt=dt)
        integrator.forces = self.active_forces
        integrator.methods = [method]
        self.operations.integrator = integrator

    def _update_integrator_forces(self):
        self.operations.integrator.forces = self.active_forces

    def _update_snapshot(self):
        snapshot = self.state.get_snapshot()
        self.compound._add_sim_data(state=snapshot)

    def update_positions(self):
        """Update compound positions from snapshot and sync back to Path if needed."""
        with self.state.cpu_local_snapshot as snap:
            particles = snap.particles.rtag[:]
            pos = snap.particles.position[particles]
            self.compound.xyz = np.copy(pos)
        _sync_back(self.compound, self._is_path, self._original_system)

    def store_current_energies(self):
        """Store energy from active forces."""
        new_energy = {}
        for force in self.active_forces:
            ffname = force.__class__.__module__ + "." + force.__class__.__name__
            new_energy[ffname] = force.energy
        if new_energy:
            self.energies.append(new_energy)

    def get_energy(self):
        """Return energies by force type across all stored frames."""
        if not self.energies:
            print(f"No energies currently stored in {self}")
            return None
        returnDict = {}
        n_frames = len(self.energies)
        for frame, energyDict in enumerate(self.energies):
            for ffkey in energyDict:
                if ffkey not in returnDict:
                    returnDict[ffkey] = np.zeros(n_frames)
                returnDict[ffkey][frame] = energyDict[ffkey]
        return returnDict

    def _log_properties(self):
        """Record properties to the logger if one is attached."""
        if self.logger is None:
            return
        for prop in self.logger.properties:
            if prop == "potential_energy":
                pe = sum(f.energy for f in self.active_forces)
                self.logger.data["potential_energy"].append(pe)
            elif prop == "kinetic_energy":
                self.logger.data["kinetic_energy"].append(None)
            elif prop == "temperature":
                self.logger.data["temperature"].append(None)

    # ─────────────────────────────────────────────────────────────────────
    # Simulation methods
    # ─────────────────────────────────────────────────────────────────────

    def nvt(
        self,
        forces_handler,
        n_steps,
        kT,
        dt,
        tau,
        thermostat=hoomd.md.methods.thermostats.Berendsen,
    ):
        """Run an NVT simulation.

        Parameters
        ----------
        forces_handler : ForcesHandler or None
        n_steps : int
        kT : float
        dt : float
        tau : float
        thermostat : hoomd thermostat class
        """
        if forces_handler is not None:
            forces_handler.scale_sim(self)
        else:
            self.active_forces = list(self.forces)

        nvt_method = hoomd.md.methods.ConstantVolume(
            filter=self.get_integrate_group(),
            thermostat=thermostat(kT=kT, tau=tau),
        )
        self.set_integrator(method=nvt_method, dt=dt)
        self.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

        if self.energies == []:
            self.run(0)
            self.store_current_energies()

        self.run(n_steps)
        self.store_current_energies()
        self._log_properties()
        self.operations.integrator = None
        self.update_positions()
        self._update_snapshot()

    def npt(
        self,
        forces_handler,
        n_steps,
        kT,
        dt,
        tau,
        pressure,
        tauS,
        couple="xyz",
        thermostat=hoomd.md.methods.thermostats.Berendsen,
    ):
        """Run an NPT simulation."""
        if forces_handler is not None:
            forces_handler.scale_sim(self)
        else:
            self.active_forces = list(self.forces)

        npt_method = hoomd.md.methods.ConstantPressure(
            filter=self.get_integrate_group(),
            thermostat=thermostat(kT=kT, tau=tau),
            S=pressure,
            tauS=tauS,
            couple=couple,
        )
        self.set_integrator(method=npt_method, dt=dt)
        self.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)

        if self.energies == []:
            self.run(0)
            self.store_current_energies()

        self.run(n_steps)
        self.store_current_energies()
        self._log_properties()
        self.operations.integrator = None
        self.update_positions()
        self._update_snapshot()

    def cap_displacement(
        self,
        forces_handler,
        n_steps,
        dt,
        max_displacement=1e-3,
    ):
        """Run a capped-displacement simulation.

        Parameters
        ----------
        forces_handler : ForcesHandler or None
        n_steps : int
        dt : float
        max_displacement : float (nm)
        """
        self.compound._kick()
        if forces_handler is not None:
            forces_handler.scale_sim(self)
        else:
            self.active_forces = list(self.forces)

        displacement_capped = hoomd.md.methods.DisplacementCapped(
            filter=self.get_integrate_group(),
            maximum_displacement=max_displacement,
        )
        self.set_integrator(method=displacement_capped, dt=dt)

        if self.energies == []:
            self.run(0)
            self.store_current_energies()

        self.run(n_steps)
        self.store_current_energies()
        self._log_properties()
        self.operations.integrator = None
        self.update_positions()
        self._update_snapshot()

    def fire(
        self,
        forces_handler,
        n_steps,
        n_iterations=1,
        dt=1e-5,
        min_steps_adapt=5,
        min_steps_conv=100,
        finc_dt=1.1,
        fdec_dt=0.5,
        alpha_start=0.1,
        fdec_alpha=0.95,
        force_tol=1e-2,
        angmom_tol=1e-2,
        energy_tol=1e-6,
    ):
        """Run FIRE energy minimization.

        Parameters
        ----------
        forces_handler : ForcesHandler or None
        n_steps : int
        n_iterations : int
        dt : float
        """
        self.compound._kick()
        if forces_handler is not None:
            forces_handler.scale_sim(self)
        else:
            self.active_forces = list(self.forces)

        nvt_method = hoomd.md.methods.ConstantVolume(filter=self.get_integrate_group())
        self.set_fire_integrator(
            dt=dt,
            force_tol=force_tol,
            angmom_tol=angmom_tol,
            energy_tol=energy_tol,
            finc_dt=finc_dt,
            fdec_dt=fdec_dt,
            alpha_start=alpha_start,
            fdec_alpha=fdec_alpha,
            min_steps_adapt=min_steps_adapt,
            min_steps_conv=min_steps_conv,
            methods=[nvt_method],
        )

        if self.energies == []:
            self.run(0)
            self.store_current_energies()

        for _ in range(n_iterations):
            self.run(n_steps)
            self.operations.integrator.reset()

        self.store_current_energies()
        self._log_properties()
        self.operations.integrator = None
        self._update_snapshot()
        self.update_positions()


# ─────────────────────────────────────────────────────────────────────────────
# OpenMM Simulation
# ─────────────────────────────────────────────────────────────────────────────


class OpenMMSimulation:
    """A wrapper around OpenMM for running simulations on mBuild Compounds or Paths.

    Accepts either an mb.Compound or an mbuild.path.Path as input.

    Parameters
    ----------
    system : mb.Compound or mbuild.path.build.Path
        The system to simulate.
    forcefield : str or None
        A foyer-compatible forcefield name or XML path.
        If None, default LJ + harmonic bond/angle parameters are applied.
    box : mb.Box or None
        Simulation box. If None, uses compound.box or a bounding box.
    constraints : str or None, default is None
        Bond constraints: None, "HBonds", "AllBonds", "HAngles".
    platform : str
        OpenMM platform: "CPU", "CUDA", "OpenCL", "Reference".
    nthreads : int
        Number of CPU threads (only relevant for CPU platform).
    logger : PropertyLogger or None
        Optional property logger.
    """

    def __init__(
        self,
        system,
        forcefield=None,
        box=None,
        constraints=None,
        platform="CPU",
        nthreads=1,
        logger=None,
    ):
        from openmm import Platform
        from openmm.app import AllBonds, HAngles, HBonds

        # Resolve input type
        compound, self._is_path, self._original_system = _resolve_input(system)

        self.compound = compound
        self.forcefield_name = forcefield
        self.logger = logger
        self.energies = []

        # Set up platform
        self._platform = Platform.getPlatformByName(platform)
        if platform == "CPU":
            self._platform.setPropertyDefaultValue("Threads", str(nthreads))

        # Parse constraints
        constraint_map = {
            None: None,
            "AllBonds": AllBonds,
            "HBonds": HBonds,
            "HAngles": HAngles,
        }
        if constraints not in constraint_map:
            raise ValueError(
                f"Invalid constraints: {constraints}. "
                f"Options: {list(constraint_map.keys())}"
            )
        self._constraints = constraint_map[constraints]

        # Build system
        if forcefield is not None:
            self._build_from_forcefield(compound, forcefield)
        else:
            self._build_default_system(compound)

    def _build_from_forcefield(self, compound, forcefield):
        """Build OpenMM system using foyer."""
        import openmm

        foyer = import_("foyer")

        to_parmed = compound.to_parmed()
        extension = os.path.splitext(forcefield)[-1]
        if extension == ".xml":
            ff = foyer.Forcefield(forcefield_files=forcefield)
        else:
            ff = foyer.Forcefield(name=forcefield)
        self._parmed = ff.apply(to_parmed)

        if self._constraints is not None:
            self.system = self._parmed.createSystem(constraints=self._constraints)
        else:
            self.system = self._parmed.createSystem()

        # For non-periodic compounds, ensure NonbondedForce uses NoCutoff
        if not compound.box:
            for i in range(self.system.getNumForces()):
                force = self.system.getForce(i)
                if isinstance(force, openmm.NonbondedForce):
                    force.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        self.topology = self._parmed.topology
        self.positions = self._parmed.positions

    def _build_default_system(self, compound):
        """Build OpenMM system with default LJ + bond + angle parameters."""
        import openmm
        import openmm.unit as u
        from openmm.app import Topology

        positions = compound.xyz  # nm
        n_particles = compound.n_particles

        self.system = openmm.System()
        for i in range(n_particles):
            self.system.addParticle(12.0 * u.amu)

        # Set periodic box vectors if compound has a box
        if compound.box:
            Lx, Ly, Lz = compound.box.lengths  # nm
            self.system.setDefaultPeriodicBoxVectors(
                openmm.Vec3(Lx, 0, 0) * u.nanometer,
                openmm.Vec3(0, Ly, 0) * u.nanometer,
                openmm.Vec3(0, 0, Lz) * u.nanometer,
            )

        # Default LJ nonbonded
        nonbonded = openmm.NonbondedForce()
        if compound.box:
            nonbonded.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
            nonbonded.setCutoffDistance(1.0 * u.nanometer)
        else:
            nonbonded.setNonbondedMethod(openmm.NonbondedForce.NoCutoff)

        sigma = DEFAULT_LJ_PARAMS["sigma"] * u.nanometer
        epsilon = DEFAULT_LJ_PARAMS["epsilon"] * u.kilojoules_per_mole
        for i in range(n_particles):
            nonbonded.addParticle(0.0, sigma, epsilon)
        self.system.addForce(nonbonded)

        # Get topology info via GMSO
        top = compound.to_gmso()
        top.identify_connections()
        sites_list = list(top.sites)

        # Default harmonic bonds
        if top.n_bonds > 0:
            bond_force = openmm.HarmonicBondForce()
            for bond in top.bonds:
                members = bond.connection_members
                i_idx = sites_list.index(members[0])
                j_idx = sites_list.index(members[1])
                bond_force.addBond(
                    i_idx,
                    j_idx,
                    DEFAULT_BOND_PARAMS["r0"] * u.nanometer,
                    DEFAULT_BOND_PARAMS["k"] * u.kilojoules_per_mole / u.nanometer**2,
                )
                nonbonded.addException(i_idx, j_idx, 0.0, 1.0, 0.0)
            self.system.addForce(bond_force)

        # Default harmonic angles
        if top.n_angles > 0:
            angle_force = openmm.HarmonicAngleForce()
            for angle in top.angles:
                members = angle.connection_members
                i_idx = sites_list.index(members[0])
                j_idx = sites_list.index(members[1])
                k_idx = sites_list.index(members[2])
                angle_force.addAngle(
                    i_idx,
                    j_idx,
                    k_idx,
                    DEFAULT_ANGLE_PARAMS["theta0"] * u.radian,
                    DEFAULT_ANGLE_PARAMS["k"] * u.kilojoules_per_mole / u.radian**2,
                )
            self.system.addForce(angle_force)

        # Build OpenMM topology
        particles_list = list(compound.particles())
        topology = Topology()
        chain = topology.addChain()
        for i, p in enumerate(particles_list):
            residue = topology.addResidue(f"R{i}", chain)
            topology.addAtom(p.name, openmm.app.Element.getByMass(12), residue)
        atoms_list = list(topology.atoms())
        if top.n_bonds > 0:
            for bond in top.bonds:
                members = bond.connection_members
                i_idx = sites_list.index(members[0])
                j_idx = sites_list.index(members[1])
                topology.addBond(atoms_list[i_idx], atoms_list[j_idx])

        self.topology = topology
        self.positions = positions * u.nanometer

    def _create_simulation(self, integrator):
        """Create an OpenMM Simulation object."""
        from openmm.app.simulation import Simulation

        self.simulation = Simulation(
            self.topology, self.system, integrator, self._platform
        )
        self.simulation.context.setPositions(self.positions)

    def _update_compound_positions(self):
        """Write positions back to compound and sync to Path if needed."""
        state = self.simulation.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True)
        pos_array = np.array(pos)  # nm
        if not np.any(np.isnan(pos_array)):
            self.compound.xyz = pos_array
            # Store updated positions so next _create_simulation uses them
            self.positions = pos
        else:
            logger.warning("NaN positions detected; compound not updated.")
            return
        _sync_back(self.compound, self._is_path, self._original_system)

    def _record_energies(self):
        """Record current potential energy."""
        import openmm.unit as u

        state = self.simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
        self.energies.append({"potential_energy": pe})
        if self.logger and "potential_energy" in self.logger.properties:
            self.logger.data["potential_energy"].append(pe)
        if self.logger and "kinetic_energy" in self.logger.properties:
            ke_state = self.simulation.context.getState(getEnergy=True)
            ke = ke_state.getKineticEnergy().value_in_unit(u.kilojoules_per_mole)
            self.logger.data["kinetic_energy"].append(ke)
        if self.logger and "temperature" in self.logger.properties:
            # Approximate from KE
            self.logger.data["temperature"].append(None)
        if self.logger and "volume" in self.logger.properties:
            box_vecs = state.getPeriodicBoxVectors(asNumpy=True)
            vol = np.dot(box_vecs[0], np.cross(box_vecs[1], box_vecs[2]))
            self.logger.data["volume"].append(float(vol))

    def get_energy(self):
        """Return recorded energies as dict of arrays."""
        if not self.energies:
            return None
        keys = self.energies[0].keys()
        return {k: np.array([e.get(k) for e in self.energies]) for k in keys}

    # ─────────────────────────────────────────────────────────────────────
    # Simulation methods
    # ─────────────────────────────────────────────────────────────────────

    def minimize(self, steps=1000, tolerance=10.0):
        """Energy minimize the system.

        Parameters
        ----------
        steps : int
        tolerance : float
            Tolerance in kJ/mol/nm.
        """
        import openmm.unit as u
        from openmm.openmm import LangevinIntegrator

        integrator = LangevinIntegrator(
            298 * u.kelvin, 1.0 / u.picosecond, 0.002 * u.picoseconds
        )
        self._create_simulation(integrator)
        self.simulation.minimizeEnergy(
            maxIterations=steps,
            tolerance=tolerance * u.kilojoules_per_mole / u.nanometer,
        )
        self._record_energies()
        self._update_compound_positions()

    def nvt(self, n_steps, kT=300, dt=0.002, friction=1.0, report_interval=None):
        """Run an NVT (Langevin) simulation.

        Parameters
        ----------
        n_steps : int
        kT : float
            Temperature in Kelvin.
        dt : float
            Timestep in picoseconds.
        friction : float
            Friction coefficient in 1/ps.
        report_interval : int or None
            If given, record energies every this many steps.
        """
        import openmm.unit as u
        from openmm.openmm import LangevinIntegrator

        integrator = LangevinIntegrator(
            kT * u.kelvin, friction / u.picosecond, dt * u.picoseconds
        )
        self._create_simulation(integrator)

        if report_interval:
            for i in range(0, n_steps, report_interval):
                chunk = min(report_interval, n_steps - i)
                self.simulation.step(chunk)
                self._record_energies()
        else:
            self.simulation.step(n_steps)
            self._record_energies()

        self._update_compound_positions()

    def npt(
        self,
        n_steps,
        kT=300,
        dt=0.002,
        friction=1.0,
        pressure=1.0,
        barostat_interval=25,
        report_interval=None,
    ):
        """Run an NPT simulation using Monte Carlo barostat.

        Parameters
        ----------
        n_steps : int
        kT : float (Kelvin)
        dt : float (ps)
        friction : float (1/ps)
        pressure : float (atm)
        barostat_interval : int
        report_interval : int or None
        """
        import openmm
        import openmm.unit as u
        from openmm.openmm import LangevinIntegrator

        barostat = openmm.MonteCarloBarostat(
            pressure * u.atmospheres, kT * u.kelvin, barostat_interval
        )
        self.system.addForce(barostat)

        integrator = LangevinIntegrator(
            kT * u.kelvin, friction / u.picosecond, dt * u.picoseconds
        )
        self._create_simulation(integrator)

        if report_interval:
            for i in range(0, n_steps, report_interval):
                chunk = min(report_interval, n_steps - i)
                self.simulation.step(chunk)
                self._record_energies()
        else:
            self.simulation.step(n_steps)
            self._record_energies()

        self._update_compound_positions()


# ─────────────────────────────────────────────────────────────────────────────
# ForcesHandler
# ─────────────────────────────────────────────────────────────────────────────


class ForcesHandler:
    """Provides a quick method for tuning HOOMD forces.

    Parameters
    ----------
    dpd : float, default 0.0
        If > 0, replaces LJ with DPDConservative.
    scale_bond : float, default 1.0
    scale_angle : float, default 1.0
    scale_lj : float, default 1.0
    scale_periodic : float, default 0.0
    scale_opls : float, default 0.0
    scale_improper : float, default 0.0
    scale_charge : float, default 0.0
    """

    def __init__(
        self,
        dpd=0.0,
        scale_angle=1,
        scale_bond=1,
        scale_lj=1,
        scale_periodic=0,
        scale_opls=0,
        scale_improper=0,
        scale_charge=0,
    ):
        self.dpd = dpd
        self.scale_forces = {
            "lj": scale_lj,
            "charge": scale_charge,
            "bond": scale_bond,
            "angle": scale_angle,
            "opls": scale_opls,
            "periodic": scale_periodic,
            "improper": scale_improper,
        }
        self.forcesDict = {}

    def scale_sim(self, sim):
        """Iterate through HOOMD force objects and apply scaling factors."""
        sim.active_forces = []
        forcesDict = {
            "lj": (hoomd.md.pair.LJ, ("epsilon")),
            "charge": (hoomd.md.special_pair.Coulomb, ("alpha")),
            "bond": (hoomd.md.bond.Harmonic, ("k")),
            "angle": (hoomd.md.angle.Harmonic, ("k")),
            "opls": (hoomd.md.dihedral.OPLS, "k1", "k2", "k3", "k4"),
            "periodic": (hoomd.md.dihedral.Periodic, ("k")),
            "improper": (hoomd.md.improper.Periodic, ("k")),
        }
        for key, scalar in self.scale_forces.items():
            if not scalar:
                continue
            if key == "lj" and self.dpd:
                dpd = sim.get_dpd_from_lj(A=self.dpd)
                sim.active_forces.append(dpd)
                self.forcesDict["lj"] = sim.get_force(forcesDict["lj"][0])
                self.forcesDict["dpd"] = dpd
                continue

            try:
                force = sim.get_force(forcesDict[key][0])
            except ValueError:
                continue
            orig_params = sim._orig_force_params.get(id(force), {})
            for param in force.params:
                for term in forcesDict[key][1:]:
                    if param in orig_params and term in orig_params[param]:
                        force.params[param][term] = orig_params[param][term] * scalar
                    else:
                        force.params[param][term] *= scalar
            sim.active_forces.append(force)
            self.forcesDict[key] = force
