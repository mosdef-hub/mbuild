"""Simulation methods that operate on mBuild compounds."""

import os
import tempfile
from warnings import warn

import gmso
import hoomd
import numpy as np
from ele.element import element_from_name, element_from_symbol
from ele.exceptions import ElementError
from gmso.parameterization import apply

from mbuild import Compound
from mbuild.exceptions import MBuildError
from mbuild.utils.io import import_


class HoomdSimulation(hoomd.simulation.Simulation):
    """A custom class to help in creating internal hoomd-based simulation methods.

    See ``hoomd_cap_displacement`` and ``hoomd_fire``.

    Parameters
    ----------
    compound : mb.Compound
        The compound to use in the simulation
    forcefield : foyer.forcefield.Forcefield or gmso.core.Forcefield
        The forcefield to apply to the system
    r_cut : float (nm)
        The cutoff distance (nm) used in the non-bonded pair neighborlist.
        Use smaller values for unstable starting conditions and faster performance.
    fixed_compounds : list of mb.Compound, default None
        If given, these compounds will be removed from the integration updates and
        held "frozen" during the hoomd simulation.
        They are still able to interact with other particles in the simulation.
        If desired, pass in a subset of children from `compound.children`.
    integrate_compounds : list of mb.Compound, default None
        If given, then only these compounds will be updated during integration
        and all other compunds in `compound` are removed from the integration group.
    run_on_gpu : bool, default False
        When `True` the HOOMD simulation uses the hoomd.device.GPU() device.
        This requires that you have a GPU compatible HOOMD install and
        a compatible GPU device.
    seed : int, default 42
        The seed passed to HOOMD
    """

    def __init__(
        self,
        compound,
        forcefield,
        r_cut,
        run_on_gpu,
        seed,
        automatic_box=False,
        box_buffer=1.5,
        integrate_compounds=None,
        fixed_compounds=None,
        gsd_file_name=None,
    ):
        if run_on_gpu:
            try:
                device = hoomd.device.GPU()
                print(f"GPU found, running on device {device.device}")
            except RuntimeError:
                print(
                    "Unable to find compatible GPU device. ",
                    "set `run_on_gpu = False` or see HOOMD documentation "
                    "for further information about GPU support.",
                )
        else:
            device = hoomd.device.CPU()
        self.compound = compound
        self.forcefield = forcefield
        self.r_cut = r_cut
        self.integrate_compounds = integrate_compounds
        self.fixed_compounds = fixed_compounds
        self.automatic_box = automatic_box
        self.box_buffer = box_buffer
        # TODO
        # self.set_box
        # Check if a hoomd sim method has been used on this compound already
        if compound._hoomd_data:
            last_snapshot, last_forces, last_forcefield = compound._get_sim_data()
            # Check if the forcefield passed this time matches last time
            if forcefield == last_forcefield:
                snapshot = last_snapshot
                self.forces = last_forces
                forcefield = last_forcefield
            else:  # New foyer/gmso forcefield has been passed, reapply
                snapshot, self.forces = self._to_hoomd_snap_forces()
                compound._add_sim_data(
                    state=snapshot, forces=self.forces, forcefield=forcefield
                )
        else:  # Sim method not used on this compound previously
            snapshot, self.forces = self._to_hoomd_snap_forces()
            compound._add_sim_data(
                state=snapshot, forces=self.forces, forcefield=forcefield
            )
        # Place holders for forces added/changed by specific methods below
        self.active_forces = []
        self.inactive_forces = []
        super(HoomdSimulation, self).__init__(device=device, seed=seed)
        self.create_state_from_snapshot(snapshot)

    def _to_hoomd_snap_forces(self):
        # Convret to GMSO, apply forcefield
        top = self.compound.to_gmso()
        top.identify_connections()
        apply(top, forcefields=self.forcefield, ignore_params=["dihedral", "improper"])
        # Get hoomd snapshot and force objects
        forces, ref = gmso.external.to_hoomd_forcefield(top, r_cut=self.r_cut)
        snap, ref = gmso.external.to_gsd_snapshot(top)
        forces = list(set().union(*forces.values()))
        return snap, forces

    def get_integrate_group(self):
        # Get indices of compounds to include in integration
        if self.integrate_compounds and not self.fixed_compounds:
            integrate_indices = []
            for comp in self.integrate_compounds:
                integrate_indices.extend(list(self.compound.get_child_indices(comp)))
            return hoomd.filter.Tags(integrate_indices)
        # Get indices of compounds to NOT include in integration
        elif self.fixed_compounds and not self.integrate_compounds:
            fix_indices = []
            for comp in self.fixed_compounds:
                fix_indices.extend(list(self.compound.get_child_indices(comp)))
            return hoomd.filter.SetDifference(
                hoomd.filter.All(), hoomd.filter.Tags(fix_indices)
            )
        # Both are passed in, not supported
        elif self.integrate_compounds and self.fixed_compounds:
            raise RuntimeError(
                "You can specify only one of integrate_compounds and fixed_compounds."
            )
        # Neither are given, include everything in integration
        else:
            return hoomd.filter.All()

    def get_force(self, instance):
        for force in set(self.forces + self.active_forces + self.inactive_forces):
            if isinstance(force, instance):
                return force

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
        # snapshot_copy = hoomd.Snapshot(snapshot)  # full copy
        self.compound._add_sim_data(state=snapshot)

    def add_gsd_writer(self, file_name, write_period):
        """"""
        pass

    def update_positions(self):
        """Update compound positions from snapshot."""
        with self.state.cpu_local_snapshot as snap:
            particles = snap.particles.rtag[:]
            pos = snap.particles.position[particles]
            self.compound.xyz = np.copy(pos)


## HOOMD METHODS ##
def hoomd_cap_displacement(
    compound,
    forcefield,
    n_steps,
    dt,
    r_cut,
    max_displacement=1e-3,
    fixed_compounds=None,
    integrate_compounds=None,
    dpd_A=None,
    bond_k_scale=1,
    angle_k_scale=1,
    n_relax_steps=0,
    run_on_gpu=False,
    seed=42,
):
    """Run a short simulation with hoomd.md.methods.DisplacementCapped

    Parameters
    ----------
    compound : mb.Compound
        The compound to use in the simulation
    forcefield : foyer.forcefield.Forcefield or gmso.core.Forcefield
        The forcefield to apply to the system
    n_steps : int
        The number of simulation time steps to run with the modified forcefield
        created by bond_k_scale, angle_k_scale, and dpd_A.
        See these parameters for more information.
    dt : float
        The simulation timestep. Note that for running capped displacement on highly unstable
        systems (e.g. overlapping particles), it is often more stable to use a larger
        value of dt and let the displacement cap limit position updates.
    r_cut : float (nm)
        The cutoff distance (nm) used in the non-bonded pair neighborlist.
        Use smaller values for unstable starting conditions and faster performance.
    max_displacement : float (nm), default = 1e-3 nm
        The maximum displacement (nm) allowed per timestep. Use smaller values for
        highly unstable systems.
    fixed_compounds : list of mb.Compound, default None
        If given, these compounds will be removed from the integration updates and
        held "frozen" during the hoomd simulation.
        They are still able to interact with other particles in the simulation.
        If desired, pass in a subset of children from `compound.children`.
    integrate_compounds : list of mb.Compound, default None
        If given, then only these compounds will be updated during integration
        and all other compunds in `compound` are removed from the integration group.
    dpd_A : float, default None
        If set to a value, then the initial simulation replaces the LJ 12-6 pair potential
        with the softer hoomd.md.pair.DPDConservative pair force from HOOMD.
        This is used for `n_steps` and replaced with the original LJ for potential which is
        used for n_relax_steps. This can be useful for highly unstable starting configutations.
        Note that the cutoff for the DPD force is set to the sigma value for each atom type
        and force cosntant (A) is set to dpd_A for all atom types.
    bond_k_scale : float, default 1
        Scale the bond-stretching force constants so that K_new = K * bond_k_scale.
        This can be useful for relaxing highly unstable starting configutations.
        The bond force constants are scaled back to their original values
        during the run of `n_relax_steps`
    angle_k_scale : float, default 1
        Scale the bond-angle force constants so that K_new = K * angle_k_scale.
        This can be useful for relaxing highly unstable starting configutations.
        The angle force constants are scaled back to their original values
        during the run of `n_relax_steps`
    n_relax_steps: int, optional
        The number of steps to run after running for `n_steps` with the modified
        forcefield. This is designed to be used when utilizing the parameters
        `dpd_A`, `bond_k_scale`, and `angle_k_scale`.
    run_on_gpu : bool, default False
        When `True` the HOOMD simulation uses the hoomd.device.GPU() device.
        This requires that you have a GPU compatible HOOMD install and
        a compatible GPU device.
    seed : int, default 42
        The seed passed to HOOMD

    Notes
    -----
    Running on GPU:
        The hoomd-based methods in ``mbuild.simulation``
        are able to utilize and run on GPUs to perform energy minimization for larger systems if needed.
        Performance improvements of using GPU over CPU may not be significant until systems reach a
        size of ~1,000 particles.

    Running multiple simulations:
        The information needed to run the HOOMD simulation is saved after the first simulation.
        Therefore, calling hoomd-based simulation methods multiple times (e.g, in a for loop or while loop)
        does not require re-running performing atom-typing, applying the forcefield, or running hoomd
        format writers again.
    """
    compound._kick()
    sim = HoomdSimulation(
        compound=compound,
        forcefield=forcefield,
        r_cut=r_cut,
        run_on_gpu=run_on_gpu,
        seed=seed,
        integrate_compounds=integrate_compounds,
        fixed_compounds=fixed_compounds,
    )
    bond = sim.get_force(hoomd.md.bond.Harmonic)
    angle = sim.get_force(hoomd.md.angle.Harmonic)
    lj = sim.get_force(hoomd.md.pair.LJ)
    # If needed scale bond K and angle K
    if bond_k_scale != 1.0:
        for param in bond.params:
            bond.params[param]["k"] *= bond_k_scale
    sim.active_forces.append(bond)
    # If angle_k_scale is zero, then "turn off" angles; skip adding the angle force
    if angle_k_scale != 0:
        if angle_k_scale != 1.0:
            for param in angle.params:
                angle.params[param]["k"] *= angle_k_scale
        sim.active_forces.append(angle)
    if dpd_A:
        dpd = sim.get_dpd_from_lj(A=dpd_A)
        sim.active_forces.append(dpd)
    else:
        sim.active_forces.append(lj)
    # Set up and run
    displacement_capped = hoomd.md.methods.DisplacementCapped(
        filter=sim.get_integrate_group(),
        maximum_displacement=max_displacement,
    )
    sim.set_integrator(method=displacement_capped, dt=dt)
    sim.run(n_steps)
    # Run with LJ pairs, original bonds and angles
    if n_relax_steps > 0:
        if dpd_A:  # Repalce DPD with original LJ for final relaxation
            sim.active_forces.remove(dpd)
            sim.active_forces.append(lj)
        if bond_k_scale != 1.0:
            for param in bond.params:
                bond.params[param]["k"] /= bond_k_scale
        if angle_k_scale != 0:  # Angle force already added, rescale
            if angle_k_scale != 1.0:
                for param in angle.params:
                    angle.params[param]["k"] /= angle_k_scale
        else:  # Don't rescale, just add the original angle force
            sim.active_forces.append(angle)
        sim._update_integrator_forces()
        sim.run(n_relax_steps)
    # Update particle positions, save latest state point snapshot
    sim.update_positions()
    sim._update_snapshot()


def hoomd_fire(
    compound,
    forcefield,
    run_on_gpu,
    fire_iteration_steps,
    num_fire_iterations=1,
    n_relax_steps=1000,
    fixed_compounds=None,
    integrate_compounds=None,
    dt=1e-5,
    dpd_A=10,
    integrate_dof=False,
    min_steps_adapt=5,
    min_steps_conv=100,
    finc_dt=1.1,
    fdec_dt=0.5,
    alpha_start=0.1,
    fdec_alpha=0.95,
    force_tol=1e-2,
    angmom_tol=1e-2,
    energy_tol=1e-6,
    seed=42,
    r_cut=1.0,
    bond_k_scale=100,
    angle_k_scale=100,
    gsd_file=None,
):
    """Run a short HOOMD-Blue simulation with the FIRE integrator.
    This method can be helpful for relaxing poor molecular
    geometries or for removing overlapping particles.

    Parameters:
    -----------
    compound : mbuild.compound.Compound; required
        The compound to perform the displacement capped simulation with.
    forcefield : foyer.focefield.ForceField or gmso.core.Forcefield; required
        The forcefield used for the simulation.
    steps : int; required
        The number of simulation steps to run
    max_displacement : float, default 1e-3
        The value of the maximum displacement (nm)
    run_on_gpu : bool, default True
        If true, attempts to run HOOMD-Blue on a GPU.
        If a GPU device isn't found, then it will run on the CPU

    Notes
    -----
    Running on GPU:
        The hoomd-based methods in ``mbuild.simulation``
        are able to utilize and run on GPUs to perform energy minimization for larger systems if needed.
        Performance improvements of using GPU over CPU may not be significant until systems reach a
        size of ~1,000 particles.

    Running multiple simulations:
        The information needed to run the HOOMD simulation is saved after the first simulation.
        Therefore, calling hoomd-based simulation methods multiple times (e.g, in a for loop or while loop)
        does not require re-running performing atom-typing, applying the forcefield, or running hoomd
        format writers again.
    """
    compound._kick()
    sim = HoomdSimulation(
        compound=compound,
        forcefield=forcefield,
        r_cut=r_cut,
        run_on_gpu=run_on_gpu,
        seed=seed,
        integrate_compounds=integrate_compounds,
        fixed_compounds=fixed_compounds,
    )
    bond = sim.get_force(hoomd.md.bond.Harmonic)
    angle = sim.get_force(hoomd.md.angle.Harmonic)
    lj = sim.get_force(hoomd.md.pair.LJ)
    # Scale bond K and angle K
    if bond_k_scale != 1.0:
        for param in bond.params:
            bond.params[param]["k"] *= bond_k_scale
    sim.active_forces.append(bond)
    # If angle_k_scale is zero, skip and don't append angle force
    if angle_k_scale != 0:
        if angle_k_scale != 1.0:
            for param in angle.params:
                angle.params[param]["k"] *= angle_k_scale
        sim.active_forces.append(angle)
    if dpd_A:
        dpd = sim.get_dpd_from_lj(A=dpd_A)
        sim.active_forces.append(dpd)
    else:
        sim.active_forces.append(lj)
    # Set up and run
    nvt = hoomd.md.methods.ConstantVolume(filter=sim.get_integrate_group())
    sim.set_fire_integrator(
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
        methods=[nvt],
    )
    # Run FIRE sims with DPD + scaled bonds and angles
    for sim_num in range(num_fire_iterations):
        sim.run(fire_iteration_steps)
        sim.operations.integrator.reset()

    # Re-scale bonds and angle force constants
    if n_relax_steps > 0:
        if dpd_A:  # Replace DPD pair force with original LJ
            sim.active_forces.remove(dpd)
            sim.active_forces.append(lj)
        if bond_k_scale != 1.0:
            for param in bond.params:
                bond.params[param]["k"] /= bond_k_scale
        if angle_k_scale != 0:  # Angle force already included, rescale.
            if angle_k_scale != 1.0:
                for param in angle.params:
                    angle.params[param]["k"] /= angle_k_scale
        else:  # Don't recale, just add angles for final relaxation step
            sim.active_forces.append(angle)
        sim._update_integrator_forces()
        sim.run(n_relax_steps)
    # Update particle positions, save latest state point snapshot
    sim._update_snapshot()
    sim.update_positions()


# Openbabel and OpenMM
def energy_minimize(
    compound,
    forcefield="UFF",
    steps=1000,
    shift_com=True,
    anchor=None,
    **kwargs,
):
    """Perform an energy minimization on a Compound.

    Default behavior utilizes `Open Babel <http://openbabel.org/docs/dev/>`_
    to perform an energy minimization/geometry optimization on a Compound by
    applying a generic force field

    Can also utilize `OpenMM <http://openmm.org/>`_ to energy minimize after
    atomtyping a Compound using
    `Foyer <https://github.com/mosdef-hub/foyer>`_ to apply a forcefield XML
    file that contains valid SMARTS strings.

    This function is primarily intended to be used on smaller components,
    with sizes on the order of 10's to 100's of particles, as the energy
    minimization scales poorly with the number of particles.

    Parameters
    ----------
    compound : mbuid.Compound, required
        The compound to perform energy minimization on.
    steps : int, optional, default=1000
        The number of optimization iterations
    forcefield : str, optional, default='UFF'
        The generic force field to apply to the Compound for minimization.
        Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', 'Ghemical'.
        Please refer to the `Open Babel documentation
        <http://open-babel.readthedocs.io/en/latest/Forcefields/Overview.html>`_
        when considering your choice of force field.
        Utilizing OpenMM for energy minimization requires a forcefield
        XML file with valid SMARTS strings. Please refer to `OpenMM docs
        <http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields>`_
        for more information.
    shift_com : bool, optional, default=True
        If True, the energy-minimized Compound is translated such that the
        center-of-mass is unchanged relative to the initial configuration.
    anchor : Compound, optional, default=None
        Translates the energy-minimized Compound such that the
        position of the anchor Compound is unchanged relative to the
        initial configuration.

    Other Parameters
    ----------------
    algorithm : str, optional, default='cg'
        The energy minimization algorithm.  Valid options are 'steep', 'cg',
        and 'md', corresponding to steepest descent, conjugate gradient, and
        equilibrium molecular dynamics respectively.
        For _energy_minimize_openbabel
    fixed_compounds : Compound, optional, default=None
        An individual Compound or list of Compounds that will have their
        position fixed during energy minimization. Note, positions are fixed
        using a restraining potential and thus may change slightly.
        Position fixing will apply to all Particles (i.e., atoms) that exist
        in the Compound and to particles in any subsequent sub-Compounds.
        By default x,y, and z position is fixed. This can be toggled by instead
        passing a list containing the Compound and an list or tuple of bool values
        corresponding to x,y and z; e.g., [Compound, (True, True, False)]
        will fix the x and y position but allow z to be free.
        For _energy_minimize_openbabel
    ignore_compounds: Compound, optional, default=None
        An individual compound or list of Compounds whose underlying particles
        will have their positions fixed and not interact with other atoms via
        the specified force field during the energy minimization process.
        Note, a restraining potential used and thus absolute position may vary
        as a result of the energy minimization process.
        Interactions of these ignored atoms can  be specified by the user,
        e.g., by explicitly setting a distance constraint.
        For _energy_minimize_openbabel
    distance_constraints: list, optional, default=None
        A list containing a pair of Compounds as a tuple or list and
        a float value specifying the target distance between the two Compounds, e.g.,:
        [(compound1, compound2), distance].
        To specify more than one constraint, pass constraints as a 2D list, e.g.,:
        [ [(compound1, compound2), distance1],  [(compound3, compound4), distance2] ].
        Note, Compounds specified here must represent individual point particles.
        For _energy_minimize_openbabel
    constraint_factor: float, optional, default=50000.0
        Harmonic springs are used to constrain distances and fix atom positions, where
        the resulting energy associated with the spring is scaled by the
        constraint_factor; the energy of this spring is considering during the minimization.
        As such, very small values of the constraint_factor may result in an energy
        minimized state that does not adequately restrain the distance/position of atoms.
        For _energy_minimize_openbabel
    scale_bonds : float, optional, default=1
        Scales the bond force constant (1 is completely on).
        For _energy_minimize_openmm
    scale_angles : float, optional, default=1
        Scales the angle force constant (1 is completely on)
        For _energy_minimize_openmm
    scale_torsions : float, optional, default=1
        Scales the torsional force constants (1 is completely on)
        For _energy_minimize_openmm
        Note: Only Ryckaert-Bellemans style torsions are currently supported
    scale_nonbonded : float, optional, default=1
        Scales epsilon (1 is completely on)
        For _energy_minimize_openmm
    constraints : str, optional, default="AllBonds"
        Specify constraints on the molecule to minimize, options are:
        None, "HBonds", "AllBonds", "HAngles"
        For _energy_minimize_openmm

    References
    ----------
    If using _energy_minimize_openmm(), please cite:

    .. [Eastman2013] P. Eastman, M. S. Friedrichs, J. D. Chodera,
       R. J. Radmer, C. M. Bruns, J. P. Ku, K. A. Beauchamp, T. J. Lane,
       L.-P. Wang, D. Shukla, T. Tye, M. Houston, T. Stich, C. Klein,
       M. R. Shirts, and V. S. Pande. "OpenMM 4: A Reusable, Extensible,
       Hardware Independent Library for High Performance Molecular
       Simulation." J. Chem. Theor. Comput. 9(1): 461-469. (2013).

    If using _energy_minimize_openbabel(), please cite:

    .. [OBoyle2011] O'Boyle, N.M.; Banck, M.; James, C.A.; Morley, C.;
       Vandermeersch, T.; Hutchison, G.R. "Open Babel: An open chemical
       toolbox." (2011) J. Cheminf. 3, 33

    .. [OpenBabel] Open Babel, version X.X.X http://openbabel.org,
       (installed Month Year)

    If using the 'MMFF94' force field please also cite the following:

    .. [Halgren1996a] T.A. Halgren, "Merck molecular force field. I. Basis,
       form, scope, parameterization, and performance of MMFF94." (1996)
       J. Comput. Chem. 17, 490-519

    .. [Halgren1996b] T.A. Halgren, "Merck molecular force field. II. MMFF94
       van der Waals and electrostatic parameters for intermolecular
       interactions." (1996) J. Comput. Chem. 17, 520-552

    .. [Halgren1996c] T.A. Halgren, "Merck molecular force field. III.
       Molecular geometries and vibrational frequencies for MMFF94." (1996)
       J. Comput. Chem. 17, 553-586

    .. [Halgren1996d] T.A. Halgren and R.B. Nachbar, "Merck molecular force
       field. IV. Conformational energies and geometries for MMFF94." (1996)
       J. Comput. Chem. 17, 587-615

    .. [Halgren1996e] T.A. Halgren, "Merck molecular force field. V.
       Extension of MMFF94 using experimental data, additional computational
       data, and empirical rules." (1996) J. Comput. Chem. 17, 616-641

    If using the 'MMFF94s' force field please cite the above along with:

    .. [Halgren1999] T.A. Halgren, "MMFF VI. MMFF94s option for energy minimization
       studies." (1999) J. Comput. Chem. 20, 720-729

    If using the 'UFF' force field please cite the following:

    .. [Rappe1992] Rappe, A.K., Casewit, C.J., Colwell, K.S., Goddard, W.A.
       III, Skiff, W.M. "UFF, a full periodic table force field for
       molecular mechanics and molecular dynamics simulations." (1992)
       J. Am. Chem. Soc. 114, 10024-10039

    If using the 'GAFF' force field please cite the following:

    .. [Wang2004] Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A.,
       Case, D.A. "Development and testing of a general AMBER force field"
       (2004) J. Comput. Chem. 25, 1157-1174

    If using the 'Ghemical' force field please cite the following:

    .. [Hassinen2001] T. Hassinen and M. Perakyla, "New energy terms for
       reduced protein models implemented in an off-lattice force field"
       (2001) J. Comput. Chem. 22, 1229-1242

    """
    # TODO: Update mbuild tutorials to provide overview of new features
    #   Preliminary tutorials: https://github.com/chrisiacovella/mbuild_energy_minimization
    com = compound.pos
    anchor_in_compound = False
    if anchor is not None:
        # check to see if the anchor exists
        # in the Compound to be energy minimized
        for succesor in compound.successors():
            if id(anchor) == id(succesor):
                anchor_in_compound = True
                anchor_pos_old = anchor.pos

        if not anchor_in_compound:
            raise MBuildError(
                f"Anchor: {anchor} is not part of the Compound: {compound}"
                "that you are trying to energy minimize."
            )
    compound._kick()
    extension = os.path.splitext(forcefield)[-1]
    openbabel_ffs = ["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"]
    if forcefield in openbabel_ffs:
        _energy_minimize_openbabel(
            compound=compound, forcefield=forcefield, steps=steps, **kwargs
        )
    else:
        tmp_dir = tempfile.mkdtemp()
        compound.save(os.path.join(tmp_dir, "un-minimized.mol2"))

        if extension == ".xml":
            _energy_minimize_openmm(
                compound=compound,
                tmp_dir=tmp_dir,
                forcefield_files=forcefield,
                forcefield_name=None,
                steps=steps,
                **kwargs,
            )
        else:
            _energy_minimize_openmm(
                compound=compound,
                tmp_dir=tmp_dir,
                forcefield_files=None,
                forcefield_name=forcefield,
                steps=steps,
                **kwargs,
            )

        compound.update_coordinates(os.path.join(tmp_dir, "minimized.pdb"))

    if shift_com:
        compound.translate_to(com)

    if anchor_in_compound:
        anchor_pos_new = anchor.pos
        delta = anchor_pos_old - anchor_pos_new
        compound.translate(delta)


def _energy_minimize_openmm(
    compound,
    tmp_dir,
    forcefield_files=None,
    forcefield_name=None,
    steps=1000,
    scale_bonds=1,
    scale_angles=1,
    scale_torsions=1,
    scale_nonbonded=1,
    constraints="AllBonds",
):
    """Perform energy minimization using OpenMM.

    Converts an mBuild Compound to a ParmEd Structure,
    applies a forcefield using Foyer, and creates an OpenMM System.

    Parameters
    ----------
    compound : mbuid.Compound, required
        The compound to perform energy minimization on.
    forcefield_files : str or list of str, optional, default=None
        Forcefield files to load
    forcefield_name : str, optional, default=None
        Apply a named forcefield to the output file using the `foyer`
        package, e.g. 'oplsaa'. `Foyer forcefields`
        <https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields>_
    steps : int, optional, default=1000
        Number of energy minimization iterations
    scale_bonds : float, optional, default=1
        Scales the bond force constant (1 is completely on)
    scale_angles : float, optiona, default=1
        Scales the angle force constant (1 is completely on)
    scale_torsions : float, optional, default=1
        Scales the torsional force constants (1 is completely on)
    scale_nonbonded : float, optional, default=1
        Scales epsilon (1 is completely on)
    constraints : str, optional, default="AllBonds"
        Specify constraints on the molecule to minimize, options are:
        None, "HBonds", "AllBonds", "HAngles"

    Notes
    -----
    Assumes a particular organization for the force groups
    (HarmonicBondForce, HarmonicAngleForce, RBTorsionForce, NonBondedForce)

    References
    ----------
    [Eastman2013]_
    """
    foyer = import_("foyer")

    to_parmed = compound.to_parmed()
    ff = foyer.Forcefield(forcefield_files=forcefield_files, name=forcefield_name)
    to_parmed = ff.apply(to_parmed)

    import openmm.unit as u
    from openmm.app import AllBonds, HAngles, HBonds
    from openmm.app.pdbreporter import PDBReporter
    from openmm.app.simulation import Simulation
    from openmm.openmm import LangevinIntegrator

    if constraints:
        if constraints == "AllBonds":
            constraints = AllBonds
        elif constraints == "HBonds":
            constraints = HBonds
        elif constraints == "HAngles":
            constraints = HAngles
        else:
            raise ValueError(
                f"Provided constraints value of: {constraints}.\n"
                f'Expected "HAngles", "AllBonds" "HBonds".'
            )
        system = to_parmed.createSystem(
            constraints=constraints
        )  # Create an OpenMM System
    else:
        system = to_parmed.createSystem()  # Create an OpenMM System
    # Create a Langenvin Integrator in OpenMM
    integrator = LangevinIntegrator(
        298 * u.kelvin, 1 / u.picosecond, 0.002 * u.picoseconds
    )
    # Create Simulation object in OpenMM
    simulation = Simulation(to_parmed.topology, system, integrator)

    # Loop through forces in OpenMM System and set parameters
    for force in system.getForces():
        if type(force).__name__ == "HarmonicBondForce":
            for bond_index in range(force.getNumBonds()):
                atom1, atom2, r0, k = force.getBondParameters(bond_index)
                force.setBondParameters(bond_index, atom1, atom2, r0, k * scale_bonds)
            force.updateParametersInContext(simulation.context)

        elif type(force).__name__ == "HarmonicAngleForce":
            for angle_index in range(force.getNumAngles()):
                atom1, atom2, atom3, r0, k = force.getAngleParameters(angle_index)
                force.setAngleParameters(
                    angle_index, atom1, atom2, atom3, r0, k * scale_angles
                )
            force.updateParametersInContext(simulation.context)

        elif type(force).__name__ == "RBTorsionForce":
            for torsion_index in range(force.getNumTorsions()):
                (
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    c0,
                    c1,
                    c2,
                    c3,
                    c4,
                    c5,
                ) = force.getTorsionParameters(torsion_index)
                force.setTorsionParameters(
                    torsion_index,
                    atom1,
                    atom2,
                    atom3,
                    atom4,
                    c0 * scale_torsions,
                    c1 * scale_torsions,
                    c2 * scale_torsions,
                    c3 * scale_torsions,
                    c4 * scale_torsions,
                    c5 * scale_torsions,
                )
            force.updateParametersInContext(simulation.context)

        elif type(force).__name__ == "NonbondedForce":
            for nb_index in range(force.getNumParticles()):
                charge, sigma, epsilon = force.getParticleParameters(nb_index)
                force.setParticleParameters(
                    nb_index, charge, sigma, epsilon * scale_nonbonded
                )
            force.updateParametersInContext(simulation.context)

        elif type(force).__name__ == "CMMotionRemover":
            pass

        else:
            warn(
                f"OpenMM Force {type(force).__name__} is "
                "not currently supported in _energy_minimize_openmm. "
                "This Force will not be updated!"
            )

    simulation.context.setPositions(to_parmed.positions)
    # Run energy minimization through OpenMM
    simulation.minimizeEnergy(maxIterations=steps)
    reporter = PDBReporter(os.path.join(tmp_dir, "minimized.pdb"), 1)
    reporter.report(simulation, simulation.context.getState(getPositions=True))


def _check_openbabel_constraints(
    compound,
    particle_list,
    successors_list,
    check_if_particle=False,
):
    """Provide routines commonly used to check constraint inputs."""
    for part in particle_list:
        if not isinstance(part, Compound):
            raise MBuildError(f"{part} is not a Compound.")
        if id(part) != id(compound) and id(part) not in successors_list:
            raise MBuildError(f"{part} is not a member of Compound {compound}.")

        if check_if_particle:
            if len(part.children) != 0:
                raise MBuildError(
                    f"{part} does not correspond to an individual particle."
                )


def _energy_minimize_openbabel(
    compound,
    steps=1000,
    algorithm="cg",
    forcefield="UFF",
    constraint_factor=50000.0,
    distance_constraints=None,
    fixed_compounds=None,
    ignore_compounds=None,
):
    """Perform an energy minimization on a Compound.

    Utilizes Open Babel (http://openbabel.org/docs/dev/) to perform an
    energy minimization/geometry optimization on a Compound by applying
    a generic force field.

    This function is primarily intended to be used on smaller components,
    with sizes on the order of 10's to 100's of particles, as the energy
    minimization scales poorly with the number of particles.

    Parameters
    ----------
    compound : mbuid.Compound, required
        The compound to perform energy minimization on.
    steps : int, optionl, default=1000
        The number of optimization iterations
    algorithm : str, optional, default='cg'
        The energy minimization algorithm.  Valid options are 'steep',
        'cg', and 'md', corresponding to steepest descent, conjugate
        gradient, and equilibrium molecular dynamics respectively.
    forcefield : str, optional, default='UFF'
        The generic force field to apply to the Compound for minimization.
        Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', 'Ghemical'.
        Please refer to the Open Babel documentation
        (http://open-babel.readthedocs.io/en/latest/Forcefields/Overview.html)
        when considering your choice of force field.
    fixed_compounds : Compound, optional, default=None
        An individual Compound or list of Compounds that will have their
        position fixed during energy minimization. Note, positions are fixed
        using a restraining potential and thus may change slightly.
        Position fixing will apply to all Particles (i.e., atoms) that exist
        in the Compound and to particles in any subsequent sub-Compounds.
        By default x,y, and z position is fixed. This can be toggled by instead
        passing a list containing the Compound and a list or tuple of bool values
        corresponding to x,y and z; e.g., [Compound, (True, True, False)]
        will fix the x and y position but allow z to be free.
    ignore_compounds: Compound, optional, default=None
        An individual compound or list of Compounds whose underlying particles
        will have their positions fixed and not interact with other atoms via
        the specified force field during the energy minimization process.
        Note, a restraining potential is used and thus absolute position may vary
        as a result of the energy minimization process.
        Interactions of these ignored atoms can  be specified by the user,
        e.g., by explicitly setting a distance constraint.
    distance_constraints: list, optional, default=None
        A list containing a pair of Compounds as a tuple or list and
        a float value specifying the target distance between the two Compounds, e.g.,:
        [(compound1, compound2), distance].
        To specify more than one constraint, pass constraints as a 2D list, e.g.,:
        [ [(compound1, compound2), distance1],  [(compound3, compound4), distance2] ].
        Note, Compounds specified here must represent individual point particles.
    constraint_factor: float, optional, default=50000.0
        Harmonic springs are used to constrain distances and fix atom positions, where
        the resulting energy associated with the spring is scaled by the
        constraint_factor; the energy of this spring is considering during the minimization.
        As such, very small values of the constraint_factor may result in an energy
        minimized state that does not adequately restrain the distance/position of atom(s)e.


    References
    ----------
    [OBoyle2011]_
    [OpenBabel]_

    If using the 'MMFF94' force field please also cite the following:
    [Halgren1996a]_
    [Halgren1996b]_
    [Halgren1996c]_
    [Halgren1996d]_
    [Halgren1996e]_

    If using the 'MMFF94s' force field please cite the above along with:
    [Halgren1999]_

    If using the 'UFF' force field please cite the following:
    [Rappe1992]_

    If using the 'GAFF' force field please cite the following:
    [Wang2001]_

    If using the 'Ghemical' force field please cite the following:
    [Hassinen2001]_
    """
    openbabel = import_("openbabel")
    for particle in compound.particles():
        if particle.element is None:
            try:
                particle._element = element_from_symbol(particle.name)
            except ElementError:
                try:
                    particle._element = element_from_name(particle.name)
                except ElementError:
                    raise MBuildError(
                        f"No element assigned to {particle}; element could not be"
                        f"inferred from particle name {particle.name}. Cannot perform"
                        "an energy minimization."
                    )
    # Create a dict containing particle id and associated index to speed up looping
    particle_idx = {
        id(particle): idx for idx, particle in enumerate(compound.particles())
    }

    # A list containing all Compounds ids contained in compound. Will be used to check if
    # compounds refered to in the constrains are actually in the Compound we are minimizing.
    successors_list = [id(comp) for comp in compound.successors()]

    # initialize constraints
    ob_constraints = openbabel.OBFFConstraints()

    if distance_constraints is not None:
        # if a user passes single constraint as a 1-D array,
        # i.e., [(p1,p2), 2.0]  rather than [[(p1,p2), 2.0]],
        # just add it to a list so we can use the same looping code
        if len(np.array(distance_constraints, dtype=object).shape) == 1:
            distance_constraints = [distance_constraints]

        for con_temp in distance_constraints:
            p1 = con_temp[0][0]
            p2 = con_temp[0][1]

            _check_openbabel_constraints(
                compound=compound,
                particle_list=[p1, p2],
                successors_list=successors_list,
                check_if_particle=True,
            )
            if id(p1) == id(p2):
                raise MBuildError(
                    f"Cannot create a constraint between a Particle and itself: {p1} {p2} ."
                )

            # openbabel indices start at 1
            pid_1 = particle_idx[id(p1)] + 1
            # openbabel indices start at 1
            pid_2 = particle_idx[id(p2)] + 1
            dist = (
                con_temp[1] * 10.0
            )  # obenbabel uses angstroms, not nm, convert to angstroms

            ob_constraints.AddDistanceConstraint(pid_1, pid_2, dist)

    if fixed_compounds is not None:
        # if we are just passed a single Compound, wrap it into
        # and array so we can just use the same looping code
        if isinstance(fixed_compounds, Compound):
            fixed_compounds = [fixed_compounds]

        # if fixed_compounds is a 1-d array and it is of length 2, we need to determine whether it is
        # a list of two Compounds or if fixed_compounds[1] should correspond to the directions to constrain
        if len(np.array(fixed_compounds, dtype=object).shape) == 1:
            if len(fixed_compounds) == 2:
                if not isinstance(fixed_compounds[1], Compound):
                    # if it is not a list of two Compounds, make a 2d array so we can use the same looping code
                    fixed_compounds = [fixed_compounds]

        for fixed_temp in fixed_compounds:
            # if an individual entry is a list, validate the input
            if isinstance(fixed_temp, list):
                if len(fixed_temp) == 2:
                    msg1 = (
                        "Expected tuple or list of length 3 to set"
                        "which dimensions to fix motion."
                    )
                    assert isinstance(fixed_temp[1], (list, tuple)), msg1

                    msg2 = (
                        "Expected tuple or list of length 3 to set"
                        "which dimensions to fix motion, "
                        f"{len(fixed_temp[1])} found."
                    )
                    assert len(fixed_temp[1]) == 3, msg2

                    dims = [dim for dim in fixed_temp[1]]
                    msg3 = (
                        "Expected bool values for which directions are fixed."
                        f"Found instead {dims}."
                    )
                    assert all(isinstance(dim, bool) for dim in dims), msg3

                    p1 = fixed_temp[0]

                # if fixed_compounds is defined as [[Compound],[Compound]],
                # fixed_temp will be a list of length 1
                elif len(fixed_temp) == 1:
                    p1 = fixed_temp[0]
                    dims = [True, True, True]

            else:
                p1 = fixed_temp
                dims = [True, True, True]

            all_true = all(dims)

            _check_openbabel_constraints(
                compound=compound, particle_list=[p1], successors_list=successors_list
            )

            if len(p1.children) == 0:
                pid = particle_idx[id(p1)] + 1  # openbabel indices start at 1

                if all_true:
                    ob_constraints.AddAtomConstraint(pid)
                else:
                    if dims[0]:
                        ob_constraints.AddAtomXConstraint(pid)
                    if dims[1]:
                        ob_constraints.AddAtomYConstraint(pid)
                    if dims[2]:
                        ob_constraints.AddAtomZConstraint(pid)
            else:
                for particle in p1.particles():
                    pid = particle_idx[id(particle)] + 1  # openbabel indices start at 1

                    if all_true:
                        ob_constraints.AddAtomConstraint(pid)
                    else:
                        if dims[0]:
                            ob_constraints.AddAtomXConstraint(pid)
                        if dims[1]:
                            ob_constraints.AddAtomYConstraint(pid)
                        if dims[2]:
                            ob_constraints.AddAtomZConstraint(pid)

    if ignore_compounds is not None:
        temp1 = np.array(ignore_compounds, dtype=object)
        if len(temp1.shape) == 2:
            ignore_compounds = list(temp1.reshape(-1))

        # Since the ignore_compounds can only be passed as a list
        # we can check the whole list at once before looping over it
        _check_openbabel_constraints(
            compound=compound,
            particle_list=ignore_compounds,
            successors_list=successors_list,
        )

        for ignore in ignore_compounds:
            p1 = ignore
            if len(p1.children) == 0:
                pid = particle_idx[id(p1)] + 1  # openbabel indices start at 1
                ob_constraints.AddIgnore(pid)

            else:
                for particle in p1.particles():
                    pid = particle_idx[id(particle)] + 1  # openbabel indices start at 1
                    ob_constraints.AddIgnore(pid)

    mol = compound.to_pybel()
    mol = mol.OBMol

    mol.PerceiveBondOrders()
    mol.SetAtomTypesPerceived()

    ff = openbabel.OBForceField.FindForceField(forcefield)
    if ff is None:
        raise MBuildError(
            f"Force field '{forcefield}' not supported for energy "
            "minimization. Valid force fields are 'MMFF94', "
            "'MMFF94s', 'UFF', 'GAFF', and 'Ghemical'."
            ""
        )
    warn(
        "Performing energy minimization using the Open Babel package. "
        "Please refer to the documentation to find the appropriate "
        f"citations for Open Babel and the {forcefield} force field"
    )

    if (
        distance_constraints is not None
        or fixed_compounds is not None
        or ignore_compounds is not None
    ):
        ob_constraints.SetFactor(constraint_factor)
        if ff.Setup(mol, ob_constraints) == 0:
            raise MBuildError("Could not setup forcefield for OpenBabel Optimization.")
    else:
        if ff.Setup(mol) == 0:
            raise MBuildError("Could not setup forcefield for OpenBabel Optimization.")

    if algorithm == "steep":
        ff.SteepestDescent(steps)
    elif algorithm == "md":
        ff.MolecularDynamicsTakeNSteps(steps, 300)
    elif algorithm == "cg":
        ff.ConjugateGradients(steps)
    else:
        raise MBuildError(
            "Invalid minimization algorithm. Valid options are 'steep', 'cg', and 'md'."
        )
    ff.UpdateCoordinates(mol)

    # update the coordinates in the Compound
    for i, obatom in enumerate(openbabel.OBMolAtomIter(mol)):
        x = obatom.GetX() / 10.0
        y = obatom.GetY() / 10.0
        z = obatom.GetZ() / 10.0
        compound[i].pos = np.array([x, y, z])
