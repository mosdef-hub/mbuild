"""HOOMD simulation format."""

import itertools
import operator
import warnings
from collections import namedtuple

import numpy as np
import parmed as pmd

import mbuild as mb
from mbuild.utils.conversion import RB_to_OPLS
from mbuild.utils.io import import_
from mbuild.utils.sorting import natural_sort
from mbuild.utils.decorators import deprecated

from .hoomd_snapshot import to_hoomdsnapshot

gsd = import_("gsd")
gsd.hoomd = import_("gsd.hoomd")
hoomd = import_("hoomd")
hoomd.md = import_("hoomd.md")


dep_msg = """
Support for Hoomd-Blue 2.x will be removed in mBuild 1.0.
See GMSO (https://github.com/mosdef-hub/gmso) for 
continued Hoomd-Blue 3.x and 4.x support.
"""
print(dep_msg)


@deprecated(dep_msg)
def create_hoomd_simulation(
    structure,
    r_cut,
    ref_distance=1.0,
    ref_mass=1.0,
    ref_energy=1.0,
    auto_scale=False,
    snapshot_kwargs={},
    pppm_kwargs={"Nx": 8, "Ny": 8, "Nz": 8, "order": 4},
    init_snap=None,
    restart=None,
    nlist=hoomd.md.nlist.cell,
):
    """Convert a parametrized pmd.Structure to hoomd.SimulationContext.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd Structure object
    r_cut : float
        Cutoff radius in simulation units
    ref_distance : float, optional, default=1.0
        Reference distance for unit conversion (from Angstrom)
    ref_mass : float, optional, default=1.0
        Reference mass for unit conversion (from Dalton)
    ref_energy : float, optional, default=1.0
        Reference energy for unit conversion (from kcal/mol)
    auto_scale : bool, optional, default=False
        Scale to reduced units by automatically using the largest sigma value
        as ref_distance, largest mass value as ref_mass, and largest epsilon
        value as ref_energy
    snapshot_kwargs : dict
        Kwargs to pass to to_hoomdsnapshot
    pppm_kwargs : dict
        Kwargs to pass to hoomd's pppm function
    init_snap : hoomd.data.SnapshotParticleData, optional, default=None
        Initial snapshot to which to add the ParmEd structure object
        (useful for rigid bodies)
    restart : str, optional, default=None
        Path to the gsd file from which to restart the simulation.
        https://hoomd-blue.readthedocs.io/en/v2.9.4/restartable-jobs.html
        Note: It is assumed that the ParmEd structure and the system in
        restart.gsd contain the same types. The ParmEd structure is still used
        to initialize the forces, but restart.gsd is used to initialize the
        system state (e.g., particle positions, momenta, etc).
    nlist : hoomd.md.nlist, default=hoomd.md.nlist.cell
        Type of neighborlist to use, see
        https://hoomd-blue.readthedocs.io/en/stable/module-md-nlist.html
        for more information.

    Returns
    -------
    hoomd_objects : list
        List of hoomd objects created during conversion
    ReferenceValues : namedtuple
        Values used in scaling

    Note
    ----
    While the hoomd objects are returned, the hoomd.SimulationContext is
    accessible via `hoomd.context.current`. If you pass a non-parametrized
    pmd.Structure, you will not have angle, dihedral, or force field
    information. You may be better off creating a hoomd.Snapshot.
    Reference units should be expected to convert parmed Structure units:

    --- angstroms, kcal/mol, and daltons
    """
    if isinstance(structure, mb.Compound):
        raise ValueError(
            "You passed mb.Compound to create_hoomd_simulation, there will be "
            "no angles, dihedrals, or force field parameters. Please use "
            "hoomd_snapshot.to_hoomdsnapshot to create a hoomd.Snapshot, then "
            "create your own hoomd context and pass your hoomd.Snapshot to "
            "hoomd.init.read_snapshot()"
        )
    elif not isinstance(structure, pmd.Structure):
        raise ValueError(
            "Please pass a parmed.Structure to create_hoomd_simulation"
        )

    _check_hoomd_version()
    version_numbers = _check_hoomd_version()
    if float(version_numbers[0]) >= 3:
        raise RuntimeError("This method does not support HOOMD-blue v3.x.")

    hoomd_objects = []  # Potential adaptation for Hoomd v3 API

    if auto_scale:
        if not all([i == 1 for i in (ref_distance, ref_energy, ref_mass)]):
            warnings.warn(
                "Autoscale option selected--provided reference values will not "
                "be used."
            )
        pair_coeffs = list(
            set((a.type, a.epsilon, a.sigma) for a in structure.atoms)
        )

        ref_mass = max([atom.mass for atom in structure.atoms])
        ref_energy = max(pair_coeffs, key=operator.itemgetter(1))[1]
        ref_distance = max(pair_coeffs, key=operator.itemgetter(2))[2]

    ReferenceValues = namedtuple("ref_values", ["distance", "mass", "energy"])
    ref_values = ReferenceValues(ref_distance, ref_mass, ref_energy)

    if not hoomd.context.current:
        hoomd.context.initialize("")
    if restart is None:
        snapshot, _ = to_hoomdsnapshot(
            structure,
            ref_distance=ref_distance,
            ref_mass=ref_mass,
            ref_energy=ref_energy,
            **snapshot_kwargs,
            hoomd_snapshot=init_snap,
        )
        hoomd_objects.append(snapshot)
        hoomd_system = hoomd.init.read_snapshot(snapshot)
        hoomd_objects.append(hoomd_system)
    else:
        with gsd.hoomd.open(restart) as f:
            snapshot = f[-1]
        hoomd_objects.append(snapshot)
        hoomd_system = hoomd.init.read_gsd(restart, restart=restart)
        hoomd_objects.append(hoomd_system)
        print("Simulation initialized from restart file")

    nl = nlist()
    nl.reset_exclusions(exclusions=["1-2", "1-3"])
    hoomd_objects.append(nl)

    if structure.atoms[0].type != "":
        print("Processing LJ and QQ")
        lj = _init_hoomd_lj(
            structure,
            nl,
            r_cut,
            ref_distance=ref_distance,
            ref_energy=ref_energy,
        )
        qq = _init_hoomd_qq(structure, nl, r_cut, **pppm_kwargs)
        hoomd_objects.append(lj)
        hoomd_objects.append(qq)
    if structure.adjusts:
        print("Processing 1-4 interactions, adjusting neighborlist exclusions")
        lj_14, qq_14 = _init_hoomd_14_pairs(
            structure,
            nl,
            r_cut,
            ref_distance=ref_distance,
            ref_energy=ref_energy,
        )
        hoomd_objects.append(lj_14)
        hoomd_objects.append(qq_14)
    if structure.bond_types:
        print("Processing harmonic bonds")
        harmonic_bond = _init_hoomd_bonds(
            structure, ref_distance=ref_distance, ref_energy=ref_energy
        )
        hoomd_objects.append(harmonic_bond)
    if structure.angle_types:
        print("Processing harmonic angles")
        harmonic_angle = _init_hoomd_angles(structure, ref_energy=ref_energy)
        hoomd_objects.append(harmonic_angle)
    if structure.dihedral_types:
        print("Processing periodic torsions")
        periodic_torsions = _init_hoomd_dihedrals(
            structure, ref_energy=ref_energy
        )
        hoomd_objects.append(periodic_torsions)
    if structure.rb_torsion_types:
        print("Processing RB torsions")
        rb_torsions = _init_hoomd_rb_torsions(structure, ref_energy=ref_energy)
        hoomd_objects.append(rb_torsions)
    print("HOOMD SimulationContext updated from ParmEd Structure")

    return hoomd_objects, ref_values


def _init_hoomd_lj(structure, nl, r_cut, ref_distance=1.0, ref_energy=1.0):
    """LJ parameters."""
    # Identify the unique atom types before setting
    atom_type_params = {}
    for atom in structure.atoms:
        if atom.type not in atom_type_params:
            atom_type_params[atom.type] = atom.atom_type

    # Set the hoomd parameters for self-interactions
    lj = hoomd.md.pair.lj(r_cut, nl)
    for name, atom_type in atom_type_params.items():
        lj.pair_coeff.set(
            name,
            name,
            sigma=atom_type.sigma / ref_distance,
            epsilon=atom_type.epsilon / ref_energy,
        )
        if atom_type.epsilon / ref_energy == 0:
            lj.pair_coeff.set(name, name, r_cut=0)

    # Cross interactions, mixing rules, NBfixes
    all_atomtypes = sorted(atom_type_params.keys())
    for a1, a2 in itertools.combinations_with_replacement(all_atomtypes, 2):
        nb_fix_info = atom_type_params[a1].nbfix.get(a2, None)
        # nb_fix_info = (rmin, eps, rmin14, eps14)
        if nb_fix_info is None:
            # No nbfix means use mixing rule to find cross-interaction
            if structure.combining_rule == "lorentz":
                sigma = (
                    atom_type_params[a1].sigma + atom_type_params[a2].sigma
                ) / (2 * ref_distance)
                epsilon = (
                    (
                        atom_type_params[a1].epsilon
                        * atom_type_params[a2].epsilon
                    )
                    / ref_energy**2
                ) ** 0.5
            elif structure.combining_rule == "geometric":
                sigma = (
                    (atom_type_params[a1].sigma * atom_type_params[a2].sigma)
                    / ref_distance**2
                ) ** 0.5
                epsilon = (
                    (
                        atom_type_params[a1].epsilon
                        * atom_type_params[a2].epsilon
                    )
                    / ref_energy**2
                ) ** 0.5
            else:
                raise ValueError(
                    f"Mixing rule {structure.combining_rule} not supported, "
                    "use lorentz"
                )
        else:
            # If we have nbfix info, use it
            sigma = nb_fix_info[0] / (ref_distance * (2 ** (1 / 6)))
            epsilon = nb_fix_info[1] / ref_energy
        lj.pair_coeff.set(a1, a2, sigma=sigma, epsilon=epsilon)
        if epsilon == 0:
            lj.pair_coeff.set(a1, a2, r_cut=0)

    return lj


def _init_hoomd_qq(structure, nl, r_cut, Nx=1, Ny=1, Nz=1, order=4):
    """Charge interactions."""
    charged = hoomd.group.charged()
    if len(charged) == 0:
        print("No charged groups found, ignoring electrostatics")
        return None
    else:
        qq = hoomd.md.charge.pppm(charged, nl)
        qq.set_params(Nx, Ny, Nz, order, r_cut)
        return qq


def _init_hoomd_14_pairs(
    structure, nl, r_cut, ref_distance=1.0, ref_energy=1.0
):
    """Special_pairs to handle 14 scalings.

    See discussion: https://groups.google.com/forum/#!topic/hoomd-users/
    iZ9WCpHczg0
    """
    # Update neighborlist to exclude 1-4 interactions,
    # but impose a special_pair force to handle these pairs
    nl.exclusions.append("1-4")

    if hoomd.context.current.system_definition.getPairData().getN() == 0:
        print("No 1,4 pairs found in hoomd snapshot")
        return None, None

    lj_14 = hoomd.md.special_pair.lj()
    qq_14 = hoomd.md.special_pair.coulomb()
    params_14 = {}
    # Identify unique 14 scalings
    for adjust in structure.adjusts:
        t1 = adjust.atom1.type
        t2 = adjust.atom2.type
        ps = "-".join(sorted([t1, t2]))
        if ps not in params_14:
            params_14[ps] = adjust.type
    for name, adjust_type in params_14.items():
        lj_14.pair_coeff.set(
            name,
            sigma=adjust_type.sigma / ref_distance,
            # The adjust epsilon already carries the scaling
            epsilon=adjust_type.epsilon / ref_energy,
            # Do NOT use hoomd's alpha to modify any LJ terms
            alpha=1,
            r_cut=r_cut,
        )
        qq_14.pair_coeff.set(name, alpha=adjust_type.chgscale, r_cut=r_cut)

    return lj_14, qq_14


def _init_hoomd_bonds(structure, ref_distance=1.0, ref_energy=1.0):
    """Harmonic bonds."""
    # Identify the unique bond types before setting
    bond_type_params = {}
    for bond in structure.bonds:
        t1, t2 = bond.atom1.type, bond.atom2.type
        t1, t2 = sorted([t1, t2], key=natural_sort)
        if t1 != "" and t2 != "":
            bond_type = "-".join((t1, t2))
            if bond_type not in bond_type_params:
                bond_type_params[bond_type] = bond.type

    # Set the hoomd parameters
    harmonic_bond = hoomd.md.bond.harmonic()
    for name, bond_type in bond_type_params.items():
        # A (paramerized) parmed structure with no bondtype
        # is because of constraints
        if bond_type is None:
            print("Bond with no bondtype detected, setting coefficients to 0")
            harmonic_bond.bond_coeff.set(name, k=0, r0=0)
        else:
            harmonic_bond.bond_coeff.set(
                name,
                k=2 * bond_type.k * ref_distance**2 / ref_energy,
                r0=bond_type.req / ref_distance,
            )

    return harmonic_bond


def _init_hoomd_angles(structure, ref_energy=1.0):
    """Harmonic angles."""
    # Identify the unique angle types before setting
    angle_type_params = {}
    for angle in structure.angles:
        t1, t2, t3 = angle.atom1.type, angle.atom2.type, angle.atom3.type
        t1, t3 = sorted([t1, t3], key=natural_sort)
        angle_type = "-".join((t1, t2, t3))
        if angle_type not in angle_type_params:
            angle_type_params[angle_type] = angle.type

    # set the hoomd parameters
    harmonic_angle = hoomd.md.angle.harmonic()
    for name, angle_type in angle_type_params.items():
        harmonic_angle.angle_coeff.set(
            name,
            t0=np.deg2rad(angle_type.theteq),
            k=2 * angle_type.k / ref_energy,
        )

    return harmonic_angle


def _init_hoomd_dihedrals(structure, ref_energy=1.0):
    """Periodic dihedrals (dubbed harmonic dihedrals in HOOMD)."""
    # Identify the unique dihedral types before setting
    # need Hoomd 2.8.0 to use proper dihedral implemtnation
    # from this PR https://github.com/glotzerlab/hoomd-blue/pull/492
    version_numbers = _check_hoomd_version()
    if float(version_numbers[0]) < 2 or float(version_numbers[1]) < 8:
        from mbuild.exceptions import MBuildError

        raise MBuildError("Please upgrade Hoomd to at least 2.8.0")

    dihedral_type_params = {}
    for dihedral in structure.dihedrals:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        if dihedral_type not in dihedral_type_params:
            if isinstance(dihedral.type, pmd.DihedralType):
                dihedral_type_params[dihedral_type] = dihedral.type
            elif isinstance(dihedral.type, pmd.DihedralTypeList):
                if len(dihedral.type) > 1:
                    warnings.warn(
                        "Multiple dihedral types detected for single dihedral, "
                        "will ignore all except  first dihedral type. First "
                        "dihedral type: {}".format(dihedral.type[0])
                    )
                dihedral_type_params[dihedral_type] = dihedral.type[0]

    # Set the hoomd parameters
    # These are periodic torsions
    periodic_torsion = hoomd.md.dihedral.harmonic()
    for name, dihedral_type in dihedral_type_params.items():
        periodic_torsion.dihedral_coeff.set(
            name,
            k=2 * dihedral_type.phi_k / ref_energy,
            d=1,
            n=dihedral_type.per,
            phi_0=np.deg2rad(dihedral_type.phase),
        )

    return periodic_torsion


def _init_hoomd_rb_torsions(structure, ref_energy=1.0):
    """RB dihedrals (implemented as OPLS dihedrals in HOOMD)."""
    # Identify the unique dihedral types before setting
    dihedral_type_params = {}
    for dihedral in structure.rb_torsions:
        t1, t2 = dihedral.atom1.type, dihedral.atom2.type
        t3, t4 = dihedral.atom3.type, dihedral.atom4.type
        if [t2, t3] == sorted([t2, t3], key=natural_sort):
            dihedral_type = "-".join((t1, t2, t3, t4))
        else:
            dihedral_type = "-".join((t4, t3, t2, t1))
        if dihedral_type not in dihedral_type_params:
            dihedral_type_params[dihedral_type] = dihedral.type

    # Set the hoomd parameter
    rb_torsion = hoomd.md.dihedral.opls()
    for name, dihedral_type in dihedral_type_params.items():
        F_coeffs = RB_to_OPLS(
            dihedral_type.c0 / ref_energy,
            dihedral_type.c1 / ref_energy,
            dihedral_type.c2 / ref_energy,
            dihedral_type.c3 / ref_energy,
            dihedral_type.c4 / ref_energy,
            dihedral_type.c5 / ref_energy,
            error_if_outside_tolerance=False,
        )
        rb_torsion.dihedral_coeff.set(
            name, k1=F_coeffs[1], k2=F_coeffs[2], k3=F_coeffs[3], k4=F_coeffs[4]
        )

    return rb_torsion


def _check_hoomd_version():
    version = hoomd.__version__
    version_numbers = version.split(".")
    return version_numbers
