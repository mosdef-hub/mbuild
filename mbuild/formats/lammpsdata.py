"""LAMMPS data format."""
import itertools as it
from collections import OrderedDict
from warnings import warn

import numpy as np
from parmed.parameters import ParameterSet
from parmed import Structure
from scipy.constants import epsilon_0

from mbuild import Box
from mbuild.utils.conversion import RB_to_OPLS
from mbuild.utils.sorting import natural_sort

__all__ = ["write_lammpsdata"]

# returns True if both mins and maxs have been defined, and each have length 3
# otherwise returns False
def _check_minsmaxs(mins, maxs):
    if mins and maxs:
        if len(mins) == 3 and len(maxs) == 3:
            return True
        else:
            warn(
                "mins and maxs passed to write_lammpsdata, but list size is "
                "incorrect. mins and maxs will be ignored."
            )
            return False
    else:
        return False


def write_lammpsdata(
    structure,
    filename,
    atom_style="full",
    unit_style="real",
    mins=None,
    maxs=None,
    pair_coeff_label=None,
    detect_forcefield_style=True,
    nbfix_in_data_file=True,
    use_urey_bradleys=False,
    use_rb_torsions=True,
    use_dihedrals=False,
    zero_dihedral_weighting_factor=False,
    moleculeID_offset=1,
):
    """Output a LAMMPS data file.

    Outputs a LAMMPS data file in the 'full' atom style format. Default units
    are 'real' units. See http://lammps.sandia.gov/doc/atom_style.html for
    more information on atom styles.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file
    atom_style: str
        Defines the style of atoms to be saved in a LAMMPS data file. The
        following atom styles are currently supported:
        'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more information
        on atom styles.
    unit_style: str
        Defines to unit style to be save in a LAMMPS data file.  Defaults to
        'real' units. Current styles are supported: 'real', 'lj'
        see https://lammps.sandia.gov/doc/99/units.html for more information
        on unit styles
    mins : list
        minimum box dimension in x, y, z directions
    maxs : list
        maximum box dimension in x, y, z directions
    pair_coeff_label : str
        Provide a custom label to the pair_coeffs section in the lammps data
        file. Defaults to None, which means a suitable default will be chosen.
    detect_forcefield_style: boolean
        If True, format lammpsdata parameters based on the contents of
        the parmed Structure
    use_urey_bradleys: boolean
        If True, will treat angles as CHARMM-style angles with urey bradley
        terms while looking for `structure.urey_bradleys`
    use_rb_torsions:
        If True, will treat dihedrals OPLS-style torsions while looking for
        `structure.rb_torsions`
    use_dihedrals:
        If True, will treat dihedrals as CHARMM-style dihedrals while looking
        for `structure.dihedrals`
    zero_dihedral_weighting_factor:
        If True, will set weighting parameter to zero in CHARMM-style dihedrals.
        This should be True if the CHARMM dihedral style is used in non-CHARMM forcefields.
    moleculeID_offset : int , optional, default=1
        Since LAMMPS treats the MoleculeID as an additional set of information
        to identify what molecule an atom belongs to, this currently
        behaves as a residue id. This value needs to start at 1 to be
        considered a real molecule.

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full
    description of the LAMMPS data format. Currently the following sections are
    supported (in addition to the header): *Masses*, *Nonbond Coeffs*,
    *Bond Coeffs*, *Angle Coeffs*, *Dihedral Coeffs*, *Atoms*, *Bonds*,
    *Angles*, *Dihedrals*, *Impropers*
    OPLS and CHARMM forcefield styles are supported, AMBER forcefield styles
    are NOT

    Some of this function has beed adopted from `mdtraj`'s support of the
    LAMMPSTRJ trajectory format. See
    https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py
    for details.

    unique_types : a sorted list of unique atomtypes for all atoms in the
        structure where atomtype = atom.type.
    unique_bond_types: an enumarated OrderedDict of unique bond types for all
        bonds in the structure.
        Defined by bond parameters and component atomtypes, in order:
        --- k : bond.type.k
        --- req : bond.type.req
        --- atomtypes : sorted((bond.atom1.type, bond.atom2.type))
    unique_angle_types: an enumerated OrderedDict of unique angle types for all
        angles in the structure.
        Defined by angle parameters and component atomtypes, in order:
        --- k : angle.type.k
        --- theteq : angle.type.theteq
        --- vertex atomtype: angle.atom2.type
        --- atomtypes: sorted((bond.atom1.type, bond.atom3.type))

    unique_dihedral_types: an enumerated OrderedDict of unique dihedrals type
        for all dihedrals in the structure.
        Defined by dihedral parameters and component atomtypes, in order:
        --- c0 : dihedral.type.c0
        --- c1 : dihedral.type.c1
        --- c2 : dihedral.type.c2
        --- c3 : dihedral.type.c3
        --- c4 : dihedral.type.c4
        --- c5 : dihedral.type.c5
        --- scee : dihedral.type.scee
        --- scnb : dihedral.type.scnb
        --- atomtype 1 : dihedral.atom1.type
        --- atomtype 2 : dihedral.atom2.type
        --- atomtype 3 : dihedral.atom3.type
        --- atomtype 4 : dihedral.atom4.type
    """
    # copy structure so the input structure isn't modified in-place
    structure = structure.copy(cls=Structure, split_dihedrals=True)
    if atom_style not in ["atomic", "charge", "molecular", "full"]:
        raise ValueError(
            'Atom style "{atom_style}" is invalid or is not currently supported'
        )

    # Check if structure is paramterized
    if unit_style == "lj":
        if any([atom.sigma for atom in structure.atoms]) is None:
            raise ValueError(
                "LJ units specified but one or more atoms has undefined LJ "
                "parameters."
            )

    xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in structure.atoms])

    forcefield = True
    if structure[0].type == "":
        forcefield = False

    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    charges = np.array([atom.charge for atom in structure.atoms])

    # Convert coordinates to LJ units
    if unit_style == "lj":
        # Get sigma, mass, and epsilon conversions by finding maximum of each
        sigma_conversion_factor = np.max([a.sigma for a in structure.atoms])
        epsilon_conversion_factor = np.max([a.epsilon for a in structure.atoms])
        mass_conversion_factor = np.max([a.mass for a in structure.atoms])

        xyz = xyz / sigma_conversion_factor
        charges = (charges * 1.6021e-19) / np.sqrt(
            4
            * np.pi
            * (sigma_conversion_factor * 1e-10)
            * (epsilon_conversion_factor * 4184)
            * epsilon_0
        )
        charges[np.isinf(charges)] = 0
    else:
        sigma_conversion_factor = 1
        epsilon_conversion_factor = 1
        mass_conversion_factor = 1
    # lammps does not require the box to be centered at any a specific origin
    # min and max dimensions are therefore needed to write the file in a
    # consistent way the parmed structure only stores the box length.  It is
    # not rigorous to assume bounds are 0 to L or -L/2 to L/2

    # NOTE: 0 to L is current default, mins and maxs should be passed by user

    if _check_minsmaxs(mins, maxs):
        box = Box.from_mins_maxs_angles(
            mins=mins, maxs=maxs, angles=structure.box[3:6]
        )
    else:
        # Internally use nm
        box = Box(
            lengths=np.array([0.1 * val for val in structure.box[0:3]]),
            angles=structure.box[3:6],
        )

        warn(
            "Explicit box bounds (i.e., mins and maxs) were not provided. Box "
            "bounds are assumed to be min = 0 and max = length in each "
            "direction. This may not produce a system with the expected "
            "spatial location and may cause non-periodic systems to fail. "
            "Bounds can be defined explicitly by passing the them to the "
            "write_lammpsdata function or by passing box info to the save "
            "function."
        )
    # Divide by conversion factor
    Lx = box.Lx * (1 / sigma_conversion_factor)
    Ly = box.Ly * (1 / sigma_conversion_factor)
    Lz = box.Lz * (1 / sigma_conversion_factor)
    box = Box(lengths=(Lx, Ly, Lz), angles=box.angles)

    # Lammps syntax depends on the functional form
    # Infer functional form based on the properties of the structure
    if detect_forcefield_style:
        # Check angles
        if len(structure.urey_bradleys) > 0:
            print("Urey bradley terms detected, will use angle_style charmm")
            use_urey_bradleys = True
        else:
            print(
                "No urey bradley terms detected, will use angle_style harmonic"
            )
            use_urey_bradleys = False

        # Check dihedrals
        if len(structure.rb_torsions) > 0:
            print("RB Torsions detected, will use dihedral_style opls")
            use_rb_torsions = True
        else:
            use_rb_torsions = False
        if len(structure.dihedrals) > 0:
            print("Charmm dihedrals detected, will use dihedral_style charmm")
            use_dihedrals = True
        else:
            use_dihedrals = False
    if use_rb_torsions and use_dihedrals:
        raise ValueError(
            "Multiple dihedral styles detected, check your "
            "Forcefield XML and structure"
        )

    bonds = [[b.atom1.idx + 1, b.atom2.idx + 1] for b in structure.bonds]
    angles = [
        [angle.atom1.idx + 1, angle.atom2.idx + 1, angle.atom3.idx + 1]
        for angle in structure.angles
    ]
    if use_rb_torsions:
        dihedrals = [
            [d.atom1.idx + 1, d.atom2.idx + 1, d.atom3.idx + 1, d.atom4.idx + 1]
            for d in structure.rb_torsions
        ]
    elif use_dihedrals:
        dihedrals = [
            [d.atom1.idx + 1, d.atom2.idx + 1, d.atom3.idx + 1, d.atom4.idx + 1]
            for d in structure.dihedrals
            if not d.improper
        ]
    else:
        dihedrals = []
    impropers = [
        [i.atom1.idx + 1, i.atom2.idx + 1, i.atom3.idx + 1, i.atom4.idx + 1]
        for i in structure.impropers
    ]
    imp_dihedrals = [
        [d.atom1.idx + 1, d.atom2.idx + 1, d.atom3.idx + 1, d.atom4.idx + 1]
        for d in structure.dihedrals
        if d.improper
    ]

    if impropers and imp_dihedrals:
        raise ValueError("Use of multiple improper styles is not supported")

    if bonds:
        if len(structure.bond_types) == 0:
            bond_types = np.ones(len(bonds), dtype=int)
        else:
            bond_types, unique_bond_types = _get_bond_types(
                structure,
                bonds,
                sigma_conversion_factor,
                epsilon_conversion_factor,
            )

    if angles:
        angle_types, unique_angle_types = _get_angle_types(
            structure,
            use_urey_bradleys,
            sigma_conversion_factor,
            epsilon_conversion_factor,
        )

    if imp_dihedrals:
        (
            imp_dihedral_types,
            unique_imp_dihedral_types,
        ) = _get_improper_dihedral_types(structure, epsilon_conversion_factor)
    if dihedrals:
        dihedral_types, unique_dihedral_types = _get_dihedral_types(
            structure,
            use_rb_torsions,
            use_dihedrals,
            epsilon_conversion_factor,
            zero_dihedral_weighting_factor,
        )

    if impropers:
        improper_types, unique_improper_types = _get_impropers(
            structure, epsilon_conversion_factor
        )

    with open(filename, "w") as data:
        data.write(f"{filename} - created by mBuild; units = {unit_style}\n\n")
        data.write("{:d} atoms\n".format(len(structure.atoms)))
        if atom_style in ["full", "molecular"]:
            data.write("{:d} bonds\n".format(len(bonds)))
            data.write("{:d} angles\n".format(len(angles)))
            data.write("{:d} dihedrals\n".format(len(dihedrals)))
            data.write(
                "{:d} impropers\n\n".format(len(impropers) + len(imp_dihedrals))
            )

        data.write("{:d} atom types\n".format(len(set(types))))
        if atom_style in ["full", "molecular"]:
            if bonds:
                data.write("{:d} bond types\n".format(len(set(bond_types))))
            if angles:
                data.write("{:d} angle types\n".format(len(set(angle_types))))
            if dihedrals:
                data.write(
                    "{:d} dihedral types\n".format(len(set(dihedral_types)))
                )
            if impropers:
                data.write(
                    "{:d} improper types\n".format(len(set(improper_types)))
                )
            elif imp_dihedrals:
                data.write(
                    "{:d} improper types\n".format(len(set(imp_dihedral_types)))
                )

        data.write("\n")
        # Box data
        # NOTE: Needs better logic handling maxs and mins of a bounding box
        # NOTE: JBG, "this should be a method/attribute of Compound?"
        if np.allclose(box.angles, 90.0, atol=1e-5) and (mins is None):
            for i, dim in enumerate(["x", "y", "z"]):
                data.write(
                    "{0:.6f} {1:.6f} {2}lo {2}hi\n".format(
                        0.0, 10.0 * box.lengths[i], dim
                    )
                )
        # NOTE:
        # currently non-orthogonal bounding box translates
        # Compound such that mins are new origin
        else:
            a = 10.0 * box.Lx
            b = 10.0 * box.Ly
            c = 10.0 * box.Lz
            alpha, beta, gamma = np.radians(box.angles)

            xy = box.xy
            xz = box.xz
            yz = box.yz

            # NOTE: using (0,0,0) as origin
            xlo, ylo, zlo = (0.0, 0.0, 0.0)
            xhi = xlo + a
            yhi = ylo + b
            zhi = zlo + c

            xlo_bound = xlo + np.min([0.0, xy, xz, xy + xz])
            xhi_bound = xhi + np.max([0.0, xy, xz, xy + xz])
            ylo_bound = ylo + np.min([0.0, yz])
            yhi_bound = yhi + np.max([0.0, yz])
            zlo_bound = zlo
            zhi_bound = zhi

            data.write("{0:.6f} {1:.6f} xlo xhi\n".format(xlo_bound, xhi_bound))
            data.write("{0:.6f} {1:.6f} ylo yhi\n".format(ylo_bound, yhi_bound))
            data.write("{0:.6f} {1:.6f} zlo zhi\n".format(zlo_bound, zhi_bound))
            data.write("{0:.6f} {1:.6f} {2:6f} xy xz yz\n".format(xy, xz, yz))

        # Mass data
        if not forcefield:
            masses = (
                np.array([atom.mass for atom in structure.atoms])
                / mass_conversion_factor
            )
        else:
            tmp_masses = list()
            for atom in structure.atoms:
                # handle case where atomtype does not contain a mass
                try:
                    tmp_masses.append(atom.atom_type.mass)
                except AttributeError:
                    warn(
                        f"No mass or defined atomtype for atom: {atom}. Using atom mass of {atom.mass / mass_conversion_factor}"
                    )
                    tmp_masses.append(atom.mass)
            masses = np.asarray(tmp_masses) / mass_conversion_factor

        mass_dict = dict(
            [
                (unique_types.index(atom_type) + 1, mass)
                for atom_type, mass in zip(types, masses)
            ]
        )
        data.write("\nMasses\n\n")
        for atom_type, mass in sorted(mass_dict.items()):
            data.write(
                "{:d}\t{:.6f}\t# {}\n".format(
                    atom_type, mass, unique_types[atom_type - 1]
                )
            )

        if forcefield:
            epsilons = (
                np.array([atom.epsilon for atom in structure.atoms])
                / epsilon_conversion_factor
            )
            sigmas = (
                np.array([atom.sigma for atom in structure.atoms])
                / sigma_conversion_factor
            )
            forcefields = [atom.type for atom in structure.atoms]
            epsilon_dict = dict(
                [
                    (unique_types.index(atom_type) + 1, epsilon)
                    for atom_type, epsilon in zip(types, epsilons)
                ]
            )
            sigma_dict = dict(
                [
                    (unique_types.index(atom_type) + 1, sigma)
                    for atom_type, sigma in zip(types, sigmas)
                ]
            )
            forcefield_dict = dict(
                [
                    (unique_types.index(atom_type) + 1, forcefield)
                    for atom_type, forcefield in zip(types, forcefields)
                ]
            )

            # Modified cross-interactions
            if structure.has_NBFIX():
                params = ParameterSet.from_structure(structure)
                # Sort keys (maybe they should be sorted in ParmEd)
                new_nbfix_types = OrderedDict()
                for key in params.nbfix_types.keys():
                    sorted_key = tuple(sorted(key))
                    if sorted_key in new_nbfix_types:
                        warn("Sorted key matches an existing key")
                        if new_nbfix_types[sorted_key]:
                            warn(
                                "nbfixes are not symmetric, overwriting old "
                                "nbfix"
                            )
                    new_nbfix_types[sorted_key] = params.nbfix_types[key]
                params.nbfix_types = new_nbfix_types
                warn(
                    "Explicitly writing cross interactions using mixing rule: "
                    "{}".format(structure.combining_rule)
                )
                coeffs = OrderedDict()
                for combo in it.combinations_with_replacement(unique_types, 2):
                    # Attempt to find pair coeffis in nbfixes
                    if combo in params.nbfix_types:
                        type1 = unique_types.index(combo[0]) + 1
                        type2 = unique_types.index(combo[1]) + 1
                        epsilon = params.nbfix_types[combo][
                            0
                        ]  # kcal OR lj units
                        rmin = params.nbfix_types[combo][
                            1
                        ]  # Angstrom OR lj units
                        sigma = rmin / 2 ** (1 / 6)
                        coeffs[(type1, type2)] = (
                            round(sigma, 8),
                            round(epsilon, 8),
                        )
                    else:
                        type1 = unique_types.index(combo[0]) + 1
                        type2 = unique_types.index(combo[1]) + 1
                        # Might not be necessary to be this explicit
                        if type1 == type2:
                            sigma = sigma_dict[type1]
                            epsilon = epsilon_dict[type1]
                        else:
                            if structure.combining_rule == "lorentz":
                                sigma = (
                                    sigma_dict[type1] + sigma_dict[type2]
                                ) * 0.5
                            elif structure.combining_rule == "geometric":
                                sigma = (
                                    sigma_dict[type1] * sigma_dict[type2]
                                ) ** 0.5
                            else:
                                raise ValueError(
                                    "Only lorentz and geometric combining "
                                    "rules are supported"
                                )
                            epsilon = (
                                epsilon_dict[type1] * epsilon_dict[type2]
                            ) ** 0.5
                        coeffs[(type1, type2)] = (
                            round(sigma, 8),
                            round(epsilon, 8),
                        )
                if nbfix_in_data_file:
                    if pair_coeff_label:
                        data.write(
                            "\nPairIJ Coeffs # {}\n".format(pair_coeff_label)
                        )
                    else:
                        data.write("\nPairIJ Coeffs # modified lj\n")

                    data.write(
                        "# type1 type2\tepsilon (kcal/mol)\tsigma (Angstrom)\n"
                    )

                    for (type1, type2), (sigma, epsilon) in coeffs.items():
                        data.write(
                            "{0} \t{1} \t{2} \t\t{3}\t\t# {4}\t{5}\n".format(
                                type1,
                                type2,
                                epsilon,
                                sigma,
                                forcefield_dict[type1],
                                forcefield_dict[type2],
                            )
                        )
                else:
                    if pair_coeff_label:
                        data.write(
                            "\nPair Coeffs # {} \n\n".format(pair_coeff_label)
                        )
                    else:
                        data.write("\nPair Coeffs # lj\n\n")

                    for idx, epsilon in sorted(epsilon_dict.items()):
                        data.write(
                            "{}\t{:.5f}\t{:.5f}\n".format(
                                idx, epsilon, sigma_dict[idx]
                            )
                        )
                    print("Copy these commands into your input script:\n")
                    print(
                        "# type1 type2\tepsilon (kcal/mol)\tsigma (Angstrom)\n"
                    )
                    for (type1, type2), (sigma, epsilon) in coeffs.items():
                        print(
                            "pair_coeff\t{0} \t{1} \t{2} \t\t{3} \t\t# {4} \t{5}".format(
                                type1,
                                type2,
                                epsilon,
                                sigma,
                                forcefield_dict[type1],
                                forcefield_dict[type2],
                            )
                        )

            # Pair coefficients
            else:
                if pair_coeff_label:
                    data.write(
                        "\nPair Coeffs # {} \n\n".format(pair_coeff_label)
                    )
                else:
                    data.write("\nPair Coeffs # lj\n\n")

                if unit_style == "real":
                    data.write("#\tepsilon (kcal/mol)\t\tsigma (Angstrom)\n")
                elif unit_style == "lj":
                    data.write("#\treduced_epsilon \t\treduced_sigma \n")
                for idx, epsilon in sorted(epsilon_dict.items()):
                    data.write(
                        "{}\t{:.5f}\t\t{:.5f}\t\t# {}\n".format(
                            idx, epsilon, sigma_dict[idx], forcefield_dict[idx]
                        )
                    )

            # Bond coefficients
            if bonds:
                data.write("\nBond Coeffs # harmonic\n")
                if unit_style == "real":
                    data.write("#\tk(kcal/mol/angstrom^2)\t\treq(angstrom)\n")
                elif unit_style == "lj":
                    data.write("#\treduced_k\t\treduced_req\n")
                sorted_bond_types = {
                    k: v
                    for k, v in sorted(
                        unique_bond_types.items(), key=lambda item: item[1]
                    )
                }
                for params, idx in sorted_bond_types.items():
                    data.write(
                        "{}\t{}\t\t{}\t\t# {}\t{}\n".format(
                            idx,
                            params[0],
                            params[1],
                            params[2][0],
                            params[2][1],
                        )
                    )

            # Angle coefficients
            if angles:
                sorted_angle_types = {
                    k: v
                    for k, v in sorted(
                        unique_angle_types.items(), key=lambda item: item[1]
                    )
                }
                if use_urey_bradleys:
                    data.write("\nAngle Coeffs # charmm\n")
                    data.write(
                        "#\tk(kcal/mol/rad^2)\t\ttheteq(deg)\tk(kcal/mol/angstrom^2)\treq(angstrom)\n"
                    )
                    for params, idx in sorted_angle_types.items():
                        data.write(
                            "{}\t{}\t{:.5f}\t{:.5f}\t{:.5f}\n".format(
                                idx, *params
                            )
                        )

                else:
                    data.write("\nAngle Coeffs # harmonic\n")
                    data.write("#\treduced_k\t\ttheteq(deg)\n")
                    for params, idx in sorted_angle_types.items():
                        data.write(
                            "{}\t{}\t\t{:.5f}\t# {}\t{}\t{}\n".format(
                                idx,
                                params[0],
                                params[1],
                                params[3][0],
                                params[2],
                                params[3][1],
                            )
                        )

            # Dihedral coefficients
            if dihedrals:
                sorted_dihedral_types = {
                    k: v
                    for k, v in sorted(
                        unique_dihedral_types.items(), key=lambda item: item[1]
                    )
                }
                if use_rb_torsions:
                    data.write("\nDihedral Coeffs # opls\n")
                    if unit_style == "real":
                        data.write(
                            "#\tf1(kcal/mol)\tf2(kcal/mol)\tf3(kcal/mol)\tf4(kcal/mol)\n"
                        )
                    elif unit_style == "lj":
                        data.write("#\tf1\tf2\tf3\tf4 (all lj reduced units)\n")
                    for params, idx in sorted_dihedral_types.items():
                        opls_coeffs = RB_to_OPLS(
                            params[0],
                            params[1],
                            params[2],
                            params[3],
                            params[4],
                            params[5],
                            error_if_outside_tolerance=False,
                        )
                        data.write(
                            "{}\t{:.5f}\t{:.5f}\t\t{:.5f}\t\t{:.5f}\t# {}\t{}\t{}\t{}\n".format(
                                idx,
                                opls_coeffs[1],
                                opls_coeffs[2],
                                opls_coeffs[3],
                                opls_coeffs[4],
                                params[8],
                                params[9],
                                params[10],
                                params[11],
                            )
                        )
                elif use_dihedrals:
                    data.write("\nDihedral Coeffs # charmm\n")
                    data.write("#k, n, phi, weight\n")
                    for params, idx in sorted_dihedral_types.items():
                        data.write(
                            "{}\t{:.5f}\t{:d}\t{:d}\t{:.5f}\t# {}\t{}\t{}\t{}\n".format(
                                idx,
                                params[0],
                                params[1],
                                params[2],
                                params[3],
                                params[6],
                                params[7],
                                params[8],
                                params[9],
                            )
                        )

            # Improper coefficients
            if impropers:
                sorted_improper_types = {
                    k: v
                    for k, v in sorted(
                        unique_improper_types.items(), key=lambda item: item[1]
                    )
                }
                data.write("\nImproper Coeffs # harmonic\n")
                data.write("#k, phi\n")
                for params, idx in sorted_improper_types.items():
                    data.write(
                        "{}\t{:.5f}\t{:.5f}\t# {}\t{}\t{}\t{}\n".format(
                            idx,
                            params[0],
                            params[1],
                            params[2],
                            params[3],
                            params[4],
                            params[5],
                        )
                    )
            elif imp_dihedrals:
                # Improper dihedral coefficients
                sorted_imp_dihedral_types = {
                    k: v
                    for k, v in sorted(
                        unique_imp_dihedral_types.items(),
                        key=lambda item: item[1],
                    )
                }
                data.write("\nImproper Coeffs # cvff\n")
                data.write("#K, d, n\n")
                for params, idx in sorted_imp_dihedral_types.items():
                    data.write(
                        "{}\t{:.5f}\t{:d}\t{:d}\t# {}\t{}\t{}\t{}\n".format(
                            idx,
                            params[0],
                            params[1],
                            params[2],
                            params[5],
                            params[6],
                            params[7],
                            params[8],
                        )
                    )

        # Atom data
        data.write("\nAtoms # {}\n\n".format(atom_style))
        if atom_style == "atomic":
            atom_line = "{index:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
        elif atom_style == "charge":
            if unit_style == "real":
                atom_line = "{index:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
            elif unit_style == "lj":
                atom_line = "{index:d}\t{type_index:d}\t{charge:.4ef}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
        elif atom_style == "molecular":
            atom_line = "{index:d}\t{zero:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
        elif atom_style == "full":
            if unit_style == "real":
                atom_line = "{index:d}\t{zero:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"
            elif unit_style == "lj":
                atom_line = "{index:d}\t{zero:d}\t{type_index:d}\t{charge:.4e}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n"

        for i, coords in enumerate(xyz):
            data.write(
                atom_line.format(
                    index=i + 1,
                    type_index=unique_types.index(types[i]) + 1,
                    zero=structure.atoms[i].residue.idx + moleculeID_offset,
                    charge=charges[i],
                    x=coords[0],
                    y=coords[1],
                    z=coords[2],
                )
            )

        if atom_style in ["full", "molecular"]:
            # Bond data
            if bonds:
                data.write("\nBonds\n\n")
                for i, bond in enumerate(bonds):
                    data.write(
                        "{:d}\t{:d}\t{:d}\t{:d}\n".format(
                            i + 1, bond_types[i], bond[0], bond[1]
                        )
                    )

            # Angle data
            if angles:
                data.write("\nAngles\n\n")
                for i, angle in enumerate(angles):
                    data.write(
                        "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                            i + 1, angle_types[i], angle[0], angle[1], angle[2]
                        )
                    )

            # Dihedral data
            if dihedrals:
                data.write("\nDihedrals\n\n")
                for i, dihedral in enumerate(dihedrals):
                    data.write(
                        "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                            i + 1,
                            dihedral_types[i],
                            dihedral[0],
                            dihedral[1],
                            dihedral[2],
                            dihedral[3],
                        )
                    )
            # Dihedral data
            if impropers:
                data.write("\nImpropers\n\n")
                for i, improper in enumerate(impropers):
                    data.write(
                        "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                            i + 1,
                            improper_types[i],
                            improper[2],
                            improper[1],
                            improper[0],
                            improper[3],
                        )
                    )
            elif imp_dihedrals:
                data.write("\nImpropers\n\n")
                for i, improper in enumerate(imp_dihedrals):
                    # The atoms are written central-atom third in LAMMPS data file.
                    # This is correct for AMBER impropers even though
                    # LAMMPS documentation implies central-atom-first.
                    data.write(
                        "{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                            i + 1,
                            imp_dihedral_types[i],
                            improper[0],
                            improper[1],
                            improper[2],
                            improper[3],
                        )
                    )


def _get_bond_types(
    structure, bonds, sigma_conversion_factor, epsilon_conversion_factor
):
    unique_bond_types = dict(
        enumerate(
            set(
                [
                    (
                        round(
                            bond.type.k
                            * (
                                sigma_conversion_factor ** 2
                                / epsilon_conversion_factor
                            ),
                            3,
                        ),
                        round(bond.type.req / sigma_conversion_factor, 3),
                        tuple(sorted((bond.atom1.type, bond.atom2.type))),
                    )
                    for bond in structure.bonds
                ]
            )
        )
    )
    unique_bond_types = OrderedDict(
        [(y, x + 1) for x, y in unique_bond_types.items()]
    )
    bond_types = [
        unique_bond_types[
            (
                round(
                    bond.type.k
                    * (
                        sigma_conversion_factor ** 2 / epsilon_conversion_factor
                    ),
                    3,
                ),
                round(bond.type.req / sigma_conversion_factor, 3),
                tuple(sorted((bond.atom1.type, bond.atom2.type))),
            )
        ]
        for bond in structure.bonds
    ]
    return bond_types, unique_bond_types


def _get_angle_types(
    structure,
    use_urey_bradleys,
    sigma_conversion_factor,
    epsilon_conversion_factor,
):
    if use_urey_bradleys:
        charmm_angle_types = []
        for angle in structure.angles:
            ub_k = 0
            ub_req = 0
            for ub in structure.urey_bradleys:
                if (angle.atom1, angle.atom3) == (ub.atom1, ub.atom2):
                    ub_k = ub.type.k
                    ub_req = ub.type.req
            charmm_angle_types.append(
                (
                    round(
                        angle.type.k
                        * (
                            sigma_conversion_factor ** 2
                            / epsilon_conversion_factor
                        ),
                        3,
                    ),
                    round(angle.type.theteq, 3),
                    round(ub_k / epsilon_conversion_factor, 3),
                    round(ub_req, 3),
                    tuple(sorted((angle.atom1.type, angle.atom3.type))),
                )
            )

        unique_angle_types = dict(enumerate(set(charmm_angle_types)))
        unique_angle_types = OrderedDict(
            [(y, x + 1) for x, y in unique_angle_types.items()]
        )
        angle_types = [
            unique_angle_types[ub_info] for ub_info in charmm_angle_types
        ]

    else:
        unique_angle_types = dict(
            enumerate(
                set(
                    [
                        (
                            round(
                                angle.type.k
                                * (
                                    sigma_conversion_factor ** 2
                                    / epsilon_conversion_factor
                                ),
                                3,
                            ),
                            round(angle.type.theteq, 3),
                            angle.atom2.type,
                            tuple(sorted((angle.atom1.type, angle.atom3.type))),
                        )
                        for angle in structure.angles
                    ]
                )
            )
        )
        unique_angle_types = OrderedDict(
            [(y, x + 1) for x, y in unique_angle_types.items()]
        )
        angle_types = [
            unique_angle_types[
                (
                    round(
                        angle.type.k
                        * (
                            sigma_conversion_factor ** 2
                            / epsilon_conversion_factor
                        ),
                        3,
                    ),
                    round(angle.type.theteq, 3),
                    angle.atom2.type,
                    tuple(sorted((angle.atom1.type, angle.atom3.type))),
                )
            ]
            for angle in structure.angles
        ]

    return angle_types, unique_angle_types


def _get_dihedral_types(
    structure,
    use_rb_torsions,
    use_dihedrals,
    epsilon_conversion_factor,
    zero_dihedral_weighting_factor,
):
    lj_unit = 1 / epsilon_conversion_factor
    if use_rb_torsions:
        unique_dihedral_types = dict(
            enumerate(
                set(
                    [
                        (
                            round(dihedral.type.c0 * lj_unit, 3),
                            round(dihedral.type.c1 * lj_unit, 3),
                            round(dihedral.type.c2 * lj_unit, 3),
                            round(dihedral.type.c3 * lj_unit, 3),
                            round(dihedral.type.c4 * lj_unit, 3),
                            round(dihedral.type.c5 * lj_unit, 3),
                            round(dihedral.type.scee, 1),
                            round(dihedral.type.scnb, 1),
                            dihedral.atom1.type,
                            dihedral.atom2.type,
                            dihedral.atom3.type,
                            dihedral.atom4.type,
                        )
                        for dihedral in structure.rb_torsions
                    ]
                )
            )
        )
        unique_dihedral_types = OrderedDict(
            [(y, x + 1) for x, y in unique_dihedral_types.items()]
        )
        dihedral_types = [
            unique_dihedral_types[
                (
                    round(dihedral.type.c0 * lj_unit, 3),
                    round(dihedral.type.c1 * lj_unit, 3),
                    round(dihedral.type.c2 * lj_unit, 3),
                    round(dihedral.type.c3 * lj_unit, 3),
                    round(dihedral.type.c4 * lj_unit, 3),
                    round(dihedral.type.c5 * lj_unit, 3),
                    round(dihedral.type.scee, 1),
                    round(dihedral.type.scnb, 1),
                    dihedral.atom1.type,
                    dihedral.atom2.type,
                    dihedral.atom3.type,
                    dihedral.atom4.type,
                )
            ]
            for dihedral in structure.rb_torsions
        ]
    elif use_dihedrals:
        charmm_dihedrals = []
        structure.join_dihedrals()
        for dihedral in structure.dihedrals:
            if not dihedral.improper:
                if zero_dihedral_weighting_factor:
                    weight = 0.0
                else:
                    weight = 1.0 / len(dihedral.type)
                for dih_type in dihedral.type:
                    charmm_dihedrals.append(
                        (
                            round(dih_type.phi_k * lj_unit, 3),
                            int(round(dih_type.per, 0)),
                            int(round(dih_type.phase, 0)),
                            round(weight, 4),
                            round(dih_type.scee, 1),
                            round(dih_type.scnb, 1),
                            dihedral.atom1.type,
                            dihedral.atom2.type,
                            dihedral.atom3.type,
                            dihedral.atom4.type,
                        )
                    )

        unique_dihedral_types = dict(enumerate(set(charmm_dihedrals)))
        unique_dihedral_types = OrderedDict(
            [(y, x + 1) for x, y in unique_dihedral_types.items()]
        )
        dihedral_types = [
            unique_dihedral_types[dihedral_info]
            for dihedral_info in charmm_dihedrals
        ]

    return dihedral_types, unique_dihedral_types


def _get_improper_dihedral_types(structure, epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    improper_dihedrals = []
    for dihedral in structure.dihedrals:
        if dihedral.improper:
            dih_type = dihedral.type
            phase = abs(int(round(dih_type.phase, 0)))
            if not (phase == 0 or phase == 180):
                raise ValueError("Improper dihedral phase must be 0 or 180")
            if phase:
                d = -1
            else:
                d = 1
            improper_dihedrals.append(
                (
                    round(dih_type.phi_k * lj_unit, 3),
                    d,
                    int(round(dih_type.per, 0)),
                    round(dih_type.scee, 1),
                    round(dih_type.scnb, 1),
                    dihedral.atom1.type,
                    dihedral.atom2.type,
                    dihedral.atom3.type,
                    dihedral.atom4.type,
                )
            )
    unique_imp_dihedral_types = dict(enumerate(set(improper_dihedrals)))
    unique_imp_dihedral_types = OrderedDict(
        [(y, x + 1) for x, y in unique_imp_dihedral_types.items()]
    )
    imp_dihedral_types = [
        unique_imp_dihedral_types[dihedral_info]
        for dihedral_info in improper_dihedrals
    ]

    return imp_dihedral_types, unique_imp_dihedral_types


def _get_impropers(structure, epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    unique_improper_types = dict(
        enumerate(
            set(
                [
                    (
                        round(improper.type.psi_k * lj_unit, 3),
                        round(improper.type.psi_eq, 3),
                        improper.atom3.type,
                        improper.atom2.type,
                        improper.atom1.type,
                        improper.atom4.type,
                    )
                    for improper in structure.impropers
                ]
            )
        )
    )
    unique_improper_types = OrderedDict(
        [(y, x + 1) for x, y in unique_improper_types.items()]
    )
    improper_types = [
        unique_improper_types[
            (
                round(improper.type.psi_k * lj_unit, 3),
                round(improper.type.psi_eq, 3),
                improper.atom3.type,
                improper.atom2.type,
                improper.atom1.type,
                improper.atom4.type,
            )
        ]
        for improper in structure.impropers
    ]

    return improper_types, unique_improper_types


def _get_box_information(
    structure,
):
    pass
