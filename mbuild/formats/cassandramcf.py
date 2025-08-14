"""Cassandra Molecular Connectivity format.

https://cassandra-mc.readthedocs.io/en/latest/guides/input_files.html#molecular-connectivity-file
"""

from __future__ import division

import logging
from math import sqrt

import networkx as nx
import parmed as pmd

__all__ = ["write_mcf"]
logger = logging.getLogger(__name__)


def write_mcf(structure, filename, angle_style, dihedral_style, lj14=None, coul14=None):
    """Output a Cassandra molecular connectivity file (MCF).

    Outputs a Cassandra MCF from a Parmed structure object.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file
    angle_style : str
        Type of angles. 'fixed' and 'harmonic' are valid choices
    dihedral_style : str
        Type of dihedrals. 'harmonic', 'OPLS', 'CHARMM', and 'none' are
        valid choices
    lj14 : float
        Scaling factor for LJ interactions on 1-4 pairs
    coul14 : float
        Scaling factor for Coulombic interactions on 1-4 pairs

    Notes
    -----
    See https://cassandra.nd.edu/index.php/documentation for
    a complete description of the MCF format.
    """
    if not isinstance(structure, pmd.Structure):
        raise ValueError("MCF writer requires parmed structure.")
    if not all([a.type for a in structure.atoms]):
        raise ValueError("MCF writing not supported without parameterized forcefield.")

    # Conversion factors
    IG_CONSTANT_KCAL = 0.00198720425864083  # kcal/mol*K
    KCAL_TO_KJ = 4.184

    # Check some things before we start writing the MCF
    # Only will write MCF for Cassandra-supported options
    if angle_style.casefold() != "fixed" and angle_style.casefold() != "harmonic":
        raise ValueError(
            "Invalid selection for angle_style. Please choose 'fixed' or 'harmonic'"
        )

    if len(structure.urey_bradleys) > 0:
        raise ValueError(
            "Urey bradley terms detected. Cassandra only currently supports "
            "fixed or harmonic angles."
        )

    if (
        dihedral_style.casefold() != "opls"
        and dihedral_style.casefold() != "charmm"
        and dihedral_style.casefold() != "none"
    ):
        raise ValueError(
            "Invalid selection for dihedral_style. Please choose 'OPLS', "
            "'CHARMM', or 'none'"
        )

    if dihedral_style.casefold() != "none":
        if len(structure.rb_torsions) > 0 and dihedral_style.casefold() != "opls":
            raise ValueError(
                "Dihedral style declared as {} but RB torsions found.".format(
                    dihedral_style
                )
            )

        if len(structure.dihedrals) > 0 and dihedral_style.casefold() != "charmm":
            raise ValueError(
                "Dihedral style declared as {} but charmm-style dihedrals "
                "found.".format(dihedral_style)
            )

        if len(structure.rb_torsions) > 0 and len(structure.dihedrals) > 0:
            raise ValueError(
                "Multiple dihedral styles detected, check your Forcefield XML "
                "and structure"
            )

    # Identify atoms in rings and Cassandra 'fragments'
    in_ring, frag_list, frag_conn = _id_rings_fragments(structure)

    # Infer 1-4 scaling if not specified
    if coul14 is None:
        if len(structure.adjusts) > 0:
            coul14 = structure.adjusts[0].type.chgscale
        else:
            coul14 = 0.0
            if len(structure.dihedrals) > 0 or len(structure.rb_torsions) > 0:
                logger.info(
                    "Unable to infer coulombic 1-4 scaling factor. Setting to "
                    "{:.1f}".format(coul14)
                )
    if lj14 is None:
        if len(structure.adjusts) > 0:
            if (
                structure.combining_rule == "geometric"
                or structure.combining_rule == "lorentz"
            ):
                combined_eps_list = [
                    sqrt(adj.atom1.epsilon * adj.atom2.epsilon)
                    for adj in structure.adjusts
                ]
                if all([c_eps == 0 for c_eps in combined_eps_list]):
                    lj14 = 0.0
                    logger.info(
                        "Unable to infer LJ 1-4 scaling factor. Setting to "
                        "{:.1f}".format(lj14)
                    )
                else:
                    scaled_eps_list = [adj.type.epsilon for adj in structure.adjusts]
                    for i_adj, combined_eps in enumerate(combined_eps_list):
                        if combined_eps != 0:
                            lj14 = scaled_eps_list[i_adj] / combined_eps
                            break
            else:
                lj14 = 0.0
                logger.info(
                    "Unable to infer LJ 1-4 scaling factor. Setting to {:.1f}".format(
                        lj14
                    )
                )
        else:
            lj14 = 0.0
            if len(structure.dihedrals) > 0 or len(structure.rb_torsions) > 0:
                logger.info(
                    "Unable to infer LJ 1-4 scaling factor. Setting to {:.1f}".format(
                        lj14
                    )
                )

    if coul14 < 0.0 or coul14 > 1.0:
        raise ValueError("Unreasonable value {} for coul14 scaling.".format(coul14))
    if lj14 < 0.0 or lj14 > 1.0:
        raise ValueError("Unreasonable value {} for lj14 scaling.".format(lj14))

    # Now we write the MCF file
    with open(filename, "w") as mcf_file:
        header = (
            "!***************************************"
            "****************************************\n"
            "!Molecular connectivity file\n"
            "!***************************************"
            "****************************************\n"
            "!" + filename + " - created by mBuild\n\n"
        )

        mcf_file.write(header)
        _write_atom_information(mcf_file, structure, in_ring, IG_CONSTANT_KCAL)
        _write_bond_information(mcf_file, structure)
        _write_angle_information(mcf_file, structure, angle_style, IG_CONSTANT_KCAL)
        _write_dihedral_information(mcf_file, structure, dihedral_style, KCAL_TO_KJ)
        _write_improper_information(mcf_file, structure, KCAL_TO_KJ)
        _write_fragment_information(mcf_file, structure, frag_list, frag_conn)
        _write_intrascaling_information(mcf_file, lj14, coul14)

        # That's all, folks!
        mcf_file.write("\n\nEND\n")


def _id_rings_fragments(structure):
    """Identify the rings and fragments of the molecule.

    Parameters
    ----------
    structure : parmed.Structure
        Parmed structure object

    Returns
    -------
    in_ring : list
        True for each atom in a ring
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments
    """
    # Identify atoms in rings
    bond_graph = nx.Graph()
    bond_graph.add_edges_from(
        [[bond.atom1.idx, bond.atom2.idx] for bond in structure.bonds]
    )

    if len(structure.bonds) == 0:
        logger.info("No bonds found. Cassandra will interpet this as a rigid species")
        in_ring = [False] * len(structure.atoms)
        frag_list = []
        frag_conn = []
        return in_ring, frag_list, frag_conn

    # Check if entire molecule is connected. Warn if not.
    if nx.is_connected(bond_graph) is False:
        raise ValueError(
            "Not all components of the molecule are connected. MCF files are "
            "for a single molecule and thus everything should be connected "
            "through bonds."
        )

    all_rings = nx.cycle_basis(bond_graph)
    in_ring = [False] * bond_graph.number_of_nodes()
    adj_to_ring = [False] * bond_graph.number_of_nodes()
    for ring in all_rings:
        for idx in ring:
            in_ring[idx] = True

    # Identify fragments
    # See Shah and Maginn, JCP, 135, 134121, 2011, doi:10.1063/1.3644939
    frag_list = []
    frag_conn = []

    # First create a neighbor list for each atom
    neigh_dict = {
        i: list(bond_graph.neighbors(i)) for i in range(bond_graph.number_of_nodes())
    }

    # Handle fused/adjoining rings
    rings_changed = True
    while rings_changed:
        rings_changed = False
        for ring1 in all_rings:
            if rings_changed:
                break
            for ring2 in all_rings:
                if ring1 == ring2:
                    continue
                if len(set(ring1) & set(ring2)) > 0:
                    all_rings.remove(ring1)
                    all_rings.remove(ring2)
                    all_rings.append(list(set(ring1 + ring2)))
                    rings_changed = True
                    break

    # ID fragments which contain a ring
    for ring in all_rings:
        adjacentatoms = []
        for idx in ring:
            if len(neigh_dict[idx]) > 2:
                adjacentatoms.append(list(set(neigh_dict[idx]) - set(ring)))
        tmp = filter(None, adjacentatoms)
        adjacentatoms = [x for sublist in tmp for x in sublist]
        frag_list.append(ring + adjacentatoms)
        for idx in adjacentatoms:
            adj_to_ring[idx] = True
    # Now ID the other fragments
    for idx in neigh_dict:
        if len(neigh_dict[idx]) > 1:
            if in_ring[idx] is True:
                continue
            else:
                frag_list.append([idx] + neigh_dict[idx])
    # Now find connectivity (shared bonds)
    for i in range(len(frag_list)):
        frag1 = frag_list[i]
        for j in range(i + 1, len(frag_list)):
            frag2 = frag_list[j]
            shared_atoms = list(set(frag1) & set(frag2))
            if len(shared_atoms) == 2:
                frag_conn.append([i, j])
            elif len(shared_atoms) > 2:
                logger.warning(
                    "Fragments share more than two atoms... something may be "
                    "going awry unless there are fused rings in your system. "
                    "See below for details."
                )
                print("Fragment 1 atoms:")
                print(frag1)
                print("Fragment 2 atoms:")
                print(frag2)

    return in_ring, frag_list, frag_conn


def _write_atom_information(mcf_file, structure, in_ring, IG_CONSTANT_KCAL):
    """Write the atoms in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    in_ring : list
        Boolean for each atom idx True if atom belongs to a ring
    IG_CONSTANT_KCAL : float
        Ideal gas constant in kcal/mol K
    """
    elements = [atom.element_name for atom in structure.atoms]
    types = [atom.type for atom in structure.atoms]
    masses = [atom.mass for atom in structure.atoms]
    charges = [atom.charge for atom in structure.atoms]
    # Convert energy to units of K
    epsilons = [atom.epsilon / IG_CONSTANT_KCAL for atom in structure.atoms]
    sigmas = [atom.sigma for atom in structure.atoms]

    # Check constraints on atom type length and element name length
    max_element_length = 6
    max_atomtype_length = 20
    n_unique_elements = len(set(elements))
    for element in elements:
        if len(element) > max_element_length:
            logger.info(
                "Element name {} will be shortened to {} characters. Please "
                "confirm your final MCF.".format(element, max_element_length)
            )

    elements = [element[:max_element_length] for element in elements]
    if len(set(elements)) < n_unique_elements:
        logger.info(
            "The number of unique elements has been reduced due to shortening "
            "the element name to {} characters.".format(max_element_length)
        )

    n_unique_types = len(set(types))
    for itype in types:
        if len(itype) > max_atomtype_length:
            logger.info(
                "Type name {} will be shortened to {} characters as {}. Please "
                "confirm your final MCF.".format(
                    itype, max_atomtype_length, itype[-max_atomtype_length:]
                )
            )
        types = [itype[-max_atomtype_length:] for itype in types]
    if len(set(types)) < n_unique_types:
        logger.info(
            "The number of unique atomtypes has been reduced due to shortening "
            "the atomtype name to {} characters.".format(max_atomtype_length)
        )

    vdw_type = "LJ"
    header = (
        "!Atom Format\n"
        "!index type element mass charge vdw_type parameters\n"
        '!vdw_type="LJ", parms=epsilon sigma\n'
        '!vdw_type="Mie", parms=epsilon sigma '
        "repulsion_exponent dispersion_exponent\n"
        "\n# Atom_Info\n"
    )

    mcf_file.write(header)
    mcf_file.write("{:d}\n".format(len(structure.atoms)))
    for i in range(len(structure.atoms)):
        mcf_file.write(
            "{:<4d}  {:<6s}  {:<2s}  {:7.3f}  {:12.8f}  {:3s}  {:8.5f}  {:8.5f}".format(
                i + 1,
                types[i],
                elements[i],
                masses[i],
                charges[i],
                vdw_type,
                epsilons[i],
                sigmas[i],
            )
        )
        if in_ring[i] is True:
            mcf_file.write("  ring")
        mcf_file.write("\n")


def _write_bond_information(mcf_file, structure):
    """Write the bonds in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    """
    bond_parms = [f"{bond.type.req:8.3f}" for bond in structure.bonds]

    mcf_file.write("\n!Bond Format\n")
    mcf_file.write("!index i j type parameters\n" + '!type="fixed", parms=bondLength\n')
    mcf_file.write("\n# Bond_Info\n")
    mcf_file.write("{:d}\n".format(len(structure.bonds)))
    for i, bond in enumerate(structure.bonds):
        mcf_file.write(
            "{:<4d}  {:<4d}  {:<4d}  {:s}  {:s}\n".format(
                i + 1,
                bond.atom1.idx + 1,
                bond.atom2.idx + 1,
                "fixed",
                bond_parms[i],
            )
        )


def _write_angle_information(mcf_file, structure, angle_style, IG_CONSTANT_KCAL):
    """Write the angles in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    angle_style : string
        Angle style for Cassandra to use
    IG_CONSTANT_KCAL : float
        Ideal gas constant in kcal/mol K
    """
    if angle_style.casefold() == "fixed":
        angle_parms = [f"{a.type.theteq:8.2f}" for a in structure.angles]
    elif angle_style.casefold() == "harmonic":
        # Convert energies to units of K
        angle_parms = [
            f"{a.type.k / IG_CONSTANT_KCAL:8.1f}  {a.type.theteq:8.2f}"
            for a in structure.angles
        ]
    else:
        raise ValueError(
            "Only 'fixed' and 'harmonic' angle styles are supported by Cassandra"
        )

    header = (
        "\n!Angle Format\n"
        "!index i j k type parameters\n"
        '!type="fixed", parms=equilibrium_angle\n'
        '!type="harmonic", parms=force_constant equilibrium_angle\n'
        "\n# Angle_Info\n"
    )

    mcf_file.write(header)
    mcf_file.write("{:d}\n".format(len(structure.angles)))

    for i, angle in enumerate(structure.angles):
        mcf_file.write(
            "{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:s}  {:s}\n".format(
                i + 1,
                angle.atom1.idx + 1,
                angle.atom2.idx + 1,
                angle.atom3.idx + 1,
                angle_style.lower(),
                angle_parms[i],
            )
        )


def _write_dihedral_information(mcf_file, structure, dihedral_style, KCAL_TO_KJ):
    """Write the dihedrals in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    dihedral_style : string
        Dihedral style for Cassandra to use
    KCAL_TO_KJ : float
        Conversion factor from kcal to kJ
    """
    # Dihedral info
    header = (
        "\n!Dihedral Format\n"
        "!index i j k l type parameters\n"
        '!type="none"\n'
        '!type="CHARMM", parms=a0 a1 delta\n'
        '!type="OPLS", parms=c0 c1 c2 c3\n'
        '!type="harmonic", parms=force_constant equilibrium_dihedral\n'
        "\n# Dihedral_Info\n"
    )

    mcf_file.write(header)

    if len(structure.dihedrals) > 0 or len(structure.rb_torsions) > 0:
        if dihedral_style.casefold() == "opls":
            dihedral_style = dihedral_style.upper()
            dihedrals = structure.rb_torsions
            # Two things happen here:
            #  (1) convert from RB form to Cassandra OPLS
            #  (2) convert units from kcal/mol to kJ/mol
            dihedral_parms = []
            for dihedral in dihedrals:
                a0 = (
                    dihedral.type.c0
                    + dihedral.type.c1
                    + dihedral.type.c2
                    + dihedral.type.c3
                )
                a1 = -dihedral.type.c1 - (3.0 / 4.0) * dihedral.type.c3
                a2 = (-1.0 / 2.0) * dihedral.type.c2
                a3 = (-1.0 / 4.0) * dihedral.type.c3
                if not dihedral.type.c4 == 0.0 and dihedral.type.c4 == 0.0:
                    raise ValueError(
                        "Can only convert Ryckaert-Bellemans dihedrals to OPLS "
                        "if c4==0 and c5==0"
                    )
                dihedral_parms.append(
                    str("{:8.3f} ".format(a0 * KCAL_TO_KJ))
                    + str("{:8.3f} ".format(a1 * KCAL_TO_KJ))
                    + str("{:8.3f} ".format(a2 * KCAL_TO_KJ))
                    + str("{:8.3f} ".format(a3 * KCAL_TO_KJ))
                )

        elif dihedral_style.casefold() == "charmm":
            dihedral_style = dihedral_style.upper()
            dihedrals = structure.dihedrals
            dihedral_parms = [
                str("{:8.3f} ".format(dihedral.type.phi_k * KCAL_TO_KJ))
                + str("{:8.3f} ".format(dihedral.type.per))
                + str("{:8.3f}".format(dihedral.type.phase))
                for dihedral in dihedrals
            ]

        elif dihedral_style.casefold() == "none":
            logger.info("Dihedral style 'none' selected. Ignoring dihedral parameters")
            dihedral_style = dihedral_style.lower()
            if structure.dihedrals:
                dihedrals = structure.dihedrals
                dihedral_parms = ["" for dihedral in dihedrals]
            elif structure.rb_torsions:
                dihedrals = structure.rb_torsions
                dihedral_parms = ["" for dihedral in dihedrals]
        else:
            raise ValueError(
                "Only 'OPLS', 'CHARMM', and 'none' dihedral styles are supported."
            )

        mcf_file.write("{:d}\n".format(len(dihedrals)))
        for i, dihedral in enumerate(dihedrals):
            # 2021 Jun 24 Ryan S. DeFever
            # Removed improper sorting for reproducibility;
            # The atom order provided in the parmed.Structure
            # is written to the MCF file.
            mcf_file.write(
                "{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:<4d}  {:s}  {:s}\n".format(
                    i + 1,
                    dihedral.atom1.idx + 1,
                    dihedral.atom2.idx + 1,
                    dihedral.atom3.idx + 1,
                    dihedral.atom4.idx + 1,
                    dihedral_style,
                    dihedral_parms[i],
                )
            )
    else:
        mcf_file.write("0\n")


def _write_improper_information(mcf_file, structure, KCAL_TO_KJ):
    """Write the impropers in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    KCAL_TO_KJ : float
        Conversion factor from kcal to kJ
    """
    header = (
        "\n!Improper Format\n"
        "!index i j k l type parameters\n"
        '!type="harmonic", parms=force_constant equilibrium_improper\n'
        "\n# Improper_Info\n"
    )

    mcf_file.write(header)
    mcf_file.write("{:d}\n".format(len(structure.impropers)))

    improper_type = "harmonic"
    for i, improper in enumerate(structure.impropers):
        mcf_file.write(
            "{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:<4d}  {:s}  {:8.3f}  {:8.3f}\n".format(
                i + 1,
                improper.atom1.idx + 1,
                improper.atom2.idx + 1,
                improper.atom3.idx + 1,
                improper.atom4.idx + 1,
                improper_type,
                improper.type.psi_k * KCAL_TO_KJ,
                improper.type.psi_eq,
            )
        )


def _write_fragment_information(mcf_file, structure, frag_list, frag_conn):
    """Write the fragments in the molecule.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments
    """
    header = (
        "\n!Fragment Format\n"
        "!index number_of_atoms_in_fragment branch_point other_atoms\n"
        "\n# Fragment_Info\n"
    )

    mcf_file.write(header)

    # Special cases first
    if len(frag_list) == 0:
        if len(structure.atoms) == 1:
            mcf_file.write("1\n")
            mcf_file.write("1 1 1\n")
        elif len(structure.atoms) == 2:
            mcf_file.write("1\n")
            mcf_file.write("1 2 1 2\n")
        else:
            logger.info("More than two atoms present but no fragments identified.")
            mcf_file.write("0\n")
    else:
        mcf_file.write("{:d}\n".format(len(frag_list)))
        for i, frag in enumerate(frag_list):
            mcf_file.write("{:d}    {:d}".format(i + 1, len(frag)))
            for idx in frag:
                mcf_file.write("    {:d}".format(idx + 1))
            mcf_file.write("\n")

    mcf_file.write("\n\n# Fragment_Connectivity\n")
    mcf_file.write("{:d}\n".format(len(frag_conn)))
    for i, conn in enumerate(frag_conn):
        mcf_file.write("{:d}    {:d}    {:d}\n".format(i + 1, conn[0] + 1, conn[1] + 1))


def _write_intrascaling_information(mcf_file, lj14, coul14):
    """Write the intramolecular scaling in the molecule.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    lj14 : float
        The 1-4 scaling parameter for LJ interactions
    coul14 : float
        The 1-4 scaling parameter for Coulombic interactions
    """
    header = (
        "\n!Intra Scaling\n"
        "!vdw_scaling    1-2 1-3 1-4 1-N\n"
        "!charge_scaling 1-2 1-3 1-4 1-N\n"
        "\n# Intra_Scaling\n"
    )

    mcf_file.write(header)
    mcf_file.write("0. 0. {:.4f} 1.\n".format(lj14))
    mcf_file.write("0. 0. {:.4f} 1.\n".format(coul14))
