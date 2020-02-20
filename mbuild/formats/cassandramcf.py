from __future__ import division

import warnings

from math import sqrt

import networkx as nx
import parmed as pmd

__all__ = ['write_mcf']


def write_mcf(structure, filename, angle_style,
                      dihedral_style, lj14=None, coul14=None):
    """Output a Cassandra molecular connectivity file (MCF).

    Outputs a Cassandra MCF from a Parmed structure object.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file
    angle_style : str
        Type of angles. 'fixed' and 'harmonic'
        are valid choices
    dihedral_style : str
        Type of dihedrals. 'harmonic', 'OPLS', 'CHARMM',
        and 'none' are valid choices
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
    if not all ([a.type for a in structure.atoms]):
        raise ValueError("MCF writing not supported without "
            "parameterized forcefield.")

    # Conversion factors
    IG_CONSTANT_KCAL = 0.00198720425864083 # kcal/mol*K
    KCAL_TO_KJ = 4.184

    # Check some things before we start writing the MCF
    # Only will write MCF for Cassandra-supported options
    if (angle_style.casefold() != 'fixed' and
            angle_style.casefold() != 'harmonic'):
        raise ValueError("Invalid selection for angle_style. "
                "Please choose 'fixed' or 'harmonic'")

    if len(structure.urey_bradleys) > 0 :
        raise ValueError("Urey bradley terms detected. Cassandra only "
                        "currently supports fixed or harmonic angles.")

    if (dihedral_style.casefold() != 'opls' and
            dihedral_style.casefold() != 'charmm' and
            dihedral_style.casefold() != 'none'):
        raise ValueError("Invalid selection for dihedral_style. "
                "Please choose 'OPLS', 'CHARMM', or 'none'")

    if dihedral_style.casefold() != 'none':
        if (len(structure.rb_torsions) > 0 and
                dihedral_style.casefold() != 'opls'):
            raise ValueError("Dihedral style declared as {} but "
                "RB torsions found.".format(dihedral_style))

        if (len(structure.dihedrals) > 0 and
                dihedral_style.casefold() != 'charmm'):
            raise ValueError("Dihedral style declared as {} but "
                "charmm-style dihedrals found.".format(dihedral_style))

        if (len(structure.rb_torsions) > 0 and
                len(structure.dihedrals) > 0):
            raise ValueError("Multiple dihedral styles detected, check your "
                             "Forcefield XML and structure")


    # Identify atoms in rings and Cassandra 'fragments'
    in_ring,frag_list,frag_conn = _id_rings_fragments(structure)

    # Infer 1-4 scaling if not specified
    if coul14 is None:
        if len(structure.adjusts) > 0:
            coul14 = structure.adjusts[0].type.chgscale
        else:
            coul14 = 1.0
            if (len(structure.dihedrals) > 0 or
                    len(structure.rb_torsions) > 0):
                warnings.warn('Unable to infer coulombic 1-4 '
                        'scaling factor. Setting to 1.0')
    if lj14 is None:
        if len(structure.adjusts) > 0:
            type1_eps = structure.adjusts[0].atom1.epsilon
            type2_eps = structure.adjusts[0].atom2.epsilon
            scaled_eps = structure.adjusts[0].type.epsilon
            if (structure.combining_rule == 'geometric' or
                    structure.combining_rule == 'lorentz'):
                combined_eps = sqrt(type1_eps*type2_eps)
                lj14 = scaled_eps/combined_eps
            else:
                lj14 = 1.0
                warnings.warn('Unable to infer LJ 1-4 scaling'
                    'factor. Setting to 1.0')
        else:
            lj14 = 1.0
            if (len(structure.dihedrals) > 0 or
                    len(structure.rb_torsions) > 0):
                warnings.warn('Unable to infer LJ 1-4 scaling'
                    'factor. Setting to 1.0')

    if coul14 < 0.0 or coul14 > 1.0:
        raise ValueError("Unreasonable value {} for "
                "coul14 scaling.".format(coul14))
    if lj14 < 0.0 or lj14 > 1.0:
        raise ValueError("Unreasonable value {} for "
                "lj14 scaling.".format(lj14))

    # Now we write the MCF file
    with open(filename, 'w') as mcf_file:

        header = ( '!***************************************'
                   '****************************************\n'
                   '!Molecular connectivity file\n'
                   '!***************************************'
                   '****************************************\n'
                   '!'+filename+' - created by mBuild\n\n'
                 )

        mcf_file.write(header)
        _write_atom_information(mcf_file, structure, in_ring,
                IG_CONSTANT_KCAL)
        _write_bond_information(mcf_file, structure)
        _write_angle_information(mcf_file, structure, angle_style,
                IG_CONSTANT_KCAL)
        _write_dihedral_information(mcf_file, structure, dihedral_style,
                KCAL_TO_KJ)
        _write_improper_information(mcf_file, structure, KCAL_TO_KJ)
        _write_fragment_information(mcf_file, structure, frag_list, frag_conn)
        _write_intrascaling_information(mcf_file, lj14, coul14)

        # That's all, folks!
        mcf_file.write('\n\nEND\n')


def _id_rings_fragments(structure):
    """Identifies the rings and fragments of the molecule

    Parameters
    ----------
    structure : parmed.Structure
        Parmed structure object

    Returns
    ---------
    in_ring : list
        True for each atom in a ring
    frag_list : list
        Atom ids belonging to each fragment
    frag_conn : list
        Fragment ids of connected fragments

    """

    # Identify atoms in rings
    bond_graph = nx.Graph()
    bond_graph.add_edges_from([ [bond.atom1.idx, bond.atom2.idx]
                                 for bond in structure.bonds ])

    if len(structure.bonds) == 0:
        warnings.warn("No bonds found. Cassandra will interpet "
                "this as a rigid species")
        in_ring = [ False ] * len(structure.atoms)
        frag_list = []
        frag_conn = []
        return in_ring, frag_list, frag_conn

    # Check if entire molecule is connected. Warn if not.
    if nx.is_connected(bond_graph) == False:
        raise ValueError("Not all components of the molecule are connected. "
                         "MCF files are for a single molecule and thus "
                         "everything should be connected through bonds.")

    all_rings = nx.cycle_basis(bond_graph)
    in_ring = [False]*bond_graph.number_of_nodes()
    adj_to_ring = [False]*bond_graph.number_of_nodes()
    for ring in all_rings:
        for idx in ring:
            in_ring[idx] = True

    # Identify fragments
    # See Shah and Maginn, JCP, 135, 134121, 2011, doi:10.1063/1.3644939
    frag_list = []
    frag_conn = []

    # First create a neighbor list for each atom
    neigh_dict = {i:list(bond_graph.neighbors(i))
                  for i in range(bond_graph.number_of_nodes())}
    # First ID fused rings
    fused_rings = []
    rings_to_remove = []
    for i in range(len(all_rings)):
        ring1 = all_rings[i]
        for j in range(i+1, len(all_rings)):
            ring2 = all_rings[j]
            shared_atoms = list(set(ring1) & set(ring2))
            if len(shared_atoms) == 2:
                fused_rings.append(list(set(ring1+ring2)))
                rings_to_remove.append(ring1)
                rings_to_remove.append(ring2)
    for ring in rings_to_remove:
        all_rings.remove(ring)
    all_rings = all_rings + fused_rings
    # ID fragments which contain a ring
    for ring in all_rings:
        adjacentatoms = []
        for idx in ring:
            if len(neigh_dict[idx]) > 2:
                adjacentatoms.append(list(set(neigh_dict[idx])-set(ring)))
        tmp=filter(None, adjacentatoms)
        adjacentatoms = [x for sublist in tmp for x in sublist]
        frag_list.append(ring+adjacentatoms)
        for idx in adjacentatoms:
            adj_to_ring[idx] = True
    # Now ID the other fragments
    for idx in neigh_dict:
        if len(neigh_dict[idx]) > 1:
            if in_ring[idx] == True:
                continue
            else:
                frag_list.append([idx]+neigh_dict[idx])
    # Now find connectivity (shared bonds)
    for i in range(len(frag_list)):
        frag1 = frag_list[i]
        for j in range(i+1, len(frag_list)):
            frag2 = frag_list[j]
            shared_atoms = list(set(frag1) & set(frag2))
            if len(shared_atoms) == 2:
                frag_conn.append([i, j])
            elif len(shared_atoms) > 2:
                warnings.warn('Fragments share more than two atoms...'
                      'something may be going awry unless there are'
                      'fused rings in your system. See below for details.')
                print('Fragment 1 atoms:')
                print(frag1)
                print('Fragment 2 atoms:')
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
    epsilons = [atom.epsilon/IG_CONSTANT_KCAL for atom in structure.atoms]
    sigmas = [atom.sigma for atom in structure.atoms]

    # Check constraints on atom type length and element name length
    # TODO: Update these following Cassandra release
    # to be more reasonable values
    n_unique_elements = len(set(elements))
    for element in elements:
        if len(element) > 2:
            warnings.warn("Warning, element name {} will be shortened "
                 "to two characters. Please confirm your final "
                 "MCF.".format(element))

    elements = [ element[:2] for element in elements ]
    if len(set(elements)) < n_unique_elements:
        warnings.warn("Warning, the number of unique elements has been "
              "reduced due to shortening the element name to two "
              "characters.")

    n_unique_types = len(set(types))
    for itype in types:
        if len(itype) > 6:
            warnings.warn("Warning, type name {} will be shortened to six "
                      "characters as {}. Please confirm your final "
                      "MCF.".format(itype,itype[-6:]))
        types = [ itype[-6:] for itype in types ]
    if len(set(types)) < n_unique_types:
        warnings.warn("Warning, the number of unique atomtypes has been "
              "reduced due to shortening the atomtype name to six "
              "characters.")

    vdw_type = 'LJ'
    header = ('!Atom Format\n'
              '!index type element mass charge vdw_type parameters\n'
              '!vdw_type="LJ", parms=epsilon sigma\n'
              '!vdw_type="Mie", parms=epsilon sigma '
              'repulsion_exponent dispersion_exponent\n'
              '\n# Atom_Info\n'
              )

    mcf_file.write(header)
    mcf_file.write('{:d}\n'.format(len(structure.atoms)))
    for i in range(len(structure.atoms)):
        mcf_file.write('{:<4d}  {:<6s}  {:<2s}  {:7.3f}  {:7.3f}  '
                '{:3s}  {:8.3f}  {:8.3f}'.format(
            i+1, types[i], elements[i], masses[i], charges[i],
            vdw_type, epsilons[i], sigmas[i]))
        if in_ring[i] == True:
            mcf_file.write('  ring')
        mcf_file.write('\n')


def _write_bond_information(mcf_file, structure):
    """Write the bonds in the system.

    Parameters
    ----------
    mcf_file : file object
        The file object of the Cassandra mcf being written
    structure : parmed.Structure
        Parmed structure object

    """

    bond_parms = [ str('{:8.3f}'.format(bond.type.req))
            for bond in structure.bonds ]

    mcf_file.write('\n!Bond Format\n')
    mcf_file.write('!index i j type parameters\n' +
              '!type="fixed", parms=bondLength\n')
    mcf_file.write('\n# Bond_Info\n')
    mcf_file.write('{:d}\n'.format(len(structure.bonds)))
    for i, bond in enumerate(structure.bonds):
        mcf_file.write('{:<4d}  {:<4d}  {:<4d}  {:s}  {:s}\n'.format(
            i+1, bond.atom1.idx+1, bond.atom2.idx+1, 'fixed', bond_parms[i]))

def _write_angle_information(mcf_file, structure, angle_style,
        IG_CONSTANT_KCAL):
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

    if angle_style.casefold() == 'fixed':
        angle_parms = [ str('{:8.2f}'.format(angle.type.theteq))
                        for angle in structure.angles ]
    elif angle_style.casefold() == 'harmonic':
        # Convert energies to units of K
        angle_parms = [ str('{:8.1f}'.format(angle.type.k/IG_CONSTANT_KCAL)) +
                        '  ' + str('{:8.2f}'.format(angle.type.theteq))
                        for angle in structure.angles ]
    else:
        raise ValueError("Only 'fixed' and 'harmonic' angle styles "
                         "are supported by Cassandra")

    header = ( '\n!Angle Format\n'
               '!index i j k type parameters\n'
               '!type="fixed", parms=equilibrium_angle\n'
               '!type="harmonic", parms=force_constant equilibrium_angle\n'
               '\n# Angle_Info\n'
             )

    mcf_file.write(header)
    mcf_file.write('{:d}\n'.format(len(structure.angles)))

    for i, angle in enumerate(structure.angles):
        mcf_file.write('{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:s}  {:s}\n'.format(
            i+1, angle.atom1.idx+1, angle.atom2.idx+1, angle.atom3.idx+1,
            angle_style.lower(), angle_parms[i]))


def _write_dihedral_information(mcf_file, structure, dihedral_style,
        KCAL_TO_KJ):
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
    header = ( '\n!Dihedral Format\n'
               '!index i j k l type parameters\n'
               '!type="none"\n'
               '!type="CHARMM", parms=a0 a1 delta\n'
               '!type="OPLS", parms=c0 c1 c2 c3\n'
               '!type="harmonic", parms=force_constant equilibrium_dihedral\n'
               '\n# Dihedral_Info\n'
             )

    mcf_file.write(header)

    if len(structure.dihedrals) > 0 or len(structure.rb_torsions) > 0:
        if dihedral_style.casefold() == 'opls':
            dihedral_style = dihedral_style.upper()
            dihedrals = structure.rb_torsions
            # Two things happen here:
            #  (1) convert from RB form to Cassandra OPLS
            #  (2) convert units from kcal/mol to kJ/mol
            dihedral_parms = []
            for dihedral in dihedrals:
                a0 = ( dihedral.type.c0 + dihedral.type.c1 +
                        dihedral.type.c2 + dihedral.type.c3 )
                a1 = -dihedral.type.c1 - (3./4.)*dihedral.type.c3
                a2 = (-1./2.)*dihedral.type.c2
                a3 = (-1./4.)*dihedral.type.c3
                if not dihedral.type.c4 == 0. and dihedral.type.c4 == 0.:
                    raise ValueError("Can only convert Ryckaert-Bellemans "
                                     "dihedrals to OPLS if c4==0 and c5==0")
                dihedral_parms.append(str('{:8.3f} '.format(a0*KCAL_TO_KJ)) +
                                      str('{:8.3f} '.format(a1*KCAL_TO_KJ)) +
                                      str('{:8.3f} '.format(a2*KCAL_TO_KJ)) +
                                      str('{:8.3f} '.format(a3*KCAL_TO_KJ)))

        elif dihedral_style.casefold() == 'charmm':
            dihedral_style = dihedral_style.upper()
            dihedrals = structure.dihedrals
            # type.per = periodicity (a1)
            # type.phase = phase offset (delta)
            dihedral_parms = [ str('{:8.3f} '.format(
                                        dihedral.type.phi_k*KCAL_TO_KJ)) +
                               str('{:8.3f} '.format(dihedral.type.per)) +
                               str('{:8.3f}'.format(dihedral.type.phase))
                               for dihedral in dihedrals ]

        elif dihedral_style.casefold() == 'none':
            warnings.warn("Dihedral style 'none' selected. "
                          "Ignoring dihedral parameters")
            dihedral_style = dihedral_style.lower()
            if structure.dihedrals:
                dihedrals = structure.dihedrals
                dihedral_parms = [ '' for dihedral in dihedrals ]
            elif structure.rb_torsions:
                dihedrals = structure.rb_torsions
                dihedral_parms = [ '' for dihedral in dihedrals ]
        else:
            raise ValueError("Only 'OPLS', 'CHARMM', and 'none' "
                             "dihedral styles are supported.")

        mcf_file.write('{:d}\n'.format(len(dihedrals)))
        for i,dihedral in enumerate(dihedrals):
            # Sorting here to match LEAP behavior
            # See https://github.com/choderalab/openmoltools/issues/24
            # If atom types are identical, too bad.
            if dihedral.improper:
                improper_atoms = [dihedral.atom1,
                                  dihedral.atom2,
                                  dihedral.atom4]
                improper_atoms.sort(key=lambda x: x.type)
                atom1 = improper_atoms[0]
                atom2 = improper_atoms[1]
                atom3 = dihedral.atom3
                atom4 = improper_atoms[2]
            else:
                atom1 = dihedral.atom1
                atom2 = dihedral.atom2
                atom3 = dihedral.atom3
                atom4 = dihedral.atom4
            mcf_file.write('{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:<4d}'
                           '  {:s}  {:s}\n'.format(
                i+1, atom1.idx+1, atom2.idx+1, atom3.idx+1, atom4.idx+1,
                dihedral_style, dihedral_parms[i]))
    else:
        mcf_file.write('0\n')

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

    header = ( '\n!Improper Format\n'
               '!index i j k l type parameters\n'
               '!type="harmonic", parms=force_constant equilibrium_improper\n'
               '\n# Improper_Info\n'
             )

    mcf_file.write(header)
    mcf_file.write('{:d}\n'.format(len(structure.impropers)))

    improper_type = 'harmonic'
    for i, improper in enumerate(structure.impropers):
        mcf_file.write('{:<4d}  {:<4d}  {:<4d}  {:<4d}  {:<4d}'
                       '  {:s}  {:8.3f}  {:8.3f}\n'.format(
            i+1, improper.atom1.idx+1, improper.atom2.idx+1,
            improper.atom3.idx+1, improper.atom4.idx+1,
            improper_type, improper.type.psi_k*KCAL_TO_KJ,
            improper.type.psi_eq))

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

    header = ('\n!Fragment Format\n'
              '!index number_of_atoms_in_fragment branch_point other_atoms\n'
              '\n# Fragment_Info\n'
             )

    mcf_file.write(header)

    # Special cases first
    if len(frag_list) == 0:
        if len(structure.atoms) == 1:
            mcf_file.write('1\n')
            mcf_file.write('1 1 1\n')
        elif len(structure.atoms) == 2:
            mcf_file.write('1\n')
            mcf_file.write('1 2 1 2\n')
        else:
            warnings.warn('More than two atoms present but '
                          'no fragments identified.')
            mcf_file.write('0\n')
    else:
        mcf_file.write('{:d}\n'.format(len(frag_list)))
        for i, frag in enumerate(frag_list):
            mcf_file.write('{:d}    {:d}'.format(i+1, len(frag)))
            for idx in frag:
                mcf_file.write('    {:d}'.format(idx+1))
            mcf_file.write('\n')

    mcf_file.write('\n\n# Fragment_Connectivity\n')
    mcf_file.write('{:d}\n'.format(len(frag_conn)))
    for i, conn in enumerate(frag_conn):
        mcf_file.write('{:d}    {:d}    {:d}\n'.format(i+1, conn[0]+1, conn[1]+1))

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

    header = ( '\n!Intra Scaling\n'
               '!vdw_scaling    1-2 1-3 1-4 1-N\n'
               '!charge_scaling 1-2 1-3 1-4 1-N\n'
               '\n# Intra_Scaling\n'
             )

    mcf_file.write(header)
    mcf_file.write('0. 0. {:.4f} 1.\n'.format(lj14))
    mcf_file.write('0. 0. {:.4f} 1.\n'.format(coul14))


