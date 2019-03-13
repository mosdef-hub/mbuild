import mbuild as mb
import parmed as pmd
import warnings
from foyer import Forcefield

__all__ = ['write_par']

def write_par(structure, filename):
    """ Write CHARMM Par file given a parametrized structureture """

    # ATOMS
    with open(filename, 'w') as f:
        f.write("ATOMS\n")
        for atom in structure.atoms:
            f.write("MASS -1 {:8s} {:8.4f}\n".format(atom.atom_type.name, atom.mass))

        f.write("\nBONDS\n")
        unique_bonds = set()
        for bond_type in structure.bond_types:
            for bond in structure.bonds:
                if bond.type == bond_type:
                    unique_bonds.add((bond.atom1.atom_type.name, 
                        bond.atom2.atom_type.name, bond.type))
                    
        for bond in unique_bonds:                
            f.write('{:8s} {:8s} {:.5f} {:.5f}\n'.format(bond[0], bond[1],
                                                bond[2].k, bond[2].req))
                
        f.write("\nANGLES\n")
        unique_angles = set()
        for angle_type in structure.angle_types:
            for angle in structure.angles:
                if angle.type == angle_type:
                    unique_angles.add((angle.atom1.atom_type.name, 
                                        angle.atom2.atom_type.name, 
                                       angle.atom3.atom_type.name, angle.type))
                    
        for angle in unique_angles:                
            f.write('{:8s} {:8s} {:8s} {:.5f} {:.5f}\n'.format(angle[0], angle[1], 
                                                angle[2],
                                                angle[3].k, angle[3].theteq))
            
        # TODO :UREY-BRADLEYS in structure
        #unique_ubs = set()
        #for ub_type in structure.ub_types:
        #    for ub in structure.ubs:
        #        if ub.type == ub_type:
        #            unique_ubs.add((ub.atom1.atom_type.name, 
        #                                ub.atom2.atom_type.name, 
        #                               ub.atom3.atom_type.name, ub.type))
        #            
        #for ub in unique_ubs:                
        #    f.write('{:8s} {:8s} {:8s} {:.5f} {:.5f}\n'.format(ub[0], ub[1], 
        #                                        ub[2],
        #                                        ub[3].k, ub[3].theteq))

        # These dihedrals need to be put in PeriodicTorsion Style (Charmm style)
        # IF there are RB (OPLS-style) torsions, we will need to convert them
        f.write("\nDIHEDRALS\n")
        unique_dihedrals = set()
        scnb = set()
        for dihedral_type in structure.dihedral_types:
            for dihedral in structure.dihedrals:
                if dihedral.type == dihedral_type:
                    unique_dihedrals.add((dihedral.atom1.atom_type.name, 
                                         dihedral.atom2.atom_type.name,
                                         dihedral.atom3.atom_type.name, 
                                         dihedral.atom4.atom_type.name,
                                         dihedral.type))
                    scnb.add(dihedral.type.scnb)
                    
        for dihedral in unique_dihedrals:                
            f.write('{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n'.format(
                dihedral[0], dihedral[1], dihedral[2], dihedral[3],
                dihedral[4].phi_k, dihedral[4].per, dihedral[4].phase))


        f.write("\nIMPROPERS\n")
        unique_impropers = set()
        for improper_type in structure.improper_types:
            for improper in structure.impropers:
                if improper.type == improper_type:
                    unique_impropers.add((improper.atom1.atom_type.name, 
                        improper.atom2.atom_type.name,
                        improper.atom3.atom_type.name, 
                        improper.atom4.atom_type.name, improper.type))
        for improper in unique_impropers:                
            f.write('{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n'.format(
                improper[2], improper[0], improper[1], improper[3],
                improper[4].psi_k, 0, improper[4].psi_eq))

        # TODO additional nonbonded parameters
        sc_nb = [a for a in scnb]
        if len(sc_nb) > 1:
            raise ValueError("Multiple 1-4 LJ scalings were detected")
        elif len(sc_nb) == 0:
            sc_nb = [1]
        sc_nb = [1]

        f.write("\nNONBONDED\n")
        unique_atypes = set()
        for atom in structure.atoms:
            unique_atypes.add(atom.atom_type)
        for atype in unique_atypes:
            # atype, 0.0, epsilon, rmin/2
            #f.write('{:8s} {:8.3f} {:8f} {:8f}\n'.format(atype.name,
                #0.0, -1 *atype.epsilon, atype.rmin/2))

            # atype, 0.0, epsilon, rmin/2, 0.0, epsilon(1-4), rmin/2 (1-4)
            f.write('{:8s} {:8.3f} {:8f} {:8f} {:8f} {:8f} {:8f}\n'.format(
                atype.name,
                0.0, -1 *atype.epsilon, atype.rmin/2, 0.0,
                -1*sc_nb[0]*atype.epsilon, atype.rmin/2))

        if structure.has_NBFIX(): 
            warnings.warn("NBFixes detected but unsupported in .par writer")

        f.write("\nEND")
