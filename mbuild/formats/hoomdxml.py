__all__ = ['write_hoomdxml']


from foyer.forcefield import apply_forcefield
import numpy as np

def write_hoomdxml(compound, filename, forcefield=None, ref_distance=1.0, ref_mass=1.0, box=None):
    """Output a HOOMD XML file.
    
    Parameters
    ----------
    compound : mb.Compound
        Compound to be saved
    filename : str
        Path of the output file.
    forcefield : str, default=None
        Name of the force field to be applied to the compound
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, default=1.0
        Reference mass for conversion to reduced units
    box : mb.Box, default=None
        Box information to save to XML file, if not specified
        a bounding box is used
    
    """

    structure = compound.to_parmed()
    if forcefield:
        structure = apply_forcefield(structure, forcefield=forcefield)

    if not box:
        box = compound.boundingbox

    with open(filename, 'w') as xml_file:
        xml_file.write('<?xml version="1.2" encoding="UTF-8"?>\n')
        xml_file.write('<hoomd_xml version="1.2">\n')
        xml_file.write('<configuration time_step="0">\n')
        xml_file.write('<box units="sigma"  Lx="{}" Ly="{}" Lz="{}"/>\n'.format(*box.lengths/ref_distance))
    
        xml_file.write('<position units="sigma" num="{}">\n'.format(len(structure.atoms)))

        xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])
        for pos in xyz:
            xml_file.write('{}\t{}\t{}\n'.format(*pos/ref_distance))
        xml_file.write('</position>\n')
        
        if forcefield:
            types = [atom.type for atom in structure.atoms]
        else:
            types = [atom.name for atom in structure.atoms]
        xml_file.write('<type>\n')
        for atom_type in types:
            xml_file.write('{}\n'.format(atom_type))
        xml_file.write('</type>\n')

        masses = [atom.mass for atom in structure.atoms]
        xml_file.write('<mass>\n')
        for mass in masses:
            xml_file.write('{}\n'.format(mass/ref_mass)) 
        xml_file.write('</mass>\n')
        
        charges = [atom.charge for atom in structure.atoms]
        xml_file.write('<charge>\n')
        for charge in charges:
            xml_file.write('{}\n'.format(charge))
        xml_file.write('</charge>\n')

        bonds = [[bond.atom1.idx, bond.atom2.idx] for bond in structure.bonds] 
        if bonds:
            if len(set([bond.type for bond in structure.bonds])) == 1:
                bond_types = np.zeros(len(bonds),dtype=int)
            else:
                all_bond_types = dict(enumerate(set([(round(bond.type.k,3),
                                                      round(bond.type.req,3)) for bond in structure.bonds])))
                all_bond_types = {y:x for x,y in all_bond_types.iteritems()}
                bond_types = [all_bond_types[(round(bond.type.k,3),
                                              round(bond.type.req,3))] for bond in structure.bonds]
            xml_file.write('<bond>\n')
            for idx,bond in enumerate(bonds):
                xml_file.write('{}\t{}\t{}\n'.format(bond_types[idx],*bond))
            xml_file.write('</bond>\n')

        angles = [[angle.atom1.idx, 
                   angle.atom2.idx, 
                   angle.atom3.idx] for angle in structure.angles]
        if angles:
            all_angle_types = dict(enumerate(set([(round(angle.type.k,3), 
                                                   round(angle.type.theteq,3)) for angle in structure.angles])))
            all_angle_types = {y:x for x,y in all_angle_types.iteritems()}
            angle_types = [all_angle_types[(round(angle.type.k,3), 
                                            round(angle.type.theteq,3))] for angle in structure.angles]
            xml_file.write('<angle>\n')
            for idx,angle in enumerate(angles):
                xml_file.write('{}\t{}\t{}\t{}\n'.format(angle_types[idx],*angle))
            xml_file.write('</angle>\n')

        dihedrals = [[dihedral.atom1.idx,
                      dihedral.atom2.idx,
                      dihedral.atom3.idx,
                      dihedral.atom4.idx] for dihedral in structure.rb_torsions]
        if dihedrals:
            all_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0,3),
                                                      round(dihedral.type.c1,3),
                                                      round(dihedral.type.c2,3),
                                                      round(dihedral.type.c3,3),
                                                      round(dihedral.type.c4,3),
                                                      round(dihedral.type.c5,3),
                                                      round(dihedral.type.scee,1),
                                                      round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
            all_dihedral_types = {y:x for x,y in all_dihedral_types.iteritems()}
            dihedral_types = [all_dihedral_types[(round(dihedral.type.c0,3),
                                                  round(dihedral.type.c1,3),
                                                  round(dihedral.type.c2,3),
                                                  round(dihedral.type.c3,3),
                                                  round(dihedral.type.c4,3),
                                                  round(dihedral.type.c5,3),
                                                  round(dihedral.type.scee,1),
                                                  round(dihedral.type.scnb,1))] for dihedral in structure.rb_torsions]
            xml_file.write('<dihedral>\n')
            for idx,dihedral in enumerate(dihedrals):
                xml_file.write('{}\t{}\t{}\t{}\t{}\n'.format(dihedral_types[idx],
                                                             *dihedral))
            xml_file.write('</dihedral>\n')

        impropers = [[improper.atom1.idx,
                      improper.atom2.idx,
                      improper.atom3.idx,
                      impropers.atom4.idx] for improper in structure.impropers]
        if impropers:
            all_improper_types = dict(enumerate(set([(round(improper.type.psi_k,3),
                                                      round(improper.type.psi_eq,3)) for improper in structure.impropers])))
            all_improper_types = {y:x for x,y in all_improper_types.iteritems()}
            improper_types = [all_improper_types[(round(improper.type.psi_k,3),
                                                  round(improper.type.psi_eq,3))] for improper in structure.impropers]
            xml_file.write('<improper>\n')
            for idx,improper in enumerate(impropers):
                xml_file.write('{}\t{}\t{}\t{}\t{}\n'.format(improper_types[idx],
                                                             *improper))
            xml_file.write('</improper>\n')
            
        '''
        # Rigid bodies
        xml_file.write('<body>\n')
        xml_file.write('</body>\n')
        '''
        
        xml_file.write('</configuration>\n')
        xml_file.write('</hoomd_xml>')
