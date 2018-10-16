import numpy as np

from mbuild import Box

__all__ = ['write_vasp']

def write_vasp(compound, filename, lattice_constant, bravais,
        sel_dev=None,coord='cartesian'):
    """Output VASP input files.

    Outputs VASP POSCAR files.  See //https://www.vasp.at for
    more information.

    Parameters
    ----------
    compound: mb.Compound
        mBuild Compound
    filename: str
        Path of the output file
    lattice_constant: float
        Scaling constant for POSCAR file
    bravais: array
        array of bravais cell
    sel_dev: yes
    coord: str, default = 'cartesian', other option = 'direct'
        Coordinate style of atom positions
    """
    structure = compound.to_parmed()
    atom_names = np.unique([atom.name for atom in structure.atoms])
    count_list = list()
    xyz_list = list()

    """
    Coordinates are broken up into a list of np.arrays to ensure
    that the coordinates of the first atom listed are written to the
    file first
    """
    if coord == 'direct':
        for atom in structure.atoms:
            atom.xx = atom.xx / lattice_constant
            atom.xy = atom.xy / lattice_constant
            atom.xz = atom.xz / lattice_constant
   
    for atom_name in atom_names:
        atom_count = np.array([atom.name for atom in
            structure.atoms].count(atom_name))
        count_list.append(atom_count)
        xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in
            structure.atoms if atom.name == atom_name])
        xyz = xyz / 10
        xyz_list.append(xyz)
    
    with open(filename, 'w') as data:
        data.write(filename+' - created by mBuild\n')
        data.write('     {0:.15f}\n'.format(lattice_constant))
        data.write('    ')
        for item in bravais[0]:
            data.write(' {0:.15f}'.format(item))
        data.write('\n')
        data.write('    ')
        for item in bravais[1]:
            data.write(' {0:.15f}'.format(item))
        data.write('\n')
        data.write('    ')
        for item in bravais[2]:
            data.write(' {0:.15f}'.format(item))
        data.write('\n')
        data.write('{}\n'.format('   '.join(map(str,atom_names))))
        data.write('{}\n'.format('   '.join(map(str,count_list))))
        if sel_dev:
            data.write('Selective Dyn\n')
        data.write(coord+'\n')
        for xyz in xyz_list:
            for pos in xyz:
                data.write('{0:.15f} {1:.15f} {2:.15f}\n'.format(
                    pos[0],pos[1],pos[2]))
