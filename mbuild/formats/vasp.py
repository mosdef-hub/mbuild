from itertools import chain

import ele
import mbuild as mb
import numpy as np
from numpy.linalg import norm


__all__ = ['write_poscar', 'read_poscar']

def write_poscar(
        compound, filename, lattice_constant,
        bravais=[[1,0,0],[0,1,0],[0,0,1]], coord='cartesian'
        ):
    """
    Outputs VASP POSCAR files.  See //https://www.vasp.at for
    more information.

    Parameters
    ----------
    compound : mbuild.Compound
               the Compound to write to the POSCAR file
    filename : str
               Path of the output file
    lattice_constant : float
                       Scaling constant for POSCAR file, used to scale all
                       lattice vectors and atomic coordinates
    bravais : array-like
              (3x3) array of bravais cell that defines unit cell of the system
              (default [[1,0,0],[0,1,0],[0,0,1]])
    coord : str
            Coordinate style of atom positions 'cartesian' or 'direct'
            (default 'cartesian')
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
            atom.xx /= lattice_constant
            atom.xy /= lattice_constant
            atom.xz /= lattice_constant

    for atom_name in atom_names:
        atom_count = np.array([atom.name for atom in
            structure.atoms].count(atom_name))
        count_list.append(atom_count)
        xyz = np.array([[atom.xx, atom.xy, atom.xz] for atom in
            structure.atoms if atom.name == atom_name])
        xyz = xyz / 10 # unit conversion from angstroms to nm
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
        data.write(coord+'\n')
        for xyz in xyz_list:
            for pos in xyz:
                data.write('{0:.15f} {1:.15f} {2:.15f}\n'.format(
                    pos[0],pos[1],pos[2]))

def read_poscar(filename, conversion=0.1):
    """
    Reads in a VASP POSCAR or CONTCAR file and returns an mbuild Compound.

    Parameters
    ----------
    filename : str
               path to the POSCAR file
    conversion : float
                 conversion factor multiplied to coordinates when
                 converting between VASP units (angstroms)
                 and mbuild units (nm) (default = 0.1)

    Returns
    -------
    mbuild.Compound
    """

    comp = mb.Compound()

    with open(filename, "r") as f:
        data = f.readlines()

    title = data.pop(0)
    scale = float(data.pop(0).strip())

    a = np.fromiter(data.pop(0).split(), dtype="float64")
    b = np.fromiter(data.pop(0).split(), dtype="float64")
    c = np.fromiter(data.pop(0).split(), dtype="float64")

    lattice_vectors = np.stack((a,b,c))

    # POSCAR files do not require atom types to be specified
    # this block handles unspecified types
    line = data.pop(0).split()
    try:
        n_types = np.fromiter(line, dtype="int")
        types = ["_"+chr(i+64) for i in range(1,len(n_types)+1)]
        # if no types exist, assign placeholder types "_A", "_B", "_C", etc
    except ValueError:
        types = line
        n_types = np.fromiter(data.pop(0).split(), dtype="int")

    total_atoms = np.sum(n_types)
    all_types = list(chain.from_iterable(
        [[itype] * n for itype, n in zip(types,n_types)]
        ))

    # handle optional argument "Selective dynamics"
    # and required arguments "Cartesian" or "Direct"
    switch = data.pop(0)[0].upper()
    selective_dynamics = False # don't know if this is necessary
    if switch == "S":
        selective_dynamics = True
        switch = data.pop(0)[0].upper()

    if switch == "C":
        cartesian = True
    else:
        cartesian = False

    # Slice is necessary to handle files using selective dynamics
    coords = np.stack([np.fromiter(
        line.split()[:3], dtype="float64"
        ) for line in data[:total_atoms]])

    if cartesian:
        coords = coords * scale
    else:
        coords = coords.dot(lattice_vectors) * scale

    alpha = np.rad2deg(np.arccos(b.dot(c)/(norm(b) * norm(c))))
    beta = np.rad2deg(np.arccos(a.dot(c)/(norm(a) * norm(c))))
    gamma = np.rad2deg(np.arccos(a.dot(b)/(norm(a) * norm(b))))

    comp.box = mb.Box(
            lengths=norm(lattice_vectors, axis=1),
            angles=[alpha, beta, gamma]
            )

    for i,xyz in enumerate(coords):
        comp.add(mb.Particle(
            name=all_types[i],
            element=ele.element_from_symbol(all_types[i]),
            pos=xyz*conversion
            ))

    return comp
