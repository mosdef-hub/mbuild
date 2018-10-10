import numpy as np
import pytest

import mbuild as mb
from mbuild.formats.vasp import write_vasp
from mbuild.tests.base_test import BaseTest

class TestLammpsData(BaseTest):
    """
    Unit tests for Vasp writer
    """

    def test_write(self):
        copper = mb.Compound(name='Cu')
        lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        spacing = [.36149, .36149, .36149]
        copper_locations = [[0., 0., 0.], [.5, .5, 0.],
                [.5, 0., .5], [0., .5, .5]]
        basis = {'Cu': copper_locations}
        copper_lattice = mb.Lattice(lattice_spacing = spacing,
                lattice_vectors=lattice_vector, lattice_points=basis)
        copper_dict = {'Cu' : copper}
        copper_pillar = copper_lattice.populate(x=3, y=3, z=1,
                compound_dict=copper_dict)
        struct = copper_pillar.to_parmed()
        write_vasp(struct, 'test.poscar',
                lattice_constant=spacing[0], bravais=lattice_vector,
                coord='cartesian')

    def test_cartesian(self):
        coord_list = list()
        poscar_coord_list = list()
        copper = mb.Compound(name='Cu')
        lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        spacing = [.36149, .36149, .36149]
        copper_locations = [[0., 0., 0.], [.5, .5, 0.],
                [.5, 0., .5], [0., .5, .5]]
        basis = {'Cu': copper_locations}
        copper_lattice = mb.Lattice(lattice_spacing = spacing,
                lattice_vectors=lattice_vector, lattice_points=basis)
        copper_dict = {'Cu' : copper}
        copper_pillar = copper_lattice.populate(x=3, y=3, z=1,
                compound_dict=copper_dict)
        for child in copper_pillar.children:
            coord_list.append(child.xyz)
        struct = copper_pillar.to_parmed()
        write_vasp(struct, 'test.poscar',
                lattice_constant=spacing[0], bravais=lattice_vector,
                coord='cartesian')
        with open('test.poscar', 'r') as f:
            f = f.readlines()[8:]
            for line in f:
                poscar_coord_list.append(
                        np.genfromtxt(line.splitlines(True)))
        
        assert all([a.all() == b.all() for a,b in zip(
            poscar_coord_list, coord_list)])

                

    #def test_direct(self):


    #def test_file_length(self):
