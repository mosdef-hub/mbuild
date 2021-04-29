import numpy as np
import pytest

import mbuild as mb
from mbuild.formats.vasp import write_poscar
from mbuild.tests.base_test import BaseTest

class TestVasp(BaseTest):
    """
    Unit tests for Vasp POSCAR writer
    """

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    def test_write(self, copper_cell):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=.4123)


    def test_write_direct(self, copper_cell):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=.4123, coord='direct')


    def test_lattice_constant(self, copper_cell):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=.4123)
        with open('test.poscar', 'r') as f:
            for i,line in enumerate(f):
                if i == 1:
                    lattice_constant = np.genfromtxt(
                            line.splitlines(True))

        assert lattice_constant == 0.4123


    def test_bravais(self, copper_cell):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=.4123)
        with open('test.poscar', 'r') as f:
            bravais = list()
            for i,line in enumerate(f):
                if i in [2,3,4]:
                    bravais.append(np.genfromtxt(
                            line.splitlines(True)))

        assert all([a.all() == b.all() for a, b in zip(bravais,
            [np.array([1,0,0]), np.array([0,1,0]),
                np.array([0,0,1])])])


    def test_num_elements(self, cscl_crystal):
        write_poscar(cscl_crystal, 'test.poscar',
                lattice_constant=.4123)
        with open('test.poscar', 'r') as f:
            for i,line in enumerate(f):
                if i == 5:
                    elements = line.split()

        assert len(elements) == 2


    def test_num_atoms(self, copper_cell):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=0.4123)
        with open('test.poscar', 'r') as f:
            for i, line in enumerate(f):
                pass

        assert i + 1 == 44


    @pytest.mark.parametrize('coord_type', ['direct', 'cartesian'])
    def test_coordinate_header(self, copper_cell, coord_type):
        write_poscar(copper_cell, 'test.poscar',
                lattice_constant=0.4123, coord=coord_type)
        with open('test.poscar', 'r') as f:
            for i, line in enumerate(f):
                if i == 7:
                   coord = line.strip() 

        assert coord == coord_type
