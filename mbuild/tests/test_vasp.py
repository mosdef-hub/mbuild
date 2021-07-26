import numpy as np
import pytest

import mbuild as mb
from mbuild.formats.vasp import read_poscar, write_poscar
from mbuild.tests.base_test import BaseTest


class TestVasp(BaseTest):
    """Unit tests for Vasp POSCAR writer"""

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    def test_read_write(self, gilmerite):
        write_poscar(gilmerite, "test.poscar")
        new_gilmerite = read_poscar("test.poscar")
        assert np.allclose(gilmerite.box.lengths, new_gilmerite.box.lengths)
        assert np.allclose(gilmerite.box.angles, new_gilmerite.box.angles)
        assert np.allclose(gilmerite.xyz, new_gilmerite.xyz)

    def test_read_write_direct(self, gilmerite):
        write_poscar(gilmerite, "test.poscar", coord_style="direct")
        new_gilmerite = read_poscar("test.poscar")
        assert np.allclose(gilmerite.box.lengths, new_gilmerite.box.lengths)
        assert np.allclose(gilmerite.box.angles, new_gilmerite.box.angles)
        assert np.allclose(gilmerite.xyz, new_gilmerite.xyz)

    def test_lattice_constant(self, copper_cell):
        write_poscar(copper_cell, "test.poscar", lattice_constant=0.4123)
        with open("test.poscar", "r") as f:
            for i, line in enumerate(f):
                if i == 1:
                    lattice_constant = np.genfromtxt(line.splitlines(True))

        assert lattice_constant == 0.4123

    def test_bravais(self, copper_cell):
        """Test that compound with no box has a lattice that is diagonal."""
        write_poscar(copper_cell, "test.poscar")
        with open("test.poscar", "r") as f:
            lines = f.readlines()

        bravais = np.stack(
            [np.fromstring(line, sep=" ") for line in lines[2:5]]
        )

        # zero the diagonal
        for i in range(3):
            bravais[i, i] = 0
        assert np.array_equal(bravais, np.zeros((3, 3)))

    def test_num_elements(self, cscl_crystal):
        write_poscar(cscl_crystal, "test.poscar")
        with open("test.poscar", "r") as f:
            for i, line in enumerate(f):
                if i == 5:
                    elements = line.split()

        assert len(elements) == 2

    def test_num_atoms(self, copper_cell):
        write_poscar(copper_cell, "test.poscar")
        with open("test.poscar", "r") as f:
            for i, line in enumerate(f):
                pass
        assert i + 1 == 44

    @pytest.mark.parametrize("coord_type", ["direct", "cartesian"])
    def test_coordinate_header(self, gilmerite, coord_type):
        write_poscar(gilmerite, "test.poscar", coord_style=coord_type)
        with open("test.poscar", "r") as f:
            for i, line in enumerate(f):
                if i == 7:
                    coord = line.strip()

        assert coord == coord_type

    def test_warning_raised(self, copper_cell):
        copper_cell.box = None
        with pytest.warns(UserWarning):
            write_poscar(copper_cell, "test.poscar", coord_style="direct")

    def test_error_raised(self, copper_cell):
        with pytest.raises(ValueError):
            write_poscar(copper_cell, "test.poscar", coord_style="heck")
