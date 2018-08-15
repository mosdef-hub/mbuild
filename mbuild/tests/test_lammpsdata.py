import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.formats.lammpsdata import write_lammpsdata


class TestLammpsData(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.lammps')

    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.lammps', forcefield_name='oplsaa')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.lammps', forcefield_name='oplsaa', box=box)

    def test_nbfix(self, ethane):
        from foyer import Forcefield
        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
	# Add nbfixes
        types = list(set([a.atom_type for a in structure.atoms]))
        types[0].add_nbfix(types[1].name, 1.2, 2.1)
        types[1].add_nbfix(types[0].name, 1.2, 2.1)
        write_lammpsdata(filename='nbfix.lammps', structure=structure)

    @pytest.mark.parametrize('atom_style, n_columns', [('full', 7), ('atomic', 5), ('molecular', 6), ('charge', 6)])
    def test_writing_atom_styles(self, ethane, atom_style, n_columns):
        ethane.save(filename='ethane.lammps', atom_style=atom_style)
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" not in line:
                    continue
                    atoms_header = next(f)
                    first_atom_line = next(f)
                    columns = first_atom_line.split("\t")
                    assert len(columns) == n_columns
