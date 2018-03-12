import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestLammpsData(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.lammps')

    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.lammps',
                forcefield_name='oplsaa')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0])) ethane.save(filename='ethane-box.lammps',
        forcefield_name='oplsaa', box=box)

    def test_full_atoms(self, ethane):
        ethane.save(filename='ethane.lammps', atom_style='full')
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" in line:
                    for index in range(2):
                        lines = next(f)
                    columns = lines.split("\t")
                    assert len(columns) == 7

    def test_atomic_atoms(self, ethane):
        ethane.save(filename='ethane.lammps', atom_style='atomic')
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" in line:
                    for index in range(2):
                        lines = next(f)
                    columns = lines.split("\t")
                    assert len(columns) == 5

    def test_molecular_atoms(self, ethane):
        ethane.save(filename='ethane.lammps', atom_style='molecular')
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" in line:
                    for index in range(2):
                        lines = next(f)
                    columns = lines.split("\t")
                    assert len(columns) == 6

    def test_charge_atoms(self, ethane):
        ethane.save(filename='ethane.lammps', atom_style='charge')
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" in line:
                    for index in range(2):
                        lines = next(f)
                    columns = lines.split("\t")
                    assert len(columns) == 6
