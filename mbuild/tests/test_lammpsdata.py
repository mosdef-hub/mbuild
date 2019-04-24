import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn
from mbuild.utils.io import has_foyer


class TestLammpsData(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.lammps')

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.lammps', forcefield_name='oplsaa')

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_charmm(self):
        cmpd = mb.load(get_fn('charmm_dihedral.mol2'))
        for i in cmpd.particles():
            i.name = "_{}".format(i.name)
        structure = cmpd.to_parmed(box=cmpd.boundingbox, 
                                    residues=set([p.parent.name for \
                                                 p in cmpd.particles()]))

        from foyer import Forcefield
        ff = Forcefield(forcefield_files=[get_fn('charmm_truncated.xml')])
        structure = ff.apply(structure, assert_dihedral_params=False)

        from mbuild.formats.lammpsdata import write_lammpsdata
        write_lammpsdata(structure, 'charmm_dihedral.lammps')

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.lammps', forcefield_name='oplsaa', box=box)

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_triclinic_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]), angles=[60, 70, 80])
        ethane.save(filename='triclinic-box.lammps', forcefield_name='oplsaa', box=box)

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
