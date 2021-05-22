import numpy as np
import pytest

import mbuild as mb
from mbuild.formats.par_writer import write_par
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestPar(BaseTest):
    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_charmm(self):
        cmpd = mb.load(get_fn("charmm_dihedral.mol2"))
        for i in cmpd.particles():
            i.name = "_{}".format(i.name)
        structure = cmpd.to_parmed(
            box=cmpd.get_boundingbox(),
            residues=set([p.parent.name for p in cmpd.particles()]),
        )

        from foyer import Forcefield

        ff = Forcefield(forcefield_files=[get_fn("charmm_truncated.xml")])
        structure = ff.apply(structure, assert_dihedral_params=False)

        write_par(structure, "charmm_dihedral.par")

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename="ethane-opls.par", forcefield_name="oplsaa")

    def test_par_parameters(self, ethane):
        ethane.save(filename="ethane-opls.par", forcefield_name="oplsaa")
        from parmed.charmm import CharmmParameterSet

        pset = CharmmParameterSet.load_set(pfile="ethane-opls.par")
        assert len(pset.bond_types) == 3
        assert len(pset.angle_types) == 3
        assert len(pset.atom_types) == 2
