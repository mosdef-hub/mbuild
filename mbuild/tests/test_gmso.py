import numpy as np
import pytest

import mbuild as mb
from mbuild.conversion import from_gmso, to_gmso
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_gmso, import_

if has_gmso:
    gmso = import_("gmso")


@pytest.mark.skipif(not has_gmso, reason="GMSO is not installed")
class TestGMSO(BaseTest):
    def test_to_gmso(self, ethane):
        gmso_eth = ethane.to_gmso()  # Equivalent to to_gmso(ethane)

        assert isinstance(gmso_eth, gmso.Topology)
        assert len(gmso_eth.bonds) == ethane.n_bonds
        for i in range(ethane.n_particles):
            assert ethane[i].name == gmso_eth.sites[i].name
            assert np.isclose(
                ethane[i].xyz, gmso_eth.sites[i].position.value
            ).all()

    def test_full_conversion(self, ethane):
        # Note: at this point, the full conversion may lose some information regarding the hierarchical,
        # especially, if the original compound has more than 3 layers.
        gmso_eth = ethane.to_gmso()
        mb_eth = from_gmso(gmso_eth)

        assert mb_eth.n_particles == ethane.n_particles
        assert mb_eth.n_bonds == ethane.n_bonds
        for i in range(mb_eth.n_particles):
            assert mb_eth[i].name == ethane[i].name
            assert np.isclose(mb_eth[i].xyz, ethane[i].xyz).all()

    def test_coords_only(self, ethane):
        gmso_eth = ethane.to_gmso()
        mb_eth = from_gmso(gmso_eth)
        # Reset coord of the mb_eth
        for particle in mb_eth.particles():
            particle.xyz = [[0, 0, 0]]

        mb_eth.from_gmso(gmso_eth, coords_only=True)
        for i in range(mb_eth.n_particles):
            assert np.isclose(mb_eth[i].xyz, ethane[i].xyz).all()

    def test_mismatch_coords_only(self, ethane):
        gmso_eth = ethane.to_gmso()
        meth = mb.load("C", smiles=True)
        with pytest.raises(ValueError):
            meth.from_gmso(gmso_eth, coords_only=True)

    def test_infer_hier(selfs, ethane):
        # Create an ethane box, should be a four structure
        eth_box = mb.packing.fill_box(compound=ethane, n_compounds=1, density=1)

        #TODO: once infer_hierarchy is implemented in gmso, changes these functions

        # infer_hierarchy=False
        #unlabeled_top = eth_box.to_gmso(infer_hierarchy=False)
        unlabeled_top = eth_box.to_gmso()
        for site in unlabeled_top.sites:
            assert not site.residue or site.molecule or site.group

        # infer_hierarchy=True
        #labeled_top = eth_box.to_gmso(infer_hierarchy=True)
        labeled_top = eth_box.to_gmso()
        for site in labeled_top.sites:
            assert site.group == "Ethane"
            assert site.molecule.name == "Ethane"
            assert site.residue.name == "CH3"
