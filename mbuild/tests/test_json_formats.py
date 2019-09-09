from mbuild.tests.base_test import BaseTest
from mbuild.formats.json_formats import compound_to_json, compound_from_json
import mbuild as mb


class TestPB2(BaseTest):

    def test_loop(self, ethane):
        compound_to_json(ethane, 'ethane.json')
        ethane_copy = compound_from_json('ethane.json')
        assert ethane.n_particles == ethane_copy.n_particles
        assert ethane.n_bonds == ethane_copy.n_bonds
        assert len(ethane.children) == len(ethane_copy.children)

    def test_loop_with_ports(self):
        from mbuild.lib.moieties import CH3
        ethane_without_overlap = mb.Compound()
        methyl1 = CH3()
        methyl2 = CH3()
        ethane_without_overlap.add(methyl1, label='methyl1')
        ethane_without_overlap.add(methyl2, label='methyl2')
        compound_to_json(ethane_without_overlap, 'ethane_without_overlap.json', include_ports=True)
        ethane_copy = compound_from_json('ethane_without_overlap.json')
        assert ethane_copy.n_particles == ethane_without_overlap.n_particles
        assert ethane_copy.n_bonds == ethane_without_overlap.n_bonds
        assert len(ethane_copy.children) == len(ethane_without_overlap.children)
        assert len(ethane_copy.all_ports()) == len(ethane_without_overlap.all_ports())
        assert ethane_copy.labels.keys() == ethane_without_overlap.labels.keys()

    def test_loop_for_propyl(self, hexane):
        compound_to_json(hexane, './hexane.json', include_ports=True)
        propyl_copy = compound_from_json('hexane.json')





