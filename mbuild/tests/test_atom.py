from mbuild.tests.base_test import BaseTest


class TestAtom(BaseTest):

    def test_neighbors(self, ethane):
        neighbors = ethane.atoms[0].neighbors
        assert len(neighbors) == 4
        neighbor_names = [neighbor.name for neighbor in neighbors]
        assert neighbor_names.count('H') == 3
        assert neighbor_names.count('C') == 1
