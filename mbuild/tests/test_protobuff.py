from mbuild.tests.base_test import BaseTest
from mbuild.formats.protobuff import write_pb3, read_pb3

class TestPB3(BaseTest):

    def test_loop(self, ethane):
        write_pb3(ethane, 'ethane.pb3')
        proto = read_pb3('ethane.pb3')
        assert ethane.n_particles == proto.n_particles
        assert ethane.n_bonds == proto.n_bonds
        assert len(ethane.children) == len(proto.children)
