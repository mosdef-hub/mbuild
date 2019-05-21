from mbuild.tests.base_test import BaseTest
from mbuild.formats.protobuf import write_pb2, read_pb2

class TestPB2(BaseTest):

    def test_loop(self, ethane):
        write_pb2(ethane, 'ethane.pb2')
        proto = read_pb2('ethane.pb2')
        assert ethane.n_particles == proto.n_particles
        assert ethane.n_bonds == proto.n_bonds
        assert len(ethane.children) == len(proto.children)
