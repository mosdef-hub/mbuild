import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_protobuf
if has_protobuf:
    from mbuild.formats.protobuf import write_pb2, read_pb2

class TestPB2(BaseTest):

    @pytest.mark.skipif(not has_protobuf, reason="Protobuf package not installed")
    def test_loop(self, ethane):
        write_pb2(ethane, 'ethane.pb2')
        proto = read_pb2('ethane.pb2')
        assert ethane.n_particles == proto.n_particles
        assert ethane.n_bonds == proto.n_bonds
        assert len(ethane.children) == len(proto.children)
