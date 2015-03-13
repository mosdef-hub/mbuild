import pytest

from mbuild.tests.base_test import BaseTest


class TestGromacs(BaseTest):

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_write(self, ethane):
        ethane.save('ethane.top')

    def test_update_from_file(self, ethane):
        pass

    def test_read_write_compare(self, ethane):
        pass
