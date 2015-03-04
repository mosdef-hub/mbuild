import pytest
from mbuild.tests.base_test import BaseTest


class TestGromacs(BaseTest):


    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    def test_save(self, ethane):
        ethane.save(filename='ethane_out.gro')


