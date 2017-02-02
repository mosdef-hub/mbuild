import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPattern(BaseTest):

    def test_port_center(self):
        port = mb.Port()
        assert [round(coord,3)==0 for coord in port.center]

    def test_port_shift(self, ethane):
        mb.translate_to(ethane, np.ones(3))
        port = mb.Port(anchor=ethane)
        assert [ethane.center[i]==coord for i,coord in enumerate(port.center)]
