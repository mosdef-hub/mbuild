import pytest
import numpy as np
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPattern(BaseTest):

    def test_port_center(self):
        port = mb.Port()
        assert [round(coord,3)==0 for coord in port.center]

    def test_port_shift(self, ethane):
        ethane.translate_to(np.ones(3))
        port = mb.Port(anchor=ethane)
        assert [ethane.center[i]==coord for i,coord in enumerate(port.center)]

    def test_port_init_shift_0(self, ethane):
        mb.translate_to(ethane, np.ones(3))
        port = mb.Port(anchor=ethane, separation=0)
        assert [ethane.center[i]==coord for i,coord in enumerate(port.center)]

    def test_port_init_shift(self, ethane):
        mb.translate_to(ethane, np.ones(3))
        separation = [1, 2, 3]
        port = mb.Port(anchor=ethane, separation=separation)
        assert [ethane.center[i]+separation[i]==coord 
                for i, (coord, sep) in enumerate(zip(port.center, separation))]

    def test_port_init_rotate_0(self, ethane):
        port1 = mb.Port(anchor=ethane, orientation=[0, 1, 0])
        port2 = mb.Port(anchor=ethane)
        assert np.allclose(port1.xyz_with_ports, port2.xyz_with_ports, atol=1e-15)

    def test_port_init_rotate(self, ethane):
        port1 = mb.Port(anchor=ethane, orientation=[1, 1, 1])
        port2 = mb.Port(anchor=ethane)
        assert not np.allclose(port1.xyz_with_ports, 
                               port2.xyz_with_ports,
                               atol=1e-15)

    def test_port_init_rotate_180(self, ethane):
        port1 = mb.Port(anchor=ethane, orientation=[0, -1, 0])
        port2 = mb.Port(anchor=ethane)
        assert not np.allclose(port1.xyz_with_ports, 
                               port2.xyz_with_ports,
                               atol=1e-15)

    def test_port_direction(self):
        port = mb.Port()
        assert(np.allclose([0, 1, 0], port.direction, atol=1e-16))
        mb.coordinate_transform.rotate(port, np.pi, [1, 0, 0])
        assert(np.allclose([0, -1, 0], port.direction, atol=1e-15))
