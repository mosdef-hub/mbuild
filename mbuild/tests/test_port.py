import pytest
import warnings
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
        ethane.translate_to(np.ones(3))
        port = mb.Port(anchor=ethane, separation=0)
        assert [ethane.center[i]==coord for i,coord in enumerate(port.center)]

    def test_port_init_shift(self, ethane):
        ethane.translate_to(np.ones(3))
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
        port.rotate(np.pi, [1, 0, 0])
        assert(np.allclose([0, -1, 0], port.direction, atol=1e-15))

    def test_port_separation(self, ethane):
        port = mb.Port(anchor=ethane, separation=0.7)
        assert(np.allclose(0.7, port.separation, atol=1e-15))

        port_no_anchor = mb.Port()
        with warnings.catch_warnings(record=True) as w:
            separation = port_no_anchor.separation
            assert w
            assert(separation is None)

    def test_update_separation(self, ethane, hexane):
        port = mb.Port(anchor=ethane, separation=0.7)
        port.update_separation(separation=0.9)
        assert(np.allclose(0.9, port.separation, atol=1e-15))

        port_used = hexane.labels['propyl2'].labels['down']
        with warnings.catch_warnings(record=True) as w:
            port_used.update_separation(0.10)
            assert w

    def test_update_orientaiton(self, ch2, hexane):
        port = ch2['up']
        port.update_orientation(orientation=(1,0,0))
        assert(np.allclose([-1,0,0], port.direction, atol=1e-15))

        port_used = hexane.labels['propyl2'].labels['down']
        with warnings.catch_warnings(record=True) as w:
            port_used.update_separation(0.10)
            assert w

    def test_access_labels(self):
        port = mb.Port()
        compound = mb.Compound()
        compound.add(port, label='foo')
        assert port.access_labels == ["['foo']"]

        compound2 = mb.Compound(name='C2')
        compound2.add(compound, label='bar')
        assert port.access_labels == ["['bar']['foo']"]

    def test_up_down_reverse_orientation_axes(self):
        for vector in [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
            port1 = mb.Port(orientation=vector)
            port2 = mb.Port(orientation=-np.array(vector))
            assert np.allclose(port1['up'].xyz_with_ports,
                               port2['down'].xyz_with_ports)
            assert np.allclose(port1['down'].xyz_with_ports,
                               port2['up'].xyz_with_ports)

    def test_up_down_reverse_orientation_random(self):
        np.random.seed(84)
        for _ in range(15):
            vector = np.random.random(3) - 0.5
            port1 = mb.Port(orientation=vector)
            port2 = mb.Port(orientation=-vector)
            assert np.allclose(port1['up'].xyz_with_ports,
                               port2['down'].xyz_with_ports)
            assert np.allclose(port1['down'].xyz_with_ports,
                               port2['up'].xyz_with_ports)
