from __future__ import division

import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilica
from mbuild.tests.base_test import BaseTest


class TestSilicaInterface(BaseTest):

    def test_silica_interface(self):
        tile_x = 1
        tile_y = 1
        thickness = 1.2

        interface = mb.SilicaInterface(bulk_silica=AmorphousSilica(),
                                       tile_x=tile_x,
                                       tile_y=tile_y,
                                       thickness=thickness)

        thickness_tolerance = 0.31
        z = [atom.pos[2] for atom in interface.particles()
             if atom.name in ['Si', 'O']]
        assert abs(max(z) - min(z) - thickness) < thickness_tolerance

        density_tolerance = 0.1
        area = interface.periodicity[0] * interface.periodicity[1]
        oh_count = len(list(interface.particles_by_name('OS')))
        assert abs((oh_count/area) - 5.0) < density_tolerance
