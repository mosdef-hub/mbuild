from __future__ import division

import numpy as np

import mbuild as mb
from mbuild.lib.bulk_materials import AmorphousSilicaBulk
from mbuild.lib.recipes import SilicaInterface
from mbuild.tests.base_test import BaseTest


class TestSilicaInterface(BaseTest):

    def test_silica_interface(self):
        tile_x = 1
        tile_y = 1
        thickness = 0.6

        interface = SilicaInterface(bulk_silica=AmorphousSilicaBulk(),
                                    tile_x=tile_x,
                                    tile_y=tile_y,
                                    thickness=thickness)

        thickness_tolerance = 0.05
        z = [atom.pos[2] for atom in interface.particles()
             if atom.name == 'Si']
        assert abs(max(z) - min(z) - thickness) < thickness_tolerance

        density_tolerance = 0.1
        area = interface.periodicity[0] * interface.periodicity[1]
        oh_count = len(list(interface.particles_by_name('O_surface')))
        assert abs((oh_count/area) - 5.0) < density_tolerance

    def test_seed(self):
        tile_x = 1
        tile_y = 1
        thickness = 0.6
        seed = 12345

        interface1 = SilicaInterface(bulk_silica=AmorphousSilicaBulk(),
                                        tile_x=tile_x,
                                        tile_y=tile_y,
                                        thickness=thickness,
                                        seed=seed)
        atom_names1 = np.array([atom.name for atom in interface1.particles()])

        interface2 = mb.recipes.SilicaInterface(bulk_silica=AmorphousSilicaBulk(),
                                        tile_x=tile_x,
                                        tile_y=tile_y,
                                        thickness=thickness,
                                        seed=seed)
        atom_names2 = np.array([atom.name for atom in interface2.particles()])

        assert np.array_equal(atom_names1, atom_names2)
        assert np.array_equal(interface1.xyz, interface2.xyz)
