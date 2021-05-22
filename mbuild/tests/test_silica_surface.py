import numpy as np
import pytest

import mbuild as mb
from mbuild.lib.surfaces import AmorphousSilicaSurface
from mbuild.tests.base_test import BaseTest


class TestAmorphousSilicaSurface(BaseTest):
    def test_create_amorphous_silica_surface(self):
        surface = AmorphousSilicaSurface()

        assert isinstance(surface, mb.Compound)
        assert surface.n_particles == 1800
        assert surface.n_bonds == 0
        assert np.array_equal(surface.periodicity, (True, True, False))
        assert np.array_equal(surface.box.lengths, [5.4366, 4.7082, 1.0])

    def test_amorphous_silica_surface_error(self):
        with pytest.raises(ValueError):
            surface = AmorphousSilicaSurface(surface_roughness=2)
