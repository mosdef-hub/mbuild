import random
from copy import deepcopy

import numpy as np

from mbuild.examples.tnp.tnp import Tnp
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x, rotate_around_y, rotate_around_z
from mbuild.tools.mask import grid_mask_3d

__author__ = 'sallai'

class TnpBox(Compound):
    def __init__(self):
        super(TnpBox, self).__init__(self)

        tnp_proto = Tnp(ball_radius=5, n_chains=5, chain_length=8)

        mask = grid_mask_3d(3, 3, 3) * 100

        rnd = random.Random()
        rnd.seed(1928)

        for pos in mask:
            tnp = deepcopy(tnp_proto)
            rotate_around_x(tnp, rnd.uniform(0, 2*np.pi))
            rotate_around_y(tnp, rnd.uniform(0, 2*np.pi))
            rotate_around_z(tnp, rnd.uniform(0, 2*np.pi))
            translate(tnp, pos)
            self.add(tnp)

if __name__ == "__main__":
    m = TnpBox()

