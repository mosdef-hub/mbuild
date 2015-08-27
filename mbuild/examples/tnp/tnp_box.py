import random
import mbuild as mb

from numpy import pi
from mbuild.examples.tnp.tnp import Tnp
from mbuild import clone

class TnpBox(mb.Compound):
    """Several tethered nanoparticles randomly dispersed in a Box. """
    def __init__(self):
        super(TnpBox, self).__init__()

        tnp_proto = Tnp(ball_radius=5, n_chains=5, chain_length=8)

        mask = mb.grid_mask_3d(3, 3, 3) * 100

        rnd = random.Random()
        rnd.seed(1928)

        for pos in mask:
            tnp = clone(tnp_proto)
            mb.rotate_around_x(tnp, rnd.uniform(0, 2 * pi))
            mb.rotate_around_y(tnp, rnd.uniform(0, 2 * pi))
            mb.rotate_around_z(tnp, rnd.uniform(0, 2 * pi))
            mb.translate(tnp, pos)
            self.add(tnp)

def main():
    box = TnpBox()
    return box

if __name__ == "__main__":
    box = main()
    box.visualize()
