from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.testing.tools import get_fn
from mbuild.tools.tiled_compound import TiledCompound


class Surface(Compound):
    """ """
    def __init__(self):
        super(Surface, self).__init__()

        self.append_from_file(get_fn('beta-cristobalite.pdb'))
        self.periodicity = np.array([4.7689, 4.13, 0.0])

        cnt = 0
        for atom in self.atoms():
            if atom.kind == 'O' and atom.pos[2] > 1:
                cnt += 1
                port = Port(anchor=atom)
                rotate_around_x(port, pi/2)
                translate(port, atom + np.array([0, 0, .1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = Surface()
    m = TiledCompound(s, n_x=2, n_y=1, n_z=1, kind="tiled")

    m.visualize()
