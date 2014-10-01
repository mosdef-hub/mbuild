import numpy as np

from mbuild.coordinate_transform import translate, rotate_around_x
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
                rotate_around_x(port, np.pi/2)
                translate(port, atom + np.array([0, 0, .1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = Surface()
    m = TiledCompound(s, n_x=3, n_y=2, n_z=1, kind="tiled")
    m.to_trajectory()

    m.save(filename='tiled_surface.xyz')

    # print s.periodicity
    # print m.periodicity

    #from mbuild.plot import Plot
    #Plot(m, bonds=True, verbose=True, periodic_bonds=True).show()
