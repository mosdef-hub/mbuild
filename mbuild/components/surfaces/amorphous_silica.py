import numpy as np

from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_x, translate
from mbuild.port import Port
from mbuild.tools.tiled_compound import TiledCompound


class AmorphousSilica(Compound):
    """ """
    def __init__(self, surface_roughness=1.0):
        super(AmorphousSilica, self).__init__()

        if surface_roughness == 1.0:
            # TODO: description of how this surface was generated/citation
            self.append_from_file('amorphous_silica_sr1.0.pdb')
            self.periodicity = np.array([5.4366, 4.7082, 0.0])
        else:
            raise ValueError('Amorphous silica input file with surface '
                             'roughness of {0:.1f} does not exist. If you have '
                             'this structure, please submit a pull request to'
                             'add it! '.format(surface_roughness))

        cnt = 0
        for atom in self.atoms:
            if atom.kind == 'OB':
                cnt += 1
                port = Port(anchor=atom)
                rotate_around_x(port, np.pi/2)
                translate(port, atom + np.array([0, 0, .1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = AmorphousSilica()
    m = TiledCompound(s, n_x=2, n_y=1, n_z=1, kind="tiled")
    m.visualize(show_ports=True)
