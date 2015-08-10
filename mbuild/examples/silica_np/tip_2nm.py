import numpy as np

from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_x, translate
from mbuild.port import Port
from mbuild.tools.tiled_compound import TiledCompound
from mbuild.coordinate_transform import translate_to, _extract_atom_positions, _write_back_atom_positions


class Tip2nm(Compound):
    """ """
    def __init__(self, radius=2.0):
        super(Tip2nm, self).__init__()

        if radius == 2.0:
            self.append_from_file('tip_2nm.pdb',
                                  relative_to_module=self.__module__)
            #self.periodicity = np.array([5.4366, 4.7082, 0.0])
        else:
            raise ValueError('Curved tip input file with radius '
                             'of {0:.1f} does not exist. If you have '
                             'this structure, please submit a pull request to'
                             'add it! '.format(radius))

        cnt = 0
        atom_positions = _extract_atom_positions(self)
        x = atom_positions[:,0]
        y = atom_positions[:,1]
        z = atom_positions[:,2]
        minx = min(x)
        maxx = max(x)
        midx = (maxx+minx)/2.0
        miny = min(y)
        maxy = max(y)
        midy = (maxy+miny)/2.0
        maxz = max(z)
        for i in range(len(atom_positions)):
            atom_positions[i][0] -= midx
            atom_positions[i][1] -= midy
            atom_positions[i][2] -= (maxz-0.4)
        _write_back_atom_positions(self,atom_positions)        

        #print self.atoms

        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            if atom.kind == 'OB':
                cnt += 1
                port = Port(anchor=atom)
                rotate_around_x(port, -np.pi/2)
                translate(port, atom + np.array([0, 0, -.1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = Tip2nm()
    m = TiledCompound(s, n_x=1, n_y=1, n_z=1, kind="tiled")
    m.visualize(show_ports=True)
