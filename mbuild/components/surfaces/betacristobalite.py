import numpy as np

from mbuild.coordinate_transform import rotate_around_x, translate
from mbuild.compound import Compound
from mbuild.port import Port
from mbuild.tools.tiled_compound import TiledCompound


class Betacristobalite(Compound):
    """The beta-cristobalite form of SiO2.

    Area per port specifies the density of attachment sites in nm^2.
        -0.2 nm^2 / port is the most physically accurate.
        -0.25 nm^2 / port is the typical density of alkane monolayers on
         SiO2 although these are actually grown on amorphous SiO2 in
         experiment.
         This particular file is the same as the 0.2 file with its
         dimensions stretched in the x- and y-dimensions.

    See http://www.wikiwand.com/en/Silicon_dioxide for more info on the various
    crystal forms.

    Note: Port sites are currently naively determined by placing them on all
    oxygens which are above 1.0 nm in the z-direction. This only holds true for
    the beta-cristobalite.pdb and beta-cristobalite-expanded.mol2 files. If you
    add a new one, please modify the file or the method of determining port
    locations.

    """
    def __init__(self, area_per_port=0.2):
        super(Betacristobalite, self).__init__()

        if area_per_port == 0.2:
            self.append_from_file('beta-cristobalite.pdb')
            self.periodicity = np.array([4.7689, 4.13, 0.0])
        elif area_per_port == 0.25:

            self.append_from_file('beta-cristobalite-expanded.mol2')
            self.periodicity = np.array([5.3888, 4.6669, 0.0])

        cnt = 0
        for atom in self.atoms:
            if atom.kind == 'O' and atom.pos[2] > 1.0:
                cnt += 1
                port = Port(anchor=atom)
                rotate_around_x(port, np.pi/2)
                translate(port, atom + np.array([0, 0, .1]))
                self.add(port, 'port_{}'.format(cnt))

if __name__ == "__main__":
    s = Betacristobalite()
    m = TiledCompound(s, n_x=2, n_y=1, n_z=1, kind="tiled")
    m.visualize(show_ports=True)
