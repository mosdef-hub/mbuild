import numpy as np

import mbuild as mb


class Betacristobalite(mb.Compound):
    """The beta-cristobalite form of SiO2.

    Area per port specifies the density of attachment sites in nm^2.
    * 0.2 nm^2 / port is the most physically accurate.
    * 0.25 nm^2 / port is the typical density of alkane monolayers on SiO2
      although these are actually grown on amorphous SiO2 in experiment.
      This particular file is the same as the 0.2 file with its dimensions
      stretched in the x- and y-dimensions.

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
            mb.load('beta-cristobalite.pdb', compound=self,
                    relative_to_module=self.__module__)
            self.periodicity = np.array([4.7689, 4.13, 0.0])
        elif area_per_port == 0.25:

            mb.load('beta-cristobalite-expanded.mol2', compound=self,
                    relative_to_module=self.__module__)
            self.periodicity = np.array([5.3888, 4.6669, 0.0])

        count = 0
        for atom in self.atoms:
            if atom.name == 'O' and atom.pos[2] > 1.0:
                count += 1
                port = mb.Port(anchor=atom)
                mb.rotate_around_x(port, np.pi/2)
                mb.translate(port, atom + np.array([0, 0, .1]))
                self.add(port, 'port_{}'.format(count))

if __name__ == "__main__":
    single = Betacristobalite()
    multiple = mb.TiledCompound(single, n_tiles=(2, 1, 1), kind="tiled")
    multiple.visualize(show_ports=True)
