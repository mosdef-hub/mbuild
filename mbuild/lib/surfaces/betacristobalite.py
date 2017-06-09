import numpy as np

import mbuild as mb


class Betacristobalite(mb.Compound):
    """The beta-cristobalite form of SiO2.

    Area per port specifies the density of attachment sites in nm^2.
    The crystal is expanded to yield an area per port of 0.25 nm^2, the
    typical density of alkane monolayers on SiO2 although these are actually
    grown on amorphous SiO2 in experiment.

    See http://www.wikiwand.com/en/Silicon_dioxide for more info on the various
    crystal forms.

    Note: Port sites are currently naively determined by placing them on all
    oxygens which are above 1.0 nm in the z-direction. This only holds true for
    the beta-cristobalite-expanded.mol2 file. If you add a new one, please modify
    the file or the method of determining port locations.

    """
    def __init__(self):
        super(Betacristobalite, self).__init__()

        mb.load('beta-cristobalite-expanded.mol2', compound=self,
                relative_to_module=self.__module__)
        self.periodicity = np.array([5.3888, 4.6669, 0.0])

        count = 0
        for particle in self.particles():
            if particle.name.startswith('O') and particle.pos[2] > 1.0:
                count += 1
                port = mb.Port(anchor=particle, orientation=[0, 0, 1], 
                               separation=0.1)
                self.add(port, 'port_{}'.format(count))
                particle.name = 'O'  # Strip numbers required in .mol2 files.
            elif particle.name.startswith('Si'):
                particle.name = 'Si'

if __name__ == "__main__":
    single = Betacristobalite()
    multiple = mb.TiledCompound(single, n_tiles=(2, 1, 1), name="tiled")
    multiple.save('betacristobalite.mol2')
