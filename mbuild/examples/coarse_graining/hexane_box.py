from mbuild.examples.alkane.alkane import Alkane
from mbuild.lib.moieties.ch2 import CH2
from mbuild.lib.moieties.ch3 import CH3
import mbuild as mb


class Propane(mb.Compound):
    def __init__(self):
        super(Propane, self).__init__()

        c = Alkane(n=3, cap_front=True, cap_end=False)
        self.add(c, 'propane')

        self.add(c['down'], 'down', containment=False)


class Hexane(mb.Compound):
    def __init__(self):
        super(Hexane, self).__init__()

        self.add(Propane(), 'propane1')
        self.add(Propane(), 'propane2')

        mb.equivalence_transform(
                self['propane1'], self['propane1']['down'], self['propane2']['down'])


if __name__ == '__main__':
    hexane_box = mb.fill_box(Hexane(), 50, box=[3, 3, 3])

    united_atom_particles = [CH2, CH3]
    three_to_one_particles = [Propane]

    united_atom = mb.coarse_grain(hexane_box, particle_classes=united_atom_particles)
    three_to_one = mb.coarse_grain(hexane_box, particle_classes=three_to_one_particles)

    # print("Proxy: {}".format(proxy))
    #
    # print("Leaves of the proxy:")
    # for i, leaf in enumerate(proxy.particles()):
    #     print("{}: {}".format(i, leaf))
