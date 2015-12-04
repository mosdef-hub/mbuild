import mbuild as mb
from mbuild.proxy import create_proxy

from mbuild.examples.alkane.alkane import Alkane


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
                self.propane1, self.propane1['down'], self.propane2['down'])


if __name__ == '__main__':

    p = Hexane()

    tier = [Propane]

    proxy = create_proxy(p, particle_classes=tier)

    print("Proxy: {}".format(proxy))

    print("Leaves of the proxy:")
    for i, leaf in enumerate(proxy.particles):
        print("{}: {}".format(i, leaf))



    proxy.visualize()