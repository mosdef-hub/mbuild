import numpy as np

import mbuild as mb


class Graphene(mb.Compound):
    """Repeat  single layer of graphene"""
    def __init__(self):
        super(Graphene, self).__init__()

        mb.load('graphene.pdb', compound=self,
                relative_to_module=self.__module__)
        self.periodicity = np.array([1.23, 0.852, 0.0])

        for count, particle in enumerate(self.particles()):
            port = mb.Port(anchor=particle)
            mb.spin_x(port, np.pi/2)
            mb.translate(port, np.array([0, 0, 0.1]))
            self.add(port, 'port_{}'.format(count))

if __name__ == "__main__":
    graphene = Graphene()
    graphene.visualize(show_ports=True)
