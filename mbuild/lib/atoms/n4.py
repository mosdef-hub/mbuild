import mbuild as mb
import numpy as np


class N4(mb.Compound):
    """An tetravalent nitrogen atom."""
    def __init__(self):
        super(N4, self).__init__()

        self.add(mb.Particle(name='N', pos=[0, 0, 0]), label='N[$]')
        self.add(mb.Port(anchor=self[0]), label='port_0')
        self.add(mb.Port(anchor=self[0]), label='port_1')
        self.add(mb.Port(anchor=self[0]), label='port_2')
        self.add(mb.Port(anchor=self[0]), label='port_3')

        mb.spin_x(self['port_2'], 2/3*np.pi)
        mb.spin_x(self['port_3'], -2/3*np.pi)
        mb.spin_z(self['port_1'], 1/3*np.pi)
        mb.spin_z(self['port_2'], -1/3*np.pi)
        mb.spin_z(self['port_3'], -1/3*np.pi)

        mb.translate(self['port_0'], [0, 0.073, 0])
        mb.translate(self['port_1'], [0.0594, -0.0243, 0])
        mb.translate(self['port_2'], [-0.042, -0.0243, 0.042])
        mb.translate(self['port_3'], [-0.042, -0.0243, -0.042])


if __name__ == '__main__':
    m = N4()
    m.visualize(show_ports=True)
