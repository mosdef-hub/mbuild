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

        self['port_2'].spin(2/3*np.pi, [1, 0, 0])
        self['port_3'].spin(-2/3*np.pi, [1, 0, 0])
        self['port_1'].spin(1/3*np.pi, [0, 0, 1])
        self['port_2'].spin(-1/3*np.pi, [0, 0, 1])
        self['port_3'].spin(-1/3*np.pi, [0, 0, 1])

        self['port_0'].translate([0, 0.073, 0])
        self['port_1'].translate([0.0594, -0.0243, 0])
        self['port_2'].translate([-0.042, -0.0243, 0.042])
        self['port_3'].translate([-0.042, -0.0243, -0.042])


if __name__ == '__main__':
    m = N4()
    m.save('n4.mol2', overwrite=True)
