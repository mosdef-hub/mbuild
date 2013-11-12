__author__ = 'sallai'
from mbuild.compound import *


class Port(Compound):
    @classmethod
    def create(cls, kind=None):
        m = super(Port, cls).create(kind="port")
        m.add(G((0, 0.2, 0)), 'top')
        m.add(G((-0.2, -0.2, 0)), 'left')
        m.add(G((0.2, -0.2, 0)), 'right')
        return m