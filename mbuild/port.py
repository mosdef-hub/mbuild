__author__ = 'sallai'
from mbuild.compound import *


class Port(Compound):
    @classmethod
    def create(cls, label=None):
        m = super(Port, cls).create(label)
        m.add(G((0, 0.02, 0)),'top')
        m.add(G((-0.02, -0.02, 0)),'left')
        m.add(G((0.02, -0.02, 0)),'right')
        return m