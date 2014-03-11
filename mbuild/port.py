__author__ = 'sallai'
from .compound import *


class Port(Compound):
    def __init__(self, ctx={}):
        super(Port, self).__init__(kind='Port', ctx=ctx)
        self.add(Atom(kind='G', pos=(0, 0, 0)), 'middle')
        self.add(Atom(kind='G', pos=(0, 0.2, 0)), 'top')
        self.add(Atom(kind='G', pos=(-0.2, -0.1, 0)), 'left')
        self.add(Atom(kind='G', pos=(0.0, -0.2, 0.1)), 'right')
