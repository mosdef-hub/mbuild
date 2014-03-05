__author__ = 'sallai'

from mbuild.compound import *
from sandbox.twodimmesh import *
from nalkanetile import *

class AlkaneMesh(Compound):
    @classmethod
    def create(cls, n, m, ctx={}):

        m = TwoDimMesh.create(NAlkaneTile.create, n, m,
                              left_port_name='left_male_port',
                              right_port_name='right_female_port',
                              top_port_name='top_male_port',
                              bottom_port_name='bottom_female_port',
                              ctx=ctx)

        return m


if __name__ == "__main__":
    # m = AlkaneMesh.create(6,8,9,{"alkane_body.port_rotation": 5*pi/6})
    m = AlkaneMesh.create(6,8,ctx={"min_alkane_length": 10, "max_alkane_length": 15})

    print [(label,atom.pos) for label, atom in self.atoms()]
    self.plot(labels=False)