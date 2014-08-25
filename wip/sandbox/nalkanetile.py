__author__ = 'sallai'
from mbuild.port import *
from n_alkane import *

class NAlkaneTile(Compound):
    @classmethod
    def create(cls, ctx={}, label=None):

        m = super(NAlkaneTile, cls).create(label)

        length = random.randint(ctx["min_alkane_length"], ctx["max_alkane_length"])

        nalkane = NAlkane.create(length, ctx=ctx)

        bbmin, bbmax = nalkane.boundingbox()

        m.add(nalkane, 'alkane')

        space = 10

        m.add(Port.create(),'top_male_port')
        m.top_male_port.transform(Translation((0,bbmax[1]+space,0)))

        m.add(Port.create(),'left_male_port')
        m.left_male_port.transform(RotationAroundZ(pi/2))
        m.left_male_port.transform(Translation((bbmin[0]-space,0,0)))

        m.add(Port.create(),'bottom_female_port')
        m.bottom_female_port.transform(Translation((0,bbmin[1]-space,0)))

        m.add(Port.create(),'right_female_port')
        m.right_female_port.transform(RotationAroundZ(pi/2))
        m.right_female_port.transform(Translation((bbmax[0]+space,0,0)))

        return m

if __name__ == "__main__":
    t = NAlkaneTile.create(ctx={"min_alkane_length": 10, "max_alkane_length": 15, "alkane_label": "10-alkane"})
    print t
    t.plot()