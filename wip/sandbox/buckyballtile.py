__author__ = 'sallai'
from mbuild.port import *
from mbuild.file_formats.xyz import *
class BuckyBallTile(Compound):
    @classmethod
    def create(cls, ctx={}):
        """

        :param label:
        :return:
        """
        m = super(BuckyBallTile, cls).create(ctx=ctx)

        buckyball = Xyz.create('mbuild/c60.xyz')

        bbmin, bbmax = buckyball.boundingbox()

        m.add(buckyball, 'c60')

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
    t = BuckyBallTile.create()
    print t
    t.plot()