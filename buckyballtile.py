__author__ = 'sallai'
from mbuild.compound import *
from mbuild.port import *
from mbuild.xyz import *
class BuckyBallTile(Compound):
    @classmethod
    def create(cls, label=None):
        """

        :param label:
        :return:
        """
        m = super(BuckyBallTile, cls).create(label)

        buckyball =Xyz.create('mbuild/c60.xyz')

        bbmin, bbmax = buckyball.boundingbox()

        m.add(buckyball, 'c60')

        space = 10

        m.add(Port.create(),'top_male_port')
        m.top_male_port.applyTransformation(CoordinateTransform.translation((0,bbmax[1]+space,0)))

        m.add(Port.create(),'left_male_port')
        m.left_male_port.applyTransformation(CoordinateTransform.rotation_around_z(pi/2))
        m.left_male_port.applyTransformation(CoordinateTransform.translation((bbmin[0]-space,0,0)))

        m.add(Port.create(),'bottom_female_port')
        m.bottom_female_port.applyTransformation(CoordinateTransform.translation((0,bbmin[1]-space,0)))

        m.add(Port.create(),'right_female_port')
        m.right_female_port.applyTransformation(CoordinateTransform.rotation_around_z(pi/2))
        m.right_female_port.applyTransformation(CoordinateTransform.translation((bbmax[0]+space,0,0)))

        return m

if __name__ == "__main__":
    t = BuckyBallTile.create()
    print t