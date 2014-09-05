__author__ = 'sallai'
from mbuild.port import *


class Box(Compound):
    @classmethod
    def create(cls, creator, ctx={}):
        m = super(Box, cls).create(ctx=ctx)

        # create the compound
        contents = creator(ctx=ctx)
        m.add(contents, contents.kind)

        m.ctx = ctx

        m.addBoxPorts()

        # bbmin, bbmax = contents.boundingbox()
        #
        # space = 3
        #
        # m.add(Port.create(), 'top_male_port')
        # m.top_male_port.transform(Translation((0, bbmax[1] + space, 0)))
        #
        # m.add(Port.create(), 'left_male_port')
        # m.left_male_port.transform(RotationAroundZ(pi / 2))
        # m.left_male_port.transform(Translation((bbmin[0] - space, 0, 0)))
        #
        # m.add(Port.create(), 'up_male_port')
        # m.up_male_port.transform(RotationAroundX(pi / 2))
        # m.up_male_port.transform(Translation((0, 0, bbmax[2] + space)))
        #
        # m.add(Port.create(), 'bottom_female_port')
        # m.bottom_female_port.transform(Translation((0, bbmin[1] - space, 0)))
        #
        # m.add(Port.create(), 'right_female_port')
        # m.right_female_port.transform(RotationAroundZ(pi / 2))
        # m.right_female_port.transform(Translation((bbmax[0] + space, 0, 0)))
        #
        # m.add(Port.create(), 'down_female_port')
        # m.down_female_port.transform(RotationAroundX(pi / 2))
        # m.down_female_port.transform(Translation((0, 0, bbmin[0] - space)))

        return m

    def addBoxPorts(self):
        space = 3

        bbmin, bbmax = self.boundingbox()

        self.add(Port.create(), 'top_male_port')
        self.top_male_port.transform(Translation((0, bbmax[1] + space, 0)))

        self.add(Port.create(), 'left_male_port')
        self.left_male_port.transform(RotationAroundZ(pi / 2))
        self.left_male_port.transform(Translation((bbmin[0] - space, 0, 0)))

        self.add(Port.create(), 'up_male_port')
        self.up_male_port.transform(RotationAroundX(pi / 2))
        self.up_male_port.transform(Translation((0, 0, bbmax[2] + space)))

        self.add(Port.create(), 'bottom_female_port')
        self.bottom_female_port.transform(Translation((0, bbmin[1] - space, 0)))

        self.add(Port.create(), 'right_female_port')
        self.right_female_port.transform(RotationAroundZ(pi / 2))
        self.right_female_port.transform(Translation((bbmax[0] + space, 0, 0)))

        self.add(Port.create(), 'down_female_port')
        self.down_female_port.transform(RotationAroundX(pi / 2))
        self.down_female_port.transform(Translation((0, 0, bbmin[0] - space)))

if __name__ == "__main__":
    from mbuild.examples.methane.methane import *
    def creator(ctx={}):
        return Methane.create(ctx)

    t = Box.create(creator,ctx={})

    print([a for a in t.atoms()])

    t.plot(verbose=True)
    print t