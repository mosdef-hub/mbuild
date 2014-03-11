__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.compound import *
from mpcalkanebody import *

class MpcAlkane(Compound):

    def __init__(self, n, ctx={}):
        super(MpcAlkane, self).__init__(kind='MpcAlkane', ctx=ctx)

        if n < 2:
            raise Exception('n must be 2 or more')

        # bot tail (CH_3)
        self.add(AlkaneTail(ctx=ctx),'bot_ch3')
        self.bot_ch3.rename([('C','CB'),('H','HB')])

        # n times the body CH_2
        last_part = self.bot_ch3

        for body_count in range(2, n):
            if body_count % 20 == 0:
                this_part = MpcAlkaneBody(direction='left', ctx=ctx)
            elif body_count % 10 == 0:
                this_part = MpcAlkaneBody(direction='right', ctx=ctx)
            elif body_count % 2 == 0:
                this_part = AlkaneBody(ctx=ctx, direction='left')
                this_part.rename([('C','CB'),('H','HB')])
            else:
                this_part = AlkaneBody(ctx=ctx, direction='right')
                this_part.rename([('C','CB'),('H','HB')])

            this_part.transform([(this_part.male_port, last_part.female_port)])
            self.add(this_part, 'body_'+str(body_count))
            last_part = this_part

        # top tail (CH_3)
        self.add(AlkaneTail(ctx=ctx),'top_ch3')
        self.top_ch3.rename([('C','CB'),('H','HB')])

        self.top_ch3.transform( [(self.top_ch3.male_port, last_part.female_port)])


if __name__ == "__main__":
    m = MpcAlkane(21)
    #print m.atoms()
    #m.plot(labels=False, verbose=False)
    Xyz.save(m, 'mpcchain.xyz')
    TreeView(m).show()
