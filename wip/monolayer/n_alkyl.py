from mbuild.plot import Plot
from mbuild.port import Port

__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.compound import *
import operator

class NAlkyl(Compound):

    def __init__(self, n, ctx={}):

        if n < 2:
            raise Exception('n must be 2 or more')

        super(NAlkyl, self).__init__(ctx=ctx)

        # top tail (CH_3)
        self.add(AlkaneTail(ctx=ctx),'top_tail')

        # n-1 times the body CH_2
        last_part = self.top_tail
        direction='left'
        for body_count in range(1, n):
            this_part = AlkaneBody(ctx=ctx, direction=direction)
            if(direction=='left'):
                direction='right'
            else:
                direction='left'
            this_part.transform([(this_part.female_port, last_part.male_port)])
            self.add(this_part, 'body_'+str(body_count))
            last_part = this_part


        port = Port()
        translateTo = map(operator.sub,last_part.male_port.top.pos, port.top.pos)
        port.transform(Translation(tuple(translateTo)))


        self.add(port, "port")

if __name__ == "__main__":
    m = NAlkyl(8)
    Plot(m, verbose=True).show()
