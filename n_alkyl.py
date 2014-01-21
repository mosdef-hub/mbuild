from mbuild.port import Port

__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.compound import *
from methane import *
import operator

class NAlkyl(Compound):
    @classmethod
    def create(cls, n, ctx={}):

        if n < 2:
            raise Exception('n must be 2 or more')

        m = super(NAlkyl, cls).create(ctx=ctx)

        # top tail (CH_3)
        m.add(AlkaneTail.create(ctx=ctx),'top_tail')

        # n-1 times the body CH_2
        last_part = m.top_tail
        direction='left'
        for body_count in range(1, n):
            this_part = AlkaneBody.create(ctx=ctx, direction=direction)
            if(direction=='left'):
                direction='right'
            else:
                direction='left'
            this_part.transform([(this_part.female_port, last_part.male_port)])
            m.add(this_part, 'body_'+str(body_count))
            last_part = this_part


        port = Port.create()
        translateTo = map(operator.sub,last_part.female_port.top.pos, port.top.pos)
        port.transform(Translation(tuple(translateTo)))


        m.add(port, "port")

        return m

if __name__ == "__main__":
    m = NAlkyl.create(8)
    print [(label,atom.pos) for label, atom in m.atoms()]
    m.plot(labels=False, verbose=True)