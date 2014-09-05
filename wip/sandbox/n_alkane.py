__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.examples.methane.methane import *

class NAlkane(Compound):

    def __init__(self, n, ctx={}):
        super(NAlkane, self).__init__(kind='NAlkane', ctx=ctx)

        if n < 2:
            raise Exception('n must be 2 or more')

        # top tail (CH_3)
        self.add(AlkaneTail(ctx=ctx),'top_tail')

        # n times the body CH_2
        last_part = self.top_tail
        for body_count in range(1, n-1):
            this_part = AlkaneBody(ctx=ctx)
            this_part.transform([(this_part.female_port, last_part.male_port)])
            # self.add(this_part, 'body_'+str(body_count))
            self.add(this_part, 'body_#') # this will generate body_0, body_1,...
            last_part = this_part

        # bottom tail (CH_3)
        self.add(AlkaneTail(ctx=ctx),'bottom_tail')
        self.bottom_tail.transform( [(self.bottom_tail.female_port, last_part.male_port)])

    @classmethod
    def ethane(cls):
        return NAlkane(2)

    @classmethod
    def propane(cls):
        return NAlkane(3)

    @classmethod
    def buthane(cls):
        return NAlkane(4)

    @classmethod
    def pentane(cls):
        return NAlkane(5)

    @classmethod
    def hexane(cls):
        return NAlkane(6)


if __name__ == "__main__":
    # m = NAlkane(8)
    m = NAlkane.pentane()
    TreeView(m).show()
