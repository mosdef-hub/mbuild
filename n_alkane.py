__author__ = 'sallai'

from alkane_tail import AlkaneTail
from alkane_body import AlkaneBody
from mbuild.compound import *
from methane import *

class NAlkane(Compound):
    @classmethod
    def create(cls, n, label=None):

        if n < 2:
            raise Exception('n must be 2 or more')

        m = super(NAlkane, cls).create(label)

        # top tail (CH_3)
        m.add(AlkaneTail.create(),'top_tail')

        # n times the body CH_2
        last_part = m.top_tail
        for body_count in range(1, n-1):
            this_part = AlkaneBody.create()
            this_part.transform([(this_part.female_port, last_part.male_port)])
            m.add(this_part, 'body_'+str(body_count))
            last_part = this_part

        # bottom tail (CH_3)
        m.add(AlkaneTail.create(),'bottom_tail')
        m.bottom_tail.transform( [(m.bottom_tail.female_port, last_part.male_port)])

        return m

    @classmethod
    def ethane(cls):
        return NAlkane.create(2)

    @classmethod
    def propane(cls):
        return NAlkane.create(3)

    @classmethod
    def buthane(cls):
        return NAlkane.create(4)

    @classmethod
    def pentane(cls):
        return NAlkane.create(5)

    @classmethod
    def hexane(cls):
        return NAlkane.create(6)


if __name__ == "__main__":
    m = NAlkane.create(8)
    # m = NAlkane.pentane()
    print [(label,atom.pos) for label, atom in m.atoms()]
    m.plot(labels=False)