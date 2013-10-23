__author__ = 'sallai'

from alkane_tail import AlkaneTail
from mbuild.compound import *


class Ethane(Compound):

    @classmethod
    def create(cls, label=None):
        m = super(Ethane, cls).create(label)
        # two tails
        m.add(AlkaneTail.create(), 'top_tail')
        m.add(AlkaneTail.create(), 'bottom_tail')
        # transform bottom_tail to top_tail's coordinate system with point equivalencies
        # m.bottom_tail.transform(
        #                     [
        #                             (m.bottom_tail.female_port.top, m.top_tail.male_port.top),
        #                             (m.bottom_tail.female_port.left, m.top_tail.male_port.left),
        #                             (m.bottom_tail.female_port.right, m.top_tail.male_port.right)
        #                     ])

        # transform bottom_tail to top_tail's coordinate system with subcomponent equivalencies
        m.bottom_tail.transform([(m.bottom_tail.female_port, m.top_tail.male_port)])

        return m


if __name__ == "__main__":
    ethane = Ethane.create()
    # print ethane
    print ethane.atoms()
    ethane.plot(labels=False, verbose=False)
    print ethane.bottom_tail.c