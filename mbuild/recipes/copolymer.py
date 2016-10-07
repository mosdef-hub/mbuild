from mbuild.compound import Compound
from mbuild.coordinate_transform import equivalence_transform
from mbuild.utils.validation import assert_port_exists
from mbuild import clone
from mbuild.recipes import polymer

from mbuild.lib.moieties import CH2, CH3
from mbuild.lib.moieties import Silane

import itertools

__all__ = ['Copolymer']


class Copolymer(Compound):
    """Connect a component to successive copies of itself.

    Parameters
    ----------
    monomers : a set of mb.Compound
        The set of compounds to replicate by the pattern
    pattern : the pattern to form a copolymer such as "AABB", "ABAB"
    port_labels : 2-tuple of strs
        The names of the two ports to use to connect copies of proto.
    n : int, optional, default=2
        The number of times to replicate proto.

    """

    def __init__(self, pattern=None, random=False):
        # parse pattern here
        # proto_dict = {}
        # pattern_list = list(pattern)
        # proto_dict[pattern_list[0]] = monomers[0]
        # proto_dict[pattern_list[2]] = monomers[1]

        """
        :param pattern: General format for all kinds of copolymers.
            [cap_front, [A, A_num, A_port_labels], [B, B_num, B_port_labels], cap_end, repeat_num]
            Alternating copolymers: A_num = 1, B_num = 1, repeat_num = any_value
            Periodic copolymers: A_num = any_value, B_num = any_value, repeat_num = any_value
            Block copolymers: A_num = any_value, B_num = any_value, repeat_num = 1
            Statistical copolymers: repeat A block and B block randomly with repeat_num times, if random is True
        :param random: option for Statistical copolymers, which are actually random copolymers.
        """

        # default value for pattern
        if pattern is None:
            ch2 = CH2()
            silane = Silane()
            pattern = [None, [ch2, 2, ("up", "down")], [silane, 2, ("up", "down")], 2, None]

        # parse pattern into real chain of copolymer
        copolymer = []
        port_labels = ("up", "down")  # we need to modify it based on the pattern
        # 9/28/2016 work here, default port_labels only for test
        for proto in pattern[1:-2]:
            copolymer.append([proto[0]] * proto[1])  # TypeError: unsupported operand type(s) for *: 'CH2' and 'int'
            for label in proto[2]:
                assert_port_exists(label, proto[0])  # assert labels in the proto

        copolymer = list(itertools.chain.from_iterable(copolymer * pattern[-2]))  # a complete chain for copolymer

        if False:  # Exception of parse pattern here
            raise Exception('n must be 1 or more')
        super(Copolymer, self).__init__()

        first_part = None
        last_part = None
        for molecule in copolymer:
            this_part = clone(molecule)
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that it's bottom port is rotated
                # and translated to the last part's top port.
                equivalence_transform(this_part, this_part.labels[port_labels[1]],
                                      last_part.labels[port_labels[0]])
            last_part = this_part

        # if pattern[0]:  # think about the way in alkane later...
        #     this_part = pattern[0]
        #     self.add(this_part, 'monomer[$]')
        #     equivalence_transform(this_part, this_part.labels[port_labels[1]],
        #                          last_part.labels[port_labels[0]])
        # else:
        #     # Hoist the first part's bottom port to be the bottom port of the copolymer.
        #     self.add(first_part.labels[port_labels[1]], port_labels[1], containment=False)
        #
        # # Hoist the last part's top port to be the top port of the copolymer.
        # self.add(last_part.labels[port_labels[0]], port_labels[0], containment=False)

        # chain = mb.Polymer(CH2(), n=n - 2, port_labels=('up', 'down'))
        # self.add(chain, 'chain')
        #
        # if pattern[0]:
        #     self.add(CH3(), "methyl_front")
        #     equivalence_transform(self.monomer, self.monomer.up, self.methyl_front.up)
        # else:
        #     # Hoist port label to Alkane level.
        #     self.add(monomer.up, 'up', containment=False)
        #
        # if pattern[-1]:
        #     self.add(CH3(), 'methyl_end')
        #     equivalence_transform(self.methyl_end, self.methyl_end.up, self.chain.down)
        # else:
        #     # Hoist port label to Alkane level.
        #     self.add(chain.down, 'down', containment=False)



if __name__ == "__main__":

    ch3 = CH3()
    ch2 = CH2()
    silane = Silane()
    # poly = Copolymer(monomers=(ch2, silane), pattern="A3B3", port_labels=("up", "down"))
    # combine monomers and pattern together
    pattern = [ch3, [ch2, 2, ("up", "down")], [ch2, 2, ("up", "down")], 2, None]
    poly = Copolymer(pattern=pattern, random=False)
    poly.visualize(show_ports=True)
