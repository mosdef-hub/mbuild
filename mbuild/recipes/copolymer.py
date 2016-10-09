from mbuild.compound import Compound
from mbuild.coordinate_transform import equivalence_transform
from mbuild.utils.validation import assert_port_exists
from mbuild import clone

from mbuild.lib.moieties import CH2, CH3
from mbuild.lib.moieties import Silane

import itertools
from numpy.random import permutation

__all__ = ['Copolymer']


class Copolymer(Compound):
    # Use MainChain() get the main_chain of Copolymer, and then add cap_front and cap_end to it.

    def __init__(self, pattern=None, random=False):
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
        super().__init__()
        if pattern is None:
            ch2 = CH2()
            silane = Silane()
            pattern = [None, [ch2, 2, ("up", "down")], [silane, 2, ("up", "down")], 2, None]

        # get main_chain
        main_chain = MainChain(pattern, random)
        self.add(main_chain, 'main_chain')

        # add caps
        if pattern[0]:
            self.add(pattern[0], "cap_front")
            equivalence_transform(self['main_chain'], self['main_chain']['up'], self['cap_front']['up'])
        else:
            # Hoist port label to Alkane level.
            self.add(main_chain['up'], 'up', containment=False)

        if pattern[-1]:
            self.add(pattern[-1], 'cap_end')
            equivalence_transform(self['cap_end'], self['cap_end']['up'], self['main_chain']['down'])
        else:
            # Hoist port label to Alkane level.
            self.add(main_chain['down'], 'down', containment=False)


class MainChain(Compound):
    # The main chain for Copolymer. Connect one or more components in a specified sequence.

    def __init__(self, pattern=None, random=False):
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
        super().__init__()
        if pattern is None:
            ch2 = CH2()
            silane = Silane()
            pattern = [None, [ch2, 2, ("up", "down")], [silane, 2, ("up", "down")], 2, None]

        # parse pattern into real chain of copolymer
        copolymer = []
        type_A = type(pattern[1][0])
        type_B = type(pattern[2][0])
        port_labels_A = pattern[1][2]
        port_labels_B = pattern[2][2]

        for proto in pattern[1:-2]:
            if proto[1] < 1:  # Exception of parse pattern here
                raise Exception('n must be 1 or more')
            copolymer.append([proto[0]] * proto[1])  # append multiple times
            for label in proto[2]:
                assert_port_exists(label, proto[0])  # assert labels in the proto

        copolymer = list(itertools.chain.from_iterable(copolymer * pattern[-2]))  # a complete chain for copolymer
        if random:
            copolymer = list(permutation(copolymer))

        first_part = None
        last_part = None
        for molecule in copolymer:
            this_part = clone(molecule)
            # print(isinstance(this_part, type(pattern[2][0])))
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that it's bottom port is rotated
                # and translated to the last part's top port.
                if isinstance(this_part, type_A):
                    equivalence_transform(this_part,
                                          this_part.labels[port_labels_A[1]],
                                          last_part.labels[port_labels_A[0]])
                elif isinstance(this_part, type_B):
                    equivalence_transform(this_part,
                                          this_part.labels[port_labels_B[1]],
                                          last_part.labels[port_labels_B[0]])
            last_part = this_part

        # Hoist the first part's bottom port to be the bottom port of the polymer.
        if isinstance(first_part, type_A):
            self.add(first_part.labels[port_labels_A[1]], port_labels_A[1], containment=False)
        elif isinstance(first_part, type_B):
            self.add(first_part.labels[port_labels_B[1]], port_labels_B[1], containment=False)

        # Hoist the last part's top port to be the top port of the polymer.
        if isinstance(last_part, type_A):
            self.add(last_part.labels[port_labels_A[0]], port_labels_A[0], containment=False)
        elif isinstance(last_part, type_B):
            self.add(last_part.labels[port_labels_B[0]], port_labels_B[0], containment=False)


if __name__ == "__main__":

    cap_front = CH3()
    cap_end = CH3()

    chain_A = CH2()
    chain_B = Silane()

    pattern = [cap_front, [chain_A, 2, ("up", "down")], [chain_B, 2, ("up", "down")], 1, cap_end]
    poly = Copolymer(pattern=pattern, random=False)
    poly.visualize(show_ports=False)
