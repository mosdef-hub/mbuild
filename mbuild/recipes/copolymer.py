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
           [cap_front, [monomer_A, monomer_B, ...]*10, cap_end]
        :param random: option for random copolymers.
        """

        # default value for pattern
        super().__init__()
        if pattern is None:
            cap_front = CH3()
            cap_end = CH3()
            monomer_A = [CH2(), ("up", "down")] * 2
            monomer_B = [Silane(), "up", "down"] * 2
            pattern = [cap_front, [monomer_A, monomer_B] * 2, cap_end]

        # get main_chain
        main_chain = MainChain(pattern, random)
        self.add(main_chain, 'main_chain')
        # get caps
        cap_front = pattern[0]
        cap_end = pattern[-1]

        # add caps
        if cap_front:
            self.add(cap_front, "cap_front")
            equivalence_transform(self['main_chain'], self['main_chain']['up'], self['cap_front']['up'])
        else:
            # Hoist port label.
            self.add(main_chain['up'], 'up', containment=False)

        if cap_end:
            self.add(cap_end, 'cap_end')
            equivalence_transform(self['cap_end'], self['cap_end']['up'], self['main_chain']['down'])
        else:
            # Hoist port label.
            self.add(main_chain['down'], 'down', containment=False)


class MainChain(Compound):
    # The main chain for Copolymer. Connect one or more components in a specified sequence.

    def __init__(self, pattern=None, random=False):
        """
        :param pattern: General format for all kinds of copolymers.
            [cap_front, [monomer_A, monomer_B, ...]*10, cap_end]
        :param random: option for random copolymers.
        """

        # default value for pattern
        super().__init__()
        if pattern is None:
            cap_front = CH3()
            cap_end = CH3()
            monomer_A = [CH2(), ("up", "down")] * 2
            monomer_B = [Silane(), "up", "down"] * 2
            pattern = [cap_front, [monomer_A, monomer_B] * 2, cap_end]

        # delete cap_front and cap_end, so we only get the main_chain
        if random:
            main_chain_pattern = list(permutation([monomer for monomers in pattern[1:-1] for monomer in monomers]))
        else:
            main_chain_pattern = [monomer for monomers in pattern[1:-1] for monomer in monomers]
        main_chain_pattern = [monomer_ for monomer in main_chain_pattern for monomer_ in monomer]

        first_part = None
        first_part_label = None
        last_part = None
        last_part_label = None

        for i in range(0, len(main_chain_pattern), 2):
            this_part = clone(main_chain_pattern[i])
            this_part_label = main_chain_pattern[i+1]
            for label in this_part_label:
                assert_port_exists(label, this_part)  # assert labels for the monomer
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
                first_part_label = this_part_label
            else:
                equivalence_transform(this_part,
                                      this_part.labels[this_part_label[1]],
                                      last_part.labels[this_part_label[0]])
            last_part = this_part
            last_part_label = this_part_label

        self.add(first_part.labels[first_part_label[1]], first_part_label[1], containment=False)
        self.add(last_part.labels[last_part_label[0]], last_part_label[0], containment=False)


if __name__ == "__main__":
    cap_front = CH3()
    cap_end = CH3()

    monomer_A = [CH2(), ("up", "down")]*2
    monomer_B = [Silane(), ("up", "down")]*2

    pattern = [cap_front, [monomer_A, monomer_B]*2, cap_end]
    poly = Copolymer(pattern=pattern, random=True)
    poly.visualize(show_ports=False)
