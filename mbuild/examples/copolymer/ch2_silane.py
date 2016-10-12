# -- ==alkane== --
import mbuild as mb

from mbuild.lib.moieties import CH2
from mbuild.lib.moieties import CH3
from mbuild.lib.moieties import Silane

from mbuild.recipes.copolymer import Copolymer


class Ch2_Sliane(mb.Compound):
    """Ch2_Sliane copolymer."""
    cap_front = CH3()  # front cap of main chain
    cap_end = CH3()  # end cap of main chain

    chain_A = CH2()  # polymer A in the main chain
    chain_B = Silane()  # polymer B in the main chain

    ch2_sliane_pattern = [cap_front, [chain_A, 2, ("up", "down")], [chain_B, 2, ("up", "down")], 1, cap_end]

    def __init__(self, pattern=ch2_sliane_pattern, random=False):
        """
        :param pattern: General format for all kinds of copolymers.
            [cap_front, [A, A_num, A_port_labels], [B, B_num, B_port_labels], cap_end, repeat_num]
            Alternating copolymers: A_num = 1, B_num = 1, repeat_num = any_value
            Periodic copolymers: A_num = any_value, B_num = any_value, repeat_num = any_value
            Block copolymers: A_num = any_value, B_num = any_value, repeat_num = 1
            Statistical copolymers: repeat A block and B block randomly with repeat_num times, if random is True
        :param random: option for Statistical copolymers, which are actually random copolymers.
        """

        super(Ch2_Sliane, self).__init__()
        ch2_sliane_copolymer = Copolymer(pattern, random)

        self.add(ch2_sliane_copolymer, 'ch2_sliane_copolymer')

# -- ==ch2_sliane== --
