from os import path

import numpy as np

import mbuild as mb
from mbuild.lib.moieties import CH2, CH3
from mbuild.lib.molecules import Methane, Ethane
from mbuild.lib.atoms import H


class Alkane(mb.Compound):
    """An alkane which may optionally end with a hydrogen or a Port."""
    def __init__(self, n=3, cap_front=True, cap_end=True):
        """Initialize an Alkane Compound.

        Args:
            n: int, optional, default 3
                Number of carbon atoms.
            cap_front: bool, optional, default True
                Add a methyl group to the beginning of chain ('down' port).
            cap_end: bool, optional, default True
                Add methyl groups to the ends of chain ('up' and 'down' port).
        """
        if n < 1:
            raise ValueError('n must be 1 or more')
        super(Alkane, self).__init__()

        # Handle the case of Methane and Ethane separately
        if n < 3:
            if n == 1:
                if cap_front and cap_end:
                    self.add(Methane(), 'chain')
                elif cap_front != cap_end:
                    chain = CH3()
                    self.add(chain, 'chain')
                    if cap_front:
                        self.add(chain['up'], 'down', containment=False)
                    else:
                        self.add(chain['up'], 'up', containment=False)
                else:
                    chain = CH2()
                    self.add(chain)
                    self.add(chain['down'], 'down', containment=False)
                    self.add(chain['up'], 'up', containment=False)
            elif n == 2:
                if cap_front and cap_end:
                    self.add(Ethane(), 'chain')
                elif cap_front != cap_end:
                    chain = CH2()
                    self.add(chain, 'chain')
                    if cap_front:
                        self.add(CH3(), 'methyl_front')
                        mb.force_overlap(move_this=self['chain'],
                                from_positions=self['chain']['up'],
                                to_positions=self['methyl_front']['up'])
                        self.add(chain['down'], 'down', containment=False)
                    else:
                        self.add(CH3(), 'methyl_end')
                        mb.force_overlap(self['methyl_end'],
                                         self['methyl_end']['up'],
                                         self['chain']['down'])
                        self.add(chain['up'], 'up', containment=False)
                else:
                    chain = CH2()
                    self.add(chain, 'chain')
                    self.add(chain['down'], 'down', containment=False)
                    self.add(chain['up'], 'up', containment=False)

        # Handle general case of n >= 3
        else:
            # Adjust length of Polmyer for absence of methyl terminations.
            if cap_front:
                n -= 1
            if cap_end:
                n -= 1
            chain = mb.recipes.Polymer()

            cwd = path.abspath(path.dirname(__file__))
            ch2_path = path.abspath(path.join(cwd, "../moieties/ch2.pdb"))
            ch2 = mb.load(ch2_path)
            ch4 = Methane()
            if cap_front and cap_end:
                chain.add_end_groups(ch4, 1, 0.15)

            chain.add_monomer(
                    ch2,
                    bonding_indices=[0, 0],
                    orientation=[[0, 1, 0], [0, -1, 0]],
                    separation=0.15,
                    replace=False,
                    port_labels=['up','down']
                    )
            chain.build(n)
            self.add(chain, "chain")
            if not cap_front or not cap_end:
                last = [p for p in self.particles()][-1]
                bond = [b for b in self.bonds() if last in b][0]
                anchor = [p for p in bond if p != last][0]
                orientation = last.pos - anchor.pos
                separation = np.linalg.norm(orientation)
                port = mb.Port(
                        anchor=anchor,
                        orientation=orientation,
                        separation=separation
                        )
                self.remove(last)
                self.add(port, 'up')

                last = [p for p in self.particles()][-1]
                bond = [b for b in self.bonds() if last in b][0]
                anchor = [p for p in bond if p != last][0]
                orientation = last.pos - anchor.pos
                separation = np.linalg.norm(orientation)
                port = mb.Port(
                        anchor=anchor,
                        orientation=orientation,
                        separation=separation
                        )
                self.remove(last)
                self.add(port, 'down')

                ch3 = CH3()
                if cap_front:
                    self.add(ch3, 'methyl1')
                    mb.force_overlap(
                            self['methyl1'],
                            self['methyl1']['up'],
                            self['up']
                            )
                #else:
                #    # Hoist port label to Alkane level.
                #    self.add(self["chain"]['up'], 'up', containment=False)
                if cap_end:
                    self.add(mb.clone(ch3), 'methyl2')
                    mb.force_overlap(
                            self['methyl2'],
                            self['methyl2']['up'],
                            self['down']
                            )
                #else:
                #    # Hoist port label to Alkane level.
                #    self.add(self["chain"]['down'], 'down', containment=False)
