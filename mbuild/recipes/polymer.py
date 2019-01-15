import itertools as it
from copy import deepcopy

from mbuild.compound import Compound
from mbuild.coordinate_transform import force_overlap
from mbuild.utils.validation import assert_port_exists
from mbuild import clone


__all__ = ['Polymer']


class Polymer(Compound):
    """Connect one or more components in a specified sequence.

    Parameters
    ----------
    monomers : mb.Compound or list of mb.Compound
        The compound(s) to replicate.
    n : int
        The number of times to replicate the sequence.
    sequence : str, optional, default='A'
        A string of characters where each unique character represents one
        repetition of a monomer. Characters in `sequence` are assigned to
        monomers in the order assigned by the built-in `sorted()`.
    port_labels : 2-tuple of strs, optional, default=('up', 'down')
        The names of the two ports to use to connect copies of proto.
    caps : Compound or 2-tuple of Compounds, optional, default=(None, None)
        The cap compounds should contain at least one open port. The first
        argument will cap the beginning of the pattern i.e. C->ABAB and
        the second argument caps the end of the pattern i.e. ABAB<-C.
        If a Compound is given it will be used to cap both ends.
    cap_ports : str 2-tuple of strs, optional, default=(None, None)
        The names of the ports from the cap compound(s) to connect to polymer.
        If a str is given that port will be used for both ends.
    """

    def __init__(self, monomers, n, sequence='A', port_labels=(
            'up', 'down'), caps=(None, None), cap_ports=(None, None)):
        if n < 1:
            raise ValueError('n must be 1 or more')
        super(Polymer, self).__init__()
        if isinstance(monomers, Compound):
            monomers = (monomers,)
        for monomer in monomers:
            for label in port_labels:
                assert_port_exists(label, monomer)

        unique_seq_ids = sorted(set(sequence))

        if len(monomers) != len(unique_seq_ids):
            raise ValueError('Number of monomers passed to `Polymer` class must'
                             ' match number of unique entries in the specified'
                             ' sequence.')

        # 'A': monomer_1, 'B': monomer_2....
        seq_map = dict(zip(unique_seq_ids, monomers))

        last_part = None
        for n_added, seq_item in enumerate(it.cycle(sequence)):
            this_part = clone(seq_map[seq_item])
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that it's bottom port is rotated
                # and translated to the last part's top port.
                force_overlap(this_part,
                              this_part.labels[port_labels[1]],
                              last_part.labels[port_labels[0]])
            last_part = this_part
            if n_added == n * len(sequence) - 1:
                break

        # Hoist the last part's top port to be the top port of the polymer.
        self.add(last_part.labels[port_labels[0]],
                 port_labels[0], containment=False)

        # Hoist the first part's bottom port to be the bottom port of the
        # polymer.
        self.add(first_part.labels[port_labels[1]],
                 port_labels[1], containment=False)

        if isinstance(caps, Compound):
            caps = (caps, deepcopy(caps))
        if isinstance(cap_ports, str):
            cap_ports = (cap_ports, cap_ports)

        if len(caps) != 2 and len(cap_ports) != 2:
            raise ValueError('Two viable caps and cap ports must be provided')

        # Check for duplicate ports and rename to avoid ambiguity if nessecary.
        if caps[0] is not None and caps[1] is not None:
            u_ports = set(caps[0].available_ports())
            for port in u_ports:
                if port in set(caps[1].available_ports()):
                    for label in caps[0].labels[port]
                        label.name = "lcap_"+label.name
                    for label in caps[1].labels[port]
                        label.name = "rcap_"+label.name

        # Add a cap to left/right or both ends of the pattern i.e. C->ABABAB
        for cap, cap_port, end_port in zip(
                caps,
                cap_ports,
                reversed(port_labels)):
            if cap is not None:
                self.add(cap)
                force_overlap(cap,
                              cap[cap_port],
                              self[end_port])

        



if __name__ == "__main__":
    from mbuild.lib.moieties import CH2
    ch2 = CH2()
    poly = Polymer(ch2, n=13, port_labels=("up", "down"))
    poly.save('polymer.mol2')
