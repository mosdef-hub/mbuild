import itertools as it

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

    """
    def __init__(self, monomers, n, sequence='A', port_labels=('up', 'down')):
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
        self.add(last_part.labels[port_labels[0]], port_labels[0], containment=False)

        # Hoist the first part's bottom port to be the bottom port of the polymer.
        self.add(first_part.labels[port_labels[1]], port_labels[1], containment=False)

if __name__ == "__main__":
    from mbuild.lib.moieties import CH2
    ch2 = CH2()
    poly = Polymer(ch2, n=13, port_labels=("up", "down"))
    poly.save('polymer.mol2')
