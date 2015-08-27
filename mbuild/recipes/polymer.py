from mbuild.compound import Compound
from mbuild.coordinate_transform import equivalence_transform
from mbuild.utils.validation import assert_port_exists
from mbuild import clone
__all__ = ['Polymer']


class Polymer(Compound):
    """Connect a component to successive copies of itself.

    Parameters
    ----------
    proto : mb.Compound
        The compound to replicate.
    port_labels : 2-tuple of strs
        The names of the two ports to use to connect copies of proto.
    n : int, optional, default=2
        The number of times to replicate proto.

    """
    def __init__(self, proto, port_labels=("up", "down"), n=2):
        if n < 1:
            raise Exception('n must be 1 or more')
        super(Polymer, self).__init__()

        for label in port_labels:
            assert_port_exists(label, proto)

        last_part = None
        for body_count in range(0, n):
            this_part = clone(proto)
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that it's bottom port is rotated
                # and translated to the last part's top port.
                equivalence_transform(this_part, this_part.labels[port_labels[1]],
                                      last_part.labels[port_labels[0]])
            last_part = this_part

        # Hoist the last part's top port to be the top port of the polymer.
        self.add(last_part.labels[port_labels[0]], port_labels[0], containment=False)

        # Hoist the first part's bottom port to be the bottom port of the polymer.
        self.add(first_part.labels[port_labels[1]], port_labels[1], containment=False)

if __name__ == "__main__":
    from mbuild.lib.moieties import CH2
    ch2 = CH2()
    poly = Polymer(ch2, n=13, port_labels=("up", "down"))
    poly.visualize(show_ports=True)
