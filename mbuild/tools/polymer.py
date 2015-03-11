from copy import deepcopy

from mbuild.core import Compound
from mbuild.coordinate_transform import equivalence_transform

__all__ = ['Polymer']


class Polymer(Compound):
    """ """
    def __init__(self, proto, port_labels=("up", "down"), n=2):
        if n < 1:
            raise Exception('n must be 1 or more')
        super(Polymer, self).__init__()

        assert(port_labels[0] in proto.labels)
        assert(port_labels[1] in proto.labels)

        last_part = None
        for body_count in range(0, n):
            this_part = deepcopy(proto)
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
    from mbuild.components.small_groups.ch2 import CH2
    ch2 = CH2()
    poly = Polymer(ch2, n=13, port_labels=("up", "down"))
    poly.visualize(show_ports=True)
