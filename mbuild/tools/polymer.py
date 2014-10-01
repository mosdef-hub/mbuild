from copy import deepcopy

from mbuild.coordinate_transform import equivalence_transform
from mbuild.compound import Compound


class Polymer(Compound):
    """ """
    def __init__(self, proto, port_labels=("up", "down"), n=2):
        """ """
        if n < 1:
            raise Exception('n must be 1 or more')
        Compound.__init__(self)

        assert(isinstance(proto, Compound))
        assert(port_labels[0] in proto.labels)
        assert(port_labels[1] in proto.labels)

        last_part = None
        for body_count in range(0, n):
            this_part = deepcopy(proto)
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # transform this part, such that it's bottom port is rotated+translated to the last part's top port
                equivalence_transform(this_part, this_part.labels[port_labels[1]],
                                      last_part.labels[port_labels[0]])
            last_part = this_part

        # hoist the last part's top port to be the top port of the polymer
        self.add(last_part.labels[port_labels[0]], port_labels[0], containment=False)

        # hoist the first part's bottom port to be the bottom port of the polymer
        self.add(first_part.labels[port_labels[1]], port_labels[1], containment=False)

if __name__ == "__main__":
    from mbuild.examples.alkane.ch2 import Ch2
    ch2 = Ch2()
    m = Polymer(ch2, n=13, port_labels=("up", "down"))
