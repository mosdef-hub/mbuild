from copy import copy


class Rings(object):
    """Find all rings of a specified length that the atom is a part of.

    Note: Finds each ring twice because the graph is traversed in both directions.
    """
    def __init__(self, atom, ring_length):
        """Initialize a ring bearer. """
        self.rings = list()
        self.current_path = list()
        self.ring_length = ring_length

        self.current_path.append(atom)
        self.step(atom)

    def step(self, atom):
        neighbors = atom.neighbors
        if len(neighbors) > 1:
            for n in neighbors:
                # Check to see if we found a ring.
                current_length = len(self.current_path)
                if current_length > 2 and n == self.current_path[0]:
                    self.rings.append(copy(self.current_path))
                # Prevent stepping backwards.
                elif n in self.current_path:
                    continue
                else:
                    if current_length < self.ring_length:
                        # Take another step.
                        self.current_path.append(n)
                        self.step(n)
                    else:
                        # Reached max length.
                        continue
            else:
                # Finished looping over all neighbors.
                del self.current_path[-1]
        else:
            # Found a dead end.
            del self.current_path[-1]


def benzene(atom):
    """Check if atom is part of a single benzene ring. """
    ring = Rings(atom, 6).rings
    # 2 rings, because we count the traversal in both directions.
    if len(ring) == 2:
        for at in ring[0]:
            if not (at.kind == 'C' and len(at.neighbors) == 3):
                break
        else:
            return ring[0]  # Only return one direction of the ring.
    return False


def dioxolan13(atom):
    """Check if the atom is part of a single 1,3-dioxolan ring. """
    ring = Rings(atom, 5).rings
    # 2 rings, because we count the traversal in both directions.
    if len(ring) == 2:
        sequence = []
        for at in ring[0]:
            if not (at.kind == 'C' or at.kind == 'O'):
                break
            sequence.append(at.kind)
        if ''.join(sequence) in ['CCOCO', 'COCOC', 'OCOCC', 'COCCO', 'OCCOC']:
            return True
    return False