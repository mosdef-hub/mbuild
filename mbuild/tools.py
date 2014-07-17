from atom import Atom
from compound import Compound

__author__ = 'sallai'

import numpy as np

def createEquivalenceTransform(equiv):
    """Compute an equivalence transformation that transforms this compound
    to another compound's coordinate system.

    :param other: the other point cloud
    :param equiv: list of equivalent points
    :returns: the coordinatetransform object that transforms this point cloud to the other point cloud's coordinates system
    """


    self_points = np.array([])
    self_points.shape = (0, 3)
    other_points = np.array([])
    other_points.shape = (0, 3)

    for pair in equiv:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise Exception('Equivalence pair not a 2-tuple')
        if not ((isinstance(pair[0], Compound) and isinstance(pair[1], Compound)) or (
                isinstance(pair[0], Atom) and isinstance(pair[1], Atom))):
            raise Exception(
                'Equivalence pair type mismatch: pair[0] is a ' + str(type(pair[0])) + ' and pair[1] is a ' + str(
                    type(pair[1])))

        if isinstance(pair[0], Atom):
            self_points = np.vstack([self_points, pair[0].pos])
            other_points = np.vstack([other_points, pair[1].pos])
        if isinstance(pair[0], Compound):
            for atom0 in pair[0].atoms():
                self_points = np.vstack([self_points, atom0.pos])
            for atom1 in pair[1].atoms():
                other_points = np.vstack([other_points, atom1.pos])

    from coordinate_transform import *
    T = RigidTransform(self_points, other_points)
    return T


if __name__ == "__main__":
    print "hello"



