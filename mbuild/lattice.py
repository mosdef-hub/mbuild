import numpy as np


class Lattice(object):
    """A building block for regular, homogeneous or heterogeneous materials.

    Lattice is the abstract building block of any object that exhibits
    long-range order. More specfically, any object that has particles arranged
    in a predictable, infinitely repeating pattern with the particles
    occupying the cell.
    These lattice sites are related to the origin by a set of translation
    vectors.

    The long range lattice can be reduced to a smaller, more manageable "unit
    cell". In other words, a smaller unit of the lattice, when repeated
    infinitely, will generate the long range lattice.

    Parameters
    ----------
    lattice_type : str, optional, default="P"
        Type of lattice centering present (P, I, F, A, B, C).
    lattice_vectors : "SOME ARRAY", optional, default [1,0,0] [0,1,0] [0,0,1]
        3 vectors in the a1, a2, a3 directions that define the dimensionality
        of the unit cell.
    lattice_spacings: float, optional, default=[1.0,1.0,1.0]
        3 spacing constants that define the side lengths of the unit cell in
        nanometers.
    """
    def __init__(self, dimension=None, lattice_vectors=None,
                 lattice_spacings=None, basis_vectors=None):
        super(Lattice, self).__init__()

        def _validate_inputs(dimension, lattice_vectors,
                             lattice_spacings, basis_vectors):
            """
            Validate all inputs and either clean up or return errors for bad
            input
            """

            # define a 3D unit cell if not specified
            if dimension is None:
                dimension = 3

            if type(dimension) is not int():
                TypeError('Incorrect type: Dimension {} is not an integer.'
                          .format(dimension))

            if dimension >= 4 or dimension <= 0:
                raise ValueError('Incorrect dimensions: {} is not an '
                                 'acceptable dimension. '
                                 'Please use 1, 2, or 3 for the dimension.'
                                 .format(dimension))

            # lattice vector cleaning
            if lattice_vectors is None:
                if dimension is 3:
                    lattice_vectors = np.asarray([1.0, 0.0, 0.0],
                                                 [0.0, 1.0, 0.0],
                                                 [0.0, 0.0, 1.0])
                elif dimension is 2:
                    lattice_vectors = np.asarray([1.0, 0.0], [0.0, 1.0])
                else:
                    lattice_vectors = np.asarray([1.0])

            # lattice_spacings cleaning
            if lattice_spacings is None:
                lattice_spacings = np.asarray([1.0, 1.0, 1.0])

            test_spacings = np.asarray(lattice_spacings)
            if np.count_nonzero(test_spacings) != dimension:
                ValueError('Lattice spacings {} does not match dimensionality '
                           '{}.' .format(np.array_str(test_spacings),
                                         dimension))

            return dimension, lattice_vectors, lattice_spacings, basis_vectors
