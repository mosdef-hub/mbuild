import numpy as np

__all__ = ['Lattice']


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

    TODO: Add in orientation support similar to lammps implementation possibly

    Parameters
    ----------
    dimension : integer, optional, default 3
        Dimension of the lattice, can be 2 or 3
    lattice_vectors : "SOME ARRAY", optional, default [1,0,0] [0,1,0] [0,0,1]
        3 vectors in the a1, a2, a3 directions that define the dimensionality
        of the unit cell.
        The vectors should NOT be scaled by their lattice spacings.
    lattice_spacings : float, optional, default=[1.0,1.0,1.0]
        3 spacing constants that define the side lengths of the unit cell in
        nanometers.
        Defined by the variables a,b,c; where a corresponds to the a1
        direction, b the a2 direction, and c the a3 direction
    basis_vectors : float, optional, default=[ID,0.0,0.0,0.0]
        Vectors that define location of basis atoms within the unit cell.
        Given as a multiple of the 3 directions 0 >= basis <= 1
        Can define multiple Compounds based on its ID.
    """
    def __init__(self, dimension=None, lattice_vectors=None,
                 lattice_spacings=None, basis_vectors=None):
        super(Lattice, self).__init__()

        self.dimension = None
        self.lattice_vectors = None
        self.lattice_spacings = None
        self.basis_vectors = None
        def _validate_inputs(dimension, lattice_vectors,
                             lattice_spacings, basis_vectors):
            """
            Validate all inputs and either clean up or return errors for bad
            input
            """

            """
            Checking for a dimension given by the user.
            If not found, default to 3D, check for incorrect formatting.
            """
            if dimension is None:
                dimension = 3
            elif not isinstance(dimension, int):
                TypeError('Incorrect type: Dimension {} is not an integer.'
                          .format(dimension))
            elif dimension >= 4 or dimension <= 1:
                raise ValueError('Incorrect dimensions: {} is not an '
                                 'acceptable dimension. '
                                 'Please use 2, or 3 for the dimension.'
                                 .format(dimension))
            else:
                dimension = int(dimension)

            """
            Cleaning up the lattice_vectors input.
            If the input is not provided, assume unit vectors along x,y,z or
            x,y if 2D.
            Check for the user to input a list or a numpy array for the
            lattice_vectors.
            If user inputted, type will be converted to float, and
            dimensionalty will be enforced. Must define all vectors for
            implemented dimension (3x3 for 3D, 2x2 for 2D).
            Will check user implemented for co-linear vectors, these cannot
            define a space filling unit cell.
            Also checking for right-handed lattice_vectors, their
            determinant must be > 0.
            """
            if lattice_vectors is None:
                if dimension is 3:
                    lattice_vectors = np.asarray(([1.0, 0.0, 0.0],
                                                  [0.0, 1.0, 0.0],
                                                  [0.0, 0.0, 1.0]))
                else:
                    lattice_vectors = np.asarray(([1.0, 0.0], [0.0, 1.0]))
            else:
                lattice_vectors = np.asarray(lattice_vectors, dtype=float)
                shape = np.shape(lattice_vectors)
                if dimension is 3:
                    if (3, 3) != shape:
                        raise ValueError('Dimensionality of lattice_vectors is'
                                         ' of {}, not (3,3).' .format(shape))
                    volume = np.linalg.det(lattice_vectors)
                    if abs(volume) == 0.0:
                        raise ValueError('Co-linear vectors {}'
                                         'have a volume of 0.0. Does not '
                                         'define a unit cell.'
                                         .format(lattice_vectors))
                    if volume < 0.0:
                        raise ValueError('Negative Determinant: the volume of'
                                         ' {} is negative, indicating a left-'
                                         'handed system.' .format(volume))
                else:
                    if (2, 2) != shape:
                        raise ValueError('Dimensionality of lattice_vectors is'
                                         ' of {}, not (2,2).' .format(shape))
                    area = np.linalg.det(lattice_vectors)
                    if abs(area) == 0.0:
                        raise ValueError('Co-linear vectors {}'
                                         'have an area of 0.0. Does not '
                                         'define a unit cell.'
                                         .format(lattice_vectors))
                    if area < 0.0:
                        raise ValueError('Negative Determinant: the area of'
                                         ' {} is negative, indicating a left-'
                                         'handed system.' .format(area))

            # lattice_spacings cleaning
            if lattice_spacings is None:
                lattice_spacings = np.asarray(([1.0, 1.0, 1.0]))

            test_spacings = np.asarray(lattice_spacings)
            if np.count_nonzero(test_spacings) != dimension:
                ValueError('Lattice spacings {} does not match dimensionality '
                           '{}.' .format(np.array_str(test_spacings),
                                         dimension))

            return dimension, lattice_vectors, lattice_spacings, basis_vectors
