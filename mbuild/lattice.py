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
    basis_vectors : float, optional, default=[ID,x,y,z] (x,y,z=0)
        Vectors that define location of basis atoms within the unit cell.
        Given as a multiple of the 3 directions 0 >= basis <= 1
        Can define multiple Compounds based on its ID.
        Input as an array either in list form or numpy array
    """
    def __init__(self, dimension=None, lattice_vectors=None,
                 lattice_spacings=None, basis_vectors=None):
        super(Lattice, self).__init__()
        
        # TODO move validate input to class method
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
                raise ValueError('No lattice spacings provided.')

            lattice_spacings = np.asarray(lattice_spacings, dtype=float)
            if np.shape(lattice_spacings) != (dimension, ):
                ValueError('Lattice spacings should be a vector of size: '
                           '{}x0. Please include lattice spacings for each'
                           ' available dimension.'.format(dimension))
            if (lattice_spacings <= 0.0).all():
                ValueError('Negative or zero lattice spacing value. One of the'
                           ' spacings {} is negative or 0, please correct.'
                           .format(lattice_spacings))

            # basis_vectors clean up
            if basis_vectors is None:
                if dimension == 3:
                    basis_vectors = {'one': ([0, 0, 0])}
                elif dimension == 2:
                    basis_vectors = {'one': ([0, 0])}
                else:
                    basis_vectors = {'one': ([0])}
            elif isinstance(basis_vectors, dict):
                temp = basis_vectors.values()
                for val in temp:
                    if not isinstance(val, list):
                        TypeError('Values for basis vector dict not lists.'
                                  ' {} is not a list, please correct.'
                                  .format(val))
                    if len(val) != dimension:
                        ValueError('Dimension and basis vector mismatch. '
                                   ' Basis vector {} does not agree with '
                                   'dimension {}.'.format(val, dimension))
                    for i in range(dimension):
                        if val[i] < 0 and val[i] >= 1:
                            ValueError('Incorrect Basis Vectors. {} should be '
                                       'multiples 0 >= xy(z) < 1.'.format(val))
            else:
                TypeError('Incorrect Type: "basis_vectors" is of type {}, '
                          'not type dict. Please review the documentation for '
                          'proper ' .format(type(basis_vectors)))

            # basis overlap check
            """
            check for two same basis values in dict already
            then check if 2 vectors then overlap
            error if that occurs
            """
            basis_items = basis_vectors.items()
            for key, val in basis_items:
                for key_test, val_test in basis_items:
                    if key != key_test:
                        if val == val_test:
                            ValueError('Duplicate Basis Vectors: Cannot have '
                                       ' 2 equivalent basis vectors.')

                        moved_val_x = 0
                        moved_val_y = 0
                        moved_val_z = 0
                        if dimension == 3:
                            for x in range(2):
                                for y in range(2):
                                    for z in range(2):
                                        moved_val_x = val[0] + x
                                        moved_val_y = val[1] + y
                                        moved_val_z = val[2] + z
                                        if moved_val_x == val_test[0] and 
                                           moved_val_y == val_test[1] and 
                                           moved_val_z == val_test[2]:
                                            ValueError('Duplicate Basis '
                                                       'Vectors: Unit cell '
                                                       'generation impossible'
                                                       '. Key {} and {} will'
                                                       ' generate two '
                                                       'equivalent lattice' 
                                                       'points, causing' 
                                                       ' overlap.' 
                                                       .format(key, key_test))
                        if dimension == 2:
                            for x in range(2):
                                for y in range(2):
                                        moved_val_x = val[0] + x
                                        moved_val_y = val[1] + y
                                        if moved_val_x == val_test[0] and 
                                           moved_val_y == val_test[1]: 
                                            ValueError('Duplicate Basis '
                                                       'Vectors: Unit cell '
                                                       'generation impossible'
                                                       '. Key {} and {} will'
                                                       ' generate two '
                                                       'equivalent lattice' 
                                                       'points, causing' 
                                                       ' overlap.' 
                                                       .format(key, key_test))
                        if dimension == 1:
                            for x in range(2):
                                        moved_val_x = val[0] + x
                                        if moved_val_x == val_test[0]: 
                                            ValueError('Duplicate Basis '
                                                       'Vectors: Unit cell '
                                                       'generation impossible'
                                                       '. Key {} and {} will'
                                                       ' generate two '
                                                       'equivalent lattice' 
                                                       'points, causing' 
                                                       ' overlap.' 
                                                       .format(key, key_test))
            # TODO add in populate function, and replicate include
                
            return {'dim': dimension, 'lat_vec': lattice_vectors,
                    'lat_space': lattice_spacings, 'basis_vec': basis_vectors}
