import numpy as np
from collections import defaultdict

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

        tmp_dict = self._validate_inputs(dimension=dimension,
                                         lattice_vectors=lattice_vectors,
                                         lattice_spacings=lattice_spacings,
                                         basis_vectors=basis_vectors)
        self.dimension = tmp_dict.get('dimension')
        self.lattice_vectors = tmp_dict.get('lattice_vectors')
        self.lattice_spacings = tmp_dict.get('lattice_spacings')
        self.basis_vectors = tmp_dict.get('basis_vectors')
        tmp_dict.clear()

    def _validate_inputs(self, dimension, lattice_vectors,
                         lattice_spacings, basis_vectors):
        """
        Validate all inputs and either clean up or
        return errors for bad input
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
        elif dimension >= 4 or dimension < 1:
            raise ValueError('Incorrect dimensions: {} is not an '
                             'acceptable dimension. '
                             'Please use 1, 2, or 3 for the dimension.'
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
            elif dimension is 2:
                lattice_vectors = np.asarray(([1.0, 0.0], [0.0, 1.0]))
            else:
                lattice_vectors = np.asarray(([1.0]))
        else:
            lattice_vectors = np.asarray(lattice_vectors, dtype=float)
            shape = np.shape(lattice_vectors)
            if dimension is 3:
                if (3, 3) != shape:
                    raise ValueError('Dimensionality of'
                                     'lattice_vectors is of {},'
                                     ' not (3,3).' .format(shape))
                volume = np.linalg.det(lattice_vectors)
                if abs(volume) == 0.0:
                    raise ValueError('Co-linear vectors {}'
                                     'have a volume of 0.0.'
                                     'Does not define a unit cell.'
                                     .format(lattice_vectors))
                if volume < 0.0:
                    raise ValueError('Negative Determinant: '
                                     'the volume of {} is '
                                     'negative, indicating a '
                                     'left-handed system.'
                                     .format(volume))
            elif dimension is 2:
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
                    raise ValueError('Negative Determinant: the area of '
                                     '{} is negative, indicating a left-'
                                     'handed system.' .format(area))
            else:
                if (1, ) != shape:
                    raise ValueError('Dimensionality of lattice_vectors is'
                                     ' of {}, not (1, ).' .format(shape))

        # lattice_spacings cleaning
        if lattice_spacings is None:
            raise ValueError('No lattice spacings provided.')

        lattice_spacings = np.asarray(lattice_spacings, dtype=float)
        if np.shape(lattice_spacings) != (dimension, ):
            ValueError('Lattice spacings should be a vector of size: '
                       '({},). Please include lattice spacings for each'
                       ' available dimension.'.format(dimension))
        if (lattice_spacings <= 0.0).all():
            ValueError('Negative or zero lattice spacing value. One of the'
                       ' spacings {} is negative or 0, please correct.'
                       .format(lattice_spacings))

        # basis_vectors clean up
        if basis_vectors is None:
            basis_vectors = defaultdict(list)
            if dimension == 3:
                basis_vectors['default'].append(([0, 0, 0]))
            elif dimension == 2:
                basis_vectors['default'].append(([0, 0]))
            else:
                basis_vectors['default'].append(([0]))
        elif (isinstance(basis_vectors, list) or
                isinstance(basis_vectors, tuple)):
            for lst in basis_vectors:
                if len(lst) != 2:
                    ValueError('Too many arguments per basis vector. '
                               '{} was provided, but the format should be: '
                               ' vector = ( (\'ID1\' , [x1,y1,z1]), '
                               '( \'ID2\' , [x2,y2,z2]) )'.format(lst))
                # input check ensure ID is a string, and atom location is list
                elif not isinstance(lst[0], str):
                    TypeError('ID Error, not a string. {} is not of '
                              'type str. Please correct in the '
                              'basis_vectors input {}.'
                              .format(lst[0], basis_vectors))
                else:
                    pass

                if len(lst[1]) == dimension:
                    if (isinstance(lst[1], tuple) or isinstance(lst[1], list)):
                        for i in range(dimension):
                            if (lst[1][i] < 0 or lst[1][i] >= 1):
                                ValueError('Incorrect basis fractions. '
                                           '{} is supposed be 0 >= and < 1.'
                                           .format(lst[1]))
                            else:
                                pass

                    else:
                        TypeError('Incorrect Format: expected <list> '
                                  'or <tuple>. {} was passed in as type: {} '
                                  'from basis_vectors.'
                                  .format(lst[1], type(lst[1])))
                else:
                    ValueError('Basis Vector does not match dimension. '
                               '{} does not have vector components equal to '
                               'the dimension specified: {}'
                               .format(lst[1], dimension))
        else:
            TypeError('Incorrect Type: "basis_vectors" is of type {}, '
                      'not type <list> or <tuple>. Please review '
                      'the documentation for proper format.'
                      .format(type(basis_vectors)))

        # basis overlap check
        """
        check for two same basis values in dict already
        then check if 2 vectors then overlap
        error if that occurs
        Achieved by expanding unit cell out 2 iterations
        """
        if isinstance(basis_vectors, defaultdict):
            pass
        else:
            # test for overlap
            overlap_dict = defaultdict(list)
            for the_id, vector in basis_vectors:
                if dimension == 3:
                    for x in range(3):
                        for y in range(3):
                            for z in range(3):
                                tmpx = vector[0] + x
                                tmpy = vector[1] + y
                                tmpz = vector[2] + z
                                tmp_tuple = tuple((tmpx, tmpy, tmpz))
                                value = tuple((the_id, vector))
                                overlap_dict[tmp_tuple].append(value)

                elif dimension == 2:
                    for x in range(3):
                        for y in range(3):
                            tmpx = vector[0] + x
                            tmpy = vector[1] + y
                            tmp_tuple = tuple((tmpx, tmpy))
                            value = tuple((the_id, vector))
                            overlap_dict[tmp_tuple].append(value)

                else:
                    for x in range(3):
                        tmpx = vector[0] + x
                        tmp_tuple = tuple((tmpx))
                        value = tuple((the_id, vector))
                        overlap_dict[tmp_tuple].append(value)

            # now review for overlap
            for key, val in overlap_dict.items():
                if len(val) > 1:
                    ValueError('Overlapping Basis Vectors: Basis vectors '
                               'overlap when the unit cell is '
                               'expanded to {}. This is an incorrect perfect '
                               'lattice. The vectors are: {}'
                               .format(key, val))

        return {'dim': dimension, 'lat_vec': lattice_vectors,
                'lat_space': lattice_spacings, 'basis_vec': basis_vectors}

#    def populate()  # TODO, determine input parameters)
