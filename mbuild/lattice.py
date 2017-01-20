import numpy as np
from collections import defaultdict
from copy import deepcopy
import mbuild as mb
import pprint

pp = pprint.PrettyPrinter(indent=4)
__all__ = ['Lattice']


class Lattice(object):
    """
    Develop crystal structure from user defined inputs.

    Lattice is the abstract building block a crystal structure composed
    of Compounds. Once defined by the user, the crystal is returned as
    a single compound that can be either replicated through its class
    methods or through a similar replicate Compound method.

    Lattice is defined through normal crystallographic means. Lattice expects
    a right handed lattice and cell edges defined by vectors all originating
    from the origin.

    Parameters
    ----------
    dimension : int, optional, default=3
        Dimension of the system of interest.
    lattice_vectors : numpy array, shape=(dimension, dimension), optional
                      default=([1,0,0], [0,1,0], [0,0,1])
        Vectors that define edges of unit cell corresponding to dimension.
    lattice_spacings : list-like, shape=(dimension,), optional, default=None
        Length of unit cell edges.
    basis_vectors : dictionary-like, shape=(['id',[dimension,]], ...) optional
                    default={('default', ([0,0,0]))}
        Location of all basis Compounds in unit cell.

    Attributes
    ----------
    dimension : int, optional, default=3
        Dimension of system of interest
    lattice_vectors : numpy array, shape=(dimension, dimension), optional
                      default=([1,0,0], [0,1,0], [0,0,1])
        Vectors that define edges of unit cell corresponding to dimension.
    lattice_spacings : list-like, shape=(dimension,), optional, default=None
        Length of unit cell edges.
    basis_vectors : dictionary-like, shape=(['id',[dimension,]], ...) optional
                    default={('default',([0,0,0]))}
        Location of all basis Compounds in unit cell.

    TODO(Justin Gilmer) : populate method with compound input
    TODO(Justin Gilmer) : write_xyz generate .xyz file
    TODO(Justin Gilmer) : write_xyz, more robust xyz file gen
    TODO(Justin Gilmer) : inheritance(Cubic, orthorhombic, hexangonal)
    TODO(Justin Gilmer) : nested for loop cleaning up
    TODO(Justin Gilmer) : orientation functionality
    """
    def __init__(self, dimension=None, lattice_vectors=None,
                 lattice_spacings=None, basis_vectors=None):
        super(Lattice, self).__init__()

        self.validate_inputs(dimension=dimension,
                             lattice_vectors=lattice_vectors,
                             lattice_spacings=lattice_spacings,
                             basis_vectors=basis_vectors)
        self.populate_xyz()
        self.write_xyz()

    def validate_inputs(self, dimension, lattice_vectors,
                        lattice_spacings, basis_vectors):
        """
        Check for proper inputs and set instance attributes.

        validate_inputs takes the data passed to the constructor by the user
        and will ensure that the data is correctly formatted and will then
        set its instance attributes.

        validate_inputs checks that dimensionality is maintained,
        the unit cell is right handed, the area or volume of the unit cell
        is positive and non-zero for 2D and 3D respectively, lattice spacings
        are provided, basis vectors do not overlap when the unit cell is
        expanded.

        Parameters
        ----------
        dimension : int, optional, default=3
            Dimension of system of interest
        lattice_vectors : numpy array, shape=(dimension, dimension), optional
                          default=([1,0,0], [0,1,0], [0,0,1])
            Vectors that define edges of unit cell corresponding to dimension.
        lattice_spacings : list-like, shape=(dimension,), default=None
            Length of unit cell edges.
        basis_vectors : dictionary-like, shape=(['id',[dimension,]], ...)
                        optional default={('default',([0,0,0]))}
            Location of all basis Compounds in unit cell.

        Exceptions Raised
        -----------------
        TypeError : incorrect typing of the input parameters.

        ValueError : values are not within restrictions.

        Call Restrictions
        -----------------
        Called by constructor, should not be callable by user.

        """
        if dimension is None:
            dimension = 3
        else:
            try:
                dimension = int(dimension)
            except ValueError:
                raise
        if dimension < 1 or dimension > 3:
            raise ValueError('Incorrect dimensions: {} is not a proper '
                             'dimension. 1, 2, or 3 are acceptable.'
                             .format(dimension))

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

            if dimension != 1:
                det = np.linalg.det(lattice_vectors)

            if dimension is 3:
                if (3, 3) != shape:
                    raise ValueError('Dimensionality of'
                                     'lattice_vectors is of {},'
                                     ' not (3,3).' .format(shape))
            elif dimension is 2:
                if (2, 2) != shape:
                    raise ValueError('Dimensionality of lattice_vectors is'
                                     ' of {}, not (2,2).' .format(shape))
            else:
                if (1, ) != shape:
                    raise ValueError('Dimensionality of lattice_vectors is'
                                     ' of {}, not (1, ).' .format(shape))

            if abs(det) == 0.0:
                raise ValueError('Co-linear vectors: {}'
                                 'have a determinant of 0.0. Does not '
                                 'define a unit cell.'
                                 .format(lattice_vectors))

            if det < 0.0:
                raise ValueError('Negative Determinant: the determinant of '
                                 '{} is negative, indicating a left-'
                                 'handed system.' .format(det))

        if lattice_spacings is None:
            lattice_spacings = np.asarray([1, 1, 1], dtype=float)

        lattice_spacings = np.asarray(lattice_spacings, dtype=float)
        if np.shape(lattice_spacings) != (dimension, ):
            raise ValueError('Lattice spacings should be a vector of size: '
                             '({},). Please include lattice spacings for each'
                             ' available dimension.'.format(dimension))
        if (lattice_spacings <= 0.0).all():
            raise ValueError('Negative or zero lattice spacing value. One of '
                             'the spacings {} is negative or 0.'
                             .format(lattice_spacings))

        if basis_vectors is None:
            basis_vectors = defaultdict(list)
            if dimension == 3:
                basis_vectors['default'].append((0, 0, 0))
            elif dimension == 2:
                basis_vectors['default'].append((0, 0))
            else:
                basis_vectors['default'].append((0))
        elif (isinstance(basis_vectors, list) or
                isinstance(basis_vectors, tuple)):
            for lst in basis_vectors:
                if len(lst) != 2:
                    raise ValueError('Too many arguments per basis vector. '
                                     '{} was provided, but the format should '
                                     'be: vector = ( (\'ID1\' , [x1,y1,z1]), '
                                     '( \'ID2\' , [x2,y2,z2]) )'.format(lst))
                elif not isinstance(lst[0], str):
                    raise TypeError('ID Error, not a string. {} is not of '
                                    'type str. Please correct in the '
                                    'basis_vectors input {}.'
                                    .format(lst[0], basis_vectors))

                if len(lst[1]) == dimension:
                    if isinstance(lst[1], tuple) or isinstance(lst[1], list):
                        for i in range(dimension):
                            if (lst[1][i] < 0 or lst[1][i] >= 1):
                                raise ValueError('Incorrect basis fractions. '
                                                 '{} is supposed be 0 >= and '
                                                 '< 1.'.format(lst[1]))

                    else:
                        raise TypeError('Incorrect Format: expected <list> '
                                        'or <tuple>. {} was passed in as '
                                        ' type: {} from basis_vectors.'
                                        .format(lst[1], type(lst[1])))
                else:
                    raise ValueError('Basis Vector does not match dimension. '
                                     '{} does not have vector components '
                                     'equal to the dimension specified: {}'
                                     .format(lst[1], dimension))
        else:
            raise TypeError('Incorrect Type: "basis_vectors" is of type {}, '
                            'not type <list> or <tuple>. Please review '
                            'the documentation for proper format.'
                            .format(type(basis_vectors)))

        if not isinstance(basis_vectors, defaultdict):
            overlap_dict = defaultdict(list)
            for the_id, vector in basis_vectors:
                num_iterations = 3
                if dimension == 3:
                    for x in range(num_iterations):
                        for y in range(num_iterations):
                            for z in range(num_iterations):
                                tmpx = vector[0] + x
                                tmpy = vector[1] + y
                                tmpz = vector[2] + z
                                tmp_tuple = tuple((tmpx, tmpy, tmpz))
                                value = tuple((the_id, vector))
                                overlap_dict[tmp_tuple].append(value)

                elif dimension == 2:
                    for x in range(num_iterations):
                        for y in range(num_iterations):
                            tmpx = vector[0] + x
                            tmpy = vector[1] + y
                            tmp_tuple = tuple((tmpx, tmpy))
                            value = tuple((the_id, vector))
                            overlap_dict[tmp_tuple].append(value)

                else:
                    for x in range(num_iterations):
                        tmpx = vector[0] + x
                        tmp_tuple = tuple((tmpx))
                        value = tuple((the_id, vector))
                        overlap_dict[tmp_tuple].append(value)

            # now review for overlap
            for key, val in overlap_dict.items():
                if len(val) > 1:
                    raise ValueError('Overlapping Basis Vectors: Basis '
                                     'vectors overlap when the unit cell is '
                                     'expanded to {}. This is an incorrect '
                                     'perfect lattice. The offending '
                                     'vectors are: {}'
                                     .format(key, val))
        basis_dict = defaultdict(list)
        if not isinstance(basis_vectors, defaultdict):
            for key, val in basis_vectors:
                basis_dict[key].append(val)
        else:
            basis_dict = deepcopy(basis_vectors)

        self.dimension = dimension
        self.lattice_vectors = lattice_vectors
        self.lattice_spacings = lattice_spacings
        self.basis_vectors = basis_dict

    def populate_xyz(self, x=None, y=None, z=None):
        """
        Expand lattice in each respective dimension.

        populate_xyz will expand the lattice based on user input.
        If the values are not provided, populate_xyz will generate a unit
        cell. This is used to generate an xyz file for viewing with commonly
        used chemical viewing software.

        Parameters
        ----------
        x : int, optional, default=None
            How many iterations in the x direction.
        y : int, optional, default=None
            How many iterations in the y direction.
        z : int, optional, default=None
            How many iterations in the z direction.

        Exceptions Raised
        -----------------
        ValueError : incorrect x,y, or z values.
        TypeError : incorrect type for basis vector

        Call Restrictions
        -----------------
        Called after constructor by user.

        """
        if self.dimension == 3:
            a = self.lattice_spacings[0]
            b = self.lattice_spacings[1]
            c = self.lattice_spacings[2]
            if x is None:
                x = 1
            if y is None:
                y = 1
            if z is None:
                z = 1
            if x < 1 or y < 1 or z < 1:
                raise ValueError('Incorrect populate value: X, Y, or Z is < 1.'
                                 ' Cannot replicate unit cell less than 1')
        elif self.dimension == 2:
            a = self.lattice_spacings[0]
            b = self.lattice_spacings[1]
            if x is None:
                x = 1
            if y is None:
                y = 1
            if z is None:
                pass
            else:
                raise ValueError('Z is defined although dimension is 2D')
            if x < 1 or y < 1:
                raise ValueError('Incorrect populate value: X or Y is < 1. '
                                 ' Cannot replicate unit cell less than 1')
        elif self.dimension == 1:
            a = self.lattice_spacings[0]
            if x is None:
                x = 1
            if y is None:
                pass
            else:
                raise ValueError('Y is defined although dimension is 1D')
            if z is None:
                pass
            if z is not None:
                raise ValueError('Z is defined although dimension is 2D')
            if x < 1:
                raise ValueError('Incorrect populate value: X < 1. '
                                 ' Cannot replicate unit cell less than 1')
        else:
            raise ValueError('Dimension not defined.')

        if isinstance(self.basis_vectors, defaultdict):
            pass
        else:
            raise TypeError('Basis vector is not of type defaultdict')

        cell = defaultdict(list)
        for key, val in self.basis_vectors.items():
            for val_item in range(len(val)):
                if self.dimension == 3:
                    for i in range(x+1):
                        for j in range(y+1):
                            for k in range(z+1):
                                tmpx = (val[val_item][0] + i) * a
                                tmpy = (val[val_item][1] + j) * b
                                tmpz = (val[val_item][2] + k) * c
                                tmp_tuple = tuple((tmpx, tmpy, tmpz))
                                cell[key].append(((tmp_tuple)))
                elif self.dimension == 2:
                    for i in range(x+1):
                        for j in range(y+1):
                            tmpx = (val[val_item][0] + i) * a
                            tmpy = (val[val_item][1] + j) * b
                            tmp_tuple = tuple((tmpx, tmpy))
                            cell[key].append(((tmp_tuple)))
                else:
                    for i in range(x+1):
                        tmpx = (val[val_item][0] + i) * a
                        tmp_tuple = tuple((tmpx))
                        cell[key].append(((tmp_tuple)))
        self.cell = cell
        pp.pprint(cell)

    def write_xyz(self, moleculeName=None, fileName=None):
        """
        Write formatted file in XYZ format of Lattice.

        Parameters
        ----------
        moleculeName : str, optional, default=None
            Name of the molecule created from the Lattice instance.
        fileName : str, optional, default=latticeOutput.xyz
            Filename for the XYZ file generated.

        """

        if moleculeName is None:
            moleculeName = ''
        else:
            try:
                moleculeName = str(moleculeName)
            except Exception as e:
                print('String Conversion Error: moleculeName {} '
                      'can not be converted to type string: {}'
                      .format(moleculeName, e))
                raise

        if fileName is None:
            fileName = 'latticeOutput.xyz'
        else:
            try:
                fileName = str(fileName)
            except Exception as e:
                print('String Conversion Error: fileName {} '
                      'can not be converted to type string: {}'
                      .format(fileName, e))
                raise
        num_atoms = 0
        outstr = ""
        for val in self.cell.values():
            num_atoms = num_atoms + len(val)
        with open(fileName, 'w') as fout:
            fout.write(str(num_atoms) + '\n')
            fout.write(moleculeName + '\n')
            for key, val in self.cell.items():
                for pos in range(len(val)):
                    outstr = ""
                    outstr = str(key) + str(' ')
                    for i in range(self.dimension):
                        outstr = outstr + str(val[pos][i] * 10) + str(' ')

                    fout.write(outstr + '\n')
        if not fout.closed:
            fout.close()



def main():
    # basis_vec = (['Cl', (0, 0, 0)], ['Cs', (.5, .5, .5)], ['He', (0, 0, 0)])
    # the_lattice = Lattice(lattice_spacings=[1, 1, 1], basis_vectors=basis_vec)
    the_lattice = Lattice(lattice_spacings=[1, 1, 1])
    the_lattice = Lattice()
    lat_space = [.4123, .4123, .4123]
    basis_vec = [('Cl', [0, 0, 0]), ('Cs', [.5, .5, .5])]
    the_lattice = Lattice(lattice_spacings=lat_space, basis_vectors=basis_vec)

if __name__ == "__main__":
    main()
