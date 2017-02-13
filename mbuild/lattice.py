__all__ = ['Lattice']


import numpy as np
from collections import defaultdict
from copy import deepcopy
from six import string_types
import itertools as it
import mbuild as mb


class Lattice(object):
    """Develop crystal structure from user defined inputs.

    Lattice, the abstract building block of a crystal cell.
    Once defined by the user, the crystal is returned as
    a single Compound that can be either replicated through its class
    methods or through a similar replicate Compound method.

    Lattice is defined through the standard bravais lattices, which have been
    accepted by the International Union of Crystallography.
    A Lattice can be fully described with its lattice vectors and lattice
    spacings. Also, the Lattice can be fully defined by its lattice parameters:
    the lattice spacings and its set of coordinate angles will then
    generate the lattice vectors. Lattice expects a right handed lattice and
    cell edges defined by vectors all originating from the origin in
    Cartesian space.

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
    angles : list-like,  shape=(dimension,), optional, default=None
        Interplanar angles describing unit cell.

    Attributes
    ----------
    dimension : int, optional, default=3
        Dimension of system of interest
    lattice_vectors : numpy array, shape=(dimension, dimension), optional
                      default=([1,0,0], [0,1,0], [0,0,1])
        Vectors that define edges of unit cell corresponding to dimension.
    lattice_spacings : list-like, shape=(dimension,), required, default=None
        Length of unit cell edges.
    basis_vectors : list-like, shape=(['id',[dimension,]], ... ,) optional
                    default={('default',([0,0,0]))}
        Location of all basis Compounds in unit cell.
    angles : list-like, optional, default=None
        Lattice angles to define Bravais Lattice.

    Examples
    --------
    Generating a triclinc lattice for cholesterol.

    >>> import mbuild as mb
    >>> from mbuild.utils.io import get_fn
    >>> # reading in the lattice parameters for crystalline cholesterol
    >>> angle_values = [94.64, 90.67, 96.32]
    >>> spacings = [1.4172, 3.4209, 1.0481]
    >>> basis_vector = ( ('cholesterol', [0,0,0]), )
    >>> cholesterol_lattice = mb.Lattice(spacings,
    ...                                  angles=angle_values,
    ...                                  basis_vectors=basis_vector,
    ...                                  dimension=3)

    The lattice based on the bravais lattice parameters of crystalline
    cholesterol was generated.

    Replicating the triclinic unit cell out 3 in x,y,z directions.
    >>> cholesterol_unit = mb.Compound()
    >>> cholesterol_unit = mb.load(get_fn('cholesterol.pdb'))
    >>> # associate basis vector with id 'cholesterol' to cholesterol Compound
    >>> basis_dictionary = {'cholesterol' : cholesterol_unit}
    >>> expanded_cell = cholesterol_lattice.populate(x=3, y=3, z=3,
    ...                              compound_dict=basis_dictionary)

    The unit cell of cholesterol was associated with a Compound that contains
    the connectivity data and spatial arrangements of a cholesterol molecule.
    The unit cell was then expanded out in x,y,z directions and cholesterol
    Compounds were populated.


    Generating BCC CsCl crystal structure
    >>> import mbuild as mb
    >>> chlorine = mb.Compound(name='Cl')
    >>> # angles not needed, when not provided, defaults to 90,90,90
    >>> cesium = mb.Compound(name='Cs')
    >>> spacings = [.4123, .4123, .4123]
    >>> basis_vector = ( ('Cl', [0,0,0]), ('Cs', [.5, .5, .5]), )
    >>> cscl_lattice = mb.Lattice(spacings, basis_vectors=basis_vector,
    ...                           dimension=3)

    Now associate id with Compounds for basis atoms and replicate 3x3x3
    >>> cscl_dict = {'Cl' : chlorine, 'Cs' : cesium}
    >>> cscl_compound = cscl_lattice.populate(x=3, y=3, z=3,
    ...                                       compound_dict=cscl_dict)

    A multi-Compound basis was created and replicated. For each unique basis
    atom position, a separate entry must be completed for the basis_vector
    input.

    Generating FCC Copper cell with lattice_vectors instead of angles
    >>> import mbuild as mb
    >>> copper = mb.Compound(name='Cu')
    >>> lattice_vector = ( [1, 0, 0], [0, 1, 0], [0, 0, 1])
    >>> spacings = [.36149, .36149, .36149]
    >>> basis_vector = ( ('Cu', [0, 0, 0]), ('Cu', [.5, .5, 0]),
    ...                 ('Cu', [.5, 0, .5] ), ('Cu', [0, .5, .5]), )
    >>> copper_lattice = mb.Lattice(spacings, dimension=3,
    ...                           lattice_vectors=lattice_vector,
    ...                           basis_vectors=basis_vector)
    >>> copper_dict = {'Cu' : copper}
    >>> copper_cell = copper_lattice.populate(x=3, y=3, z=20,
    ...                                       compound_dict=copper_dict)

    TODO(Justin Gilmer) : migrate data cleaning to separate functions
    TODO(Justin Gilmer) : Print function to display info about Lattice (repr)
    TODO(Justin Gilmer) : inheritance(Cubic, orthorhombic, hexangonal)
    TODO(Justin Gilmer) : nested for loop cleaning up
    TODO(Justin Gilmer) : orientation functionality
    TODO(Justin Gilmer) : conversion to idiomatic python
    """
    def __init__(self, lattice_spacings, dimension=None,
                 lattice_vectors=None, basis_vectors=None,
                 angles=None):
        super(Lattice, self).__init__()
        self._validate_inputs(lattice_vectors=lattice_vectors,
                              dimension=dimension,
                              lattice_spacings=lattice_spacings,
                              basis_vectors=basis_vectors,
                              angles=angles)

    def _validate_inputs(self, lattice_vectors, dimension,
                         lattice_spacings, basis_vectors, angles):
        """Check for proper inputs and set instance attributes.

        validate_inputs takes the data passed to the constructor by the user
        and will ensure that the data is correctly formatted and will then
        set its instance attributes.

        validate_inputs checks that dimensionality is maintained,
        the unit cell is right handed, the area or volume of the unit cell
        is positive and non-zero for 2D and 3D respectively, lattice spacings
        are provided, basis vectors do not overlap when the unit cell is
        expanded.

        Exceptions Raised
        -----------------
        TypeError : incorrect typing of the input parameters.

        ValueError : values are not within restrictions.
        """
        if dimension is None:
            dimension = 3
        else:
            dimension = int(dimension)
        if dimension < 1 or dimension > 3:
            raise ValueError('Incorrect dimensions: {} is not a proper '
                             'dimension. 1, 2, or 3 are acceptable.'
                             .format(dimension))

        if lattice_spacings is not None:
            lattice_spacings = np.asarray(lattice_spacings, dtype=float)
            if np.shape(lattice_spacings) != (dimension, ):
                raise ValueError('Lattice spacings should be a vector of '
                                 'size:({},). Please include lattice spacings '
                                 'for each available dimension.'
                                 .format(dimension))
        else:
            raise ValueError('Lattice Spacing Issue: None provided, '
                             'must provide lattice spacings matching '
                             'the dimension ({}) of the system.'
                             .format(dimension))
        if np.any(lattice_spacings <= 0.0):
            raise ValueError('Negative or zero lattice spacing value. One of '
                             'the spacings {} is negative or 0.'
                             .format(lattice_spacings))

        if angles is not None:
            if (len(angles), dimension) == (3, 3):
                if all(isinstance(theAngle, float) for theAngle in angles):
                    if sum(angles) < 360.0 or sum(angles) > -360.0:
                        for theAngle in angles:
                            if(theAngle != 180.0 and theAngle != 0.0):
                                pass
                            else:
                                raise ValueError('Angles cannot be 180.0 or '
                                                 '0.0.')
                    else:
                        raise ValueError('Angles sum to a value greater than '
                                         '360.0 or less than -360.0.')
                else:
                    raise TypeError('Angles are not type float.')

                for subset in it.permutations(angles, 3):
                    if not subset[0] < sum(angles) - subset[0]:
                        raise ValueError('Each angle provided must be less '
                                         'than the sum of the other two '
                                         'angles. {} is greater.'
                                         .format(subset[0]))

            elif len(angles) == 1 and dimension == 2:
                if all(isinstance(theAngle, float) for theAngle in angles):
                    for theAngle in angles:
                        if (theAngle != 180.0 and theAngle != 0.0 and
                                theAngle < 180.0 and theAngle > -180.0):
                            pass
                        else:
                            raise ValueError('Angle incorrectly defined. {} '
                                             'does not follow the proper '
                                             'guidelines for a bravais angle. '
                                             'Refer to documentation.'
                                             .format(theAngle))
            else:
                raise ValueError('Incorrect amount of angles provided for '
                                 'dimension {}. Recieved {} angles.'
                                 .format(dimension, len(angles)))

        if angles is not None and lattice_vectors is not None:
            return ValueError('Over defined system. Lattice vectors and '
                              'angles passed to constructor. Only one of '
                              'these sets are required.')
        if angles is not None and dimension != 1 and lattice_vectors is None:
            lattice_vectors = self._from_lattice_parameters(angles, dimension)

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

            if dimension in (3, 2):
                if (dimension, dimension) != shape:
                    raise ValueError('Dimensionality of lattice_vectors is '
                                     ' of shape {} not {}.'
                                     .format(shape, (dimension, dimension)))
            else:
                if (1, ) != shape:
                    raise ValueError('Dimensionality of lattice_vectors is'
                                     ' of {}, not (1, ).' .format(shape))
            if dimension != 1:
                det = np.linalg.det(lattice_vectors)
                if abs(det) == 0.0:
                    raise ValueError('Co-linear vectors: {}'
                                     'have a determinant of 0.0. Does not '
                                     'define a unit cell.'
                                     .format(lattice_vectors))

                if det <= 0.0:
                    raise ValueError('Negative Determinant: the determinant '
                                     'of {} is negative, indicating a left-'
                                     'handed system.' .format(det))

        if basis_vectors is None:
            basis_vectors = defaultdict(list)
            basis_vectors['default'].append((0,) * dimension)
        elif isinstance(basis_vectors, tuple):
            for lst in basis_vectors:
                if len(lst) != 2:
                    raise ValueError('Too many arguments per basis vector. '
                                     '{} was provided, but the format should '
                                     'be: vector = ( (\'ID1\' , [x1,y1,z1]), '
                                     '( \'ID2\' , [x2,y2,z2]) )'.format(lst))
                elif not isinstance(lst[0], string_types):
                    raise TypeError('ID Error, not a string. {} is not of '
                                    'type string_types. Please correct in the '
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
                            'not <tuple>. Please review '
                            'the documentation and examples for proper format.'
                            .format(type(basis_vectors)))

        if not isinstance(basis_vectors, defaultdict):
            overlap_dict = defaultdict(list)
            for the_id, vector in basis_vectors:
                num_iter = 3
                for offsets in it.product(range(num_iter), repeat=dimension):
                    offset_vector = tuple((v + offset for v, offset in zip(vector, offsets)))
                    overlap_dict[offset_vector].append((the_id, vector))

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

    def _from_lattice_parameters(self, angles, dimension):
        """Convert Bravais lattice parameters to lattice vectors.

        _from_lattice_parameters will generate the lattice vectors based on
        the parameters necessary to build a Bravais Lattice.

        This was adapted from the ASE triclinic.py lattice parameter code.

        S. R. Bahn and K. W. Jacobsen
        An object-oriented scripting interface to a
        legacy electronic structure code Comput. Sci. Eng., Vol. 4, 56-66, 2002

        Parameters
        ----------
        angles : list-like, required
            Angles of bravais lattice.
        """
        if dimension is 3:
            (alpha, beta, gamma) = angles

            degree = np.pi / 180.0
            cosa = np.cos(alpha*degree)
            cosb = np.cos(beta*degree)
            sinb = np.sin(beta*degree)
            cosg = np.cos(gamma*degree)
            sing = np.sin(gamma*degree)
            lattice_vec = ([1, 0, 0],
                           [cosg, sing, 0],
                           [cosb, (cosa-cosb*cosg)/sing,
                            np.sqrt(sinb**2 - ((cosa - cosb*cosg)/sing)**2)])
        else:
            alpha = angles
            degree = np.pi / 180.0
            cosa = np.cos(alpha*degree)
            sina = np.sin(alpha*degree)
            lattice_vec = ([1, 0], [cosa, sina])

        return lattice_vec

    def populate(self, compound_dict=None, x=None, y=None, z=None):
        """Expand lattice and create compound from lattice.

        populate will expand lattice based on user input. The user must also
        pass in a dictionary that contains the keys that exist in the
        basis_dict. The corresponding Compound will be the full lattice
        returned to the user.

        If no dictionary is passed to the user, Dummy Compounds will be used.

        Parameters
        ----------
        x : int, optional, default=None
            How many iterations in the x direction.
        y : int, optional, default=None
            How many iterations in the y direction.
        z : int, optional, default=None
            How many iterations in the z direction.
        compound_dict : dictionary, optional, default=None
            Link between basis_dict and Compounds.

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

        cell = defaultdict(list)
        for key, val in self.basis_vectors.items():
            for val_item in range(len(val)):
                if self.dimension == 3:
                    for i in range(x):
                        for j in range(y):
                            for k in range(z):
                                tmpx = (val[val_item][0] + i) * a
                                tmpy = (val[val_item][1] + j) * b
                                tmpz = (val[val_item][2] + k) * c
                                tmp_tuple = tuple((tmpx, tmpy, tmpz))
                                cell[key].append(((tmp_tuple)))
                elif self.dimension == 2:
                    for i in range(x):
                        for j in range(y):
                            tmpx = (val[val_item][0] + i) * a
                            tmpy = (val[val_item][1] + j) * b
                            tmp_tuple = tuple((tmpx, tmpy))
                            cell[key].append(((tmp_tuple)))
                else:
                    for i in range(x):
                        tmpx = (val[val_item][0] + i) * a
                        tmp_tuple = tuple((tmpx))
                        cell[key].append(((tmp_tuple)))

        ret_lattice = mb.Compound()
        if compound_dict is None:
            for key_id, all_pos in cell.items():
                particle = mb.Particle(name=key_id, pos=[0, 0, 0])
                for pos in all_pos:
                    particle_to_add = mb.clone(particle)
                    mb.translate(particle_to_add, list(pos))
                    ret_lattice.add(particle_to_add)
        else:
            for key_id, all_pos in cell.items():
                if isinstance(compound_dict[key_id], mb.Compound):
                    compound_to_move = compound_dict[key_id]
                    for pos in all_pos:
                        tmp_comp = mb.clone(compound_to_move)
                        mb.translate(tmp_comp, list(pos))
                        ret_lattice.add(tmp_comp)
                else:
                    err_type = type(compound_dict.get(key_id))
                    TypeError('Invalid type in provided Compound Dictionary. '
                              'For key {}, type: {} was provided, '
                              'not Compound.'.format(key_id, err_type))
        return ret_lattice
