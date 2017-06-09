from collections import defaultdict
from copy import deepcopy
import itertools as it

import numpy as np
from six import string_types
from six.moves import zip_longest

import mbuild as mb

__all__ = ['Lattice']


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
    basis_atoms : dictionary, shape={'id':[nested list of coordinate pairs]}
                    default={'default':[[0., 0., 0.]]
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
    basis_atoms : list-like, shape=(['id',[dimension,]], ... ,) optional
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
    >>> basis = {'cholesterol':[[0., 0., 0.]]}
    >>> cholesterol_lattice = mb.Lattice(spacings,
    ...                                  angles=angle_values,
    ...                                  basis_atoms=basis,
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
    >>> basis = {'Cl' : [[0., 0., 0.]], 'Cs' : [[.5, .5, .5]]}
    >>> cscl_lattice = mb.Lattice(spacings, basis_atoms=basis,
    ...                           dimension=3)

    Now associate id with Compounds for basis atoms and replicate 3x3x3
    >>> cscl_dict = {'Cl' : chlorine, 'Cs' : cesium}
    >>> cscl_compound = cscl_lattice.populate(x=3, y=3, z=3,
    ...                                       compound_dict=cscl_dict)

    A multi-Compound basis was created and replicated. For each unique basis
    atom position, a separate entry must be completed for the basis_atom
    input.

    Generating FCC Copper cell with lattice_vectors instead of angles
    >>> import mbuild as mb
    >>> copper = mb.Compound(name='Cu')
    >>> lattice_vector = ( [1, 0, 0], [0, 1, 0], [0, 0, 1])
    >>> spacings = [.36149, .36149, .36149]
    >>> copper_locations = [[0., 0., 0.], [.5, .5, 0.],
    ...                     [.5, 0., .5], [0., .5, .5]]
    >>> basis = {'Cu' : copper_locations}
    >>> copper_lattice = mb.Lattice(spacings, dimension=3,
    ...                           lattice_vectors=lattice_vector,
    ...                           basis_atoms=basis)
    >>> copper_dict = {'Cu' : copper}
    >>> copper_cell = copper_lattice.populate(x=3, y=3, z=20,
    ...                                       compound_dict=copper_dict)

    TODO(Justin Gilmer) : Print function to display info about Lattice (repr)
    TODO(Justin Gilmer) : inheritance(Cubic, orthorhombic, hexangonal)
    TODO(Justin Gilmer) : orientation functionality
    """

    def __init__(self, lattice_spacings, dimension=None,
                 lattice_vectors=None, basis_atoms=None,
                 angles=None):
        super(Lattice, self).__init__()
        self.lattice_spacings = None
        self.dimension = None
        self.lattice_vectors = None
        self.basis_atoms = dict()
        self.angles = None
        self._sanitize_inputs(lattice_vectors=lattice_vectors,
                              dimension=dimension,
                              lattice_spacings=lattice_spacings,
                              basis_atoms=basis_atoms,
                              angles=angles)

    def _sanitize_inputs(self, lattice_vectors, dimension,
                         lattice_spacings, basis_atoms, angles):
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

        self._validate_dimension(dimension)
        self._validate_lattice_spacing(lattice_spacings, self.dimension)

        if angles is not None and lattice_vectors is not None:
            raise ValueError('Overdefined system: angles and lattice_vectors '
                             'provided. Only one of these should be passed.')
        if angles is not None:
            self._validate_angles(angles, self.dimension)
            self.lattice_vectors = self._from_lattice_parameters(
                self.angles, self.dimension)
        else:
            self._validate_lattice_vectors(lattice_vectors, self.dimension)

        self._validate_basis_atoms(basis_atoms, self.dimension)

    def _validate_dimension(self, dimension):
        """Ensure that dimension input is correct.

        _validate_dimension will check for that the dimensionality
        passed to the constructor is a proper input.

        If the dimensionality is None, the default value is 3,
        or the user can specify 1D or 2D.

        If _validate_dimension cannot convert the passed in value to an int,
        or if the dimension is <1 or >3, a ValueError will be raised.

        Exceptions Raised
        -----------------
        ValueError : Incorrect typing of the input parameter.
        """
        if dimension is None:
            dimension = 3
        else:
            dimension = int(dimension)
        if dimension < 1 or dimension > 3:
            raise ValueError('Incorrect dimensions: {} is not a proper '
                             'dimension. 1, 2, or 3 are acceptable.'
                             .format(dimension))
        self.dimension = dimension

    def _validate_lattice_spacing(self, lattice_spacings, dimension):
        """Ensure that lattice spacing is provided and correct.

        _validate_lattice_spacing will ensure that the lattice spacings
        provided are acceptable values and dimensionally constant.

        Exceptions Raised
        -----------------
        ValueError : Incorrect lattice_vectors input
        """
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
        self.lattice_spacings = lattice_spacings

    def _validate_angles(self, angles, dimension):
        if angles is not None:
            for index, value in enumerate(angles):
                angles[index] = float(value)
            if (len(angles), dimension) == (3, 3):
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

                for subset in it.permutations(angles, 3):
                    if not subset[0] < sum(angles) - subset[0]:
                        raise ValueError('Each angle provided must be less '
                                         'than the sum of the other two '
                                         'angles. {} is greater.'
                                         .format(subset[0]))
                self.angles = angles

            elif len(angles) == 1 and dimension == 2:
                for theAngle in angles:
                    if (theAngle != 180.0 and theAngle != 0.0 and
                            theAngle < 180.0 and theAngle > -180.0):
                        pass
                    else:
                        raise ValueError('Angle incorrectly defined. {} '
                                         'does not follow the proper '
                                         'guidelines for a bravais angle. '
                                         .format(theAngle))
                self.angles = angles
            else:
                raise ValueError('Incorrect amount of angles provided for '
                                 'dimension {}. Recieved {} angles.'
                                 .format(dimension, len(angles)))

    def _validate_lattice_vectors(self, lattice_vectors, dimension):
        """Ensure that the lattice_vectors are reasonable inputs.

        """
        if lattice_vectors is None:
                lattice_vectors = np.identity(dimension, dtype=float)
        else:
            lattice_vectors = np.asarray(lattice_vectors, dtype=float)
            shape = np.shape(lattice_vectors)

            if (dimension, dimension) != shape:
                raise ValueError('Dimensionality of lattice_vectors is '
                                 ' of shape {} not {}.'
                                 .format(shape, (dimension, dimension)))
            if dimension > 1:
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
        self.lattice_vectors = lattice_vectors

    def _validate_basis_atoms(self, basis_atoms, dimension):
        if basis_atoms is None:
            basis_atoms = {}
            basis_atoms = {'default': [[0. for x in range(dimension)]]}
        elif isinstance(basis_atoms, dict):
            pass
        else:
            raise TypeError('Incorrect type, basis_atoms is of type {}, '
                            'Expected dict.'.format(type(basis_atoms)))

        for name in basis_atoms.keys():
            positions = basis_atoms[name]
            for pos in positions:
                location_check = []
                if len(pos) != dimension:
                    raise ValueError("Incorrect basis atom position size. "
                                     "Basis atom {} was passed with location "
                                     "{}, which is inconsistent with the "
                                     "dimension {}.".format(name, pos,
                                                            dimension))
                if pos is None:
                    raise ValueError("NoneType passed, expected float. "
                                     "None was passed in as position for {}."
                                     .format(name))

                location_check = [coord for coord in pos if coord is None or coord >= 1. or coord < 0.]
                if len(location_check) != 0:
                    raise ValueError("Incorrect coordinate value for basis. "
                                     "Basis {}, was passed coordinates {}. "
                                     "The coordinates {}, were either < 0, or"
                                     " > 1.".format(name, pos, location_check))

        self.basis_atoms = self._check_for_overlap(basis_atoms, dimension)

    def _check_for_overlap(self, basis_atoms, dimension):

        overlap_dict = defaultdict(list)
        num_iter = 3
        for name in basis_atoms.keys():
            positions = basis_atoms[name]
            for pos in positions:
                for offsets in it.product(range(num_iter), repeat=dimension):
                    offset_vector = tuple((v + offset for v, offset in zip(pos, offsets)))
                    overlap_dict[offset_vector].append((pos))

        for key, val in overlap_dict.items():
            if len(val) > 1:
                raise ValueError('Overlapping Basis Vectors: Basis '
                                 'vectors overlap when the unit cell is '
                                 'expanded to {}. This is an incorrect '
                                 'perfect lattice. The offending '
                                 'vectors are: {}'
                                 .format(key, val))
        return basis_atoms

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
        dimension : integer, required
            Dimensionality of system, can only be 2 or 3.
        """
        if dimension is 3:
            (alpha, beta, gamma) = angles

            degree = np.pi / 180.0
            cosa = np.cos(alpha * degree)
            cosb = np.cos(beta * degree)
            sinb = np.sin(beta * degree)
            cosg = np.cos(gamma * degree)
            sing = np.sin(gamma * degree)
            lattice_vec = ([1, 0, 0],
                           [cosg, sing, 0],
                           [cosb, (cosa - cosb * cosg) / sing,
                            np.sqrt(sinb**2 - ((cosa - cosb * cosg) / sing)**2)])
        else:
            alpha = angles
            degree = np.pi / 180.0
            cosa = np.cos(alpha * degree)
            sina = np.sin(alpha * degree)
            lattice_vec = ([1, 0], [cosa, sina])

        return lattice_vec

    def populate(self, compound_dict=None, x=1, y=1, z=1):
        """Expand lattice and create compound from lattice.

        populate will expand lattice based on user input. The user must also
        pass in a dictionary that contains the keys that exist in the
        basis_dict. The corresponding Compound will be the full lattice
        returned to the user.

        If no dictionary is passed to the user, Dummy Compounds will be used.

        Parameters
        ----------
        x : int, optional, default=1
            How many iterations in the x direction.
        y : int, optional, default=1
            How many iterations in the y direction.
        z : int, optional, default=1
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
        error_dict = {0:'X', 1:'Y', 2:'Z'}

        # padded for Compound compatibility
        cell_edges = [edge[0] for edge in zip_longest(self.lattice_spacings, range(3), fillvalue=0.0)]

        for replication_amount in x, y, z:
            if replication_amount is None:
                raise ValueError('Attempt to replicate None times. '
                                 'None is not an acceptable replication amount, '
                                 '1 is the default.')

        for replication_amount, index in zip([x, y, z], range(3)):
            if replication_amount < 1:
                raise ValueError('Incorrect populate value: {} : {} is < 1. '
                                 .format(error_dict[index], replication_amount))

        if self.dimension == 2:
            if z > 1:
                raise ValueError('Attempting to replicate in Z. '
                                 'A non-default value for Z is being '
                                 'passed. 1 is the default value, not {}.'
                                 .format(z))
        elif self.dimension == 1:
            if (y > 1) or (z > 1):
                raise ValueError('Attempting to replicate in Y or Z. '
                                 'A non-default value for Y or Z is being '
                                 'passed. 1 is the default value.')
        else:
            pass

        if ((isinstance(compound_dict, dict)) or (compound_dict is None)):
            pass
        else:
            raise TypeError('Compound dictionary is not of type dict. '
                            '{} was passed.'.format(type(compound_dict)))

        cell = defaultdict(list)
        [a, b, c] = cell_edges
        for key, locations in self.basis_atoms.items():
            for coords in range(len(locations)):
                for replication in it.product(range(x), range(y), range(z)):
                    tmpx = (locations[coords][0] + replication[0]) * a

                    try:
                        tmpy = (locations[coords][1] + replication[1]) * b
                    except IndexError:
                        tmpy = 0.0

                    try:
                        tmpz = (locations[coords][2] + replication[2]) * c
                    except IndexError:
                        tmpz = 0.0

                    tmp_tuple = tuple((tmpx, tmpy, tmpz))
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
                    raise TypeError('Invalid type in provided Compound dictionary. '
                              'For key {}, type: {} was provided, '
                              'not mbuild.Compound.'.format(key_id, err_type))
        return ret_lattice
