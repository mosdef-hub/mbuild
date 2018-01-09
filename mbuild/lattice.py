from collections import defaultdict
import itertools as it
from warnings import warn

import numpy as np


import mbuild as mb

__all__ = ['Lattice']


class Lattice(object):
    """Develop crystal structure from user defined inputs.

    Lattice, the abstract building block of a crystal cell.
    Once defined by the user, the lattice can then be populated with
    Compounds and replicated as many cell lengths desired in 3D space.

    A Lattice is defined through the Bravais lattice definitions. With edge
    vectors a1, a2, a3; lattice spacing a,b,c; and lattice points at unique
    fractional positions between 0-1 in 3 dimensions. This encapsulates
    distance, area, volume, depending on the parameters defined.


    Parameters
    ----------
    lattice_spacing : array-like, shape=(3,), required, dtype=float
        Array of lattice spacings a,b,c for the cell.
    lattice_vectors : array-like, shape=(3, 3), optional
                      default=[[1,0,0], [0,1,0], [0,0,1]]
        Vectors that encase the unit cell corresponding to dimension. Will
        only default to these values if no angles were defined as well.
    lattice_points : dictionary, shape={'id': [[nested list of positions]]
        optional, default={'default': [[0.,0.,0.]]}
        Locations of all lattice points in cell using fractional coordinates.
    angles : array-like, shape=(3,), optional, dtype=float
        Array of inter-planar Bravais angles in degrees.

    Attributes
    ----------
    dimension : int, 3
        Default dimensionality within mBuild. If choosing a lower dimension,
        pad the relevant arrays with zeroes.
    lattice_spacing : numpy array, shape=(3,), required, dtype=float
        Array of lattice spacings a,b,c for the cell.
    lattice_vectors : numpy array, shape=(3, 3), optional
                      default=[[1,0,0], [0,1,0], [0,0,1]]
        Vectors that encase the unit cell corresponding to dimension. Will
        only default to these values if no angles were defined as well.
    lattice_points : dictionary, shape={'id': [[nested list of positions]]
        optional, default={'default': [[0.,0.,0.]]}
        Locations of all lattice points in cell using fractional coordinates.
    angles : numpy array, shape=(3,), optional, dtype=float
        Array of inter-planar Bravais angles

    Examples
    --------
    Generating a triclinic lattice for cholesterol.

    >>> import mbuild as mb
    >>> from mbuild.utils.io import get_fn
    >>> # reading in the lattice parameters for crystalline cholesterol
    >>> angle_values = [94.64, 90.67, 96.32]
    >>> spacing = [1.4172, 3.4209, 1.0481]
    >>> basis = {'cholesterol':[[0., 0., 0.]]}
    >>> cholesterol_lattice = mb.Lattice(spacing,
    ...                                  angles=angle_values,
    ...                                  lattice_points=basis)

    >>> # The lattice based on the bravais lattice parameters of crystalline
    >>> # cholesterol was generated.

    >>> # Replicating the triclinic unit cell out 3 replications
    >>> # in x,y,z directions.

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
    >>> spacing = [.4123, .4123, .4123]
    >>> basis = {'Cl' : [[0., 0., 0.]], 'Cs' : [[.5, .5, .5]]}
    >>> cscl_lattice = mb.Lattice(spacing, lattice_points=basis)

    >>> # Now associate id with Compounds for lattice points and replicate 3x

    >>> cscl_dict = {'Cl' : chlorine, 'Cs' : cesium}
    >>> cscl_compound = cscl_lattice.populate(x=3, y=3, z=3,
    ...                                       compound_dict=cscl_dict)

    A multi-Compound basis was created and replicated. For each unique basis
    atom position, a separate entry must be completed for the basis_atom
    input.

    Generating FCC Copper cell with lattice_vectors instead of angles

    >>> import mbuild as mb
    >>> copper = mb.Compound(name='Cu')
    >>> lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    >>> spacing = [.36149, .36149, .36149]
    >>> copper_locations = [[0., 0., 0.], [.5, .5, 0.],
    ...                     [.5, 0., .5], [0., .5, .5]]
    >>> basis = {'Cu' : copper_locations}
    >>> copper_lattice = mb.Lattice(lattice_spacing = spacing,
    ...                             lattice_vectors=lattice_vector,
    ...                             lattice_points=basis)
    >>> copper_dict = {'Cu' : copper}
    >>> copper_pillar = copper_lattice.populate(x=3, y=3, z=20,
    ...                                       compound_dict=copper_dict)

    Generating the 2d Structure Graphene carbon backbone

    >>> import mbuild as mb
    >>> carbon = mb.Compound(name='C')
    >>> angles = [90, 90, 120]
    >>> carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
    >>> basis = {'C' : carbon_locations}
    >>> graphene = mb.Lattice(lattice_spacing=[.2456, .2456, 0],
    ...                        angles=angles, lattice_points=basis)
    >>> carbon_dict = {'C' : carbon}
    >>> graphene_cell = graphene.populate(compound_dict=carbon_dict,
    ...                                   x=3, y=3, z=1)

    """

    def __init__(self, lattice_spacing=None, lattice_vectors=None,
                 lattice_points=None, angles=None):
        super(Lattice, self).__init__()
        self.dimension = 3
        self.lattice_spacing = None
        self.lattice_vectors = None
        self.lattice_points = dict()
        self.angles = None
        self._sanitize_inputs(lattice_spacing=lattice_spacing,
                              lattice_vectors=lattice_vectors,
                              lattice_points=lattice_points,
                              angles=angles)

    def _sanitize_inputs(self, lattice_spacing, lattice_vectors,
                         lattice_points, angles):
        """Check for proper inputs and set instance attributes.

        validate_inputs takes the data passed to the constructor by the user
        and will ensure that the data is correctly formatted and will then
        set its instance attributes.

        validate_inputs checks that dimensionality is maintained,
        the unit cell is right handed, the area or volume of the unit cell
        is positive and non-zero for 2D and 3D respectively, lattice spacing
        are provided, basis vectors do not overlap when the unit cell is
        expanded.

        Exceptions Raised
        -----------------
        TypeError : incorrect typing of the input parameters.

        ValueError : values are not within restrictions.
        """

        if angles is not None and lattice_vectors is not None:
            raise ValueError('Overdefined system: angles and lattice_vectors '
                             'provided. Only one of these should be passed.')

        self._validate_lattice_spacing(lattice_spacing)

        if angles is not None:
            self._validate_angles(angles)
            self.lattice_vectors = self._from_lattice_parameters(self.angles)
        else:
            self._validate_lattice_vectors(lattice_vectors)
            self.angles = self._from_lattice_vectors()

        self._validate_lattice_points(lattice_points)

    def _validate_lattice_spacing(self, lattice_spacing):
        """Ensure that lattice spacing is provided and correct.

        _validate_lattice_spacing will ensure that the lattice spacing
        provided are acceptable values. Additional Numpy errors can also occur
        due to the conversion to a Numpy array.

        Exceptions Raised
        -----------------
        ValueError : Incorrect lattice_spacing input
        """

        dataType = np.float64

        if lattice_spacing is not None:
            lattice_spacing = np.asarray(lattice_spacing, dtype=dataType)
            lattice_spacing = lattice_spacing.reshape((3,))
            if np.shape(lattice_spacing) != (self.dimension,):
                raise ValueError('Lattice spacing should be a vector of '
                                 'size:({},). Please include lattice spacing '
                                 'of size >= 0 depending on desired '
                                 'dimensionality.'
                                 .format(self.dimension))
        else:
            raise ValueError('No lattice_spacing provided. Please provide '
                             'lattice spacing\'s that are >= 0. with size ({},)'
                             .format((self.dimension)))

        if np.any(np.isnan(lattice_spacing)):
            raise ValueError('None type or NaN type values present in '
                             'lattice_spacing: {}.'.format(lattice_spacing))
        elif np.any(lattice_spacing < 0.0):
            raise ValueError('Negative lattice spacing value. One of '
                             'the spacing: {} is negative.'
                             .format(lattice_spacing))

        self.lattice_spacing = lattice_spacing

    def _validate_angles(self, angles):
        """Ensure that the angles between the lattice_vectors are correct"""

        dataType = np.float64
        tempAngles = np.asarray(angles, dtype=dataType)
        tempAngles = tempAngles.reshape((3,))

        if np.shape(tempAngles) == (self.dimension,):
            if np.sum(tempAngles) < 360.0 or np.sum(tempAngles) > -360.0:
                if (np.all(tempAngles != 180.0)
                        and np.all(tempAngles != 0.0)):
                    pass
                else:
                    raise ValueError('Angles cannot be 180.0 or 0.0')
            else:
                raise ValueError('Angles sum: {} is either greater than '
                                 '360.0 or less than -360.0'
                                 .format(np.sum(tempAngles)))

            for subset in it.permutations(tempAngles, r=self.dimension):
                if not subset[0] < np.sum(tempAngles) - subset[0]:
                    raise ValueError('Each angle provided must be less'
                                     'than the sum of the other angles. '
                                     '{} is greater.'.format(subset[0]))
        else:
            raise ValueError('Incorrect array size. When converted to a '
                             'Numpy array, the shape is: {}, expected {}.'
                             .format(np.shape(tempAngles),
                                     (3,)))
        self.angles = tempAngles

    def _validate_lattice_vectors(self, lattice_vectors):
        """Ensure that the lattice_vectors are reasonable inputs.

        """
        dataType = np.float64
        if lattice_vectors is None:
                lattice_vectors = np.identity(self.dimension, dtype=dataType)
        else:
            lattice_vectors = np.asarray(lattice_vectors, dtype=dataType)

            if (self.dimension, self.dimension) != np.shape(lattice_vectors):
                raise ValueError('Dimensionality of lattice_vectors is '
                                 ' of shape {} not {}.'
                                 .format(np.shape(lattice_vectors),
                                         (self.dimension, self.dimension)))

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

    def _validate_lattice_points(self, lattice_points):
        if lattice_points is None:
            lattice_points = {}
            lattice_points = {'id': [[0. for x in range(self.dimension)]]}
        elif isinstance(lattice_points, dict):
            pass
        else:
            raise TypeError('Incorrect type, lattice_points is of type {}, '
                            'Expected dict.'.format(type(lattice_points)))

        for name, positions in lattice_points.items():
            for pos in positions:
                if len(pos) != self.dimension:
                    raise ValueError("Incorrect lattice point position size. "
                                     "lattice point {} has location "
                                     "{}, which is inconsistent with the "
                                     "dimension {}.".format(name, pos,
                                                            self.dimension))
                if pos is None:
                    raise ValueError("NoneType passed, expected float. "
                                     "None was passed in as position for {}."
                                     .format(name))
                for coord in pos:
                    if (coord is None) or (0 > coord) or (coord >=1):
                        raise ValueError('Incorrect lattice point fractional '
                                         'coordinates. Coordinates cannot be '
                                         '{}, {}, or {}. You passed {}.'
                                         .format('None', '>= 1', '< 0', coord))

        self.lattice_points = self._check_for_overlap(lattice_points)

    def _check_for_overlap(self, lattice_points):

        overlap_dict = defaultdict(list)
        num_iter = 3
        dim = self.dimension
        for name, positions in lattice_points.items():
            for pos in positions:
                for offsets in it.product(range(num_iter), repeat=dim):
                    offset_vector = tuple((v + offset for v, offset in zip(pos, offsets)))
                    overlap_dict[offset_vector].append((pos))

        for key, val in overlap_dict.items():
            if len(val) > 1:
                raise ValueError('Overlapping lattice points: Lattice '
                                 'points overlap when the unit cell is '
                                 'expanded to {}. This is an incorrect '
                                 'perfect lattice. The offending '
                                 'points are: {}'
                                 .format(key, val))
        return lattice_points

    def _from_lattice_parameters(self, angles):
        """Convert Bravais lattice parameters to lattice vectors.

        _from_lattice_parameters will generate the lattice vectors based on
        the parameters necessary to build a Bravais Lattice. The lattice
        vectors are in the lower diagonal matrix form.

        This was adapted from the ASE triclinic.py lattice parameter code.

        S. R. Bahn and K. W. Jacobsen
        An object-oriented scripting interface to a
        legacy electronic structure code Comput. Sci. Eng., Vol. 4, 56-66, 2002

        Parameters
        ----------
        angles : list-like, required
            Angles of bravais lattice.
        """

        dataType = np.float64
        (alpha, beta, gamma) = angles

        radianConversion = np.pi / 180.0
        cosa = np.cos(alpha * radianConversion)
        cosb = np.cos(beta * radianConversion)
        sinb = np.sin(beta * radianConversion)
        cosg = np.cos(gamma * radianConversion)
        sing = np.sin(gamma * radianConversion)
        matCoef_y = (cosa - cosb * cosg) / sing
        matCoef_z = np.power(sinb, 2, dtype=dataType) - \
            np.power(matCoef_y, 2, dtype=dataType)

        if matCoef_z > 0.:
            matCoef_z = np.sqrt(matCoef_z)
        else:
            raise ValueError('Incorrect lattice vector coefficients.'
                             'Lattice parameters chosen return a non-positive '
                             'z vector.')

        lattice_vec = [[1, 0, 0],
                       [cosg, sing, 0],
                       [cosb, matCoef_y, matCoef_z]]

        return np.asarray(lattice_vec, dtype=np.float64)

    def _from_lattice_vectors(self):
        """Calculate the angles between the vectors that define the lattice.

        _from_lattice_vectors will calculate the angles alpha, beta, and
        gamma from the Lattice object attribute lattice_vectors.
        """

        degreeConvsersion = 180.0 / np.pi
        vector_magnitudes = np.linalg.norm(self.lattice_vectors, axis=1)

        a_dot_b = np.dot(self.lattice_vectors[0], self.lattice_vectors[1])
        b_dot_c = np.dot(self.lattice_vectors[1], self.lattice_vectors[2])
        a_dot_c = np.dot(self.lattice_vectors[0], self.lattice_vectors[2])

        alpha_raw = a_dot_c / (vector_magnitudes[0] * vector_magnitudes[2])
        beta_raw = b_dot_c / (vector_magnitudes[1] * vector_magnitudes[2])
        gamma_raw = a_dot_b / (vector_magnitudes[0] * vector_magnitudes[1])

        alpha = np.arccos(np.clip(alpha_raw, -1.0, 1.0)) * degreeConvsersion
        beta = np.arccos(np.clip(beta_raw, -1.0, 1.0)) * degreeConvsersion
        gamma = np.arccos(np.clip(gamma_raw, -1.0, 1.0)) * degreeConvsersion

        return np.asarray([alpha, beta, gamma], dtype=np.float64)

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
        error_dict = {0: 'X', 1: 'Y', 2: 'Z'}
        try:
            x = int(x)
            y = int(y)
            z = int(z)
        except (ValueError, TypeError):
            raise ValueError('Cannot convert replication amounts into '
                             'integers. x= {}, y= {}, z= {} needs to '
                             'be an int.'.format(x, y, z))

        for replication_amount in x, y, z:
            if replication_amount is None:
                raise ValueError('Attempt to replicate None times. '
                                 'None is not an acceptable replication '
                                 'amount, 1 is the default.')

        for replication_amount, index in zip([x, y, z], range(3)):
            if replication_amount < 1:
                raise ValueError('Incorrect populate value: {} : {} is < 1. '
                                 .format(error_dict[index],
                                         replication_amount))

        if ((isinstance(compound_dict, dict)) or (compound_dict is None)):
            pass
        else:
            raise TypeError('Compound dictionary is not of type dict. '
                            '{} was passed.'.format(type(compound_dict)))

        cell = defaultdict(list)
        [a, b, c] = self.lattice_spacing

        transform_mat = self.lattice_vectors
        # unit vectors
        transform_mat = np.asarray(transform_mat, dtype=np.float64)
        transform_mat = np.reshape(transform_mat, newshape=(3,3))
        norms = np.linalg.norm(transform_mat, axis=1)

        # normalized vectors for change of basis
        unit_vecs = np.divide(transform_mat.transpose(), norms)

        for key, locations in self.lattice_points.items():
            for coords in locations:
                for replication in it.product(range(x), range(y), range(z)):
                    temp_location = list()

                    new_coords = np.asarray(coords, dtype=np.float64)
                    new_coords = np.reshape(new_coords, (1, 3), order='C')

                    new_coords[0][0] = new_coords[0][0] + replication[0]
                    new_coords[0][1] = new_coords[0][1] + replication[1]
                    new_coords[0][2] = new_coords[0][2] + replication[2]

                    # change of basis to cartesian
                    new_coords = np.dot(unit_vecs, new_coords.transpose())

                    new_coords[0] = new_coords[0] * a
                    new_coords[1] = new_coords[1] * b
                    new_coords[2] = new_coords[2] * c
                    new_coords = np.reshape(new_coords, (1, 3), order='C')

                    tuple_of_coords = tuple(new_coords.flatten())
                    cell[key].append(tuple_of_coords)

        ret_lattice = mb.Compound()

        if compound_dict is None:
            for key_id, all_pos in cell.items():
                particle = mb.Compound(name=key_id, pos=[0, 0, 0])
                for pos in all_pos:
                    particle_to_add = mb.clone(particle)
                    particle_to_add.translate_to(list(pos))
                    ret_lattice.add(particle_to_add)
        else:
            for key_id, all_pos in cell.items():
                if isinstance(compound_dict[key_id], mb.Compound):
                    compound_to_move = compound_dict[key_id]
                    for pos in all_pos:
                        tmp_comp = mb.clone(compound_to_move)
                        tmp_comp.translate_to(list(pos))
                        ret_lattice.add(tmp_comp)
                else:
                    err_type = type(compound_dict.get(key_id))
                    raise TypeError('Invalid type in provided Compound '
                                    'dictionary. For key {}, type: {} was '
                                    'provided, not mbuild.Compound.'
                                    .format(key_id, err_type))
        # set periodicity
        ret_lattice.periodicity = np.asarray([a * x, b * y, c * z], dtype=np.float64)
        warn('Periodicity of non-rectangular lattices are not valid with '
                    'default boxes. Only rectangular lattices are valid '
                    'at this time.')

        # if coordinates are below a certain threshold, set to 0
        tolerance = 1e-12
        ret_lattice.xyz_with_ports[ret_lattice.xyz_with_ports <= tolerance] = 0.

        return ret_lattice
