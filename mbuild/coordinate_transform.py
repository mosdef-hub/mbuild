__author__ = 'sallai'

# from mbuild.tools import createEquivalenceTransform
from itertools import *
from numpy import *
import numpy as np
from numpy.linalg import *


class CoordinateTransform(object):
    def __init__(self, T=None):
        if (T == None):
            T = eye(4);

        self.T = T
        self.Tinv = inv(T)

    def applyTo(self, A):
        """
        Apply the coordinate transformation to points in A
        :param A: list of points (nx3)
        :return: transformed list of points (nx3)
        """
        (rows, cols) = A.shape
        A_new = zeros((rows, 4))
        A_new = hstack([A, ones((rows, 1))])

        A_new = transpose(self.T.dot(transpose(A_new)))
        return A_new[:, 0:cols]

    def applyInverseTo(self, A):
        """
        Apply the inverse coordinate transformation to points in A
        :param A: list of points (nx3)
        :return: transformed list of points (nx3)
        """
        (rows, cols) = A.shape
        A_new = zeros((rows, 4))
        A_new = hstack([A, ones((rows, 1))])

        A_new = transpose(self.Tinv.dot(transpose(A_new)))
        return A_new[:, 0:cols]


class Translation(CoordinateTransform):
    def __init__(self, P):
        T = eye(4)
        T[0, 3] = P[0]
        T[1, 3] = P[1]
        T[2, 3] = P[2]
        super(Translation, self).__init__(T)


class RotationAroundZ(CoordinateTransform):
    def __init__(self, theta):
        T = eye(4)
        T[0, 0] = cos(theta)
        T[0, 1] = -sin(theta)
        T[1, 0] = sin(theta)
        T[1, 1] = cos(theta)
        super(RotationAroundZ, self).__init__(T)


class RotationAroundY(CoordinateTransform):
    def __init__(self, theta):
        T = eye(4)
        T[0, 0] = cos(theta)
        T[0, 2] = sin(theta)
        T[2, 0] = -sin(theta)
        T[2, 2] = cos(theta)
        super(RotationAroundY, self).__init__(T)


class RotationAroundX(CoordinateTransform):
    def __init__(self, theta):
        T = eye(4)
        T[1, 1] = cos(theta)
        T[1, 2] = -sin(theta)
        T[2, 1] = sin(theta)
        T[2, 2] = cos(theta)
        super(RotationAroundX, self).__init__(T)


class Rotation(CoordinateTransform):
    def __init__(self, theta, around):
        # rotation around vector around by angle theta
        assert around.size == 3

        T = eye(4)

        s = sin(theta)
        c = cos(theta)
        t = 1 - c

        n = around / linalg.norm(around)

        x = n[0]
        y = n[1]
        z = n[2]
        m = array([
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c]])
        T[0:3, 0:3] = m
        super(Rotation, self).__init__(T)


class ChangeOfBasis(CoordinateTransform):
    def __init__(self, basis, origin=array([0.0, 0.0, 0.0])):
        assert shape(basis) == (3, 3)

        T = eye(4)

        T[0:3, 0:3] = basis
        T = inv(T)

        T[0:3, 3:4] = -array([origin]).transpose()
        super(ChangeOfBasis, self).__init__(T)


class AxisTransform(CoordinateTransform):
    def __init__(self, new_origin=array([0.0, 0.0, 0.0]),
                 point_on_x_axis=array([1.0, 0.0, 0.0]),
                 point_on_xy_plane=array([1.0, 1.0, 0.0])):
        # change the basis such that p1 is the origin, p2 is on the x axis and p3 is in the xy plane
        p1 = new_origin
        p2 = point_on_x_axis  # positive x axis
        p3 = point_on_xy_plane  # positive y part of the x axis

        # the direction vector of our new x axis
        newx = unit_vector(p2 - p1)
        p3_u = unit_vector(p3 - p1)
        newz = unit_vector(cross(newx, p3_u))
        newy = cross(newz, newx)

        # print "newx=" +str(newx)
        # print "newy=" +str(newy)
        # print "newz=" +str(newz)

        # translation that moves new_origin to the origin
        T_tr = eye(4)
        T_tr[0:3, 3:4] = -array([p1]).transpose()

        # rotation that moves newx to the x axis, newy to the y axis, newz to the z axis
        B = eye(4)
        B[0:3, 0:3] = vstack((newx, newy, newz))

        # the concatentaion of translation and rotation
        B_tr = dot(B, T_tr)

        super(AxisTransform, self).__init__(B_tr)


class RigidTransform(CoordinateTransform):
    """
    Computes the rigid transformation that maps points A to points B
    see http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.173.2196&rep=rep1&type=pdf
    :param A: list of points (nx3) in source coordinate system
    :param B: list of points (nx3) in destination coordinate system
    :return: the CoordinateTransformation object
    """

    def __init__(self, A, B):
        (rows, cols) = shape(A)
        centroid_A = mean(A, axis=0)
        centroid_B = mean(B, axis=0)
        centroid_A.shape = (1, 3)
        centroid_B.shape = (1, 3)

        H = zeros((3, 3), dtype=float)

        for i in range(0, rows):
            H = H + transpose(A[i, :] - centroid_A).dot((B[i, :] - centroid_B))

        U, s, V = svd(H)

        # print "centroid_A:" + str(centroid_A)
        # print "H:" + str(H)

        # print "U:" + str(U)
        V = transpose(V)
        # print "V:" + str(V)

        R = V.dot(transpose(U))
        # print "R:" + str(R)

        C_A = eye(3)
        C_A = vstack(
            [hstack([C_A, transpose(centroid_A) * -1.0]), array([0, 0, 0, 1])])
        # print "C_A:" + str(C_A)

        R_new = eye(3)
        R_new = vstack(
            [hstack([R, array([[0], [0], [0]])]), array([0, 0, 0, 1])])
        # print "R_new:" + str(R_new)

        C_B = eye(3)

        C_B = vstack(
            [hstack([C_B, transpose(centroid_B)]), array([0, 0, 0, 1])])
        # print "C_B:" + str(C_B)

        T = C_B.dot(R_new).dot(C_A)

        super(RigidTransform, self).__init__(T)


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / linalg.norm(vector)


def vec_angle(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = CoordinateTransform.unit_vector(v1)
    v2_u = CoordinateTransform.unit_vector(v2)

    d = dot(v1_u, v2_u)
    if abs(d - 1.0) < 0.000000001:
        angle = 0.0
    else:
        angle = np.arccos(d)
    if isnan(angle):
        if (v1_u == v2_u).all():
            return 0.0
        else:
            return np.pi
    return angle


def _extract_atom_positions(compound):
    arr = np.fromiter(
        chain.from_iterable(atom.pos for atom in compound.atoms()),
        dtype=np.float64)
    arrnx3 = arr.reshape((-1, 3))
    return arrnx3


def _write_back_atom_positions(compound, arrnx3):
    arr = arrnx3.reshape((-1))
    for i, atom in enumerate(compound.atoms()):
        atom.pos = np.array([arr[3 * i], arr[3 * i + 1], arr[3 * i + 2]])


def _createEquivalenceTransform(equiv):
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

    from compound import Compound
    from atom import Atom

    for pair in equiv:
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise Exception('Equivalence pair not a 2-tuple')
        if not ((isinstance(pair[0], Compound) and isinstance(pair[1],
                                                              Compound)) or (
                    isinstance(pair[0], Atom) and isinstance(pair[1], Atom))):
            raise Exception(
                'Equivalence pair type mismatch: pair[0] is a ' + str(
                    type(pair[0])) + ' and pair[1] is a ' + str(
                    type(pair[1])))

        if isinstance(pair[0], Atom):
            self_points = np.vstack([self_points, pair[0].pos])
            other_points = np.vstack([other_points, pair[1].pos])
        if isinstance(pair[0], Compound):
            for atom0 in pair[0].atoms():
                self_points = np.vstack([self_points, atom0.pos])
            for atom1 in pair[1].atoms():
                other_points = np.vstack([other_points, atom1.pos])

    T = RigidTransform(self_points, other_points)
    return T


def equivalence_transform(compound, from_positions=None, to_positions=None,
                          add_bond=True):
    """Computes an affine transformation that maps the from_positions to the respective
    to_positions, and applies this transformation to the compound.
    :param equivalence_pairs: list of equivalence pairs (tuples)
    """

    if isinstance(from_positions, (list, tuple)) and isinstance(to_positions,
                                                                (list, tuple)):
        equivalence_pairs = zip(from_positions, to_positions)
    else:
        equivalence_pairs = [(from_positions, to_positions)]

    T = _createEquivalenceTransform(equivalence_pairs)
    atom_positions = _extract_atom_positions(compound)
    atom_positions = T.applyTo(atom_positions)
    _write_back_atom_positions(compound, atom_positions)

    if add_bond:
        from mbuild.port import Port
        from mbuild.bond import Bond

        if isinstance(from_positions, Port) and isinstance(to_positions, Port):
            compound.add(Bond(from_positions, to_positions))


def translate(compound, v):
    atom_positions = _extract_atom_positions(compound)
    atom_positions = Translation(v).applyTo(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_z(compound, theta):
    atom_positions = _extract_atom_positions(compound)
    atom_positions = RotationAroundZ(theta).applyTo(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_y(compound, theta):
    atom_positions = _extract_atom_positions(compound)
    atom_positions = RotationAroundY(theta).applyTo(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def rotate_around_x(compound, theta):
    atom_positions = _extract_atom_positions(compound)
    atom_positions = RotationAroundX(theta).applyTo(atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def x_axis_transform(compound, new_origin=array([0.0, 0.0, 0.0]),
                     point_on_x_axis=array([1.0, 0.0, 0.0]),
                     point_on_xy_plane=array([1.0, 1.0, 0.0])):
    from atom import Atom

    if isinstance(new_origin, Atom):
        new_origin = new_origin.pos
    if isinstance(point_on_x_axis, Atom):
        point_on_x_axis = point_on_x_axis.pos
    if isinstance(point_on_xy_plane, Atom):
        point_on_xy_plane = point_on_xy_plane.pos

    atom_positions = _extract_atom_positions(compound)
    atom_positions = AxisTransform(new_origin=new_origin,
                                   point_on_x_axis=point_on_x_axis,
                                   point_on_xy_plane=point_on_xy_plane).applyTo(
        atom_positions)
    _write_back_atom_positions(compound, atom_positions)


def y_axis_transform(compound, new_origin=array([0.0, 0.0, 0.0]),
                     point_on_y_axis=array([1.0, 0.0, 0.0]),
                     point_on_xy_plane=array([1.0, 1.0, 0.0])):
    x_axis_transform(compound, new_origin=new_origin,
                     point_on_x_axis=point_on_y_axis,
                     point_on_xy_plane=point_on_xy_plane)
    rotate_around_z(compound, pi / 2)


if __name__ == "__main__":
    # # matrix with n rows containing x, y, z coordinates(N >= 3) of points in coordinate system 1
    #
    # # Test 1
    # A1 = array([
    # [0.1239, 0.2085, 0.9479],
    #     [0.4904, 0.5650, 0.0821],
    #     [0.8530, 0.6403, 0.1057],
    #     [0.8739, 0.4170, 0.1420],
    #     [0.2703, 0.2060, 0.1665]])
    #
    # # matrix with n rows containing x, y, z coordinates(N >= 3) of the same points in coordinate system 2
    # B1 = array([
    #     [-0.4477, 0.4830, 0.6862],
    #     [0.2384, -0.1611, 0.3321],
    #     [0.0970, -0.3850, 0.0721],
    #     [0.0666, -0.1926, -0.0449],
    #     [0.2457, 0.2672, 0.3625]])
    #
    #
    # # Test 2
    # A2 = array([
    #     [0.0, -0.2, 0.0],
    #     [-1.0, -0.5, 0.0],
    #     [1.0, -0.5, 0.0]])
    #
    # # matrix with n rows containing x, y, z coordinates(N >= 3) of the same points in coordinate system 2
    # B2 = array([
    #     [0.0, 0.8, 0.0],
    #     [-1.0, 0.5, 0.0],
    #     [1.0, 0.5, 0.0]])
    #
    #
    # # find out the coordinate transform
    #
    #
    # A = A1
    # B = B1
    #
    # tAB = RigidTransform(A,B)
    #
    # print "T:" + str(tAB.T)
    #
    # # test if it works:
    # A2 = tAB.apply(A)
    #
    # print "B:" + str(B)
    # print "A in B:" + str(A2)
    # print "Diff:" + str(A2-B)
    #
    #
    # translation = Translation((10,10,10))
    # print translation.apply(array([[1,1,1]]))
    #
    # rotation_around_z = RotationAroundZ(pi/4)
    # print rotation_around_z.apply(array([[1,1,1]]))
    #

    # print "Change of basis"
    # CT = ChangeOfBasis(array([[sqrt(3)/2,0.5,0.0],
    #                           [-0.5,sqrt(3)/2,0.0],
    #                           [0.0,0.0,1.0]
    #                         ]))
    #
    # A = array([
    #     [1.0, 0.0, 0.0],
    #     [0.5, sqrt(3)/2, 0.0]
    # ])
    # # A = array([
    # #     [1.0, 0.0, 0.0],
    # #     [0.0, 1.0, 0.0],
    # #     [0.0, 0.0, 1.0]])
    #
    # A_prime = CT.apply(A)
    #
    # print "A_prime=" + str(A_prime)


    #
    #
    # print "Axis Transform"
    # new_origin=array([1,0,0])
    # point_on_x_axis=new_origin + array([.5, 0.0+sqrt(3)/2, 0.0])
    # CT = AxisTransform(new_origin=new_origin, point_on_x_axis=point_on_x_axis, point_on_xy_plane=new_origin+array([1.0,1.0,0.0]))
    #
    # A = array([
    #     new_origin,
    #     point_on_x_axis,
    #     [1.0, 0.0, 0.0],
    #     [0.5, sqrt(3)/2, 0.0]
    # ])
    # # A = array([
    # #     [1.0, 0.0, 0.0],
    # #     [0.0, 1.0, 0.0],
    # #     [0.0, 0.0, 1.0]])
    #
    # A_prime = CT.apply(A)
    #
    # print "A=" + str(A)
    # print "A_prime=" + str(A_prime)


    # print "Axis Transform"
    # CT = AxisTransform(array([0.0,0.0,0.0]), array([1.0, 0.0, 0.0]), array([0.0,1.0,0.0]))
    #
    # A = array([
    #     [2.0, 0.0, 0.0],
    #     [0.5, sqrt(3)/2, 1.0]
    # ])
    # # A = array([
    # #     [1.0, 0.0, 0.0],
    # #     [0.0, 1.0, 0.0],
    # #     [0.0, 0.0, 1.0]])
    #
    # A_prime = CT.apply(A)
    #
    # print "A=" + str(A)
    # print "A_prime=" + str(A_prime)


    # B = array([
    #     [0.5, sqrt(3)/2, 0.0, 0.0],
    #     [-sqrt(3)/2, 0.5, 0.0, 0.0],
    #     [0.0, 0.0, 1.0, 0.0],
    #     [0.0, 0.0, 0.0, 1.0]])
    #
    # Binv = inv(B)
    #
    # v1 = array([[1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0]])
    #
    # v1_prime = dot(Binv,v1.transpose()).transpose()
    # print "the image of v1="+str(v1)+" is v1_prime=" + str(v1_prime)
    #
    #
    # v_tr = array([2.0, 3.0, 4.0])
    # # v_tr = array([0.0, 0.0, 0.0])
    # Binv_tr = Binv
    # Binv_tr[0:3, 3:4] = array([v_tr]).transpose()
    # print "Binv_tr="+str(Binv_tr)
    #
    # v1_prime_tr = dot(Binv_tr, v1.transpose()).transpose()
    # print "the translated image of v1="+str(v1)+" is v1_prime_tr=" + str(v1_prime_tr)
    #
    # Basis = B[0:3,0:3]
    # T = ChangeOfBasis(Basis, -v_tr)
    # v1_prime_tr_T = T.apply(array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]]))
    # print "the translated image of v1="+str(v1)+" is v1_prime_tr_T=" + str(v1_prime_tr_T)
    #
    # AT = AxisTransform(new_origin=-v_tr, point_on_x_axis=B[0,0:3]-v_tr, point_on_xy_plane=B[1,0:3]-v_tr)
    # v1_prime_tr_AT = AT.apply(array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]]))
    # print "the translated image of v1="+str(v1)+" is v1_prime_tr_AT=" + str(v1_prime_tr_AT)

    v_tr = array([2.0, 3.0, 4.0])

    # translate by -v_tr
    T_tr = eye(4)
    T_tr[0:3, 3:4] = -array([v_tr]).transpose()

    # rotate 60 degrees
    B = array([
        [0.5, sqrt(3) / 2, 0.0, 0.0],
        [-sqrt(3) / 2, 0.5, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]])

    B_tr = dot(B, T_tr)

    v1 = array([[2.0, 3.0, 4.0, 1.0], [2.0 + 0.5, 3.0 + sqrt(3) / 2, 4.0, 1.0],
                [1.0, 0.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0]])

    v1_prime = dot(B_tr, v1.transpose()).transpose()
    print "the image of v1=" + str(v1) + " is v1_prime=" + str(v1_prime)

    AT = AxisTransform(new_origin=v_tr,
                       point_on_x_axis=[2.0 + 0.5, 3.0 + sqrt(3) / 2, 4.0],
                       point_on_xy_plane=[2.0, 4.0, 4.0])
    v1_prime_tr_AT = AT.apply(v1[:, 0:3])
    print "the image of v1=" + str(v1[:, 0:3]) + " is v1_prime_tr_AT=" + str(
        v1_prime_tr_AT)
