__author__ = 'sallai'

from numpy import *
from numpy.linalg import *

class CoordinateTransform:
    def __init__(self, T=None):
        if(T==None):
            T = eye(4);

        self.T = T
        self.Tinv = inv(T)

    @classmethod
    def translation(cls, P):
        T = eye(4)
        T[0,3] = P[0]
        T[1,3] = P[1]
        T[2,3] = P[2]
        return cls(T)

    @classmethod
    def rotation_around_z(cls, theta):
        T = eye(4)
        T[0,0] = cos(theta)
        T[0,1] = -sin(theta)
        T[1,0] = sin(theta)
        T[1,1] = cos(theta)
        # T[0,3] = 1.0
        # T[1,3] = 1.0
        # T[2,3] = 1.0
        return cls(T)

    @classmethod
    def rotation_around_y(cls, theta):
        T = eye(4)
        T[0,0] = cos(theta)
        T[0,2] = sin(theta)
        T[2,0] = -sin(theta)
        T[2,2] = cos(theta)
        # T[0,3] = 1.0
        # T[1,3] = 1.0
        # T[2,3] = 1.0
        return cls(T)

    @classmethod
    def rotation_around_x(cls, theta):
        T = eye(4)
        T[1,1] = cos(theta)
        T[1,2] = -sin(theta)
        T[2,1] = sin(theta)
        T[2,2] = cos(theta)
        # T[0,3] = 1.0
        # T[1,3] = 1.0
        # T[2,3] = 1.0
        return cls(T)


    @classmethod
    def fromMatrix(cls, T):
        return cls(T)

    @classmethod
    def compute(cls, A, B):
        (rows, cols) = shape(A)

        centroid_A = mean(A, axis=0)
        centroid_B = mean(B, axis=0)
        centroid_A.shape = (1,3)
        centroid_B.shape = (1,3)

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
        C_A = vstack([hstack([C_A, transpose(centroid_A) * -1.0]),array([0,0,0,1])])
        # print "C_A:" + str(C_A)

        R_new = eye(3)
        R_new = vstack([hstack([R,array([[0],[0],[0]])]),array([0,0,0,1])])
        # print "R_new:" + str(R_new)

        C_B = eye(3)

        C_B = vstack([hstack([C_B, transpose(centroid_B)]), array([0,0,0,1])])
        # print "C_B:" + str(C_B)

        T = C_B.dot(R_new).dot(C_A)

        return cls(T)

    def transform(self, A):
        (rows, cols) = A.shape
        A_new = zeros((rows, 4))
        A_new = hstack([A, ones((rows,1))])

        # print "A_new_prime" + str(transpose(A_new))
        # print "T*A_new_prime" + str(self.T.dot(transpose(A_new)))

        A_new = transpose(self.T.dot(transpose(A_new)))
        return A_new[:, 0:cols]


if __name__ == "__main__":
    # matrix with n rows containing x, y, z coordinates(N >= 3) of points in coordinate system 1

    # Test 1
    A1 = array([
        [0.1239, 0.2085, 0.9479],
        [0.4904, 0.5650, 0.0821],
        [0.8530, 0.6403, 0.1057],
        [0.8739, 0.4170, 0.1420],
        [0.2703, 0.2060, 0.1665]])

    # matrix with n rows containing x, y, z coordinates(N >= 3) of the same points in coordinate system 2
    B1 = array([
        [-0.4477, 0.4830, 0.6862],
        [0.2384, -0.1611, 0.3321],
        [0.0970, -0.3850, 0.0721],
        [0.0666, -0.1926, -0.0449],
        [0.2457, 0.2672, 0.3625]])


    # Test 2
    A2 = array([
        [0.0, -0.2, 0.0],
        [-1.0, -0.5, 0.0],
        [1.0, -0.5, 0.0]])

    # matrix with n rows containing x, y, z coordinates(N >= 3) of the same points in coordinate system 2
    B2 = array([
        [0.0, 0.8, 0.0],
        [-1.0, 0.5, 0.0],
        [1.0, 0.5, 0.0]])


    # find out the coordinate transform


    A = A1
    B = B1

    tAB = CoordinateTransform.compute(A,B)

    print "T:" + str(tAB.T)

    # test if it works:
    A2 = tAB.transform(A)

    print "B:" + str(B)
    print "A in B:" + str(A2)
    print "Diff:" + str(A2-B)


    translation = CoordinateTransform.translation((10,10,10))
    print translation.transform(array([[1,1,1]]))

    rotation_around_z = CoordinateTransform.rotation_around_z(pi/4)
    print rotation_around_z.transform(array([[1,1,1]]))
