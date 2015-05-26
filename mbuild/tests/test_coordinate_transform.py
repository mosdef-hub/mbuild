import numpy as np

import mbuild as mb
from mbuild.coordinate_transform import Translation, CoordinateTransform, RotationAroundZ, RotationAroundY, RotationAroundX, Rotation, ChangeOfBasis, AxisTransform, RigidTransform
from mbuild.tests.base_test import BaseTest


class TestCoordinateTransform(BaseTest):

    def test_apply_to(self):
        double = CoordinateTransform(T=np.eye(4)*2)
        A = np.array([[1,2,3],[4,5,6],[7,8,9]])
        assert (double.apply_to(A) == np.array([[2,4,6],[8,10,12],[14,16,18]])).all()

    def test_translation(self):
        translation = Translation((10,10,10))
        assert (translation.apply_to(np.array([[1,1,1]])) == np.array([11,11,11])).all()

    def test_rotation_around_z(self):
        Z_rotation = RotationAroundZ(np.pi)
        a = Z_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([[-2.,-3.,4.]])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation_around_y(self):
        Y_rotation = RotationAroundY(np.pi)
        a = Y_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([-2,3,-4])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation_around_x(self):
        X_rotation = RotationAroundX(np.pi)
        a = X_rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([2,-3,-4])
        assert np.allclose(a, b, atol=1.e-16)

    def test_rotation(self):
        rotation = Rotation(np.pi*2/3, np.array([1,1,1]))
        a = rotation.apply_to(np.array([[2,3,4]]))
        b = np.array([4,2,3])
        assert np.allclose(a, b, atol=1.e-16)

    def test_change_of_basis(self):
        change_basis = ChangeOfBasis(np.array([[-2,0,0],[0,-2,0],[0,0,-2]]))
        assert (change_basis.apply_to(np.array([[2,3,4]])) == np.array([[-1.,-1.5,-2.]])).all()

    def test_axis_transform(self):
        origin_transform = AxisTransform(new_origin=np.array([1,1,1]))
        assert (origin_transform.apply_to(np.array([[1,1,1]])) == np.array([[0,0,0]])).all()
        orientation_transform = AxisTransform(point_on_x_axis=np.array([0,0,1]),point_on_xy_plane=np.array([0,1,1]))
        assert (orientation_transform.apply_to(np.array([[2,3,4]])) == np.array([[4,3,-2]])).all()
        axis_transform = AxisTransform(np.array([1,1,1]),np.array([1,1,2]),np.array([1,2,1]))
        assert (axis_transform.apply_to(np.array([[2,3,4]])) == np.array([3,2,-1])).all()

    def test_rigid_transform(self):
        A = np.array([[2,3,4]])
        B = np.array([[3,-4,9]])
        rigid_transform = RigidTransform(A, B)
        assert (rigid_transform.apply_to(np.array([[2,3,4]])) == B).all()

if __name__=='__main__':
    TestCoordinateTransform().test_apply_to()
    TestCoordinateTransform().test_translation()
    TestCoordinateTransform().test_rotation_around_z()
    TestCoordinateTransform().test_rotation_around_y()
    TestCoordinateTransform().test_rotation_around_x()
    TestCoordinateTransform().test_rotation()
    TestCoordinateTransform().test_change_of_basis()
    TestCoordinateTransform().test_axis_transform()
    TestCoordinateTransform().test_rigid_transform()