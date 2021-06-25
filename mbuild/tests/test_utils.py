import difflib
from warnings import warn

import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.conversion import OPLS_to_RB, RB_to_OPLS
from mbuild.utils.exceptions import RemovedFuncError
from mbuild.utils.geometry import wrap_coords
from mbuild.utils.io import get_fn, import_, run_from_ipython
from mbuild.utils.jsutils import overwrite_nglview_default
from mbuild.utils.orderedset import OrderedSet
from mbuild.utils.validation import assert_port_exists


class TestUtils(BaseTest):
    def test_orderedset(self):
        oset_empty = OrderedSet()

        assert isinstance(oset_empty, OrderedSet)
        assert len(oset_empty) == 0

        oset = OrderedSet("heck", 0, 1)

        assert isinstance(oset, OrderedSet)
        assert len(oset) == 3
        assert "heck" in oset
        assert oset[1] == 0

        oset.add("this")

        assert "this" in oset

        oset.discard("heck")

        assert "heck" not in oset

        assert [i for i in oset] == [0, 1, "this"]

    def test_orderedset_setmethods(self):
        oset = OrderedSet(1, 2, 3, 4)
        oset2 = OrderedSet(2, 4, 6, 8)

        union = oset.union(oset2)
        union2 = oset2.union(oset)

        assert union == union2
        assert union == set([1, 2, 3, 4, 6, 8])

        inter = oset.intersection(oset2)
        inter2 = oset2.intersection(oset)

        assert inter == inter2
        assert inter == set([2, 4])

        diff = oset.difference(oset2)
        diff2 = oset2.difference(oset)

        assert diff == set([1, 3])
        assert diff2 == set([6, 8])

    def test_assert_port_exists(self, ch2):
        assert_port_exists("up", ch2)
        with pytest.raises(ValueError):
            assert_port_exists("dog", ch2)

    @pytest.mark.xfail(strict=False)
    def test_structure_reproducibility(self):
        from mbuild.lib.recipes import Alkane

        filename = "decane-tmp.xyz"
        decane = Alkane(10)
        decane.save(filename)
        with open(get_fn("decane.xyz")) as file1:
            with open(filename) as file2:
                diff = difflib.ndiff(file1.readlines(), file2.readlines())
        changes = [l for l in diff if l.startswith("+ ") or l.startswith("- ")]
        assert not changes

    def test_fn(self):
        get_fn("benzene.mol2")

        with pytest.raises((IOError, OSError)):
            get_fn("garbage_file_name.foo")

    def test_import(self):
        assert np == import_("numpy")

        with pytest.raises(ImportError):
            import_("garbagepackagename")

    def test_js_utils(self):
        nglview = import_("nglview")
        with pytest.raises(TypeError):
            overwrite_nglview_default(object())
        test_widget = nglview.NGLWidget()
        overwrite_nglview_default(test_widget)
        assert hasattr(test_widget, "stage")
        assert isinstance(test_widget._ngl_msg_archive, list)
        assert len(test_widget._ngl_msg_archive) == 2
        assert isinstance(test_widget._ngl_msg_archive[0], dict)
        message_dict = test_widget._ngl_msg_archive[0]
        assert message_dict["target"] == "Widget"
        assert message_dict["type"] == "call_method"
        assert message_dict["methodName"] == "executeCode"
        assert message_dict["args"] == [
            """
                    this.stage.mouseControls.add('hoverPick', (stage, pickingProxy) => {
                        let tooltip = this.stage.tooltip;
                        if(pickingProxy && pickingProxy.atom && !pickingProxy.bond){
                            let atom = pickingProxy.atom;
                            tooltip.innerText = "ATOM: " + atom.qualifiedName() + ", Index: " + atom.index;
                        }
                    });
                 """
        ]
        message_dict = test_widget._ngl_msg_archive[1]
        assert message_dict["target"] == "Widget"
        assert message_dict["type"] == "call_method"
        assert message_dict["methodName"] == "executeCode"
        assert message_dict["args"] == [
            """
                    this.stage.signals.clicked.removeAll();
                    this.stage.signals.clicked.add((pickingProxy) => {
                            if(pickingProxy){
                               let pickingText = null;
                               this.model.set('picked', {});
                               this.touch();
                               let currentPick = {};
                               if(pickingProxy.atom){
                                    currentPick.atom1 = pickingProxy.atom.toObject();
                                    currentPick.atom1.name = pickingProxy.atom.qualifiedName();
                                    pickingText = "Atom: " + currentPick.atom1.name + ", Index: "
                                                  + pickingProxy.atom.index;
                               }
                               else if(pickingProxy.bond){
                                    currentPick.bond = pickingProxy.bond.toObject();
                                    currentPick.atom1 = pickingProxy.bond.atom1.toObject();
                                    currentPick.atom1.name = pickingProxy.bond.atom1.qualifiedName();
                                    currentPick.atom2 = pickingProxy.bond.atom2.toObject();
                                    currentPick.atom2.name = pickingProxy.bond.atom2.qualifiedName();
                                    pickingText = "Bond: " + currentPick.atom1.name +
                                                    `(${pickingProxy.bond.atom1.index})` +
                                                    " - " + currentPick.atom2.name    +
                                                    `(${pickingProxy.bond.atom2.index})`;
                               }

                               if(pickingProxy.instance){
                                    currentPick.instance = pickingProxy.instance;
                               }
                               var nComponents = this.stage.compList.length;
                               for(let i = 0; i < nComponents; i++){
                                    let comp = this.stage.compList[i];
                                    if(comp.uuid == pickingProxy.component.uuid){
                                        currentPick.component = i;
                                    }
                               }
                               this.model.set('picked', currentPick);
                               this.touch();
                               this.$pickingInfo.text(pickingText);
                            }
                    });
                """
        ]

    def test_coord_wrap(self):
        xyz = np.array([[3, 3, 1], [1, 1, 0]])
        box = [2, 2, 2]
        new_xyz = wrap_coords(xyz, box)
        assert (new_xyz[0, :] == np.array([1, 1, 1])).all()
        assert (new_xyz[1, :] == xyz[1, :]).all()

    def test_coord_wrap_box(self):
        xyz = np.array([[3, 3, 1], [1, 1, 0]])

        box = mb.Box(lengths=[2.0, 2.0, 2.0], angles=[90, 90, 90])

        new_xyz = wrap_coords(xyz, box, mins=np.min(xyz, axis=0))

        assert (new_xyz[0, :] == np.array([1, 1, 1])).all()
        assert (new_xyz[1, :] == xyz[1, :]).all()

    def test_has_ipython(self):
        __IPYTHON__ = None
        assert run_from_ipython() is False

    def test_removed_func_error(self):
        RemovedFuncError("old_function", "new_function", "0.0.0", "0.1.0")


class TestUtilsConversion(BaseTest):
    def test_RB_to_OPLS_c5_not_0(self):
        with pytest.raises(
            ValueError,
            match=r"c5 must equal zero, so this conversion is not possible.",
        ):
            c0 = 0.1
            c1 = 0.1
            c2 = -0.2
            c3 = -0.1
            c4 = -0.2
            c5 = 0.3
            RB_to_OPLS(c0, c1, c2, c3, c4, c5)

    def test_RB_to_OPLS_f0_not_0_within_tolerance_warn(self):
        # should throw a warning that f0 is not zero
        text_for_error_tolerance = (
            "f0 \= 2 \* \( c0 \+ c1 \+ c2 \+ c3 \+ c4 \+ c5 \) is not zero. "
            "The f0/2 term is the constant for the OPLS dihedral. "
            "Since the f0 term is not zero, the dihedral is not an "
            "exact conversion; since this constant does not contribute "
            "to the force equation, this should provide matching results "
            "for MD, but the energy for each dihedral will be shifted "
            "by the f0/2 value."
        )

        with pytest.warns(UserWarning, match=f"{text_for_error_tolerance}"):
            c0 = 0.4
            c1 = 0.4
            c2 = -0.1
            c3 = 0.4
            c4 = -0.2
            c5 = 0
            RB_to_OPLS(c0, c1, c2, c3, c4, c5, error_if_outside_tolerance=False)

    def test_RB_to_OPLS_f0_not_0_within_tolerance_error(self):
        text_for_error_tolerance = (
            "f0 \= 2 \* \( c0 \+ c1 \+ c2 \+ c3 \+ c4 \+ c5 \) is not zero. "
            "The f0/2 term is the constant for the OPLS dihedral. "
            "Since the f0 term is not zero, the dihedral is not an "
            "exact conversion; since this constant does not contribute "
            "to the force equation, this should provide matching results "
            "for MD, but the energy for each dihedral will be shifted "
            "by the f0/2 value."
        )

        with pytest.raises(
            TypeError,
            match=f"{text_for_error_tolerance}",
        ):
            c0 = 0.4
            c1 = 0.4
            c2 = -0.1
            c3 = 0.4
            c4 = -0.2
            c5 = 0
            RB_to_OPLS(c0, c1, c2, c3, c4, c5)

    def test_RB_to_OPLS_f0_not_0_within_tolerance_error(self):
        with pytest.raises(
            TypeError,
            match=f"The error_tolerance variable must be a float, is type {type('s')}.",
        ):
            c0 = 0.1
            c1 = 0.1
            c2 = -0.2
            c3 = -0.1
            c4 = -0.2
            c5 = 0.0
            RB_to_OPLS(
                c0,
                c1,
                c2,
                c3,
                c4,
                c5,
                error_tolerance="s",
                error_if_outside_tolerance=False,
            )

    def test_RB_to_OPLS_text_for_error_tolerance_not_bool(self):
        with pytest.raises(
            TypeError,
            match=f"The text_for_error_tolerance variable must be a bool, is type {type('s')}.",
        ):
            c0 = 0.1
            c1 = 0.1
            c2 = -0.2
            c3 = -0.1
            c4 = -0.2
            c5 = 0.0
            RB_to_OPLS(c0, c1, c2, c3, c4, c5, error_if_outside_tolerance="s")

    @pytest.mark.parametrize(
        "c0, c1, c2, c3, c4, c5",
        [
            (-8.2723, 0.2263, 10.22, -3.208, 1.034, 0.0),
            (-0.363, 2.726, 2.849, 7.373, -12.585, 0),
            (0, -2.3927, 0, 9.17, -6.7773, 0),
            (10, 2.37, 0, 0, -12.37, 0),
            (2.10, -8.22, 1.22, 7.20, -2.3, 0),
            (-12, 1.24, -22.816, 11.28654, 22.28946, 0),
            (-0.3936, -8.283, 3.183, -3.28, 8.7736, 0),
            (-8.7, 0, 0, 33.12, -24.42, 0),
            (-12, 1.393, -1.234, 0.2323, 11.6087, 0),
            (10.22, 3.293, -34.32, 2.596, 18.211, 0),
            (3.28629, 7.44211, 1.85995, -14.67569, 2.08734, 0),
            (5.77183, -2.67148, 0.95814, -4.05848, -0.00001, 0),
            (5.77183, -2.67148, 0.95814, -4.05848, 20, 0),
            (10.1, 10.1, 10.1, 10.1, 10.1, 0),
            (-5.5, -5.5, -5.5, -5.5, -5.5, 0),
        ],
    )
    def test_RB_to_OPLS_and_back_random_values(self, c0, c1, c2, c3, c4, c5):
        # However, this may not be true for real dihedrals.
        test_error_tolerance = 1e-10

        opls_coeffs = RB_to_OPLS(
            c0, c1, c2, c3, c4, c5, error_if_outside_tolerance=False
        )
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[0],
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert np.all(
            np.isclose(
                c0, reversed_RB_coeffs[0], atol=test_error_tolerance, rtol=0
            )
        )
        assert np.all(
            np.isclose(
                c1, reversed_RB_coeffs[1], atol=test_error_tolerance, rtol=0
            )
        )
        assert np.all(
            np.isclose(
                c2, reversed_RB_coeffs[2], atol=test_error_tolerance, rtol=0
            )
        )
        assert np.all(
            np.isclose(
                c3, reversed_RB_coeffs[3], atol=test_error_tolerance, rtol=0
            )
        )
        assert np.all(
            np.isclose(
                c4, reversed_RB_coeffs[4], atol=test_error_tolerance, rtol=0
            )
        )
        assert np.all(
            np.isclose(
                c5, reversed_RB_coeffs[5], atol=test_error_tolerance, rtol=0
            )
        )

    def test_OPLS_to_RB_error_tolerance_not_float(self):
        with pytest.raises(
            TypeError,
            match=f"The error_tolerance variable must be a float, is type {type('s')}.",
        ):
            f0 = 0.1
            f1 = 0.1
            f2 = -0.2
            f3 = -0.1
            f4 = 0.2
            OPLS_to_RB(f0, f1, f2, f3, f4, error_tolerance="s")

    def test_OPLS_to_RB_f0_is_zero(self):
        # should throw a warning that f0 is zero
        text_for_error_tolerance = (
            "The f0/2 term is the constant for the OPLS dihedral equation, "
            "which is added to a constant for the RB torsions equation via the c0 coefficient. "
            "The f0 term is zero in the OPLS dihedral form or is force set to zero in this equation, "
            "so the dihedral is may not an exact conversion; "
            "since this constant does not contribute to the force equation, "
            "this should provide matching results for MD, but the energy for each"
            "dihedral will be shifted by the real f0/2 value."
        )
        with pytest.warns(UserWarning, match=text_for_error_tolerance):
            f0 = 0
            f1 = 0.1
            f2 = -0.2
            f3 = -0.1
            f4 = 0.2
            OPLS_to_RB(f0, f1, f2, f3, f4)
