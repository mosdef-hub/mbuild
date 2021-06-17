import difflib

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
            match=r"ERROR: c5 must equal zero, so this conversion is not possible.",
        ):
            c0 = 0.1
            c1 = 0.1
            c2 = -0.2
            c3 = -0.1
            c4 = -0.2
            c5 = 0.3
            RB_to_OPLS(c0, c1, c2, c3, c4, c5)

    def test_RB_to_OPLS_f0_not_0(self):
        with pytest.raises(
            ValueError,
            match=r"ERROR: f0 \= 2 \* \(c0 \+ c1 \+ c2 \+ c3 \+ c4 \+ c5\) must equal zero, "
            "so this conversion is not possible.",
        ):
            c0 = 0.4
            c1 = 0.4
            c2 = -0.1
            c3 = 0.4
            c4 = -0.2
            c5 = 0
            RB_to_OPLS(c0, c1, c2, c3, c4, c5)

    def test_RB_to_OPLS_and_back_random_values(self):
        # Note the sum of c0 to c5 must be zero for these test.
        # However, this may not be true for real dihedrals.

        # Test Number 1
        c0 = -8.2723
        c1 = 0.2263
        c2 = 10.22
        c3 = -3.208
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 2
        c0 = -0.363
        c1 = 2.726
        c2 = 2.849
        c3 = 7.373
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 3
        c0 = 0
        c1 = -2.3927
        c2 = 0
        c3 = 9.17
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 4
        c0 = 10
        c1 = 2.37
        c2 = 0
        c3 = 0
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 5
        c0 = 2.10
        c1 = -8.22
        c2 = 1.22
        c3 = 7.20
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 6
        c0 = -12
        c1 = 1.24
        c2 = -22.816
        c3 = 11.28654
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 7
        c0 = -0.3936
        c1 = -8.283
        c2 = 3.183
        c3 = -3.28
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 8
        c0 = -8.7
        c1 = 0
        c2 = 0
        c3 = 33.12
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 9
        c0 = -12
        c1 = 1.393
        c2 = -1.234
        c3 = 0.2323
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 10
        c0 = 10.22
        c1 = 3.293
        c2 = -34.32
        c3 = 2.596
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 11
        c0 = 3.28629
        c1 = 7.44211
        c2 = 1.85995
        c3 = -14.67569
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))

        # Test Number 12
        c0 = 5.77183
        c1 = -2.67148
        c2 = 0.95814
        c3 = -4.05848
        c4 = -(c0 + c1 + c2 + c3)
        c5 = 0

        opls_coeffs = RB_to_OPLS(c0, c1, c2, c3, c4, c5)
        reversed_RB_coeffs = OPLS_to_RB(
            opls_coeffs[1],
            opls_coeffs[2],
            opls_coeffs[3],
            opls_coeffs[4],
        )

        assert bool(np.isclose(c0, reversed_RB_coeffs[0], atol=1e-12, rtol=0))
        assert bool(np.isclose(c1, reversed_RB_coeffs[1], atol=1e-12, rtol=0))
        assert bool(np.isclose(c2, reversed_RB_coeffs[2], atol=1e-12, rtol=0))
        assert bool(np.isclose(c3, reversed_RB_coeffs[3], atol=1e-12, rtol=0))
        assert bool(np.isclose(c4, reversed_RB_coeffs[4], atol=1e-12, rtol=0))
        assert bool(np.isclose(c5, reversed_RB_coeffs[5], atol=1e-12, rtol=0))
