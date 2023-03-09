import numpy as np

import mbuild as mb
from mbuild.formats.json_formats import compound_from_json, compound_to_json
from mbuild.tests.base_test import BaseTest


class TestJSONFormats(BaseTest):
    def test_loop(self, ethane):
        for part in ethane:
            part.element = part.name
        compound_to_json(ethane, "ethane.json")
        ethane_copy = compound_from_json("ethane.json")
        for part_orig, part_copy in zip(ethane, ethane_copy):
            assert part_orig.element.symbol == part_copy.element.symbol
        assert ethane.n_particles == ethane_copy.n_particles
        assert ethane.n_bonds == ethane_copy.n_bonds
        assert len(ethane.children) == len(ethane_copy.children)

    def test_box(self):
        meth = mb.load("C", smiles=True)

        meth.save("meth_no_box.json")
        meth_no_box = mb.load("meth_no_box.json")
        assert meth_no_box.box is None

        meth.box = mb.Box(lengths=(3, 3, 3), angles=(45, 45, 45))
        meth.save("meth_with_box.json")
        meth_with_box = mb.load("meth_with_box.json")
        assert (
            meth.n_particles
            == meth_with_box.n_particles
            == meth_no_box.n_particles
        )
        assert meth.n_bonds == meth_with_box.n_bonds == meth_no_box.n_bonds
        assert meth.box.lengths == meth_with_box.box.lengths == (3, 3, 3)
        assert meth.box.angles == meth_with_box.box.angles

    def test_loop_with_ports(self):
        from mbuild.lib.moieties import CH3

        ethane_without_overlap = mb.Compound()
        methyl1 = CH3()
        methyl2 = CH3()
        ethane_without_overlap.add(methyl1, label="methyl1")
        ethane_without_overlap.add(methyl2, label="methyl2")
        compound_to_json(
            ethane_without_overlap,
            "ethane_without_overlap.json",
            include_ports=True,
        )
        ethane_copy = compound_from_json("ethane_without_overlap.json")
        assert ethane_copy.n_particles == ethane_without_overlap.n_particles
        assert ethane_copy.n_bonds == ethane_without_overlap.n_bonds
        assert len(ethane_copy.children) == len(ethane_without_overlap.children)
        assert len(ethane_copy.all_ports()) == len(
            ethane_without_overlap.all_ports()
        )
        assert ethane_copy.labels.keys() == ethane_without_overlap.labels.keys()
        assert (
            ethane_without_overlap["methyl2"].labels.keys()
            == ethane_copy["methyl2"].labels.keys()
        )

    def test_loop_for_propyl(self, hexane):
        compound_to_json(hexane, "hexane.json", include_ports=True)
        hexane_copy = compound_from_json("hexane.json")
        assert hexane.n_particles == hexane_copy.n_particles
        assert hexane.n_bonds == hexane_copy.n_bonds
        assert len(hexane.children) == len(hexane_copy.children)
        assert len(hexane.all_ports()) == len(hexane_copy.all_ports())
        assert hexane.labels.keys() == hexane_copy.labels.keys()

    def test_nested_compound(self):
        num_chidren = 10
        num_grand_children = 10
        num_ports = 2
        ancestor = mb.Compound(name="Ancestor")
        for i in range(num_chidren):
            this_child = mb.Compound(name="Child{}".format(i + 1))
            ancestor.add(this_child, label="Ancestor'sChild{}".format(i + 1))
            for j in range(num_ports):
                port1 = mb.Port(anchor=this_child)
                this_child.add(port1, label="port{}".format(j + 1))
            for k in range(num_grand_children):
                this_grand_child = mb.Compound(
                    name="GrandChild{}".format(k + 1)
                )
                this_child.add(
                    this_grand_child,
                    label="Child{0}GrandChild{1}".format(i + 1, k + 1),
                )
        compound_to_json(ancestor, "large_compound.json", include_ports=True)
        ancestor_copy = compound_from_json("large_compound.json")
        assert ancestor.n_particles == ancestor_copy.n_particles
        assert ancestor.n_bonds == ancestor_copy.n_bonds
        assert len(ancestor.children) == len(ancestor_copy.children)
        assert len(ancestor.all_ports()) == len(ancestor_copy.all_ports())
        assert ancestor.labels.keys() == ancestor_copy.labels.keys()

    def test_label_consistency(self):
        from mbuild.lib.moieties import CH2, CH3

        parent = mb.Compound(name="Hierarchy1")
        for i in range(10):
            parent.add(CH2())
            parent.add(CH3())
        compound_to_json(parent, "parent.json", include_ports=True)
        parent_copy = compound_from_json("parent.json")
        assert len(parent_copy["CH2"]) == len(parent["CH2"])
        assert parent_copy.labels.keys() == parent.labels.keys()
        for child, child_copy in zip(
            parent.successors(), parent_copy.successors()
        ):
            assert child.labels.keys() == child_copy.labels.keys()
        assert parent_copy.available_ports() == parent.available_ports()

    def test_float_64_position(self):
        ethane = mb.lib.molecules.Ethane()
        ethane.xyz = np.asarray(ethane.xyz, dtype=np.float64)
        compound_to_json(ethane, "ethane.json", include_ports=True)
        ethane_copy = compound_from_json("ethane.json")
        assert np.allclose(ethane.xyz, ethane_copy.xyz, atol=10**-6)
