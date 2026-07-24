import networkx as nx
import numpy as np
import pytest

from mbuild.lib.molecules import (
    FattyAcid,
    WaterOPC,
    WaterOPC3,
    WaterSPC,
    WaterTIP3P,
    WaterTIP4P,
    WaterTIP4P2005,
    WaterTIP4PIce,
)
from mbuild.tests.base_test import BaseTest


class TestWater(BaseTest):
    @pytest.mark.parametrize(
        "model, bond_length, angle",
        [
            (WaterTIP3P, 0.09572, 104.52),
            (WaterSPC, 0.1, 109.47),
            (WaterOPC3, 0.09789, 109.47),
        ],
    )
    def test_water_3site(self, model, bond_length, angle):
        water = model()
        o1 = [p for p in water.particles_by_name("OW")]
        assert len(o1) == 1
        o1 = o1[0]
        h1 = [p for p in water.particles_by_name("HW1")]
        assert len(h1) == 1
        h1 = h1[0]
        h2 = [p for p in water.particles_by_name("HW2")]
        assert len(h2) == 1
        h2 = h2[0]
        v1 = h1.xyz[0] - o1.xyz[0]
        v2 = h2.xyz[0] - o1.xyz[0]
        assert np.allclose(np.linalg.norm(v1), bond_length)
        assert np.allclose(np.linalg.norm(v2), bond_length)
        assert np.allclose(
            np.degrees(
                np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            ),
            angle,
        )

    @pytest.mark.parametrize(
        "model, bond_length, vsite_length, angle",
        [
            (WaterTIP4P, 0.09572, 0.015, 104.52),
            (WaterTIP4PIce, 0.09572, 0.01577, 104.52),
            (WaterTIP4P2005, 0.09572, 0.01546, 104.52),
            (WaterOPC, 0.08724, 0.01594, 103.6),
        ],
    )
    def test_water_4site(self, model, bond_length, vsite_length, angle):
        water = model()
        o1 = [p for p in water.particles_by_name("OW")]
        assert len(o1) == 1
        o1 = o1[0]
        h1 = [p for p in water.particles_by_name("HW1")]
        assert len(h1) == 1
        h1 = h1[0]
        h2 = [p for p in water.particles_by_name("HW2")]
        assert len(h2) == 1
        h2 = h2[0]
        m1 = [p for p in water.particles_by_name("MW")]
        assert len(m1) == 1
        m1 = m1[0]
        v1 = h1.xyz[0] - o1.xyz[0]
        v2 = h2.xyz[0] - o1.xyz[0]
        v3 = m1.xyz[0] - o1.xyz[0]
        assert np.allclose(np.linalg.norm(v1), bond_length)
        assert np.allclose(np.linalg.norm(v2), bond_length)
        assert np.allclose(np.linalg.norm(v3), vsite_length)
        assert np.allclose(
            np.degrees(
                np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            ),
            angle,
        )


class TestFattyAcid(BaseTest):
    @pytest.mark.parametrize("chain_length", [4, 12, 18])
    def test_chain_length(self, chain_length):
        fatty_acid = FattyAcid(chain_length=chain_length)

        assert len(list(fatty_acid.particles_by_element("C"))) == chain_length

    @pytest.mark.parametrize(
        "double_bonds, expected_count",
        [
            (None, 0),
            ([(9, "cis")], 1),
            ([(9, "trans"), (12, "cis")], 2),
        ],
    )
    def test_number_of_double_bonds(self, double_bonds, expected_count):
        fatty_acid = FattyAcid(double_bonds=double_bonds)

        tail_double_bonds = [
            bond
            for bond in fatty_acid.bonds(return_bond_order=True)
            if bond[2]["bond_order"] == 2.0
            and bond[0].element.symbol == "C"
            and bond[1].element.symbol == "C"
        ]
        assert len(tail_double_bonds) == expected_count

    @pytest.mark.parametrize(
        "chain_length, double_bonds",
        [(12, [(4, "trans")]), (18, [(9, "cis")]), (18, [(9, "cis"), (12, "cis")])],
    )
    def test_double_bond_positions(self, chain_length, double_bonds):
        fatty_acid = FattyAcid(
            chain_length=chain_length,
            double_bonds=double_bonds,
        )
        bonds = list(fatty_acid.bonds(return_bond_order=True))

        carbonyl_carbon = next(
            atom1 if atom1.element.symbol == "C" else atom2
            for atom1, atom2, bond_data in bonds
            if bond_data["bond_order"] == 2.0
            and {atom1.element.symbol, atom2.element.symbol} == {"C", "O"}
        )

        carbon_double_bonds = [
            bond
            for bond in bonds
            if bond[2]["bond_order"] == 2.0
            and bond[0].element.symbol == "C"
            and bond[1].element.symbol == "C"
        ]

        actual_positions_in_chain = sorted(
            min(
                nx.shortest_path_length(
                    fatty_acid.bond_graph,
                    carbonyl_carbon,
                    atom,
                )
                for atom in bond[:2]
            )
            + 1  # since we count distances to reference, add 1 to get
            # back into carbon numbering
            for bond in carbon_double_bonds
        )

        expected_positions = sorted(position for position, _ in double_bonds)

        assert actual_positions_in_chain == expected_positions
