import math

import numpy as np
import pytest

from mbuild.path.build import Path, lamellar
from mbuild.path.namers import CyclicNamer
from mbuild.tests.base_test import BaseTest

# Fractional centering translations expected of each crystal lattice, written
# out here so the tests do not depend on the table in the module under test.
TRANSLATIONS = {
    "primitive": ((0.0, 0.0, 0.0),),
    "base_centered": ((0.0, 0.0, 0.0), (0.5, 0.5, 0.0)),
    "body_centered": ((0.0, 0.0, 0.0), (0.5, 0.5, 0.5)),
    "face_centered": (
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    ),
}
FACE_TRANSLATIONS = {
    "ab": (0.5, 0.5, 0.0),
    "ac": (0.5, 0.0, 0.5),
    "bc": (0.0, 0.5, 0.5),
}
CELL = (0.5, 0.5, 0.2)
REPEATS = (2, 2, 4)


def lattice_sites(cell, repeat_units, translations, initial_point=(0, 0, 0)):
    """Every site of the decorated cell, as a set of rounded coordinate tuples."""
    a, b, c = cell
    num_layers, num_stacks, num_sites = repeat_units
    sites = set()
    for offset_a, offset_b, offset_c in translations:
        for layer in range(num_layers):
            for stack in range(num_stacks):
                for site in range(num_sites):
                    sites.add(
                        (
                            round(initial_point[0] + (layer + offset_a) * a, 9),
                            round(initial_point[1] + (stack + offset_b) * b, 9),
                            round(initial_point[2] + (site + offset_c) * c, 9),
                        )
                    )
    return sites


def column_runs(coordinates, sites):
    """Split lattice beads into consecutive runs sharing an (x, y) column."""
    runs = []
    for point in coordinates:
        rounded = tuple(np.round(point, 9))
        if rounded not in sites:
            continue
        if not runs or runs[-1][0] != rounded[:2]:
            runs.append((rounded[:2], [rounded[2]]))
        else:
            runs[-1][1].append(rounded[2])
    return runs


class TestLamellar(BaseTest):
    def test_lamellar_direction_start(self):
        path_left_to_right = lamellar(
            cell=(1, 3, 1),
            repeat_units=(3, 3, 3),
        )

        path_right_to_left = lamellar(
            cell=(1, 3, 1),
            repeat_units=(3, 3, 3),
            left_to_right=False,
        )

        assert np.array_equal(
            path_left_to_right.coordinates[0], path_right_to_left.coordinates[0]
        )
        assert path_right_to_left.coordinates[1][2] < 0
        assert path_left_to_right.coordinates[1][2] > 0

    def test_lamellar(self):
        path = lamellar(cell=CELL, repeat_units=(3, 3, 5))
        assert path.bond_graph.number_of_edges() == len(path.coordinates) - 1
        compound = path.to_compound()
        Lx, Ly, Lz = compound.get_boundingbox().lengths
        # a and b span the layer and stack directions exactly; the folds bulge
        # past the ends of the columns, so z is longer than (num_sites - 1) * c.
        assert np.allclose(Lx, 2 * CELL[0])
        assert np.allclose(Ly, 2 * CELL[1])
        assert 4 * CELL[2] < Lz <= 4 * CELL[2] + CELL[0]

    def test_lamellar_bonding(self):
        path = lamellar(cell=CELL, repeat_units=REPEATS)
        assert path.bond_graph.number_of_nodes() == len(path.coordinates)
        assert path.bond_graph.number_of_edges() == len(path.coordinates) - 1
        compound = path.to_compound()
        assert compound.n_bonds == compound.n_particles - 1

    def test_lamellar_single_cell(self):
        path = lamellar(cell=CELL, repeat_units=(1, 1, 1))
        assert len(path.coordinates) == 1
        assert path.bond_graph.number_of_edges() == 0
        assert np.allclose(path.coordinates[0], (0, 0, 0))

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_visits_every_site(self, lattice):
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        sites = lattice_sites(CELL, REPEATS, TRANSLATIONS[lattice])
        visited = {tuple(np.round(point, 9)) for point in path.coordinates}
        assert sites <= visited
        assert len(sites) == len(TRANSLATIONS[lattice]) * math.prod(REPEATS)

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_visits_every_site_once(self, lattice):
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        sites = lattice_sites(CELL, REPEATS, TRANSLATIONS[lattice])
        rounded = [tuple(np.round(point, 9)) for point in path.coordinates]
        assert len([point for point in rounded if point in sites]) == len(sites)
        assert len(set(rounded)) == len(rounded)

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_sites_on_half_cell_fractions(self, lattice):
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        sites = lattice_sites(CELL, REPEATS, TRANSLATIONS[lattice])
        for point in path.coordinates:
            if tuple(np.round(point, 9)) in sites:
                doubled = 2.0 * point / np.asarray(CELL)
                assert np.allclose(doubled, np.round(doubled))

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_columns_entered_once(self, lattice):
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        sites = lattice_sites(CELL, REPEATS, TRANSLATIONS[lattice])
        keys = [key for key, _ in column_runs(path.coordinates, sites)]
        assert len(keys) == len(set(keys))
        assert len(keys) == len(TRANSLATIONS[lattice]) * REPEATS[0] * REPEATS[1]

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_column_runs_monotonic(self, lattice):
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        sites = lattice_sites(CELL, REPEATS, TRANSLATIONS[lattice])
        for _, heights in column_runs(path.coordinates, sites):
            assert len(heights) == REPEATS[2]
            steps = np.diff(heights)
            assert np.allclose(np.abs(steps), CELL[2])
            assert np.all(steps > 0) or np.all(steps < 0)

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_bond_lengths(self, lattice):
        """Verify that bonds lengths are reasonable."""
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        bonds = np.linalg.norm(np.diff(path.coordinates, axis=0), axis=1)
        assert np.all(bonds <= CELL[2] + 1e-8)

    @pytest.mark.parametrize(
        "lattice, num_beads", [("primitive", 25), ("base-centered", 50)]
    )
    def test_lamellar_bead_count(self, lattice, num_beads):
        # 4 columns of 4 sites joined by 3 arcs of 3 beads, doubled and
        # re-threaded for the C-centered lattice.
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        assert len(path.coordinates) == num_beads

    @pytest.mark.parametrize(
        "lattice", ["base_centered", "body_centered", "face_centered"]
    )
    def test_lamellar_centering_adds_beads(self, lattice):
        primitive = lamellar(cell=CELL, repeat_units=REPEATS)
        centered = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)
        assert len(centered.coordinates) > len(primitive.coordinates)
        # Centering only adds sites; the primitive ones stay where they were.
        visited = {tuple(np.round(point, 9)) for point in centered.coordinates}
        assert lattice_sites(CELL, REPEATS, TRANSLATIONS["primitive"]) <= visited

    @pytest.mark.parametrize("face", list(FACE_TRANSLATIONS))
    def test_lamellar_base_centered_face(self, face):
        path = lamellar(
            cell=CELL,
            repeat_units=REPEATS,
            lattice_centerings="base_centered",
            base_centered_face=face,
        )
        visited = {tuple(np.round(point, 9)) for point in path.coordinates}
        expected = lattice_sites(
            CELL, REPEATS, ((0.0, 0.0, 0.0), FACE_TRANSLATIONS[face])
        )
        assert expected <= visited
        for other, translation in FACE_TRANSLATIONS.items():
            if other != face:
                unique = lattice_sites(CELL, REPEATS, (translation,)) - expected
                assert not unique <= visited

    def test_lamellar_default_base_centered_face(self):
        default = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings="base_centered"
        )
        explicit = lamellar(
            cell=CELL,
            repeat_units=REPEATS,
            lattice_centerings="base_centered",
            base_centered_face="ab",
        )
        assert np.allclose(default.coordinates, explicit.coordinates)

    @pytest.mark.parametrize("lattice", ["primitive", "body_centered", "face_centered"])
    def test_lamellar_base_centered_face_ignored(self, lattice):
        reference = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice
        )
        for face in FACE_TRANSLATIONS:
            path = lamellar(
                cell=CELL,
                repeat_units=REPEATS,
                lattice_centerings=lattice,
                base_centered_face=face,
            )
            assert np.allclose(path.coordinates, reference.coordinates)

    @pytest.mark.parametrize("lattice", list(TRANSLATIONS))
    def test_lamellar_initial_point(self, lattice):
        shift = np.array([1.5, -2.0, 0.25])
        at_origin = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice
        ).coordinates
        shifted = lamellar(
            cell=CELL,
            repeat_units=REPEATS,
            lattice_centerings=lattice,
            initial_point=shift,
        ).coordinates
        assert np.allclose(shifted, at_origin + shift)

    def test_lamellar_direction(self):
        left_to_right = lamellar(cell=CELL, repeat_units=REPEATS)
        right_to_left = lamellar(cell=CELL, repeat_units=REPEATS, left_to_right=False)
        # The first column is traversed up +z, or down from the far end of it.
        assert np.allclose(left_to_right.coordinates[0], (0, 0, 0))
        assert left_to_right.coordinates[1][2] > left_to_right.coordinates[0][2]
        assert np.allclose(right_to_left.coordinates[0], (0, 0, 0))
        assert right_to_left.coordinates[1][2] < 0

    def test_lamellar_appends_to_path(self):
        path = lamellar(cell=CELL, repeat_units=REPEATS)
        count = len(path.coordinates)
        lamellar(
            path=path,
            cell=CELL,
            repeat_units=REPEATS,
            initial_point=(10.0, 0.0, 0.0),
        )
        assert len(path.coordinates) == 2 * count
        # The two segments are bonded separately, not to each other.
        assert path.bond_graph.number_of_edges() == 2 * (count - 1)
        assert np.allclose(
            path.coordinates[count:] - np.array([10.0, 0.0, 0.0]),
            path.coordinates[:count],
        )

    @pytest.mark.parametrize(
        "spelling",
        ["face_centered", "face-centered", "FACE CENTERED", " Face-Centered "],
    )
    def test_lamellar_lattice_centerings_spellings(self, spelling):
        reference = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings="face_centered"
        )
        path = lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=spelling)
        assert np.allclose(path.coordinates, reference.coordinates)

    @pytest.mark.parametrize("repeats", [[2, 2, 4], np.array([2, 2, 4])])
    def test_lamellar_repeat_units_types(self, repeats):
        reference = lamellar(cell=CELL, repeat_units=REPEATS)
        path = lamellar(cell=CELL, repeat_units=repeats)
        assert np.allclose(path.coordinates, reference.coordinates)

    def test_lamellar_is_deterministic(self):
        first = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings="face_centered"
        )
        second = lamellar(
            cell=CELL, repeat_units=REPEATS, lattice_centerings="face_centered"
        )
        assert np.allclose(first.coordinates, second.coordinates)

    def test_lamellar_cyclic_namer(self):
        path = Path()
        lamellar(
            path=path,
            cell=CELL,
            repeat_units=(1, 1, 4),
            bead_name=CyclicNamer(["_A", "_B"]),
        )
        assert list(path.beads) == ["_A", "_B", "_A", "_B"]

    def test_lamellar_string_bead_name(self):
        path = lamellar(cell=CELL, repeat_units=REPEATS, bead_name="_X")
        assert set(path.beads) == {"_X"}
        node_names = [d["name"] for _, d in path.bond_graph.nodes(data=True)]
        assert node_names == ["_X"] * len(path.coordinates)

    def test_lamellar_no_cell(self):
        path = lamellar(repeat_units=REPEATS)
        assert np.linalg.norm(path.coordinates[0] - path.coordinates[1]) == 1.0

    @pytest.mark.parametrize("cell", [(0.5, 0.5), (0.5, 0.5, 0.2, 0.1)])
    def test_lamellar_bad_cell_shape(self, cell):
        with pytest.raises(ValueError):
            lamellar(cell=cell, repeat_units=REPEATS)

    @pytest.mark.parametrize("cell", [(0.0, 0.5, 0.2), (0.5, -0.5, 0.2)])
    def test_lamellar_nonpositive_cell(self, cell):
        with pytest.raises(ValueError):
            lamellar(cell=cell, repeat_units=REPEATS)

    @pytest.mark.parametrize("repeats", [(2, 2), (2, 2, 4, 4), (2, 2, 0), (2, -1, 4)])
    def test_lamellar_bad_repeat_units(self, repeats):
        with pytest.raises(ValueError):
            lamellar(cell=CELL, repeat_units=repeats)

    @pytest.mark.parametrize("repeats", [(2, 2, 4.0), (2.5, 2, 4)])
    def test_lamellar_noninteger_repeat_units(self, repeats):
        with pytest.raises(TypeError):
            lamellar(cell=CELL, repeat_units=repeats)

    @pytest.mark.parametrize("lattice", ["hexagonal", "side_centered", "primitve", ""])
    def test_lamellar_bad_lattice_centerings(self, lattice):
        with pytest.raises(ValueError):
            lamellar(cell=CELL, repeat_units=REPEATS, lattice_centerings=lattice)

    @pytest.mark.parametrize("face", ["ba", "abc", "a", ""])
    def test_lamellar_bad_base_centered_face(self, face):
        with pytest.raises(ValueError):
            lamellar(
                cell=CELL,
                repeat_units=REPEATS,
                lattice_centerings="base_centered",
                base_centered_face=face,
            )
