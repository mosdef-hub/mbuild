import networkx as nx
import numpy as np
import pytest

from mbuild.path.build import (
    Path,
    hard_sphere_random_walk,
    straight_line,
)
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_rdkit


def compound_graph(compound):
    """Build a particle-index graph from a compound's bonds."""
    particles = list(compound.particles())
    indices = {p: i for i, p in enumerate(particles)}
    graph = nx.Graph()
    graph.add_nodes_from(range(len(particles)))
    for p1, p2 in compound.bonds():
        graph.add_edge(indices[p1], indices[p2])
    return graph


@pytest.fixture
def branched_path():
    """A random walk backbone of A beads with a B bead branch on node 4."""
    path = hard_sphere_random_walk(
        termination=10, radius=0.2, bond_length=0.25, bead_name="A", seed=7
    )
    path = hard_sphere_random_walk(
        path=path,
        termination=15,
        bead_name="B",
        radius=0.2,
        bond_length=0.25,
        initial_point=4,
        connectivity="link-linear",
        seed=7,
    )
    return path


class TestCoarseGraining(BaseTest):
    def test_to_cgsmiles_graph(self):
        path = straight_line(spacing=0.25, N=4, bead_name="PEO")
        meta_graph = path.to_cgsmiles_graph()
        assert len(meta_graph) == 4
        assert meta_graph.number_of_edges() == 3
        for node in meta_graph.nodes:
            assert meta_graph.nodes[node]["fragname"] == "PEO"
            assert np.allclose(
                meta_graph.nodes[node]["position"], path.coordinates[node]
            )
        for edge in meta_graph.edges:
            assert meta_graph.edges[edge]["order"] == 1

    def test_to_cgsmiles_graph_fragname_map(self):
        path = straight_line(spacing=0.25, N=4)  # default bead name "_A"
        meta_graph = path.to_cgsmiles_graph(fragname_map={"_A": "PE"})
        fragnames = nx.get_node_attributes(meta_graph, "fragname")
        assert set(fragnames.values()) == {"PE"}

    def test_to_cgsmiles_string_linear(self):
        path = straight_line(spacing=0.25, N=3, bead_name="PEO")
        assert path.to_cgsmiles() == "{[#PEO][#PEO][#PEO]}"

    def test_to_cgsmiles_string_branched(self, branched_path):
        cgsmiles_str = branched_path.to_cgsmiles()
        # branches are written in parentheses
        assert "(" in cgsmiles_str and ")" in cgsmiles_str
        assert "[#A]" in cgsmiles_str and "[#B]" in cgsmiles_str
        # round trip through the cgsmiles reader gives an isomorphic graph
        from cgsmiles import read_cgsmiles

        read_back = read_cgsmiles(cgsmiles_str)
        assert nx.is_isomorphic(
            read_back,
            branched_path.bond_graph,
        )

    def test_resolve(self):
        from mbuild.coarse_graining import resolve

        path = straight_line(spacing=0.25, N=3, bead_name="PEO")
        meta_graph, molecule, node_to_beads = resolve(path, "{#PEO=[>]COC[<]}")
        elements = nx.get_node_attributes(molecule, "element")
        assert sum(1 for e in elements.values() if e == "C") == 6
        assert sum(1 for e in elements.values() if e == "O") == 3
        # interior monomer has 4 H, the two ends have 5 H each
        assert sum(1 for e in elements.values() if e == "H") == 14
        assert nx.is_connected(molecule)
        # every atom maps back to exactly one bead
        for node in molecule.nodes:
            assert len(node_to_beads[node]) == 1
            assert node_to_beads[node][0] in path.bond_graph.nodes

    def test_backmap_linear(self):
        path = straight_line(spacing=0.35, N=5, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        assert compound.n_particles == 37
        assert compound.n_bonds == 36
        assert nx.is_connected(compound_graph(compound))
        # hierarchy mirrors the path: one child compound per bead
        assert [c.name for c in compound.children] == ["PEO"] * 5

    def test_backmap_branched(self, branched_path):
        n_beads = len(branched_path.coordinates)
        degrees = dict(branched_path.bond_graph.degree())
        assert max(degrees.values()) == 3  # the case build_from_path can't handle
        compound = branched_path.backmap("{#A=[$]C([$])C[$],#B=[$]CC[$]}")
        assert nx.is_connected(compound_graph(compound))
        assert len(compound.children) == n_beads
        assert [c.name for c in compound.children] == list(branched_path.beads)

    def test_backmap_positions_centered(self):
        path = straight_line(spacing=0.35, N=4, bead_name="PE")
        compound = path.backmap("{#PE=[>]CC[<]}")
        for i, child in enumerate(compound.children):
            centroid = np.mean([p.pos for p in child.particles()], axis=0)
            assert np.linalg.norm(centroid - path.coordinates[i]) < 0.05

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_rdkit_positions(self, branched_path):
        compound = branched_path.backmap(
            "{#A=[$]C([$])C[$],#B=[$]CC[$]}", placement="rdkit"
        )
        # each bead's atoms are centered on the bead coordinate
        for i, child in enumerate(compound.children):
            centroid = np.mean([p.pos for p in child.particles()], axis=0)
            assert np.linalg.norm(centroid - branched_path.coordinates[i]) < 0.35

    def test_backmap_fragname_map(self):
        path = straight_line(spacing=0.35, N=4)  # default bead name "_A"
        compound = path.backmap("{#PE=[>]CC[<]}", fragname_map={"_A": "PE"})
        assert compound.n_particles == 26
        assert compound.n_bonds == 25
        # children keep the original bead names
        assert [c.name for c in compound.children] == ["_A"] * 4

    def test_backmap_elements(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        for particle in compound.particles():
            assert particle.element is not None
            assert particle.element.symbol == particle.name

    def test_backmap_bad_placement(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PE")
        with pytest.raises(ValueError):
            path.backmap("{#PE=[>]CC[<]}", placement="bad-option")

    def test_validate_missing_descriptors(self, branched_path):
        # A only has two descriptors but the branch point has degree 3;
        # cgsmiles silently skips the branch bond without validation
        fragments = "{#A=[>]CC[<],#B=[>]CC[<]}"
        with pytest.raises(ValueError, match="bonding descriptors"):
            branched_path.backmap(fragments)
        compound = branched_path.backmap(fragments, validate=False)
        assert not nx.is_connected(compound_graph(compound))

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=5, bead_name="PEO")
        template = mb.load("COC", smiles=True)
        compound = path.backmap(
            "{#PEO=[>]COC[<]}",
            templates={"PEO": template},
            placement="template",
        )
        assert compound.n_particles == 37
        assert compound.n_bonds == 36
        assert nx.is_connected(compound_graph(compound))
        for i, child in enumerate(compound.children):
            centroid = np.mean([p.pos for p in child.particles()], axis=0)
            assert np.linalg.norm(centroid - path.coordinates[i]) < 0.35

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_branched(self, branched_path):
        import mbuild as mb

        compound = branched_path.backmap(
            "{#A=[$]C([$])C[$],#B=[$]CC[$]}",
            templates={
                "A": mb.load("CC", smiles=True),
                "B": mb.load("CC", smiles=True),
            },
            placement="template",
        )
        assert nx.is_connected(compound_graph(compound))
        assert [c.name for c in compound.children] == list(branched_path.beads)

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_partial(self, branched_path):
        import mbuild as mb

        # A beads use the template, B beads fall back to rdkit embedding
        compound = branched_path.backmap(
            "{#A=[$]C([$])C[$],#B=[$]CC[$]}",
            templates={"A": mb.load("CC", smiles=True)},
        )
        assert nx.is_connected(compound_graph(compound))

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_mismatch(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=3, bead_name="PEO")
        # wrong heavy atom count: template CC vs fragment COC
        with pytest.raises(ValueError, match="heavy atoms"):
            path.backmap(
                "{#PEO=[>]COC[<]}", templates={"PEO": mb.load("CC", smiles=True)}
            )
        # right count, wrong element order: CCO vs COC
        with pytest.raises(ValueError, match="same order"):
            path.backmap(
                "{#PEO=[>]COC[<]}", templates={"PEO": mb.load("CCO", smiles=True)}
            )

    def test_backmap_template_placement_requires_templates(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PE")
        with pytest.raises(ValueError, match="requires"):
            path.backmap("{#PE=[>]CC[<]}", placement="template")

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_tagged(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        # template SMILES written in a different atom order than the
        # fragment (COC); tags pin each heavy atom to its fragment index
        tagged = mb.load("O{1}(C{0})C{2}", smiles=True)
        compound = path.backmap(
            "{#PEO=[>]COC[<]}", templates={"PEO": tagged}, placement="template"
        )
        assert compound.n_particles == 30
        assert nx.is_connected(compound_graph(compound))
        # the same reordered template without tags fails the order check
        with pytest.raises(ValueError, match="same order"):
            path.backmap(
                "{#PEO=[>]COC[<]}",
                templates={"PEO": mb.load("O(C)C", smiles=True)},
                placement="template",
            )

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_bad_tags(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=3, bead_name="PEO")
        with pytest.raises(ValueError, match="all heavy atoms or none"):
            path.backmap(
                "{#PEO=[>]COC[<]}",
                templates={"PEO": mb.load("C{0}OC", smiles=True)},
            )
        with pytest.raises(ValueError, match="duplicate"):
            path.backmap(
                "{#PEO=[>]COC[<]}",
                templates={"PEO": mb.load("C{0}O{1}C{1}", smiles=True)},
            )

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_compound_to_fragment_graph(self):
        import mbuild as mb
        from mbuild.coarse_graining import compound_to_fragment_graph

        template = mb.load("C{<}C{>,$br}", smiles=True)
        graph = compound_to_fragment_graph(template, "PE")
        assert len(graph) == 2  # hydrogens folded into hcount
        assert graph.nodes[0]["element"] == "C"
        assert graph.nodes[0]["hcount"] == 3
        assert graph.nodes[0]["bonding"] == ["<1"]
        # comma separates multiple descriptors on one atom
        assert graph.nodes[1]["bonding"] == [">1", "$br1"]
        assert graph.edges[(0, 1)]["order"] == 1

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_compound_defined_fragment(self):
        import mbuild as mb

        # no CGsmiles fragment string at all: the tagged compound
        # defines the fragment (equivalent to "{#PE=[<]CC[>]}")
        path = straight_line(spacing=0.35, N=4, bead_name="PE")
        template = mb.load("C{<}C{>}", smiles=True)
        compound = path.backmap(templates={"PE": template})
        assert compound.n_particles == 26
        assert compound.n_bonds == 25
        assert nx.is_connected(compound_graph(compound))
        # template coordinates are used for placement
        for i, child in enumerate(compound.children):
            centroid = np.mean([p.pos for p in child.particles()], axis=0)
            assert np.linalg.norm(centroid - path.coordinates[i]) < 1e-6

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_compound_defined_branched(self, branched_path):
        import mbuild as mb

        templates = {
            "A": mb.load("C{$,$}C{$}", smiles=True),
            "B": mb.load("C{$}C{$}", smiles=True),
        }
        compound = branched_path.backmap(templates=templates)
        assert nx.is_connected(compound_graph(compound))
        assert [c.name for c in compound.children] == list(branched_path.beads)

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_compound_defined_aromatic(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=3, bead_name="PS")
        template = mb.load("C{>}C{<}(c1ccccc1)", smiles=True)
        compound = path.backmap(templates={"PS": template})
        assert nx.is_connected(compound_graph(compound))
        carbons = sum(1 for p in compound.particles() if p.name == "C")
        assert carbons == 3 * 8

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_mixed_string_and_compound_definitions(self, branched_path):
        import mbuild as mb

        # A defined by the string, B defined by a tagged compound
        compound = branched_path.backmap(
            "{#A=[$]C([$])C[$]}",
            templates={"B": mb.load("C{$}C{$}", smiles=True)},
        )
        assert nx.is_connected(compound_graph(compound))

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_compound_defined_untagged(self):
        import mbuild as mb

        # a compound-defined fragment without descriptor tags cannot
        # bond to its neighbors; validation reports the unrealized bonds
        path = straight_line(spacing=0.35, N=3, bead_name="PE")
        with pytest.raises(ValueError, match="bonding descriptors"):
            path.backmap(templates={"PE": mb.load("CC", smiles=True)})

    def test_backmap_requires_fragments_or_templates(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PE")
        with pytest.raises(ValueError, match="fragment string"):
            path.backmap()

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_backmap_template_non_integer_tags_ignored(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=3, bead_name="PE")
        # head/tail tags from the Polymer workflow fall back to order matching
        template = mb.load("C{head}C{tail}", smiles=True)
        compound = path.backmap(
            "{#PE=[>]CC[<]}", templates={"PE": template}, placement="template"
        )
        assert nx.is_connected(compound_graph(compound))

    def test_backmap_missing_fragment(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PE")
        with pytest.raises(SyntaxError):
            path.backmap("{#PEO=[>]COC[<]}")

    def test_backmap_from_compound_path(self, ethane):
        # Paths created from compounds can also be used as meta graphs
        path = Path.from_compound(ethane)
        meta_graph = path.to_cgsmiles_graph()
        assert len(meta_graph) == ethane.n_particles

    @pytest.fixture
    def multi_chain_path(self):
        """Three disconnected chains in one Path, the third one branched."""
        path = hard_sphere_random_walk(
            termination=8, radius=0.2, bond_length=0.25, bead_name="A", seed=1
        )
        path = hard_sphere_random_walk(
            path=path,
            termination=16,
            radius=0.2,
            bond_length=0.25,
            bead_name="A",
            seed=2,
        )
        path = hard_sphere_random_walk(
            path=path,
            termination=24,
            radius=0.2,
            bond_length=0.25,
            bead_name="A",
            seed=3,
        )
        path = hard_sphere_random_walk(
            path=path,
            termination=28,
            radius=0.2,
            bond_length=0.25,
            bead_name="B",
            initial_point=20,
            connectivity="link-linear",
            seed=4,
        )
        assert nx.number_connected_components(path.bond_graph) == 3
        return path

    def test_to_cgsmiles_multi_chain(self, multi_chain_path):
        from cgsmiles import read_cgsmiles

        cgsmiles_str = multi_chain_path.to_cgsmiles()
        n_beads = len(multi_chain_path.coordinates)
        assert cgsmiles_str.count("[#") == n_beads
        assert "." in cgsmiles_str
        # '.' reads back as zero-order (non-bonded) edges
        read_back = read_cgsmiles(cgsmiles_str)
        assert read_back.number_of_nodes() == n_beads
        orders = nx.get_edge_attributes(read_back, "order")
        bonded = read_back.edge_subgraph(
            [edge for edge, order in orders.items() if order != 0]
        )
        assert nx.number_connected_components(bonded) == 3

    def test_backmap_multi_chain(self, multi_chain_path):
        compound = multi_chain_path.backmap("{#A=[$]C([$])C[$],#B=[$]CC[$]}")
        # root -> Molecule compounds -> bead compounds
        assert len(compound.children) == 3
        assert [c.name for c in compound.children] == ["Molecule"] * 3
        n_beads = sum(len(mol.children) for mol in compound.children)
        assert n_beads == len(multi_chain_path.coordinates)
        # one atomistic molecule per CG component
        assert nx.number_connected_components(compound_graph(compound)) == 3
        # bead compounds keep bead names, grouped by their own chain
        for mol, component in zip(
            compound.children,
            sorted(
                (
                    sorted(c)
                    for c in nx.connected_components(multi_chain_path.bond_graph)
                ),
                key=lambda c: c[0],
            ),
        ):
            assert [b.name for b in mol.children] == [
                str(multi_chain_path.beads[i]) for i in component
            ]

    # --- CG Compound input -------------------------------------------------

    @staticmethod
    def cg_compound_from_path(path):
        """Build a CG Compound with one leaf particle per path bead."""
        import mbuild as mb

        compound = mb.Compound(name="CG")
        beads = [
            mb.Compound(name=str(name), pos=xyz)
            for name, xyz in zip(path.beads, path.coordinates)
        ]
        compound.add(beads)
        for u, v in path.bond_graph.edges:
            compound.add_bond((beads[int(u)], beads[int(v)]))
        return compound

    def test_compound_to_cgsmiles_graph(self, branched_path):
        from mbuild.coarse_graining import to_cgsmiles_graph

        cg_compound = self.cg_compound_from_path(branched_path)
        from_compound = to_cgsmiles_graph(cg_compound)
        from_path = branched_path.to_cgsmiles_graph()
        assert nx.utils.graphs_equal(
            nx.Graph(from_path.edges), nx.Graph(from_compound.edges)
        )
        for node in from_path.nodes:
            assert (
                from_compound.nodes[node]["fragname"]
                == from_path.nodes[node]["fragname"]
            )
            assert np.allclose(
                from_compound.nodes[node]["position"],
                from_path.nodes[node]["position"],
            )

    def test_compound_to_cgsmiles_string(self, branched_path):
        cg_compound = self.cg_compound_from_path(branched_path)
        assert cg_compound.to_cgsmiles() == branched_path.to_cgsmiles()

    def test_backmap_cg_compound_matches_path(self, branched_path):
        cg_compound = self.cg_compound_from_path(branched_path)
        via_path = branched_path.backmap("{#A=[>bb]C([<br])C[<bb],#B=[>br]CC[<br]}")
        via_compound = cg_compound.backmap("{#A=[>bb]C([<br])C[<bb],#B=[>br]CC[<br]}")
        assert via_compound.n_particles == via_path.n_particles
        assert via_compound.n_bonds == via_path.n_bonds
        assert [c.name for c in via_compound.children] == [
            c.name for c in via_path.children
        ]
        assert np.allclose(via_compound.xyz, via_path.xyz)

    def test_backmap_cg_compound_fragname_map(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=4)  # default bead name "_A"
        cg_compound = self.cg_compound_from_path(path)
        compound = mb.backmap(
            cg_compound,
            "{#PE=[>]CC[<]}",
            fragname_map={"_A": "PE"},
        )
        assert compound.n_particles == 26  # C8H18
        # bead compounds keep the original particle names
        assert [c.name for c in compound.children] == ["_A"] * 4

    def test_to_cgsmiles_graph_bad_type(self):
        from mbuild.coarse_graining import to_cgsmiles_graph

        with pytest.raises(TypeError):
            to_cgsmiles_graph([0, 1, 2])

    def test_to_cgsmiles_graph_single_particle(self):
        # a childless Compound is itself a particle -> a one-bead graph
        import mbuild as mb
        from mbuild.coarse_graining import to_cgsmiles_graph

        graph = to_cgsmiles_graph(mb.Compound(name="A", pos=[1.0, 2.0, 3.0]))
        assert len(graph) == 1
        assert graph.nodes[0]["fragname"] == "A"
        assert np.allclose(graph.nodes[0]["position"], [1.0, 2.0, 3.0])

    def test_coarse_grain_round_trip_linear(self):
        from mbuild.coarse_graining import coarse_grain

        path = straight_line(spacing=0.35, N=6, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        cg_path, mapping = coarse_grain(compound, "{#PEO=[>]COC[<]}")
        assert nx.is_isomorphic(cg_path.bond_graph, path.bond_graph)
        assert list(cg_path.beads) == list(path.beads)
        assert np.allclose(cg_path.coordinates, path.coordinates, atol=0.01)
        assert mapping.n_beads == 6
        assert len(mapping.atom_to_bead) == compound.n_particles
        # bead_to_atoms is the inverse of atom_to_bead
        for bead, atoms in mapping.bead_to_atoms.items():
            assert all(mapping.atom_to_bead[a] == bead for a in atoms)

    def test_coarse_grain_round_trip_branched(self, branched_path):
        # chemically distinct fragments -> bead labels are recovered exactly
        fragments = "{#A=[>bb]CC(c1ccccc1)([<bb])[<br],#B=[>br]C[<br]}"
        compound = branched_path.backmap(
            fragments,
        )
        cg_path, _ = compound.coarse_grain(fragments)
        assert nx.is_isomorphic(cg_path.bond_graph, branched_path.bond_graph)
        assert list(cg_path.beads) == list(branched_path.beads)
        assert np.allclose(cg_path.coordinates, branched_path.coordinates, atol=0.01)

    def test_coarse_grain_identical_fragments_ambiguous_names(self, branched_path):
        # A and B are both CC (only descriptor labels differ): the tiling
        # frame is recovered, but the A/B labeling is inherently ambiguous
        fragments = "{#A=[>bb]C([<br])C[<bb],#B=[>br]CC[<br]}"
        compound = branched_path.backmap(
            fragments,
        )
        cg_path, _ = compound.coarse_grain(fragments)
        assert nx.is_isomorphic(cg_path.bond_graph, branched_path.bond_graph)
        assert np.allclose(cg_path.coordinates, branched_path.coordinates, atol=0.01)

    def test_coarse_grain_multi_molecule(self):
        import mbuild as mb

        water = mb.load("O", smiles=True)
        box = mb.Compound()
        for i in range(3):
            clone = mb.clone(water)
            clone.translate([float(i), 0, 0])
            box.add(clone)
        cg_path, mapping = mb.coarse_grain(box, "{#W=O}")
        assert len(cg_path.coordinates) == 3
        assert cg_path.bond_graph.number_of_edges() == 0
        assert set(cg_path.beads) == {"W"}
        assert mapping.n_beads == 3

    def test_coarse_grain_hierarchy(self):
        path = straight_line(spacing=0.35, N=5, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        cg_path, mapping = compound.coarse_grain(beads=["PEO"])
        assert nx.is_isomorphic(cg_path.bond_graph, path.bond_graph)
        assert list(cg_path.beads) == ["PEO"] * 5
        assert mapping.n_beads == 5

    def test_coarse_grain_hierarchy_unmapped_raises(self):
        path = straight_line(spacing=0.35, N=5, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        with pytest.raises(ValueError, match="not contained in any"):
            compound.coarse_grain(beads=["NotABead"])

    def test_coarse_grain_explicit_mapping(self):
        import mbuild as mb

        ethane = mb.load("CC", smiles=True)
        carbons = [p for p in ethane.particles() if p.name == "C"]
        mapping = {}
        for particle in ethane.particles():
            if particle.name == "C":
                mapping[particle] = carbons.index(particle)
        for p1, p2 in ethane.bonds():
            if p1.name == "H":
                mapping[p1] = mapping[p2]
            elif p2.name == "H":
                mapping[p2] = mapping[p1]
        cg_path, cg_mapping = ethane.coarse_grain(
            mapping=mapping, bead_names=["CH3", "CH3"]
        )
        assert len(cg_path.coordinates) == 2
        assert cg_path.bond_graph.number_of_edges() == 1
        assert list(cg_path.beads) == ["CH3", "CH3"]
        assert sorted(len(v) for v in cg_mapping.bead_to_atoms.values()) == [4, 4]

    def test_coarse_grain_wrong_fragment_raises(self):
        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        with pytest.raises(ValueError, match="matches"):
            compound.coarse_grain("{#PE=[>]CC[<]}")

    def test_coarse_grain_missing_descriptors_raises(self):
        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        with pytest.raises(ValueError):
            compound.coarse_grain("{#PEO=COC}")

    @pytest.mark.skipif(not has_rdkit, reason="rdkit is not installed")
    def test_coarse_grain_template_fragment(self):
        import mbuild as mb

        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        template = mb.load("C{>}O{}C{<}", smiles=True)
        cg_path, _ = compound.coarse_grain(templates={"PEO": template})
        assert nx.is_isomorphic(cg_path.bond_graph, path.bond_graph)
        assert list(cg_path.beads) == ["PEO"] * 4

    def test_coarse_grain_mass_center(self):
        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        geo, _ = compound.coarse_grain("{#PEO=[>]COC[<]}", center="geometry")
        mass, _ = compound.coarse_grain("{#PEO=[>]COC[<]}", center="mass")
        assert not np.allclose(geo.coordinates, mass.coordinates)

    def test_coarse_grain_mode_validation(self):
        path = straight_line(spacing=0.35, N=3, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        with pytest.raises(ValueError, match="exactly one"):
            compound.coarse_grain()
        with pytest.raises(ValueError, match="exactly one"):
            compound.coarse_grain("{#PEO=[>]COC[<]}", beads=["PEO"])
        with pytest.raises(ValueError, match="invalid"):
            compound.coarse_grain("{#PEO=[>]COC[<]}", center="bad")

    def test_coarse_grain_then_backmap(self):
        path = straight_line(spacing=0.35, N=5, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        cg_path, _ = compound.coarse_grain("{#PEO=[>]COC[<]}")
        again = cg_path.backmap("{#PEO=[>]COC[<]}")
        assert again.n_particles == compound.n_particles
        assert again.n_bonds == compound.n_bonds

    def test_coarse_grain_path_to_compound(self):
        # Paths built by the constructor (e.g. coarse_grain output) have no
        # node attributes on their bond graph; to_compound must use the
        # beads array for names
        path = straight_line(spacing=0.35, N=4, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}")
        cg_path, _ = compound.coarse_grain("{#PEO=[>]COC[<]}")
        cg_compound = cg_path.to_compound()
        assert [p.name for p in cg_compound.particles()] == ["PEO"] * 4
        assert cg_compound.n_bonds == 3
        assert np.allclose(
            np.array([p.pos for p in cg_compound.particles()]),
            cg_path.coordinates,
        )

    def test_backmap_hand_built_template_manual_tags(self):
        # a template built by hand (no SMILES): descriptor tags set via the
        # particle_tag attribute, bonds added with add_bond (which stores
        # bond_order=0.0 -> must be read as an unspecified single bond)
        import mbuild as mb
        from mbuild.coarse_graining import compound_to_fragment_graph

        frag = mb.Compound(name="Dimer")
        c1 = mb.Compound(name="C", pos=[0.00, 0, 0])
        c2 = mb.Compound(name="C", pos=[0.15, 0, 0])
        frag.add([c1, c2])
        frag.add_bond((c1, c2))
        for parent, dx in [(c1, 0.0), (c2, 0.15)]:
            for dy, dz in [(0.1, 0.05), (-0.1, 0.05), (0.0, -0.1)]:
                h = mb.Compound(name="H", pos=[dx, dy, dz])
                frag.add(h)
                frag.add_bond((parent, h))
        c1.particle_tag = "<"
        c2.particle_tag = ">"

        graph = compound_to_fragment_graph(frag, "A")
        assert list(graph.edges(data="order")) == [(0, 1, 1)]

        path = straight_line(spacing=0.3, N=4, bead_name="A")
        compound = path.backmap(templates={"A": frag})
        assert compound.n_particles == 26  # C8H18
        assert compound.n_bonds == 25
        assert nx.is_connected(compound_graph(compound))

    def test_backmap_bead_to_dimer_to_water(self):
        # multi-resolution to a molecule: each W bead resolves to a two-bead
        # dimer (CG), then to a water molecule (all-atom). The last block is
        # SMILES, so an all-atom endpoint is inferred.
        import mbuild as mb

        box = mb.fill_box(compound=mb.Compound(name="W"), n_compounds=5, box=[3, 3, 3])
        water = box.backmap(["{#W=[#HH][#OH]}", "{#HH=[$][H],#OH=[$]O}"])
        # 5 W -> 5 H2O: O with two H, two bonds each
        assert water.n_particles == 15
        assert water.n_bonds == 10
        names = sorted(p.name for p in water.particles())
        assert names.count("O") == 5 and names.count("H") == 10
        # the waters are separate molecules (the box beads are unbonded)
        assert nx.number_connected_components(compound_graph(water)) == 5

    def test_backmap_cg_chain_refine_then_atomistic(self):
        # a CG chain of 4 E2 beads refines to 8 E beads (CG endpoint), then
        # all the way to an atomistic polyethylene chain of 8 CC units.
        path = straight_line(spacing=0.5, N=4, bead_name="E2")

        # CG -> CG: 4 E2 beads -> 8 E beads (all_atom=False inferred)
        cg_chain = path.backmap("{#E2=[>][#E][#E][<]}")
        assert cg_chain.n_particles == 8
        assert cg_chain.n_bonds == 7
        assert {p.name for p in cg_chain.particles()} == {"E"}
        assert all(p.element is None for p in cg_chain.particles())
        assert nx.is_connected(compound_graph(cg_chain))

        # CG -> CG -> atomistic in one multi-resolution call
        pe = path.backmap(["{#E2=[>][#E][#E][<]}", "{#E=[>]CC[<]}"])
        assert pe.n_particles == 50  # 16 C + 34 H
        assert pe.n_bonds == 49
        names = [p.name for p in pe.particles()]
        assert names.count("C") == 16 and names.count("H") == 34
        assert nx.is_connected(compound_graph(pe))

    def test_backmap_cg_chain_direct_to_atomistic(self):
        # the same polyethylene chain reached directly: each E2 bead resolves
        # to four carbons in one step, skipping the intermediate CG level.
        path = straight_line(spacing=0.5, N=4, bead_name="E2")

        direct = path.backmap("{#E2=[>]CCCC[<]}")
        assert direct.n_particles == 50
        assert direct.n_bonds == 49
        assert nx.is_connected(compound_graph(direct))

        # the direct and refined routes give the same atomistic chain
        refined = path.backmap(["{#E2=[>][#E][#E][<]}", "{#E=[>]CC[<]}"])
        assert direct.n_particles == refined.n_particles
        assert direct.n_bonds == refined.n_bonds

    def test_coarse_grain_zero_mass_raises(self):
        import mbuild as mb

        chain = mb.Compound(name="CG")
        beads = [mb.Compound(name="A", pos=[0.3 * i, 0, 0]) for i in range(3)]
        chain.add(beads)
        for b1, b2 in zip(beads[:-1], beads[1:]):
            chain.add_bond((b1, b2))
        # CG beads carry no element, so their mass is 0/None, not usable
        with pytest.raises(ValueError, match="mass"):
            chain.coarse_grain(beads=["A"], center="mass")
