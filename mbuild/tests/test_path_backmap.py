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


class TestPathBackmap(BaseTest):
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

    def test_resolve_path(self):
        from mbuild.path.backmap import resolve_path

        path = straight_line(spacing=0.25, N=3, bead_name="PEO")
        meta_graph, molecule, node_to_beads = resolve_path(path, "{#PEO=[>]COC[<]}")
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
        compound = path.backmap("{#PEO=[>]COC[<]}", placement="jitter")
        assert compound.n_particles == 37
        assert compound.n_bonds == 36
        assert nx.is_connected(compound_graph(compound))
        # hierarchy mirrors the path: one child compound per bead
        assert [c.name for c in compound.children] == ["PEO"] * 5

    def test_backmap_branched(self, branched_path):
        n_beads = len(branched_path.coordinates)
        degrees = dict(branched_path.bond_graph.degree())
        assert max(degrees.values()) == 3  # the case build_from_path can't handle
        compound = branched_path.backmap(
            "{#A=[$]C([$])C[$],#B=[$]CC[$]}", placement="jitter"
        )
        assert nx.is_connected(compound_graph(compound))
        assert len(compound.children) == n_beads
        assert [c.name for c in compound.children] == list(branched_path.beads)

    def test_backmap_jitter_positions(self):
        path = straight_line(spacing=0.35, N=4, bead_name="PE")
        compound = path.backmap("{#PE=[>]CC[<]}", placement="jitter", jitter=0.01)
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
        compound = path.backmap(
            "{#PE=[>]CC[<]}", fragname_map={"_A": "PE"}, placement="jitter"
        )
        assert compound.n_particles == 26
        assert compound.n_bonds == 25
        # children keep the original bead names
        assert [c.name for c in compound.children] == ["_A"] * 4

    def test_backmap_elements(self):
        path = straight_line(spacing=0.35, N=2, bead_name="PEO")
        compound = path.backmap("{#PEO=[>]COC[<]}", placement="jitter")
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
            branched_path.backmap(fragments, placement="jitter")
        compound = branched_path.backmap(
            fragments, placement="jitter", validate=False
        )
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
            templates={"A": mb.load("CC", smiles=True), "B": mb.load("CC", smiles=True)},
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

        from mbuild.path.backmap import compound_to_fragment_graph

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
            path=path, termination=16, radius=0.2, bond_length=0.25,
            bead_name="A", seed=2,
        )
        path = hard_sphere_random_walk(
            path=path, termination=24, radius=0.2, bond_length=0.25,
            bead_name="A", seed=3,
        )
        path = hard_sphere_random_walk(
            path=path, termination=28, radius=0.2, bond_length=0.25,
            bead_name="B", initial_point=20, connectivity="link-linear", seed=4,
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
        compound = multi_chain_path.backmap(
            "{#A=[$]C([$])C[$],#B=[$]CC[$]}", placement="jitter"
        )
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
