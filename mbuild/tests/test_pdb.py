import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.formats.pdb import _decode_index, _encode_index
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn


def _records(filename, kind):
    with open(filename) as f:
        return [line for line in f if line.startswith(kind)]


def _residue(name, number, atom_names, x):
    residue = mb.Compound(name=name)
    for i, atom_name in enumerate(atom_names):
        residue.add(
            mb.Particle(name=atom_name, element=atom_name[0], pos=[x, 0.1 * i, 0])
        )
    residue.residue_number = number
    return residue


class TestPDB(BaseTest):
    def test_save_writes_bonds(self, ethane):
        ethane.save("ethane.pdb")
        assert len(_records("ethane.pdb", "HETATM")) == 8
        assert len(_records("ethane.pdb", "CONECT")) == 8
        assert mb.load("ethane.pdb").n_bonds == 7

    def test_default_residue(self, ethane):
        ethane.save("ethane.pdb")
        assert {line[17:21] for line in _records("ethane.pdb", "HETATM")} == {"RES "}

    def test_residues_by_name(self, ethane, h2o):
        system = mb.Compound([ethane, h2o])
        system.save("system.pdb", residues=["Ethane", "H2O"])
        atoms = _records("system.pdb", "HETATM")
        assert [line[17:21] for line in atoms] == ["Etha"] * 8 + ["H2O "] * 3
        assert [int(line[22:26]) for line in atoms] == [1] * 8 + [2] * 3

    def test_residue_and_chain_fields(self):
        residue = _residue("LIG", 52, ["C1"], 0.0)
        residue.insertion_code = "A"
        chain = mb.Compound(subcompounds=residue)
        chain.chain_id = "B"
        chain.save("chain.pdb")
        atom = _records("chain.pdb", "HETATM")[0]
        assert atom[17:21] == "LIG "
        assert atom[21] == "B"
        assert atom[22:27] == "  52A"
        ter = _records("chain.pdb", "TER")[0]
        assert ter[21:27] == "B  52A"

    def test_standard_residues_skip_implied_bonds(self):
        ala = _residue("ALA", 1, ["N", "CA", "C"], 0.0)
        gly = _residue("GLY", 2, ["N", "CA", "C"], 0.4)
        cys = _residue("CYS", 3, ["N", "CA", "SG"], 0.8)
        chain = mb.Compound([ala, gly, cys])
        chain.chain_id = "A"
        for res in (ala, gly, cys):
            chain.add_bond((res[0], res[1]))
        chain.add_bond((ala[2], gly[0]))  # peptide bond, implied
        chain.add_bond((ala[1], cys[2]))  # not implied by the residues
        chain.save("peptide.pdb")
        assert len(_records("peptide.pdb", "ATOM")) == 9
        conect = _records("peptide.pdb", "CONECT")
        assert [line.split() for line in conect] == [
            ["CONECT", "2", "9"],
            ["CONECT", "9", "2"],
        ]

    def test_formal_charge_column(self):
        sodium = mb.Particle(name="NA", element="Na", formal_charge=1)
        chloride = mb.Particle(
            name="CL", element="Cl", pos=[0.3, 0, 0], formal_charge=-1
        )
        mb.Compound([sodium, chloride]).save("ions.pdb")
        atoms = _records("ions.pdb", "HETATM")
        assert [line[76:80] for line in atoms] == ["NA1+", "CL1-"]
        # two-letter element symbols start the atom name in column 13
        assert [line[12:16] for line in atoms] == ["NA  ", "CL  "]

    def test_bond_order_in_conect(self):
        carbonyl = mb.Compound(name="FOR")
        carbon = mb.Particle(name="C", element="C")
        oxygen = mb.Particle(name="O", element="O", pos=[0.12, 0, 0])
        carbonyl.add([carbon, oxygen])
        carbonyl.add_bond((carbon, oxygen), bond_order=2.0)
        carbonyl.save("carbonyl.pdb")
        assert [line.split() for line in _records("carbonyl.pdb", "CONECT")] == [
            ["CONECT", "1", "2", "2"],
            ["CONECT", "2", "1", "1"],
        ]

    def test_box(self, ethane):
        ethane.box = mb.Box([2.0, 3.0, 4.0])
        ethane.save("box.pdb")
        cryst1 = _records("box.pdb", "CRYST1")[0]
        assert np.allclose([float(v) for v in cryst1.split()[1:4]], [20, 30, 40])

    def test_include_ports(self, ch2):
        ch2.save("ports.pdb", include_ports=True)
        names = [line[17:20] for line in _records("ports.pdb", "HETATM")]
        n_port_particles = len(list(ch2.particles(include_ports=True))) - (
            ch2.n_particles
        )
        assert names.count("PRT") == n_port_particles > 0

    def test_long_chain_id(self, ethane):
        ethane.chain_id = "AB"
        with pytest.raises(MBuildError):
            ethane.save("long.pdb")

    def test_unknown_keyword(self, ethane):
        with pytest.raises(TypeError):
            ethane.save("ethane.pdb", renumber=False)

    def test_large_serials(self):
        assert _encode_index(99999, 5) == "99999"
        assert _encode_index(100000, 5) == "A0000"
        with pytest.raises(MBuildError):
            _encode_index(10**5 + 6 * 16**4, 5)


class TestReadPDB(BaseTest):
    def test_round_trip(self, ethane):
        ethane.save("ethane.pdb")
        loaded = mb.load("ethane.pdb", backend="internal")
        assert loaded.n_particles == 8
        assert loaded.n_bonds == 7
        assert np.allclose(loaded.xyz, ethane.xyz, atol=1e-4)
        assert [p.element.symbol for p in loaded.particles()] == [
            p.element.symbol for p in ethane.particles()
        ]

    def test_bond_orders_and_charges(self):
        carbonyl = mb.Compound(name="FOR")
        carbon = mb.Particle(name="C", element="C")
        oxygen = mb.Particle(name="O", element="O", pos=[0.12, 0, 0], formal_charge=-1)
        carbonyl.add([carbon, oxygen])
        carbonyl.add_bond((carbon, oxygen), bond_order=2.0)
        carbonyl.save("carbonyl.pdb")
        loaded = mb.load("carbonyl.pdb", backend="internal")
        assert [b[2]["bond_order"] for b in loaded.bonds(return_bond_order=True)] == [
            2.0
        ]
        assert [p.formal_charge for p in loaded.particles()] == [None, -1]

    def test_residues_and_chains(self):
        chains = []
        for chain_id, x in (("A", 0.0), ("B", 1.0)):
            residue = _residue("LIG", 52, ["C1", "C2"], x)
            residue.insertion_code = "A"
            chain = mb.Compound(subcompounds=residue)
            chain.chain_id = chain_id
            chains.append(chain)
        mb.Compound(chains).save("chains.pdb")
        loaded = mb.load("chains.pdb", backend="internal")
        assert list(loaded.labels) == ["all-chains", "chain[0]", "chain[1]"]
        assert [chain.chain_id for chain in loaded.children] == ["A", "B"]
        residue = loaded.children[1].children[0]
        assert (residue.name, residue.residue_number, residue.insertion_code) == (
            "LIG",
            52,
            "A",
        )
        assert [p.name for p in residue.particles()] == ["C1", "C2"]

    def test_protein(self):
        protein = mb.load(get_fn("6m03.pdb"), backend="internal")
        reference = mb.load(get_fn("6m03.pdb"))
        assert protein.n_particles == reference.n_particles
        assert [c.chain_id for c in protein.children] == ["A", "A"]
        residues = protein.children[0].children
        assert (residues[0].name, residues[0].residue_number) == ("SER", 1)
        assert [r.name for r in residues] == [
            r.name for r in reference.children[0].children
        ]

    def test_alternate_locations(self, caplog):
        with open("altloc.pdb", "w") as f:
            f.write(
                "ATOM      1  CA ASER A   1       0.000   0.000   0.000  0.60  0.00           C\n"
                "ATOM      2  CA BSER A   1       1.000   0.000   0.000  0.40  0.00           C\n"
                "ATOM      3  CB ASER A   1       0.000   1.000   0.000  0.60  0.00           C\n"
                "ATOM      4  CB BSER A   1       1.000   1.000   0.000  0.40  0.00           C\n"
                "END\n"
            )
        loaded = mb.load("altloc.pdb", backend="internal")
        assert loaded.n_particles == 2
        assert np.allclose(loaded.xyz[:, 0], 0.0)
        assert "first alternate location" in caplog.text

    def test_infer_hierarchy_false(self, ethane):
        ethane.save("ethane.pdb")
        loaded = mb.load("ethane.pdb", backend="internal", infer_hierarchy=False)
        assert len(loaded.children) == loaded.n_particles == 8

    def test_coords_only(self, ethane):
        moved = mb.clone(ethane)
        moved.translate([1.0, 0, 0])
        moved.save("moved.pdb")
        mb.load("moved.pdb", compound=ethane, coords_only=True, backend="internal")
        assert np.allclose(ethane.xyz, moved.xyz, atol=1e-4)

    def test_large_serials(self):
        assert _decode_index("A0000", width=5) == 100000
        assert _decode_index(_encode_index(123456, 5), width=5) == 123456
        assert _decode_index("  42", width=4) == 42
