import numpy as np
import pytest

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.formats.pdb import _encode_index, pdb_residues
from mbuild.tests.base_test import BaseTest


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

    def test_conect_lists_every_bond(self):
        ala = _residue("ALA", 1, ["N", "CA", "C"], 0.0)
        gly = _residue("GLY", 2, ["N", "CA", "C"], 0.4)
        cys = _residue("CYS", 3, ["N", "CA", "SG"], 0.8)
        chain = mb.Compound([ala, gly, cys])
        chain.chain_id = "A"
        for res in (ala, gly, cys):
            chain.add_bond((res[0], res[1]))
        chain.add_bond((ala[2], gly[0]))
        chain.add_bond((ala[1], cys[2]))
        chain.save("peptide.pdb")
        assert len(_records("peptide.pdb", "ATOM")) == 9
        pairs = set()
        for line in _records("peptide.pdb", "CONECT"):
            serial, *partners = (int(field) for field in line.split()[1:])
            pairs.update(tuple(sorted((serial, p))) for p in partners)
        assert pairs == {(1, 2), (4, 5), (7, 8), (3, 4), (2, 9)}

    def test_residue_named_like_a_nucleotide_keeps_bonds(self):
        monomer = _residue("A", 1, ["C1", "C2", "O1"], 0.0)
        monomer.add_bond((monomer[0], monomer[1]))
        monomer.add_bond((monomer[1], monomer[2]))
        monomer.save("named_a.pdb")
        assert len(_records("named_a.pdb", "CONECT")) == 3

    def test_pdb_residues_match_the_file(self, ethane):
        ethane.name = "ETHANE"
        ethane.residue_number = 7
        ethane.save("ethane.pdb")
        (residue,) = pdb_residues(ethane)
        atoms = _records("ethane.pdb", "HETATM")
        assert residue.name == atoms[0][17:21].strip() == "ETHA"
        assert residue.original_name == "ETHANE"
        assert residue.number == int(atoms[0][22:26]) == 7
        assert residue.insertion_code == residue.chain_id == ""
        assert [residue.atom_names[p] for p in residue.particles] == [
            line[12:16].strip() for line in atoms
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

    def test_auto_numbers_skip_explicit_numbers(self):
        ala = mb.Compound(name="ALA")
        ala.residue_number = 2
        lig = mb.Compound(name="LIG")
        lig.add(mb.Particle(name="C1", element="C", pos=[0.3, 0, 0]))
        chain = mb.Compound([ala, lig])
        chain.chain_id = "A"
        chain.save("numbers.pdb", residues=["LIG"])
        atoms = _records("numbers.pdb", "ATOM") + _records("numbers.pdb", "HETATM")
        assert [(line[17:20], int(line[22:26])) for line in atoms] == [
            ("ALA", 2),
            ("LIG", 3),
        ]

    def test_auto_numbers_restart_in_each_chain(self):
        chains = []
        for chain_id in "AB":
            chain = mb.Compound([_residue("LIG", None, ["C1"], 0.0)])
            chain.chain_id = chain_id
            chains.append(chain)
        mb.Compound(chains).save("chains.pdb", residues=["LIG"])
        atoms = _records("chains.pdb", "HETATM")
        assert [(line[21], int(line[22:26])) for line in atoms] == [
            ("A", 1),
            ("B", 1),
        ]

    def test_repeated_explicit_numbers_warn(self, caplog):
        chain = mb.Compound(
            [_residue("LIG", 5, ["C1"], 0.0), _residue("LIG", 5, ["C1"], 0.5)]
        )
        chain.save("repeated.pdb")
        assert "share residue numbers [5]" in caplog.text

    def test_chains_with_the_same_id_stay_separate(self):
        chains = []
        for x in (0.0, 0.5):
            chain = mb.Compound([_residue("LIG", 1, ["C1"], x)])
            chain.chain_id = "A"
            chains.append(chain)
        mb.Compound(chains).save("same_id.pdb")
        assert len(_records("same_id.pdb", "TER")) == 2

    def test_repeated_atom_names_are_numbered(self, ethane, caplog):
        ethane.save("ethane.pdb")
        names = [line[12:16].strip() for line in _records("ethane.pdb", "HETATM")]
        assert names == ["C1", "H1", "H2", "H3", "C2", "H4", "H5", "H6"]
        assert {p.name for p in ethane.particles()} == {"C", "H"}
        assert "numbered in the file" in caplog.text

    def test_long_residue_name_is_cut_with_warning(self, caplog):
        residue = _residue("LIGAND", 1, ["C1"], 0.0)
        residue.save("long.pdb")
        assert _records("long.pdb", "HETATM")[0][17:21] == "LIGA"
        assert "longer than three characters" in caplog.text

    def test_long_atom_names_are_cut_with_warning(self, caplog):
        residue = _residue("LIG", 1, ["C1234", "O5678"], 0.0)
        residue.save("long.pdb")
        names = [line[12:16].strip() for line in _records("long.pdb", "HETATM")]
        assert names == ["C123", "O567"]
        assert "cut to four" in caplog.text

    def test_too_many_repeated_names_raise(self):
        residue = _residue("RES", 1, ["CA"] * 99, 0.0)
        names = set(pdb_residues(residue)[0].atom_names.values())
        assert len(names) == 99 and "CA99" in names
        residue.add(mb.Particle(name="CA", element="C"))
        with pytest.raises(MBuildError, match="residue_number"):
            residue.save("too_many.pdb")

    def test_values_too_wide_raise(self):
        far = mb.Particle(name="C", element="C", pos=[1000.0, 0, 0])
        with pytest.raises(MBuildError, match="do not fit"):
            far.save("far.pdb")
        charged = mb.Particle(name="C", element="C", formal_charge=10)
        with pytest.raises(MBuildError, match="charge columns"):
            charged.save("charged.pdb")

    def test_no_cryst1_without_box(self, ethane):
        ethane.save("ethane.pdb")
        assert _records("ethane.pdb", "CRYST1") == []

    def test_large_residue_numbers_warn(self, caplog):
        _residue("LIG", 10000, ["C1"], 0.0).save("big.pdb")
        assert _records("big.pdb", "HETATM")[0][22:26] == "A000"
        assert "hexadecimal" in caplog.text

    def test_atom_name_starting_with_digit(self):
        residue = _residue("LIG", 1, ["C"], 0.0)
        residue.add(mb.Particle(name="1HB", element="H", pos=[0.1, 0, 0]))
        residue.save("digit.pdb")
        names = [line[12:16] for line in _records("digit.pdb", "HETATM")]
        assert names == [" C  ", "1HB "]

    def test_large_serials(self):
        assert _encode_index(99999, 5) == "99999"
        assert _encode_index(100000, 5) == "A0000"
        with pytest.raises(MBuildError):
            _encode_index(10**5 + 6 * 16**4, 5)
