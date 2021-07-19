import numpy as np
import pytest

import mbuild as mb
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestLammpsData(BaseTest):
    def test_save(self, ethane):
        ethane.save(filename="ethane.lammps")

    @pytest.mark.parametrize("unit_style", ["real", "lj"])
    def test_save_forcefield(self, ethane, unit_style):
        ethane.save(
            filename="ethane-opls.lammps",
            forcefield_name="oplsaa",
            unit_style=unit_style,
        )

    @pytest.mark.parametrize("pair_coeff_label", [None, "lj", "lj/coul/cut"])
    def test_save_pair_coeff_label(self, ethane, pair_coeff_label):
        ethane.save(
            filename="ethane-opls.lammps",
            forcefield_name="oplsaa",
            unit_style="real",
            pair_coeff_label=pair_coeff_label,
        )
        with open("ethane-opls.lammps", "r") as fp:
            for line in fp:
                if "Pair Coeffs" in line:
                    if pair_coeff_label is None:
                        assert "# lj" in line
                    else:
                        assert "# {}".format(pair_coeff_label) in line

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_charmm(self):
        from foyer import Forcefield

        from mbuild.formats.lammpsdata import write_lammpsdata

        cmpd = mb.load(get_fn("charmm_dihedral.mol2"))
        for i in cmpd.particles():
            i.name = "_{}".format(i.name)
        structure = cmpd.to_parmed(
            box=cmpd.get_boundingbox(),
            residues=set([p.parent.name for p in cmpd.particles()]),
        )
        ff = Forcefield(forcefield_files=[get_fn("charmm_truncated.xml")])
        structure = ff.apply(structure, assert_dihedral_params=False)
        write_lammpsdata(structure, "charmm_dihedral.lammps")
        out_lammps = open("charmm_dihedral.lammps", "r").readlines()
        found_angles = False
        found_dihedrals = False
        for i, line in enumerate(out_lammps):
            if "Angle Coeffs" in line:
                assert "# charmm" in line
                assert (
                    "#\tk(kcal/mol/rad^2)\t\ttheteq(deg)\tk(kcal/mol/angstrom^2)\treq(angstrom)\n"
                    in out_lammps[i + 1]
                )
                assert len(out_lammps[i + 2].split("#")[0].split()) == 5
                found_angles = True
            elif "Dihedral Coeffs" in line:
                assert "# charmm" in line
                assert "#k, n, phi, weight" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 5
                found_dihedrals = True
            else:
                pass
        assert found_angles
        assert found_dihedrals

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_singleterm_charmm(self):
        from foyer import Forcefield

        from mbuild.formats.lammpsdata import write_lammpsdata

        cmpd = mb.load(get_fn("charmm_dihedral.mol2"))
        for i in cmpd.particles():
            i.name = "_{}".format(i.name)
        structure = cmpd.to_parmed(
            box=cmpd.get_boundingbox(),
            residues=set([p.parent.name for p in cmpd.particles()]),
        )
        ff = Forcefield(
            forcefield_files=[get_fn("charmm_truncated_singleterm.xml")]
        )
        structure = ff.apply(structure, assert_dihedral_params=False)
        write_lammpsdata(structure, "charmm_dihedral_singleterm.lammps")
        out_lammps = open("charmm_dihedral_singleterm.lammps", "r").readlines()
        found_dihedrals = False
        for i, line in enumerate(out_lammps):
            if "Dihedral Coeffs" in line:
                assert "# charmm" in line
                assert "#k, n, phi, weight" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 5
                assert float(
                    out_lammps[i + 2].split("#")[0].split()[4]
                ) == float("1.0")
                found_dihedrals = True
            else:
                pass
        assert found_dihedrals

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_charmm_improper(self):
        from foyer import Forcefield

        import mbuild as mb
        from mbuild.formats.lammpsdata import write_lammpsdata

        system = mb.Compound()
        first = mb.Particle(name="_CTL2", pos=[-1, 0, 0])
        second = mb.Particle(name="_CL", pos=[0, 0, 0])
        third = mb.Particle(name="_OBL", pos=[1, 0, 0])
        fourth = mb.Particle(name="_OHL", pos=[0, 1, 0])

        system.add([first, second, third, fourth])

        system.add_bond((first, second))
        system.add_bond((second, third))
        system.add_bond((second, fourth))

        ff = Forcefield(forcefield_files=[get_fn("charmm36_cooh.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )

        write_lammpsdata(struc, "charmm_improper.lammps")
        out_lammps = open("charmm_improper.lammps", "r").readlines()
        found_impropers = False
        for i, line in enumerate(out_lammps):
            if "Improper Coeffs" in line:
                assert "# harmonic" in line
                assert "k, phi" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 3
                found_impropers = True
        assert found_impropers

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_amber(self):
        from foyer import Forcefield

        from mbuild.formats.lammpsdata import write_lammpsdata

        cmpd = mb.load("C1(=CC=CC=C1)F", smiles=True)

        ff = Forcefield(forcefield_files=[get_fn("gaff_test.xml")])
        structure = ff.apply(cmpd)

        write_lammpsdata(
            structure, "amber.lammps", zero_dihedral_weighting_factor=True
        )
        out_lammps = open("amber.lammps", "r").readlines()
        found_angles = False
        found_dihedrals = False
        found_impropers = False
        for i, line in enumerate(out_lammps):
            if "Angle Coeffs" in line:
                assert "# harmonic" in line
                assert "#\treduced_k\t\ttheteq(deg)" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 3
                found_angles = True
            elif "Dihedral Coeffs" in line:
                assert "# charmm" in line
                assert "#k, n, phi, weight" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 5
                assert float(out_lammps[i + 2].split("#")[0].split()[4]) == 0.0
                found_dihedrals = True
            elif "Improper Coeffs" in line:
                assert "# cvff" in line
                assert "#K, d, n" in out_lammps[i + 1]
                assert len(out_lammps[i + 2].split("#")[0].split()) == 4
                assert out_lammps[i + 2].split("#")[0].split()[2] == "-1"
                found_impropers = True
            else:
                pass
        assert found_angles
        assert found_dihedrals
        assert found_impropers

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    @pytest.mark.parametrize("unit_style", ["real", "lj"])
    def test_save_box(self, ethane, unit_style):
        box = mb.Box(
            lengths=np.array([2.0, 2.0, 2.0]), angles=[90.0, 90.0, 90.0]
        )
        ethane.save(
            filename="ethane-box.lammps",
            forcefield_name="oplsaa",
            box=box,
            unit_style=unit_style,
        )

    def test_nbfix(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        # Add nbfixes
        types = list(set([a.atom_type for a in structure.atoms]))
        types[0].add_nbfix(types[1].name, 1.2, 2.1)
        types[1].add_nbfix(types[0].name, 1.2, 2.1)
        write_lammpsdata(filename="nbfix.lammps", structure=structure)

        checked_section = False
        with open("nbfix.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "PairIJ Coeffs" in line:
                    fi.readline()
                    line = fi.readline().partition("#")[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float),
                        [1, 1, 0.066, 3.5],
                    )
                    line = fi.readline().partition("#")[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float),
                        [1, 2, 2.1, 1.06907846],
                    )
                    line = fi.readline().partition("#")[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float), [2, 2, 0.03, 2.5]
                    )
                    line = fi.readline()
                    checked_section = True
                # Break if PairIJ Coeffs is not found
                if "Atoms" in line:
                    break

    def test_save_triclinic_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]), angles=[60, 70, 80])
        ethane.save(
            filename="triclinic-box.lammps", forcefield_name="oplsaa", box=box
        )

    def test_save_orthorhombic_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(
            filename="ortho-box.lammps", forcefield_name="oplsaa", box=box
        )
        with open("ortho-box.lammps", "r") as fi:
            for line in fi:
                assert "xy" not in line
                assert "xz" not in line
                assert "yz" not in line

    @pytest.mark.parametrize(
        "atom_style, n_columns",
        [("full", 7), ("atomic", 5), ("molecular", 6), ("charge", 6)],
    )
    def test_writing_atom_styles(self, ethane, atom_style, n_columns):
        ethane.save(filename="ethane.lammps", atom_style=atom_style)
        with open("ethane.lammps", "r") as f:
            for line in f:
                if "Atoms" not in line:
                    continue
                    atoms_header = next(f)
                    first_atom_line = next(f)
                    columns = first_atom_line.split("\t")
                    assert len(columns) == n_columns
                else:
                    assert "# " + atom_style in line

    def test_resid(self, ethane, methane):
        structure = ethane.to_parmed() + methane.to_parmed()
        n_atoms = len(structure.atoms)
        write_lammpsdata(structure, "compound.lammps")
        res_list = list()
        with open("compound.lammps", "r") as f:
            for i, line in enumerate(f):
                if "Atoms" in line:
                    break
        atom_lines = open("compound.lammps", "r").readlines()[
            i + 2 : i + n_atoms + 2
        ]
        for line in atom_lines:
            res_list.append(line.rstrip().split()[1])

        assert set(res_list) == set(["1", "0"])

    def test_box_bounds(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        box = mb.Box.from_mins_maxs_angles(
            mins=np.array([-1.0, -2.0, -3.0]),
            maxs=np.array([3.0, 2.0, 1.0]),
            angles=[90.0, 90.0, 90.0],
        )
        # box = ethane.get_boundingbox()

        write_lammpsdata(
            filename="box.lammps",
            structure=structure,
            unit_style="real",
            mins=[0.0, 0.0, 0.0],
            maxs=[m for m in box.lengths],
        )

        checked_section = False
        # NOTE, need to figure out how to handle lo and hi of a compound
        with open("box.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "xlo" in line:
                    xlo = float(line.split()[0])
                    xhi = float(line.split()[1])
                    assert np.isclose(xlo, 0.0)
                    assert np.isclose(xhi, 40.0)

                    line = fi.readline()
                    ylo = float(line.split()[0])
                    yhi = float(line.split()[1])
                    assert np.isclose(ylo, 0.0)
                    assert np.isclose(yhi, 40.0)

                    line = fi.readline()
                    zlo = float(line.split()[0])
                    zhi = float(line.split()[1])
                    assert np.isclose(zlo, 0.0)
                    assert np.isclose(zhi, 40.0)

                    checked_section = True

        write_lammpsdata(
            filename="box.lammps", structure=structure, unit_style="real"
        )

        checked_section = False
        with open("box.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "xlo" in line:
                    xlo = float(line.split()[0])
                    xhi = float(line.split()[1])
                    assert np.isclose(xlo, 0.0)
                    assert np.isclose(xhi, 7.13999987)

                    line = fi.readline()
                    ylo = float(line.split()[0])
                    yhi = float(line.split()[1])
                    assert np.isclose(ylo, 0.0)
                    assert np.isclose(yhi, 7.93800011)

                    line = fi.readline()
                    zlo = float(line.split()[0])
                    zhi = float(line.split()[1])
                    assert np.isclose(zlo, 0.0)
                    assert np.isclose(zhi, 6.646)

                    checked_section = True

    def test_lj_box(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "dihedral types" in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(float(line), 2.04)
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 2.268)
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 1.898857)
                    checked_section = True

    def test_lj_masses(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Masses" in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 1.00)
                    checked_section = True

    def test_lj_pairs(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Pair Coeffs" in line:
                    fi.readline()
                    fi.readline()
                    line = fi.readline().split()
                    epsilon = float(line[1])
                    sigma = float(line[2])
                    assert np.isclose(epsilon, 1.00)
                    assert np.isclose(sigma, 1.00)
                    checked_section = True

    def test_lj_bonds(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Bond Coeffs" in line:
                    fi.readline()
                    bonds = list()
                    bonds.append(float(fi.readline().split()[1]))
                    bonds.append(float(fi.readline().split()[1]))
                    assert np.allclose(sorted(bonds), [49742.424, 63106.06])
                    checked_section = True

    def test_lj_angles(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Angle Coeffs" in line:
                    fi.readline()
                    angles = list()
                    angles.append(float(fi.readline().split()[1]))
                    angles.append(float(fi.readline().split()[1]))

                    assert np.allclose(sorted(angles), [6125.0, 6960.227])
                    checked_section = True

    def test_lj_dihedrals(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name="oplsaa")
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(
            filename="lj.lammps", structure=structure, unit_style="lj"
        )

        checked_section = False
        with open("lj.lammps", "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Dihedral Coeffs" in line:
                    fi.readline()
                    dihedrals = fi.readline().split()[1:5]
                    dihedrals = [float(i) for i in dihedrals]
                    assert np.allclose(dihedrals, [0.0005, 0.0, 4.5455, -0.0])
                    checked_section = True
