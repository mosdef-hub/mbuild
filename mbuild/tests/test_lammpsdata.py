from pathlib import Path

import numpy as np
import pytest
from pytest import FixtureRequest
from pathlib import Path

import mbuild as mb
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer

KCAL_TO_KJ = 4.184
ANG_TO_NM = 0.10


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestLammpsData(BaseTest):
    def test_save(self, ethane):
        ethane.save(filename="ethane.lammps")

    @pytest.fixture(scope="session")
    def lj_save(self, tmpdir_factory, sigma=None, epsilon=None, mass=None):
        tmp = tmpdir_factory

        def _create_lammps(ethane, tmp, sigma=sigma, epsilon=epsilon, mass=mass):
            from foyer import Forcefield

            OPLSAA = Forcefield(name="oplsaa")
            structure = OPLSAA.apply(ethane)
            fn = tmpdir_factory.mktemp("data").join("lj.lammps")
            write_lammpsdata(
                filename=str(fn),
                structure=structure,
                unit_style="lj",
                sigma_conversion_factor=sigma,
                epsilon_conversion_factor=epsilon,
                mass_conversion_factor=mass,
            )
            return str(fn)

        return _create_lammps

    @pytest.fixture(scope="session")
    def real_save(self, tmpdir_factory):
        tmp = tmpdir_factory

        def _create_lammps(ethane, tmp):
            from foyer import Forcefield

            OPLSAA = Forcefield(name="oplsaa")
            structure = OPLSAA.apply(ethane)
            fn = tmpdir_factory.mktemp("data").join("lj.lammps")
            write_lammpsdata(filename=str(fn), structure=structure, unit_style="real")
            return str(fn)

        return _create_lammps

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
                assert (
                    "#\tk(kcal/mol/rad^2)\t\ttheteq(deg)" in out_lammps[i + 1]
                )
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
    def test_save_forcefield_with_same_struct(self):
        from foyer import Forcefield

        from mbuild.formats.lammpsdata import write_lammpsdata

        system = mb.load("C1(=CC=CC=C1)F", smiles=True)

        ff = Forcefield(forcefield_files=[get_fn("gaff_test.xml")])
        struc = ff.apply(
            system,
            assert_angle_params=False,
            assert_dihedral_params=False,
            assert_improper_params=False,
        )
        write_lammpsdata(
            struc, "charmm_improper.lammps", zero_dihedral_weighting_factor=True
        )
        for i in range(3):
            xyz = struc.coordinates
            xyz = xyz + np.array([1, 1, 1])
            struc.coordinates = xyz
            write_lammpsdata(
                struc,
                f"charmm_improper{i}.lammps",
                zero_dihedral_weighting_factor=True,
            )

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

    @pytest.mark.parametrize(
        "offset, expected_value", [(0, ("0", "1")), (1, ("1", "2"))]
    )
    def test_resid(self, offset, expected_value, ethane, methane):
        structure = ethane.to_parmed() + methane.to_parmed()
        n_atoms = len(structure.atoms)
        write_lammpsdata(structure, "compound.lammps", moleculeID_offset=offset)
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
        assert set(res_list) == set(expected_value)

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

    def test_lj_box(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        ethane_sigma = 0.35  # nm
        box_length = np.array(ethane.get_boundingbox().lengths) + [
            0.5,
            0.5,
            0.5,
        ]  # 0.5nm buffer around box
        box_length /= ethane_sigma  # Unitless
        with open(str(fn), "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "dihedral types" in line:  # line before box info
                    fi.readline()
                    for i in range(len(box_length)):
                        line = float(fi.readline().split()[1])
                        assert np.isclose(float(line), box_length[i])
                    checked_section = True

    def test_lj_masses(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Masses" in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 1.00)
                    checked_section = True

    def test_real_masses(self, ethane, real_save):
        fn = real_save(ethane, Path.cwd())
        checked_section = False
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Masses" in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 12.010780)
                    checked_section = True

    def test_lj_pairs(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Pair Coeffs" in line:
                    fi.readline()
                    line = fi.readline().split()
                    epsilon = float(line[1])
                    sigma = float(line[2])
                    assert np.isclose(epsilon, 1.00)
                    assert np.isclose(sigma, 1.00)
                    checked_section = True

    def test_lj_passed_pairs(self, ethane, lj_save):
        fn = lj_save(
            ethane,
            Path.cwd(),
            sigma=1,
            epsilon=1,
        )
        checked_section = False
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Pair Coeffs" in line:
                    fi.readline()
                    line = fi.readline().split()
                    epsilon = float(line[1])
                    sigma = float(line[2])
                    assert np.isclose(epsilon, 0.276144, atol=1e-5)
                    assert np.isclose(sigma, 0.35)
                    checked_section = True

    def test_real_pairs(self, ethane, real_save):
        fn = real_save(ethane, Path.cwd())
        checked_section = False
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Pair Coeffs" in line:
                    fi.readline()
                    line = fi.readline().split()
                    epsilon = float(line[1])
                    sigma = float(line[2])
                    assert np.isclose(epsilon, 0.066)
                    assert np.isclose(sigma, 3.5)
                    checked_section = True

    def test_lj_bonds(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        ethane_bondsk = (
            np.array([224262.4, 284512.0]) * (0.35) ** 2 / (0.276144) / 2
        )  # kj/mol/nm**2 to lammps
        ethane_bondsreq = np.array([0.1529, 0.109]) * 0.35 ** -1  # unitless
        ethane_lj_bonds = zip(ethane_bondsk, ethane_bondsreq)
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Bond Coeffs" in line:
                    fi.readline()
                    for bond_params in ethane_lj_bonds:
                        print(bond_params)
                        line = fi.readline()
                        assert np.allclose(
                            (float(line.split()[1]), float(line.split()[2])),
                            bond_params,
                            atol=1e-3,
                        )
                    checked_section = True

    def test_real_bonds(self, ethane, real_save):
        fn = real_save(ethane, Path.cwd())
        checked_section = False
        ethane_bondsk = (
            np.array([224262.4, 284512.0]) / KCAL_TO_KJ * ANG_TO_NM ** 2 / 2
        )  # kj/mol/nm**2 to lammps
        ethane_bondsreq = np.array([0.1529, 0.109]) / ANG_TO_NM
        ethane_real_bonds = zip(ethane_bondsk, ethane_bondsreq)
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Bond Coeffs" in line:
                    fi.readline()
                    for bond_params in ethane_real_bonds:
                        print(bond_params)
                        line = fi.readline()
                        assert np.allclose(
                            (float(line.split()[1]), float(line.split()[2])),
                            bond_params,
                            atol=1e-3,
                        )
                    checked_section = True

    def test_lj_angles(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        ethane_anglesk = (
            np.array([313.8, 276.144]) / 0.276144 / 2
        )  # kj/mol/nm**2 to lammps
        ethane_anglesreq = (
            np.array([1.93207948196, 1.88146493365]) * 180 / np.pi
        )  # radians to lammps
        ethane_lj_angles = zip(ethane_anglesk, ethane_anglesreq)
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Angle Coeffs" in line:
                    fi.readline()
                    for angle_params in ethane_lj_angles:
                        line = fi.readline()
                        assert np.allclose(
                            (float(line.split()[1]), float(line.split()[2])),
                            angle_params,
                        )
                    checked_section = True

    def test_real_angles(self, ethane, real_save):
        fn = real_save(ethane, Path.cwd())
        checked_section = False
        ethane_anglesk = (
            np.array([313.8, 276.144]) / 2 / KCAL_TO_KJ
        )  # kj/mol to lammps kcal/mol
        ethane_anglesreq = (
            np.array([1.93207948196, 1.88146493365]) * 180 / np.pi
        )  # radians to lammps
        ethane_real_angles = zip(ethane_anglesk, ethane_anglesreq)
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Angle Coeffs" in line:
                    fi.readline()
                    for angle_params in ethane_real_angles:
                        line = fi.readline()
                        assert np.allclose(
                            (float(line.split()[1]), float(line.split()[2])),
                            angle_params,
                        )
                    checked_section = True

    def test_lj_dihedrals(self, ethane, lj_save):
        fn = lj_save(ethane, Path.cwd())
        checked_section = False
        ethane_diheds = (
            np.array([-0.00000, -0.00000, 0.30000, -0.00000]) / 0.066
        )  # kcal/mol to lammps
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Dihedral Coeffs" in line:
                    fi.readline()
                    line = fi.readline()
                    assert np.allclose(
                        (
                            float(line.split()[1]),
                            float(line.split()[2]),
                            float(line.split()[3]),
                            float(line.split()[4]),
                        ),
                        ethane_diheds,
                    )
                    checked_section = True

    def test_real_dihedrals(self, ethane, real_save):
        fn = real_save(ethane, Path.cwd())
        checked_section = False
        ethane_diheds = np.array([-0.00000, -0.00000, 0.30000, -0.00000])
        with open(fn, "r") as fi:
            while not checked_section:
                line = fi.readline()
                if "Dihedral Coeffs" in line:
                    fi.readline()
                    line = fi.readline()
                    assert np.allclose(
                        (
                            float(line.split()[1]),
                            float(line.split()[2]),
                            float(line.split()[3]),
                            float(line.split()[4]),
                        ),
                        ethane_diheds,
                    )
                    checked_section = True
