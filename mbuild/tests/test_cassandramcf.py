import pytest
from numpy import isclose

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestCassandraMCF(BaseTest):
    def test_not_parameterized(self, ethane):
        with pytest.raises(
            ValueError, match=r"MCF writing not supported without"
        ):
            ethane.save(
                filename="ethane.mcf",
                angle_style="harmonic",
                dihedral_style="opls",
            )

    def test_invalid_structure(self, ethane):
        with pytest.raises(ValueError, match=r"requires parmed structure"):
            mb.formats.cassandramcf.write_mcf(
                ethane,
                "ethane.mcf",
                angle_style="harmonic",
                dihedral_style="opls",
            )

    def test_invalid_angle_style(self, ethane):
        with pytest.raises(
            ValueError, match=r"Invalid selection for angle_style"
        ):
            ethane.save(
                filename="ethane-opls.mcf",
                forcefield_name="oplsaa",
                angle_style="harm",
                dihedral_style="opls",
            )

    def test_invalid_dihedral_style(self, ethane):
        with pytest.raises(
            ValueError, match=r"Invalid selection for dihedral_style"
        ):
            ethane.save(
                filename="ethane-opls.mcf",
                forcefield_name="oplsaa",
                angle_style="harmonic",
                dihedral_style="op",
            )

    def test_dihedral_style_none(self, ethane):
        ethane.save(
            filename="ethane-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="harmonic",
            dihedral_style="none",
        )

        mcf_data = []
        with open("ethane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Dihedral_Info":
                    dihedral_section_start = idx

        assert mcf_data[dihedral_section_start + 2][5] == "none"

    def test_no_dihedrals(self, methane):
        methane.save(
            filename="methane-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="harmonic",
            dihedral_style="none",
        )

        mcf_data = []
        with open("methane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Dihedral_Info":
                    dihedral_section_start = idx

        assert mcf_data[dihedral_section_start + 1][0] == "0"

    def test_no_14(self, methane):
        methane.save(
            filename="methane-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="harmonic",
            dihedral_style="none",
        )

        mcf_data = []
        with open("methane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Intra_Scaling":
                    intrascaling_section_start = idx

        assert isclose(float(mcf_data[intrascaling_section_start + 1][2]), 0.0)
        assert isclose(float(mcf_data[intrascaling_section_start + 2][2]), 0.0)

    def test_infer_14(self, ethane):
        ethane.save(
            filename="ethane-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="harmonic",
            dihedral_style="none",
        )

        mcf_data = []
        with open("ethane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Intra_Scaling":
                    intrascaling_section_start = idx

        assert isclose(float(mcf_data[intrascaling_section_start + 1][2]), 0.5)
        assert isclose(float(mcf_data[intrascaling_section_start + 2][2]), 0.5)

    def test_unmatched_dihedral_style(self, ethane):
        with pytest.raises(ValueError, match=r"but RB torsions found"):
            ethane.save(
                filename="ethane-opls.mcf",
                forcefield_name="oplsaa",
                angle_style="harmonic",
                dihedral_style="charmm",
            )

    def test_unreasonable_lj14(self, ethane):
        import foyer

        oplsaa = foyer.Forcefield(name="oplsaa")
        ethane = oplsaa.apply(ethane)
        with pytest.raises(ValueError, match=r"Unreasonable value"):
            mb.formats.cassandramcf.write_mcf(
                ethane,
                "ethane.mcf",
                angle_style="harmonic",
                dihedral_style="opls",
                lj14=2.0,
            )

    def test_unreasonable_coul14(self, ethane):
        import foyer

        oplsaa = foyer.Forcefield(name="oplsaa")
        ethane = oplsaa.apply(ethane)
        with pytest.raises(ValueError, match=r"Unreasonable value"):
            mb.formats.cassandramcf.write_mcf(
                ethane,
                "ethane.mcf",
                angle_style="harmonic",
                dihedral_style="opls",
                coul14=-1.0,
            )

    def test_multiple_molecules(self, ethane):
        n_ethane = 2
        ethane.name = "Ethane"
        filled = mb.fill_box(
            ethane, n_compounds=n_ethane, box=[0, 0, 0, 4, 4, 4]
        )
        with pytest.raises(
            ValueError, match=r"Not all components of the molecule"
        ):
            filled.save(
                filename="box-ethane-opls.mcf",
                forcefield_name="oplsaa",
                angle_style="harmonic",
                dihedral_style="opls",
            )

    def test_save_forcefield(self, ethane):
        ethane.save(
            filename="ethane-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="harmonic",
            dihedral_style="opls",
        )

        mcf_data = []
        with open("ethane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Atom_Info":
                    atom_section_start = idx
                elif line[1] == "Bond_Info":
                    bond_section_start = idx
                elif line[1] == "Angle_Info":
                    angle_section_start = idx
                elif line[1] == "Dihedral_Info":
                    dihedral_section_start = idx
                elif line[1] == "Improper_Info":
                    improper_section_start = idx
                elif line[1] == "Fragment_Info":
                    fragment_section_start = idx
                elif line[1] == "Fragment_Connectivity":
                    fragment_conn_start = idx

        # Check a some atom info
        assert mcf_data[atom_section_start + 1][0] == "8"
        assert mcf_data[atom_section_start + 2][1] == "opls_135"
        assert isclose(float(mcf_data[atom_section_start + 2][3]), 12.011)
        assert isclose(float(mcf_data[atom_section_start + 2][4]), -0.18)
        assert mcf_data[atom_section_start + 2][5] == "LJ"
        assert isclose(float(mcf_data[atom_section_start + 2][6]), 33.212)
        assert isclose(float(mcf_data[atom_section_start + 2][7]), 3.500)

        # Bond info
        assert mcf_data[bond_section_start + 1][0] == "7"
        passed_test = False
        for line in mcf_data[bond_section_start + 2 : angle_section_start - 6]:
            a1 = line[1]
            a2 = line[2]
            if (a1 == "1" and a2 == "2") or (a2 == "1" and a1 == "2"):
                assert line[3] == "fixed"
                assert isclose(float(line[4]), 1.090)
                passed_test = True
                break
        assert passed_test

        # Angle info
        assert mcf_data[angle_section_start + 1][0] == "12"
        passed_test = False
        for line in mcf_data[
            angle_section_start + 2 : dihedral_section_start - 8
        ]:
            if line[2] == "5":
                a1 = line[1]
                a3 = line[3]
                if (a1 == "1" and a3 == "6") or (a3 == "1" and a1 == "6"):
                    assert line[4] == "harmonic"
                    assert isclose(float(line[5]), 18870.7)
                    assert isclose(float(line[6]), 110.7)
                    passed_test = True
                    break
        assert passed_test

        # Dihedral info
        assert mcf_data[dihedral_section_start + 1][0] == "9"
        passed_test = False
        for line in mcf_data[
            dihedral_section_start + 2 : improper_section_start - 5
        ]:
            a1 = line[1]
            a2 = line[2]
            a3 = line[3]
            a4 = line[4]
            if (a1 == "2" and a2 == "1" and a3 == "5" and a4 == "6") or (
                a4 == "2" and a3 == "1" and a2 == "5" and a1 == "6"
            ):
                assert line[5] == "OPLS"
                assert isclose(float(line[6]), 0.000)
                assert isclose(float(line[7]), 0.000)
                assert isclose(float(line[8]), 0.000)
                assert isclose(float(line[9]), 0.628)

        assert mcf_data[improper_section_start + 1][0] == "0"
        assert mcf_data[fragment_section_start + 1][0] == "2"

        # Check fragment connectivity
        assert mcf_data[fragment_conn_start + 1][0] == "1"
        assert mcf_data[fragment_conn_start + 2][0] == "1"
        assert mcf_data[fragment_conn_start + 2][1] == "1"
        assert mcf_data[fragment_conn_start + 2][2] == "2"

    def test_save_ring_forcefield(self, benzene):
        benzene.save(
            filename="benzene-opls.mcf",
            forcefield_name="oplsaa",
            angle_style="fixed",
            dihedral_style="opls",
        )

        mcf_data = []
        with open("benzene-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Atom_Info":
                    atom_section_start = idx
                elif line[1] == "Bond_Info":
                    bond_section_start = idx
                elif line[1] == "Angle_Info":
                    angle_section_start = idx
                elif line[1] == "Dihedral_Info":
                    dihedral_section_start = idx
                elif line[1] == "Improper_Info":
                    improper_section_start = idx
                elif line[1] == "Fragment_Info":
                    fragment_section_start = idx
                elif line[1] == "Fragment_Connectivity":
                    fragment_conn_start = idx

        # Check a some atom info
        assert mcf_data[atom_section_start + 1][0] == "12"
        assert mcf_data[atom_section_start + 2][1] == "opls_145"
        assert isclose(float(mcf_data[atom_section_start + 2][3]), 12.011)
        assert isclose(float(mcf_data[atom_section_start + 2][4]), -0.115)
        assert mcf_data[atom_section_start + 2][5] == "LJ"
        assert isclose(float(mcf_data[atom_section_start + 2][6]), 35.225)
        assert isclose(float(mcf_data[atom_section_start + 2][7]), 3.550)
        assert mcf_data[atom_section_start + 2][8] == "ring"

        # Bond info
        assert mcf_data[bond_section_start + 1][0] == "12"
        passed_test = False
        for line in mcf_data[bond_section_start + 2 : angle_section_start - 6]:
            a1 = line[1]
            a2 = line[2]
            if (a1 == "1" and a2 == "2") or (a2 == "1" and a1 == "2"):
                assert line[3] == "fixed"
                assert isclose(float(line[4]), 1.400)
                passed_test = True
                break
        assert passed_test

        # Angle info
        assert mcf_data[angle_section_start + 1][0] == "18"
        passed_test = False
        for line in mcf_data[
            angle_section_start + 2 : dihedral_section_start - 8
        ]:
            if line[2] == "2":
                a1 = line[1]
                a3 = line[3]
                if (a1 == "1" and a3 == "3") or (a3 == "1" and a1 == "3"):
                    assert line[4] == "fixed"
                    assert isclose(float(line[5]), 120.00)
                    passed_test = True
                    break
        assert passed_test

        # Dihedral info
        assert mcf_data[dihedral_section_start + 1][0] == "24"
        passed_test = False
        for line in mcf_data[
            dihedral_section_start + 2 : improper_section_start - 5
        ]:
            a1 = line[1]
            a2 = line[2]
            a3 = line[3]
            a4 = line[4]
            if (a1 == "1" and a2 == "2" and a3 == "3" and a4 == "4") or (
                a4 == "1" and a3 == "2" and a2 == "3" and a1 == "4"
            ):
                assert line[5] == "OPLS"
                assert isclose(float(line[6]), 0.000)
                assert isclose(float(line[7]), -0.000)
                assert isclose(float(line[8]), 15.167)
                assert isclose(float(line[9]), -0.000)

        assert mcf_data[improper_section_start + 1][0] == "0"
        assert mcf_data[fragment_section_start + 1][0] == "1"

        # Check fragment connectivity
        assert mcf_data[fragment_conn_start + 1][0] == "0"

    def test_shorten_atomname(self, ethane):
        import foyer

        from mbuild.formats.cassandramcf import write_mcf

        typed_ethane = foyer.forcefields.load_OPLSAA().apply(ethane)
        typed_ethane[0].type = "C_very_very_very_extended"
        write_mcf(
            typed_ethane,
            "ethane-opls.mcf",
            angle_style="harmonic",
            dihedral_style="opls",
        )

        mcf_data = []
        with open("ethane-opls.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Atom_Info":
                    atom_section_start = idx

        assert mcf_data[atom_section_start + 2][1] == "y_very_very_extended"

    def test_fused_rings(self):
        import foyer

        import mbuild
        from mbuild.formats.cassandramcf import write_mcf

        naph = mbuild.load("C1=CC=C2C=CC=CC2=C1", smiles=True)
        # Note the atomtyping is wrong -- doesn't matter for test though
        naph_ff = foyer.forcefields.load_OPLSAA().apply(naph)
        write_mcf(
            naph_ff, "naph.mcf", angle_style="harmonic", dihedral_style="opls"
        )

        mcf_data = []
        with open("naph.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Fragment_Info":
                    frag_section_start = idx

        assert int(mcf_data[frag_section_start + 1][0]) == 1
        assert int(mcf_data[frag_section_start + 2][1]) == 18

        tmada = mbuild.load("C[C](C)(C)C12CC3CC(C1)CC(C3)C2", smiles=True)
        tmada_ff = foyer.forcefields.load_OPLSAA().apply(tmada)
        write_mcf(
            tmada_ff, "tmada.mcf", angle_style="harmonic", dihedral_style="opls"
        )

        mcf_data = []
        with open("tmada.mcf") as f:
            for line in f:
                mcf_data.append(line.strip().split())

        for idx, line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == "Fragment_Info":
                    frag_section_start = idx

        assert int(mcf_data[frag_section_start + 1][0]) == 5
        assert int(mcf_data[frag_section_start + 2][1]) == 26
        assert int(mcf_data[frag_section_start + 3][1]) == 5
        assert int(mcf_data[frag_section_start + 4][1]) == 5
        assert int(mcf_data[frag_section_start + 5][1]) == 5
        assert int(mcf_data[frag_section_start + 6][1]) == 5

    def test_infer_14_scaling_zero_eps(self):
        import foyer

        import mbuild
        from mbuild.formats.cassandramcf import write_mcf

        mol = mbuild.load("CO", smiles=True)
        mol_ff = foyer.forcefields.load_OPLSAA().apply(mol)
        with pytest.warns(UserWarning):
            write_mcf(
                mol_ff,
                "test.mcf",
                angle_style="harmonic",
                dihedral_style="opls",
            )
