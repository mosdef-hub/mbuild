from collections import OrderedDict

import numpy as np
import pytest
from foyer.forcefields import forcefields

import mbuild as mb
from mbuild import Box, Compound
from mbuild.formats import charmm_writer
from mbuild.formats.charmm_writer import Charmm
from mbuild.lattice import load_cif
from mbuild.tests.base_test import BaseTest
from mbuild.utils.conversion import (
    base10_to_base16_alph_num,
    base10_to_base26_alph,
    base10_to_base52_alph,
    base10_to_base62_alph_num,
)
from mbuild.utils.io import get_fn, has_foyer
from mbuild.utils.specific_ff_to_residue import specific_ff_to_residue


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestCharmmWriterData(BaseTest):
    def test_save(self, ethane_gomc):
        Charmm(
            ethane_gomc,
            "ethane",
            ff_filename="ethane",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

    def test_save_charmm_gomc_ff(self, ethane_gomc):
        charmm = Charmm(
            ethane_gomc,
            "charmm_data",
            ff_filename="charmm_data",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()

        with open("charmm_data.inp", "r") as fp:
            masses_read = False
            bonds_read = False
            angles_read = False
            dihedrals_read = False
            nonbondeds_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):

                if (
                    "!atom_types" in line
                    and "mass" in line
                    and "atomTypeForceFieldName_ResidueName" in line
                    and "(i.e., atoms_type_per_utilized_FF)" in line
                ):
                    masses_read = True
                    assert len(out_gomc[i + 1].split("!")[0].split()) == 3
                    assert out_gomc[i + 1].split("!")[0].split()[0:3] == [
                        "*",
                        "A",
                        "12.010780",
                    ]
                    assert len(out_gomc[i + 2].split("!")[0].split()) == 3
                    assert out_gomc[i + 2].split("!")[0].split()[0:3] == [
                        "*",
                        "B",
                        "1.007947",
                    ]
                    assert out_gomc[i + 1].split()[4:5] == ["opls_135_ETH"]
                    assert out_gomc[i + 2].split()[4:5] == ["opls_140_ETH"]

                elif (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["A", "B", "340.0", "1.09"],
                        ["A", "A", "268.0", "1.529"],
                    ]
                    assert len(out_gomc[i + 1].split("!")[0].split()) == 4
                    assert len(out_gomc[i + 2].split("!")[0].split()) == 4
                    if (
                        out_gomc[i + 1].split("!")[0].split()[0:4]
                        == bond_types[0]
                    ):
                        assert (
                            out_gomc[i + 1].split("!")[0].split()[0:4]
                            == bond_types[0]
                        )
                        assert (
                            out_gomc[i + 2].split("!")[0].split()[0:4]
                            == bond_types[1]
                        )
                    elif (
                        out_gomc[i + 1].split("!")[0].split()[0:4]
                        == bond_types[1]
                    ):
                        assert (
                            out_gomc[i + 1].split("!")[0].split()[0:4]
                            == bond_types[1]
                        )
                        assert (
                            out_gomc[i + 2].split("!")[0].split()[0:4]
                            == bond_types[0]
                        )

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    angles_read = True
                    angle_types = [
                        ["A", "A", "B", "37.5", "110.70000"],
                        ["B", "A", "B", "33.0", "107.80000"],
                    ]
                    assert len(out_gomc[i + 1].split("!")[0].split()) == 5
                    assert len(out_gomc[i + 2].split("!")[0].split()) == 5
                    if (
                        out_gomc[i + 1].split("!")[0].split()[0:5]
                        == angle_types[0]
                    ):
                        assert (
                            out_gomc[i + 1].split("!")[0].split()[0:5]
                            == angle_types[0]
                        )
                        assert (
                            out_gomc[i + 2].split("!")[0].split()[0:5]
                            == angle_types[1]
                        )
                    elif (
                        out_gomc[i + 1].split("!")[0].split()[0:4]
                        == angle_types[1]
                    ):
                        assert (
                            out_gomc[i + 1].split("!")[0].split()[0:5]
                            == angle_types[1]
                        )
                        assert (
                            out_gomc[i + 2].split("!")[0].split()[0:5]
                            == angle_types[0]
                        )

                elif (
                    "!atom_types" in line
                    and "Kchi" in line
                    and "n" in line
                    and "delta" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    dihedrals_read = True
                    dihed_types = [
                        ["B", "A", "A", "B", "0.000000", "1", "180.0"],
                        ["B", "A", "A", "B", "0.000000", "2", "0.0"],
                        ["B", "A", "A", "B", "-0.150000", "3", "180.0"],
                        ["B", "A", "A", "B", "0.000000", "4", "0.0"],
                        ["B", "A", "A", "B", "0.000000", "5", "180.0"],
                    ]
                    for j in range(0, len(dihed_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 7
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:7]
                            == dihed_types[j]
                        )

                elif (
                    "!atype" in line
                    and "ignored	epsilon" in line
                    and "Rmin/2" in line
                    and "ignored" in line
                    and "eps,1-4" in line
                    and "Rmin/2,1-4" in line
                    and "atom_type_per_utilized_FF" in line
                ):
                    nonbondeds_read = True
                    nb_types = [
                        [
                            "A",
                            "0.00",
                            "-0.066000000",
                            "1.96430858454",
                            "0.00",
                            "-0.033000000",
                            "1.96430858454",
                        ],
                        [
                            "B",
                            "0.00",
                            "-0.030000000",
                            "1.40307756039",
                            "0.00",
                            "-0.015000000",
                            "1.40307756039",
                        ],
                    ]

                    for j in range(0, len(nb_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 7
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:7]
                            == nb_types[j]
                        )

                else:
                    pass

        assert masses_read
        assert bonds_read
        assert angles_read
        assert dihedrals_read
        assert nonbondeds_read

    def test_save_charmm_psf(self, ethane_gomc):
        charmm = Charmm(
            ethane_gomc,
            "charmm_data",
            ff_filename="charmm_data",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_psf()

        with open("charmm_data.psf", "r") as fp:
            charges_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "8 !NATOM" in line:
                    charges_read = True
                    atom_type_charge_etc_list = [
                        [
                            "1",
                            "SYS",
                            "1",
                            "ETH",
                            "C1",
                            "A",
                            "-0.180000",
                            "12.0108",
                        ],
                        [
                            "2",
                            "SYS",
                            "1",
                            "ETH",
                            "C2",
                            "A",
                            "-0.180000",
                            "12.0108",
                        ],
                        [
                            "3",
                            "SYS",
                            "1",
                            "ETH",
                            "H1",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                        [
                            "4",
                            "SYS",
                            "1",
                            "ETH",
                            "H2",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                        [
                            "5",
                            "SYS",
                            "1",
                            "ETH",
                            "H3",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                        [
                            "6",
                            "SYS",
                            "1",
                            "ETH",
                            "H4",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                        [
                            "7",
                            "SYS",
                            "1",
                            "ETH",
                            "H5",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                        [
                            "8",
                            "SYS",
                            "1",
                            "ETH",
                            "H6",
                            "B",
                            "0.060000",
                            "1.0079",
                        ],
                    ]
                    for j in range(0, len(atom_type_charge_etc_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:8]
                            == atom_type_charge_etc_list[j]
                        )

                else:
                    pass

        assert charges_read

    def test_save_charmm_pdb(self, ethane_gomc):
        charmm = Charmm(
            ethane_gomc,
            "charmm_data",
            ff_filename="charmm_data",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_pdb()

        with open("charmm_data.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "ETH", "A", "1"],
                        ["ATOM", "2", "C2", "ETH", "A", "1"],
                        ["ATOM", "3", "H1", "ETH", "A", "1"],
                        ["ATOM", "4", "H2", "ETH", "A", "1"],
                        ["ATOM", "5", "H3", "ETH", "A", "1"],
                        ["ATOM", "6", "H4", "ETH", "A", "1"],
                        ["ATOM", "7", "H5", "ETH", "A", "1"],
                        ["ATOM", "8", "H6", "ETH", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read

    def test_save_charmm_ua_gomc_ff(self, two_propanol_ua):
        charmm = Charmm(
            two_propanol_ua,
            "charmm_data_UA",
            ff_filename="charmm_data_UA",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_inp()

        with open("charmm_data_UA.inp", "r") as fp:
            masses_read = False
            bonds_read = False
            angles_read = False
            dihedrals_read = False
            nonbondeds_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "mass" in line
                    and "atomTypeForceFieldName_ResidueName" in line
                    and "(i.e., atoms_type_per_utilized_FF)" in line
                ):
                    masses_read = True
                    atom_types_1 = [
                        ["*", "A", "15.035000"],
                        ["*", "B", "13.019000"],
                        ["*", "D", "15.999430"],
                        ["*", "C", "1.007947"],
                    ]
                    atom_types_2 = [
                        ["CH3_sp3_POL"],
                        ["CH_O_POL"],
                        ["O_POL"],
                        ["H_POL"],
                    ]

                    for j in range(0, len(atom_types_1)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 3
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:3]
                            == atom_types_1[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[4:5] == atom_types_2[j]
                        )

                elif (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["C", "D", "600.40152964", "0.945"],
                        ["B", "D", "600.40152964", "1.43"],
                        ["A", "B", "600.40152964", "1.54"],
                    ]
                    total_bonds_evaluated = []
                    total_bonds_evaluated_reorg = []
                    for j in range(0, len(bond_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 4
                        )
                        if (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            == bond_types[0]
                            or bond_types[1]
                            or bond_types[2]
                        ):
                            total_bonds_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    for k in range(0, len(bond_types)):
                        if bond_types[k] in total_bonds_evaluated:
                            total_bonds_evaluated_reorg.append(bond_types[k])
                    assert total_bonds_evaluated_reorg == bond_types

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line in line
                ):
                    angles_read = True
                    angle_types = [
                        ["A", "B", "A", "62.10013026", "112.00007"],
                        ["A", "B", "D", "50.07754422", "109.46989"],
                        ["B", "D", "C", "55.04555449", "108.49987"],
                    ]
                    total_angles_evaluated = []
                    total_angles_evaluated_reorg = []
                    for j in range(0, len(angle_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 5
                        )
                        if (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:5]
                            == angle_types[0]
                            or angle_types[1]
                            or angle_types[2]
                        ):
                            total_angles_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:5]
                            )
                    for k in range(0, len(angle_types)):
                        if angle_types[k] in total_angles_evaluated:
                            total_angles_evaluated_reorg.append(angle_types[k])
                    assert total_angles_evaluated_reorg == angle_types

                elif (
                    "!atom_types" in line
                    and "Kchi" in line
                    and "n" in line
                    and "delta" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    dihedrals_read = True
                    dihedral_types = [
                        ["A", "B", "D", "C", "-0.392135", "1", "180.0"],
                        ["A", "B", "D", "C", "-0.062518", "2", "0.0"],
                        ["A", "B", "D", "C", "0.345615", "3", "180.0"],
                        ["A", "B", "D", "C", "0.000000", "4", "0.0"],
                        ["A", "B", "D", "C", "0.000000", "5", "180.0"],
                    ]
                    for j in range(0, len(dihedral_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 7
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:7]
                            == dihedral_types[j]
                        )

                elif (
                    "!atype" in line
                    and "ignored	epsilon" in line
                    and "Rmin/2" in line
                    and "ignored" in line
                    and "eps,1-4" in line
                    and "Rmin/2,1-4" in line
                    and "atom_type_per_utilized_FF" in line
                ):
                    nonbondeds_read = True
                    nb_types = [
                        [
                            "A",
                            "0.00",
                            "-0.194745937",
                            "2.10461634058",
                            "0.00",
                            "-0.000000000",
                            "2.10461634058",
                        ],
                        [
                            "B",
                            "0.00",
                            "-0.019872012",
                            "2.43013033459",
                            "0.00",
                            "-0.000000000",
                            "2.43013033459",
                        ],
                        [
                            "D",
                            "0.00",
                            "-0.184809990",
                            "1.69491769295",
                            "0.00",
                            "-0.000000000",
                            "1.69491769295",
                        ],
                        [
                            "C",
                            "0.00",
                            "-0.000000000",
                            "5.61231024155",
                            "0.00",
                            "-0.000000000",
                            "5.61231024155",
                        ],
                    ]

                    for j in range(0, len(nb_types)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 7
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:7]
                            == nb_types[j]
                        )

                else:
                    pass

        assert masses_read
        assert bonds_read
        assert angles_read
        assert dihedrals_read
        assert nonbondeds_read

    def test_save_charmm_ua_psf(self, two_propanol_ua):
        charmm = Charmm(
            two_propanol_ua,
            "charmm_data_UA",
            ff_filename="charmm_data_UA",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_psf()

        with open("charmm_data_UA.psf", "r") as fp:
            read_psf = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "5 !NATOM" in line:
                    read_psf = True
                    atom_type_charge_etc_list = [
                        [
                            "1",
                            "SYS",
                            "1",
                            "POL",
                            "C1",
                            "A",
                            "0.000000",
                            "15.0350",
                        ],
                        [
                            "2",
                            "SYS",
                            "1",
                            "POL",
                            "BD1",
                            "B",
                            "0.265000",
                            "13.0190",
                        ],
                        [
                            "3",
                            "SYS",
                            "1",
                            "POL",
                            "O1",
                            "D",
                            "-0.700000",
                            "15.9994",
                        ],
                        [
                            "4",
                            "SYS",
                            "1",
                            "POL",
                            "H1",
                            "C",
                            "0.435000",
                            "1.0079",
                        ],
                        [
                            "5",
                            "SYS",
                            "1",
                            "POL",
                            "C2",
                            "A",
                            "0.000000",
                            "15.0350",
                        ],
                    ]

                    for j in range(0, len(atom_type_charge_etc_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:8]
                            == atom_type_charge_etc_list[j]
                        )

                else:
                    pass

        assert read_psf

    def test_save_charmm_ua_pdb(self, two_propanol_ua):
        charmm = Charmm(
            two_propanol_ua,
            "charmm_data_UA",
            ff_filename="charmm_data_UA",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()

        with open("charmm_data_UA.pdb", "r") as fp:
            read_pdb = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    read_pdb = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "POL", "A", "1"],
                        ["ATOM", "2", "BD1", "POL", "A", "1"],
                        ["ATOM", "3", "O1", "POL", "A", "1"],
                        ["ATOM", "4", "H1", "POL", "A", "1"],
                        ["ATOM", "5", "C2", "POL", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "EP"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert read_pdb

    def test_charmm_pdb_fix_angle_bond_fix_atoms(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_propane = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )
        charmm = Charmm(
            test_box_ethane_propane,
            "Test_fixes_angle_bond_atoms",
            ff_filename="Test_fixes_angle_bond_atoms",
            residues=[ethanol_gomc.name, ethane_gomc.name],
            forcefield_selection="oplsaa",
            fix_residue=[ethane_gomc.name],
            fix_residue_in_box=[ethanol_gomc.name],
            gomc_fix_bonds_angles=[ethane_gomc.name],
        )
        charmm.write_inp()
        charmm.write_pdb()

        with open("Test_fixes_angle_bond_atoms.inp", "r") as fp:
            masses_read = False
            bonds_read = False
            angles_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "mass" in line
                    and "atomTypeForceFieldName_ResidueName" in line
                    and "(i.e., atoms_type_per_utilized_FF)" in line
                ):
                    masses_read = True
                    mass_type_1 = [
                        ["*", "A", "12.010780"],
                        ["*", "C", "1.007947"],
                        ["*", "B", "12.010780"],
                        ["*", "G", "12.010780"],
                        ["*", "E", "15.999430"],
                        ["*", "D", "1.007947"],
                        ["*", "F", "1.007947"],
                    ]
                    mass_type_2 = [
                        ["opls_135_ETH"],
                        ["opls_140_ETH"],
                        ["opls_135_ETO"],
                        ["opls_157_ETO"],
                        ["opls_154_ETO"],
                        ["opls_140_ETO"],
                        ["opls_155_ETO"],
                    ]

                    for j in range(0, len(mass_type_1)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 3
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:3]
                            == mass_type_1[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[4:5] == mass_type_2[j]
                        )

                elif (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["D", "G", "340.0", "1.09"],
                        ["E", "G", "320.0", "1.41"],
                        ["E", "F", "553.0", "0.945"],
                        ["A", "C", "999999999999", "1.09"],
                        ["B", "D", "340.0", "1.09"],
                        ["A", "A", "999999999999", "1.529"],
                        ["B", "G", "268.0", "1.529"],
                    ]
                    total_bonds_evaluated = []
                    total_fixed_bonds = []
                    for j in range(0, 7):
                        total_bonds_evaluated.append(
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                        )
                        if out_gomc[i + 1 + j].split("!")[0].split()[2:3] == [
                            "999999999999"
                        ]:
                            total_fixed_bonds.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert total_bonds_evaluated.sort() == bond_types.sort()
                    assert len(total_fixed_bonds) == 2

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line in line
                ):
                    angles_read = True
                    fixed_angle_types = [
                        ["A", "A", "C", "999999999999", "110.70000"],
                        ["C", "A", "C", "999999999999", "107.80000"],
                    ]
                    total_angles_evaluated = []
                    total_fixed_angles = []
                    for j in range(0, 9):
                        if out_gomc[i + 1 + j].split("!")[0].split()[0:4] == (
                            fixed_angle_types[0] or fixed_angle_types[1]
                        ):
                            total_angles_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                        if out_gomc[i + 1 + j].split("!")[0].split()[3:4] == [
                            "999999999999"
                        ]:
                            total_fixed_angles.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert (
                        fixed_angle_types.sort()
                        == total_angles_evaluated.sort()
                    )
                    assert len(total_fixed_angles) == len(fixed_angle_types)

                else:
                    pass

        assert masses_read
        assert bonds_read
        assert angles_read

        with open("Test_fixes_angle_bond_atoms.pdb", "r") as fp:
            read_pdb_part_1 = False
            read_pdb_part_2 = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    read_pdb_part_1 = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "20.000",
                        "20.000",
                        "20.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]

                if "CRYST1" in line:
                    read_pdb_part_2 = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "ETH", "A", "1"],
                        ["ATOM", "2", "C2", "ETH", "A", "1"],
                        ["ATOM", "3", "H1", "ETH", "A", "1"],
                        ["ATOM", "4", "H2", "ETH", "A", "1"],
                        ["ATOM", "5", "H3", "ETH", "A", "1"],
                        ["ATOM", "6", "H4", "ETH", "A", "1"],
                        ["ATOM", "7", "H5", "ETH", "A", "1"],
                        ["ATOM", "8", "H6", "ETH", "A", "1"],
                        ["ATOM", "9", "C1", "ETO", "A", "2"],
                        ["ATOM", "10", "C2", "ETO", "A", "2"],
                        ["ATOM", "11", "O1", "ETO", "A", "2"],
                        ["ATOM", "12", "H1", "ETO", "A", "2"],
                        ["ATOM", "13", "H2", "ETO", "A", "2"],
                        ["ATOM", "14", "H3", "ETO", "A", "2"],
                        ["ATOM", "15", "H4", "ETO", "A", "2"],
                        ["ATOM", "16", "H5", "ETO", "A", "2"],
                        ["ATOM", "17", "H6", "ETO", "A", "2"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "1.00", "C"],
                        ["1.00", "1.00", "C"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "1.00", "H"],
                        ["1.00", "2.00", "C"],
                        ["1.00", "2.00", "C"],
                        ["1.00", "2.00", "O"],
                        ["1.00", "2.00", "H"],
                        ["1.00", "2.00", "H"],
                        ["1.00", "2.00", "H"],
                        ["1.00", "2.00", "H"],
                        ["1.00", "2.00", "H"],
                        ["1.00", "2.00", "H"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert read_pdb_part_1
        assert read_pdb_part_2

    def test_charmm_pdb_fix_bonds_only(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_propane = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )
        charmm = Charmm(
            test_box_ethane_propane,
            "Test_fixes_bonds_only",
            ff_filename="Test_fixes_bonds_only",
            residues=[ethanol_gomc.name, ethane_gomc.name],
            forcefield_selection="oplsaa",
            gomc_fix_bonds=[ethane_gomc.name],
        )
        charmm.write_inp()

        with open("Test_fixes_bonds_only.inp", "r") as fp:
            bonds_read = False
            angles_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["D", "G", "340.0", "1.09"],
                        ["E", "G", "320.0", "1.41"],
                        ["E", "F", "553.0", "0.945"],
                        ["A", "C", "999999999999", "1.09"],
                        ["B", "D", "340.0", "1.09"],
                        ["A", "A", "999999999999", "1.529"],
                        ["B", "G", "268.0", "1.529"],
                    ]
                    total_bonds_evaluated = []
                    total_fixed_bonds = []
                    for j in range(0, 7):
                        total_bonds_evaluated.append(
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                        )
                        if out_gomc[i + 1 + j].split("!")[0].split()[2:3] == [
                            "999999999999"
                        ]:
                            total_fixed_bonds.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert total_bonds_evaluated.sort() == bond_types.sort()
                    assert len(total_fixed_bonds) == 2

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    angles_read = True
                    fixed_angle_types = []
                    total_angles_evaluated = []
                    total_fixed_angles = []
                    for j in range(0, 9):
                        if len(fixed_angle_types) > 0:
                            if out_gomc[i + 1 + j].split("!")[0].split()[
                                0:4
                            ] == (fixed_angle_types[0] or fixed_angle_types[1]):
                                total_angles_evaluated.append(
                                    out_gomc[i + 1 + j]
                                    .split("!")[0]
                                    .split()[0:4]
                                )
                        if out_gomc[i + 1 + j].split("!")[0].split()[3:4] == [
                            "999999999999"
                        ]:
                            total_fixed_angles.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert (
                        fixed_angle_types.sort()
                        == total_angles_evaluated.sort()
                    )
                    assert len(total_fixed_angles) == len(fixed_angle_types)

                else:
                    pass

        assert bonds_read
        assert angles_read

    def test_charmm_pdb_fix_bonds_only_and_fix_bonds_angles(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_propane = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )
        charmm = Charmm(
            test_box_ethane_propane,
            "Test_fixes_bonds_only_and_fix_bonds_angles",
            ff_filename="Test_fixes_bonds_only_and_fix_bonds_angles",
            residues=[ethanol_gomc.name, ethane_gomc.name],
            forcefield_selection="oplsaa",
            gomc_fix_bonds=[ethane_gomc.name],
            gomc_fix_bonds_angles=[ethane_gomc.name],
        )
        charmm.write_inp()

        with open("Test_fixes_bonds_only_and_fix_bonds_angles.inp", "r") as fp:
            bonds_read = False
            angles_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["D", "G", "340.0", "1.09"],
                        ["E", "G", "320.0", "1.41"],
                        ["E", "F", "553.0", "0.945"],
                        ["A", "C", "999999999999", "1.09"],
                        ["B", "D", "340.0", "1.09"],
                        ["A", "A", "999999999999", "1.529"],
                        ["B", "G", "268.0", "1.529"],
                    ]
                    total_bonds_evaluated = []
                    total_fixed_bonds = []
                    for j in range(0, 7):
                        total_bonds_evaluated.append(
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                        )
                        if out_gomc[i + 1 + j].split("!")[0].split()[2:3] == [
                            "999999999999"
                        ]:
                            total_fixed_bonds.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert total_bonds_evaluated.sort() == bond_types.sort()
                    assert len(total_fixed_bonds) == 2

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    angles_read = True
                    fixed_angle_types = [
                        ["A", "A", "C", "999999999999", "110.70000"],
                        ["C", "A", "C", "999999999999", "107.80000"],
                    ]
                    total_angles_evaluated = []
                    total_fixed_angles = []
                    for j in range(0, 9):
                        if out_gomc[i + 1 + j].split("!")[0].split()[0:4] == (
                            fixed_angle_types[0] or fixed_angle_types[1]
                        ):
                            total_angles_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                        if out_gomc[i + 1 + j].split("!")[0].split()[3:4] == [
                            "999999999999"
                        ]:
                            total_fixed_angles.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert (
                        fixed_angle_types.sort()
                        == total_angles_evaluated.sort()
                    )
                    assert len(total_fixed_angles) == len(fixed_angle_types)

                else:
                    pass

        assert bonds_read
        assert angles_read

    def test_charmm_pdb_fix_angles_only(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_propane = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )
        charmm = Charmm(
            test_box_ethane_propane,
            "Test_fixes_angles_only",
            ff_filename="Test_fixes_angles_only",
            residues=[ethanol_gomc.name, ethane_gomc.name],
            forcefield_selection="oplsaa",
            gomc_fix_angles=[ethane_gomc.name],
        )
        charmm.write_inp()

        with open("Test_fixes_angles_only.inp", "r") as fp:
            bonds_read = False
            angles_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["D", "G", "340.0", "1.09"],
                        ["E", "G", "320.0", "1.41"],
                        ["E", "F", "553.0", "0.945"],
                        ["A", "C", "340.0", "1.09"],
                        ["B", "D", "340.0", "1.09"],
                        ["A", "A", "268.0", "1.529"],
                        ["B", "G", "268.0", "1.529"],
                    ]
                    total_bonds_evaluated = []
                    total_fixed_bonds = []
                    for j in range(0, 7):
                        total_bonds_evaluated.append(
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                        )
                        if out_gomc[i + 1 + j].split("!")[0].split()[2:3] == [
                            "999999999999"
                        ]:
                            total_fixed_bonds.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert total_bonds_evaluated.sort() == bond_types.sort()
                    assert len(total_fixed_bonds) == 0

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    angles_read = True
                    fixed_angle_types = [
                        ["A", "A", "C", "999999999999", "110.70000"],
                        ["C", "A", "C", "999999999999", "107.80000"],
                    ]
                    total_angles_evaluated = []
                    total_fixed_angles = []
                    for j in range(0, 9):
                        if out_gomc[i + 1 + j].split("!")[0].split()[0:4] == (
                            fixed_angle_types[0] or fixed_angle_types[1]
                        ):
                            total_angles_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                        if out_gomc[i + 1 + j].split("!")[0].split()[3:4] == [
                            "999999999999"
                        ]:
                            total_fixed_angles.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert (
                        fixed_angle_types.sort()
                        == total_angles_evaluated.sort()
                    )
                    assert len(total_fixed_angles) == len(fixed_angle_types)

                else:
                    pass

        assert bonds_read
        assert angles_read

    def test_charmm_pdb_fix_angles_only_and_fix_bonds_angles(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_propane = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )
        charmm = Charmm(
            test_box_ethane_propane,
            "Test_fixes_angles_only_and_fix_bonds_angles",
            ff_filename="Test_fixes_angles_only_and_fix_bonds_angles",
            residues=[ethanol_gomc.name, ethane_gomc.name],
            forcefield_selection="oplsaa",
            gomc_fix_angles=[ethane_gomc.name],
            gomc_fix_bonds_angles=[ethane_gomc.name],
        )
        charmm.write_inp()

        with open("Test_fixes_angles_only_and_fix_bonds_angles.inp", "r") as fp:
            bonds_read = False
            angles_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "Kb" in line
                    and "b0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    bonds_read = True
                    bond_types = [
                        ["D", "G", "340.0", "1.09"],
                        ["E", "G", "320.0", "1.41"],
                        ["E", "F", "553.0", "0.945"],
                        ["A", "C", "999999999999", "1.09"],
                        ["B", "D", "340.0", "1.09"],
                        ["A", "A", "999999999999", "1.529"],
                        ["B", "G", "268.0", "1.529"],
                    ]
                    total_bonds_evaluated = []
                    total_fixed_bonds = []
                    for j in range(0, 7):
                        total_bonds_evaluated.append(
                            out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                        )
                        if out_gomc[i + 1 + j].split("!")[0].split()[2:3] == [
                            "999999999999"
                        ]:
                            total_fixed_bonds.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert total_bonds_evaluated.sort() == bond_types.sort()
                    assert len(total_fixed_bonds) == 2

                elif (
                    "!atom_types" in line
                    and "Ktheta" in line
                    and "Theta0" in line
                    and "atoms_types_per_utilized_FF" in line
                ):
                    angles_read = True
                    fixed_angle_types = [
                        ["A", "A", "C", "999999999999", "110.70000"],
                        ["C", "A", "C", "999999999999", "107.80000"],
                    ]
                    total_angles_evaluated = []
                    total_fixed_angles = []
                    for j in range(0, 9):
                        if out_gomc[i + 1 + j].split("!")[0].split()[0:4] == (
                            fixed_angle_types[0] or fixed_angle_types[1]
                        ):
                            total_angles_evaluated.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                        if out_gomc[i + 1 + j].split("!")[0].split()[3:4] == [
                            "999999999999"
                        ]:
                            total_fixed_angles.append(
                                out_gomc[i + 1 + j].split("!")[0].split()[0:4]
                            )
                    assert (
                        fixed_angle_types.sort()
                        == total_angles_evaluated.sort()
                    )
                    assert len(total_fixed_angles) == len(fixed_angle_types)

                else:
                    pass

        assert bonds_read
        assert angles_read

    def test_charmm_pdb_no_differenc_1_4_coul_scalars(
        self, two_propanol_ua, ethane_gomc
    ):
        test_box_ethane_two_propanol_ua = mb.fill_box(
            compound=[two_propanol_ua, ethane_gomc],
            n_compounds=[1, 1],
            box=[2.0, 2.0, 2.0],
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: There are multiple 1,4-coulombic scaling factors "
            "GOMC will only accept a singular input for the 1,4-coulombic "
            "scaling factors",
        ):
            Charmm(
                test_box_ethane_two_propanol_ua,
                "residue_reorder_box_sizing_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="residue_reorder_box_sizing_box_1",
                ff_filename="residue_reorder_box",
                residues=[two_propanol_ua.name, ethane_gomc.name],
                forcefield_selection={
                    two_propanol_ua.name: "trappe-ua",
                    ethane_gomc.name: "oplsaa",
                },
                fix_residue=None,
                fix_residue_in_box=None,
                gomc_fix_bonds_angles=None,
                reorder_res_in_pdb_psf=False,
                bead_to_atom_name_dict={"_CH3": "C"},
            )

    def test_charmm_pdb_residue_reorder_and_ff_filename_box_sizing(
        self, ethanol_gomc, ethane_gomc
    ):
        test_box_ethane_ethanol_gomc = mb.fill_box(
            compound=[ethanol_gomc, ethane_gomc],
            n_compounds=[1, 1],
            box=[3, 3, 3],
        )

        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[4, 4, 4]
        )

        charmm = Charmm(
            test_box_ethane_ethanol_gomc,
            "residue_reorder_box_sizing_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="residue_reorder_box_sizing_box_1",
            ff_filename=None,
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection=str(forcefields.get_ff_path()[0])
            + "/xml/"
            + "oplsaa.xml",
            fix_residue=None,
            fix_residue_in_box=None,
            gomc_fix_bonds_angles=None,
            reorder_res_in_pdb_psf=True,
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()

        with open("residue_reorder_box_sizing_box_0.pdb", "r") as fp:
            pdb_box_0_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_box_0_read = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "30.000",
                        "30.000",
                        "30.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]

                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "ETH", "A", "1"],
                        ["ATOM", "2", "C2", "ETH", "A", "1"],
                        ["ATOM", "3", "H1", "ETH", "A", "1"],
                        ["ATOM", "4", "H2", "ETH", "A", "1"],
                        ["ATOM", "5", "H3", "ETH", "A", "1"],
                        ["ATOM", "6", "H4", "ETH", "A", "1"],
                        ["ATOM", "7", "H5", "ETH", "A", "1"],
                        ["ATOM", "8", "H6", "ETH", "A", "1"],
                        ["ATOM", "9", "C1", "ETO", "A", "2"],
                        ["ATOM", "10", "C2", "ETO", "A", "2"],
                        ["ATOM", "11", "O1", "ETO", "A", "2"],
                        ["ATOM", "12", "H1", "ETO", "A", "2"],
                        ["ATOM", "13", "H2", "ETO", "A", "2"],
                        ["ATOM", "14", "H3", "ETO", "A", "2"],
                        ["ATOM", "15", "H4", "ETO", "A", "2"],
                        ["ATOM", "16", "H5", "ETO", "A", "2"],
                        ["ATOM", "17", "H6", "ETO", "A", "2"],
                    ]

                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "C"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "H"],
                    ]
                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_box_0_read

        with open("residue_reorder_box_sizing_box_1.pdb", "r") as fp:
            pdb_box_1_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_box_1_read = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "40.000",
                        "40.000",
                        "40.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]
                else:
                    pass

        assert pdb_box_1_read

    # test utils base 10 to base 16 converter
    def test_base_10_to_base_16(self):
        list_base_10_and_16 = [
            [15, "f"],
            [16, "10"],
            [17, "11"],
            [200, "c8"],
            [1000, "3e8"],
            [5000, "1388"],
            [int(16**3 - 1), "fff"],
            [int(16**3), "1000"],
        ]

        for test_base_16_iter in range(0, len(list_base_10_and_16)):
            test_10_iter = list_base_10_and_16[test_base_16_iter][0]
            test_16_iter = list_base_10_and_16[test_base_16_iter][1]
            assert str(base10_to_base16_alph_num(test_10_iter)) == str(
                test_16_iter
            )

        unique_entries_base_16_list = []
        for test_unique_base_16 in range(0, 16**2):
            unique_entries_base_16_list.append(
                base10_to_base16_alph_num(test_unique_base_16)
            )

        verified_unique_entries_base_16_list = np.unique(
            unique_entries_base_16_list
        )
        assert len(verified_unique_entries_base_16_list) == len(
            unique_entries_base_16_list
        )

        add_same_values_list = ["1", "a"]
        for add_same_base_16 in range(0, len(add_same_values_list)):
            verified_unique_entries_base_16_list = np.append(
                verified_unique_entries_base_16_list,
                add_same_values_list[add_same_base_16],
            )
        assert len(verified_unique_entries_base_16_list) - len(
            add_same_values_list
        ) == len(unique_entries_base_16_list)

    # test utils base 10 to base 26 converter
    def test_base_10_to_base_26(self):
        list_base_10_and_26 = [
            [0, "A"],
            [5, "F"],
            [25, "Z"],
            [26, "BA"],
            [200, "HS"],
            [1000, "BMM"],
            [5000, "HKI"],
            [int(26**3 - 1), "ZZZ"],
            [int(26**3), "BAAA"],
        ]

        for test_base_26_iter in range(0, len(list_base_10_and_26)):
            test_10_iter = list_base_10_and_26[test_base_26_iter][0]
            test_26_iter = list_base_10_and_26[test_base_26_iter][1]
            assert str(base10_to_base26_alph(test_10_iter)) == str(test_26_iter)

        unique_entries_base_26_list = []
        for test_unique_base_26 in range(0, 26**2):
            unique_entries_base_26_list.append(
                base10_to_base26_alph(test_unique_base_26)
            )

        verified_unique_entries_base_26_list = np.unique(
            unique_entries_base_26_list
        )
        assert len(verified_unique_entries_base_26_list) == len(
            unique_entries_base_26_list
        )

        add_same_values_list = ["1", "a"]
        for add_same_base_26 in range(0, len(add_same_values_list)):
            verified_unique_entries_base_26_list = np.append(
                verified_unique_entries_base_26_list,
                add_same_values_list[add_same_base_26],
            )
        assert len(verified_unique_entries_base_26_list) - len(
            add_same_values_list
        ) == len(unique_entries_base_26_list)

    # test utils base 10 to base 52 converter
    def test_base_10_to_base_52(self):
        list_base_10_and_52 = [
            [17, "R"],
            [51, "z"],
            [52, "BA"],
            [53, "BB"],
            [200, "Ds"],
            [1000, "TM"],
            [5000, "BsI"],
            [int(52**3 - 1), "zzz"],
            [int(52**3), "BAAA"],
        ]

        for test_base_52_iter in range(0, len(list_base_10_and_52)):
            test_10_iter = list_base_10_and_52[test_base_52_iter][0]
            test_52_iter = list_base_10_and_52[test_base_52_iter][1]
            assert str(base10_to_base52_alph(test_10_iter)) == str(test_52_iter)

        unique_entries_base_52_list = []
        for test_unique_base_52 in range(0, 52**2):
            unique_entries_base_52_list.append(
                base10_to_base52_alph(test_unique_base_52)
            )

        verified_unique_entries_base_52_list = np.unique(
            unique_entries_base_52_list
        )
        assert len(verified_unique_entries_base_52_list) == len(
            unique_entries_base_52_list
        )

        add_same_values_list = ["1", "a"]
        for add_same_base_52 in range(0, len(add_same_values_list)):
            verified_unique_entries_base_52_list = np.append(
                verified_unique_entries_base_52_list,
                add_same_values_list[add_same_base_52],
            )
        assert len(verified_unique_entries_base_52_list) - len(
            add_same_values_list
        ) == len(unique_entries_base_52_list)

    # test utils base 10 to base 62 converter
    def test_base_10_to_base_62(self):
        list_base_10_and_62 = [
            [17, "H"],
            [61, "z"],
            [62, "10"],
            [63, "11"],
            [200, "3E"],
            [1000, "G8"],
            [5000, "1Ie"],
            [int(62**3 - 1), "zzz"],
            [int(62**3), "1000"],
        ]

        for test_base_62_iter in range(0, len(list_base_10_and_62)):
            test_10_iter = list_base_10_and_62[test_base_62_iter][0]
            test_62_iter = list_base_10_and_62[test_base_62_iter][1]
            assert str(base10_to_base62_alph_num(test_10_iter)) == str(
                test_62_iter
            )

        unique_entries_base_62_list = []
        for test_unique_base_62 in range(0, 62**2):
            unique_entries_base_62_list.append(
                base10_to_base62_alph_num(test_unique_base_62)
            )

        verified_unique_entries_base_62_list = np.unique(
            unique_entries_base_62_list
        )
        assert len(verified_unique_entries_base_62_list) == len(
            unique_entries_base_62_list
        )

        add_same_values_list = ["1", "a"]
        for add_same_base_62 in range(0, len(add_same_values_list)):
            verified_unique_entries_base_62_list = np.append(
                verified_unique_entries_base_62_list,
                add_same_values_list[add_same_base_62],
            )
        assert len(verified_unique_entries_base_62_list) - len(
            add_same_values_list
        ) == len(unique_entries_base_62_list)

    # Tests for the mbuild.utils.specific_FF_to_residue.Specific_FF_to_residue() function
    def test_specific_ff_ff_is_none(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"Please the force field selection \(forcefield_selection\) as a "
            r"dictionary with all the residues specified to a force field "
            '-> Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection=None,
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_ff_extention(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct "
            r"foyer FF name and not a path to a FF file. "
            r"If you are entering a path to a FF file, "
            r"please use the forcefield_files variable with the "
            r"proper XML extension \(.xml\).",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa.pdb"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_all_residue_not_input(self, ethane_gomc, ethanol_gomc):
        with pytest.raises(
            ValueError,
            match=r"All the residues are not specified, or the residues "
            r"entered does not match the residues that were found "
            r"and built for structure.",
        ):
            box = mb.fill_box(
                compound=[ethane_gomc, ethanol_gomc],
                box=[1, 1, 1],
                n_compounds=[1, 1],
            )

            specific_ff_to_residue(
                box,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=2,
            )

    def test_specific_ff_to_residue_ff_selection_not_dict(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"The force field selection \(forcefield_selection\) "
            "is not a dictionary. Please enter a dictionary "
            "with all the residues specified to a force field "
            '-> Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file "
            "or by using the standard force field name provided the `foyer` package.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection="oplsaa",
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_is_none(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"Please enter the residues in the Specific_FF_to_residue function.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=None,
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_reorder_not_true_or_false(
        self, ethane_gomc
    ):
        with pytest.raises(
            TypeError,
            match=r"Please enter the reorder_res_in_pdb_psf "
            r"in the Specific_FF_to_residue function \(i.e., True or False\).",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=None,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_simulation_boxes_not_1_or_2(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please enter boxes_for_simulation equal the integer 1 or 2.",
        ):
            test_box_ethane_gomc = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[2, 3, 4]
            )

            specific_ff_to_residue(
                test_box_ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=3,
            )

    def test_specific_ff_to_residue_ffselection_wrong_path(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct foyer FF path, "
            r"including the FF file name.xml "
            r"If you are using the pre-build FF files in foyer, "
            r"only use the string name without any extension.",
        ):
            test_box_ethane_gomc = mb.fill_box(
                compound=[ethane_gomc], n_compounds=[1], box=[4, 5, 6]
            )

            specific_ff_to_residue(
                test_box_ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa.xml"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_path(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct foyer FF path, "
            r"including the FF file name.xml "
            r"If you are using the pre-build FF files in foyer, "
            r"only use the string name without any extension.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa.xml"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_input_string_as_compound(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The structure expected to be of type: "
            r"<class 'mbuild.compound.Compound'> or <class 'mbuild.box.Box'>, "
            r"received: <class 'str'>",
        ):
            specific_ff_to_residue(
                "ethane_gomc",
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_boxes_for_simulation_not_int(
        self, ethane_gomc
    ):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter boxes_for_simulation equal "
            "the integer 1 or 2.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1.1,
            )

    def test_specific_ff_to_residues_no_ff(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"The forcefield_selection variable are not provided, "
            r"but there are residues provided.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_no_residues(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"The residues variable is an empty list but there are "
            "forcefield_selection variables provided.",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
                residues=[],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_wrong_foyer_name(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"Please make sure you are entering the correct foyer FF name, "
            r"or the correct file extension \(i.e., .xml, if required\).",
        ):
            specific_ff_to_residue(
                ethane_gomc,
                forcefield_selection={ethane_gomc.name: "xxx"},
                residues=[ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_specific_ff_to_residue_ffselection_run(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[4, 5, 6]
        )

        [
            test_value_0,
            test_value_1,
            test_value_2,
            test_value_3,
        ] = specific_ff_to_residue(
            test_box_ethane_gomc,
            forcefield_selection={
                ethane_gomc.name: forcefields.get_ff_path()[0]
                + "/xml/"
                + "oplsaa.xml"
            },
            residues=[ethane_gomc.name],
            reorder_res_in_pdb_psf=False,
            boxes_for_simulation=1,
        )
        assert test_value_1 == {"ETH": 0.5}
        assert test_value_2 == {"ETH": 0.5}
        assert test_value_3 == ["ETH"]

    def test_specific_ff_to_no_atoms_in_residue(self):
        with pytest.raises(
            ValueError,
            match=r"The residues variable is an empty list but there "
            r"are forcefield_selection variables provided.",
        ):
            empty_compound = mb.Compound()

            specific_ff_to_residue(
                empty_compound,
                forcefield_selection={"empty_compound": "oplsaa"},
                residues=[],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_charmm_methane_test_no_children(self, methane_ua_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: If you are not providing an empty box, "
            r"you need to specify the atoms/beads as children in the mb.Compound. "
            r"If you are providing and empty box, please do so by specifying and "
            r"mbuild Box \({}\)".format(type(Box(lengths=[1, 1, 1]))),
        ):
            specific_ff_to_residue(
                methane_ua_gomc,
                forcefield_selection={methane_ua_gomc.name: "trappe-ua"},
                residues=[methane_ua_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_charmm_a_few_mbuild_layers(self, ethane_gomc, ethanol_gomc):
        box_reservior_1 = mb.fill_box(
            compound=[ethane_gomc], box=[1, 1, 1], n_compounds=[1]
        )
        box_reservior_1.periodicity = (True, True, True)
        box_reservior_2 = mb.fill_box(
            compound=[ethanol_gomc], box=[1, 1, 1], n_compounds=[1]
        )
        box_reservior_2.translate([0, 0, 1])
        box_reservior_1.add(box_reservior_2, inherit_periodicity=False)

        [
            test_value_0,
            test_value_1,
            test_value_2,
            test_value_3,
        ] = specific_ff_to_residue(
            box_reservior_1,
            forcefield_selection={
                ethanol_gomc.name: "oplsaa",
                ethane_gomc.name: "oplsaa",
            },
            residues=[ethanol_gomc.name, ethane_gomc.name],
            reorder_res_in_pdb_psf=False,
            boxes_for_simulation=1,
        )

        assert (
            str(test_value_0)
            == "<Structure 17 atoms; 2 residues; 15 bonds; PBC (orthogonal); parametrized>"
        )
        assert test_value_1 == {"ETO": 0.5, "ETH": 0.5}
        assert test_value_2 == {"ETO": 0.5, "ETH": 0.5}
        assert test_value_3 == ["ETH", "ETO"]

    def test_charmm_all_residues_not_in_dict(self, ethane_gomc, ethanol_gomc):
        with pytest.raises(
            ValueError,
            match=r"All the residues were not used from the forcefield_selection "
            r"string or dictionary. There may be residues below other "
            r"specified residues in the mbuild.Compound hierarchy. "
            r"If so, all the highest listed residues pass down the force "
            r"fields through the hierarchy. Alternatively, residues that "
            r"are not in the structure may have been specified. ",
        ):
            box_reservior_1 = mb.fill_box(
                compound=[ethane_gomc], box=[1, 1, 1], n_compounds=[1]
            )
            specific_ff_to_residue(
                box_reservior_1,
                forcefield_selection={ethanol_gomc.name: "oplsaa"},
                residues=[ethanol_gomc.name, ethane_gomc.name],
                reorder_res_in_pdb_psf=False,
                boxes_for_simulation=1,
            )

    def test_charmm_correct_residue_format(self, ethane_gomc):
        test_value = Charmm(
            ethane_gomc,
            "box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename=None,
            residues=[ethane_gomc.name],
            forcefield_selection={ethane_gomc.name: "oplsaa"},
        )

        assert test_value.input_error is False

    def test_charmm_residue_not_list(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the residues list \(residues\) in a list format.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename=None,
                residues=ethane_gomc.name,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_residue_string(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the residues list \(residues\) in a list format.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename=None,
                residues="ethane_gomc.name",
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_residue_is_none(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the residues list \(residues\)",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename=None,
                residues=None,
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_filename_0_is_not_string(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the filename_box_0 as a string.",
        ):
            Charmm(
                ethane_gomc,
                0,
                structure_box_1=None,
                filename_box_1=None,
                ff_filename=None,
                residues=[ethane_gomc.name],
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_filename_box_1_is_not_string(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the filename_box_1 as a string.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=ethane_gomc,
                filename_box_1=["box_0"],
                ff_filename=None,
                residues=[ethane_gomc.name],
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_gomc_filename_not_string(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter GOMC force field name \(ff_filename\) as a string.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename=0,
                residues=[ethane_gomc.name],
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_gomc_filename_ext_not_dot_inp(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter GOMC force field name without an "
            "extention or the .inp extension.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename="box.test",
                residues=[ethane_gomc.name],
                forcefield_selection={ethane_gomc.name: "oplsaa"},
            )

    def test_charmm_ffselection_not_dict(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The force field selection \(forcefield_selection\) "
            "is not a string or a dictionary with all the residues specified "
            'to a force field. -> String Ex: "path/trappe-ua.xml" or Ex: "trappe-ua" '
            "Otherise provided a dictionary with all the residues specified "
            "to a force field "
            '->Dictionary Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file if "
            "a standard foyer force field is not used.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename="box_0",
                residues=[ethane_gomc.name],
                forcefield_selection=["oplsaa", "oplsaa"],
            )

    def test_charmm_ffselection_string(self, ethane_gomc):
        test_value = Charmm(
            ethane_gomc,
            "box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="box_0",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        assert test_value.input_error is False

    def test_charmm_residue_name_not_in_residues(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the residues are not specified, or "
            "the residues entered does not match the residues that "
            "were found and built for structure.",
        ):
            Charmm(
                ethane_gomc,
                "box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename="box_0",
                residues=["XXX"],
                forcefield_selection="oplsaa",
            )

    def test_ffselection_string(self, two_propanol_ua):
        charmm = Charmm(
            two_propanol_ua,
            "ffselection_string",
            ff_filename="ffselection_string",
            residues=[two_propanol_ua.name],
            forcefield_selection=forcefields.get_ff_path()[0]
            + "/xml/"
            + "trappe-ua.xml",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()

        with open("ffselection_string.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "POL", "A", "1"],
                        ["ATOM", "2", "BD1", "POL", "A", "1"],
                        ["ATOM", "3", "O1", "POL", "A", "1"],
                        ["ATOM", "4", "H1", "POL", "A", "1"],
                        ["ATOM", "5", "C2", "POL", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "EP"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read

    def test_ff_selection_list(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The force field selection \(forcefield_selection\) "
            "is not a string or a dictionary with all the residues specified "
            'to a force field. -> String Ex: "path/trappe-ua.xml" or Ex: "trappe-ua" '
            "Otherise provided a dictionary with all the residues specified "
            "to a force field "
            '->Dictionary Ex: {"Water" : "oplsaa", "OCT": "path/trappe-ua.xml"}, '
            "Note: the file path must be specified the force field file if "
            "a standard foyer force field is not used.",
        ):
            Charmm(
                two_propanol_ua,
                "S",
                ff_filename="S",
                residues=[two_propanol_ua.name],
                forcefield_selection=[
                    str(forcefields.get_ff_path()[0])
                    + "/xml/"
                    + "trappe-ua.xml"
                ],
                bead_to_atom_name_dict={"_CH3": "C"},
            )

    def test_residues_not_a_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter a residues list "
            r"\(residues\) with only string values.",
        ):
            Charmm(
                two_propanol_ua,
                "box_0",
                ff_filename="box_0",
                residues=[2],
                forcefield_selection={two_propanol_ua.name: "trappe-ua"},
                bead_to_atom_name_dict={"_CH3": "C"},
            )

    # charmm writer sub-function testing
    def test_charmm_bond_reorder_angle_urey_bradleys(
        self, two_propanol_gomc, ethanol_gomc
    ):
        box_reservior_0 = mb.fill_box(
            compound=[two_propanol_gomc, ethanol_gomc],
            box=[2, 2, 2],
            n_compounds=[2, 2],
        )

        [
            structure_ff,
            coulomb14scalar_dict,
            lj14_scalar_dict,
            residues_applied_list,
        ] = specific_ff_to_residue(
            box_reservior_0,
            forcefield_selection={
                two_propanol_gomc.name: "oplsaa",
                ethanol_gomc.name: "oplsaa",
            },
            residues=[ethanol_gomc.name, two_propanol_gomc.name],
            reorder_res_in_pdb_psf=False,
            boxes_for_simulation=1,
        )

        sigma_conversion_factor = 1
        epsilon_conversion_factor = 1
        # reversed the bond order so it fixes itself
        bonds_1 = [
            [bond.atom1.idx + 1, bond.atom2.idx + 1]
            for bond in structure_ff.bonds
        ]
        bond_types_1, unique_bond_types_1 = charmm_writer._get_bond_types(
            structure_ff, sigma_conversion_factor, epsilon_conversion_factor
        )

        bonds_2 = [
            [bond.atom2.idx + 1, bond.atom1.idx + 1]
            for bond in structure_ff.bonds
        ]
        bond_types_2, unique_bond_types_2 = charmm_writer._get_bond_types(
            structure_ff, sigma_conversion_factor, epsilon_conversion_factor
        )

        assert bonds_1 != bonds_2
        assert bond_types_1 == bond_types_2
        assert unique_bond_types_1 == unique_bond_types_2

        # test for error if trying to use urey_bradleys in th angles
        use_urey_bradleys = True
        angle_types_1, unique_angle_types_1 = charmm_writer._get_angle_types(
            structure_ff,
            sigma_conversion_factor,
            epsilon_conversion_factor,
            use_urey_bradleys=use_urey_bradleys,
        )

        assert angle_types_1 is None
        assert unique_angle_types_1 is None

    # test for error if trying to use  use_dihedrals and impropers in the dihedrals (i.e. only RB torsion allowed)
    def test_charmm_dihedral_reorder(self, ethyl_ether_gomc, methyl_ether_gomc):

        box_reservior_0 = mb.fill_box(
            compound=[ethyl_ether_gomc, methyl_ether_gomc],
            box=[10, 10, 10],
            n_compounds=[10, 10],
        )

        [
            structure_ff,
            coulomb14scalar_dict,
            lj14_scalar_dict,
            residues_applied_list,
        ] = specific_ff_to_residue(
            box_reservior_0,
            forcefield_selection={
                ethyl_ether_gomc.name: "oplsaa",
                methyl_ether_gomc.name: "oplsaa",
            },
            residues=[ethyl_ether_gomc.name, methyl_ether_gomc.name],
            reorder_res_in_pdb_psf=False,
            boxes_for_simulation=1,
        )

        use_rb_torsions_1 = False
        use_dihedrals_1 = True
        epsilon_conversion_factor = 1
        lj_unit = 1 / epsilon_conversion_factor
        (
            dihedral_types_1,
            unique_dihedral_types_1,
        ) = charmm_writer._get_dihedral_types(
            structure_ff,
            use_rb_torsions_1,
            use_dihedrals_1,
            epsilon_conversion_factor,
        )
        assert dihedral_types_1 is None
        assert unique_dihedral_types_1 is None

        use_rb_torsions_2 = True
        use_dihedrals_2 = False

        (
            dihedral_types_2,
            unique_dihedral_types_2,
        ) = charmm_writer._get_dihedral_types(
            structure_ff,
            use_rb_torsions_2,
            use_dihedrals_2,
            epsilon_conversion_factor,
        )

        unique_dih_typ_unsorted_2 = dict(
            enumerate(
                set(
                    [
                        (
                            round(dihedral.type.c0 * lj_unit, 3),
                            round(dihedral.type.c1 * lj_unit, 3),
                            round(dihedral.type.c2 * lj_unit, 3),
                            round(dihedral.type.c3 * lj_unit, 3),
                            round(dihedral.type.c4 * lj_unit, 3),
                            round(dihedral.type.c5 * lj_unit, 3),
                            round(dihedral.type.scee, 1),
                            round(dihedral.type.scnb, 1),
                            dihedral.atom1.type,
                            dihedral.atom2.type,
                            dihedral.atom3.type,
                            dihedral.atom4.type,
                            dihedral.atom1.residue.name,
                            dihedral.atom2.residue.name,
                            dihedral.atom3.residue.name,
                            dihedral.atom4.residue.name,
                        )
                        for dihedral in structure_ff.rb_torsions
                    ]
                )
            )
        )

        unique_dih_typ_unsorted_2 = OrderedDict(
            [(y, x + 1) for x, y in unique_dih_typ_unsorted_2.items()]
        )

        assert len(unique_dih_typ_unsorted_2) == 7
        assert len(unique_dihedral_types_2) == 5

        # test for error if trying to use impropers in the dihedrals (currently impropers are found but not used in
        # the output)
        # ******** NOTE*************************
        # ******** NOTE*************************
        # These impropers are blank and will need filled in upon adding the improper functionallity.
        # They are kept in the code to identify if there are any impropers in the system and count them
        # ******** NOTE*************************
        # ******** NOTE*************************
        # ******** NOTE*************************
        (
            improper_types_1,
            unique_improper_types_1,
        ) = charmm_writer._get_impropers(
            structure_ff, epsilon_conversion_factor
        )

        assert str(improper_types_1) == "[]"
        assert str(unique_improper_types_1) == "OrderedDict()"

    def test_charmm_angle_reorder(self, ethyl_ether_gomc, methyl_ether_gomc):
        box_reservior_0 = mb.fill_box(
            compound=[ethyl_ether_gomc, methyl_ether_gomc],
            box=[10, 10, 10],
            n_compounds=[10, 10],
        )

        [
            structure_ff,
            coulomb14scalar_dict,
            lj14_scalar_dict,
            residues_applied_list,
        ] = specific_ff_to_residue(
            box_reservior_0,
            forcefield_selection={
                ethyl_ether_gomc.name: "oplsaa",
                methyl_ether_gomc.name: "oplsaa",
            },
            residues=[ethyl_ether_gomc.name, methyl_ether_gomc.name],
            reorder_res_in_pdb_psf=False,
            boxes_for_simulation=1,
        )

        sigma_conversion_factor = 1
        epsilon_conversion_factor = 1
        use_urey_bradleys = False
        angle_types_1, unique_angle_types_1 = charmm_writer._get_angle_types(
            structure_ff,
            sigma_conversion_factor,
            epsilon_conversion_factor,
            use_urey_bradleys,
        )

        # note this sorts all the possible combinations, so this should be the same as the double check (i.e, both 10)
        unique_angle_types_1_unsorted = dict(
            enumerate(
                set(
                    [
                        (
                            round(
                                angle.type.k
                                * (
                                    sigma_conversion_factor**2
                                    / epsilon_conversion_factor
                                ),
                                3,
                            ),
                            round(angle.type.theteq, 3),
                            angle.atom2.type,
                            tuple(sorted((angle.atom1.type, angle.atom3.type))),
                            angle.atom1.residue.name,
                            angle.atom2.residue.name,
                            angle.atom3.residue.name,
                        )
                        for angle in structure_ff.angles
                    ]
                )
            )
        )
        unique_angle_types_1_unsorted = OrderedDict(
            [(y, x + 1) for x, y in unique_angle_types_1_unsorted.items()]
        )

        assert len(unique_angle_types_1_unsorted) == 10
        assert len(unique_angle_types_1) == 10

    def test_bead_atomname_equal_3(self, two_propanol_ua):
        # testing def unique_atom_naming in charmm_writer, expecting when failing
        with pytest.raises(
            ValueError,
            match=r"ERROR: The unique_atom_naming function failed while "
            "running the charmm_writer function. Ensure the proper inputs are "
            "in the bead_to_atom_name_dict.",
        ):
            box_reservior_0 = mb.fill_box(
                compound=[two_propanol_ua], box=[10, 10, 10], n_compounds=[10]
            )

            value_0 = Charmm(
                box_reservior_0,
                "test_bead_atomname_equal_3",
                ff_filename="test_bead_atomname_equal_3",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "Cx", "_HC": "Cxx"},
            )
            value_0.write_inp()
            value_0.write_pdb()
            value_0.write_psf()

    def test_gomc_fix_bonds_angles_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please ensure the residue names in the \({}\) variable "
            r"are in a list.".format("gomc_fix_bonds_angles"),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_bonds_angles="two_propanol_ua.name",
            )

    def test_gomc_fix_bonds_angles_residue_not_in_system(self, two_propanol_ua):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please ensure that all the residue names in the "
            r"{} list are also in the residues list.".format(
                "gomc_fix_bonds_angles"
            ),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_bonds_angles=["WNG"],
            )

    def test_gomc_fix_bonds_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please ensure the residue names in the \({}\) variable "
            r"are in a list.".format("gomc_fix_bonds"),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_bonds="two_propanol_ua.name",
            )

    def test_gomc_fix_bonds_residue_not_in_system(self, two_propanol_ua):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please ensure that all the residue names in the "
            r"{} list are also in the residues list.".format("gomc_fix_bonds"),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_bonds=["WNG"],
            )

    def test_gomc_fix_angles_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please ensure the residue names in the \({}\) variable "
            r"are in a list.".format("gomc_fix_angles"),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_angles="two_propanol_ua.name",
            )

    def test_gomc_fix_angles_residue_not_in_system(self, two_propanol_ua):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please ensure that all the residue names in the "
            r"{} list are also in the residues list.".format("gomc_fix_angles"),
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                gomc_fix_angles=["WNG"],
            )

    def test_fix_residue_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the fix_residue in a list format",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                fix_residue="two_propanol_ua.name",
            )

    def test_fix_residue_string_residue_not_in_system(self, two_propanol_ua):
        with pytest.raises(
            ValueError,
            match=r"Error: Please ensure that all the residue names in the fix_residue "
            r"list are also in the residues list.",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                fix_residue=["WNG"],
            )

    def test_fix_residue_in_box_string(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the fix_residue_in_box in a list format.",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                fix_residue_in_box="two_propanol_ua.name",
            )

    def test_fix_residue_in_box_string_residue_not_in_system(
        self, two_propanol_ua
    ):
        with pytest.raises(
            ValueError,
            match=r"Error: Please ensure that all the residue names in the "
            r"fix_residue_in_box list are also in the residues list.",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
                fix_residue_in_box=["WNG"],
            )

    def test_bead_to_atom_name_dict_list(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the a bead type to atom in the dictionary "
            r"\(bead_to_atom_name_dict\) so GOMC can properly evaluate the "
            r"unique atom names",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict=["_CH3", "C"],
            )

    def test_bead_to_atom_name_dict_not_string_0(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the bead_to_atom_name_dict with only "
            r"string inputs.",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": 0},
            )

    def test_bead_to_atom_name_dict_not_string_1(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the bead_to_atom_name_dict with only "
            r"string inputs.",
        ):
            Charmm(
                two_propanol_ua,
                "charmm_data_UA",
                ff_filename="charmm_data_UA",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={0: "C"},
            )

    def test_1_box_residues_not_all_listed_box_0(
        self, ethane_gomc, ethanol_gomc
    ):
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the residues are not specified, or the residues "
            r"entered does not match the residues that were found and "
            r"built for structure.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename="charmm_data",
                residues=[ethanol_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_2_box_residues_not_all_listed_box_0(
        self, ethane_gomc, ethanol_gomc
    ):
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the residues are not specified, or the residues "
            r"entered does not match the residues that were found and "
            r"built for structure.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethanol_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=["XXX", ethanol_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_2_box_residues_not_all_listed_box_1(
        self, ethane_gomc, ethanol_gomc
    ):
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the residues are not specified, or the residues "
            r"entered does not match the residues that were found and "
            r"built for structure.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethanol_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=["XXX", ethane_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_2_box_residues_listed_2x(self, ethane_gomc, ethanol_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Please enter the residues list \(residues\) that has "
            r"only unique residue names.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethanol_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethanol_gomc.name, ethanol_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_all_residues_are_listed(self, ethane_gomc, ethanol_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the residues are not specified, or the residues "
            r"entered does not match the residues that were found and "
            r"built for structure.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethanol_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethanol_gomc.name],
                forcefield_selection="oplsaa",
            )

    # Test that an empty box (psf and pdb files) can be created to start a simulation
    def test_box_1_empty_test_1(self, two_propanol_ua):
        empty_compound = Box(lengths=[2, 2, 2])

        charmm = Charmm(
            two_propanol_ua,
            "charmm_filled_box",
            structure_box_1=empty_compound,
            filename_box_1="charmm_empty_box",
            ff_filename="charmm_empty_box.inp",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()
        charmm.write_psf()

        with open("charmm_empty_box.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "20.000",
                        "20.000",
                        "20.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]
                    assert out_gomc[i + 1].split() == ["END"]

                else:
                    pass

        assert pdb_read

        with open("charmm_filled_box.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "POL", "A", "1"],
                        ["ATOM", "2", "BD1", "POL", "A", "1"],
                        ["ATOM", "3", "O1", "POL", "A", "1"],
                        ["ATOM", "4", "H1", "POL", "A", "1"],
                        ["ATOM", "5", "C2", "POL", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "EP"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read

    def test_box_1_empty_test_2(self, two_propanol_ua):
        empty_compound = Box(lengths=[3, 3, 3], angles=[90, 90, 90])

        charmm = Charmm(
            two_propanol_ua,
            "charmm_filled_box",
            structure_box_1=empty_compound,
            filename_box_1="charmm_empty_box",
            ff_filename="charmm_empty_box.inp",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()
        charmm.write_psf()

        with open("charmm_empty_box.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "30.000",
                        "30.000",
                        "30.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]
                    assert out_gomc[i + 1].split() == ["END"]

                else:
                    pass

        assert pdb_read

        with open("charmm_filled_box.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "POL", "A", "1"],
                        ["ATOM", "2", "BD1", "POL", "A", "1"],
                        ["ATOM", "3", "O1", "POL", "A", "1"],
                        ["ATOM", "4", "H1", "POL", "A", "1"],
                        ["ATOM", "5", "C2", "POL", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "EP"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read

    def test_box_1_empty_test_3(self, two_propanol_ua):
        empty_compound = Box(lengths=[4, 5, 6])

        test_box_two_propanol_ua_gomc = mb.fill_box(
            compound=[two_propanol_ua], n_compounds=[1], box=[3, 4, 5]
        )

        charmm = Charmm(
            empty_compound,
            "charmm_empty_box",
            structure_box_1=test_box_two_propanol_ua_gomc,
            filename_box_1="charmm_filled_box",
            ff_filename="charmm_empty_box",
            residues=[two_propanol_ua.name],
            forcefield_selection="trappe-ua",
            bead_to_atom_name_dict={"_CH3": "C"},
        )
        charmm.write_pdb()
        charmm.write_psf()

        with open("charmm_empty_box.pdb", "r") as fp:
            pdb_part_1_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_part_1_read = True
                    assert out_gomc[i].split()[0:7] == [
                        "CRYST1",
                        "40.000",
                        "50.000",
                        "60.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]
                    assert out_gomc[i + 1].split() == ["END"]

                else:
                    pass

        assert pdb_part_1_read

        with open("charmm_filled_box.pdb", "r") as fp:
            pdb_part_2_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_part_2_read = True
                    atom_type_res_part_1_list = [
                        ["ATOM", "1", "C1", "POL", "A", "1"],
                        ["ATOM", "2", "BD1", "POL", "A", "1"],
                        ["ATOM", "3", "O1", "POL", "A", "1"],
                        ["ATOM", "4", "H1", "POL", "A", "1"],
                        ["ATOM", "5", "C2", "POL", "A", "1"],
                    ]
                    atom_type_res_part_2_list = [
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "EP"],
                        ["1.00", "0.00", "O"],
                        ["1.00", "0.00", "H"],
                        ["1.00", "0.00", "EP"],
                    ]

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_part_2_read

    def test_box_1_empty_test_4(self):
        empty_compound_box_0 = Box(lengths=[2, 2, 2])
        empty_compound_box_1 = Box(lengths=[3, 3, 3])
        with pytest.raises(
            TypeError,
            match=r"ERROR: Both structure_box_0 and structure_box_0 are empty Boxes {}. "
            "At least 1 structure must be an mbuild compound {} with 1 "
            "or more atoms in it".format(
                type(Box(lengths=[1, 1, 1])), type(Compound())
            ),
        ):
            Charmm(
                empty_compound_box_0,
                "charmm_data_box_0",
                structure_box_1=empty_compound_box_1,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[],
                forcefield_selection="oplsaa",
            )

    def test_box_1_empty_test_5(self):
        empty_compound_box_0 = Box(lengths=[2, 2, 2])
        with pytest.raises(
            TypeError,
            match=r"ERROR: Only 1 structure is provided and it can not be an empty "
            r"mbuild Box {}. "
            "it must be an mbuild compound {} with at least 1 "
            "or more atoms in it.".format(
                type(Box(lengths=[1, 1, 1])), type(Compound())
            ),
        ):
            Charmm(
                empty_compound_box_0,
                "charmm_data_box_0",
                structure_box_1=None,
                filename_box_1=None,
                ff_filename="charmm_data",
                residues=[],
                forcefield_selection="oplsaa",
            )

    def test_box_1_empty_test_6(self, two_propanol_ua):
        with pytest.raises(
            TypeError,
            match=r"ERROR: If you are not providing an empty box, "
            r"you need to specify the atoms/beads as children in the mb.Compound. "
            r"If you are providing and empty box, please do so by specifying and "
            r"mbuild Box \({}\)".format(type(Box(lengths=[1, 1, 1]))),
        ):
            test_box_two_propanol_ua_gomc = mb.fill_box(
                compound=[two_propanol_ua], n_compounds=[1], box=[3, 4, 5]
            )

            empty_compound = mb.Compound()
            Charmm(
                empty_compound,
                "charmm_empty_box",
                structure_box_1=test_box_two_propanol_ua_gomc,
                filename_box_1="charmm_filled_box",
                ff_filename="charmm_empty_box",
                residues=[two_propanol_ua.name],
                forcefield_selection="trappe-ua",
                bead_to_atom_name_dict={"_CH3": "C"},
            )

    def test_structure_box_0_not_mb_compound(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The structure_box_0 expected to be of type: "
            r"{} or {}, received: {}".format(
                type(Compound()),
                type(Box(lengths=[1, 1, 1])),
                type("ethane_gomc"),
            ),
        ):
            Charmm(
                "ethane_gomc",
                "charmm_data_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_structure_box_1_not_mb_compound(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The structure_box_1 expected to be of type: "
            "{} or {}, received: {}".format(
                type(Compound()), type(Box(lengths=[1, 1, 1])), type(0)
            ),
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=0,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name],
                forcefield_selection="oplsaa",
            )

    def test_ff_dict_not_entered(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: Please enter the forcefield_selection as it was not provided.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name],
                forcefield_selection=None,
            )

    def test_mie_non_bonded_type(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Currently the Mie potential \(non_bonded_type\) is not "
            r"supported in this MoSDeF GOMC parameter writer.",
        ):
            charmm = Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name],
                forcefield_selection="oplsaa",
                non_bonded_type="Mie",
            )
            charmm.write_inp()

    def test_other_non_bonded_type(self, ethane_gomc):
        with pytest.raises(
            ValueError,
            match=r"ERROR: Currently this potential \(non_bonded_type\) is not "
            r"supported in this MoSDeF GOMC parameter writer.",
        ):
            charmm = Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name],
                forcefield_selection="oplsaa",
                non_bonded_type="XXX",
            )
            charmm.write_inp()

    def test_diff_1_4_coul_scalars(self, ethane_gomc, two_propanol_ua):
        with pytest.raises(
            ValueError,
            match=r"ERROR: There are multiple 1,4-coulombic scaling factors "
            "GOMC will only accept a singular input for the 1,4-coulombic "
            "scaling factors.",
        ):
            Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=two_propanol_ua,
                filename_box_1="charmm_data_box_1",
                ff_filename="charmm_data",
                residues=[ethane_gomc.name, two_propanol_ua.name],
                forcefield_selection={
                    ethane_gomc.name: "oplsaa",
                    two_propanol_ua.name: "trappe-ua",
                },
            )

    def test_write_inp_wo_ff_filename(self, ethane_gomc):
        with pytest.raises(
            TypeError,
            match=r"ERROR: The force field file name was not specified and in the "
            r"Charmm object. "
            r"Therefore, the force field file \(.inp\) can not be written. "
            r"Please use the force field file name when building the Charmm object, "
            r"then use the write_inp function.",
        ):
            charmm = Charmm(
                ethane_gomc,
                "charmm_data_box_0",
                structure_box_1=ethane_gomc,
                filename_box_1="charmm_data_box_1",
                ff_filename=None,
                forcefield_selection="oplsaa",
                residues=[ethane_gomc.name],
            )
            charmm.write_inp()

    def test_write_inp_with_2_boxes(self, ethane_gomc):
        charmm = Charmm(
            ethane_gomc,
            "charmm_data_box_0",
            structure_box_1=ethane_gomc,
            filename_box_1="charmm_data_box_1",
            ff_filename="charmm_data",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()

        with open("charmm_data.inp", "r") as fp:
            masses_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if (
                    "!atom_types" in line
                    and "mass" in line
                    and "atomTypeForceFieldName_ResidueName" in line
                    and "(i.e., atoms_type_per_utilized_FF)" in line
                ):
                    masses_read = True
                    mass_type_1 = [
                        ["*", "A", "12.010780"],
                        ["*", "B", "1.007947"],
                    ]
                    mass_type_2 = [["opls_135_ETH"], ["opls_140_ETH"]]
                    for j in range(0, len(mass_type_1)):
                        assert (
                            len(out_gomc[i + 1 + j].split("!")[0].split()) == 3
                        )
                        assert (
                            out_gomc[i + 1 + j].split("!")[0].split()[0:3]
                            == mass_type_1[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[4:5] == mass_type_2[j]
                        )

        assert masses_read

    # test cif reader ETA psf writer outputs correct atom and residue numbering using non-orthoganol box
    def test_save_non_othoganol_box_psf(self):
        lattice_cif_ETV_triclinic = load_cif(
            file_or_path=get_fn("ETV_triclinic.cif")
        )
        ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=1, y=1, z=1)
        ETV_triclinic.name = "ETV"

        charmm = Charmm(
            ETV_triclinic,
            "ETV_triclinic",
            ff_filename="ETV_triclinic_FF",
            forcefield_selection={
                ETV_triclinic.name: get_fn(
                    "Charmm_writer_testing_only_zeolite.xml"
                )
            },
            residues=[ETV_triclinic.name],
            bead_to_atom_name_dict=None,
            fix_residue=[ETV_triclinic.name],
        )

        charmm.write_psf()

        with open("ETV_triclinic.psf", "r") as fp:
            psf_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "42 !NATOM" in line:
                    psf_read = True
                    no_O_atoms = 28
                    no_Si_atoms = 14
                    atom_type_charge_etc_list = []
                    for f_i in range(0, no_O_atoms):
                        atom_type_charge_etc_list.append(
                            [
                                str(f_i + 1),
                                "SYS",
                                str(f_i + 1),
                                "ETV",
                                "O1",
                                "A",
                                "-0.400000",
                                "15.9994",
                            ],
                        )
                    for f_i in range(no_O_atoms, no_O_atoms + no_Si_atoms):
                        atom_type_charge_etc_list.append(
                            [
                                str(f_i + 1),
                                "SYS",
                                str(f_i + 1),
                                "ETV",
                                "Si1",
                                "B",
                                "0.800000",
                                "28.0855",
                            ],
                        )

                    for j in range(0, len(atom_type_charge_etc_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:8]
                            == atom_type_charge_etc_list[j]
                        )

                else:
                    pass

        assert psf_read

    # test cif reader ETA pdb writer outputs correct atom and residue numbering using non-orthoganol box
    def test_save_non_othoganol_box_pdb(self):
        lattice_cif_ETV_triclinic = load_cif(
            file_or_path=get_fn("ETV_triclinic.cif")
        )
        ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=1, y=1, z=1)
        ETV_triclinic.name = "ETV"

        charmm = Charmm(
            ETV_triclinic,
            "ETV_triclinic",
            ff_filename="ETV_triclinic_FF",
            forcefield_selection={
                ETV_triclinic.name: get_fn(
                    "Charmm_writer_testing_only_zeolite.xml"
                )
            },
            residues=[ETV_triclinic.name],
            bead_to_atom_name_dict=None,
            fix_residue=[ETV_triclinic.name],
        )

        charmm.write_pdb()

        with open("ETV_triclinic.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):

                if "CRYST1" in line:
                    pdb_read = True
                    crystal_box_length_angles = [
                        "CRYST1",
                        "8.750",
                        "9.648",
                        "10.272",
                        "105.72",
                        "100.19",
                        "97.02",
                    ]

                    no_O_atoms = 28
                    no_Si_atoms = 14
                    atom_type_res_part_1_list = []
                    for f_i in range(0, no_O_atoms):
                        atom_type_res_part_1_list.append(
                            [
                                "ATOM",
                                str(f_i + 1),
                                "O1",
                                "ETV",
                                "A",
                                str(f_i + 1),
                            ]
                        )
                    for f_i in range(no_O_atoms, no_O_atoms + no_Si_atoms):
                        atom_type_res_part_1_list.append(
                            [
                                "ATOM",
                                str(f_i + 1),
                                "Si1",
                                "ETV",
                                "A",
                                str(f_i + 1),
                            ]
                        )

                    atom_type_res_part_2_list = []
                    for f_i in range(0, no_O_atoms):
                        atom_type_res_part_2_list.append(["1.00", "1.00", "O"])
                    for f_i in range(no_O_atoms, no_O_atoms + no_Si_atoms):
                        atom_type_res_part_2_list.append(["1.00", "1.00", "SI"])

                    assert out_gomc[i].split()[0:7] == crystal_box_length_angles

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read

    # test methane UA psf writer outputs correct atom and residue numbering using orthoganol box
    def test_save_othoganol_methane_ua_psf(self):
        methane = mb.Compound(name="MET")
        methane_child_bead = mb.Compound(name="_CH4")
        methane.add(methane_child_bead, inherit_periodicity=False)

        methane_box = mb.fill_box(
            compound=methane, n_compounds=4, box=[1, 1, 1]
        )

        charmm = Charmm(
            methane_box,
            "methane_box",
            ff_filename="methane_box_FF",
            forcefield_selection={methane.name: "trappe-ua"},
            residues=[methane.name],
            bead_to_atom_name_dict={"_CH4": "C"},
        )

        charmm.write_psf()

        with open("methane_box.psf", "r") as fp:
            psf_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "4 !NATOM" in line:
                    psf_read = True
                    no_methane_atoms = 4
                    atom_type_charge_etc_list = []
                    for f_i in range(0, no_methane_atoms):
                        atom_type_charge_etc_list.append(
                            [
                                str(f_i + 1),
                                "SYS",
                                str(f_i + 1),
                                "MET",
                                "C1",
                                "A",
                                "0.000000",
                                "16.0430",
                            ],
                        )

                    for j in range(0, len(atom_type_charge_etc_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:8]
                            == atom_type_charge_etc_list[j]
                        )

                else:
                    pass

        assert psf_read

    # test methane UA pdb writer outputs correct atom and residue numbering using orthoganol box
    def test_save_othoganol_methane_ua_pdb(self):
        methane = mb.Compound(name="MET")
        methane_child_bead = mb.Compound(name="_CH4")
        methane.add(methane_child_bead, inherit_periodicity=False)

        methane_box = mb.fill_box(
            compound=methane, n_compounds=10, box=[1, 2, 3]
        )

        charmm = Charmm(
            methane_box,
            "methane_box",
            ff_filename="methane_box_FF",
            forcefield_selection={methane.name: "trappe-ua"},
            residues=[methane.name],
            bead_to_atom_name_dict={"_CH4": "C"},
        )

        charmm.write_pdb()

        with open("methane_box.pdb", "r") as fp:
            pdb_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if "CRYST1" in line:
                    pdb_read = True
                    crystal_box_length_angles = [
                        "CRYST1",
                        "10.000",
                        "20.000",
                        "30.000",
                        "90.00",
                        "90.00",
                        "90.00",
                    ]

                    no_methane_atoms = 4
                    atom_type_res_part_1_list = []
                    for f_i in range(0, no_methane_atoms):
                        atom_type_res_part_1_list.append(
                            [
                                "ATOM",
                                str(f_i + 1),
                                "C1",
                                "MET",
                                "A",
                                str(f_i + 1),
                            ]
                        )

                    atom_type_res_part_2_list = []
                    for f_i in range(0, no_methane_atoms):
                        atom_type_res_part_2_list.append(["1.00", "0.00", "EP"])

                    assert out_gomc[i].split()[0:7] == crystal_box_length_angles

                    for j in range(0, len(atom_type_res_part_1_list)):
                        assert (
                            out_gomc[i + 1 + j].split()[0:6]
                            == atom_type_res_part_1_list[j]
                        )
                        assert (
                            out_gomc[i + 1 + j].split()[9:12]
                            == atom_type_res_part_2_list[j]
                        )

                else:
                    pass

        assert pdb_read
