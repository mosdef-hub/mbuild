import os

import pytest

import mbuild as mb
import mbuild.formats.gomc_conf_writer as gomc_control
from mbuild.formats.charmm_writer import Charmm
from mbuild.lattice import load_cif
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestGOMCControlFileWriter(BaseTest):
    def test_dict_keys_to_list(
        self,
    ):
        dict = {"a": "1", "b": "2", "c": "3"}
        keys = gomc_control.dict_keys_to_list(dict)

        assert keys == ["a", "b", "c"]

    def test_get_required_data(self):
        value = gomc_control._get_required_data(description=False)
        assert (
            value.sort()
            == [
                "charmm_object",
                "ensemble_type",
                "RunSteps",
                "Temperature",
                "ff_psf_pdb_file_directory",
                "check_input_files_exist",
                "Restart",
                "RestartCheckpoint",
                "Parameters",
                "Coordinates_box_0",
                "override_psf_box_0",
                "Coordinates_box_1",
                "Structure_box_1",
                "binCoordinates_box_0",
                "extendedSystem_box_0",
                "binVelocities_box_0",
                "binCoordinates_box_1",
                "extendedSystem_box_1",
                "binVelocities_box_1",
            ].sort()
        )

        value = gomc_control._get_required_data(description=True)
        assert (
            gomc_control.dict_keys_to_list(value).sort()
            == [
                "charmm_object",
                "ensemble_type",
                "RunSteps",
                "Temperature",
                "ff_psf_pdb_file_directory",
                "Restart",
                "RestartCheckpoint",
                "ExpertMode",
                "check_input_files_exist",
                "Parameters",
                "Coordinate_box_0",
                "Structure_box_0",
                "Coordinate_box_1",
                "Structure_box_1",
                "binCoordinates_box_0",
                "extendedSystem_box_0",
                "binVelocities_box_0",
                "binCoordinates_box_1",
                "extendedSystem_box_1",
                "binVelocities_box_1",
            ].sort()
        )

    def test_get_all_possible_input_variable(self):
        value = gomc_control._get_all_possible_input_variables(
            description=False
        )
        assert (
            value.sort()
            == [
                "PRNG",
                "ParaTypeCHARMM",
                "ParaTypeMie",
                "ParaTypeMARTINI",
                "RcutCoulomb_box_0",
                "RcutCoulomb_box_1",
                "Pressure",
                "Rcut",
                "RcutLow",
                "LRC",
                "IPC",
                "Exclude",
                "Potential",
                "Rswitch",
                "ElectroStatic",
                "Ewald",
                "CachedFourier",
                "Tolerance",
                "Dielectric",
                "PressureCalc",
                "EqSteps",
                "AdjSteps",
                "VDWGeometricSigma",
                "useConstantArea",
                "FixVolBox0",
                "ChemPot",
                "Fugacity",
                "CBMC_First",
                "CBMC_Nth",
                "CBMC_Ang",
                "CBMC_Dih",
                "OutputName",
                "CoordinatesFreq",
                "DCDFreq",
                "RestartFreq",
                "CheckpointFreq",
                "ConsoleFreq",
                "BlockAverageFreq",
                "HistogramFreq",
                "DistName",
                "HistName",
                "RunNumber",
                "RunLetter",
                "SampleFreq",
                "OutEnergy",
                "OutPressure",
                "OutMolNum",
                "OutDensity",
                "OutVolume",
                "OutSurfaceTension",
                "FreeEnergyCalc",
                "MoleculeType",
                "InitialState",
                "LambdaVDW",
                "LambdaCoulomb",
                "ScaleCoulomb",
                "ScalePower",
                "ScaleAlpha",
                "MinSigma",
                "DisFreq",
                "RotFreq",
                "IntraSwapFreq",
                "SwapFreq",
                "RegrowthFreq",
                "CrankShaftFreq",
                "VolFreq",
                "MultiParticleFreq",
                "IntraMEMC-1Freq",
                "MEMC-1Freq",
                "IntraMEMC-2Freq",
                "MEMC-2Freq",
                "IntraMEMC-3Freq",
                "MEMC-3Freq",
                "ExchangeVolumeDim",
                "MEMC_DataInput",
                "TargetedSwap",
            ].sort()
        )

        value = gomc_control._get_all_possible_input_variables(description=True)
        assert (
            gomc_control.dict_keys_to_list(value).sort()
            == [
                "PRNG",
                "ParaTypeCHARMM",
                "ParaTypeMie",
                "ParaTypeMARTINI",
                "RcutCoulomb_box_0",
                "RcutCoulomb_box_1",
                "Pressure",
                "Rcut",
                "RcutLow",
                "LRC",
                "IPC",
                "Exclude",
                "Potential",
                "Rswitch",
                "ElectroStatic",
                "Ewald",
                "CachedFourier",
                "Tolerance",
                "Dielectric",
                "PressureCalc",
                "EqSteps",
                "AdjSteps",
                "VDWGeometricSigma",
                "useConstantArea",
                "FixVolBox0",
                "ChemPot",
                "Fugacity",
                "CBMC_First",
                "CBMC_Nth",
                "CBMC_Ang",
                "CBMC_Dih",
                "OutputName",
                "CoordinatesFreq",
                "DCDFreq",
                "RestartFreq",
                "CheckpointFreq",
                "ConsoleFreq",
                "BlockAverageFreq",
                "HistogramFreq",
                "DistName",
                "HistName",
                "RunNumber",
                "RunLetter",
                "SampleFreq",
                "OutEnergy",
                "OutPressure",
                "OutMolNum",
                "OutDensity",
                "OutVolume",
                "OutSurfaceTension",
                "FreeEnergyCalc",
                "MoleculeType",
                "InitialState",
                "LambdaVDW",
                "LambdaCoulomb",
                "ScaleCoulomb",
                "ScalePower",
                "ScaleAlpha",
                "MinSigma",
                "DisFreq",
                "RotFreq",
                "IntraSwapFreq",
                "SwapFreq",
                "RegrowthFreq",
                "CrankShaftFreq",
                "VolFreq",
                "MultiParticleFreq",
                "IntraMEMC-1Freq",
                "MEMC-1Freq",
                "IntraMEMC-2Freq",
                "MEMC-2Freq",
                "IntraMEMC-3Freq",
                "MEMC-3Freq",
                "ExchangeVolumeDim",
                "MEMC_DataInput",
                "TargetedSwap",
            ].sort()
        )

    def test_get_default_variables_dict(self):
        value = gomc_control._get_default_variables_dict()
        assert (
            gomc_control.dict_keys_to_list(value).sort()
            == [
                "PRNG",
                "ParaTypeCHARMM",
                "ParaTypeMie",
                "ParaTypeMARTINI",
                "RcutCoulomb_box_0",
                "RcutCoulomb_box_1",
                "Pressure",
                "Rcut",
                "RcutLow",
                "LRC",
                "IPC",
                "Exclude",
                "coul_1_4_scaling",
                "Potential",
                "Rswitch",
                "ElectroStatic",
                "Ewald",
                "CachedFourier",
                "Tolerance",
                "Dielectric",
                "PressureCalc",
                "EqSteps",
                "AdjSteps",
                "VDWGeometricSigma",
                "useConstantArea",
                "FixVolBox0",
                "ChemPot",
                "Fugacity",
                "CBMC_First",
                "CBMC_Nth",
                "CBMC_Ang",
                "CBMC_Dih",
                "OutputName",
                "CoordinatesFreq",
                "DCDFreq",
                "RestartFreq",
                "CheckpointFreq",
                "ConsoleFreq",
                "BlockAverageFreq",
                "HistogramFreq",
                "DistName",
                "HistName",
                "RunNumber",
                "RunLetter",
                "SampleFreq",
                "OutEnergy",
                "OutPressure",
                "OutMolNum",
                "OutDensity",
                "OutVolume",
                "OutSurfaceTension",
                "FreeEnergyCalc",
                "MoleculeType",
                "InitialState",
                "LambdaVDW",
                "LambdaCoulomb",
                "ScaleCoulomb",
                "ScalePower",
                "ScaleAlpha",
                "MinSigma",
                "ExchangeVolumeDim",
                "MEMC_DataInput",
                "DisFreq",
                "RotFreq",
                "IntraSwapFreq",
                "SwapFreq",
                "RegrowthFreq",
                "CrankShaftFreq",
                "VolFreq",
                "MultiParticleFreq",
                "IntraMEMC-1Freq",
                "MEMC-1Freq",
                "IntraMEMC-2Freq",
                "MEMC-2Freq",
                "IntraMEMC-3Freq",
                "MEMC-3Freq",
                "TargetedSwap",
            ].sort()
        )

    def test_print_ensemble_info(self):

        try:
            gomc_control.print_required_input(description=True)
            gomc_control.print_required_input(description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "NVT", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "NVT", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "NPT", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "NPT", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "GEMC_NVT", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "GEMC_NVT", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "GEMC_NPT", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "GEMC_NPT", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "GCMC", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "GCMC", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables(
                "XXXXX", description=True
            )
            gomc_control.print_valid_ensemble_input_variables(
                "XXXXX", description=False
            )
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "FAILED"

    def test_get_possible_ensemble_input_variables(self):
        with pytest.warns(
            UserWarning,
            match="WARNING: The ensemble_type selected for "
            "the _get_possible_ensemble_input_variables "
            "function is not valid.",
        ):
            gomc_control._get_possible_ensemble_input_variables("XXX")

    def test_wrong_ensemble_gomccontrol(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        ensemble_input = "XXXXX"
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of '{}' is not a valid ensemble option. "
            r"Please choose the 'NPT', 'NVT', 'GEMC_NVT', 'GEMC_NPT', or 'GCMC' "
            "ensembles".format(ensemble_input),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_wrong_ensemble_gomccontrol",
                ensemble_input,
                100,
                300,
                check_input_files_exist=False,
            )

    def test_charmm_ff_name_is_none(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename=None,
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"The force field file name was not specified and in the Charmm object \({}\)."
            r"Therefore, the force field file \(.inp\) can not be written, and thus, the "
            r"GOMC control file \(.conf\) can not be created. Please use the force field file "
            r"name when building the Charmm object".format(type(Charmm)),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_charmm_ff_name_is_none",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
            )

    def test_input_variables_dict_wrong_value(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The input_variables_dict variable is not None or a dictionary.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_input_variables_dict_wrong_value",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                input_variables_dict="XXXXX",
            )

    def test_not_entered_charmm_object(self):
        not_charmm_object = "XXXXX"
        with pytest.raises(
            TypeError,
            match=r"ERROR: The variable supplied is a \({}\), not a charmm_object \({}\)"
            r"".format(type(not_charmm_object), type(Charmm)),
        ):
            gomc_control.write_gomc_control_file(
                not_charmm_object,
                "test_not_charmm_object",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
            )

    def test_save_basic_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane",
            ff_filename="ethane",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()
        charmm.write_psf()
        charmm.write_pdb()

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_NVT.conf",
            "NVT",
            10,
            300,
            check_input_files_exist=True,
            Restart=False,
        )

        with open("test_save_basic_NVT.conf", "r") as fp:
            variables_read_dict = {
                "Restart": False,
                "ExpertMode": False,
                "PRNG": False,
                "ParaTypeCHARMM": False,
                "Parameters": False,
                "Coordinates": False,
                "Structure": False,
                "Temperature": False,
                "Potential": False,
                "LRC": False,
                "IPC": False,
                "Rcut": False,
                "RcutLow": False,
                "VDWGeometricSigma": False,
                "Exclude": False,
                "Ewald": False,
                "ElectroStatic": False,
                "CachedFourier": False,
                "Tolerance": False,
                "1-4scaling": False,
                "PressureCalc": False,
                "RunSteps": False,
                "EqSteps": False,
                "AdjSteps": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
                "CellBasisVector1": False,
                "CellBasisVector2": False,
                "CellBasisVector3": False,
                "CBMC_First": False,
                "CBMC_Nth": False,
                "CBMC_Ang": False,
                "CBMC_Dih": False,
                "OutputName": False,
                "RestartFreq": False,
                "CheckpointFreq": False,
                "CoordinatesFreq": False,
                "ConsoleFreq": False,
                "BlockAverageFreq": False,
                "HistogramFreq": False,
                "DistName": False,
                "HistName": False,
                "RunNumber": False,
                "RunLetter": False,
                "SampleFreq": False,
                "OutEnergy": False,
                "OutPressure": False,
                "OutMolNum": False,
                "OutDensity": False,
                "OutVolume": False,
                "OutSurfaceTension": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("PRNG "):
                    variables_read_dict["PRNG"] = True
                    split_line = line.split()
                    assert split_line[1] == "RANDOM"

                elif line.startswith("ParaTypeCHARMM "):
                    variables_read_dict["ParaTypeCHARMM"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "ethane.inp"

                elif line.startswith("Coordinates "):
                    variables_read_dict["Coordinates"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane.pdb"

                elif line.startswith("Structure "):
                    variables_read_dict["Structure"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane.psf"

                elif line.startswith("Temperature "):
                    variables_read_dict["Temperature"] = True
                    split_line = line.split()
                    assert split_line[1] == "300"

                elif line.startswith("Potential "):
                    variables_read_dict["Potential"] = True
                    split_line = line.split()
                    assert split_line[1] == "VDW"

                elif line.startswith("LRC "):
                    variables_read_dict["LRC"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("IPC "):
                    variables_read_dict["IPC"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("Rcut "):
                    variables_read_dict["Rcut"] = True
                    split_line = line.split()
                    assert split_line[1] == "10"

                elif line.startswith("RcutLow "):
                    variables_read_dict["RcutLow"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("Exclude "):
                    variables_read_dict["Exclude"] = True
                    split_line = line.split()
                    assert split_line[1] == "1-3"

                elif line.startswith("Ewald "):
                    variables_read_dict["Ewald"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("ElectroStatic "):
                    variables_read_dict["ElectroStatic"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("CachedFourier "):
                    variables_read_dict["CachedFourier"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("Tolerance "):
                    variables_read_dict["Tolerance"] = True
                    split_line = line.split()
                    assert split_line[1] == "1e-05"

                elif line.startswith("1-4scaling "):
                    variables_read_dict["1-4scaling"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.5"

                elif line.startswith("PressureCalc "):
                    variables_read_dict["PressureCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("RunSteps "):
                    variables_read_dict["RunSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "10"

                elif line.startswith("EqSteps "):
                    variables_read_dict["EqSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"

                elif line.startswith("AdjSteps "):
                    variables_read_dict["AdjSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.3"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.3"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("CellBasisVector1 "):
                    variables_read_dict["CellBasisVector1"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "10.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 "):
                    variables_read_dict["CellBasisVector2"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "10.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 "):
                    variables_read_dict["CellBasisVector3"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "10.0"

                elif line.startswith("CBMC_First "):
                    variables_read_dict["CBMC_First"] = True
                    split_line = line.split()
                    assert split_line[1] == "12"

                elif line.startswith("CBMC_Nth"):
                    variables_read_dict["CBMC_Nth"] = True
                    split_line = line.split()
                    assert split_line[1] == "10"

                elif line.startswith("CBMC_Ang "):
                    variables_read_dict["CBMC_Ang"] = True
                    split_line = line.split()
                    assert split_line[1] == "50"

                elif line.startswith("CBMC_Dih "):
                    variables_read_dict["CBMC_Dih"] = True
                    split_line = line.split()
                    assert split_line[1] == "50"

                elif line.startswith("OutputName "):
                    variables_read_dict["OutputName"] = True
                    split_line = line.split()
                    assert split_line[1] == "Output_data"

                elif line.startswith("RestartFreq "):
                    variables_read_dict["RestartFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("CheckpointFreq "):
                    variables_read_dict["CheckpointFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("CoordinatesFreq "):
                    variables_read_dict["CoordinatesFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("ConsoleFreq "):
                    variables_read_dict["ConsoleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("BlockAverageFreq "):
                    variables_read_dict["BlockAverageFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("HistogramFreq "):
                    variables_read_dict["HistogramFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("DistName "):
                    variables_read_dict["DistName"] = True
                    split_line = line.split()
                    assert split_line[1] == "dis"

                elif line.startswith("HistName "):
                    variables_read_dict["HistName"] = True
                    split_line = line.split()
                    assert split_line[1] == "his"

                elif line.startswith("RunNumber "):
                    variables_read_dict["RunNumber"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"

                elif line.startswith("RunLetter "):
                    variables_read_dict["RunLetter"] = True
                    split_line = line.split()
                    assert split_line[1] == "a"

                elif line.startswith("SampleFreq "):
                    variables_read_dict["SampleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"

                elif line.startswith("OutEnergy "):
                    variables_read_dict["OutEnergy"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                elif line.startswith("OutPressure "):
                    variables_read_dict["OutPressure"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                elif line.startswith("OutMolNum "):
                    variables_read_dict["OutMolNum"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                elif line.startswith("OutDensity "):
                    variables_read_dict["OutDensity"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                elif line.startswith("OutVolume "):
                    variables_read_dict["OutVolume"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                elif line.startswith("OutSurfaceTension "):
                    variables_read_dict["OutSurfaceTension"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                else:
                    pass

        assert variables_read_dict == {
            "Restart": True,
            "ExpertMode": True,
            "PRNG": True,
            "ParaTypeCHARMM": True,
            "Parameters": True,
            "Coordinates": True,
            "Structure": True,
            "Temperature": True,
            "Potential": True,
            "LRC": True,
            "IPC": True,
            "Rcut": True,
            "RcutLow": True,
            "VDWGeometricSigma": True,
            "Exclude": True,
            "Ewald": True,
            "ElectroStatic": True,
            "CachedFourier": True,
            "Tolerance": True,
            "1-4scaling": True,
            "PressureCalc": True,
            "RunSteps": True,
            "EqSteps": True,
            "AdjSteps": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
            "CellBasisVector1": True,
            "CellBasisVector2": True,
            "CellBasisVector3": True,
            "CBMC_First": True,
            "CBMC_Nth": True,
            "CBMC_Ang": True,
            "CBMC_Dih": True,
            "OutputName": True,
            "RestartFreq": True,
            "CheckpointFreq": True,
            "CoordinatesFreq": True,
            "ConsoleFreq": True,
            "BlockAverageFreq": True,
            "HistogramFreq": True,
            "DistName": True,
            "HistName": True,
            "RunNumber": True,
            "RunLetter": True,
            "SampleFreq": True,
            "OutEnergy": True,
            "OutPressure": True,
            "OutMolNum": True,
            "OutDensity": True,
            "OutVolume": True,
            "OutSurfaceTension": True,
        }

    def test_save_basic_NPT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane",
            ff_filename="ethane",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()
        charmm.write_psf()
        charmm.write_pdb()

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_NPT.conf",
            "NPT",
            1000,
            500,
            check_input_files_exist=True,
            input_variables_dict={"Potential": "SWITCH"},
        )

        with open("test_save_basic_NPT.conf", "r") as fp:
            variables_read_dict = {
                "Potential": False,
                "Rswitch": False,
                "Rcut": False,
                "Pressure": False,
                "Temperature": False,
                "PressureCalc": False,
                "RunSteps": False,
                "EqSteps": False,
                "AdjSteps": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
                "CellBasisVector1": False,
                "CellBasisVector2": False,
                "CellBasisVector3": False,
                "RestartFreq": False,
                "CheckpointFreq": False,
                "CoordinatesFreq": False,
                "ConsoleFreq": False,
                "BlockAverageFreq": False,
                "HistogramFreq": False,
                "SampleFreq": False,
                "VDWGeometricSigma": False,
                "useConstantArea": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Potential "):
                    variables_read_dict["Potential"] = True
                    split_line = line.split()
                    assert split_line[1] == "SWITCH"

                elif line.startswith("Rswitch "):
                    variables_read_dict["Rswitch"] = True
                    split_line = line.split()
                    assert split_line[1] == "9"

                elif line.startswith("Rcut "):
                    variables_read_dict["Rcut"] = True
                    split_line = line.split()
                    assert split_line[1] == "10"

                elif line.startswith("Pressure "):
                    variables_read_dict["Pressure"] = True
                    split_line = line.split()
                    assert split_line[1] == "1.01325"

                elif line.startswith("Temperature "):
                    variables_read_dict["Temperature"] = True
                    split_line = line.split()
                    assert split_line[1] == "500"

                elif line.startswith("PressureCalc "):
                    variables_read_dict["PressureCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("RunSteps "):
                    variables_read_dict["RunSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "1000"

                elif line.startswith("EqSteps "):
                    variables_read_dict["EqSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "100"

                elif line.startswith("AdjSteps "):
                    variables_read_dict["AdjSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "100"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.29"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.3"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("CellBasisVector1 "):
                    variables_read_dict["CellBasisVector1"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "20.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 "):
                    variables_read_dict["CellBasisVector2"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "20.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 "):
                    variables_read_dict["CellBasisVector3"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "20.0"

                elif line.startswith("RestartFreq "):
                    variables_read_dict["RestartFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("CheckpointFreq "):
                    variables_read_dict["CheckpointFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("CoordinatesFreq "):
                    variables_read_dict["CoordinatesFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("ConsoleFreq "):
                    variables_read_dict["ConsoleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("BlockAverageFreq "):
                    variables_read_dict["BlockAverageFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("HistogramFreq "):
                    variables_read_dict["HistogramFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "100"

                elif line.startswith("SampleFreq "):
                    variables_read_dict["SampleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "100"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("useConstantArea "):
                    variables_read_dict["useConstantArea"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                else:
                    pass

        assert variables_read_dict == {
            "Potential": True,
            "Rswitch": True,
            "Rcut": True,
            "Pressure": True,
            "Temperature": True,
            "PressureCalc": True,
            "RunSteps": True,
            "EqSteps": True,
            "AdjSteps": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
            "CellBasisVector1": True,
            "CellBasisVector2": True,
            "CellBasisVector3": True,
            "RestartFreq": True,
            "CheckpointFreq": True,
            "CoordinatesFreq": True,
            "ConsoleFreq": True,
            "BlockAverageFreq": True,
            "HistogramFreq": True,
            "SampleFreq": True,
            "VDWGeometricSigma": True,
            "useConstantArea": True,
        }

    def test_save_basic_GCMC(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()
        charmm.write_psf()
        charmm.write_pdb()

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GCMC.conf",
            "GCMC",
            100000,
            500,
            check_input_files_exist=True,
            input_variables_dict={
                "ChemPot": {"ETH": -4000},
                "VDWGeometricSigma": True,
            },
        )

        with open("test_save_basic_GCMC.conf", "r") as fp:
            variables_read_dict = {
                "Parameters": False,
                "Coordinates 0": False,
                "Coordinates 1": False,
                "Structure 0": False,
                "Structure 1": False,
                "Temperature": False,
                "ChemPot": False,
                "PressureCalc": False,
                "RunSteps": False,
                "EqSteps": False,
                "AdjSteps": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "SwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
                "CellBasisVector1 0": False,
                "CellBasisVector2 0": False,
                "CellBasisVector3 0": False,
                "CellBasisVector1 1": False,
                "CellBasisVector2 1": False,
                "CellBasisVector3 1": False,
                "RestartFreq": False,
                "CheckpointFreq": False,
                "CoordinatesFreq": False,
                "ConsoleFreq": False,
                "BlockAverageFreq": False,
                "HistogramFreq": False,
                "SampleFreq": False,
                "VDWGeometricSigma": False,
                "ExpertMode": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "ethane_FF.inp"

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_box_0.pdb"

                elif line.startswith("Coordinates 1"):
                    variables_read_dict["Coordinates 1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "ethane_box_1.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_box_0.psf"

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure 1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "ethane_box_1.psf"

                elif line.startswith("Temperature "):
                    variables_read_dict["Temperature"] = True
                    split_line = line.split()
                    assert split_line[1] == "500"

                elif line.startswith("ChemPot "):
                    variables_read_dict["ChemPot"] = True
                    split_line = line.split()
                    assert split_line[1] == "ETH"
                    assert split_line[2] == "-4000"

                elif line.startswith("PressureCalc "):
                    variables_read_dict["PressureCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("RunSteps "):
                    variables_read_dict["RunSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "100000"

                elif line.startswith("EqSteps "):
                    variables_read_dict["EqSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "10000"

                elif line.startswith("AdjSteps "):
                    variables_read_dict["AdjSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "1000"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("SwapFreq "):
                    variables_read_dict["SwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.35"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("CellBasisVector1 0"):
                    variables_read_dict["CellBasisVector1 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "20.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 0"):
                    variables_read_dict["CellBasisVector2 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "20.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 0"):
                    variables_read_dict["CellBasisVector3 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "20.0"

                elif line.startswith("CellBasisVector1 1"):
                    variables_read_dict["CellBasisVector1 1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "20.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 1"):
                    variables_read_dict["CellBasisVector2 1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "20.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 1"):
                    variables_read_dict["CellBasisVector3 1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "20.0"

                elif line.startswith("RestartFreq "):
                    variables_read_dict["RestartFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("CheckpointFreq "):
                    variables_read_dict["CheckpointFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("CoordinatesFreq "):
                    variables_read_dict["CoordinatesFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("ConsoleFreq "):
                    variables_read_dict["ConsoleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("BlockAverageFreq "):
                    variables_read_dict["BlockAverageFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("HistogramFreq "):
                    variables_read_dict["HistogramFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "10000"

                elif line.startswith("SampleFreq "):
                    variables_read_dict["SampleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "500"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                else:
                    pass

        assert variables_read_dict == {
            "Parameters": True,
            "Coordinates 0": True,
            "Coordinates 1": True,
            "Structure 0": True,
            "Structure 1": True,
            "Temperature": True,
            "ChemPot": True,
            "PressureCalc": True,
            "RunSteps": True,
            "EqSteps": True,
            "AdjSteps": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "SwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
            "CellBasisVector1 0": True,
            "CellBasisVector2 0": True,
            "CellBasisVector3 0": True,
            "CellBasisVector1 1": True,
            "CellBasisVector2 1": True,
            "CellBasisVector3 1": True,
            "RestartFreq": True,
            "CheckpointFreq": True,
            "CoordinatesFreq": True,
            "ConsoleFreq": True,
            "BlockAverageFreq": True,
            "HistogramFreq": True,
            "SampleFreq": True,
            "VDWGeometricSigma": True,
            "ExpertMode": True,
        }

    def test_save_basic_GEMC_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()
        charmm.write_psf()
        charmm.write_pdb()

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GEMC_NVT.conf",
            "GEMC_NVT",
            1000000,
            500,
            check_input_files_exist=True,
        )

        with open("test_save_basic_GEMC_NVT.conf", "r") as fp:
            variables_read_dict = {
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "SwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("SwapFreq "):
                    variables_read_dict["SwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                else:
                    pass

        assert variables_read_dict == {
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "SwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
        }

    def test_save_basic_GEMC_NPT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GEMC_NPT.conf",
            "GEMC-NPT",
            1000000,
            500,
            check_input_files_exist=False,
            input_variables_dict={
                "Pressure": 10,
                "useConstantArea": True,
                "FixVolBox0": True,
                "RcutCoulomb_box_0": 14,
                "RcutCoulomb_box_1": 14,
            },
        )

        with open("test_save_basic_GEMC_NPT.conf", "r") as fp:
            variables_read_dict = {
                "Pressure": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "SwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
                "VolFreq": False,
                "useConstantArea": False,
                "FixVolBox0": False,
                "RcutCoulomb_box_0": False,
                "RcutCoulomb_box_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Pressure "):
                    variables_read_dict["Pressure"] = True
                    split_line = line.split()
                    assert split_line[1] == "10"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.19"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("SwapFreq "):
                    variables_read_dict["SwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("VolFreq "):
                    variables_read_dict["VolFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.01"

                elif line.startswith("useConstantArea "):
                    variables_read_dict["useConstantArea"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("FixVolBox0 "):
                    variables_read_dict["FixVolBox0"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("RcutCoulomb 0 "):
                    variables_read_dict["RcutCoulomb_box_0"] = True
                    split_line = line.split()
                    assert split_line[2] == "14"

                elif line.startswith("RcutCoulomb 1 "):
                    variables_read_dict["RcutCoulomb_box_1"] = True
                    split_line = line.split()
                    assert split_line[2] == "14"

                else:
                    pass

        assert variables_read_dict == {
            "Pressure": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "SwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
            "VolFreq": True,
            "useConstantArea": True,
            "FixVolBox0": True,
            "RcutCoulomb_box_0": True,
            "RcutCoulomb_box_1": True,
        }

    def test_save_change_most_variable_NVT(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[4.0, 4.0, 4.0],
        )
        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_change_most_variable_NVT.conf",
            "NVT",
            100000,
            300,
            check_input_files_exist=False,
            Restart=False,
            input_variables_dict={
                "PRNG": 123,
                "ParaTypeCHARMM": True,
                "ParaTypeMARTINI": False,
                "ParaTypeMie": False,
                "LRC": False,
                "IPC": True,
                "Rcut": 12,
                "RcutLow": 8,
                "Exclude": "1-4",
                "Ewald": False,
                "ElectroStatic": False,
                "CachedFourier": True,
                "RcutCoulomb_box_0": 14,
                "PressureCalc": [False, 4],
                "Tolerance": 0.01,
                "DisFreq": 0.2,
                "RotFreq": 0.2,
                "IntraSwapFreq": 0.1,
                "RegrowthFreq": 0.1,
                "CrankShaftFreq": 0.2,
                "MultiParticleFreq": 0.05,
                "IntraMEMC-1Freq": 0.05,
                "IntraMEMC-2Freq": 0.05,
                "IntraMEMC-3Freq": 0.05,
                "TargetedSwapFreq": 0.00,
                "CBMC_First": 55,
                "CBMC_Nth": 66,
                "CBMC_Ang": 33,
                "CBMC_Dih": 22,
                "OutputName": "test_out",
                "RestartFreq": [False, 50],
                "CheckpointFreq": [False, 50],
                "CoordinatesFreq": [False, 50],
                "ConsoleFreq": [False, 500],
                "BlockAverageFreq": [False, 50],
                "HistogramFreq": [False, 50],
                "DistName": "dist",
                "HistName": "hist",
                "RunNumber": 4,
                "RunLetter": "c",
                "SampleFreq": 25,
                "FreeEnergyCalc": [True, 50],
                "MoleculeType": ["ETH", 1],
                "InitialState": 3,
                "LambdaVDW": [0.0, 0.1, 0.1, 0.2, 1.0],
                "LambdaCoulomb": [0.0, 0.0, 0.3, 0.8, 1.0],
                "MEMC_DataInput": [
                    [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                ],
                "OutEnergy": [False, False],
                "OutPressure": [False, False],
                "OutMolNum": [False, False],
                "OutDensity": [False, False],
                "OutVolume": [False, False],
                "OutSurfaceTension": [True, True],
            },
        )

        with open("test_save_change_most_variable_NVT.conf", "r") as fp:
            variables_read_dict = {
                "Restart": False,
                "ExpertMode": False,
                "PRNG": False,
                "Random_Seed": False,
                "ParaTypeCHARMM": False,
                "Parameters": False,
                "Coordinates": False,
                "Structure": False,
                "Temperature": False,
                "Potential": False,
                "LRC": False,
                "IPC": False,
                "Rcut": False,
                "RcutLow": False,
                "Exclude": False,
                "Ewald": False,
                "ElectroStatic": False,
                "CachedFourier": False,
                "Tolerance": False,
                "1-4scaling": False,
                "RcutCoulomb 0": False,
                "PressureCalc": False,
                "RunSteps": False,
                "EqSteps": False,
                "AdjSteps": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
                "CrankShaftFreq": False,
                "MultiParticleFreq": False,
                "IntraMEMC-1Freq": False,
                "IntraMEMC-2Freq": False,
                "IntraMEMC-3Freq": False,
                "CellBasisVector1": False,
                "CellBasisVector2": False,
                "CellBasisVector3": False,
                "FreeEnergyCalc": False,
                "MoleculeType": False,
                "InitialState": False,
                "ScalePower": False,
                "ScaleAlpha": False,
                "MinSigma": False,
                "ScaleCoulomb": False,
                "# States": False,
                "LambdaVDW": False,
                "LambdaCoulomb": False,
                "CBMC_First": False,
                "CBMC_Nth": False,
                "CBMC_Ang": False,
                "CBMC_Dih": False,
                "OutputName": False,
                "RestartFreq": False,
                "CheckpointFreq": False,
                "CoordinatesFreq": False,
                "ConsoleFreq": False,
                "BlockAverageFreq": False,
                "HistogramFreq": False,
                "DistName": False,
                "HistName": False,
                "RunNumber": False,
                "RunLetter": False,
                "SampleFreq": False,
                "OutEnergy": False,
                "OutPressure": False,
                "OutMolNum": False,
                "OutDensity": False,
                "OutVolume": False,
                "OutSurfaceTension": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("PRNG "):
                    variables_read_dict["PRNG"] = True
                    split_line = line.split()
                    assert split_line[1] == "INTSEED"

                elif line.startswith("Random_Seed "):
                    variables_read_dict["Random_Seed"] = True
                    split_line = line.split()
                    assert split_line[1] == "123"

                elif line.startswith("ParaTypeCHARMM "):
                    variables_read_dict["ParaTypeCHARMM"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "ethane_ethanol.inp"

                elif line.startswith("Coordinates "):
                    variables_read_dict["Coordinates"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_ethanol.pdb"

                elif line.startswith("Structure "):
                    variables_read_dict["Structure"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_ethanol.psf"

                elif line.startswith("Temperature "):
                    variables_read_dict["Temperature"] = True
                    split_line = line.split()
                    assert split_line[1] == "300"

                elif line.startswith("Potential "):
                    variables_read_dict["Potential"] = True
                    split_line = line.split()
                    assert split_line[1] == "VDW"

                elif line.startswith("LRC "):
                    variables_read_dict["LRC"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("IPC "):
                    variables_read_dict["IPC"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("Rcut "):
                    variables_read_dict["Rcut"] = True
                    split_line = line.split()
                    assert split_line[1] == "12"

                elif line.startswith("RcutLow "):
                    variables_read_dict["RcutLow"] = True
                    split_line = line.split()
                    assert split_line[1] == "8"

                elif line.startswith("Exclude "):
                    variables_read_dict["Exclude"] = True
                    split_line = line.split()
                    assert split_line[1] == "1-4"

                elif line.startswith("Ewald "):
                    variables_read_dict["Ewald"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("ElectroStatic "):
                    variables_read_dict["ElectroStatic"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("CachedFourier "):
                    variables_read_dict["CachedFourier"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("Tolerance "):
                    variables_read_dict["Tolerance"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.01"

                elif line.startswith("1-4scaling "):
                    variables_read_dict["1-4scaling"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.5"

                elif line.startswith("RcutCoulomb 0 "):
                    variables_read_dict["RcutCoulomb 0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "14"

                elif line.startswith("PressureCalc "):
                    variables_read_dict["PressureCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "4"

                elif line.startswith("RunSteps "):
                    variables_read_dict["RunSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "100000"

                elif line.startswith("EqSteps "):
                    variables_read_dict["EqSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "10000"

                elif line.startswith("AdjSteps "):
                    variables_read_dict["AdjSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "1000"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("CrankShaftFreq "):
                    variables_read_dict["CrankShaftFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.2"

                elif line.startswith("MultiParticleFreq "):
                    variables_read_dict["MultiParticleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.05"

                elif line.startswith("IntraMEMC-1Freq "):
                    variables_read_dict["IntraMEMC-1Freq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.05"

                elif line.startswith("IntraMEMC-2Freq "):
                    variables_read_dict["IntraMEMC-2Freq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.05"

                elif line.startswith("IntraMEMC-3Freq "):
                    variables_read_dict["IntraMEMC-3Freq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.05"

                elif line.startswith("CellBasisVector1 "):
                    variables_read_dict["CellBasisVector1"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "40.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 "):
                    variables_read_dict["CellBasisVector2"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "40.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 "):
                    variables_read_dict["CellBasisVector3"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "40.0"

                elif line.startswith("FreeEnergyCalc "):
                    variables_read_dict["FreeEnergyCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "50"

                elif line.startswith("MoleculeType "):
                    variables_read_dict["MoleculeType"] = True
                    split_line = line.split()
                    assert split_line[1] == "ETH"
                    assert split_line[2] == "1"

                elif line.startswith("InitialState "):
                    variables_read_dict["InitialState"] = True
                    split_line = line.split()
                    assert split_line[1] == "3"

                elif line.startswith("ScalePower "):
                    variables_read_dict["ScalePower"] = True
                    split_line = line.split()
                    assert split_line[1] == "2"

                elif line.startswith("ScaleAlpha "):
                    variables_read_dict["ScaleAlpha"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.5"

                elif line.startswith("MinSigma "):
                    variables_read_dict["MinSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "3"

                elif line.startswith("ScaleCoulomb "):
                    variables_read_dict["ScaleCoulomb"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"

                elif line.startswith("# States "):
                    variables_read_dict["# States"] = True
                    split_line = line.split()
                    assert split_line[2] == "0"
                    assert split_line[3] == "1"
                    assert split_line[4] == "2"
                    assert split_line[5] == "3"

                elif line.startswith("LambdaVDW "):
                    variables_read_dict["LambdaVDW"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.0"
                    assert split_line[2] == "0.1"
                    assert split_line[3] == "0.1"
                    assert split_line[4] == "0.2"
                    assert split_line[5] == "1.0"

                elif line.startswith("LambdaCoulomb "):
                    variables_read_dict["LambdaCoulomb"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.0"
                    assert split_line[2] == "0.0"
                    assert split_line[3] == "0.3"
                    assert split_line[4] == "0.8"
                    assert split_line[5] == "1.0"

                elif line.startswith("CBMC_First "):
                    variables_read_dict["CBMC_First"] = True
                    split_line = line.split()
                    assert split_line[1] == "55"

                elif line.startswith("CBMC_Nth "):
                    variables_read_dict["CBMC_Nth"] = True
                    split_line = line.split()
                    assert split_line[1] == "66"

                elif line.startswith("CBMC_Ang "):
                    variables_read_dict["CBMC_Ang"] = True
                    split_line = line.split()
                    assert split_line[1] == "33"

                elif line.startswith("CBMC_Dih "):
                    variables_read_dict["CBMC_Dih"] = True
                    split_line = line.split()
                    assert split_line[1] == "22"

                elif line.startswith("OutputName "):
                    variables_read_dict["OutputName"] = True
                    split_line = line.split()
                    assert split_line[1] == "test_out"

                elif line.startswith("RestartFreq "):
                    variables_read_dict["RestartFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "50"

                elif line.startswith("CheckpointFreq "):
                    variables_read_dict["CheckpointFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "50"

                elif line.startswith("CoordinatesFreq "):
                    variables_read_dict["CoordinatesFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "50"

                elif line.startswith("ConsoleFreq "):
                    variables_read_dict["ConsoleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "500"

                elif line.startswith("BlockAverageFreq "):
                    variables_read_dict["BlockAverageFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "50"

                elif line.startswith("HistogramFreq "):
                    variables_read_dict["HistogramFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "50"

                elif line.startswith("DistName "):
                    variables_read_dict["DistName"] = True
                    split_line = line.split()
                    assert split_line[1] == "dist"

                elif line.startswith("HistName "):
                    variables_read_dict["HistName"] = True
                    split_line = line.split()
                    assert split_line[1] == "hist"

                elif line.startswith("RunNumber "):
                    variables_read_dict["RunNumber"] = True
                    split_line = line.split()
                    assert split_line[1] == "4"

                elif line.startswith("RunLetter "):
                    variables_read_dict["RunLetter"] = True
                    split_line = line.split()
                    assert split_line[1] == "c"

                elif line.startswith("SampleFreq "):
                    variables_read_dict["SampleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "25"

                elif line.startswith("OutEnergy "):
                    variables_read_dict["OutEnergy"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                elif line.startswith("OutPressure "):
                    variables_read_dict["OutPressure"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                elif line.startswith("OutMolNum "):
                    variables_read_dict["OutMolNum"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                elif line.startswith("OutDensity "):
                    variables_read_dict["OutDensity"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                elif line.startswith("OutVolume "):
                    variables_read_dict["OutVolume"] = True
                    split_line = line.split()
                    assert split_line[1] == "False"
                    assert split_line[2] == "False"

                elif line.startswith("OutSurfaceTension "):
                    variables_read_dict["OutSurfaceTension"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "True"

                else:
                    pass

        assert variables_read_dict == {
            "Restart": True,
            "ExpertMode": True,
            "PRNG": True,
            "Random_Seed": True,
            "ParaTypeCHARMM": True,
            "Parameters": True,
            "Coordinates": True,
            "Structure": True,
            "Temperature": True,
            "Potential": True,
            "LRC": True,
            "IPC": True,
            "Rcut": True,
            "RcutLow": True,
            "Exclude": True,
            "Ewald": True,
            "ElectroStatic": True,
            "CachedFourier": True,
            "Tolerance": True,
            "1-4scaling": True,
            "RcutCoulomb 0": True,
            "PressureCalc": True,
            "RunSteps": True,
            "EqSteps": True,
            "AdjSteps": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
            "CrankShaftFreq": True,
            "MultiParticleFreq": True,
            "IntraMEMC-1Freq": True,
            "IntraMEMC-2Freq": True,
            "IntraMEMC-3Freq": True,
            "CellBasisVector1": True,
            "CellBasisVector2": True,
            "CellBasisVector3": True,
            "FreeEnergyCalc": True,
            "MoleculeType": True,
            "InitialState": True,
            "ScalePower": True,
            "ScaleAlpha": True,
            "MinSigma": True,
            "ScaleCoulomb": True,
            "# States": True,
            "LambdaVDW": True,
            "LambdaCoulomb": True,
            "CBMC_First": True,
            "CBMC_Nth": True,
            "CBMC_Ang": True,
            "CBMC_Dih": True,
            "OutputName": True,
            "RestartFreq": True,
            "CheckpointFreq": True,
            "CoordinatesFreq": True,
            "ConsoleFreq": True,
            "BlockAverageFreq": True,
            "HistogramFreq": True,
            "DistName": True,
            "HistName": True,
            "RunNumber": True,
            "RunLetter": True,
            "SampleFreq": True,
            "OutEnergy": True,
            "OutPressure": True,
            "OutMolNum": True,
            "OutDensity": True,
            "OutVolume": True,
            "OutSurfaceTension": True,
        }

    def test_save_NVT_bad_lamda_value(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The last value in the LambdaCoulomb variable list must be a 1.0",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 3,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 0.9],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The first value in the LambdaCoulomb variable list must be a 0.0",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 3,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.1, 0.3, 0.8, 0.9, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The last value in the LambdaVDW variable list must be a 1.0",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 3,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 0.9],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The first value in the LambdaVDW variable list must be a 0.0",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 3,
                    "LambdaVDW": [0.1, 0.2, 0.4, 0.9, 1.0],
                },
            )

    def test_save_NVT_bad_variables_part_1(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PRNG'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PRNG": [1]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeCHARMM'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeCHARMM": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeMie'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeMie": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeMARTINI'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeMARTINI": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RcutCoulomb_box_0'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutCoulomb_box_0": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['RcutCoulomb_box_1'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutCoulomb_box_1": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Pressure'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Pressure": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Rcut'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Rcut": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RcutLow'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutLow": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LRC'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"LRC": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Exclude'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Exclude": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Potential'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Potential": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Rswitch'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Rswitch": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ElectroStatic'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ElectroStatic": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Ewald'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Ewald": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CachedFourier'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CachedFourier": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Tolerance'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Tolerance": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Dielectric'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Dielectric": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['EqSteps'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"EqSteps": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['EqSteps'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"EqSteps": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['useConstantArea'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"useConstantArea": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ChemPot": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Fugacity": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_First'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_First": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Nth'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Nth": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Ang'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Ang": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Dih'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Dih": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutputName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutputName": 1},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CoordinatesFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CoordinatesFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RestartFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RestartFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CheckpointFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CheckpointFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ConsoleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ConsoleFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['BlockAverageFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"BlockAverageFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['HistogramFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"HistogramFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['DistName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"DistName": 1},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['HistName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"HistName": 1},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RunNumber'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RunNumber": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RunLetter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RunLetter": 1},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['SampleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"SampleFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutPressure'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutPressure": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutMolNum'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutMolNum": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutDensity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutDensity": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutVolume'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutVolume": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutSurfaceTension'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutSurfaceTension": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": "s",
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", "s"],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [["ETH"], 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [{"ETH": "1"}, 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": "s",
                    "LambdaVDW": [0.0, 0.1, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaVDW'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": "s",
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": "s",
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the "
            r"following variables need to be set, and not equal to "
            r"None: FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"FreeEnergyCalc": [True, 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScaleCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScaleCoulomb": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScalePower'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScalePower": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScaleAlpha'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScaleAlpha": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MinSigma'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MinSigma": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ExchangeVolumeDim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC_DataInput": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['DisFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"DisFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RotFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RotFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraSwapFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraSwapFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['SwapFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"SwapFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RegrowthFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RegrowthFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CrankShaftFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CrankShaftFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['VolFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"VolFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MultiParticleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MultiParticleFreq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-1Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-1Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-1Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-1Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-2Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-2Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-2Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-2Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-3Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-3Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-3Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-3Freq": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['XXXXXX'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"XXXXXX": "s"},
            )

    def test_save_NVT_bad_variables_part_2(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PRNG'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PRNG": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeCHARMM'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeCHARMM": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeMie'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeMie": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ParaTypeMARTINI'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ParaTypeMARTINI": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RcutCoulomb_box_0'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutCoulomb_box_0": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['RcutCoulomb_box_1'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutCoulomb_box_1": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Pressure'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Pressure": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Rcut'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Rcut": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RcutLow'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutLow": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LRC'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"LRC": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Exclude'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Exclude": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Potential'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Potential": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Rswitch'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Rswitch": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ElectroStatic'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ElectroStatic": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Ewald'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Ewald": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CachedFourier'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CachedFourier": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Tolerance'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Tolerance": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Dielectric'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Dielectric": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['EqSteps'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"EqSteps": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['AdjSteps'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"AdjSteps": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['useConstantArea'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"useConstantArea": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['FixVolBox0'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"FixVolBox0": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ChemPot": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Fugacity": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_First'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_First": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Nth'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Nth": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Ang'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Ang": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CBMC_Dih'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CBMC_Dih": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutputName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutputName": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CoordinatesFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CoordinatesFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RestartFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RestartFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CheckpointFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CheckpointFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ConsoleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ConsoleFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['BlockAverageFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"BlockAverageFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['HistogramFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"HistogramFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['DistName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"DistName": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['HistName'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"HistName": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RunNumber'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RunNumber": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RunLetter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RunLetter": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['SampleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"SampleFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutPressure'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutPressure": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutMolNum'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutMolNum": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutDensity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutDensity": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutVolume'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutVolume": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutSurfaceTension'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutSurfaceTension": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", []],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [["ETH"], 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [{"ETH": "1"}, 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": [],
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaVDW'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the "
            r"following variables need to be set, and not equal to "
            r"None: FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"FreeEnergyCalc": [True, 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScaleCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScaleCoulomb": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScalePower'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScalePower": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ScaleAlpha'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ScaleAlpha": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MinSigma'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MinSigma": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ExchangeVolumeDim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC_DataInput": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['DisFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"DisFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['DisFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"DisFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraSwapFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraSwapFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraSwapFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraSwapFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RegrowthFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RegrowthFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['CrankShaftFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"CrankShaftFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['VolFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"VolFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MultiParticleFreq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MultiParticleFreq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-1Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-1Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-1Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-1Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-2Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-2Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-2Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-2Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['IntraMEMC-3Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"IntraMEMC-3Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC-3Freq'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-3Freq": []},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the correct input variables where not provided for "
            r"the NVT ensemble. Please be sure to check that the keys in the "
            r"input variables dictionary \(input_variables_dict\) is correct, and "
            r"be aware that added spaces before or after the variable in any keys "
            r"will also give this warning. The bad variable inputs ensemble "
            r"inputs = \['XXXXXX'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_2.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"XXXXXX": []},
            )

    def test_save_NVT_bad_variables_part_5(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, 10000]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [False, 10000]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [1, 10000]},
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, 10000]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [False, 10000]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [1, 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": ["", 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [["x"], 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [{"s": 1}, 10000]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, 1.0]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, "x"]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, ["x"]]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [True, {"s": 1}]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['PressureCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_5.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"PressureCalc": [1, True]},
            )

    def test_save_NVT_bad_variables_part_6(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, True]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [False, True]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [False, False]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, True]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [False, True]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [False, False]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [1, True]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": ["", True]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [["x"], True]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [{"s": 1}, True]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, 1.0]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, "x"]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, ["x"]]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['OutEnergy'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_6.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"OutEnergy": [True, {"s": 1}]},
            )

    def test_save_NVT_bad_variables_part_7(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[1, 1, 1],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 0.9, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 0.8, 1.0],
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETH", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the following "
            r"variables need to be set, and not equal to None: FreeEnergyCalc, "
            r"MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the following "
            r"variables need to be set, and not equal to None: FreeEnergyCalc, "
            r"MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the following "
            r"variables need to be set, and not equal to None: FreeEnergyCalc, "
            r"MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: To utilize the free energy calculations all the following "
            r"variables need to be set, and not equal to None: FreeEnergyCalc, "
            r"MoleculeType, InitialState, LambdaVDW.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaCoulomb": [0.1, 0.3, 0.8, 1.0],
                },
            )

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # starting bad inputs for the Free engergy calcs side from not using all required variables
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [1, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": ["1", 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [["1"], 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [{"a": "1"}, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [False, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # starting bad inputs for the Free engergy calcs side from not using all required variables
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 1.0],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, "1"],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, ["1"]],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, {"a": "1"}],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FreeEnergyCalc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000, "s"],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        # start checking the MoleculeType variable for errors
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [1, 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [[1], 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": [{"a": "1"}, 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", "1"],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", ["1"]],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", {"a": "1"}],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MoleculeType'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETOa", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        # start checking the initial state variable
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": "s",
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": ["s"],
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": {"a": "1"},
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['InitialState'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1.0,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        # start checking the LamdaVDW variable
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaVDW'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": ["x", 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaVDW'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [[0.0], 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaVDW'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [{"a": "1"}, 0.2, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 1.0],
                },
            )

        # start testing the LambdaCoulomb
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": ["x", 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [[0.0], 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['LambdaCoulomb'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 1.0],
                    "LambdaCoulomb": [{"a": "1"}, 0.3, 1.0],
                },
            )
        with pytest.raises(
            ValueError,
            match=r"ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.1, 0.3, 0.8, 1.0],
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_7.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "FreeEnergyCalc": [True, 10000],
                    "MoleculeType": ["ETO", 1],
                    "InitialState": 1,
                    "LambdaVDW": [0.0, 0.1, 0.2, 0.4, 1.0],
                    "LambdaCoulomb": [0.0, 0.3, 0.8, 1.0],
                },
            )

    def test_save_NVT_bad_variables_part_8(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[4.0, 4.0, 4.0],
        )

        charmm = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol_box_0",
            structure_box_1=test_box_ethane_ethanol,
            filename_box_1="ethane_ethanol_box_1",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        test_box_ethane_ethanol = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[4.0, 4.0, 4.0],
        )

        charmm_NPT_NVT = Charmm(
            test_box_ethane_ethanol,
            "ethane_ethanol_box_0",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['RcutCoulomb_box_1'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"RcutCoulomb_box_1": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['FixVolBox0'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"FixVolBox0": "s"},
            )

        # test ExchangeVolumeDim for errors
        with pytest.raises(
            ValueError,
            match=r"The MEMC_DataInput variable is equal to None, but at least one "
            r"of the MEMC move ratios are all non-zero \(IntraMEMC-1Freq, "
            r"MEMC-1Freq, IntraMEMC-2Freq, MEMC-2Freq, IntraMEMC-3Freq, "
            r"and MEMC-3Freq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"MEMC-1Freq": 1},
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": [1.0, 1.0, 1.0]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": [1, 1, 1]},
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ExchangeVolumeDim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": ["s", 1.0, 1.0]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ExchangeVolumeDim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ExchangeVolumeDim": [1.0, [1.0], 1.0]},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ExchangeVolumeDim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "ExchangeVolumeDim": [1.0, 1.0, {"a": 1.0}]
                },
            )

        # testing failures and passes for MEMC_DataInput
        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "O1"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C2", "C1"], "ETO", ["O1", "C1"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "O1"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["O1", "C1"], "ETO", ["C2", "C1"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1.0, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        ["s", "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [[1], "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [{"a": "1"}, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETHa", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, 1, ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, [1], ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", [1, "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", [[1], "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", 1], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", [1]], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], 1, ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], [1], ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", [1, "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", [[1], "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", 1]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['MEMC_DataInput'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", [1]]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        # test the MEMC move ratios cant be set without specifying the MEMC move paramters ("MEMC_DataInput")
        with pytest.raises(
            ValueError,
            match=r"ERROR: The MEMC_DataInput variable is equal to None, but at least "
            r"one of the MEMC move ratios are all non-zero "
            r"\(IntraMEMC-1Freq, MEMC-1Freq, IntraMEMC-2Freq, MEMC-2Freq, "
            r"IntraMEMC-3Freq, and MEMC-3Freq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IntraMEMC-1Freq": 0.20,
                    "MEMC-1Freq": 0.20,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.20,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )

        # test some GCMC variable errors with Chempot and fugacity
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"ChemPot": "s"},
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_1.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={"Fugacity": "s"},
            )

        # testing the move frequency sum to 1 for all ensembles
        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
            r"Note: The sum that was manually entered may equal 1, but some "
            r"moves may not be valid for the provided ensemble. The moves that "
            r"are invalid for a given ensemble are set to zero. If the default "
            r"moves are not being used, all the move frequencies which do not have "
            r"default values of zero will need to be set manually so the sum equals "
            r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
            r"CrankShaftFreq, and VolFreq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.20,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.20,
                },
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.1,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
            r"Note: The sum that was manually entered may equal 1, but some "
            r"moves may not be valid for the provided ensemble. The moves that "
            r"are invalid for a given ensemble are set to zero. If the default "
            r"moves are not being used, all the move frequencies which do not have "
            r"default values of zero will need to be set manually so the sum equals "
            r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
            r"CrankShaftFreq, and VolFreq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GEMC_NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.1,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.20,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.20,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.20,
                },
            )

        try:
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.1,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                    "ChemPot": {"ETH": -4000, "ETO": 8000},
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
            r"Note: The sum that was manually entered may equal 1, but some "
            r"moves may not be valid for the provided ensemble. The moves that "
            r"are invalid for a given ensemble are set to zero. If the default "
            r"moves are not being used, all the move frequencies which do not have "
            r"default values of zero will need to be set manually so the sum equals "
            r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
            r"CrankShaftFreq, and VolFreq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.05,
                    "RegrowthFreq": 0.05,
                    "CrankShaftFreq": 0.1,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.20,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.20,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.20,
                    "Fugacity": {"ETH": 0, "ETO": 1.0},
                },
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.20,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
            r"Note: The sum that was manually entered may equal 1, but some "
            r"moves may not be valid for the provided ensemble. The moves that "
            r"are invalid for a given ensemble are set to zero. If the default "
            r"moves are not being used, all the move frequencies which do not have "
            r"default values of zero will need to be set manually so the sum equals "
            r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
            r"CrankShaftFreq, and VolFreq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.05,
                    "MultiParticleFreq": 0.05,
                    "IntraMEMC-1Freq": 0.20,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.19,
                    "MEMC-3Freq": 0.00,
                },
            )

        try:
            value = gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.20,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(
            ValueError,
            match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
            r"Note: The sum that was manually entered may equal 1, but some "
            r"moves may not be valid for the provided ensemble. The moves that "
            r"are invalid for a given ensemble are set to zero. If the default "
            r"moves are not being used, all the move frequencies which do not have "
            r"default values of zero will need to be set manually so the sum equals "
            r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
            r"CrankShaftFreq, and VolFreq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.20,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.21,
                    "MEMC-3Freq": 0.00,
                },
            )

        # test good values of Volume for NVT, and GCMC if set to zero
        try:
            value = gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.20,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.20,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "ChemPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.10,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.10,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.10,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # test come MEMC with GCMC
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 1,
                    "Fugacity": {1: 0, "ETO": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 1,
                    "Fugacity": {"ETH": -1, "ETO": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 1,
                    "Fugacity": {"ETH": "1", "ETO": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The MEMC_DataInput variable is not equal to None, "
            r"but all the MEMC move ratios are zero \(IntraMEMC-1Freq, MEMC-1Freq, "
            r"IntraMEMC-2Freq, MEMC-2Freq, IntraMEMC-3Freq, and MEMC-3Freq\).",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 1,
                    "Fugacity": {"ETH": 2, "ETO": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "Fugacity": {"ETH": 0, "XXX": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['Fugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "Fugacity": {"XXX": 0, "ETO": 1.0},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "ChemPot": {1: -4000, "ETO": -8000},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "ChemPot": {"XXX": -4000, "ETO": -8000},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "ChemPot": {"ETH": -4000, "XXX": -8000},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "ChemPot": {"ETH": "40", "ETO": -8000},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['ChemPot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 1,
                    "ChemPot": {"ETH": ["40"], "ETO": -8000},
                },
            )

        # test bad values of Volume for NVT, and GCMC
        with pytest.raises(
            ValueError,
            match=r"ERROR: The input variable VolFreq is non-zero \(0\). "
            r'VolFreq must be zero \(0\) for the "NVT", "GEMC_NVT", '
            r'and "GCMC" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.1,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The input variable VolFreq is non-zero \(0\). "
            r'VolFreq must be zero \(0\) for the "NVT", "GEMC_NVT", '
            r'and "GCMC" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "ChemPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.1,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        # test bad values of MEMC  for NVT, NPT
        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NVT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: All the MC move input variables must be non-zero \(0\) "
            r"for the SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq. "
            r"The SwapFreq, MEMC-1Freq, MEMC-2Freq, and MEMC-3Freq need to be zero "
            r'\(0\) for the "NVT" and "NPT" ensembles.',
        ):
            gomc_control.write_gomc_control_file(
                charmm_NPT_NVT,
                "test_save_NVT_bad_variables_part_8.conf",
                "NPT",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.10,
                },
            )

        # test good values of MEMC  with GCMC
        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "ChemPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.10,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.10,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "ChemPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.10,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.00,
                    "TargetedSwapFreq": 0.00,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]]
                    ],
                    "ChemPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreq": 0.05,
                    "RotFreq": 0.05,
                    "IntraSwapFreq": 0.05,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.10,
                    "CrankShaftFreq": 0.05,
                    "VolFreq": 0.0,
                    "MultiParticleFreq": 0.1,
                    "IntraMEMC-1Freq": 0.10,
                    "MEMC-1Freq": 0.00,
                    "IntraMEMC-2Freq": 0.20,
                    "MEMC-2Freq": 0.00,
                    "IntraMEMC-3Freq": 0.20,
                    "MEMC-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # try all case unspecific values
        try:
            value = gomc_control.write_gomc_control_file(
                charmm,
                "test_save_NVT_bad_variables_part_8.conf",
                "GCMC",
                10,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "MEMC_DataInput": [
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "C2"]],
                        [1, "ETH", ["C1", "C2"], "ETO", ["C1", "O1"]],
                    ],
                    "ChEmPot": {"ETH": -4000, "ETO": -8000},
                    "DisFreQ": 0.05,
                    "RotFreq": 0.05,
                    "InTraSwapFreq": 0.05,
                    "SwaPFreq": 0.00,
                    "ReGrowthFreq": 0.10,
                    "crankshaftfreq": 0.05,
                    "VolFrEQ": 0.0,
                    "MULtiParticleFreq": 0.1,
                    "IntRAMEMC-1Freq": 0.10,
                    "MEMC-1FREQ": 0.00,
                    "IntrAMEMC-2Freq": 0.20,
                    "MEMC-2FReq": 0.00,
                    "intramemc-3Freq": 0.20,
                    "memc-3Freq": 0.10,
                },
            )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

    def test_charmm_object_has_proper_no_boxes_for_ensemble_part_9(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_ethanol_liq = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[4, 4],
            box=[4.0, 4.0, 4.0],
        )

        test_box_ethane_ethanol_vap = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[1, 1],
            box=[8.0, 8.0, 8.0],
        )

        charmm_one_box = Charmm(
            test_box_ethane_ethanol_liq,
            "ethane_ethanol_1_box_liq",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        charmm_two_boxes = Charmm(
            test_box_ethane_ethanol_liq,
            "ethane_ethanol_2_boxes_liq",
            structure_box_1=test_box_ethane_ethanol_vap,
            filename_box_1="ethane_box_2_boxes_vap",
            ff_filename="ethane_ethanol",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        # test that it fails with the GEMC_NVT with only 1 box in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of {} is using a Charmm "
            r"object with one simulation boxes, and the {} ensemble only accepts "
            r"two boxes \(box 0 and box 1\).".format("GEMC_NVT", "GEMC_NVT"),
        ):

            gomc_control.write_gomc_control_file(
                charmm_one_box,
                "test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box",
                "GEMC_NVT",
                100,
                300,
                check_input_files_exist=False,
            )

        # test that it fails with the GEMC_NPT with only 1 box in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of {} is using a Charmm "
            r"object with one simulation boxes, and the {} ensemble only accepts "
            r"two boxes \(box 0 and box 1\).".format("GEMC_NPT", "GEMC_NPT"),
        ):

            gomc_control.write_gomc_control_file(
                charmm_one_box,
                "test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box",
                "GEMC_NPT",
                100,
                300,
                check_input_files_exist=False,
            )

        # test that it fails with the GCMC with only 1 box in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of {} is using a Charmm "
            r"object with one simulation boxes, and the {} ensemble only accepts "
            r"two boxes \(box 0 and box 1\).".format("GCMC", "GCMC"),
        ):

            gomc_control.write_gomc_control_file(
                charmm_one_box,
                "test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box",
                "GCMC",
                100,
                300,
                check_input_files_exist=False,
            )

        # test that it fails with the NVT with 2 boxes in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of {} is using a Charmm "
            r"object with two simulation boxes, and the {} ensemble only accepts "
            r"one box \(box 0\).".format("NVT", "NVT"),
        ):

            gomc_control.write_gomc_control_file(
                charmm_two_boxes,
                "test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
            )

        # test that it fails with the NPT with 2 boxes in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: The ensemble type selection of {} is using a Charmm "
            r"object with two simulation boxes, and the {} ensemble only accepts "
            r"one box \(box 0\).".format("NPT", "NPT"),
        ):
            gomc_control.write_gomc_control_file(
                charmm_two_boxes,
                "test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box",
                "NPT",
                100,
                300,
                check_input_files_exist=False,
            )

    def test_save_non_othoganol_writer(self):
        lattice_cif_ETV_triclinic = load_cif(
            file_or_path=get_fn("ETV_triclinic.cif")
        )
        ETV_triclinic_1_cell = lattice_cif_ETV_triclinic.populate(x=1, y=1, z=1)
        ETV_triclinic_1_cell.name = "ETV_1"
        ETV_triclinic_3_cell = lattice_cif_ETV_triclinic.populate(x=3, y=3, z=3)
        ETV_triclinic_3_cell.name = "ETV_3"

        charmm = Charmm(
            ETV_triclinic_1_cell,
            "ETV_triclinic_1_cell_box_0",
            structure_box_1=ETV_triclinic_3_cell,
            filename_box_1="ETV_triclinic_3_cell_box_1",
            ff_filename="ETV_triclinic_FF",
            forcefield_selection={
                ETV_triclinic_1_cell.name: get_fn(
                    "Charmm_writer_testing_only_zeolite.xml"
                ),
                ETV_triclinic_3_cell.name: get_fn(
                    "Charmm_writer_testing_only_zeolite.xml"
                ),
            },
            residues=[ETV_triclinic_1_cell.name, ETV_triclinic_3_cell.name],
            bead_to_atom_name_dict=None,
            fix_residue=[ETV_triclinic_1_cell.name, ETV_triclinic_3_cell.name],
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_non_othoganol_writer.conf",
            "GEMC_NVT",
            100000,
            300,
            check_input_files_exist=False,
        )

        with open("test_save_non_othoganol_writer.conf", "r") as fp:
            cell_vector_box_0_1_read = False
            cell_vector_box_0_2_read = False
            cell_vector_box_0_3_read = False
            cell_vector_box_1_1_read = False
            cell_vector_box_1_2_read = False
            cell_vector_box_1_3_read = False
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("CellBasisVector1 0"):
                    cell_vector_box_0_1_read = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "8.7503"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 0"):
                    cell_vector_box_0_2_read = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "-1.179131"
                    assert split_line[3] == "9.575585"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 0"):
                    cell_vector_box_0_3_read = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "-1.817231"
                    assert split_line[3] == "-3.027821"
                    assert split_line[4] == "9.645823"

                if line.startswith("CellBasisVector1 1"):
                    cell_vector_box_1_1_read = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "26.2509"
                    assert split_line[3] == "0.0"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector2 1"):
                    cell_vector_box_1_2_read = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "-3.537381"
                    assert split_line[3] == "28.726735"
                    assert split_line[4] == "0.0"

                elif line.startswith("CellBasisVector3 1"):
                    cell_vector_box_1_3_read = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "-5.451699"
                    assert split_line[3] == "-9.083469"
                    assert split_line[4] == "28.937455"

                else:
                    pass

        assert cell_vector_box_0_1_read == True
        assert cell_vector_box_0_2_read == True
        assert cell_vector_box_0_3_read == True
        assert cell_vector_box_1_1_read == True
        assert cell_vector_box_1_2_read == True
        assert cell_vector_box_1_3_read == True

    def test_box_vector_too_many_char(self):

        methane = mb.Compound(name="MET")
        methane_child_bead = mb.Compound(name="_CH4")
        methane.add(methane_child_bead, inherit_periodicity=False)

        methane_box_orth = mb.fill_box(
            compound=methane, n_compounds=10, box=[1, 2, 3]
        )

        charmm_bad_box_0 = Charmm(
            methane_box_orth,
            "methane_box_0_orth",
            ff_filename="methane_box_orth_bad_box_0_non_orth",
            residues=[methane.name],
            forcefield_selection="trappe-ua",
        )

        # set the vectors all too long
        charmm_bad_box_0.box_0_vectors[0][0] = -0.45678901234561
        charmm_bad_box_0.box_0_vectors[0][1] = -0.45678901234562
        charmm_bad_box_0.box_0_vectors[0][2] = -0.45678901234563
        charmm_bad_box_0.box_0_vectors[1][0] = -0.45678901234564
        charmm_bad_box_0.box_0_vectors[1][1] = -0.45678901234565
        charmm_bad_box_0.box_0_vectors[1][2] = -0.45678901234566
        charmm_bad_box_0.box_0_vectors[2][0] = -0.45678901234567
        charmm_bad_box_0.box_0_vectors[2][1] = -0.45678901234568
        charmm_bad_box_0.box_0_vectors[2][2] = -0.45678901234569

        charmm_bad_box_1 = Charmm(
            methane_box_orth,
            "methane_box_0_orth",
            structure_box_1=methane_box_orth,
            filename_box_1="methane_box_1_orth",
            ff_filename="methane_box_orth_bad_box_1_non_orth",
            residues=[methane.name],
            forcefield_selection="trappe-ua",
        )

        # set the vectors all too long
        charmm_bad_box_1.box_1_vectors[0][0] = -0.45678901234561
        charmm_bad_box_1.box_1_vectors[0][1] = -0.45678901234562
        charmm_bad_box_1.box_1_vectors[0][2] = -0.45678901234563
        charmm_bad_box_1.box_1_vectors[1][0] = -0.45678901234564
        charmm_bad_box_1.box_1_vectors[1][1] = -0.45678901234565
        charmm_bad_box_1.box_1_vectors[1][2] = -0.45678901234566
        charmm_bad_box_1.box_1_vectors[2][0] = -0.45678901234567
        charmm_bad_box_1.box_1_vectors[2][1] = -0.45678901234568
        charmm_bad_box_1.box_1_vectors[2][2] = -0.45678901234569

        # test that it fails with the GEMC_NVT with only 1 box in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: At lease one of the individual box {} vectors are too large "
            "or greater than {} characters."
            "".format(0, 16),
        ):

            gomc_control.write_gomc_control_file(
                charmm_bad_box_0,
                "test_box_vector_too_many_char_box_0",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
            )

        # test that it fails with the GEMC_NPT with only 1 box in the Charmm object
        with pytest.raises(
            ValueError,
            match=r"ERROR: At lease one of the individual box {} vectors are too large "
            "or greater than {} characters."
            "".format(1, 16),
        ):

            gomc_control.write_gomc_control_file(
                charmm_bad_box_1,
                "test_box_vector_too_many_char_box_1",
                "GCMC",
                100,
                300,
                check_input_files_exist=False,
            )

    def test_adjustment_steps_and_ff_psf_pdb_file_directory(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        changed_file_path = "../files"
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        # test the failure of the ff_psf_pdb_file_directory variable is not None or a string
        with pytest.raises(
            TypeError,
            match=f"ERROR: The {'ff_psf_pdb_file_directory'} variable for directly entering the "
            f"{'force field, pdb, and psf'} file directory and name is a {type(['x'])} and not a string.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_charmm_object_ff_psf_pdb_file_directory_not_a_sting",
                "GEMC_NPT",
                100,
                300,
                check_input_files_exist=False,
                ff_psf_pdb_file_directory=["x"],
            )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_adjustment_steps.conf",
            "GEMC_NVT",
            100000,
            500,
            check_input_files_exist=False,
            ff_psf_pdb_file_directory=changed_file_path,
            input_variables_dict={
                "PressureCalc": [True, 1],
                "AdjSteps": 2,
                "EqSteps": 3,
                "CoordinatesFreq": [True, 4],
                "RestartFreq": [True, 5],
                "CheckpointFreq": [True, 6],
                "ConsoleFreq": [True, 7],
                "BlockAverageFreq": [True, 8],
                "HistogramFreq": [True, 9],
                "SampleFreq": 11,
                "VDWGeometricSigma": True,
            },
        )

        with open("test_adjustment_steps.conf", "r") as fp:
            variables_read_dict = {
                "PressureCalc": False,
                "AdjSteps": False,
                "EqSteps": False,
                "CoordinatesFreq": False,
                "RestartFreq": False,
                "CheckpointFreq": False,
                "ConsoleFreq": False,
                "BlockAverageFreq": False,
                "HistogramFreq": False,
                "SampleFreq": False,
                "VDWGeometricSigma": False,
                "Parameters": False,
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "Coordinates_box_1": False,
                "Structure_box_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("PressureCalc "):
                    variables_read_dict["PressureCalc"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1"

                elif line.startswith("AdjSteps "):
                    variables_read_dict["AdjSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "2"

                elif line.startswith("EqSteps "):
                    variables_read_dict["EqSteps"] = True
                    split_line = line.split()
                    assert split_line[1] == "3"

                elif line.startswith("CoordinatesFreq "):
                    variables_read_dict["CoordinatesFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "4"

                elif line.startswith("RestartFreq "):
                    variables_read_dict["RestartFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "5"

                elif line.startswith("CheckpointFreq "):
                    variables_read_dict["CheckpointFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "6"

                elif line.startswith("ConsoleFreq "):
                    variables_read_dict["ConsoleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "7"

                elif line.startswith("BlockAverageFreq "):
                    variables_read_dict["BlockAverageFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "8"

                elif line.startswith("HistogramFreq "):
                    variables_read_dict["HistogramFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "9"

                elif line.startswith("SampleFreq "):
                    variables_read_dict["SampleFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "11"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "{}/ethane_FF.inp".format(
                        changed_file_path
                    )

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[2] == "{}/ethane_box_0.pdb".format(
                        changed_file_path
                    )

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[2] == "{}/ethane_box_0.psf".format(
                        changed_file_path
                    )

                elif line.startswith("Coordinates 1 "):
                    variables_read_dict["Coordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[2] == "{}/ethane_box_1.pdb".format(
                        changed_file_path
                    )

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure_box_1"] = True
                    split_line = line.split()
                    assert split_line[2] == "{}/ethane_box_1.psf".format(
                        changed_file_path
                    )

        assert variables_read_dict == {
            "PressureCalc": True,
            "AdjSteps": True,
            "EqSteps": True,
            "CoordinatesFreq": True,
            "RestartFreq": True,
            "CheckpointFreq": True,
            "ConsoleFreq": True,
            "BlockAverageFreq": True,
            "HistogramFreq": True,
            "SampleFreq": True,
            "VDWGeometricSigma": True,
            "Parameters": True,
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "Coordinates_box_1": True,
            "Structure_box_1": True,
        }

    def test_check_required_gomc_files_inp_exist_GEMC_NPT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"The {} with the file directory and name {}, "
            "does not exist.".format(
                "force field file or parameter file",
                "{}".format("ethane_FF.inp"),
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_gomc_ff_exists_GEMC_NPT",
                "GEMC_NPT",
                100,
                300,
                check_input_files_exist=True,
            )

    def test_check_required_gomc_files_pdb_exist_GEMC_NPT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        charmm.write_inp()

        with pytest.raises(
            ValueError,
            match=r"The {} with the file directory and name {}, "
            "does not exist.".format(
                "box 0 pdb file", "{}".format("ethane_box_0.pdb")
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_gomc_pdb_exists_GEMC_NPT",
                "GEMC_NPT",
                100,
                300,
                check_input_files_exist=True,
            )

    def test_check_required_gomc_files_psf_exist_GEMC_NPT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        charmm.write_inp()
        charmm.write_pdb()

        with pytest.raises(
            ValueError,
            match=r"The {} with the file directory and name {}, "
            "does not exist.".format(
                "box 0 psf file", "{}".format("ethane_box_0.psf")
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_gomc_psf_exists_GEMC_NPT",
                "GEMC_NPT",
                100,
                300,
                check_input_files_exist=True,
            )

    def test_check_required_gomc_files_pdb_exist_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        charmm.write_inp()

        with pytest.raises(
            ValueError,
            match=r"The {} with the file directory and name {}, "
            "does not exist.".format(
                "box 0 pdb file", "{}".format("ethane_box_0.pdb")
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_gomc_pdb_exists_NVT",
                "NVT",
                100,
                300,
                check_input_files_exist=True,
            )

    def test_check_required_gomc_files_psf_exist_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        charmm.write_inp()
        charmm.write_pdb()

        with pytest.raises(
            ValueError,
            match=r"The {} with the file directory and name {}, "
            "does not exist.".format(
                "box 0 psf file", "{}".format("ethane_box_0.psf")
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_gomc_psf_exists_NVT",
                "NVT",
                100,
                300,
                check_input_files_exist=True,
            )

    def test_check_restart_bool(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        restart_input = "XXXXX"
        with pytest.raises(
            TypeError,
            match=r"ERROR: The {} input is {} and needs to be a boolean \(i.e., True or False\)."
            "".format("Restart", type(restart_input)),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_checkpoint_error",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                Restart=restart_input,
            )

        restart_checkpoint_input = "XXXXX"
        with pytest.raises(
            TypeError,
            match=r"ERROR: The {} input is {} and needs to be a boolean \(i.e., True or False\)."
            "".format("RestartCheckpoint", type(restart_checkpoint_input)),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_checkpoint_error",
                "NVT",
                100,
                300,
                check_input_files_exist=False,
                RestartCheckpoint="XXXXX",
            )

        check_input_files_exist_input = "XXXXX"
        with pytest.raises(
            TypeError,
            match=r"ERROR: The {} input is {} and needs to be a boolean \(i.e., True or False\)."
            "".format(
                "check_input_files_exist", type(check_input_files_exist_input)
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "check_input_files_exist_error",
                "NVT",
                100,
                300,
                check_input_files_exist="XXXXX",
            )

    def test_restarting_dcd_and_binary_files_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_dcd_and_binary_files_NVT",
            "NVT",
            1000,
            300,
            Restart=True,
            check_input_files_exist=False,
            Coordinates_box_0="../test_files/NVT_toluene_box_0.pdb",
            Structure_box_0="../test_files/NVT_toluene_box_0.psf",
            binCoordinates_box_0="../test_files/NVT_toluene_box_0.coor",
            extendedSystem_box_0="../test_files/NVT_toluene_box_0.xsc",
            binVelocities_box_0="../test_files/NVT_toluene_box_0.vel",
            input_variables_dict={
                "VDWGeometricSigma": True,
                "DCDFreq": [True, 1000],
            },
        )

        with open("test_restarting_dcd_and_binary_files_NVT.conf", "r") as fp:
            variables_read_dict = {
                "VDWGeometricSigma": False,
                "DCDFreq": False,
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "binCoordinates_box_0": False,
                "extendedSystem_box_0": False,
                "binVelocities_box_0": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DCDFreq "):
                    variables_read_dict["DCDFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1000"

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_toluene_box_0.pdb"
                    )

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_toluene_box_0.psf"
                    )

                elif line.startswith("binCoordinates   0 "):
                    variables_read_dict["binCoordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_toluene_box_0.coor"
                    )

                elif line.startswith("extendedSystem 	0 "):
                    variables_read_dict["extendedSystem_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_toluene_box_0.xsc"
                    )

                elif line.startswith("binVelocities   	0"):
                    variables_read_dict["binVelocities_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_toluene_box_0.vel"
                    )

        assert variables_read_dict == {
            "Restart": True,
            "VDWGeometricSigma": True,
            "DCDFreq": True,
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "binCoordinates_box_0": True,
            "extendedSystem_box_0": True,
            "binVelocities_box_0": True,
        }

    def test_restarting_dcd_and_binary_files_GEMC_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_dcd_and_binary_files_GEMC_NVT",
            "GEMC_NVT",
            1000,
            300,
            Restart=True,
            check_input_files_exist=False,
            Coordinates_box_0="../test_files/NVT_ethane_box_0.pdb",
            Structure_box_0="../test_files/NVT_ethane_box_0.psf",
            binCoordinates_box_0="../test_files/NVT_ethane_box_0.coor",
            extendedSystem_box_0="../test_files/NVT_ethane_box_0.xsc",
            binVelocities_box_0="../test_files/NVT_ethane_box_0.vel",
            Coordinates_box_1="../test_files/NVT_ethane_box_1.pdb",
            Structure_box_1="../test_files/NVT_ethane_box_1.psf",
            binCoordinates_box_1="../test_files/NVT_ethane_box_1.coor",
            extendedSystem_box_1="../test_files/NVT_ethane_box_1.xsc",
            binVelocities_box_1="../test_files/NVT_ethane_box_1.vel",
            input_variables_dict={
                "VDWGeometricSigma": True,
                "DCDFreq": [True, 1000],
            },
        )

        with open(
            "test_restarting_dcd_and_binary_files_GEMC_NVT.conf", "r"
        ) as fp:
            variables_read_dict = {
                "VDWGeometricSigma": False,
                "DCDFreq": False,
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "binCoordinates_box_0": False,
                "extendedSystem_box_0": False,
                "binVelocities_box_0": False,
                "Coordinates_box_1": False,
                "Structure_box_1": False,
                "binCoordinates_box_1": False,
                "extendedSystem_box_1": False,
                "binVelocities_box_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("VDWGeometricSigma "):
                    variables_read_dict["VDWGeometricSigma"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DCDFreq "):
                    variables_read_dict["DCDFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"
                    assert split_line[2] == "1000"

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../test_files/NVT_ethane_box_0.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../test_files/NVT_ethane_box_0.psf"

                elif line.startswith("binCoordinates   0 "):
                    variables_read_dict["binCoordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files/NVT_ethane_box_0.coor"
                    )

                elif line.startswith("extendedSystem 	0 "):
                    variables_read_dict["extendedSystem_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../test_files/NVT_ethane_box_0.xsc"

                elif line.startswith("binVelocities   	0"):
                    variables_read_dict["binVelocities_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../test_files/NVT_ethane_box_0.vel"

                elif line.startswith("Coordinates 1 "):
                    variables_read_dict["Coordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../test_files/NVT_ethane_box_1.pdb"

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../test_files/NVT_ethane_box_1.psf"

                elif line.startswith("binCoordinates   1 "):
                    variables_read_dict["binCoordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert (
                        split_line[2] == "../test_files/NVT_ethane_box_1.coor"
                    )

                elif line.startswith("extendedSystem 	1 "):
                    variables_read_dict["extendedSystem_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../test_files/NVT_ethane_box_1.xsc"

                elif line.startswith("binVelocities   	1"):
                    variables_read_dict["binVelocities_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../test_files/NVT_ethane_box_1.vel"

        assert variables_read_dict == {
            "Restart": True,
            "VDWGeometricSigma": True,
            "DCDFreq": True,
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "binCoordinates_box_0": True,
            "extendedSystem_box_0": True,
            "binVelocities_box_0": True,
            "Coordinates_box_1": True,
            "Structure_box_1": True,
            "binCoordinates_box_1": True,
            "extendedSystem_box_1": True,
            "binVelocities_box_1": True,
        }

    def test_restarting_pdb_psf_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_NVT",
            "NVT",
            1000,
            300,
            ff_psf_pdb_file_directory="../Test",
            Parameters="../test_folder/new.par",
            Restart=True,
            RestartCheckpoint=True,
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open("test_restarting_pdb_psf_NVT.conf", "r") as fp:
            variables_read_dict = {
                "Parameters": False,
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "Restart": False,
                "RestartCheckpoint": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "../test_folder/new.par"

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.psf"

                elif line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("RestartCheckpoint "):
                    variables_read_dict["RestartCheckpoint"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

        assert variables_read_dict == {
            "Parameters": True,
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "Restart": True,
            "RestartCheckpoint": True,
        }

    def test_restarting_pdb_psf_NVT_only_rename_coordinates(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_NVT_only_rename_coordinates",
            "NVT",
            1000,
            300,
            ff_psf_pdb_file_directory="../Test",
            Restart=True,
            RestartCheckpoint=True,
            Coordinates_box_0="../test_files_1/NVT_toluene_box_0.pdb",
            Structure_box_0=None,
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open(
            "test_restarting_pdb_psf_NVT_only_rename_coordinates.conf", "r"
        ) as fp:
            variables_read_dict = {
                "Coordinates_box_0": False,
                "Structure_box_0": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files_1/NVT_toluene_box_0.pdb"
                    )

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.psf"

        assert variables_read_dict == {
            "Coordinates_box_0": True,
            "Structure_box_0": True,
        }

    def test_restarting_pdb_psf_NVT_only_rename_structure(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_NVT_only_rename_structure",
            "NVT",
            1000,
            300,
            ff_psf_pdb_file_directory=None,
            Restart=True,
            RestartCheckpoint=True,
            Coordinates_box_0=None,
            Structure_box_0="../test_files_2/NVT_toluene_box_0.psf",
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open(
            "test_restarting_pdb_psf_NVT_only_rename_structure.conf", "r"
        ) as fp:
            variables_read_dict = {
                "Coordinates_box_0": False,
                "Structure_box_0": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_box_0.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files_2/NVT_toluene_box_0.psf"
                    )

        assert variables_read_dict == {
            "Coordinates_box_0": True,
            "Structure_box_0": True,
        }

    def test_restarting_pdb_psf_GEMC_NVT_only_rename_coordinates(
        self, ethane_gomc
    ):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_GEMC_NVT_only_rename_coordinates",
            "GEMC_NVT",
            1000,
            300,
            ff_psf_pdb_file_directory=None,
            Restart=True,
            RestartCheckpoint=True,
            Coordinates_box_0="../test_files_1/NVT_toluene_box_0.pdb",
            Structure_box_0=None,
            Coordinates_box_1="../test_files_2/NVT_toluene_box_1.pdb",
            Structure_box_1=None,
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open(
            "test_restarting_pdb_psf_GEMC_NVT_only_rename_coordinates.conf", "r"
        ) as fp:
            variables_read_dict = {
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "Coordinates_box_1": False,
                "Structure_box_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files_1/NVT_toluene_box_0.pdb"
                    )

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "ethane_box_0.psf"

                elif line.startswith("Coordinates 1 "):
                    variables_read_dict["Coordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert (
                        split_line[2] == "../test_files_2/NVT_toluene_box_1.pdb"
                    )

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "ethane_box_1.psf"

        assert variables_read_dict == {
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "Coordinates_box_1": True,
            "Structure_box_1": True,
        }

    def test_restarting_pdb_psf_GEMC_NVT_only_rename_structure(
        self, ethane_gomc
    ):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_GEMC_NVT_only_rename_structure",
            "GEMC_NVT",
            1000,
            300,
            ff_psf_pdb_file_directory="../Test",
            Restart=True,
            RestartCheckpoint=True,
            Coordinates_box_0=None,
            Structure_box_0="../test_files_1/NVT_toluene_box_0.psf",
            Coordinates_box_1=None,
            Structure_box_1="../test_files_2/NVT_toluene_box_1.psf",
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open(
            "test_restarting_pdb_psf_GEMC_NVT_only_rename_structure.conf", "r"
        ) as fp:
            variables_read_dict = {
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "Coordinates_box_1": False,
                "Structure_box_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert (
                        split_line[2] == "../test_files_1/NVT_toluene_box_0.psf"
                    )

                if line.startswith("Coordinates 1 "):
                    variables_read_dict["Coordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../Test/ethane_box_1.pdb"

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert (
                        split_line[2] == "../test_files_2/NVT_toluene_box_1.psf"
                    )

        assert variables_read_dict == {
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "Coordinates_box_1": True,
            "Structure_box_1": True,
        }

    def test_restarting_pdb_psf_GEMC_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_restarting_pdb_psf_GEMC_NVT",
            "GEMC-NVT",
            1000,
            300,
            ff_psf_pdb_file_directory="../Test",
            Parameters="../test_folder/new.inp",
            Restart=True,
            RestartCheckpoint=True,
            check_input_files_exist=False,
            input_variables_dict={},
        )

        with open("test_restarting_pdb_psf_GEMC_NVT.conf", "r") as fp:
            variables_read_dict = {
                "Parameters": False,
                "Coordinates_box_0": False,
                "Structure_box_0": False,
                "Coordinates_box_1": False,
                "Structure_box_1": False,
                "Restart": False,
                "RestartCheckpoint": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("Parameters "):
                    variables_read_dict["Parameters"] = True
                    split_line = line.split()
                    assert split_line[1] == "../test_folder/new.inp"

                elif line.startswith("Coordinates 0 "):
                    variables_read_dict["Coordinates_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.pdb"

                elif line.startswith("Structure 0 "):
                    variables_read_dict["Structure_box_0"] = True
                    split_line = line.split()
                    assert split_line[1] == "0"
                    assert split_line[2] == "../Test/ethane_box_0.psf"

                elif line.startswith("Coordinates 1 "):
                    variables_read_dict["Coordinates_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../Test/ethane_box_1.pdb"

                elif line.startswith("Structure 1 "):
                    variables_read_dict["Structure_box_1"] = True
                    split_line = line.split()
                    assert split_line[1] == "1"
                    assert split_line[2] == "../Test/ethane_box_1.psf"

                elif line.startswith("Restart "):
                    variables_read_dict["Restart"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("RestartCheckpoint "):
                    variables_read_dict["RestartCheckpoint"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

        assert variables_read_dict == {
            "Parameters": True,
            "Coordinates_box_0": True,
            "Structure_box_0": True,
            "Coordinates_box_1": True,
            "Structure_box_1": True,
            "Restart": True,
            "RestartCheckpoint": True,
        }

    def test_failures_restarting_dcd_and_binary_files_NVT(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        with pytest.raises(
            ValueError,
            match="ERROR: To restart a simulation with the binary files both the coor "
            "and xsc files for box 0 must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0="../test_files/NVT_ethane_box_0.coor",
                extendedSystem_box_0=None,
                binVelocities_box_0=None,
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

        with pytest.raises(
            ValueError,
            match="ERROR: To restart a simulation with the binary files both the coor and "
            "xsc files for box 0 must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0=None,
                extendedSystem_box_0="../test_files/NVT_ethane_box_0.xsc",
                binVelocities_box_0=None,
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

        with pytest.raises(
            ValueError,
            match='ERROR: To restart a "NVT", "NPT" simulation with the '
            "velocity binary files, the velocity files for box 0 "
            "must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0=None,
                extendedSystem_box_0=None,
                binVelocities_box_0="../test_files/NVT_ethane_box_0.vel",
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

    def test_failures_restarting_dcd_and_binary_files_GEMC_NVT(
        self, ethane_gomc
    ):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )

        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )
        charmm.write_inp()
        charmm.write_pdb()
        charmm.write_psf()

        with pytest.raises(
            ValueError,
            match="ERROR: To restart a simulation with the binary files both the coor and "
            "xsc files for box 0 and box 1 must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0="../test_files/NVT_ethane_box_0.coor",
                extendedSystem_box_0="../test_files/NVT_ethane_box_0.xsc",
                binVelocities_box_0=None,
                binCoordinates_box_1="../test_files/NVT_ethane_box_1.coor",
                extendedSystem_box_1=None,
                binVelocities_box_1=None,
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

        with pytest.raises(
            ValueError,
            match="ERROR: To restart a simulation with the binary files both the coor and "
            "xsc files for box 0 and box 1 must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0="../test_files/NVT_ethane_box_0.coor",
                extendedSystem_box_0="../test_files/NVT_ethane_box_0.xsc",
                binVelocities_box_0=None,
                binCoordinates_box_1=None,
                extendedSystem_box_1="../test_files/NVT_ethane_box_0.xsc",
                binVelocities_box_1=None,
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

        with pytest.raises(
            ValueError,
            match='ERROR: To restart a "GEMC_NPT", "GEMC_NVT", "GCMC" simulation with the '
            "velocity binary files, both the velocity files for box 0 and box 1 "
            "must be provided.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                Restart=True,
                check_input_files_exist=False,
                binCoordinates_box_0="../test_files/NVT_ethane_box_0.coor",
                extendedSystem_box_0="../test_files/NVT_ethane_box_0.xsc",
                binVelocities_box_0="../test_files/NVT_ethane_box_0.vel",
                binCoordinates_box_1="../test_files/NVT_ethane_box_1.coor",
                extendedSystem_box_1="../test_files/NVT_ethane_box_1.xsc",
                binVelocities_box_1=None,
                input_variables_dict={
                    "VDWGeometricSigma": True,
                    "DCDFreq": [True, 1000],
                },
            )

        test_box_0_pdb = "XXXX"
        with pytest.raises(
            TypeError,
            match=r'ERROR: The {} variable expects a file extension of {}, but the actual file extension is "{}". '
            r"".format(
                "Coordinates_box_0",
                r"\['.pdb'\]",
                os.path.splitext(test_box_0_pdb)[-1],
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                check_input_files_exist=True,
                Coordinates_box_0=test_box_0_pdb,
                Structure_box_0="ethane_box_0.psf",
                Coordinates_box_1="ethane_box_0.pdb",
                Structure_box_1="ethane_box_1.psf",
            )

        test_box_1_pdb = "XXXX"
        with pytest.raises(
            TypeError,
            match=r'ERROR: The {} variable expects a file extension of {}, but the actual file extension is "{}". '
            r"".format(
                "Coordinates_box_1",
                r"\['.pdb'\]",
                os.path.splitext(test_box_1_pdb)[-1],
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                check_input_files_exist=True,
                Coordinates_box_0="ethane_box_0.pdb",
                Structure_box_0="ethane_box_0.psf",
                Coordinates_box_1=test_box_1_pdb,
                Structure_box_1="ethane_box_1.psf",
            )

        test_box_0_psf = "XXXX"
        with pytest.raises(
            TypeError,
            match=r'ERROR: The {} variable expects a file extension of {}, but the actual file extension is "{}". '
            r"".format(
                "Structure_box_0",
                r"\['.psf'\]",
                os.path.splitext(test_box_0_psf)[-1],
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                check_input_files_exist=True,
                Coordinates_box_0="ethane_box_0.pdb",
                Structure_box_0=test_box_0_psf,
                Coordinates_box_1="ethane_box_1.pdb",
                Structure_box_1="ethane_box_1.psf",
            )

        test_box_1_psf = "XXXX"
        with pytest.raises(
            TypeError,
            match=r'ERROR: The {} variable expects a file extension of {}, but the actual file extension is "{}". '
            r"".format(
                "Structure_box_1",
                r"\['.psf'\]",
                os.path.splitext(test_box_1_psf)[-1],
            ),
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_restart_inputs",
                "GEMC_NVT",
                1000,
                300,
                check_input_files_exist=True,
                Coordinates_box_0="ethane_box_0.pdb",
                Structure_box_0="ethane_box_0.psf",
                Coordinates_box_1="ethane_box_1.pdb",
                Structure_box_1=test_box_1_psf,
            )

            test_parameters = ["XXXX"]
            with pytest.raises(
                TypeError,
                match=r"ERROR: The {} variable for directly entering the "
                r"{} file directory and name is a {} and not a string.".format(
                    "Parameters", "force field", type(test_parameters)
                ),
            ):
                gomc_control.write_gomc_control_file(
                    charmm,
                    "test_restart_inputs",
                    "GEMC_NVT",
                    1000,
                    300,
                    check_input_files_exist=True,
                    Parameters=test_parameters,
                    Coordinates_box_0="ethane_box_0.pdb",
                    Structure_box_0="ethane_box_0.psf",
                    Coordinates_box_1="ethane_box_1.pdb",
                    Structure_box_1=test_box_1_psf,
                )

    def test_save_basic_NPT_use_ExpertMode(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane",
            ff_filename="ethane",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_NPT_ExpertMode.conf",
            "NPT",
            1000,
            500,
            ExpertMode=True,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.25,
                "RotFreq": 0.25,
                "IntraSwapFreq": 0.25,
                "SwapFreq": 0.0,
                "RegrowthFreq": 0.25,
                "VolFreq": 0.0,
            },
        )

        with open("test_save_basic_NPT_ExpertMode.conf", "r") as fp:
            variables_read_dict = {
                "ExpertMode": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                else:
                    pass

        assert variables_read_dict == {
            "ExpertMode": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
        }

    def test_save_basic_GCMC_use_ExpertMode(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GCMC_ExpertMode.conf",
            "GCMC",
            1000,
            500,
            ExpertMode=True,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.25,
                "RotFreq": 0.25,
                "IntraSwapFreq": 0.25,
                "SwapFreq": 0.0,
                "RegrowthFreq": 0.25,
                "VolFreq": 0.0,
                "Fugacity": {"ETH": 1.0},
            },
        )

        with open("test_save_basic_GCMC_ExpertMode.conf", "r") as fp:
            variables_read_dict = {
                "ExpertMode": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                else:
                    pass

        assert variables_read_dict == {
            "ExpertMode": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
        }

    def test_save_basic_GEMC_NVT_use_ExpertMode(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GEMC_NVT_ExpertMode.conf",
            "GEMC-NVT",
            1000,
            500,
            ExpertMode=True,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.25,
                "RotFreq": 0.25,
                "IntraSwapFreq": 0.25,
                "SwapFreq": 0.0,
                "RegrowthFreq": 0.25,
                "VolFreq": 0.0,
            },
        )

        with open("test_save_basic_GEMC_NVT_ExpertMode.conf", "r") as fp:
            variables_read_dict = {
                "ExpertMode": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "RegrowthFreq": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                else:
                    pass

        assert variables_read_dict == {
            "ExpertMode": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "RegrowthFreq": True,
        }

    def test_save_basic_GEMC_NPT_use_ExpertMode(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[2, 2, 2]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GEMC_NPT_ExpertMode.conf",
            "GEMC-NPT",
            1000,
            500,
            ExpertMode=True,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.25,
                "RotFreq": 0.25,
                "IntraSwapFreq": 0.15,
                "SwapFreq": 0.10,
                "RegrowthFreq": 0.25,
                "VolFreq": 0.0,
            },
        )

        with open("test_save_basic_GEMC_NPT_ExpertMode.conf", "r") as fp:
            variables_read_dict = {
                "ExpertMode": False,
                "DisFreq": False,
                "RotFreq": False,
                "IntraSwapFreq": False,
                "SwapFreq": False,
                "RegrowthFreq": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("ExpertMode"):
                    variables_read_dict["ExpertMode"] = True
                    split_line = line.split()
                    assert split_line[1] == "True"

                elif line.startswith("DisFreq "):
                    variables_read_dict["DisFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("RotFreq "):
                    variables_read_dict["RotFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                elif line.startswith("IntraSwapFreq "):
                    variables_read_dict["IntraSwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.15"

                elif line.startswith("SwapFreq "):
                    variables_read_dict["SwapFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.1"

                elif line.startswith("RegrowthFreq "):
                    variables_read_dict["RegrowthFreq"] = True
                    split_line = line.split()
                    assert split_line[1] == "0.25"

                else:
                    pass

        assert variables_read_dict == {
            "ExpertMode": True,
            "DisFreq": True,
            "RotFreq": True,
            "IntraSwapFreq": True,
            "SwapFreq": True,
            "RegrowthFreq": True,
        }

    def test_save_basic_GCMC_use_targetedswap_chempot(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )
        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc_ethanol_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GCMC_use_targetedswap_chempot.conf",
            "GCMC",
            1000,
            500,
            ExpertMode=False,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.20,
                "RotFreq": 0.20,
                "IntraSwapFreq": 0.10,
                "SwapFreq": 0.10,
                "RegrowthFreq": 0.20,
                "VolFreq": 0.0,
                "TargetedSwapFreq": 0.10,
                "IntraTargetedSwapFreq": 0.10,
                "ChemPot": {"ETH": -444, "ETO": -555},
                "TargetedSwap_DataInput": {
                    0: {
                        "SubVolumeType": "dynamic",
                        "SubVolumeBox": 0,
                        "SubVolumeCenterList": ["1-2", 3, 4],
                        "SubVolumeDim": [3, 2, 1],
                        "SubVolumeResidueKind": ["all"],
                        "SubVolumeRigidSwap": True,
                        "SubVolumePBC": "XY",
                        "SubVolumeChemPot": {"ETH": -222, "ETO": -22},
                    },
                    1: {
                        "SubVolumeType": "static",
                        "SubVolumeBox": 0,
                        "SubVolumeCenter": [2, 3, 4],
                        "SubVolumeDim": [4, 3, 2],
                        "SubVolumeResidueKind": "All",
                        "SubVolumeRigidSwap": False,
                        "SubVolumePBC": "XYZ",
                        "SubVolumeChemPot": {"ETH": -333, "ETO": -33},
                    },
                },
            },
        )

        with open(
            "test_save_basic_GCMC_use_targetedswap_chempot.conf", "r"
        ) as fp:
            variables_read_dict = {
                "TargetedSwapFreq": False,
                "IntraTargetedSwapFreq": False,
                "SubVolumeBox_number_0": False,
                "SubVolumeCenterList_number_0": False,
                "SubVolumeDim_number_0": False,
                "SubVolumeResidueKind_number_0": False,
                "SubVolumeRigidSwap_number_0": False,
                "SubVolumePBC_number_0": False,
                "SubVolumeChemPot_number_0_value_0": False,
                "SubVolumeChemPot_number_0_value_1": False,
                "SubVolumeBox_number_1": False,
                "SubVolumeCenter_number_1": False,
                "SubVolumeDim_number_1": False,
                "SubVolumeResidueKind_number_1": False,
                "SubVolumeRigidSwap_number_1": False,
                "SubVolumePBC_number_1": False,
                "SubVolumeChemPot_number_1_value_0": False,
                "SubVolumeChemPot_number_1_value_1": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("TargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["TargetedSwapFreq"] = True

                elif line.startswith("IntraTargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["IntraTargetedSwapFreq"] = True

                elif line.startswith("SubVolumeBox "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_1"] = True

                elif line.startswith("SubVolumeCenterList "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "1-2"
                        assert split_line[3] == "3"
                        assert split_line[4] == "4"
                        variables_read_dict[
                            "SubVolumeCenterList_number_0"
                        ] = True

                elif line.startswith("SubVolumeCenter "):
                    split_line = line.split()
                    if split_line[1] == "1":
                        assert split_line[2] == "2"
                        assert split_line[3] == "3"
                        assert split_line[4] == "4"
                        variables_read_dict["SubVolumeCenter_number_1"] = True

                elif line.startswith("SubVolumeDim "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "3"
                        assert split_line[3] == "2"
                        assert split_line[4] == "1"
                        variables_read_dict["SubVolumeDim_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "4"
                        assert split_line[3] == "3"
                        assert split_line[4] == "2"
                        variables_read_dict["SubVolumeDim_number_1"] = True

                elif line.startswith("SubVolumeResidueKind "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "ALL"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "ALL"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_1"
                        ] = True

                elif line.startswith("SubVolumeRigidSwap "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "True"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "False"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_1"
                        ] = True

                elif line.startswith("SubVolumePBC "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "XY"
                        variables_read_dict["SubVolumePBC_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "XYZ"
                        variables_read_dict["SubVolumePBC_number_1"] = True

                elif line.startswith("SubVolumeChemPot "):
                    split_line = line.split()
                    if split_line[1] == "0" and split_line[2] == "ETH":
                        assert split_line[3] == "-222"
                        variables_read_dict[
                            "SubVolumeChemPot_number_0_value_0"
                        ] = True
                    if split_line[1] == "0" and split_line[2] == "ETO":
                        assert split_line[3] == "-22"
                        variables_read_dict[
                            "SubVolumeChemPot_number_0_value_1"
                        ] = True
                    if split_line[1] == "1" and split_line[2] == "ETH":
                        assert split_line[3] == "-333"
                        variables_read_dict[
                            "SubVolumeChemPot_number_1_value_0"
                        ] = True
                    if split_line[1] == "1" and split_line[2] == "ETO":
                        assert split_line[3] == "-33"
                        variables_read_dict[
                            "SubVolumeChemPot_number_1_value_1"
                        ] = True

                else:
                    pass

        assert variables_read_dict == {
            "TargetedSwapFreq": True,
            "IntraTargetedSwapFreq": True,
            "SubVolumeBox_number_0": True,
            "SubVolumeCenterList_number_0": True,
            "SubVolumeDim_number_0": True,
            "SubVolumeResidueKind_number_0": True,
            "SubVolumeRigidSwap_number_0": True,
            "SubVolumePBC_number_0": True,
            "SubVolumeChemPot_number_0_value_0": True,
            "SubVolumeChemPot_number_0_value_1": True,
            "SubVolumeBox_number_1": True,
            "SubVolumeCenter_number_1": True,
            "SubVolumeDim_number_1": True,
            "SubVolumeResidueKind_number_1": True,
            "SubVolumeRigidSwap_number_1": True,
            "SubVolumePBC_number_1": True,
            "SubVolumeChemPot_number_1_value_0": True,
            "SubVolumeChemPot_number_1_value_1": True,
        }

    def test_save_basic_GCMC_use_targetedswap_fugacity(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )
        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc_ethanol_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GCMC_use_targetedswap_fugacity.conf",
            "GCMC",
            1000,
            500,
            ExpertMode=False,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.20,
                "RotFreq": 0.20,
                "IntraSwapFreq": 0.10,
                "SwapFreq": 0.10,
                "RegrowthFreq": 0.20,
                "VolFreq": 0.0,
                "TargetedSwapFreq": 0.10,
                "IntraTargetedSwapFreq": 0.10,
                "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                "TargetedSwap_DataInput": {
                    0: {
                        "SubVolumeType": "dynamic",
                        "SubVolumeBox": 0,
                        "SubVolumeCenterList": ["1-6", 7, 8],
                        "SubVolumeDim": [3, 2, 1],
                        "SubVolumeResidueKind": ["ETH", "ETO"],
                        "SubVolumeRigidSwap": False,
                        "SubVolumePBC": "XY",
                        "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                    },
                    1: {
                        "SubVolumeType": "static",
                        "SubVolumeBox": 0,
                        "SubVolumeCenter": [2, 3, 4],
                        "SubVolumeDim": [4, 3, 2],
                        "SubVolumeResidueKind": "ETH",
                        "SubVolumeFugacity": {"ETH": 3.33},
                    },
                },
            },
        )

        with open(
            "test_save_basic_GCMC_use_targetedswap_fugacity.conf", "r"
        ) as fp:
            variables_read_dict = {
                "TargetedSwapFreq": False,
                "IntraTargetedSwapFreq": False,
                "SubVolumeBox_number_0": False,
                "SubVolumeCenterList_number_0": False,
                "SubVolumeDim_number_0": False,
                "SubVolumeResidueKind_number_0": False,
                "SubVolumeRigidSwap_number_0": False,
                "SubVolumePBC_number_0": False,
                "SubVolumeFugacity_number_0_value_0": False,
                "SubVolumeFugacity_number_0_value_1": False,
                "SubVolumeBox_number_1": False,
                "SubVolumeCenter_number_1": False,
                "SubVolumeDim_number_1": False,
                "SubVolumeResidueKind_number_1": False,
                "SubVolumeRigidSwap_number_1": False,
                "SubVolumePBC_number_1": False,
                "SubVolumeFugacity_number_1_value_0": False,
            }
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("TargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["TargetedSwapFreq"] = True

                elif line.startswith("IntraTargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["IntraTargetedSwapFreq"] = True

                elif line.startswith("SubVolumeBox "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_1"] = True

                elif line.startswith("SubVolumeCenterList "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "1-6"
                        assert split_line[3] == "7"
                        assert split_line[4] == "8"
                        variables_read_dict[
                            "SubVolumeCenterList_number_0"
                        ] = True

                elif line.startswith("SubVolumeCenter "):
                    split_line = line.split()
                    if split_line[1] == "1":
                        assert split_line[2] == "2"
                        assert split_line[3] == "3"
                        assert split_line[4] == "4"
                        variables_read_dict["SubVolumeCenter_number_1"] = True

                elif line.startswith("SubVolumeDim "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "3"
                        assert split_line[3] == "2"
                        assert split_line[4] == "1"
                        variables_read_dict["SubVolumeDim_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "4"
                        assert split_line[3] == "3"
                        assert split_line[4] == "2"
                        variables_read_dict["SubVolumeDim_number_1"] = True

                elif line.startswith("SubVolumeResidueKind "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "ETH"
                        assert split_line[3] == "ETO"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "ETH"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_1"
                        ] = True

                elif line.startswith("SubVolumeRigidSwap "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "False"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "True"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_1"
                        ] = True

                elif line.startswith("SubVolumePBC "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "XY"
                        variables_read_dict["SubVolumePBC_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "XYZ"
                        variables_read_dict["SubVolumePBC_number_1"] = True

                elif line.startswith("SubVolumeFugacity "):
                    split_line = line.split()
                    if split_line[1] == "0" and split_line[2] == "ETH":
                        assert split_line[3] == "2.22"
                        variables_read_dict[
                            "SubVolumeFugacity_number_0_value_0"
                        ] = True
                    if split_line[1] == "0" and split_line[2] == "ETO":
                        assert split_line[3] == "0.22"
                        variables_read_dict[
                            "SubVolumeFugacity_number_0_value_1"
                        ] = True
                    if split_line[1] == "1" and split_line[2] == "ETH":
                        assert split_line[3] == "3.33"
                        variables_read_dict[
                            "SubVolumeFugacity_number_1_value_0"
                        ] = True

                else:
                    pass

        assert variables_read_dict == {
            "TargetedSwapFreq": True,
            "IntraTargetedSwapFreq": True,
            "SubVolumeBox_number_0": True,
            "SubVolumeCenterList_number_0": True,
            "SubVolumeDim_number_0": True,
            "SubVolumeResidueKind_number_0": True,
            "SubVolumeRigidSwap_number_0": True,
            "SubVolumePBC_number_0": True,
            "SubVolumeFugacity_number_0_value_0": True,
            "SubVolumeFugacity_number_0_value_1": True,
            "SubVolumeBox_number_1": True,
            "SubVolumeCenter_number_1": True,
            "SubVolumeDim_number_1": True,
            "SubVolumeResidueKind_number_1": True,
            "SubVolumeRigidSwap_number_1": True,
            "SubVolumePBC_number_1": True,
            "SubVolumeFugacity_number_1_value_0": True,
        }

    def test_save_basic_GEMC_NVT_use_targetedswap(
        self, ethane_gomc, ethanol_gomc
    ):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )
        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc_ethanol_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_GEMC_NVT_use_targetedswap.conf",
            "GEMC_NVT",
            1000,
            500,
            ExpertMode=False,
            check_input_files_exist=False,
            input_variables_dict={
                "IPC": True,
                "LRC": False,
                "DisFreq": 0.20,
                "RotFreq": 0.20,
                "IntraSwapFreq": 0.10,
                "SwapFreq": 0.10,
                "RegrowthFreq": 0.20,
                "VolFreq": 0.0,
                "TargetedSwapFreq": 0.10,
                "IntraTargetedSwapFreq": 0.10,
                "TargetedSwap_DataInput": {
                    0: {
                        "SubVolumeType": "dynamic",
                        "SubVolumeBox": 0,
                        "SubVolumeCenterList": ["1-6", 7, 8],
                        "SubVolumeDim": [3, 2, 1],
                        "SubVolumeResidueKind": ["ETH", "ETO"],
                        "SubVolumeRigidSwap": True,
                        "SubVolumePBC": "XY",
                    },
                    1: {
                        "SubVolumeType": "static",
                        "SubVolumeBox": 1,
                        "SubVolumeCenter": [2, 3, 4],
                        "SubVolumeDim": [4, 3, 2],
                        "SubVolumeResidueKind": "ETH",
                        "SubVolumeRigidSwap": False,
                        "SubVolumePBC": "XYZ",
                    },
                },
            },
        )

        with open("test_save_basic_GEMC_NVT_use_targetedswap.conf", "r") as fp:
            variables_read_dict = {
                "IPC": False,
                "LRC": False,
                "TargetedSwapFreq": False,
                "IntraTargetedSwapFreq": False,
                "SubVolumeBox_number_0": False,
                "SubVolumeCenterList_number_0": False,
                "SubVolumeDim_number_0": False,
                "SubVolumeResidueKind_number_0": False,
                "SubVolumeRigidSwap_number_0": False,
                "SubVolumePBC_number_0": False,
                "SubVolumeBox_number_1": False,
                "SubVolumeCenter_number_1": False,
                "SubVolumeDim_number_1": False,
                "SubVolumeResidueKind_number_1": False,
                "SubVolumeRigidSwap_number_1": False,
                "SubVolumePBC_number_1": False,
            }

            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("IPC "):
                    split_line = line.split()
                    assert split_line[1] == "True"
                    variables_read_dict["IPC"] = True

                elif line.startswith("LRC "):
                    split_line = line.split()
                    assert split_line[1] == "False"
                    variables_read_dict["LRC"] = True

                elif line.startswith("TargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["TargetedSwapFreq"] = True

                elif line.startswith("IntraTargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.1"
                    variables_read_dict["IntraTargetedSwapFreq"] = True

                elif line.startswith("SubVolumeBox "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "1"
                        variables_read_dict["SubVolumeBox_number_1"] = True

                elif line.startswith("SubVolumeCenterList "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "1-6"
                        assert split_line[3] == "7"
                        assert split_line[4] == "8"
                        variables_read_dict[
                            "SubVolumeCenterList_number_0"
                        ] = True

                elif line.startswith("SubVolumeCenter "):
                    split_line = line.split()
                    if split_line[1] == "1":
                        assert split_line[2] == "2"
                        assert split_line[3] == "3"
                        assert split_line[4] == "4"
                        variables_read_dict["SubVolumeCenter_number_1"] = True

                elif line.startswith("SubVolumeDim "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "3"
                        assert split_line[3] == "2"
                        assert split_line[4] == "1"
                        variables_read_dict["SubVolumeDim_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "4"
                        assert split_line[3] == "3"
                        assert split_line[4] == "2"
                        variables_read_dict["SubVolumeDim_number_1"] = True

                elif line.startswith("SubVolumeResidueKind "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "ETH"
                        assert split_line[3] == "ETO"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "ETH"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_1"
                        ] = True

                elif line.startswith("SubVolumeRigidSwap "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "True"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "False"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_1"
                        ] = True

                elif line.startswith("SubVolumePBC "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "XY"
                        variables_read_dict["SubVolumePBC_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "XYZ"
                        variables_read_dict["SubVolumePBC_number_1"] = True

                else:
                    pass

        assert variables_read_dict == {
            "IPC": True,
            "LRC": True,
            "TargetedSwapFreq": True,
            "IntraTargetedSwapFreq": True,
            "SubVolumeBox_number_0": True,
            "SubVolumeCenterList_number_0": True,
            "SubVolumeDim_number_0": True,
            "SubVolumeResidueKind_number_0": True,
            "SubVolumeRigidSwap_number_0": True,
            "SubVolumePBC_number_0": True,
            "SubVolumeBox_number_1": True,
            "SubVolumeCenter_number_1": True,
            "SubVolumeDim_number_1": True,
            "SubVolumeResidueKind_number_1": True,
            "SubVolumeRigidSwap_number_1": True,
            "SubVolumePBC_number_1": True,
        }

    def test_save_basic_NVT_use_targetedswap(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )
        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=None,
            filename_box_1=None,
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        gomc_control.write_gomc_control_file(
            charmm,
            "test_save_basic_NVT_use_targetedswap.conf",
            "NVT",
            1000,
            500,
            ExpertMode=False,
            check_input_files_exist=False,
            input_variables_dict={
                "DisFreq": 0.20,
                "RotFreq": 0.20,
                "IntraSwapFreq": 0.10,
                "SwapFreq": 0.00,
                "RegrowthFreq": 0.20,
                "VolFreq": 0.0,
                "TargetedSwapFreq": 0.00,
                "IntraTargetedSwapFreq": 0.30,
                "TargetedSwap_DataInput": {
                    0: {
                        "SubVolumeType": "dynamic",
                        "SubVolumeBox": 0,
                        "SubVolumeCenterList": ["1-6", 7, 8],
                        "SubVolumeDim": [3, 2, 1],
                        "SubVolumeResidueKind": ["ETH", "ETO"],
                        "SubVolumeRigidSwap": True,
                        "SubVolumePBC": "XY",
                    },
                    1: {
                        "SubVolumeType": "static",
                        "SubVolumeBox": 0,
                        "SubVolumeCenter": [2, 3, 4],
                        "SubVolumeDim": [4, 3, 2],
                        "SubVolumeResidueKind": "ETH",
                        "SubVolumeRigidSwap": False,
                        "SubVolumePBC": "XYZ",
                    },
                },
            },
        )

        with open("test_save_basic_NVT_use_targetedswap.conf", "r") as fp:
            variables_read_dict = {
                "IntraTargetedSwapFreq": False,
                "SubVolumeBox_number_0": False,
                "SubVolumeCenterList_number_0": False,
                "SubVolumeDim_number_0": False,
                "SubVolumeResidueKind_number_0": False,
                "SubVolumeRigidSwap_number_0": False,
                "SubVolumePBC_number_0": False,
                "SubVolumeBox_number_1": False,
                "SubVolumeCenter_number_1": False,
                "SubVolumeDim_number_1": False,
                "SubVolumeResidueKind_number_1": False,
                "SubVolumeRigidSwap_number_1": False,
                "SubVolumePBC_number_1": False,
            }

            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
                if line.startswith("IntraTargetedSwapFreq "):
                    split_line = line.split()
                    assert split_line[1] == "0.3"
                    variables_read_dict["IntraTargetedSwapFreq"] = True

                elif line.startswith("SubVolumeBox "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "0"
                        variables_read_dict["SubVolumeBox_number_1"] = True

                elif line.startswith("SubVolumeCenterList "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "1-6"
                        assert split_line[3] == "7"
                        assert split_line[4] == "8"
                        variables_read_dict[
                            "SubVolumeCenterList_number_0"
                        ] = True

                elif line.startswith("SubVolumeCenter "):
                    split_line = line.split()
                    if split_line[1] == "1":
                        assert split_line[2] == "2"
                        assert split_line[3] == "3"
                        assert split_line[4] == "4"
                        variables_read_dict["SubVolumeCenter_number_1"] = True

                elif line.startswith("SubVolumeDim "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "3"
                        assert split_line[3] == "2"
                        assert split_line[4] == "1"
                        variables_read_dict["SubVolumeDim_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "4"
                        assert split_line[3] == "3"
                        assert split_line[4] == "2"
                        variables_read_dict["SubVolumeDim_number_1"] = True

                elif line.startswith("SubVolumeResidueKind "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "ETH"
                        assert split_line[3] == "ETO"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "ETH"
                        variables_read_dict[
                            "SubVolumeResidueKind_number_1"
                        ] = True

                elif line.startswith("SubVolumeRigidSwap "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "True"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_0"
                        ] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "False"
                        variables_read_dict[
                            "SubVolumeRigidSwap_number_1"
                        ] = True

                elif line.startswith("SubVolumePBC "):
                    split_line = line.split()
                    if split_line[1] == "0":
                        assert split_line[2] == "XY"
                        variables_read_dict["SubVolumePBC_number_0"] = True
                    if split_line[1] == "1":
                        assert split_line[2] == "XYZ"
                        variables_read_dict["SubVolumePBC_number_1"] = True

                else:
                    pass

        assert variables_read_dict == {
            "IntraTargetedSwapFreq": True,
            "SubVolumeBox_number_0": True,
            "SubVolumeCenterList_number_0": True,
            "SubVolumeDim_number_0": True,
            "SubVolumeResidueKind_number_0": True,
            "SubVolumeRigidSwap_number_0": True,
            "SubVolumePBC_number_0": True,
            "SubVolumeBox_number_1": True,
            "SubVolumeCenter_number_1": True,
            "SubVolumeDim_number_1": True,
            "SubVolumeResidueKind_number_1": True,
            "SubVolumeRigidSwap_number_1": True,
            "SubVolumePBC_number_1": True,
        }

    def test_failures_targetedswap_GCMC(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )

        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc_ethanol_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The TargetedSwap_DataInput variable is equal to None, "
            r"but the TargetedSwapFreq or IntraTargetedSwapFreq move ratio is non-zero.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.00,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.20,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                },
            )

        with pytest.raises(
            ValueError,
            match=r"The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            r" or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        "x": {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumefugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumetype'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "s",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumetype'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": 0,
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": "x",
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 2,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": ["s"],
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6.0", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": 3,
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": "x",
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": "x",
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": -2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumedim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [-3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumedim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": ["s", 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumedim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": "s",
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumedim'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 2, 3, 4],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": [0, "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": "x",
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": 0,
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["x"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumerigidswap'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": 0,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumerigidswap'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": "s",
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumerigidswap'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": [True],
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumepbc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "s",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumepbc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": 0,
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumepbc'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": [0],
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"X": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ET": -2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ET": 2.22, "ETO": "x"},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": ["s"],
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": "s",
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": 1,
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match="Both ChemPot and Fugacity were used in the "
            "TargetedSwap_DataInput dictionaries "
            "and in the standard GOMC swap inputs. "
            "However, only ChemPot or Fugacity may be used, not both.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match="Both ChemPot and Fugacity were used in the "
            "TargetedSwap_DataInput dictionaries "
            "and in the standard GOMC swap inputs. "
            "However, only ChemPot or Fugacity may be used, not both.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match="Both ChemPot and Fugacity were used in the "
            "TargetedSwap_DataInput dictionaries. "
            "However, only ChemPot or Fugacity may be used, not both.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 2.22, "ETO": 0.22},
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 0,
                            "SubVolumeCenter": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "subVolumetype": "static",
                            "subVolumeBox": 0,
                            "SubVolumecenter": 0,
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumechemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 0,
                            "SubVolumeCenter": "0",
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        # missing the subvolumecenter variable
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 0,
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        # missing the subvolumecenterlist variable
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "Swapfreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumepbc": "XY",
                            "SubVolumechemPot": {"ETH": 4.44},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenter'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 0,
                            "SubVolumeCenter": [1, 1, 1, 1],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": "All",
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            r" or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            0: "static",
                            "SubVolumeBox": 0,
                            "SubVolumeCenter": [1, 1, 1],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": "all",
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            r" or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 0,
                            "SubVolumeCenter": [1, 1, 1],
                            "x": [3, 2, 1],
                            "SubVolumeResidueKind": "ETH",
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            r" or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidswap": True,
                            "Subvolume": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": ["ETO"],
                            "SubVolumeRigidswap": True,
                            "SubvolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumeresiduekind'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": "ETH",
                            "SubVolumeRigidswap": True,
                            "SubvolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 1,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": "all",
                            "SubVolumeRigidswap": True,
                            "SubvolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.10,
                    "IntraTargetedSwapFreq": 0.10,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 1,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": "all",
                            "SubVolumeRigidswap": True,
                            "SubVolumeFugacity": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        # "SubVolumeFugacity" not a list
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumefugacity'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.10,
                    "IntraTargetedSwapFreq": 0.10,
                    "Fugacity": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": "all",
                            "SubVolumeRigidswap": True,
                            "SubVolumeFugacity": ["ETH"],
                        },
                    },
                },
            )

        # "SubVolumeFugacity" not a list
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumechempot'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GCMC",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.10,
                    "IntraTargetedSwapFreq": 0.10,
                    "ChemPot": {"ETH": 4.44, "ETO": 5.55},
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "Subvolumebox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [1, 1, 1],
                            "SubVolumeresidueKind": "all",
                            "SubVolumeRigidswap": True,
                            "SubVolumeChemPot": ["ETH"],
                        },
                    },
                },
            )

    def test_failures_targetedswap_GEMC_NVT(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )

        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=test_box_ethane_gomc_ethanol_gomc,
            filename_box_1="ethane_box_1",
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The TargetedSwap_DataInput variable is equal to None, "
            r"but the TargetedSwapFreq or IntraTargetedSwapFreq move ratio is non-zero.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GEMC_NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.00,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.20,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.10,
                    "IntraTargetedSwapFreq": 0.10,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            r" or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GEMC_NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        "x": {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match="Either the ChemPot and Fugacity were used in the "
            "TargetedSwap_DataInput dictionaries, "
            "which can not be used for the 'NPT', 'NVT', 'GEMC_NVT', "
            "or 'GEMC_NPT' ensembles.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GEMC_NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeChemPot": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        # This is missing the "SubVolumeDim" in the  "TargetedSwap_DataInput"
        with pytest.raises(
            ValueError,
            match=f"The TargetedSwap_DataInput dictionaries do not have all the required "
            f"keys or inputs per the subvolumetype and specified ensemble. "
            f"Remember that the 'subvolumetype' values are 'static' and 'dynamic', "
            f"which must only have the cooresponding "
            f"'subvolumecenter' and 'subvolumecenterlist' values, respectively,"
            f" in each subvolume.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GEMC_NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 1,
                            "SubVolumeCenter": [1, 7, 8],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        # This box "SubVolumeBox" is not 0 or 1
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "GEMC_NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.20,
                    "SwapFreq": 0.20,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "static",
                            "SubVolumeBox": 4,
                            "SubVolumeCenter": [1, 7, 8],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

    def test_failures_targetedswap_NVT(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_gomc_ethanol_gomc = mb.fill_box(
            compound=[ethane_gomc, ethanol_gomc],
            n_compounds=[10, 10],
            box=[2, 2, 2],
        )

        charmm = Charmm(
            test_box_ethane_gomc_ethanol_gomc,
            "ethane_box_0",
            ff_filename="ethane",
            structure_box_1=None,
            filename_box_1=None,
            residues=[ethane_gomc.name, ethanol_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The TargetedSwap_DataInput variable is equal to None, "
            r"but the TargetedSwapFreq or IntraTargetedSwapFreq move ratio is non-zero.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.00,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.20,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumebox'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 1,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        with pytest.raises(
            ValueError,
            match="Either the ChemPot and Fugacity were used in the "
            "TargetedSwap_DataInput dictionaries, "
            "which can not be used for the 'NPT', 'NVT', 'GEMC_NVT', "
            "or 'GEMC_NPT' ensembles.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                            "SubVolumeFugacity": {"ETH": 4.44, "ETO": 5.55},
                        },
                    },
                },
            )

        # This is missing the "SubVolumeDim" in the  "TargetedSwap_DataInput"
        with pytest.raises(
            ValueError,
            match=f"The TargetedSwap_DataInput dictionaries do not have all the required "
            f"keys or inputs per the subvolumetype and specified ensemble. "
            f"Remember that the 'subvolumetype' values are 'static' and 'dynamic', "
            f"which must only have the cooresponding "
            f"'subvolumecenter' and 'subvolumecenterlist' values, respectively,"
            f" in each subvolume.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": ["1-6", 7, 8],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        # TargetedSwap_DataInput is a list not a dict
        with pytest.raises(
            ValueError,
            match="The TargetedSwap_DataInput is not formatted correctly as a dictionary"
            " or has the wrong input keys, values, or types.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": ["s"],
                },
            )

        # This "SubVolumeCenterList" is negative
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": [-7, 8, "1-6"],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        # This "SubVolumeCenterList" is missing a value in the - separator
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": [1, 8, "1-"],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

        # This "SubVolumeCenterList" is same value on both sides of - symbol
        with pytest.raises(
            ValueError,
            match=r"ERROR: The following input variables have "
            r"bad values \(check spelling and for empty spaces in the keys or that "
            r"the values are in the correct form with the acceptable values\)"
            r": \['subvolumecenterlist'\]",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_failures_targetedswap",
                "NVT",
                1000,
                500,
                ExpertMode=False,
                check_input_files_exist=False,
                input_variables_dict={
                    "DisFreq": 0.20,
                    "RotFreq": 0.20,
                    "IntraSwapFreq": 0.40,
                    "SwapFreq": 0.00,
                    "RegrowthFreq": 0.00,
                    "VolFreq": 0.0,
                    "TargetedSwapFreq": 0.00,
                    "IntraTargetedSwapFreq": 0.20,
                    "TargetedSwap_DataInput": {
                        0: {
                            "SubVolumeType": "dynamic",
                            "SubVolumeBox": 0,
                            "SubVolumeCenterList": [1, 8, "1-1"],
                            "SubVolumeDim": [3, 2, 1],
                            "SubVolumeResidueKind": ["ETH", "ETO"],
                            "SubVolumeRigidSwap": True,
                            "SubVolumePBC": "XY",
                        },
                    },
                },
            )

    def test_IPC_input_errors(self, ethane_gomc):
        test_box_ethane_gomc = mb.fill_box(
            compound=[ethane_gomc], n_compounds=[1], box=[1, 1, 1]
        )
        charmm = Charmm(
            test_box_ethane_gomc,
            "ethane_box_0",
            structure_box_1=None,
            filename_box_1=None,
            ff_filename="ethane_FF",
            residues=[ethane_gomc.name],
            forcefield_selection="oplsaa",
        )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The impulse correction term \(IPC\) can not be set as True "
            r"if the LRC=True or the Potential is SHIFT or SWITCH.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_IPC_input_errors",
                "NVT",
                1000,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IPC": True,
                    "LRC": True,
                    "Potential": "VDW",
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The impulse correction term \(IPC\) can not be set as True "
            r"if the LRC=True or the Potential is SHIFT or SWITCH.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_IPC_input_errors",
                "NVT",
                1000,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IPC": True,
                    "LRC": False,
                    "Potential": "SHIFT",
                },
            )

        with pytest.raises(
            ValueError,
            match=r"ERROR: The impulse correction term \(IPC\) can not be set as True "
            r"if the LRC=True or the Potential is SHIFT or SWITCH.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_IPC_input_errors",
                "NVT",
                1000,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IPC": True,
                    "LRC": False,
                    "Potential": "SWITCH",
                },
            )

        with pytest.warns(
            UserWarning,
            match=r"WARNING: The impulse correction term \(IPC\) is False, but likely needs to be True, "
            r"as the LRC=False when the Potential is VDW or EXP6.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_IPC_input_errors",
                "NVT",
                1000,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IPC": False,
                    "LRC": False,
                    "Potential": "VDW",
                },
            )

        with pytest.warns(
            UserWarning,
            match=r"WARNING: The impulse correction term \(IPC\) is False, but likely needs to be True, "
            r"as the LRC=False when the Potential is VDW or EXP6.",
        ):
            gomc_control.write_gomc_control_file(
                charmm,
                "test_IPC_input_errors",
                "NVT",
                1000,
                300,
                check_input_files_exist=False,
                input_variables_dict={
                    "IPC": False,
                    "LRC": False,
                    "Potential": "EXP6",
                },
            )
