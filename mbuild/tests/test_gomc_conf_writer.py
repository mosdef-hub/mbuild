import pytest
import mbuild as mb
import mbuild.formats.gomc_conf_writer as gomc_control

from mbuild.formats.charmm_writer import Charmm
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer



@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestGOMCControlFileWriter(BaseTest):

    def test_dict_keys_to_list(self,):
        dict = {'a': '1', 'b': '2', 'c': '3'}
        keys = gomc_control.dict_keys_to_list(dict)

        assert keys == ['a', 'b', 'c']


    def test_get_required_data(self):
        value = gomc_control._get_required_data(description=False)
        assert value == ['charmm_object', 'ensemble_type', 'RunSteps', 'Temperature']

        value = gomc_control._get_required_data(description=True)
        assert gomc_control.dict_keys_to_list(value) == ['charmm_object', 'ensemble_type', 'RunSteps', 'Temperature']


    def test_get_all_possible_input_variable(self):
        value = gomc_control._get_all_possible_input_variables(description=False)
        assert value == ['Restart', 'RestartCheckpoint', 'PRNG', 'ParaTypeCHARMM', 'ParaTypeMie',
                         'ParaTypeMARTINI', 'RcutCoulomb_box_0', 'RcutCoulomb_box_1', 'Pressure',
                         'Rcut', 'RcutLow', 'LRC', 'Exclude', 'Potential', 'Rswitch', 'VDWGeometricSigma',
                         'ElectroStatic', 'Ewald', 'CachedFourier', 'Tolerance', 'Dielectric', 'PressureCalc',
                         'EqSteps', 'AdjSteps', 'useConstantArea', 'FixVolBox0', 'ChemPot', 'Fugacity',
                         'CBMC_First', 'CBMC_Nth', 'CBMC_Ang', 'CBMC_Dih', 'OutputName', 'CoordinatesFreq',
                         'RestartFreq', 'CheckpointFreq', 'ConsoleFreq', 'BlockAverageFreq', 'HistogramFreq',
                         'DistName', 'HistName', 'RunNumber', 'RunLetter', 'SampleFreq', 'OutEnergy',
                         'OutPressure', 'OutMolNumber', 'OutDensity', 'OutVolume', 'OutSurfaceTension',
                         'FreeEnergyCalc', 'MoleculeType', 'InitialState', 'LambdaVDW', 'LambdaCoulomb',
                         'ScaleCoulomb', 'ScalePower', 'ScaleAlpha', 'MinSigma', 'DisFreq', 'RotFreq',
                         'IntraSwapFreq', 'SwapFreq', 'RegrowthFreq', 'CrankShaftFreq', 'VolFreq',
                         'MultiParticleFreq', 'IntraMEMC-1Freq', 'MEMC-1Freq', 'IntraMEMC-2Freq', 'MEMC-2Freq',
                         'IntraMEMC-3Freq', 'MEMC-3Freq', 'ExchangeVolumeDim', 'MEMC_DataInput']

        value = gomc_control._get_all_possible_input_variables(description=True)
        assert gomc_control.dict_keys_to_list(value) == ['Restart', 'RestartCheckpoint', 'PRNG', 'ParaTypeCHARMM',
                                                         'ParaTypeMie', 'ParaTypeMARTINI', 'RcutCoulomb_box_0',
                                                         'RcutCoulomb_box_1', 'Pressure', 'Rcut', 'RcutLow',
                                                         'LRC', 'Exclude', 'Potential', 'Rswitch', 'VDWGeometricSigma',
                                                         'ElectroStatic', 'Ewald', 'CachedFourier', 'Tolerance',
                                                         'Dielectric', 'PressureCalc', 'EqSteps', 'AdjSteps',
                                                         'useConstantArea', 'FixVolBox0', 'ChemPot', 'Fugacity',
                                                         'CBMC_First', 'CBMC_Nth', 'CBMC_Ang', 'CBMC_Dih',
                                                         'OutputName', 'CoordinatesFreq', 'RestartFreq',
                                                         'CheckpointFreq', 'ConsoleFreq', 'BlockAverageFreq',
                                                         'HistogramFreq', 'DistName', 'HistName', 'RunNumber',
                                                         'RunLetter', 'SampleFreq', 'OutEnergy', 'OutPressure',
                                                         'OutMolNumber', 'OutDensity', 'OutVolume',
                                                         'OutSurfaceTension', 'FreeEnergyCalc', 'MoleculeType',
                                                         'InitialState', 'LambdaVDW', 'LambdaCoulomb','ScaleCoulomb',
                                                         'ScalePower', 'ScaleAlpha', 'MinSigma', 'DisFreq', 'RotFreq',
                                                         'IntraSwapFreq', 'SwapFreq', 'RegrowthFreq', 'CrankShaftFreq',
                                                         'VolFreq', 'MultiParticleFreq', 'IntraMEMC-1Freq',
                                                         'MEMC-1Freq', 'IntraMEMC-2Freq', 'MEMC-2Freq',
                                                         'IntraMEMC-3Freq', 'MEMC-3Freq', 'ExchangeVolumeDim',
                                                         'MEMC_DataInput']


    def test_get_default_variables_dict(self):
        value = gomc_control._get_default_variables_dict()
        assert gomc_control.dict_keys_to_list(value) == ['Restart', 'RestartCheckpoint', 'PRNG', 'ParaTypeCHARMM',
                                                         'ParaTypeMie', 'ParaTypeMARTINI', 'RcutCoulomb_box_0',
                                                         'RcutCoulomb_box_1', 'Pressure', 'Rcut', 'RcutLow', 'LRC',
                                                         'Exclude', 'coul_1_4_scaling', 'Potential', 'Rswitch',
                                                         'ElectroStatic', 'Ewald', 'CachedFourier', 'Tolerance',
                                                         'Dielectric', 'PressureCalc', 'EqSteps', 'AdjSteps',
                                                         'useConstantArea', 'FixVolBox0', 'ChemPot', 'Fugacity',
                                                         'CBMC_First', 'CBMC_Nth', 'CBMC_Ang', 'CBMC_Dih',
                                                         'OutputName', 'CoordinatesFreq', 'RestartFreq',
                                                         'CheckpointFreq', 'ConsoleFreq', 'BlockAverageFreq',
                                                         'HistogramFreq', 'DistName', 'HistName', 'RunNumber',
                                                         'RunLetter', 'SampleFreq', 'OutEnergy', 'OutPressure',
                                                         'OutMolNumber', 'OutDensity', 'OutVolume',
                                                         'OutSurfaceTension', 'FreeEnergyCalc', 'MoleculeType',
                                                         'InitialState', 'LambdaVDW', 'LambdaCoulomb',
                                                         'ScaleCoulomb', 'ScalePower', 'ScaleAlpha', 'MinSigma',
                                                         'ExchangeVolumeDim', 'MEMC_DataInput', 'DisFreq',
                                                         'RotFreq', 'IntraSwapFreq', 'SwapFreq', 'RegrowthFreq',
                                                         'CrankShaftFreq', 'VolFreq', 'MultiParticleFreq',
                                                         'IntraMEMC-1Freq', 'MEMC-1Freq', 'IntraMEMC-2Freq',
                                                         'MEMC-2Freq', 'IntraMEMC-3Freq', 'MEMC-3Freq']






    def test_print_ensemble_info(self):
        try:
            gomc_control.print_valid_ensemble_input_variables('NVT', description=True)
            gomc_control.print_required_ensemble_files('NVT', description=True)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('NVT', description=False)
            gomc_control.print_required_ensemble_files('NVT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('NPT', description=True)
            gomc_control.print_required_ensemble_files('NPT', description=True)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('NPT', description=False)
            gomc_control.print_required_ensemble_files('NPT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NVT', description=True)
            gomc_control.print_required_ensemble_files('GEMC_NVT', description=True)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NVT', description=False)
            gomc_control.print_required_ensemble_files('GEMC_NVT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NPT', description=True)
            gomc_control.print_required_ensemble_files('GEMC_NPT', description=True)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NPT', description=False)
            gomc_control.print_required_ensemble_files('GEMC_NPT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GCMC', description=True)
            gomc_control.print_required_ensemble_files('GCMC', description=True)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GCMC', description=False)
            gomc_control.print_required_ensemble_files('GCMC', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"


    def test_save_basic_NVT(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'ethane', FF_filename='ethane',
                        residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_NVT.conf', 'NVT', 10, 300,
                                             )

        out_GOMC = open('test_save_basic_NVT.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if line.startswith('Restart '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('PRNG '):
                split_line = line.split()
                assert split_line[1] == 'RANDOM'

            elif line.startswith('ParaTypeCHARMM '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('Parameters '):
                split_line = line.split()
                assert split_line[1] == 'ethane.inp'

            elif line.startswith('Coordinates '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane.pdb'

            elif line.startswith('Structure '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane.psf'

            elif line.startswith('Temperature '):
                split_line = line.split()
                assert split_line[1] == '300'

            elif line.startswith('Potential '):
                split_line = line.split()
                assert split_line[1] == 'VDW'

            elif line.startswith('LRC '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('Rcut '):
                split_line = line.split()
                assert split_line[1] == '10'

            elif line.startswith('RcutLow '):
                split_line = line.split()
                assert split_line[1] == '1'

            elif line.startswith('Exclude '):
                split_line = line.split()
                assert split_line[1] == '1-3'

            elif line.startswith('Ewald '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('ElectroStatic '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('CachedFourier '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('Tolerance '):
                split_line = line.split()
                assert split_line[1] == '0.000010000000'

            elif line.startswith('1-4scaling '):
                split_line = line.split()
                assert split_line[1] == '0.5'

            elif line.startswith('PressureCalc '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('RunSteps '):
                split_line = line.split()
                assert split_line[1] == '10'

            elif line.startswith('EqSteps '):
                split_line = line.split()
                assert split_line[1] == '1'

            elif line.startswith('AdjSteps '):
                split_line = line.split()
                assert split_line[1] == '1'

            elif line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.3'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.3'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('CellBasisVector1 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '1.0'
                assert split_line[3] == '0.00'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector2 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '1.0'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector3 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '0.00'
                assert split_line[4] == '1.0'

            elif line.startswith('CBMC_First '):
                split_line = line.split()
                assert split_line[1] == '12'

            elif line.startswith('CBMC_Nth'):
                split_line = line.split()
                assert split_line[1] == '10'

            elif line.startswith('CBMC_Ang '):
                split_line = line.split()
                assert split_line[1] == '50'

            elif line.startswith('CBMC_Dih '):
                split_line = line.split()
                assert split_line[1] == '50'

            elif line.startswith('OutputName '):
                split_line = line.split()
                assert split_line[1] == 'Output_data'

            elif line.startswith('RestartFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('CheckpointFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('CoordinatesFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('ConsoleFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('BlockAverageFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('HistogramFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '1'

            elif line.startswith('DistName '):
                split_line = line.split()
                assert split_line[1] == 'dis'

            elif line.startswith('HistName '):
                split_line = line.split()
                assert split_line[1] == 'his'

            elif line.startswith('RunNumber '):
                split_line = line.split()
                assert split_line[1] == '1'

            elif line.startswith('RunLetter '):
                split_line = line.split()
                assert split_line[1] == 'a'

            elif line.startswith('SampleFreq '):
                split_line = line.split()
                assert split_line[1] == '1'

            elif line.startswith('OutEnergy '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            elif line.startswith('OutPressure '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            elif line.startswith('OutMolNumber '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            elif line.startswith('OutDensity '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            elif line.startswith('OutVolume '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            elif line.startswith('OutSurfaceTension '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            else:
                pass


    def test_save_basic_NPT(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'ethane', FF_filename='ethane',
                        residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                        box_0=[2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_NPT.conf', 'NPT', 1000, 500)

        out_GOMC = open('test_save_basic_NPT.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):

            if line.startswith('Pressure '):
                split_line = line.split()
                assert split_line[1] == '1.01325'

            elif line.startswith('Temperature '):
                split_line = line.split()
                assert split_line[1] == '500'

            elif line.startswith('PressureCalc '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('RunSteps '):
                split_line = line.split()
                assert split_line[1] == '1000'

            elif line.startswith('EqSteps '):
                split_line = line.split()
                assert split_line[1] == '100'

            elif line.startswith('AdjSteps '):
                split_line = line.split()
                assert split_line[1] == '100'

            elif line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.29'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.3'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.01'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('CellBasisVector1 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '2.0'
                assert split_line[3] == '0.00'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector2 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '2.0'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector3 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '0.00'
                assert split_line[4] == '2.0'

            elif line.startswith('RestartFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('CheckpointFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('CoordinatesFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('ConsoleFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('BlockAverageFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('HistogramFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '100'

            elif line.startswith('SampleFreq '):
                split_line = line.split()
                assert split_line[1] == '100'

            else:
                pass


    def test_save_basic_GCMC(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'ethane_box_0',
                        structure_box_1= EthaneGOMC, filename_box_1= 'ethane_box_1',
                        FF_filename='ethane_FF',
                        residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                        box_0= [2, 2, 2], box_1 = [2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GCMC.conf', 'GCMC', 100000, 500,
                                             input_variables_dict={"ChemPot" : {'ETH': -4000}
                                                                   }
                                             )

        out_GOMC = open('test_save_basic_GCMC.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):

            if line.startswith('Parameters '):
                split_line = line.split()
                assert split_line[1] == 'ethane_FF.inp'

            elif line.startswith('Coordinates 0'):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane_box_0.pdb'

            elif line.startswith('Coordinates 1'):
                split_line = line.split()
                assert split_line[1] == '1'
                assert split_line[2] == 'ethane_box_1.pdb'

            elif line.startswith('Structure 0'):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane_box_0.psf'

            elif line.startswith('Structure 1'):
                split_line = line.split()
                assert split_line[1] == '1'
                assert split_line[2] == 'ethane_box_1.psf'

            elif line.startswith('Temperature '):
                split_line = line.split()
                assert split_line[1] == '500'

            elif line.startswith('ChemPot '):
                split_line = line.split()
                assert split_line[1] == 'ETH'
                assert split_line[2] == '-4000'

            elif line.startswith('PressureCalc '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('RunSteps '):
                split_line = line.split()
                assert split_line[1] == '100000'

            elif line.startswith('EqSteps '):
                split_line = line.split()
                assert split_line[1] == '10000'

            elif line.startswith('AdjSteps '):
                split_line = line.split()
                assert split_line[1] == '1000'

            elif line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.35'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.15'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('CellBasisVector1 0'):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '2.0'
                assert split_line[3] == '0.00'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector2 0'):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '2.0'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector3 0'):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '0.00'
                assert split_line[4] == '2.0'

            elif line.startswith('CellBasisVector1 1'):
                split_line = line.split()
                assert split_line[1] == '1'
                assert split_line[2] == '2.0'
                assert split_line[3] == '0.00'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector2 1'):
                split_line = line.split()
                assert split_line[1] == '1'
                assert split_line[2] == '0.00'
                assert split_line[3] == '2.0'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector3 1'):
                split_line = line.split()
                assert split_line[1] == '1'
                assert split_line[2] == '0.00'
                assert split_line[3] == '0.00'
                assert split_line[4] == '2.0'

            elif line.startswith('RestartFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('CheckpointFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('CoordinatesFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('ConsoleFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('BlockAverageFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('HistogramFreq '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '10000'

            elif line.startswith('SampleFreq '):
                split_line = line.split()
                assert split_line[1] == '500'

            else:
                pass


    def test_save_basic_GEMC_NVT(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'ethane_box_0',
                        structure_box_1= EthaneGOMC, filename_box_1= 'ethane_box_1',
                        FF_filename='ethane_FF',
                        residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                        box_0= [2, 2, 2], box_1 = [2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GEMC_NVT.conf', 'GEMC_NVT', 1000000, 500,
                                             )

        out_GOMC = open('test_save_basic_GEMC_NVT.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):

            if line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'


    def test_save_basic_GEMC_NPT(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'ethane_box_0',
                        structure_box_1= EthaneGOMC, filename_box_1= 'ethane_box_1',
                        FF_filename='ethane_FF',
                        residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                        box_0= [2, 2, 2], box_1 = [2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GEMC_NPT.conf', 'GEMC_NPT', 1000000, 500,
                                             input_variables_dict={"Pressure": 10
                                                                   }
                                             )

        out_GOMC = open('test_save_basic_GEMC_NPT.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if line.startswith('Pressure '):
                split_line = line.split()
                assert split_line[1] == '10'

            elif line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.19'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.01'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'


    def test_save_change_most_variable_NVT(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])
        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol',
                        FF_filename='ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name],
                        forcefield_selection='oplsaa'
                        )

        gomc_control.write_gomc_control_file(charmm, 'test_save_change_most_variable_NVT.conf', 'NVT', 100000, 300,
                                             input_variables_dict={'Restart': True, 'PRNG': 123,
                                                                   'ParaTypeCHARMM': True,
                                                                   'ParaTypeMARTINI': False,
                                                                   'ParaTypeMie': False,
                                                                   'LRC' : False,
                                                                   'Rcut': 12,
                                                                   'RcutLow': 8,
                                                                   'Exclude' : '1-4',
                                                                   'Ewald' : False,
                                                                   'ElectroStatic' : False,
                                                                   'CachedFourier' : True,
                                                                   "RcutCoulomb_box_0": 14,
                                                                   "PressureCalc": [False, 4],
                                                                   "Tolerance": 0.01,
                                                                   "DisFreq": 0.2,
                                                                   "RotFreq": 0.2,
                                                                   "IntraSwapFreq": 0.1,
                                                                   "RegrowthFreq": 0.1,
                                                                   "CrankShaftFreq": 0.2,
                                                                   "MultiParticleFreq": 0.05,
                                                                   'IntraMEMC-1Freq': 0.05,
                                                                   'IntraMEMC-2Freq': 0.05,
                                                                   'IntraMEMC-3Freq': 0.05,
                                                                   "CBMC_First": 55,
                                                                   "CBMC_Nth": 66,
                                                                   "CBMC_Ang": 33,
                                                                   "CBMC_Dih": 22,
                                                                   "OutputName": 'test_out',
                                                                   "RestartFreq": [False, 50],
                                                                   "CheckpointFreq": [False, 50],
                                                                   "CoordinatesFreq": [False, 50],
                                                                   "ConsoleFreq": [False, 500],
                                                                   "BlockAverageFreq": [False, 50],
                                                                   "HistogramFreq": [False, 50],
                                                                   "DistName": 'dist',
                                                                   "HistName": 'hist',
                                                                   "RunNumber": 4,
                                                                   "RunLetter": 'c',
                                                                   "SampleFreq": 25,
                                                                   "FreeEnergyCalc": [True, 50],
                                                                   "MoleculeType": ['ETH', 1],
                                                                   "InitialState": 3,
                                                                   "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                   "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                   "MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]],
                                                                   "OutEnergy" : [False, False],
                                                                   "OutPressure": [False, False],
                                                                   "OutMolNumber": [False, False],
                                                                   "OutDensity": [False, False],
                                                                   "OutVolume": [False, False],
                                                                   "OutSurfaceTension": [True, True],
                                                                   }
                                             )

        out_GOMC = open('test_save_change_most_variable_NVT.conf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if line.startswith('Restart '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('PRNG '):
                split_line = line.split()
                assert split_line[1] == 'INTSEED'

            elif line.startswith('Random_Seed '):
                split_line = line.split()
                assert split_line[1] == '123'

            elif line.startswith('ParaTypeCHARMM '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('Parameters '):
                split_line = line.split()
                assert split_line[1] == 'ethane_ethanol.inp'

            elif line.startswith('Coordinates '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane_ethanol.pdb'

            elif line.startswith('Structure '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == 'ethane_ethanol.psf'

            elif line.startswith('Temperature '):
                split_line = line.split()
                assert split_line[1] == '300'

            elif line.startswith('Potential '):
                split_line = line.split()
                assert split_line[1] == 'VDW'

            elif line.startswith('LRC '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('Rcut '):
                split_line = line.split()
                assert split_line[1] == '12'

            elif line.startswith('RcutLow '):
                split_line = line.split()
                assert split_line[1] == '8'

            elif line.startswith('Exclude '):
                split_line = line.split()
                assert split_line[1] == '1-4'

            elif line.startswith('Ewald '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('ElectroStatic '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('CachedFourier '):
                split_line = line.split()
                assert split_line[1] == 'True'

            elif line.startswith('Tolerance '):
                split_line = line.split()
                assert split_line[1] == '0.010000000000'

            elif line.startswith('1-4scaling '):
                split_line = line.split()
                assert split_line[1] == '0.5'

            elif line.startswith('RcutCoulomb 0 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '14'

            elif line.startswith('PressureCalc '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '4'

            elif line.startswith('RunSteps '):
                split_line = line.split()
                assert split_line[1] == '100000'

            elif line.startswith('EqSteps '):
                split_line = line.split()
                assert split_line[1] == '10000'

            elif line.startswith('AdjSteps '):
                split_line = line.split()
                assert split_line[1] == '1000'

            elif line.startswith('DisFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('RotFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('IntraSwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('SwapFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('RegrowthFreq '):
                split_line = line.split()
                assert split_line[1] == '0.1'

            elif line.startswith('CrankShaftFreq '):
                split_line = line.split()
                assert split_line[1] == '0.2'

            elif line.startswith('VolFreq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('MultiParticleFreq '):
                split_line = line.split()
                assert split_line[1] == '0.05'

            elif line.startswith('IntraMEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.05'

            elif line.startswith('MEMC-1Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.05'

            elif line.startswith('MEMC-2Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('IntraMEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.05'

            elif line.startswith('MEMC-3Freq '):
                split_line = line.split()
                assert split_line[1] == '0.0'

            elif line.startswith('CellBasisVector1 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '4.0'
                assert split_line[3] == '0.00'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector2 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '4.0'
                assert split_line[4] == '0.00'

            elif line.startswith('CellBasisVector3 '):
                split_line = line.split()
                assert split_line[1] == '0'
                assert split_line[2] == '0.00'
                assert split_line[3] == '0.00'
                assert split_line[4] == '4.0'

            elif line.startswith('FreeEnergyCalc '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == '50'

            elif line.startswith('MoleculeType '):
                split_line = line.split()
                assert split_line[1] == 'ETH'
                assert split_line[2] == '1'

            elif line.startswith('InitialState '):
                split_line = line.split()
                assert split_line[1] == '3'

            elif line.startswith('ScalePower '):
                split_line = line.split()
                assert split_line[1] == '2'

            elif line.startswith('ScaleAlpha '):
                split_line = line.split()
                assert split_line[1] == '0.5'

            elif line.startswith('MinSigma '):
                split_line = line.split()
                assert split_line[1] == '3'

            elif line.startswith('ScaleCoulomb '):
                split_line = line.split()
                assert split_line[1] == 'False'

            elif line.startswith('# States '):
                split_line = line.split()
                assert split_line[2] == '0'
                assert split_line[3] == '1'
                assert split_line[4] == '2'
                assert split_line[5] == '3'

            elif line.startswith('LambdaVDW '):
                split_line = line.split()
                assert split_line[1] == '0.1'
                assert split_line[2] == '0.2'
                assert split_line[3] == '0.4'
                assert split_line[4] == '0.9'

            elif line.startswith('LambdaCoulomb '):
                split_line = line.split()
                assert split_line[1] == '0.1'
                assert split_line[2] == '0.3'
                assert split_line[3] == '0.8'
                assert split_line[4] == '0.8'

            elif line.startswith('CBMC_First '):
                split_line = line.split()
                assert split_line[1] == '55'

            elif line.startswith('CBMC_Nth '):
                split_line = line.split()
                assert split_line[1] == '66'

            elif line.startswith('CBMC_Ang '):
                split_line = line.split()
                assert split_line[1] == '33'

            elif line.startswith('CBMC_Dih '):
                split_line = line.split()
                assert split_line[1] == '22'

            elif line.startswith('OutputName '):
                split_line = line.split()
                assert split_line[1] == 'test_out'

            elif line.startswith('RestartFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '50'

            elif line.startswith('CheckpointFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '50'

            elif line.startswith('CoordinatesFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '50'

            elif line.startswith('ConsoleFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '500'

            elif line.startswith('BlockAverageFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '50'

            elif line.startswith('HistogramFreq '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == '50'

            elif line.startswith('DistName '):
                split_line = line.split()
                assert split_line[1] == 'dist'

            elif line.startswith('HistName '):
                split_line = line.split()
                assert split_line[1] == 'hist'

            elif line.startswith('RunNumber '):
                split_line = line.split()
                assert split_line[1] == '4'

            elif line.startswith('RunLetter '):
                split_line = line.split()
                assert split_line[1] == 'c'

            elif line.startswith('SampleFreq '):
                split_line = line.split()
                assert split_line[1] == '25'

            elif line.startswith('OutEnergy '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            elif line.startswith('OutPressure '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            elif line.startswith('OutMolNumber '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            elif line.startswith('OutDensity '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            elif line.startswith('OutVolume '):
                split_line = line.split()
                assert split_line[1] == 'False'
                assert split_line[2] == 'False'

            elif line.startswith('OutSurfaceTension '):
                split_line = line.split()
                assert split_line[1] == 'True'
                assert split_line[2] == 'True'

            else:
                pass

    def test_save_NVT_bad_variables_part_1(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Restart': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartCheckpoint' : 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PRNG' : [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeCHARMM': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMie': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMARTINI': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_0': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_1': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Pressure': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rcut': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutLow': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LRC': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Exclude': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'coul_1_4_scaling': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Potential': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rswitch': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ElectroStatic': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Ewald': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CachedFourier': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Tolerance': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Dielectric': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'EqSteps': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'AdjSteps': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'useConstantArea': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FixVolBox0': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ChemPot': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Fugacity': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_First': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Nth': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Ang': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Dih': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutputName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CoordinatesFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CheckpointFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ConsoleFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'BlockAverageFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistogramFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DistName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunNumber': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunLetter': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SampleFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutPressure': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutMolNumber': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutDensity': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutVolume': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutSurfaceTension': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FreeEnergyCalc': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MoleculeType': ['s'], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'InitialState': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaVDW': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaCoulomb': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleCoulomb': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScalePower': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleAlpha': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MinSigma': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ExchangeVolumeDim': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC_DataInput': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DisFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RotFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraSwapFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SwapFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RegrowthFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CrankShaftFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'VolFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MultiParticleFreq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-1Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-1Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-2Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-2Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-3Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-3Freq': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'XXXXXX': 's', }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_2(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Restart': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartCheckpoint' : [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PRNG' : [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeCHARMM': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMie': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMARTINI': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_0': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_1': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Pressure': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rcut': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutLow': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LRC': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Exclude': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'coul_1_4_scaling': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Potential': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rswitch': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ElectroStatic': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Ewald': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CachedFourier': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Tolerance': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Dielectric': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'EqSteps': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'AdjSteps': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'useConstantArea': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FixVolBox0': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ChemPot': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Fugacity': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_First': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Nth': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Ang': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Dih': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutputName': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CoordinatesFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CheckpointFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ConsoleFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'BlockAverageFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistogramFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DistName': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistName': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunNumber': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunLetter': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SampleFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutPressure': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutMolNumber': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutDensity': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutVolume': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutSurfaceTension': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FreeEnergyCalc': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MoleculeType': [[]], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'InitialState': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaVDW': [], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaCoulomb': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleCoulomb': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScalePower': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleAlpha': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MinSigma': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ExchangeVolumeDim': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC_DataInput': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DisFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RotFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraSwapFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SwapFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RegrowthFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CrankShaftFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'VolFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MultiParticleFreq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-1Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-1Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-2Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-2Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-3Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-3Freq': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'XXXXXX': [], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_3(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Restart': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartCheckpoint' : {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PRNG' : {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeCHARMM': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMie': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMARTINI': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_0': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_1': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Pressure': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rcut': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutLow': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LRC': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Exclude': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'coul_1_4_scaling': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Potential': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rswitch': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ElectroStatic': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Ewald': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CachedFourier': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Tolerance': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Dielectric': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'EqSteps': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'AdjSteps': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'useConstantArea': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FixVolBox0': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ChemPot': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Fugacity': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_First': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Nth': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Ang': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Dih': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutputName': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CoordinatesFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CheckpointFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ConsoleFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'BlockAverageFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistogramFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DistName': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistName': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunNumber': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunLetter': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SampleFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutPressure': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutMolNumber': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutDensity': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutVolume': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutSurfaceTension': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FreeEnergyCalc': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MoleculeType': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'InitialState': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaVDW': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaCoulomb': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleCoulomb': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScalePower': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleAlpha': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MinSigma': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ExchangeVolumeDim': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC_DataInput': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DisFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RotFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraSwapFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SwapFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RegrowthFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CrankShaftFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'VolFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MultiParticleFreq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-1Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-1Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-2Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-2Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-3Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-3Freq': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_3.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'XXXXXX': {}, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_4(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Restart': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartCheckpoint' : 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PRNG' : [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeCHARMM': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMie': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ParaTypeMARTINI': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_0': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutCoulomb_box_1': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Pressure': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rcut': [3], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RcutLow': 20, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LRC': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Exclude': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'coul_1_4_scaling': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Potential': 1 }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Rswitch': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ElectroStatic': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Ewald': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CachedFourier': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Tolerance': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Dielectric': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'EqSteps': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'AdjSteps': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'useConstantArea': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FixVolBox0': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ChemPot': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Fugacity': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_First': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Nth': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Ang': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CBMC_Dih': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutputName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CoordinatesFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RestartFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CheckpointFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ConsoleFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'BlockAverageFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistogramFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DistName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'HistName': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunNumber': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RunLetter': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SampleFreq': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutPressure': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutMolNumber': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutDensity': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutVolume': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutSurfaceTension': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'FreeEnergyCalc': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MoleculeType': 1, }
                                                         )


        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'InitialState': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaVDW': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'LambdaCoulomb': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleCoulomb': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScalePower': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleAlpha': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MinSigma': [1], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ExchangeVolumeDim': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC_DataInput': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DisFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'DisFreq': 1,
                                                                               'RotFreq': 0.01}
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RotFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraSwapFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'SwapFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'RegrowthFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'CrankShaftFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'VolFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MultiParticleFreq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-1Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-1Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-2Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-2Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'IntraMEMC-3Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'MEMC-3Freq': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_4.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'XXXXXX': 1, }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_5(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [True, 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [False, 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [1, 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [True , 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [False , 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [1 , 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : ['' , 10000], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [['x'] , 10000], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [{'s' : 1} , 10000], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [True , 1.0], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [True , 'x'], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [True , ['x']], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc' : [True , {'s' : 1}], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"



        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [1 , True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_6(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [True, True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [False, True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [False, False], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [True , True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [False , True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [False , False], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy': [1 , True], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : ['' , True], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [['x'] , True], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [{'s' : 1} , True], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [True , 1.0], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [True , 'x'], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [True , ['x']], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'OutEnergy' : [True , {'s' : 1}], }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_7(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', FF_filename = 'ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                        box_0 = [1,1,1]
                        )

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETH', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9, 0.99],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8,
                                                                                                 0.8, 0.99],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETH', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )

        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],

                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],

                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,

                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # this passes as the "LambdaCoulomb" is default set to "LambdaVDW" if not used
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],

                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # starting bad inputs for the Free engergy calcs side from not using all required variables
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [1, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": ['1', 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [['1'], 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [{'a' : '1'}, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],

                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        # starting bad inputs for the Free engergy calcs side from not using all required variables
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 1.0],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, '1'],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, ['1']],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, {'a' : '1'}],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000, 's'],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # start checking the MoleculeType variable for errors
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": [1, 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": [[1], 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": [{'a' : '1'}, 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000, 's'],
                                                                               "MoleculeType": ['ETO', '1'],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000, 's'],
                                                                               "MoleculeType": ['ETO', ['1']],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000, 's'],
                                                                               "MoleculeType": ['ETO', {'a' : '1'}],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # start checking the initial state variable
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 's',
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": ['s'],
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": {'a' : '1'},
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1.0,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # start checking the LamdaVDW variable
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": ["x", 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [[0.1], 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [{'a' : '1'}, 0.2, 0.4],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # start testing the LambdaCoulomb
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": ["x", 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [[0.1], 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4],
                                                                               "LambdaCoulomb": [{'a': '1'}, 0.3, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # different LambdaVDW and LambdaCoulomb list lengths
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [ 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                               "MoleculeType": ['ETO', 1],
                                                                               "InitialState": 1,
                                                                               "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                               "LambdaCoulomb": [0.3, 0.8, 0.8],
                                                                               }
                                                        )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


    def test_save_NVT_bad_variables_part_8(self, EthaneGOMC, EthanolGOMC):
        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol_box_0',
                        structure_box_1 = test_box_ethane_ethanol, filename_box_1 = 'ethane_ethanol_box_1',
                        FF_filename='ethane_ethanol',
                        residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa'
                        )

        test_box_ethane_ethanol = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm_NPT_NVT = Charmm(test_box_ethane_ethanol, 'ethane_ethanol_box_0',
                                FF_filename='ethane_ethanol',
                                residues=[EthaneGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa'
                                )

        # test ExchangeVolumeDim for errors
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"MEMC-1Freq": 1,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"ExchangeVolumeDim": [1.0, 1.0, 1.0],
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"ExchangeVolumeDim": [1, 1, 1] }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"ExchangeVolumeDim": ['s', 1.0, 1.0]  }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"ExchangeVolumeDim": [1.0, [1.0], 1.0] }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"ExchangeVolumeDim": [1.0, 1.0, {1.0}]}
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # testing failures and passes for MEMC_DataInput
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1.0, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ ['s', 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [[1], 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [[{'a' : '1'}, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [[1, 'ETHaa', ['C1', 'C2'], 'ETO', ['C1', 'C2']]],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 1, ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, [1], ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', [1, 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', [[1], 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 1], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', [1]], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 1, ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], [1], ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', [1, 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', [[1], 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 1]]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', [1]]]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # test the MEMC move ratios cant be set without specifying the MEMC move paramters ("MEMC_DataInput")
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={
                                                                               "IntraMEMC-1Freq": 0.20 ,
                                                                               "MEMC-1Freq": 0.20 ,
                                                                               "IntraMEMC-2Freq": 0.20 ,
                                                                               "MEMC-2Freq": 0.20 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        #testing the move frequency sum to 1 for all ensembles
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.10 ,
                                                                               "MEMC-2Freq": 0.10 ,
                                                                               "IntraMEMC-3Freq": 0.10 ,
                                                                               "MEMC-3Freq": 0.10 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                             [ [1, 'ETH', ['C1', 'C2'], 'ETO', ['C1', 'C2']]  ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq":   0.05 ,
                                                                               "CrankShaftFreq":  0.05 ,
                                                                               "VolFreq":   0.05 ,
                                                                               "MultiParticleFreq": 0.05 ,
                                                                               "IntraMEMC-1Freq": 0.10 ,
                                                                               "MEMC-1Freq": 0.10 ,
                                                                               "IntraMEMC-2Freq": 0.20 ,
                                                                               "MEMC-2Freq": 0.20 ,
                                                                               "IntraMEMC-3Freq": 0.20 ,
                                                                               "MEMC-3Freq": 0.20 ,
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # check the  input_variables_dict ChemPot and Fugacity for errors
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "DisFreq": 1,
                                                                               "Fugacity": {1: 0, "ETO": 1.0},
                                                                               }
                                                        )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"


        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                               "Fugacity": {"ETH": -1, "ETO": 1.0},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                               "Fugacity": {"ETH": "1", "ETO": 1.0},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                               "Fugacity": {"ETH": ["1"], "ETO": 1.0},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                               "Fugacity": {"ETH": 0, "ETO": 1.0},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                                "ChemPot": {1: -4000, "ETO": -8000},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                                "ChemPot": {'ETH': '40', "ETO": -8000},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                                "DisFreq": 1,
                                                                                "ChemPot": {'ETH': ['40'], "ETO": -8000},
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # test good values of Volume for NVT, and GCMC if set to zero
        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "ChemPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"



        # test bad values of Volume for NVT, and GCMC
        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "ChemPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # test bad values of MEMC  for NVT, NPT
        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT' , 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value== "TEST_FAILED"

        # test good values of MEMC  with GCMC
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC' , 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "ChemPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "ChemPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'], 'ETO',
                                                                                     ['C1', 'C2']]],
                                                                               "ChemPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"