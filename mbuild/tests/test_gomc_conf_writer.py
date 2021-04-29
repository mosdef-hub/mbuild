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
                         'Rcut', 'RcutLow', 'LRC', 'Exclude', 'Potential', 'Rswitch',
                         'ElectroStatic', 'Ewald', 'CachedFourier', 'Tolerance', 'Dielectric', 'PressureCalc',
                         'EqSteps', 'AdjSteps', 'VDWGeometricSigma', 'useConstantArea', 'FixVolBox0',
                         'ChemPot', 'Fugacity',
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
                                                         'LRC', 'Exclude', 'Potential', 'Rswitch',
                                                         'ElectroStatic', 'Ewald', 'CachedFourier', 'Tolerance',
                                                         'Dielectric', 'PressureCalc', 'EqSteps', 'AdjSteps',
                                                         'VDWGeometricSigma', 'useConstantArea', 'FixVolBox0',
                                                         'ChemPot', 'Fugacity',
                                                         'CBMC_First', 'CBMC_Nth', 'CBMC_Ang', 'CBMC_Dih',
                                                         'OutputName', 'CoordinatesFreq', 'RestartFreq',
                                                         'CheckpointFreq', 'ConsoleFreq', 'BlockAverageFreq',
                                                         'HistogramFreq', 'DistName', 'HistName', 'RunNumber',
                                                         'RunLetter', 'SampleFreq', 'OutEnergy', 'OutPressure',
                                                         'OutMolNumber', 'OutDensity', 'OutVolume',
                                                         'OutSurfaceTension', 'FreeEnergyCalc', 'MoleculeType',
                                                         'InitialState', 'LambdaVDW', 'LambdaCoulomb', 'ScaleCoulomb',
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
                                                         'VDWGeometricSigma', 'useConstantArea', 'FixVolBox0',
                                                         'ChemPot', 'Fugacity',
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
            gomc_control.print_required_input(description=True)
            gomc_control.print_required_input(description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('NVT', description=True)
            gomc_control.print_valid_ensemble_input_variables('NVT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('NPT', description=True)
            gomc_control.print_valid_ensemble_input_variables('NPT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NVT', description=True)
            gomc_control.print_valid_ensemble_input_variables('GEMC_NVT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GEMC_NPT', description=True)
            gomc_control.print_valid_ensemble_input_variables('GEMC_NPT', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

        try:
            gomc_control.print_valid_ensemble_input_variables('GCMC', description=True)
            gomc_control.print_valid_ensemble_input_variables('GCMC', description=False)
            test_status = "PASSED"
        except:
            test_status = "FAILED"
        assert test_status == "PASSED"

    def test_save_basic_NVT(self, ethane_gomc):
        charmm = Charmm(ethane_gomc, 'ethane', ff_filename='ethane',
                        residues=[ethane_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_NVT.conf', 'NVT', 10, 300,
                                             )

        with open('test_save_basic_NVT.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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

                elif line.startswith('VDWGeometricSigma '):
                    split_line = line.split()
                    assert split_line[1] == 'False'

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
                    print('split_line[2] = '+str(split_line[2]))
                    assert split_line[2] == '10.0'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector2 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '10.0'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector3 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '10.0'

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

    def test_save_basic_NPT(self, ethane_gomc):
        charmm = Charmm(ethane_gomc, 'ethane', ff_filename='ethane',
                        residues=[ethane_gomc.name], forcefield_selection='oplsaa',
                        box_0=[2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_NPT.conf', 'NPT', 1000, 500)

        with open('test_save_basic_NPT.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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
                    assert split_line[2] == '20.0'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector2 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '20.0'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector3 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '20.0'

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

                elif line.startswith('VDWGeometricSigma '):
                    split_line = line.split()
                    assert split_line[1] == 'False'

                elif line.startswith('useConstantArea '):
                    split_line = line.split()
                    assert split_line[1] == 'False'

                else:
                    pass

    def test_save_basic_GCMC(self, ethane_gomc):
        charmm = Charmm(ethane_gomc, 'ethane_box_0',
                        structure_box_1=ethane_gomc, filename_box_1='ethane_box_1',
                        ff_filename='ethane_FF',
                        residues=[ethane_gomc.name], forcefield_selection='oplsaa',
                        box_0=[2, 2, 2], box_1=[2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GCMC.conf', 'GCMC', 100000, 500,
                                             input_variables_dict={"ChemPot": {'ETH': -4000},
                                                                   "VDWGeometricSigma": True
                                                                   }
                                             )

        with open('test_save_basic_GCMC.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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
                    assert split_line[2] == '20.0'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector2 0'):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '20.0'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector3 0'):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '20.0'

                elif line.startswith('CellBasisVector1 1'):
                    split_line = line.split()
                    assert split_line[1] == '1'
                    assert split_line[2] == '20.0'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector2 1'):
                    split_line = line.split()
                    assert split_line[1] == '1'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '20.0'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector3 1'):
                    split_line = line.split()
                    assert split_line[1] == '1'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '20.0'

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

                elif line.startswith('VDWGeometricSigma '):
                    split_line = line.split()
                    assert split_line[1] == 'True'

                else:
                    pass

    def test_save_basic_GEMC_NVT(self, ethane_gomc):
        charmm = Charmm(ethane_gomc, 'ethane_box_0',
                        structure_box_1=ethane_gomc, filename_box_1='ethane_box_1',
                        ff_filename='ethane_FF',
                        residues=[ethane_gomc.name], forcefield_selection='oplsaa',
                        box_0=[2, 2, 2], box_1=[2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GEMC_NVT.conf', 'GEMC_NVT', 1000000, 500,
                                             )

        with open('test_save_basic_GEMC_NVT.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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

                else:
                    pass

    def test_save_basic_GEMC_NPT(self, ethane_gomc):
        charmm = Charmm(ethane_gomc, 'ethane_box_0',
                        structure_box_1=ethane_gomc, filename_box_1='ethane_box_1',
                        ff_filename='ethane_FF',
                        residues=[ethane_gomc.name], forcefield_selection='oplsaa',
                        box_0=[2, 2, 2], box_1=[2, 2, 2]
                        )
        gomc_control.write_gomc_control_file(charmm, 'test_save_basic_GEMC_NPT.conf', 'GEMC_NPT', 1000000, 500,
                                             input_variables_dict={"Pressure": 10,
                                                                   "useConstantArea": True,
                                                                   "FixVolBox0": True
                                                                   }
                                             )

        with open('test_save_basic_GEMC_NPT.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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

                elif line.startswith('useConstantArea '):
                    split_line = line.split()
                    assert split_line[1] == 'True'

                elif line.startswith('FixVolBox0 '):
                    split_line = line.split()
                    assert split_line[1] == 'True'

                else:
                    pass

    def test_save_change_most_variable_NVT(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])
        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol',
                        ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name],
                        forcefield_selection='oplsaa'
                        )

        gomc_control.write_gomc_control_file(charmm, 'test_save_change_most_variable_NVT.conf', 'NVT', 100000, 300,
                                             input_variables_dict={'Restart': True, 'PRNG': 123,
                                                                   'ParaTypeCHARMM': True,
                                                                   'ParaTypeMARTINI': False,
                                                                   'ParaTypeMie': False,
                                                                   'LRC': False,
                                                                   'Rcut': 12,
                                                                   'RcutLow': 8,
                                                                   'Exclude': '1-4',
                                                                   'Ewald': False,
                                                                   'ElectroStatic': False,
                                                                   'CachedFourier': True,
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
                                                                   "OutEnergy": [False, False],
                                                                   "OutPressure": [False, False],
                                                                   "OutMolNumber": [False, False],
                                                                   "OutDensity": [False, False],
                                                                   "OutVolume": [False, False],
                                                                   "OutSurfaceTension": [True, True],
                                                                   }
                                             )

        with open('test_save_change_most_variable_NVT.conf', 'r') as fp:
            out_gomc = fp.readlines()
            for i, line in enumerate(out_gomc):
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
                    assert split_line[2] == '40.0'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector2 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '40.0'
                    assert split_line[4] == '0.00'

                elif line.startswith('CellBasisVector3 '):
                    split_line = line.split()
                    assert split_line[1] == '0'
                    assert split_line[2] == '0.00'
                    assert split_line[3] == '0.00'
                    assert split_line[4] == '40.0'

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

    def test_save_NVT_bad_variables_part_1(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
                        )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Restart'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Restart': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RestartCheckpoint'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RestartCheckpoint': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PRNG'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PRNG': [1], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeCHARMM'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeCHARMM': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeMie'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeMie': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeMARTINI'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeMARTINI': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RcutCoulomb_box_0'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutCoulomb_box_0': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['RcutCoulomb_box_1'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutCoulomb_box_1': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Pressure'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Pressure': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Rcut'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Rcut': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RcutLow'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutLow': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LRC'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'LRC': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Exclude'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Exclude': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Potential'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'Potential': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Rswitch'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Rswitch': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ElectroStatic'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ElectroStatic': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Ewald'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Ewald': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CachedFourier'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CachedFourier': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Tolerance'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Tolerance': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Dielectric'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Dielectric': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['EqSteps'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'EqSteps': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['EqSteps'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'EqSteps': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['useConstantArea'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={'useConstantArea': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ChemPot': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Fugacity': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_First'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_First': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Nth'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Nth': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Ang'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Ang': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Dih'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Dih': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutputName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutputName': 1, }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CoordinatesFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CoordinatesFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RestartFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RestartFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CheckpointFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CheckpointFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ConsoleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ConsoleFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['BlockAverageFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'BlockAverageFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['HistogramFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'HistogramFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['DistName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'DistName': 1, }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['HistName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'HistName': 1, }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RunNumber'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RunNumber': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RunLetter'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RunLetter': 1, }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['SampleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'SampleFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutPressure'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutPressure': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutMolNumber'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutMolNumber': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutDensity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutDensity': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutVolume'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutVolume': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutSurfaceTension'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutSurfaceTension': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc":  's',
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 's'],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [['ETH'], 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [{'ETH': "1"}, 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 's',
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaVDW'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": 's',
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": 's'
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the "
                                             r"following variables need to be set, and not equal to "
                                             r"None: FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScaleCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'ScaleCoulomb': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScalePower'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ScalePower': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScaleAlpha'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ScaleAlpha': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MinSigma'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MinSigma': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ExchangeVolumeDim'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ExchangeVolumeDim': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC_DataInput': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['DisFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'DisFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RotFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RotFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraSwapFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraSwapFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['SwapFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'SwapFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RegrowthFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RegrowthFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CrankShaftFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CrankShaftFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['VolFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'VolFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MultiParticleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MultiParticleFreq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-1Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-1Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-1Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-1Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-2Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-2Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-2Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-2Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-3Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-3Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-3Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-3Freq': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['XXXXXX'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'XXXXXX': 's', }
                                                 )

    def test_save_NVT_bad_variables_part_2(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
                        )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Restart'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Restart': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RestartCheckpoint'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RestartCheckpoint': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PRNG'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PRNG': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeCHARMM'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeCHARMM': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeMie'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeMie': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ParaTypeMARTINI'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ParaTypeMARTINI': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RcutCoulomb_box_0'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutCoulomb_box_0': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['RcutCoulomb_box_1'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutCoulomb_box_1': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Pressure'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Pressure': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Rcut'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Rcut': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RcutLow'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RcutLow': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LRC'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'LRC': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Exclude'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Exclude': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Potential'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Potential': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Rswitch'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Rswitch': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ElectroStatic'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ElectroStatic': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Ewald'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Ewald': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CachedFourier'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CachedFourier': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Tolerance'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Tolerance': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Dielectric'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Dielectric': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['EqSteps'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'EqSteps': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['AdjSteps'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'AdjSteps': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['useConstantArea'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'useConstantArea': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['FixVolBox0'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'FixVolBox0': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ChemPot': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'Fugacity': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_First'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_First': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Nth'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Nth': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Ang'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Ang': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CBMC_Dih'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CBMC_Dih': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutputName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutputName': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CoordinatesFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CoordinatesFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RestartFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RestartFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CheckpointFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CheckpointFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ConsoleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ConsoleFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['BlockAverageFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'BlockAverageFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['HistogramFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'HistogramFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['DistName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'DistName': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['HistName'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'HistName': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RunNumber'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RunNumber': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RunLetter'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RunLetter': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['SampleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'SampleFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutPressure'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutPressure': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutMolNumber'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutMolNumber': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutDensity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutDensity': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutVolume'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutVolume': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutSurfaceTension'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutSurfaceTension': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc":  [],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', []],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [['ETH'], 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [{'ETH': "1"}, 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": [],
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaVDW'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETH', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": []}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the "
                                             r"following variables need to be set, and not equal to "
                                             r"None: FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScaleCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ScaleCoulomb': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScalePower'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ScalePower': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ScaleAlpha'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ScaleAlpha': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MinSigma'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MinSigma': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ExchangeVolumeDim'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'ExchangeVolumeDim': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC_DataInput': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['DisFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'DisFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['DisFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'DisFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraSwapFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraSwapFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraSwapFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraSwapFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RegrowthFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'RegrowthFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['CrankShaftFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'CrankShaftFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['VolFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'VolFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MultiParticleFreq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MultiParticleFreq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-1Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-1Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-1Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-1Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-2Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-2Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-2Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-2Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['IntraMEMC-3Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'IntraMEMC-3Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC-3Freq'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'MEMC-3Freq': [], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: All the correct input variables where not provided for "
                                             r"the NVT ensemble. Please be sure to check that the keys in the "
                                             r"input variables dictionary \(input_variables_dict\) is correct, and "
                                             r"be aware that added spaces before or after the variable in any keys "
                                             r"will also give this warning. The bad variable inputs ensemble "
                                             r"inputs = \['XXXXXX'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_2.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'XXXXXX': [], }
                                                 )

    def test_save_NVT_bad_variables_part_5(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
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

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [1, 10000], }
                                                 )

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [True, 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={'PressureCalc': [False, 10000], }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [1, 10000], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': ['', 10000], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [['x'], 10000], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [{'s': 1}, 10000], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [True, 1.0], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [True, 'x'], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [True, ['x']], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [True, {'s': 1}], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['PressureCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_5.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'PressureCalc': [1, True], }
                                                 )

    def test_save_NVT_bad_variables_part_6(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
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
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [False, False], }
                                                 )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [True, True], }
                                                 )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [False, True], }
                                                 )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [False, False], }
                                                 )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [1, True], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': ['', True], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [['x'], True], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [{'s': 1}, True], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [True, 1.0], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [True, 'x'], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [True, ['x']], }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['OutEnergy'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_6.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={'OutEnergy': [True, {'s': 1}], }
                                                 )

    def test_save_NVT_bad_variables_part_7(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol', ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                        box_0=[1, 1, 1]
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

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the following "
                                             r"variables need to be set, and not equal to None: FreeEnergyCalc, "
                                             r"MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={
                                                     "MoleculeType": ['ETO', 1],
                                                     "InitialState": 1,
                                                     "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                     "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                 }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the following "
                                             r"variables need to be set, and not equal to None: FreeEnergyCalc, "
                                             r"MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [False, 10000],

                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the following "
                                             r"variables need to be set, and not equal to None: FreeEnergyCalc, "
                                             r"MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                       "MoleculeType": ['ETO', 1],

                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: To utilize the free energy calculations all the following "
                                             r"variables need to be set, and not equal to None: FreeEnergyCalc, "
                                             r"MoleculeType, InitialState, LambdaVDW."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [False, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,

                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        try:
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
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
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [1, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": ['1', 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [['1'], 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [{'a': '1'}, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

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
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 1.0],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, '1'],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, ['1']],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, {'a': '1'}],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FreeEnergyCalc'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000, 's'],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        # start checking the MoleculeType variable for errors
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [1, 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [[1], 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": [{'a': '1'}, 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', '1'],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', ['1']],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', {'a': '1'}],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MoleculeType'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETOa', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        # start checking the initial state variable
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 's',
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": ['s'],
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": {'a': '1'},
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['InitialState'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1.0,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        # start checking the LamdaVDW variable
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaVDW'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": ["x", 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaVDW'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [[0.1], 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaVDW'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [{'a': '1'}, 0.2, 0.4],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8],
                                                                       }
                                                 )

        # start testing the LambdaCoulomb
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": ["x", 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [[0.1], 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['LambdaCoulomb'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4],
                                                                       "LambdaCoulomb": [{'a': '1'}, 0.3, 0.8],
                                                                       }
                                                 )

        # different LambdaVDW and LambdaCoulomb list lengths
        with pytest.raises(ValueError, match=r"ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_7.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"FreeEnergyCalc": [True, 10000],
                                                                       "MoleculeType": ['ETO', 1],
                                                                       "InitialState": 1,
                                                                       "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                                       "LambdaCoulomb": [0.3, 0.8, 0.8],
                                                                       }
                                                 )

    def test_save_NVT_bad_variables_part_8(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm = Charmm(test_box_ethane_ethanol, 'ethane_ethanol_box_0',
                        structure_box_1=test_box_ethane_ethanol, filename_box_1='ethane_ethanol_box_1',
                        ff_filename='ethane_ethanol',
                        residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa'
                        )

        test_box_ethane_ethanol = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                              n_compounds=[1, 1],
                                              box=[4.0, 4.0, 4.0])

        charmm_NPT_NVT = Charmm(test_box_ethane_ethanol, 'ethane_ethanol_box_0',
                                ff_filename='ethane_ethanol',
                                residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa'
                                )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['RcutCoulomb_box_1'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={'RcutCoulomb_box_1': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['FixVolBox0'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={'FixVolBox0': 's', }
                                                 )

        # test ExchangeVolumeDim for errors
        with pytest.raises(ValueError, match=r"The MEMC_DataInput variable is equal to None, but at least one "
                                             r"of the MEMC move ratios are all non-zero \(IntraMEMC_1Freq, "
                                             r"MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, IntraMEMC_3Freq, "
                                             r"and MEMC_3Freq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={"MEMC-1Freq": 1,
                                                                       }
                                                 )

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
                                                         input_variables_dict={"ExchangeVolumeDim": [1, 1, 1]}
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ExchangeVolumeDim'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={"ExchangeVolumeDim": ['s', 1.0, 1.0]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ExchangeVolumeDim'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={"ExchangeVolumeDim": [1.0, [1.0], 1.0]}
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ExchangeVolumeDim'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={"ExchangeVolumeDim": [1.0, 1.0, {'a': 1.0}]}
                                                 )

        # testing failures and passes for MEMC_DataInput
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
                                                                                ],
                                                                               "DisFreq": 0.05,
                                                                               "RotFreq": 0.05,
                                                                               "IntraSwapFreq":  0.05,
                                                                               "SwapFreq": 0.05,
                                                                               "RegrowthFreq": 0.05,
                                                                               "CrankShaftFreq":  0.05,
                                                                               "VolFreq": 0.05,
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
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'O1']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C2', 'C1'],
                                                                                 'ETO', ['O1', 'C1']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'O1'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['O1', 'C1'],
                                                                         'ETO', ['C2', 'C1']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1.0, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [['s', 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[[1], 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[{'a': '1'}, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETHa', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 1, ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, [1], ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', [1, 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', [[1], 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 1],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', [1]],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         1, ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         [1], ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', [1, 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', [[1], 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 1]]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['MEMC_DataInput'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', [1]]]
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
                                                                       }
                                                 )

        # test the MEMC move ratios cant be set without specifying the MEMC move paramters ("MEMC_DataInput")
        with pytest.raises(ValueError, match=r"ERROR: The MEMC_DataInput variable is equal to None, but at least "
                                             r"one of the MEMC move ratios are all non-zero "
                                             r"\(IntraMEMC_1Freq, MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, "
                                             r"IntraMEMC_3Freq, and MEMC_3Freq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={
                                                     "IntraMEMC-1Freq": 0.20,
                                                     "MEMC-1Freq": 0.20,
                                                     "IntraMEMC-2Freq": 0.20,
                                                     "MEMC-2Freq": 0.20,
                                                     "IntraMEMC-3Freq": 0.10,
                                                     "MEMC-3Freq": 0.10,
                                                 }
                                                 )

        # test some GCMC variable errors with Chempot and fugacity
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={'ChemPot': 's', }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_1.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={'Fugacity': 's', }
                                                 )

        # testing the move frequency sum to 1 for all ensembles
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
                                             r"Note: The sum that was manually entered may equal 1, but some "
                                             r"moves may not be valid for the provided ensemble. The moves that "
                                             r"are invalid for a given ensemble are set to zero. If the default "
                                             r"moves are not being used, all the move frequencies which do not have "
                                             r"default values of zero will need to be set manually so the sum equals "
                                             r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
                                             r"CrankShaftFreq, and VolFreq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GEMC_NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
                                             r"Note: The sum that was manually entered may equal 1, but some "
                                             r"moves may not be valid for the provided ensemble. The moves that "
                                             r"are invalid for a given ensemble are set to zero. If the default "
                                             r"moves are not being used, all the move frequencies which do not have "
                                             r"default values of zero will need to be set manually so the sum equals "
                                             r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
                                             r"CrankShaftFreq, and VolFreq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GEMC_NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        try:
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       "ChemPot": {'ETH': -4000, 'ETO': 8000},
                                                                       }
                                                 )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
                                             r"Note: The sum that was manually entered may equal 1, but some "
                                             r"moves may not be valid for the provided ensemble. The moves that "
                                             r"are invalid for a given ensemble are set to zero. If the default "
                                             r"moves are not being used, all the move frequencies which do not have "
                                             r"default values of zero will need to be set manually so the sum equals "
                                             r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
                                             r"CrankShaftFreq, and VolFreq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       "Fugacity": {'ETH': 0, 'ETO': 1.0},
                                                                       }
                                                 )

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
                                             r"Note: The sum that was manually entered may equal 1, but some "
                                             r"moves may not be valid for the provided ensemble. The moves that "
                                             r"are invalid for a given ensemble are set to zero. If the default "
                                             r"moves are not being used, all the move frequencies which do not have "
                                             r"default values of zero will need to be set manually so the sum equals "
                                             r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
                                             r"CrankShaftFreq, and VolFreq\)."):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        with pytest.raises(ValueError, match=r"ERROR: The sum of the Monte Carlo move ratios does not equal 1. "
                                             r"Note: The sum that was manually entered may equal 1, but some "
                                             r"moves may not be valid for the provided ensemble. The moves that "
                                             r"are invalid for a given ensemble are set to zero. If the default "
                                             r"moves are not being used, all the move frequencies which do not have "
                                             r"default values of zero will need to be set manually so the sum equals "
                                             r"\(DisFreq, RotFreq, IntraSwapFreq, SwapFreq, RegrowthFreq, "
                                             r"CrankShaftFreq, and VolFreq\)."):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        # test good values of Volume for NVT, and GCMC if set to zero
        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NVT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'NPT', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"

        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
                                                                                ],
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

        # test come MEMC with GCMC
        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
                                                                        ],
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {1: 0, "ETO": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
                                                                        ],
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {"ETH": -1, "ETO": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
                                                                        ],
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {"ETH": "1", "ETO": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The MEMC_DataInput variable is not equal to None, "
                                             r"but all the MEMC move ratios are zero \(IntraMEMC_1Freq, MEMC_1Freq, "
                                             r"IntraMEMC_2Freq, MEMC_2Freq, IntraMEMC_3Freq, and MEMC_3Freq\)."):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
                                                                        ],
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {"ETH": 2, "ETO": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {"ETH": 0, "XXX": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['Fugacity'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "Fugacity": {"XXX": 0, "ETO": 1.0},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "ChemPot": {1: -4000, "ETO": -8000},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "ChemPot": {"XXX": -4000, "ETO": -8000},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "ChemPot": {"ETH": -4000, "XXX": -8000},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={
                                                                       "DisFreq": 1,
                                                                       "ChemPot": {'ETH': '40', "ETO": -8000},
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r"ERROR: The following input variables have "
                                             r"bad values \(check spelling and for empty spaces in the keys or that "
                                             r"the values are in the correct form with the acceptable values\)"
                                             r": \['ChemPot'\]"):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"DisFreq": 1,
                                                                       "ChemPot": {'ETH': ['40'], "ETO": -8000},
                                                                       }
                                                 )

        # test bad values of Volume for NVT, and GCMC
        with pytest.raises(ValueError, match=r'ERROR: The input variable VolFreq is non-zero \(0\). '
                                             r'VolFreq must be zero \(0\) for the "NVT", "GEMC_NVT", '
                                             r'and "GCMC" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: The input variable VolFreq is non-zero \(0\). '
                                             r'VolFreq must be zero \(0\) for the "NVT", "GEMC_NVT", '
                                             r'and "GCMC" ensembles.'):
            gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'GCMC', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
                                                                        ],
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

        # test bad values of MEMC  for NVT, NPT
        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NVT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        with pytest.raises(ValueError, match=r'ERROR: All the MC move input variables must be non-zero \(0\) '
                                             r'for the SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. '
                                             r'The SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq need to be zero '
                                             r'\(0\) for the "NVT" and "NPT" ensembles.'):
            gomc_control.write_gomc_control_file(charmm_NPT_NVT, 'test_save_NVT_bad_variables_part_8.conf',
                                                 'NPT', 10, 300,
                                                 input_variables_dict={"MEMC_DataInput":
                                                                       [[1, 'ETH', ['C1', 'C2'],
                                                                         'ETO', ['C1', 'C2']]
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
                                                                       }
                                                 )

        # test good values of MEMC  with GCMC
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC' , 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
                                                                                ],
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
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
                                                                                ],
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
                                                                               [[1, 'ETH', ['C1', 'C2'],
                                                                                 'ETO', ['C1', 'C2']]
                                                                                ],
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

        # try all case unspecific values
        try:
            value = gomc_control.write_gomc_control_file(charmm, 'test_save_NVT_bad_variables_part_8.conf',
                                                         'GCMC', 10, 300,
                                                         input_variables_dict={"MEMC_DataInput":
                                                                                   [[1, 'ETH', ['C1', 'C2'],
                                                                                     'ETO', ['C1', 'C2']]
                                                                                    ],
                                                                               "ChEmPot": {'ETH': -4000, "ETO": -8000},
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
                                                                               }
                                                         )
        except:
            value = "TEST_FAILED"

        assert value == "GOMC_CONTROL_FILE_WRITTEN"


    def test_charmm_object_has_proper_no_boxes_for_ensemble_part_9(self, ethane_gomc, ethanol_gomc):
        test_box_ethane_ethanol_liq = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                                  n_compounds=[4, 4],
                                                  box=[4.0, 4.0, 4.0])

        test_box_ethane_ethanol_vap = mb.fill_box(compound=[ethane_gomc, ethanol_gomc],
                                                  n_compounds=[1, 1],
                                                  box=[8.0, 8.0, 8.0])

        charmm_one_box = Charmm(test_box_ethane_ethanol_liq, 'ethane_ethanol_1_box_liq', ff_filename='ethane_ethanol',
                                residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                                )

        charmm_two_boxes = Charmm(test_box_ethane_ethanol_liq, 'ethane_ethanol_2_boxes_liq',
                                  structure_box_1=test_box_ethane_ethanol_vap, filename_box_1='ethane_box_2_boxes_vap',
                                  ff_filename='ethane_ethanol',
                                  residues=[ethane_gomc.name, ethanol_gomc.name], forcefield_selection='oplsaa',
                                  )

        # test that it fails with the GEMC_NVT with only 1 box in the Charmm object
        with pytest.raises(ValueError, match=r"ERROR: The ensemble type selection of {} is using a Charmm " 
                                             r"object with one simulation boxes, and the {} ensemble only accepts "
                                             r"two boxes \(box 0 and box 1\).".format('GEMC_NVT', 'GEMC_NVT')
                           ):

            gomc_control.write_gomc_control_file(charmm_one_box,
                                                 'test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box',
                                                 'GEMC_NVT', 100, 300
                                                 )

        # test that it fails with the GEMC_NPT with only 1 box in the Charmm object
        with pytest.raises(ValueError, match=r"ERROR: The ensemble type selection of {} is using a Charmm " 
                                             r"object with one simulation boxes, and the {} ensemble only accepts "
                                             r"two boxes \(box 0 and box 1\).".format('GEMC_NPT', 'GEMC_NPT')
                           ):

            gomc_control.write_gomc_control_file(charmm_one_box,
                                                 'test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box',
                                                 'GEMC_NPT', 100, 300
                                                 )

        # test that it fails with the GCMC with only 1 box in the Charmm object
        with pytest.raises(ValueError, match=r"ERROR: The ensemble type selection of {} is using a Charmm " 
                                             r"object with one simulation boxes, and the {} ensemble only accepts "
                                             r"two boxes \(box 0 and box 1\).".format('GCMC', 'GCMC')
                           ):

            gomc_control.write_gomc_control_file(charmm_one_box,
                                                 'test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box',
                                                 'GCMC', 100, 300
                                                 )

        # test that it fails with the NVT with 2 boxes in the Charmm object
        with pytest.raises(ValueError, match=r"ERROR: The ensemble type selection of {} is using a Charmm "
                                             r"object with two simulation boxes, and the {} ensemble only accepts " 
                                             r"one box \(box 0\).".format('NVT', 'NVT')
                           ):

            gomc_control.write_gomc_control_file(charmm_two_boxes,
                                                 'test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box',
                                                 'NVT', 100, 300
                                                 )

        # test that it fails with the NPT with 2 boxes in the Charmm object
        with pytest.raises(ValueError, match=r"ERROR: The ensemble type selection of {} is using a Charmm "
                                             r"object with two simulation boxes, and the {} ensemble only accepts "
                                             r"one box \(box 0\).".format('NPT', 'NPT')
                           ):
            gomc_control.write_gomc_control_file(charmm_two_boxes,
                                                 'test_charmm_object_has_proper_no_boxes_for_ensemble_part_9_1_box',
                                                 'NPT', 100, 300
                                                 )