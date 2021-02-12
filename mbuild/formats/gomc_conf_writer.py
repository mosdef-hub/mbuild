import datetime
import os
import mbuild.formats.charmm_writer as mf_charmm
from warnings import warn


def dict_keys_to_list(dict):
    """
    Parameters
    ----------
    dict : dict, a provided dictionary

    Outputs
    ---------
    list : list of keys from the provided dictionary.
    """

    list = []
    for key in dict.keys():
        list.append(key)
    return list



def print_valid_required_input_variables(description=False):
    """
    Parameters
    ----------
    description =  bool, default = False
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Outputs
    ---------
    Prints out the valid input variables (user optional) on the screen
        , which can be entered in the GOMC writer. These are the valid input
        variables for all ensembles.
    """

    valid_args = _get_all_possible_input_variables(description=description)
    for arg, description in valid_args.items():
        print("{:10s}:    {}".format(arg, description))


def _get_required_data(description=False):
    """
    Provides a list of the required inputs for all possible ensembles.

    Parameters
    ----------
    description =  bool, default = False.
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Outputs
    ---------
    required_data = dict or list, default = list.
        If the description = True then a dict is provided with the key and value.
        if the description = False then a list of the dict keys is provided.
    """


    required_data = {
        "charmm_object": 'Charmm object; ' \
                      + 'A Charmm object, which by definition has been parameterized ' \
                      + 'from the selected force field.',
        "ensemble_type": "Required files or System Info (all ensembles): str " \
                         + "(valid strings are 'NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', or 'GCMC'); " \
                      + 'the ensemble type for the simulation.',
        "RunSteps" : "Required files or System Info (all ensembles): int (> 0); " \
                      + "The number or run steps for the simulation.",
        "Temperature" :  "Required files or System Info (all ensembles): float or integer (> 0); " \
                      + "Temperature of system in Kelvin (K)",

    }

    if description:
        return required_data
    else:
        return list(required_data.keys())


def _get_all_possible_input_variables(description=False):
    """
    Provides a list of the variables inputs (user optional) for all possible ensembles.

    Parameters
    ----------
    description =  bool, default = False.
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Outputs
    ---------
    valid_input_variables = dict or list, default = list.
        If the description = True then a dict is provided with the key and value.
        if the description = False then a list of the dict keys is provided.
    """

    valid_input_variables = {

        # ******************************************************************************************************
        # Definitions in this function are copied to a large extent from the GOMC manual release version 2.50 (start)
        # insert citation here:
        # ******************************************************************************************************
        "Restart": 'Simulation info (all ensembles): boolean; default = False. ' \
                      + 'Determines whether to restart the simulation ' \
                      + 'from restart file (*_restart.pdb and *_restart.psf) or not.',
        "RestartCheckpoint": 'Simulation info (all ensembles): boolean; default = False. ' \
                      + 'Determines whether to restart the ' \
                      + 'simulation with the checkpoint file (checkpoint.dat) or not. Restarting the ' \
                      + 'simulation with checkpoint.dat would result in an identitcal outcome, as if ' \
                      + 'previous simulation was continued.' ,
        "PRNG" : 'Simulation info (all ensembles): string or int (>= 0) ("RANDOM" or integer); default = "RANDOM". ' \
                      + 'Note PRNG = Pseudo-Random Number Generator (PRNG). ' \
                      + 'The first options are to enter the string \n' \
                      + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "RANDOM", which selects a random seed number. ' \
                      + 'This will enter the line "PRNG RANDOM" in the gomc configuration file. \n'\
                      + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- integer, which defines the integer seed number ' \
                      + 'for the simulation. ' \
                      + 'This is equivelent to entering the following two lines in the configuration file: ' \
                      + 'line 1 = PRNG INTSEED, ' \
                      + 'line 2 = Random_Seed user_selected_integer. ' \
                      + 'Example 1: for Random enter the string "RANDOM. ' \
                      + 'Example 2: for a specific seed number enter a integer of your choosing. ',
        "ParaTypeCHARMM": 'Simulation info (all ensembles): boolean; default = True. ' \
                      + 'True if a CHARMM forcefield, False otherwise.',
        "ParaTypeMie": 'Simulation info (all ensembles): boolean; default = False. ' \
                      + 'True if a Mie forcefield type, False otherwise.',
        "ParaTypeMARTINI": 'Simulation info (all ensembles): boolean; default = False. ' \
                      + 'True if a MARTINI forcefield, False otherwise.',

        "RcutCoulomb_box_0": 'Simulation info (all ensembles): int or float (>= 0); '\
                             + 'default = None (Note: if None, GOMC will default to the Rcut value). '\
                             + 'Sets a specific radius in box 0 where the short range ' \
                             + 'electrostatic energy will be calculated (i.e., The distance to truncate the ' \
                             + 'short range electrostatic energy in box 0.)',
        "RcutCoulomb_box_1": 'Simulation info (only GEMC_NPT, GEMC_NVT, and GCMC): int, or float (>= 0); ' \
                             + 'default = None (Note: if None, GOMC will default to the Rcut value). '\
                             + 'Sets a specific radius in box 1 where the short range  ' \
                             + 'electrostatic energy will be calculated. (i.e., The distance to truncate the ' \
                             + 'short range electrostatic energy in box 1.)',
        "Pressure": 'Simulation info (only GEMC_NPT and NPT): int (>= 0); default = 1.01325. ' \
                    + 'The pressure in bar utilized for the NPT ' \
                    + 'and GEMC_NPT simulations.',
        "Rcut": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut); default = 10. '
                    + 'Sets a specific radius in Angstroms that non-bonded interaction ' \
                    + 'energy and force will be considered and calculated using defined potential function. ' \
                    + 'The distance in Angstoms to truncate the LJ, Mie, or other VDW type potential at. '\
                    + 'Note: Rswitch is only used when the "Potential" = SWITCH. ',
        "RcutLow": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut); default = 1. '\
                    + 'Sets a specific minimum possible distance in Angstroms that reject ' \
                    + 'any move that places any atom closer than specified distance. The minimum possible ' \
                    + 'distance between any atoms. '  \
                    + 'Sets a specific radius in Angstroms that non-bonded interaction '\
                    + 'Note: Rswitch is only used when the "Potential" = SWITCH. ',
        "LRC": 'Simulation info (all ensembles): boolean; default = True. ' \
                    + 'If True, the simulation considers the long range tail corrections for the non-bonded VDW or '\
                    + 'dispersion interactions. ' \
                    + 'Note: In case of using SHIFT or SWITCH potential functions, LRC will be ignored.',
        "Exclude": 'Simulation info (all ensembles): str ' \
                    + '(The string inputs are "1-2", "1-3", or "1-4"); default = "1-3". ' \
                    + 'Note: In CHARMM force field, the 1-4 interaction needs to be considered. ' \
                    + 'Choosing "Excude 1-3", will modify 1-4 interaction based on 1-4 parameters ' \
                    + 'in parameter file. If a kind force field is used, where ' \
                    + '1-4 interaction needs to be ignored, such as TraPPE, either "Excude 1-4" needs to be ' \
                    + 'chosen or 1-4 parameter needs to be assigned to zero in the parameter file. \n'
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "1-2": All interaction pairs of bonded atoms, ' \
                    + 'except the ones that separated with one bond, ' \
                    + 'will be considered and modified using 1-4 parameters defined in parameter file. \n' \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "1-3": All interaction pairs of bonded atoms, ' \
                    + 'except the ones that separated with one or two ' \
                    + 'bonds, will be considered and modified using 1-4 parameters defined in parameter file. \n' \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "1-4": All interaction pairs of bonded atoms, ' \
                    + 'except the ones that separated with one, ' \
                    + 'two or three bonds, will be considered using non-bonded parameters defined in parameter file',
        "Potential": 'Simulation info (all ensembles): str ' \
                    + '(The string inputs are "VDW", "EXP6", "SHIFT" and "SWITCH"); default = "VDW". ' \
                    + 'Defines the potential function type to calculate non-bonded dispersion interaction ' \
                    + 'energy and force between atoms. \n' \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "VDW":    Nonbonded dispersion interaction energy and force ' \
                    + 'calculated based on n-6(Lennard - Jones) equation. This function will be discussed ' \
                    + 'further in the Intermolecular energy and ' \
                    + 'Virial calculation section. \n'  \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "EXP6":   Nonbonded dispersion interaction energy and force ' \
                    + 'calculated based on exp-6(Buckingham potential) equation. \n' \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "SHIFT":  This option forces the potential energy to be ' \
                    + 'zero at Rcut distance.  \n' \
                    + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- "SWITCH": This option smoothly forces the potential ' \
                    + 'energy to be zero at Rcut distance and starts modifying the potential at Rswitch ' \
                    + 'distance. Depending on force field type, specific potential function will be applied. ' ,
        "Rswitch": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut); default = 12. '\
                    + 'Note: Rswitch is only used when the SWITCH function is used (i.e., "Potential" = SWITCH). '\
                    + 'The Rswitch distance is in Angstrom. If the “SWITCH” function is chosen, ' \
                    + 'Rswitch needs to be defined; otherwise, the program will be terminated. When using ' \
                    + 'choosing "SWITCH" as potential function, the Rswitch distance defines where the' \
                    + 'non-bonded interaction energy modification is started, which is eventually truncated ' \
                    + 'smoothly at Rcut distance.',
        "VDWGeometricSigma": 'Simulation info (all ensembles): boolean; default = False. ' \
                    + 'Use geometric mean, as required by OPLS force field, ' \
                    + 'to combining Lennard-Jones sigma parameters for different atom types. If set to True, GOMC ' \
                    + 'uses geometric mean to combine Lennard-Jones or VDW sigmas. Note: The default setting of  ' \
                    + 'VDWGeometricSigma is false to use arithmetic mean when combining Lennard-Jones or VDW ' \
                    + 'sigma parameters for different atom types.',
        "ElectroStatic": 'Simulation info (all ensembles): boolean; default = True. ' \
                    + 'Considers the coulomb interactions or not. '\
                    + 'If True, coulomb interactions are considered and false if not.  ' \
                    + 'Note: To simulate the polar molecule in MARTINI force field, ElectroStatic needs to be  ' \
                    + 'turn on. MARTINI force field uses short-range coulomb interaction with constant  ' \
                    + 'Dielectric 15.0.',
        "Ewald":  'Simulation info (all ensembles): boolean; default = True. ' \
                    + 'Considers the standard Ewald summation method for electrostatic calculations. ' \
                    + 'If True, Ewald summation calculation needs to be considered and false if not. '\
                    + 'Note: By default, GOMC will set ElectroStatic to True if Ewald summation  ' \
                    + 'method was used to calculate coulomb interaction.',
        "CachedFourier": 'Simulation info (all ensembles): boolean; default = False. ' \
                    + 'Considers storing the reciprocal terms for Ewald summation ' \
                    + 'calculation in order to improve the code performance. This option would increase the code ' \
                    + 'performance with the cost of memory usage. If True, to store reciprocal terms of Ewald ' \
                    + 'summation calculation and False if not. ' \
                    + 'Warning: Monte Carlo moves, such as MEMC-1, MEMC-2, MEMC-3, ' \
                    + 'IntraMEMC-1, IntraMEMC-2, and IntraMEMC-3 are not support with CachedFourier.',
        "Tolerance": 'Simulation info (all ensembles): float (0.0 < float < 1.0); default = 0.00001. ' \
                    + 'Sets the accuracy in Ewald summation calculation. Ewald separation parameter and number ' \
                    + 'of reciprocal vectors for the Ewald summation are determined based on the accuracy parameter',
        "Dielectric": 'Simulation info (all ensembles): int or float (>= 0.0); default = 15. ' \
                    + 'Sets dielectric value used in coulomb interaction when the Martini ' \
                    + 'force field is used. Note: In MARTINI force field, Dielectric needs to be set to 15.0.',
        "PressureCalc": 'Simulation info (all ensembles): list [bool , int (> 0)] or [bool , step_frequency]; ' \
                    + 'default = [True , set via formula or 10,000 max]. ' \
                    + 'Considers to calculate the pressure or not. bool = True, enables the pressure calculation ' \
                    + 'during the simulation, false disables the calculation. The int/step frequency sets the ' \
                    + 'frequency of calculating the pressure.',
        "EqSteps": 'Simulation info (all ensembles): int (> 0); default = set via formula or 1M max. ' \
                    + 'Sets the number of steps necessary to equilibrate the system ' \
                    + 'Averaging will begin at this step. ' \
                    + 'Note: In GCMC simulations, the Histogram files will be outputed at EqSteps.',
        "AdjSteps": 'Simulation info (all ensembles): int (> 0); default = set via formula or 1,000 max. ' \
                    + 'Number of steps per move adjustment. ' \
                    + 'Sets the number of steps per adjustment of the parameter associated with each move ' \
                    + '(e.g. maximum translate distance, maximum rotation, maximum volume exchange, etc.)',
        "useConstantArea": 'Simulation info (only GEMC_NPT and NPT): boolean: default = False. ' \
                    + 'Considers to change the volume of the simulation box by fixing the cross-sectional ' \
                    + 'area (x-y plane). If true, the volume will change only in z axis, If false, the volume ' \
                    + 'will change with constant axis ratio. ',
        "FixVolBox0": 'Simulation info (only GEMC_NPT): boolean; default = False . ' \
                    + 'Changing the volume of fluid phase (Box 1) to maintain the constant imposed pressure and ' \
                    + 'temperature, while keeping the volume of adsorbed phase (Box 0) fixed. Note: By default, ' \
                    + 'GOMC will set useConstantArea to False if no value was set. It means, the volume of the ' \
                    + 'box will change in a way to maintain the constant axis ratio.',
        # GCMC only properties
        "ChemPot": 'Simulation info (only GCMC): dict {str (4 dig limit) , int or float}; '  \
                    + 'default = None (i.e., user must set as there is no working default).' \
                    + 'The chemical potentials in GOMC units of energy, K. ' \
                    + 'There is a 4 character limit for the string/residue name since the PDB/PSF ' \
                    +  'Note: These strings must match the residue in the psf and psb files or it will fail. ' \
                    + 'files have a 4 character limitation and require and exact match in the conf file. ' \
                    + 'The name of the residues and their corresponding chemical potential must specified ' \
                    + 'for ever residue in the system (i.e., {"residue_name" : chemical_potential}). ' \
                    + 'Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY ' \
                    + 'OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER. ' \
                    + 'Example 1 (system with only water):  {"H2O" : -4000} . ' \
                    + 'Example 2 (system with water and ethanol):  {"H2O" : -4000, "ETH" : -8000} ',
        "Fugacity": 'Simulation info (only GCMC): dict {str , int or float (>= 0)}; '  \
                    + 'default = None (i.e., user must set as there is no working default). ' \
                    + 'The fugacity in GOMC units of pressure, bar. ' \
                    + 'There is a 4 character limit for the string/residue name since the PDB/PSF' \
                    + 'Note: These strings must match the residue in the psf and psb files or it will fail. ' \
                    + 'files have a 4 character limitation and require and exact match in the conf file. ' \
                   + 'The name of the residues and their corresponding fugacity must specified ' \
                   + 'for ever residue in the system (i.e., {"residue_name" : fugacity}). ' \
                   + 'Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY ' \
                   + 'OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER. ' \
                    + 'Example 1 (system with only water):  {"H2O" : 1} . ' \
                    + 'Example 2 (system with water and ethanol):  {"H2O" : 0.5, "ETH" : 10} ' ,

        # CBMC inputs
        "CBMC_First": 'CBMC inputs (all ensembles): int (>= 0); default = 12, ' \
                    + 'The Number of CD-CBMC trials to choose the first atom position' \
                    + '(Lennard-Jones trials for first seed growth).',
        "CBMC_Nth": 'CBMC inputs (all ensembles): int (>= 0); default = 10,  ' \
                    + 'The Number of CD-CBMC trials to choose the later atom positions ' \
                    + '(Lennard-Jones trials for first seed growth).',
        "CBMC_Ang": 'CBMC inputs (all ensembles): int (>= 0); default = 50, ' \
                    + 'The Number of CD-CBMC bending angle trials to perform for geometry ' \
                    + '(per the coupled-decoupled CBMC scheme).',
        "CBMC_Dih": 'CBMC inputs (all ensembles): int (>= 0); default = 50, ' \
                    + 'The Number of CD-CBMC dihedral angle trials to perform for geometry ' \
                    + '(per the coupled-decoupled CBMC scheme).',


        # Control file (.conf file ) output controls/parameters
        "OutputName": 'Output Frequency (all ensembles): str (NO SPACES); default = "Output_data". ' \
                      + 'UNIQUE STRING NAME WITH NO SPACES for simulation used to name the block average, ' \
                      + 'PDB, and PSF output files.',
        "CoordinatesFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                      + '[Generate_data_bool , steps_per_data_output_int]; ' \
                      + 'default = [True , set via formula or 1M max]. ' \
                      + 'Controls output of PDB file (coordinates). ' \
                      + 'If bool is True, this enables outputing the coordinate files at the ' \
                      + 'integer frequency (set steps_per_data_ouput_int), ' \
                      + 'while "False" disables outputing the coordinates.',
        "RestartFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                       + '[Generate_data_bool , steps_per_data_output_int]; ' \
                       + 'default = [True , set via formula or 1M max], ' \
                       + 'This creates the PDB and PSF (coordinate and topology) files for restarting the system ' \
                       + 'at the set steps_per_data_ouput_int (frequency) '\
                       + 'If bool is True, this enables outputing the PDB/PSF restart files at the ' \
                       + 'integer frequency (set steps_per_data_ouput_int), ' \
                       + 'while “false” disables outputing the PDB/PSF restart files. ',
        "CheckpointFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                      + '[Generate_data_bool , steps_per_data_output_int]; ' \
                      + 'default = [True , set via formula or 1M max], ' \
                      + 'Controls the output of the last state of simulation at a specified step, in a ' \
                      + 'binary file format (checkpoint.dat). Checkpoint file contains the following ' \
                      + 'information in full precision: ' \
                      + '(1) Last simulation step that saved into checkpoint file. '\
                      + '(2) Simulation cell dimensions and angles. ' \
                      + '(3) Maximum amount of displacement (Å), rotation (δ), and volume (Å^3) that is used in ' \
                      + 'Displacement, Rotation, MultiParticle, and Volume move. ' \
                      + '(4) Number of Monte Carlo move trial and acceptance. ' \
                      + '(5) All molecule’s coordinates. ' \
                      + '(6) Random number sequence. ' \
                      + 'If bool is True, this enables outputing the checkpoint file at the ' \
                      + 'integer frequency (set steps_per_data_ouput_int), ' \
                      + 'while "False" disables outputing the checkpoint file.',
          "ConsoleFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                      + '[Generate_data_bool , steps_per_data_output_int]; ' \
                      + 'default = [True , set via formular or 10,000 max]. ' \
                      + 'Controls the output to STDIO (“the console” or log file) of messages such as ' \
                      + 'acceptance statistics, and run timing info. In addition, instantaneously-selected ' \
                      + 'thermodynamic properties will be output to this file.' \
                      + 'If bool is True, this enables outputing the consol data at the ' \
                      + 'integer frequency (set steps_per_data_ouput_int), ' \
                      + 'while "False" disables outputing the consol data file.',
        "BlockAverageFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                      + '[Generate_data_bool , steps_per_data_output_int]; ' \
                      + 'default = [True , set via formula or 10,000 max]. ' \
                      + 'Controls the block averages output of selected thermodynamic properties. ' \
                      + 'Block averages are averages of thermodynamic values of interest for chunks of the ' \
                      + 'simulation (for post-processing of averages or std. dev. in those values).' \
                      + 'If bool is True, this enables outputing the block averaging data/file at the ' \
                      + 'integer frequency (set steps_per_data_ouput_int), ' \
                      + 'while "False" disables outputing the block averaging data/file.',
        "HistogramFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' \
                      + '[Generate_data_bool , steps_per_data_output_int]; ' \
                      + 'default = [True , set via formula or 10,000 max]. ' \
                      + 'Controls the histograms. Histograms are a binned listing of observation frequency ' \
                      + 'for a specific thermodynamic variable. In this code, they also control the output ' \
                      + 'of a file containing energy/molecule samples; ' \
                      + 'it only will be used in "GCMC" ensemble simulations for histogram reweighting purposes.' \
                      + 'If bool is True, this enables outputing the data to the histogram data at the ' \
                      + 'integer frequency (set steps_per_data_ouput_int), ' \
                      + 'while "False" disables outputing the histogram data.',

        # Histogram data
        "DistName": 'Histogram Output (all ensembles): str (NO SPACES); default = "dis". '  \
                    + 'Short phrase which will be combined with RunNumber and RunLetter ' \
                    + 'to use in the name of the binned histogram for molecule distribution.' \
                    + 'Sets short phrase to naming molecule distribution file.',
        "HistName": 'Histogram Output (all ensembles): str (NO SPACES); default = "his". '  \
                    + 'Short phrase, which will be combined with RunNumber and RunLetter, ' \
                    + 'to use in the name of the energy/molecule count sample file.' \
                    + 'Sets short phrase to naming energy sample file.',
        "RunNumber": 'Histogram Output (all ensembles): int  ( > 0 ); default = 1. '  \
                    + 'Run number to be used in the above file names.  ' \
                    + 'Sets a number, which is a part of DistName and HistName file name.',
        "RunLetter": 'Histogram Output (all ensembles): str (1 alphabetic character only); default = "a". '  \
                    + 'Run letter to be used in above file names.' \
                    + 'Sets a letter, which is a part of DistName and HistName file name.',
        "SampleFreq": 'Histogram Output (all ensembles): int ( > 0 ); default = 500. '  \
                    + 'The number of steps per histogram sample.' \
                    + 'Controls histogram sampling frequency.',

        # Data output for the consol and bulk properties calculations
        "OutEnergy": 'Output Data (all ensembles): [bool, bool]; default = [True, True].   '\
                    + 'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the energy data into the block averages and consol output/log files.',
        "OutPressure": 'Output Data (all ensembles): [bool, bool]; default = [True, True].   '\
                    + 'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the pressure data into the block averages and consol output/log files.',
        "OutMolNumber": 'Output Data (all ensembles): [bool, bool]; default = [True, True].   '\
                    + 'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the number of molecules data into the block averages and consol output/log files.',
        "OutDensity": 'Output Data (all ensembles): [bool, bool]; default = [True, True].   '\
                    + 'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the density data into the block averages and consol output/log files.',
        "OutVolume": 'Output Data (all ensembles): [bool, bool]; default = [True, True].   '\
                    + 'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the volume data into the block averages and consol output/log files.',
        "OutSurfaceTension": 'Output Data (all ensembles): [bool, bool]; default = [False, False]. ' \
                    +  'The list provides the booleans to [block_averages_bool, consol_output_bool]. '\
                    + 'This ouputs the surface tension data into the block averages and consol output/log files.',

        # free energy calculation in NVT and NPT ensembles.
        "FreeEnergyCalc": 'Free Energy Calcs (NVT and NPT only): list [bool , int (> 0)] or '  \
                    + '[Generate_data_bool , steps_per_data_output_int]; default = None. ' \
                    + 'bool = True enabling free energy calculation during the simulation, false disables '  \
                    + 'the calculation. The int/step frequency sets the frequency of calculating the free energy.',
        "MoleculeType": 'Free Energy Calcs (NVT and NPT only): list [str , int (> 0)] or '  \
                    + '["residue_name" , residue_ID]; ' \
                    + 'user must set as there is no working default (default = None). ' \
                    + 'Note: ONLY 4 characters can be used for the string (i.e., "residue_name"). ' \
                    + 'Sets the solute molecule kind (residue name) and molecule number (residue ID), '  \
                    + 'which absolute solvation free will be calculated for.',
        "InitialState": 'Free Energy Calcs (NVT and NPT only): int (>= 0); user must set as there is no ' \
                    + 'usable default (default = None). ' \
                    + 'The index of LambdaCoulomb and LambdaVDW vectors. Sets the index of the' \
                    + 'LambdaCoulomb and LambdaVDW vectors, to determine the simulation lambda value for' \
                    + 'VDW and Coulomb interactions. ' \
                    + 'WARNRING : This must an integer within the vector count of the LambdaVDW and LambdaCoulomb, ' \
                    + 'in which the counting starts at 0.  ',
        "LambdaVDW": 'Free Energy Calcs (NVT and NPT only): list of floats (0 <= floats <= 1); ' \
                    + 'user must set as there is no usable default (default = None). ' \
                    + 'Lambda values for VDW interaction in ascending order. Sets the intermediate ' \
                    + 'lambda states to which solute-solvent VDW interactions are scaled. ' \
                    + 'WARNRING : All lambda values must be stated in the ascending order, otherwise the program ' \
                    + 'will terminate. ' \
                    + 'WARNRING : This list must be the same length as the LambdaCoulomb list length.' \
                    + 'Example 1: [0.1, 1.0,] . ' \
                    + 'Example 2: [0.1, 0.2, 0.4, 0.9] . ',
        "LambdaCoulomb": 'Free Energy Calcs (NVT and NPT only):  list of floats (0 <= floats <= 1); ' \
                    + 'user must set as there is no usable default (default = None). ' \
                    + 'Lambda values for Coulombic interaction in ascending order. Sets the intermediate '\
                    + 'lambda states to which solute-solvent Coulombic interactions are scaled.' \
                    + 'WARNRING : All lambda values must be stated in the ascending order, otherwise the program ' \
                    + 'will terminate.  ' \
                    + 'WARNRING : This list must be the same length as the LambdaVDW list.' \
                    + 'NOTE: By default (i.e., LambdaCoulomb = None or default = None),' \
                    + 'the lambda values for Coulombic interaction will be set to zero if ' \
                    + "ElectroStatic or Ewald is deactivated. By GOMC's default, or the lambda values for Coulombic " \
                    + 'interaction will be set to Lambda values for VDW interaction if ElectroStatic or ' \
                    + 'Ewald is activated.' \
                    + 'Example 1: [0.1, 1.0,] . ' \
                    + 'Example 2: [0.1, 0.2, 0.4, 0.9] . ',
        "ScaleCoulomb": 'Free Energy Calcs (NVT and NPT only): bool; default = False, '\
                    + 'True if coulombic interaction needs to be scaled non-linearly, ' \
                    + 'False if coulombic interaction needs to be scaled linearly. Determines to scale the ' \
                    + 'Coulombic interaction non-linearly (soft-core scheme) or not.',
        "ScalePower": 'Free Energy Calcs (NVT and NPT only): int (>= 0); default = 2, '\
                    + 'The p value in the soft-core scaling scheme.  Sets the p value in ' \
                    + 'soft-core scaling scheme, where the distance between solute and solvent is scaled' \
                    + 'non-linearly.',
        "ScaleAlpha": 'Free Energy Calcs (NVT and NPT only): int or float (>= 0), default = 0.5, '\
                    + 'alpha vaule in the soft-core scaling scheme. Sets the α value' \
                    + 'in soft-core scaling scheme, where the distance between solute and solvent is scaled' \
                    + 'non-linearly.',
        "MinSigma": 'Free Energy Calcs (NVT and NPT only): int or float (>= 0), default = 3, '\
                    + 'Minimum sigma value in the soft-core scaling scheme.' \
                    + 'Sets the minimum σ value in soft-core scaling scheme, where the distance between' \
                    + 'solute and solvent is scaled non-linearly.',


        # moves without MEMC
        "DisFreq": 'Std. MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.15, NPT = 0.15, GEMC_NVT = 0.20, GEMC_NPT = 0.19, GCMC = 0.15). ' \
                    + 'Fractional percentage at which the displacement move will occur ' \
                    + '(i.e., fraction of displacement moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "RotFreq": 'Std. MC moves (all ensembles)                     : ' \
                     + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.15, NPT = 0.15, GEMC_NVT = 0.20, GEMC_NPT = 0.20, GCMC = 0.15). ' \
                    + 'Fractional percentage at which the rotation move will occur ' \
                   + '(i.e., fraction of rotation moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "IntraSwapFreq": 'Std. MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.30, NPT = 0.29, GEMC_NVT = 0.10, GEMC_NPT = 0.10, GCMC = 0.10). ' \
                    + 'Fractional percentage at which the molecule will be removed from a ' \
                    + 'box and inserted into the same box using coupled-decoupled configurational-bias' \
                    + 'algorithm. (i.e., fraction of intra-molecule swap moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "SwapFreq": 'Std. MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.20, GEMC_NPT = 0.20, GCMC = 0.35). ' \
                    + 'For Gibbs and Grand Canonical (GC) ensemble runs only: Fractional ' \
                    + 'percentage at which molecule swap move will occur using coupled-decoupled ' \
                    + 'configurational-bias. (i.e., fraction of molecule swaps moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "RegrowthFreq": 'Std. MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.30, NPT = 0.30, GEMC_NVT = 0.20, GEMC_NPT = 0.20, GCMC = 0.15). ' \
                    + 'Fractional percentage at which part of the molecule will be ' \
                    + 'deleted and then regrown using coupled- decoupled configurational-bias algorithm ' \
                    + '(i.e., fraction of molecular growth moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "CrankShaftFreq": 'Std. MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.10, NPT = 0.10, GEMC_NVT = 0.10, GEMC_NPT = 0.10, GCMC = 0.10). ' \
                    + 'Fractional percentage at which crankshaft move will occur. ' \
                    + 'In this move, two atoms that are forming angle or dihedral are selected randomly and ' \
                    + 'form a shaft. Then any atoms or group that are within these two selected atoms, will ' \
                    + 'rotate around the shaft to sample intra-molecular degree of freedom ' \
                    + '(i.e., fraction of crankshaft moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "VolFreq": 'Std. MC moves (only  GEMC_NPT  and  NPT )         : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.00, NPT = 0.01, GEMC_NVT = 0.00, GEMC_NPT = 0.01, GCMC = 0.00). ' \
                    + 'Fractional percentage at which molecule will be removed from one box and inserted into ' \
                    + 'the other box using configurational bias algorithm ' \
                    + '(i.e., fraction of Volume swaps moves) Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',
        "MultiParticleFreq": 'Std. MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
                    + 'Fractional percentage at which multi-particle move will ' \
                    + 'occur. In this move, all molecules in the selected simulation box will be rigidly ' \
                    + 'rotated or displaced simultaneously, along the calculated torque or force ' \
                    + 'respectively (i.e., fraction of multi-particle moves). Note: all of the move types'  \
                    + 'are not available in for every ensemble. Note: all of the move fractions'  \
                    + 'must sum to 1, or the control file writer will fail.  ',

        # MEMC moves
        "IntraMEMC-1Freq": 'MEMC MC moves (all ensembles)                     : ' \
                    + 'int or float (0 <= value <= 1); default are specific for each ' \
                    + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
                    + 'Fractional percentage at which specified number of small molecule kind will be ' \
                    + 'exchanged with a specified large molecule kind in defined sub-volume within ' \
                    + 'same simulation box.  This move need additional information such as ' \
                    + 'ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.' \
                    + 'Note: all of the move types are not available in for every ensemble.' \
                    + 'Note: all of the move fractions must sum to 1, or the control file writer will fail.  ',
        "MEMC-1Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' \
            + 'int or float (0 <= value <= 1); default are specific for each ' \
            + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
            + 'Fractional percentage at which specified number of small molecule kind will be exchanged with ' \
            + 'a specified large molecule kind in defined sub-volume in dense simulation box. This move need ' \
            + 'additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, ' \
            + 'and ExchangeLargeKind.' \
            + 'Note: all of the move types are not available in for every ensemble.' \
            + 'Note: all of the move fractions must sum to 1, or the control file writer will fail.  ',
        "IntraMEMC-2Freq": 'MEMC MC moves (all ensembles)                     : ' \
            + 'int or float (0 <= value <= 1); default are specific for each ' \
            + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
            + 'Fractional percentage at which specified number of small molecule kind will be exchanged with ' \
            + 'a specified large molecule kind in defined sub-volume within same simulation box. ' \
            + 'Backbone of small and large molecule kind will be used to insert the large molecule more efficiently. '\
            + 'This move need additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, ' \
            + 'ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            + 'Note: all of the move types are not available in for every ensemble.' \
            + 'Note: all of the move fractions must sum to 1, or the control file writer will fail.  ',
        "MEMC-2Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' \
            + 'int or float (0 <= value <= 1); default are specific for each ' \
            + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
            + 'Fractional percentage at which specified number of small molecule kind will be exchanged with ' \
            + 'a specified large molecule kind in defined sub-volume in dense simulation box. Backbone of small ' \
            + 'and large molecule kind will be used to insert the large molecule more efficiently. ' \
            + 'This move need additional information such as ExchangeVolumeDim, ExchangeRatio, ' \
            + 'ExchangeSmallKind, ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            + 'Note: all of the move types are not available in for every ensemble.' \
            + 'Note: all of the move fractions must sum to 1, or the control file writer will fail.  ',
        "IntraMEMC-3Freq": 'MEMC MC moves (all ensembles)                     : ' \
            + 'int or float (0 <= value <= 1); default are specific for each ' \
            + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
            + 'Fractional percentage at which specified number of small molecule kind will be exchanged with ' \
            + 'a specified large molecule kind in defined sub-volume within same simulation box. Specified atom ' \
            + 'of the large molecule kind will be used to insert the large molecule using coupled-decoupled ' \
            + 'configurational-bias. This move need additional information such as ExchangeVolumeDim, ' \
            + 'ExchangeRatio, ExchangeSmallKind, ExchangeLargeKind, and LargeKindBackBone. ' \
            + 'Note: all of the move types are not available in for every ensemble.' \
            + 'Note: all of the move fractions must sum to 1, or the control file writer will fail. ',
        "MEMC-3Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' \
            + 'int or float (0 <= value <= 1); default are specific for each ' \
            + 'ensemble (NVT = 0.00, NPT = 0.00, GEMC_NVT = 0.00, GEMC_NPT = 0.00, GCMC = 0.00). ' \
            + 'Fractional percentage at which specified number of small molecule kind will be exchanged with ' \
            + 'a specified large molecule kind in defined sub-volume in dense simulation box. Specified atom ' \
            + 'of the large molecule kind will be used to insert the large molecule using coupled-decoupled ' \
            + 'configurational-bias. This move need additional information such as ExchangeVolumeDim, ' \
            + 'ExchangeRatio, ExchangeSmallKind, ExchangeLargeKind, and LargeKindBackBone. ' \
            + 'Note: all of the move types are not available in for every ensemble.' \
            + 'Note: all of the move fractions must sum to 1, or the control file writer will fail.  ',

        # MEMC move parameters
        "ExchangeVolumeDim": 'MEMC parameters (all ensembles)                   : ' \
                      + 'list of 3 floats or integers ' \
                      +  '[int or float (> 0), int or float (> 0), int or float (> 0)]'\
                      +  ' or [X-dimension, Y-dimension, Z-dimension)]; '\
                      + 'default is [1.0, 1.0, 1.0]. ' \
                      + 'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, the exchange ' \
                      + 'subvolume must be defined. The exchange sub-volume is defined as an orthogonal box ' \
                      + 'with x, y, and z-dimensions, where small molecule/molecules kind will be selected ' \
                      + 'from to be exchanged with a large molecule kind. ' \
                      + 'Note: Currently, the X and Y dimension cannot be set independently (X = Y = max(X, Y)). ' \
                      + 'Note: A heuristic for setting good values of the x, y, and z-dimensions is to use' \
                      + 'the geometric size of the large molecule plus 1-2 Å in each dimension. ' \
                      + 'Note: In case of exchanging 1 small molecule kind with 1 large molecule kind in ' \
                      + 'IntraMEMC-2, IntraMEMC-3, MEMC-2, MEMC-3 Monte Carlo moves, the sub-volume ' \
                      + 'dimension has no effect on acceptance rate. ',
        "MEMC_DataInput": 'MEMC parameters (availablity based on selelection): nested lists; default = None.  ' \
                         + 'Enter data as a list with some sub-lists as follows: ' \
                      + '[[ExchangeRatio_int (> 0), ExchangeLargeKind_str, ' \
                      + '[LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ], ' \
                      + 'ExchangeSmallKind_str, ' \
                      + '[SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]], ..., ' \
                      + '[ExchangeRatio_int (> 0), ExchangeLargeKind_str, ' \
                      + '[LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ], ' \
                      + 'ExchangeSmallKind_str, '
                      + '[SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]. ' \
                      + 'NOTE: CURRENTLY ALL THESE INPUTS NEED TO BE SPECIFIED, REGARDLESS OF THE MEMC TYPE ' \
                      + 'SELECTION. IF THE SmallKindBackBone or LargeKindBackBone IS NOT REQUIRED FOR THE MEMC TYPE, '\
                      + 'None CAN BE USED IN PLACE OF A STRING. ' \
                      + 'Note: These strings must match the residue in the psf and psb files or it will fail. ' \
                      + 'It is recommended that the user print the Charmm object psf and pdb files '
                      + 'and review the residue names that match the atom name before using the in '
                      + 'the  MEMC_DataInput variable input'
                      + 'Note: see the below data explanations for the ExchangeRatio, ExchangeSmallKind, '\
                      + 'ExchangeLargeKind, LargeKindBackBone, SmallKindBackBone. ' \
                      + "Example 1 (MEMC-1) : [ [1, 'WAT', [None, None], 'wat', [None, None]] , "
                      + "[1, 'WAT', [None, None], 'wat', [None, None]] . "
                      + "Example 2 (MEMC-2): [ [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1' ]] , "
                      + " [1, 'WAT', ['H1', 'H2'], 'wat', ['H1', 'H2' ]] . "
                      + "Example 3 (MEMC-3) : [ [2, 'WAT', 'O1', 'H1'], 'wat', [None, None]] , "
                      + "[2, 'WAT', ['H1', 'H2'], 'wat', [None, None]] .\n"
                +'\t\t\t\t\t\t\t\t\t\t\t\t\t -- ExchangeRatio     = MEMC parameters (all ensembles): ' \
                                     + 'int (> 0); default = None. The Ratio of exchanging ' \
                                 + 'small molecule/molecules with 1 large molecule. ' \
                                 + 'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, ' \
                                 + 'the exchange ratio must be defined. ' \
                                 + 'The exchange ratio defines how many small molecule will be ' \
                                 + 'exchanged with 1 large molecule. For each large-small molecule pairs, ' \
                                 + 'one exchange ratio must be defined. \n' \
                + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- ExchangeSmallKind = MEMC parameters (all ensembles): ' \
                                  + 'str; default = None. The small molecule ' \
                                  + 'kind (resname) to be exchanged. ' \
                                  + 'Note: ONLY 4 characters can be used for the strings. ' \
                                  + 'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, ' \
                                 + 'the small molecule kind to be exchanged with a large molecule ' \
                                 + 'kind must be defined. Multiple small molecule kind can be specified.  \n' \
                + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- ExchangeLargeKind = MEMC parameters (all ensembles): ' \
                                 + 'str; default = None. The large molecule ' \
                                + 'kind (resname) to be exchanged. ' \
                                + 'Note: ONLY 4 characters can be used for the strings. ' \
                                + 'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, ' \
                                + 'the large molecule kind to be exchanged with small molecule ' \
                                 + 'kind must be defined. Multiple large molecule kind can be specified. \n' \
                + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- LargeKindBackBone = MEMC parameters (all ensembles): ' \
                                + '2 strings in a list [str, str] or 2 Nones in a list [None, None]; ' \
                                + 'default = None. ' \
                                + 'Note: ONLY 4 characters can be used for the strings. ' \
                                + 'The [None, None] values can only be used if that MEMC type does not require them. ' \
                                + 'The strings for the the atom name 1 and atom name 2 that belong to the large ' \
                                + 'molecule’s backbone (i.e., [str_for_atom_name_1, str_for_atom_name_2]) '  \
                                + 'To use MEMC-2, MEMC-3, IntraMEMC-2, and IntraMEMC-3 Monte Carlo moves, the ' \
                                + 'large molecule backbone must be defined. The backbone of the molecule is defined ' \
                                + 'as a vector that connects two atoms belong to the large molecule. The large ' \
                                + 'molecule backbone will be used to align the sub-volume in MEMC-2 and IntraMEMC-2 ' \
                                + 'moves, while in MEMC-3 and IntraMEMC-3 moves, it uses the atom name to start ' \
                                + 'growing the large molecule using coupled-decoupled configurational-bias. For ' \
                                + 'each large-small molecule pairs, two atom names must be defined. ' \
                                + 'Note: all atom names in the molecule must be unique. ' \
                                + 'Note: In MEMC-3 and IntraMEMC-3 Monte Carlo moves, both atom names must be same, ' \
                                + 'otherwise program will be terminated. ' \
                                + 'Note: If the large molecule has only one atom (mono atomic molecules), ' \
                                + 'same atom name must be used for str_for_atom_name_1 and str_for_atom_name_2 ' \
                                + 'of the LargeKindBackBone.  \n' \
                + '\t\t\t\t\t\t\t\t\t\t\t\t\t -- SmallKindBackBone = MEMC parameters (all ensembles): ' \
                                + '2 strings in a list [str, str] or 2 Nones in a list [None, None]; default = None. ' \
                                + 'Note: ONLY 4 characters can be used for the strings. ' \
                                + 'The [None, None] values can only be used if that MEMC type does not require them.' \
                                + 'The strings for the the atom name 1 and atom name 2 that belong to the small ' \
                                + 'molecule’s backbone (i.e., [str_for_atom_name_1, str_for_atom_name_2]) ' \
                                + 'To use MEMC-2, and IntraMEMC-2 Monte Carlo moves, the small molecule backbone ' \
                                + 'must be defined. The backbone of the molecule is defined as a vector that ' \
                                + 'connects two atoms belong to the small molecule and will be used to align the ' \
                                + 'sub-volume. For each large-small molecule pairs, two atom names must be defined. ' \
                                + 'Note: all atom names in the molecule must be unique. '  \
                                + 'Note: If the small molecule has only one atom (mono atomic molecules), same atom ' \
                                + 'name must be used str_for_atom_name_1 and str_for_atom_name_2 ' \
                                + 'of the SmallKindBackBone. ',

        # ******************************************************************************************************
        # Definitions in this function are copied to a large extent from the GOMC manual release version 2.50 (end)
        # insert citation here:
        # ******************************************************************************************************
    }
    if description:
        return valid_input_variables
    else:
        return list(valid_input_variables.keys())


def _get_default_variables_dict():
    """
    Provides a dict of the default variables inputs (user optional)

    Parameters
    ----------
    None

    Outputs
    ---------
    default_input_variables_dict =  Provides a dict of the default variables inputs (user optional)
    """

    default_input_variables_dict = {
        "Restart": False ,
        "RestartCheckpoint" : False  ,
        "PRNG" : "RANDOM",
        "ParaTypeCHARMM": True ,
        "ParaTypeMie": False,
        "ParaTypeMARTINI": False ,
        "RcutCoulomb_box_0" : None,
        "RcutCoulomb_box_1" : None,
        "Pressure": 1.01325,
        "Rcut": 10,
        "RcutLow":  1 ,
        "LRC": True,
        "Exclude": '1-3' ,
        "coul_1_4_scaling" : None,
        "Potential": 'VDW',
        "Rswitch": 9,
        "ElectroStatic": True,
        "Ewald":  True,
        "CachedFourier": False,
        "Tolerance": 0.00001,
        "Dielectric" : 15,
        "PressureCalc": [True , 10000],
        "EqSteps": 1000000,
        "AdjSteps": 1000,
        "useConstantArea": False,
        "FixVolBox0": False,
        # GCMC only properties
        "ChemPot": None,
        "Fugacity": None,

        # CBMC inputs
        "CBMC_First": 12,
        "CBMC_Nth": 10,
        "CBMC_Ang": 50,
        "CBMC_Dih": 50,


        # Control file (.conf file ) output controls/parameters
        "OutputName": "Output_data",
        "CoordinatesFreq": [True , 1000000],
        "RestartFreq": [True , 1000000],
        "CheckpointFreq": [True , 1000000],
        "ConsoleFreq": [True , 10000],
        "BlockAverageFreq": [True , 10000],
        "HistogramFreq": [True , 10000],

        # Histogram data
        "DistName": "dis",
        "HistName": "his",
        "RunNumber": 1,
        "RunLetter": "a",
        "SampleFreq": 500,

        # Data output for the consol and bulk properties calculations
        "OutEnergy": [True, True],
        "OutPressure": [True, True],
        "OutMolNumber": [True, True],
        "OutDensity": [True, True],
        "OutVolume": [True, True],
        "OutSurfaceTension": [False, False],


        # free energy calculation in NVT and NPT ensembles.
        "FreeEnergyCalc": None,
        "MoleculeType": None,
        "InitialState": None,
        "LambdaVDW": None,
        "LambdaCoulomb": None,
        "ScaleCoulomb": False,
        "ScalePower": 2,
        "ScaleAlpha": 0.5,
        "MinSigma": 3,

        # MEMC move info
        "ExchangeVolumeDim": [1.0, 1.0, 1.0],
        "MEMC_DataInput": None,


        # moves without MEMC
        "DisFreq":           {'NVT' : 0.15, 'NPT' : 0.15, 'GEMC_NVT' : 0.20, 'GEMC_NPT' : 0.19, 'GCMC' : 0.15},
        "RotFreq":           {'NVT' : 0.15, 'NPT' : 0.15, 'GEMC_NVT' : 0.20, 'GEMC_NPT' : 0.20, 'GCMC' : 0.15},
        "IntraSwapFreq":     {'NVT' : 0.30, 'NPT' : 0.29, 'GEMC_NVT' : 0.10, 'GEMC_NPT' : 0.10, 'GCMC' : 0.10},
        "SwapFreq":          {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.20, 'GEMC_NPT' : 0.20, 'GCMC' : 0.35},
        "RegrowthFreq":      {'NVT' : 0.30, 'NPT' : 0.30, 'GEMC_NVT' : 0.20, 'GEMC_NPT' : 0.20, 'GCMC' : 0.15},
        "CrankShaftFreq":    {'NVT' : 0.10, 'NPT' : 0.10, 'GEMC_NVT' : 0.10, 'GEMC_NPT' : 0.10, 'GCMC' : 0.10},
        "VolFreq":           {'NVT' : 0.00, 'NPT' : 0.01, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.01, 'GCMC' : 0.00},
        "MultiParticleFreq": {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        # MEMC moves
        "IntraMEMC-1Freq":   {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        "MEMC-1Freq":        {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        "IntraMEMC-2Freq":   {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        "MEMC-2Freq":        {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        "IntraMEMC-3Freq":   {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
        "MEMC-3Freq":        {'NVT' : 0.00, 'NPT' : 0.00, 'GEMC_NVT' : 0.00, 'GEMC_NPT' : 0.00, 'GCMC' : 0.00},
    }

    return default_input_variables_dict



def check_valid_ensemble_files(ensemble_type, testing_ensemble_files_List):
    """
    Checks if all the required ensemble inputs are provided,
        and provides a list of the bad variables in the printed output.

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
    testing_ensemble_files_List =  list, a list containing the required ensemble
        files variables, which will be tested for to see if they are valid.

    Outputs
    ---------
    bool, returns a bool (True or False) depending on if all variables
        are valid or not
    """

    bad_key_inputs_List = []

    required_ensemble_files_List =  _get_required_ensemble_files(ensemble_type)

    ensemble_has_required_ensemble_files_List = True
    for iter in range(0, len(testing_ensemble_files_List)):
        if testing_ensemble_files_List[iter] not in required_ensemble_files_List:
            bad_key_inputs_List.append(testing_ensemble_files_List[iter])
            ensemble_has_required_ensemble_files_List = False

    for iter in range(0, len(required_ensemble_files_List)):
        if required_ensemble_files_List[iter] not in testing_ensemble_files_List:
            bad_key_inputs_List.append(testing_ensemble_files_List[iter])
            ensemble_has_required_ensemble_files_List = False

    if ensemble_has_required_ensemble_files_List:
        return True

    else:
        warn("ERROR: checking the valid ensemble file (check_valid_ensemble_files function), "
                   "shows that all the required files are not included for the " +str(ensemble_type) + " " 
                   "ensemble. \n " + "The bad required ensemble inputs = " +str(bad_key_inputs_List))
        return False


def print_required_ensemble_files(ensemble_type, description=False):
    """
    Prints the required ensemble arguments with an optional description based on the ensemble type

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
    description =  bool, default = False.
        If True, it prints the descriptions of the required ensemble inputs (i.e. dict),
        If False, it only prints the required ensemble inputs without the descriptions (i.e. list)

    Outputs
    ---------
    Prints the required ensemble arguments with an optional description based on the ensemble type

    """


    required_data_dict =  _get_required_data(description=True)
    ensemble_has_all_valid_required_data = True
    required_data =  _get_required_data()
    required_data_List =  _get_required_ensemble_files(ensemble_type)

    for iter in range(0, len(required_data_List)):
        if required_data_List[iter] not in required_data:
            ensemble_has_all_valid_required_data = False
    if ensemble_has_all_valid_required_data and description == False:
        for iter_2 in range(0, len(required_data_List)):
            required_data_iter = required_data_List[iter_2]
            print("{:10s}:    {}".format(str(iter_2), str(required_data_iter)))

    elif ensemble_has_all_valid_required_data and description == True:
        for iter_2 in range(0, len(required_data_List)):
            required_data_iter = required_data_List[iter_2]
            print("{:10s}:    {:30s}    {}".format(str(iter_2), str(required_data_iter),
                                                      str(required_data_dict[required_data_iter])))
    else:
        print("ERROR: Some files in this ensemble are not in the required file list")




def check_valid_ensemble_input_variables(ensemble_type, testing_input_variables_List):
    """
    Checks if all the input variables (user optional) inputs are valid for the given
        ensemble, and provides a list of the bad variables in the printed output.

    Outputs
    ---------
    bool, returns a bool (True or False) depending on if all variables
        are valid or not

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
    testing_input_variables_List =  list, list containing the optional ensemble
        input variables which will be tested for to see if they are valid.
        """

    bad_key_inputs_List = []

    valid_input_variables_List  = _get_possible_ensemble_input_variables(ensemble_type)
    ensemble_has_valid_input_variables_List = True
    for iter in range(0, len(testing_input_variables_List)):
        if testing_input_variables_List[iter] not in valid_input_variables_List:
            bad_key_inputs_List.append(testing_input_variables_List[iter])
            ensemble_has_valid_input_variables_List = False

    if ensemble_has_valid_input_variables_List:
        return [True, bad_key_inputs_List]


    else:
        return [False, bad_key_inputs_List]


def print_valid_ensemble_input_variables(ensemble_type, description=False):
    """
    Prints the arguments for optional variables brief description based on the ensemble type

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
    description =  bool, default = False.
        If True, it prints the descriptions of the optional variable ensemble inputs (i.e. dict),
        If False, it only prints the  optional variable ensemble inputs without the
        descriptions (i.e. list)

    Outputs
    ---------
    Prints the arguments for optional variables brief description based on the ensemble type
    """

    valid_input_variables_dict = _get_all_possible_input_variables(description=True)
    ensemble_has_all_valid_input_variables = True
    all_valid_input_variables = _get_all_possible_input_variables()

    valid_input_variables_List  = _get_possible_ensemble_input_variables(ensemble_type)

    for iter in range(0, len(valid_input_variables_List)):
        if valid_input_variables_List[iter] not in all_valid_input_variables:

            ensemble_has_all_valid_input_variables = False
    if ensemble_has_all_valid_input_variables and description == False:
        for iter_2 in range(0, len(valid_input_variables_List)):
            ensemble_kwarg_iter = valid_input_variables_List[iter_2]
            print("{:10s}:    {}".format(str(iter_2), str(ensemble_kwarg_iter)))

    elif ensemble_has_all_valid_input_variables and description == True:
        for iter_2 in range(0, len(valid_input_variables_List)):
            ensemble_kwarg_iter = valid_input_variables_List[iter_2]
            print("{:10s}:    {:30s}    {}".format(str(iter_2), str(ensemble_kwarg_iter),
                                                      str(valid_input_variables_dict[ensemble_kwarg_iter])))
    else:
        print("ERROR: Some input_variables in the ensemble are not in the main input_variables list")



def _get_possible_ensemble_input_variables(ensemble_type):
    """
    Provides list of the possible optional input variables based on the ensemble type

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'

    Outputs
    ---------
    valid_input_variables_List = list; a list possible optional input variables
        for the provided ensemble type.

    """


    if ensemble_type in ['NPT', 'NVT']:
        sim_info_variables_List = ["Restart", "RestartCheckpoint", "PRNG",
                                      "ParaTypeCHARMM" , "ParaTypeMie",  "ParaTypeMARTINI",
                                      "RcutCoulomb_box_0",
                                      "Pressure",
                                      "Rcut", "RcutLow", "LRC", "Exclude", "Potential", "Rswitch",
                                      "VDWGeometricSigma", "ElectroStatic", "Ewald", "CachedFourier",  "Tolerance",
                                      "Dielectric", "PressureCalc", "EqSteps", "AdjSteps",
                                      "useConstantArea",]

        CBMC_variables_List = ["CBMC_First", "CBMC_Nth", "CBMC_Ang", "CBMC_Dih",]
        output_freq_variables_List = ["OutputName",  "CoordinatesFreq", "RestartFreq", "CheckpointFreq",
                                      "ConsoleFreq", "BlockAverageFreq", "HistogramFreq",]

        histogram_output_variables_List = ["DistName", "HistName", "RunNumber", "RunLetter", "SampleFreq",]

        output_data_variables_List = ["OutEnergy", "OutPressure", "OutMolNumber", "OutDensity",
                                      "OutVolume", "OutSurfaceTension",]

        free_energy_variables_List = ["FreeEnergyCalc", "MoleculeType", "InitialState", "LambdaVDW",
                                      "LambdaCoulomb", "ScaleCoulomb", "ScalePower", "ScaleAlpha", "MinSigma",]

        Std_MC_moves_variables_List = [ "DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                                      "CrankShaftFreq", "VolFreq", "MultiParticleFreq",]

        MEMC_MC_moves_variables_List = ["IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                                        "IntraMEMC-3Freq", "MEMC-3Freq",
                                        "ExchangeVolumeDim", "MEMC_DataInput"]

        valid_input_variables_List = sim_info_variables_List + CBMC_variables_List \
                                     + output_freq_variables_List + histogram_output_variables_List \
                                     + output_data_variables_List + free_energy_variables_List \
                                     + Std_MC_moves_variables_List + MEMC_MC_moves_variables_List


    elif ensemble_type in ['GEMC_NPT', 'GEMC_NVT']:
        sim_info_variables_List = ["Restart", "RestartCheckpoint", "PRNG",
                                      "ParaTypeCHARMM" , "ParaTypeMie",  "ParaTypeMARTINI",
                                      "RcutCoulomb_box_0", "RcutCoulomb_box_1",
                                      "Pressure",
                                      "Rcut", "RcutLow", "LRC", "Exclude", "Potential", "Rswitch",
                                      "VDWGeometricSigma", "ElectroStatic", "Ewald", "CachedFourier",  "Tolerance",
                                      "Dielectric", "PressureCalc", "EqSteps", "AdjSteps",
                                      "useConstantArea", "FixVolBox0",]

        CBMC_variables_List = [ "CBMC_First", "CBMC_Nth", "CBMC_Ang", "CBMC_Dih",]

        output_freq_variables_List = [  "OutputName",  "CoordinatesFreq", "RestartFreq", "CheckpointFreq",
                                      "ConsoleFreq", "BlockAverageFreq", "HistogramFreq",]

        free_energy_variables_List = []  # always empty for GEMC

        histogram_output_variables_List = ["DistName", "HistName", "RunNumber", "RunLetter", "SampleFreq",]

        output_data_variables_List = [ "OutEnergy", "OutPressure", "OutMolNumber", "OutDensity",
                                      "OutVolume", "OutSurfaceTension",]

        Std_MC_moves_variables_List = [ "DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                                      "CrankShaftFreq", "VolFreq", "MultiParticleFreq",]

        MEMC_MC_moves_variables_List = ["IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                                        "IntraMEMC-3Freq", "MEMC-3Freq",
                                        "ExchangeVolumeDim", "MEMC_DataInput"]

        valid_input_variables_List = sim_info_variables_List + CBMC_variables_List \
                                     + output_freq_variables_List + histogram_output_variables_List \
                                     + output_data_variables_List + free_energy_variables_List \
                                     + Std_MC_moves_variables_List + MEMC_MC_moves_variables_List


    elif ensemble_type == 'GCMC':
        sim_info_variables_List = ["Restart", "RestartCheckpoint", "PRNG",
                                      "ParaTypeCHARMM" , "ParaTypeMie",  "ParaTypeMARTINI",
                                      "RcutCoulomb_box_0", "RcutCoulomb_box_1",
                                      "Rcut", "RcutLow", "LRC", "Exclude", "Potential", "Rswitch",
                                      "VDWGeometricSigma", "ElectroStatic", "Ewald", "CachedFourier",  "Tolerance",
                                      "Dielectric",  "PressureCalc",  "EqSteps", "AdjSteps",  "ChemPot", "Fugacity"]

        CBMC_variables_List = ["CBMC_First", "CBMC_Nth", "CBMC_Ang", "CBMC_Dih",]

        output_freq_variables_List = [  "OutputName",  "CoordinatesFreq", "RestartFreq", "CheckpointFreq",
                                      "ConsoleFreq", "BlockAverageFreq", "HistogramFreq",]

        histogram_output_variables_List = [ "DistName", "HistName", "RunNumber", "RunLetter", "SampleFreq",]

        output_data_variables_List = [ "OutEnergy", "OutPressure", "OutMolNumber", "OutDensity",
                                      "OutVolume", "OutSurfaceTension",]

        free_energy_variables_List = [] # always empty for GCMC

        Std_MC_moves_variables_List = [ "DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                                      "CrankShaftFreq", "VolFreq", "MultiParticleFreq",]

        MEMC_MC_moves_variables_List = ["IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                                        "IntraMEMC-3Freq", "MEMC-3Freq",
                                        "ExchangeVolumeDim", "MEMC_DataInput"]

        valid_input_variables_List = sim_info_variables_List + CBMC_variables_List \
                                     + output_freq_variables_List + histogram_output_variables_List \
                                     + output_data_variables_List + free_energy_variables_List \
                                     + Std_MC_moves_variables_List + MEMC_MC_moves_variables_List


    else:
        print('ERROR: The ensemble_type selected for the _get_possible_ensemble_input_variables function is not valid.')
        valid_input_variables_List = None


    return  valid_input_variables_List




def _get_required_ensemble_files(ensemble_type):
    """
    Provides list of the possible optional input variables based on the ensemble type

    Parameters
    ----------
    ensemble_type = str; valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'

    Outputs
    ---------
    required_ensemble_files_List = list; a list of the required ensemble input variables

    """
    if ensemble_type in ['NVT', 'NPT']:

        simulation_settings = ["charmm_object", "ensemble_type", "RunSteps", "Temperature"]

        required_ensemble_files_List =  simulation_settings

    elif ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:

        simulation_settings = ["charmm_object", "ensemble_type", "RunSteps", "Temperature"]

        required_ensemble_files_List =   simulation_settings


    else:
        print('ERROR: The ensemble_type selected for the _get_required_ensemble_files function is not valid.')
        required_ensemble_files_List = None

    return required_ensemble_files_List




class GOMCControl():
    def __init__(self, charmm_object, ensemble_type, RunSteps, Temperature, input_variables_dict = None
                 ):

        """

        Construct the GOMC configuration input file with the defaults,
        or adding additional data in the input_variable section
        Default setting for the GOMC configuraion files are based upon the
        a educated guess which should result in reasonable sampling for a
        given ensemble/simulation type. However, there is no guarantee that
        the default setting will provide the best or adequate sampling for
        the selected system. The user has the option to modify the
        configuration/contorl files based on the simulation specifics or in to
        optimize the system beyond the standard settings.  These override
        options are available via the keyword arguments in input_variable_dict.
        Parameters
        ----------
        charmm_object :  Charmm object, which by definition has been parameterized ' \
                      + 'from the selected force field.',
        ensemble_typ : str, only accepts 'NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', 'GCMC'
        RunSteps : int; must be an integer greater than zero.
            Sets the total number of simulation steps.
        Temperature : Temperature of system in Kelvin (K)

        input_variables_dict: dict, default = None
            These input variables are optional and override the default settings.
            Changing these variables likely required for more advanced systems.

            The details of the acceptable input variables for the selected
            ensembles can be found by running this python workbook,

                print_valid_ensemble_input_variables('GCMC', description = True)

            which prints the input_variables with their subsection description
            for the selected 'GCMC' ensemble (other ensembles can be set as well).

            Example : input_variables_dict = {'Restart' : False, 'PRNG' : 123,
                                              'ParaTypeCHARMM' : True }

        Returns
        -------
        None

        Notes
        -------
        The details of the required inputs for the selected
        ensembles can be found by running this python workbook,

             print_valid_required_input_variables('NVT', description = True)

        which prints the required inputs with their subsection description
        for the selected 'NVT' ensemble (other ensembles can be set as well).

        The box units imported are in nm.  They need and are to be converted to Ang for GOMC or NAMD
        """

        # set this to check and see if all the input pass
        self.all_inputs_pass = True

        # set this to check and see if all the input pass
        self.all_failed_input_List = []

        # Check if charmm_object is really a Charmm() object
        if isinstance(charmm_object, mf_charmm.Charmm) == False:
            self.all_inputs_pass = False
            print_error_message = 'The variable supplied as a charmm_object ({}}) is not a ' \
                                  'charmm_object ({}})'.format(type(mf_charmm.Charmm), type(mf_charmm.Charmm))
            raise ValueError(print_error_message)

        # check ensemble is a correct type
        print('INFO: ensemble_type = ' + str(ensemble_type))
        if ensemble_type in ['NPT', 'NVT', 'GCMC', 'GEMC_NVT', 'GEMC_NPT']:
            self.ensemble_type = ensemble_type
        else:
            self.all_inputs_pass = False
            print_error_message = "The ensemble_type is not a valid entry. "\
                                  "Only 'NPT', 'NVT', 'GCMC', 'GEMC_NVT', 'GEMC_NPT' are valid entries."
            raise ValueError(print_error_message)

        self.RunSteps = RunSteps
        self.Temperature = Temperature
        if charmm_object.FF_filename != None and isinstance(charmm_object.FF_filename, str) == True:
            self.FF_filename = charmm_object.FF_filename
        elif charmm_object.FF_filename is None or isinstance(charmm_object.FF_filename, str) == False:
            self.all_inputs_pass = False
            print_error_message = "The force field file name was not specified and in the Charmm object ({}})." \
                                  "Therefore, the force field file (.inp) can not be written, and thus, the " \
                                  "GOMC control file (.conf) can not be created. Please use the force field file " \
                                  "name when building the Charmm object ({}})".format(type(mf_charmm.Charmm),
                                                                                      type(mf_charmm.Charmm))
            raise ValueError(print_error_message)

        if charmm_object.filename_box_0 != None and isinstance(charmm_object.filename_box_0, str) == True:
            self.Coordinates_box_0 = str(charmm_object.filename_box_0) + '.pdb'
            self.Structures_box_0 = str(charmm_object.filename_box_0) + '.psf'
        if charmm_object.filename_box_1 != None and isinstance(charmm_object.filename_box_1, str) == True:
            self.Coordinates_box_1 = str(charmm_object.filename_box_1) + '.pdb'
            self.Structures_box_1 = str(charmm_object.filename_box_1) + '.psf'
        else:
            self.Coordinates_box_1 = None
            self.Structures_box_1 = None

        self.coul_1_4_scaling = charmm_object.coul_1_4
        self.input_variables_dict = input_variables_dict
        self.residues_List = charmm_object.residues
        self.all_atom_names_and_res_pairs_dict = charmm_object.all_atom_name_res_pairs_dict
        self.all_atom_names_and_res_pairs_List = list(self.all_atom_names_and_res_pairs_dict.keys())

        self.x_dim_box_0 = charmm_object.box_0.maxs[0] * 10   # times 10 to convert from nm to Angstroms
        self.y_dim_box_0 = charmm_object.box_0.maxs[1] * 10   # times 10 to convert from nm to Angstroms
        self.z_dim_box_0 = charmm_object.box_0.maxs[2] * 10   # times 10 to convert from nm to Angstroms
        if charmm_object.filename_box_1 != None and isinstance(charmm_object.filename_box_1, str) == True:
            self.x_dim_box_1 = charmm_object.box_1.maxs[0] * 10   # times 10 to convert from nm to Angstroms
            self.y_dim_box_1 = charmm_object.box_1.maxs[1] * 10   # times 10 to convert from nm to Angstroms
            self.z_dim_box_1 = charmm_object.box_1.maxs[2] * 10   # times 10 to convert from nm to Angstroms
        else:
            self.x_dim_box_1 = None
            self.y_dim_box_1 = None
            self.z_dim_box_1 = None

        # the future control file name is entered now as None
        self.conf_filename = None

        # list of bad variable inputs
        bad_input_variables_values_List = []

        # set all the other variable initally to None (they are corrected to their set or default values later)
        default_input_variables_dict = _get_default_variables_dict()

        self.Restart = default_input_variables_dict['Restart']
        self.RestartCheckpoint = default_input_variables_dict['RestartCheckpoint']
        self.PRNG = default_input_variables_dict['PRNG']
        self.ParaTypeCHARMM = default_input_variables_dict['ParaTypeCHARMM']
        self.ParaTypeMie = default_input_variables_dict['ParaTypeMie']
        self.ParaTypeMARTINI = default_input_variables_dict['ParaTypeMARTINI']
        self.RcutCoulomb_box_0 = default_input_variables_dict['RcutCoulomb_box_0']
        self.RcutCoulomb_box_1 = default_input_variables_dict['RcutCoulomb_box_1']
        self.Pressure = default_input_variables_dict['Pressure']
        self.Rcut = default_input_variables_dict['Rcut']
        self.RcutLow = default_input_variables_dict['RcutLow']
        self.LRC = default_input_variables_dict['LRC']
        self.Exclude = default_input_variables_dict['Exclude']
        self.Potential = default_input_variables_dict['Potential']
        self.Rswitch = default_input_variables_dict['Rswitch']
        self.ElectroStatic = default_input_variables_dict['ElectroStatic']
        self.Ewald = default_input_variables_dict['Ewald']
        self.CachedFourier = default_input_variables_dict['CachedFourier']
        self.Tolerance = default_input_variables_dict['Tolerance']
        self.Dielectric = default_input_variables_dict['Dielectric']
        self.PressureCalc = default_input_variables_dict['PressureCalc']
        self.EqSteps = default_input_variables_dict['EqSteps']
        self.AdjSteps = default_input_variables_dict['AdjSteps']
        self.useConstantArea = default_input_variables_dict['useConstantArea']
        self.FixVolBox0 = default_input_variables_dict['FixVolBox0']
        self.ChemPot = default_input_variables_dict['ChemPot']
        self.Fugacity = default_input_variables_dict['Fugacity']
        self.CBMC_First = default_input_variables_dict['CBMC_First']
        self.CBMC_Nth = default_input_variables_dict['CBMC_Nth']
        self.CBMC_Ang = default_input_variables_dict['CBMC_Ang']
        self.CBMC_Dih = default_input_variables_dict['CBMC_Dih']
        self.OutputName = default_input_variables_dict['OutputName']
        self.CoordinatesFreq = default_input_variables_dict['CoordinatesFreq']
        self.RestartFreq = default_input_variables_dict['RestartFreq']
        self.CheckpointFreq = default_input_variables_dict['CheckpointFreq']
        self.ConsoleFreq = default_input_variables_dict['ConsoleFreq']
        self.BlockAverageFreq = default_input_variables_dict['BlockAverageFreq']
        self.HistogramFreq = default_input_variables_dict['HistogramFreq']
        self.DistName = default_input_variables_dict['DistName']
        self.HistName = default_input_variables_dict['HistName']
        self.RunNumber = default_input_variables_dict['RunNumber']
        self.RunLetter = default_input_variables_dict['RunLetter']
        self.SampleFreq = default_input_variables_dict['SampleFreq']
        self.OutEnergy = default_input_variables_dict['OutEnergy']
        self.OutPressure = default_input_variables_dict['OutPressure']
        self.OutMolNumber = default_input_variables_dict['OutMolNumber']
        self.OutDensity = default_input_variables_dict['OutDensity']
        self.OutVolume = default_input_variables_dict['OutVolume']
        self.OutSurfaceTension = default_input_variables_dict['OutSurfaceTension']
        self.FreeEnergyCalc = default_input_variables_dict['FreeEnergyCalc']
        self.MoleculeType = default_input_variables_dict['MoleculeType']
        self.InitialState = default_input_variables_dict['InitialState']
        self.LambdaVDW = default_input_variables_dict['LambdaVDW']
        self.LambdaCoulomb = default_input_variables_dict['LambdaCoulomb']
        self.ScaleCoulomb = default_input_variables_dict['ScaleCoulomb']
        self.ScalePower = default_input_variables_dict['ScalePower']
        self.ScaleAlpha = default_input_variables_dict['ScaleAlpha']
        self.MinSigma = default_input_variables_dict['MinSigma']
        self.DisFreq = default_input_variables_dict['DisFreq']
        self.RotFreq = default_input_variables_dict['RotFreq']
        self.IntraSwapFreq = default_input_variables_dict['IntraSwapFreq']
        self.SwapFreq = default_input_variables_dict['SwapFreq']
        self.RegrowthFreq = default_input_variables_dict['RegrowthFreq']
        self.CrankShaftFreq = default_input_variables_dict['CrankShaftFreq']
        self.VolFreq = default_input_variables_dict['VolFreq']
        self.MultiParticleFreq = default_input_variables_dict['MultiParticleFreq']
        # standard moves
        self.DisFreq = default_input_variables_dict['DisFreq'][self.ensemble_type]
        self.RotFreq = default_input_variables_dict['RotFreq'][self.ensemble_type]
        self.IntraSwapFreq = default_input_variables_dict['IntraSwapFreq'][self.ensemble_type]
        self.SwapFreq = default_input_variables_dict['SwapFreq'][self.ensemble_type]
        self.RegrowthFreq = default_input_variables_dict['RegrowthFreq'][self.ensemble_type]
        self.CrankShaftFreq = default_input_variables_dict['CrankShaftFreq'][self.ensemble_type]
        self.VolFreq = default_input_variables_dict['VolFreq'][self.ensemble_type]
        self.MultiParticleFreq = default_input_variables_dict['MultiParticleFreq'][self.ensemble_type]

        self.IntraMEMC_1Freq = default_input_variables_dict['IntraMEMC-1Freq'][self.ensemble_type]
        self.MEMC_1Freq = default_input_variables_dict['MEMC-1Freq'][self.ensemble_type]
        self.IntraMEMC_2Freq = default_input_variables_dict['IntraMEMC-2Freq'][self.ensemble_type]
        self.MEMC_2Freq = default_input_variables_dict['MEMC-2Freq'][self.ensemble_type]
        self.IntraMEMC_3Freq = default_input_variables_dict['IntraMEMC-3Freq'][self.ensemble_type]
        self.MEMC_3Freq = default_input_variables_dict['MEMC-3Freq'][self.ensemble_type]

        self.ExchangeVolumeDim = default_input_variables_dict['ExchangeVolumeDim']
        self.MEMC_DataInput = default_input_variables_dict['MEMC_DataInput']

        # auto calculate the best EqSteps (number of Equilbrium Steps) and Adj_Steps (number of AdjSteps Steps)
        set_max_steps_equib = self.EqSteps #1 * 10 ** 6

        if RunSteps / 10 >= set_max_steps_equib and RunSteps / 10 >= 1:
            self.EqSteps = int(set_max_steps_equib)
        elif RunSteps / 10  >= 1:
            self.EqSteps = int(RunSteps / 10)
        else:
            self.EqSteps = int( 1)

        set_max_steps_adj = self.AdjSteps #1000
        if RunSteps / 10 >= set_max_steps_adj and RunSteps / 10 >= 1:
            self.AdjSteps = int(set_max_steps_adj)
        elif RunSteps / 10 >= 1:
            self.AdjSteps = int(RunSteps / 10)
        else:
            self.AdjSteps = int(1)



        # auto calculate the best RestartFreq  for the number of RunSteps
        set_max_steps_RestartFreq  = self.RestartFreq[1] #1 * 10 ** 6

        if RunSteps / 10 >= set_max_steps_RestartFreq and RunSteps / 10 >= 1:
            self.RestartFreq[1] = int(set_max_steps_RestartFreq)
        elif RunSteps / 10  >= 1:
            self.RestartFreq[1] = int(RunSteps / 10)
        else:
            self.RestartFreq[1] = int( 1)

        # auto calculate the best CheckpointFreq  for the number of RunSteps
        set_max_steps_CheckpointFreq  = self.CheckpointFreq[1] # 1 * 10 ** 6

        if RunSteps / 10 >= set_max_steps_CheckpointFreq and RunSteps / 10 >= 1:
            self.CheckpointFreq[1] = int(set_max_steps_CheckpointFreq)
        elif RunSteps / 10  >= 1:
            self.CheckpointFreq[1] = int(RunSteps / 10)
        else:
            self.CheckpointFreq[1] = int( 1)

        # auto calculate the best CoordinatesFreq  for the number of RunSteps
        set_max_steps_CoordinatesFreq  = self.CoordinatesFreq[1] # 1 * 10 ** 6

        if RunSteps / 10 >= set_max_steps_CoordinatesFreq and RunSteps / 10 >= 1:
            self.CoordinatesFreq[1] = int(set_max_steps_CoordinatesFreq)
        elif RunSteps / 10  >= 1:
            self.CoordinatesFreq[1] = int(RunSteps / 10)
        else:
            self.CoordinatesFreq[1] = int( 1)

        # auto calculate the best ConsoleFreq  for the number of RunSteps
        set_max_steps_ConsoleFreq  = self.ConsoleFreq[1] #1 * 10 ** 4

        if RunSteps / 10 >= set_max_steps_ConsoleFreq and RunSteps / 10 >= 1:
            self.ConsoleFreq[1] = int(set_max_steps_ConsoleFreq)
        elif RunSteps / 10  >= 1:
            self.ConsoleFreq[1] = int(RunSteps / 10)
        else:
            self.ConsoleFreq[1] = int( 1)

        # auto calculate the best PressureCalc  for the number of RunSteps
        set_max_steps_PressureCalc  = self.PressureCalc[1]  #1 * 10 ** 4

        if RunSteps / 10 >= set_max_steps_PressureCalc and RunSteps / 10 >= 1:
            self.PressureCalc[1] = int(set_max_steps_PressureCalc)
        elif RunSteps / 10  >= 1:
            self.PressureCalc[1] = int(RunSteps / 10)
        else:
            self.PressureCalc[1] = int( 1)

        # auto calculate the best BlockAverageFreq  for the number of RunSteps
        set_max_steps_BlockAverageFreq  = self.BlockAverageFreq[1] #1 * 10 ** 4

        if RunSteps / 10 >= set_max_steps_BlockAverageFreq and RunSteps / 10 >= 1:
            self.BlockAverageFreq[1] = int(set_max_steps_BlockAverageFreq)
        elif RunSteps / 10  >= 1:
            self.BlockAverageFreq[1] = int(RunSteps / 10)
        else:
            self.BlockAverageFreq[1] = int( 1)

        # auto calculate the best HistogramFreq  for the number of RunSteps
        set_max_steps_HistogramFreq  = self.HistogramFreq[1] #1 * 10 ** 4

        if RunSteps / 10 >= set_max_steps_HistogramFreq and RunSteps / 10 >= 1:
            self.HistogramFreq[1] = int(set_max_steps_HistogramFreq)
        elif RunSteps / 10  >= 1:
            self.HistogramFreq[1] = int(RunSteps / 10)
        else:
            self.HistogramFreq[1] = int( 1)

        # auto calculate the best SampleFreq  for the number of RunSteps
        set_max_steps_SampleFreq  = self.SampleFreq # 500

        if RunSteps / 10 >= set_max_steps_SampleFreq and RunSteps / 10 >= 1:
            self.SampleFreq = int(set_max_steps_SampleFreq)
        elif RunSteps / 10  >= 1:
            self.SampleFreq = int(RunSteps / 10)
        else:
            self.SampleFreq = int( 1)



        # check the box dimensions
        if (isinstance( self.x_dim_box_0, int) == False \
            and isinstance( self.x_dim_box_0, float) == False) \
                or self.x_dim_box_0 <= 0:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The x-dimension for box 0 is, {} , not an integer, float, " \
                                  "or is <= 0.".format(self.x_dim_box_0)
            raise ValueError(print_error_message)

        if (isinstance( self.y_dim_box_0, int) == False \
            and isinstance( self.y_dim_box_0, float) == False) \
                or self.y_dim_box_0 <= 0:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The y-dimension for box 0 is, {} , not an integer, float, " \
                                  "or is <= 0.".format(self.y_dim_box_0)
            raise ValueError(print_error_message)

        if (isinstance( self.z_dim_box_0, int) == False \
            and isinstance( self.z_dim_box_0, float) == False) \
                or self.z_dim_box_0 <= 0:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The z-dimension for box 0 is, {} , not an integer, float, " \
                                  "or is <= 0.".format(self.z_dim_box_0)
            raise ValueError(print_error_message)

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            if (isinstance(self.x_dim_box_1, int) == False \
                and isinstance(self.x_dim_box_1, float) == False) \
                    or self.x_dim_box_1 <= 0:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The x-dimension for box 1 is, {} , not an integer, float, " \
                                      "or is <= 0.".format(self.x_dim_box_1)
                raise ValueError(print_error_message)

            elif self.x_dim_box_1 is None:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The x-dimension for box 1 was not provided.  The x-dimension for box 1, "\
                                      " is required to be an an integer, float, and be > 0."
                raise ValueError(print_error_message)

            if (isinstance(self.y_dim_box_1, int) == False \
                and isinstance(self.y_dim_box_1, float) == False) \
                    or self.y_dim_box_1 <= 0:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The y-dimension for box 1 is, {} , not an integer, float, " \
                                      "or is <= 0.".format(self.y_dim_box_1)
                raise ValueError(print_error_message)

            elif self.y_dim_box_1 is None:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The y-dimension for box 1 was not provided.  The y-dimension for box 1, "\
                                      "is required to be an an integer, float, and be > 0."
                raise ValueError(print_error_message)

            if (isinstance(self.z_dim_box_1, int) == False \
                and isinstance(self.z_dim_box_1, float) == False) \
                    or self.z_dim_box_1 <= 0:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The z-dimension for box 1 is, {} , not an integer, float, " \
                                      "or is <= 0.".format(self.z_dim_box_1)
                raise ValueError(print_error_message)

            elif self.z_dim_box_1 is None:
                self.all_inputs_pass = False
                print_error_message = "ERROR: The z-dimension for box 1 was not provided.  The z-dimension for box 1, "\
                                      "is required to be an an integer, float, and be > 0."
                raise ValueError(print_error_message)
        else:
            self.x_dim_box_1 = None
            self.y_dim_box_1 = None
            self.z_dim_box_1 = None



        # Checking for a valid ensemble type
        if self.ensemble_type not in ['NPT', 'NVT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The ensemble type selection of {}  is not a valid ensemble option. " \
                                  "Please choose the 'NPT', 'NVT', 'GEMC_NVT','GEMC_NPT', or 'GCMC' " \
                                  "ensembles".format(ensemble_type)
            raise ValueError(print_error_message)
        else:
            print("INFO: All the ensemble (ensemble_type) input passed the intial error checking")

        # check that the coulombic 1-4 scalar is : 0 =< 1-4 scalar <=1
        if (isinstance(self.coul_1_4_scaling, int) == False \
                and isinstance(self.coul_1_4_scaling, float) == False ) \
                or self.coul_1_4_scaling < 0 or self.coul_1_4_scaling > 1:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The selected 1-4 Coulombic scalar ({{) is not correct. "\
                                  "The 1-4 Coulombic scalar need to be an integer or float from " \
                                  "0 to 1.".format(self.coul_1_4_scaling)
            raise ValueError(print_error_message)

        # check that the Temperature is valid
        if self.Temperature <= 1:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The selected temperature ({}) is equal to or less than 1 Kelvin. "\
                                  "Please select a valid temperature".format(self.Temperature)
            raise ValueError(print_error_message)
        else:
            print("INFO: All the temperature  (Temperature) input passed the initial error checking")

        # RunSteps
        if not isinstance(self.RunSteps, int) or self.RunSteps <= 0:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The selected run steps (RunSteps variable = {}) is not "\
                                  "an integer or is less than or equal to 0.".format(self.RunSteps)
            raise ValueError(print_error_message)

        # create a list of the possible required files and check them based on the ensemble
        required_data_List = [self.FF_filename, self.Coordinates_box_0,
                               self.Structures_box_0]

        if self.Coordinates_box_1 != None:
            required_data_List.append(self.Coordinates_box_1)
        if self.Structures_box_1 != None:
            required_data_List.append(self.Structures_box_1)

        if self.ensemble_type in ['NVT', 'NPT']:
            if len(required_data_List) != 3 \
                    or os.path.splitext(self.FF_filename)[1] != '.inp' \
                    or os.path.splitext(self.Coordinates_box_0)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_0)[1] != '.psf':
                self.all_inputs_pass = False
                print_error_message = 'ERROR: The proper force field, PDB, and psf files were not provided, '\
                                      'or at least their extentions are not correct '\
                                      '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were '\
                                      'provided for the NVT or NPT simulations, which is not allowed'
                raise ValueError(print_error_message)

        else:
            print("INFO: All the required force field, pdb, and psf files for box 0 (.inp, .pdb, and .psf) all " \
                  + "passed the intial error checking. Note: the file names and their existance is not confirmed.")

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            if len(required_data_List) != 5 \
                    or os.path.splitext(self.FF_filename)[1] != '.inp' \
                    or os.path.splitext(self.Coordinates_box_0)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_0)[1] != '.psf' \
                    or os.path.splitext(self.Coordinates_box_1)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_1)[1] != '.psf':
                warn('ERROR: The proper force field, PDB, and psf files were not provided, '
                           'or at least their extentions are not correct '
                           '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were not provided '
                           'for the GEMC_NVT, GEMC_NPT or GCMC simulations, which is not allowed')
                self.all_inputs_pass = False
                print_error_message = 'ERROR: The proper force field, PDB, and psf files were not provided, '\
                                      'or at least their extentions are not correct '\
                                      '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were not provided '\
                                      'for the GEMC_NVT, GEMC_NPT or GCMC simulations, which is not allowed'
                raise ValueError(print_error_message)
        else:
            print("INFO: All the required force field, pdb, and psf files for box 0 and 1 (.inp, .pdb, and .psf) all "
                  "passed the intial error checking. Note: the file names and their existance is not confirmed.")

        if input_variables_dict is None:
            self.input_variables_dict = {}
        elif isinstance(input_variables_dict, dict) == True:
            self.input_variables_dict = input_variables_dict
        else:
            self.all_inputs_pass = False
            print_error_message = "ERROR: The input_variables_dict variable is not None or a dictionary. "
            raise ValueError(print_error_message)

        # verify all input variables keys are valid
        input_variables_dict_keys_list = dict_keys_to_list(self.input_variables_dict)
        if check_valid_ensemble_input_variables(self.ensemble_type, input_variables_dict_keys_list)[0] == False:
            returned_ck_bad_inputs_List = check_valid_ensemble_input_variables(self.ensemble_type,
                                                                               input_variables_dict_keys_list)[1]
            self.all_inputs_pass = False
            print_error_message = "ERROR: All the correct input variables where not provided for the {} " \
                                  "ensemble. Please be sure to check that the keys in the " \
                                  "input variables dictionary (input_variables_dict) is correct, and be aware "\
                                  "that added spaces before or after the variable in any keys " \
                                  "will also give this warning. The bad variable inputs " \
                                  "ensemble inputs = {}".format(self.ensemble_type, returned_ck_bad_inputs_List)
            raise ValueError(print_error_message)

        # verify all input variable values are valid, for their keys
        input_var_keys_list = dict_keys_to_list(self.input_variables_dict)


        possible_ensemble_variables_List = _get_possible_ensemble_input_variables(self.ensemble_type)


        # check to make sure the VDW FF (ParaTypeCHARMM) is set true  for multiple ones by the user
        # (i.e., ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI)
        VDW_ck_list = []
        for VDW_ck_inter in range(0,len(input_var_keys_list)):
            if input_var_keys_list[VDW_ck_inter] in ['ParaTypeCHARMM', 'ParaTypeMie', 'ParaTypeMARTINI'] \
                    and self.input_variables_dict[ input_var_keys_list[VDW_ck_inter]] == True:
                VDW_ck_list.append(True)

        if sum(VDW_ck_list) > 1:
            self.all_inputs_pass = False
            print_error_message = 'ERROR: there can only be 1 VDW set to true.  Please set only one of the '\
                                  'ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI types to True in the ' \
                                  'user variable input'
            raise ValueError(print_error_message)


        # check for MC move ratios and zero all of them if any are in the input variables
        for var_iter in range(0, len(input_var_keys_list)):
            # standard MC moves
            key_move_List = ["DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                             "CrankShaftFreq", "VolFreq", "MultiParticleFreq",
                             "IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                             "IntraMEMC-3Freq", "MEMC-3Freq"]


            if input_var_keys_list[var_iter] in key_move_List:
                self.DisFreq = 0.00
                self.RotFreq = 0.00
                self.IntraSwapFreq = 0.00
                self.SwapFreq = 0.00
                self.RegrowthFreq = 0.00
                self.CrankShaftFreq = 0.00
                self.VolFreq = 0.00
                self.MultiParticleFreq = 0.00
                self.IntraMEMC_1Freq = 0.00
                self.MEMC_1Freq = 0.00
                self.IntraMEMC_2Freq = 0.00
                self.MEMC_2Freq = 0.00
                self.IntraMEMC_3Freq = 0.00
                self.MEMC_3Freq = 0.00

            # set all the "RcutLow", "Rcut", "Rswitch" variable ahead of time so they can check the values
            # relative to each other in the next interation, regardless of their user entered order
            if input_var_keys_list[var_iter] == "Rcut":
                self.Rcut = self.input_variables_dict["Rcut"]

            if input_var_keys_list[var_iter] == "RcutLow":
                self.RcutLow = self.input_variables_dict["RcutLow"]

            if input_var_keys_list[var_iter] == "Rswitch":
                self.Rswitch = self.input_variables_dict["Rswitch"]


        # check for bad input variables and list the bad ones
        for var_iter in range(0, len(input_var_keys_list)):
            key = "Restart"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Restart = self.input_variables_dict[key]

            key = "RestartCheckpoint"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RestartCheckpoint = self.input_variables_dict[key]

            key = "PRNG"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != "RANDOM" \
                        and isinstance(self.input_variables_dict[key], int) != True:
                    bad_input_variables_values_List.append(key)
                if isinstance(self.input_variables_dict[key], int) == True:
                    if self.input_variables_dict[key] < 0 :
                        bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.PRNG = self.input_variables_dict[key]

            key = "ParaTypeCHARMM"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ParaTypeCHARMM = self.input_variables_dict[key]
                    if self.input_variables_dict[key] == True:
                        self.ParaTypeMie = False
                        self.ParaTypeMARTINI = False

            key = "ParaTypeMie"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ParaTypeMie = self.input_variables_dict[key]
                    if self.input_variables_dict[key] == True:
                        self.ParaTypeCHARMM = False

            key = "ParaTypeMARTINI"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ParaTypeMARTINI = self.input_variables_dict[key]
                    if self.input_variables_dict[key] == True:
                        self.ParaTypeCHARMM = False

            key = "RcutCoulomb_box_0"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RcutCoulomb_box_0 = self.input_variables_dict[key]

            key = "RcutCoulomb_box_1"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RcutCoulomb_box_1 = self.input_variables_dict[key]

            key = "Pressure"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Pressure = self.input_variables_dict[key]


            key = "Rcut"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Rcut = self.input_variables_dict[key]



            key = "RcutLow"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)
                if isinstance(self.input_variables_dict[key], float) \
                    or isinstance(self.input_variables_dict[key], int):
                    if self.input_variables_dict[key] >  self.Rcut:
                        bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RcutLow = self.input_variables_dict[key]


            key = "LRC"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.LRC = self.input_variables_dict[key]

            key = "Exclude"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != '1-2' and self.input_variables_dict[
                    key] != '1-3' \
                        and self.input_variables_dict[key] != '1-4':
                    bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Exclude = self.input_variables_dict[key]

            key = "Potential"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != 'VDW' and self.input_variables_dict[
                    key] != 'EXP6' \
                        and self.input_variables_dict[key] != 'SHIFT' and \
                        self.input_variables_dict[key] != 'SWITCH':
                    bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Potential = self.input_variables_dict[key]

            key = "Rswitch"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)
                if (isinstance(self.input_variables_dict[key], float) \
                        or isinstance(self.input_variables_dict[key], int)) \
                        and (isinstance(self.RcutLow, float) \
                        or isinstance(self.RcutLow, int)) \
                        and (isinstance(self.Rcut, float) \
                        or isinstance(self.Rcut, int)):
                    if self.input_variables_dict[key] <= self.RcutLow \
                            or self.input_variables_dict[key] >= self.Rcut:
                        bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Rswitch = self.input_variables_dict[key]

            key = "ElectroStatic"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ElectroStatic = self.input_variables_dict[key]

            key = "Ewald"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Ewald = self.input_variables_dict[key]

            key = "CachedFourier"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CachedFourier = self.input_variables_dict[key]

            key = "Tolerance"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_float_greater_zero_less_1(self.input_variables_dict,
                                                                            key,
                                                                            bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Tolerance = self.input_variables_dict[key]

            key = "Dielectric"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Dielectric = self.input_variables_dict[key]

            key = "PressureCalc"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_zero_or_greater(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.PressureCalc = self.input_variables_dict[key]

            key = "EqSteps"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                                   key,
                                                                   bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.EqSteps = self.input_variables_dict[key]

            key = "AdjSteps"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                                   key,
                                                                   bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.AdjSteps = self.input_variables_dict[key]

            key = "useConstantArea"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.useConstantArea = self.input_variables_dict[key]

            key = "FixVolBox0"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.FixVolBox0 = self.input_variables_dict[key]

            # ChemPot and Fugacity are only for GCMC
            key = "ChemPot"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_GCMC_dict_str_int_or_float(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ChemPot = self.input_variables_dict[key]

            key = "Fugacity"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_GCMC_dict_str_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                                  key,
                                                                                  bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.Fugacity = self.input_variables_dict[key]

            key = "CBMC_First"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CBMC_First = self.input_variables_dict[key]

            key = "CBMC_Nth"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CBMC_Nth = self.input_variables_dict[key]

            key = "CBMC_Ang"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CBMC_Ang = self.input_variables_dict[key]

            key = "CBMC_Dih"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CBMC_Dih = self.input_variables_dict[key]

            key = "OutputName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                                          key,
                                                                          bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutputName = self.input_variables_dict[key]

            key = "CoordinatesFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CoordinatesFreq = self.input_variables_dict[key]

            key = "RestartFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RestartFreq = self.input_variables_dict[key]

            key = "CheckpointFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CheckpointFreq = self.input_variables_dict[key]

            key = "ConsoleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ConsoleFreq = self.input_variables_dict[key]

            key = "BlockAverageFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.BlockAverageFreq = self.input_variables_dict[key]

            key = "HistogramFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                                key,
                                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.HistogramFreq = self.input_variables_dict[key]

            key = "DistName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.DistName = self.input_variables_dict[key]

            key = "HistName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.HistName = self.input_variables_dict[key]

            key = "RunNumber"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                                   key,
                                                                   bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RunNumber = self.input_variables_dict[key]

            key = "RunLetter"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], str) != True:
                    bad_input_variables_values_List.append(key)
                if isinstance(self.input_variables_dict[key], str) == True:
                    if len(self.input_variables_dict[key]) != 1 :
                        bad_input_variables_values_List.append(key)
                    elif len(self.input_variables_dict[key]) == 1 :
                        Is_Run_Letter_alphabet_char = self.input_variables_dict[key].isalpha()
                        if Is_Run_Letter_alphabet_char == False:
                            bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RunLetter = self.input_variables_dict[key]

            key = "SampleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                                   key,
                                                                   bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.SampleFreq = self.input_variables_dict[key]


            key = "OutEnergy"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutEnergy = self.input_variables_dict[key]

            key = "OutPressure"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutPressure = self.input_variables_dict[key]

            key = "OutMolNumber"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)
                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutMolNumber = self.input_variables_dict[key]

            key = "OutDensity"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutDensity = self.input_variables_dict[key]

            key = "OutVolume"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutVolume = self.input_variables_dict[key]

            key = "OutSurfaceTension"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.OutSurfaceTension = self.input_variables_dict[key]

            key = "FreeEnergyCalc"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.FreeEnergyCalc = self.input_variables_dict[key]

            key = "MoleculeType"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_residue_str_int_greater_zero(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)


                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MoleculeType = self.input_variables_dict[key]

            key = "InitialState"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.InitialState = self.input_variables_dict[key]

            key = "LambdaVDW"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_of_floats_zero_to_1(self.input_variables_dict,
                                                                           key,
                                                                           bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.LambdaVDW = self.input_variables_dict[key]

            key = "LambdaCoulomb"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_of_floats_zero_to_1(self.input_variables_dict,
                                                                           key,
                                                                           bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.LambdaCoulomb = self.input_variables_dict[key]

            key = "ScaleCoulomb"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_True_or_False(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ScaleCoulomb = self.input_variables_dict[key]

            key = "ScalePower"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                                      key,
                                                                      bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ScalePower = self.input_variables_dict[key]

            key = "ScaleAlpha"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ScaleAlpha = self.input_variables_dict[key]

            key = "MinSigma"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                               key,
                                                                               bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MinSigma = self.input_variables_dict[key]


            # standard MC moves
            key = "DisFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.DisFreq = self.input_variables_dict[key]

            key = "RotFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RotFreq = self.input_variables_dict[key]

            key = "IntraSwapFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.IntraSwapFreq = self.input_variables_dict[key]

            key = "SwapFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.SwapFreq = self.input_variables_dict[key]
                else:
                    self.VolFreq = 0.00

            key = "RegrowthFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.RegrowthFreq = self.input_variables_dict[key]


            key = "CrankShaftFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.CrankShaftFreq = self.input_variables_dict[key]

            key = "VolFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.VolFreq = self.input_variables_dict[key]
                else:
                    self.VolFreq = 0.00

            key = "MultiParticleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MultiParticleFreq = self.input_variables_dict[key]

            # MEMC moves freqencies
            key = "IntraMEMC-1Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.IntraMEMC_1Freq  = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_1Freq = 0.00

            key = "MEMC-1Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MEMC_1Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_1Freq = 0.00


            key = "IntraMEMC-2Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.IntraMEMC_2Freq = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_2Freq = 0.00


            key = "MEMC-2Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MEMC_2Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_2Freq = 0.00


            key = "IntraMEMC-3Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.IntraMEMC_3Freq = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_3Freq = 0.00

            key = "MEMC-3Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_List)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MEMC_3Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_3Freq = 0.00


            key = "ExchangeVolumeDim"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], list) == False:
                    bad_input_variables_values_List.append(key)
                elif isinstance(self.input_variables_dict[key], list) == True:
                    if len(self.input_variables_dict[key]) != 3 \
                            or (isinstance(self.input_variables_dict[key][0], float) != True \
                                and isinstance(self.input_variables_dict[key][0], int) != True) \
                            or (isinstance(self.input_variables_dict[key][1], float) != True \
                                and isinstance(self.input_variables_dict[key][1], int) != True) \
                            or (isinstance(self.input_variables_dict[key][2], float) != True \
                                and isinstance(self.input_variables_dict[key][2], int) != True) \
                            or str(self.input_variables_dict[key][0]) == str(True) \
                                or str(self.input_variables_dict[key][0]) == str(False) \
                            or str(self.input_variables_dict[key][1]) == str(True) \
                                or str(self.input_variables_dict[key][1]) == str(False) \
                            or str(self.input_variables_dict[key][2]) == str(True) \
                                or str(self.input_variables_dict[key][2]) == str(False) :
                        bad_input_variables_values_List.append(key)
                    elif len(self.input_variables_dict[key]) == 3 :
                        if (isinstance(self.input_variables_dict[key][0], float) == True \
                            or isinstance(self.input_variables_dict[key][0], int) == True) \
                                and (isinstance(self.input_variables_dict[key][1], float) == True \
                                    or isinstance(self.input_variables_dict[key][1], int) == True) \
                                and (isinstance(self.input_variables_dict[key][2], float) == True \
                                    or isinstance(self.input_variables_dict[key][2], int) == True):
                                if self.input_variables_dict[key][0] <= 0 \
                                    or self.input_variables_dict[key][1] <= 0  \
                                    or self.input_variables_dict[key][2] <= 0 :
                                    bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.ExchangeVolumeDim = self.input_variables_dict[key]

            key = "MEMC_DataInput"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], list) == False:
                    bad_input_variables_values_List.append(key)
                elif isinstance(self.input_variables_dict[key], list) == True:
                    if len(self.input_variables_dict[key]) == 0 :
                        bad_input_variables_values_List.append(key)
                    elif len(self.input_variables_dict[key]) > 0:
                        No_MEMC_combos = len(self.input_variables_dict[key])
                        for MEMC_iter in range(0, No_MEMC_combos):
                            if isinstance(self.input_variables_dict[key][MEMC_iter], list) == False:
                                bad_input_variables_values_List.append(key)
                            elif isinstance(self.input_variables_dict[key][MEMC_iter], list) == True \
                                    and len(self.input_variables_dict[key][MEMC_iter]) == 5 :
                                if  isinstance(self.input_variables_dict[key][MEMC_iter][2], list) == False \
                                        or isinstance(self.input_variables_dict[key][MEMC_iter][4], list) == False:
                                    bad_input_variables_values_List.append(key)
                                elif isinstance(self.input_variables_dict[key][MEMC_iter][2], list) == True \
                                        and isinstance(self.input_variables_dict[key][MEMC_iter][4], list) == True:
                                    if (len(self.input_variables_dict[key][MEMC_iter][2]) != 2 \
                                            or len(self.input_variables_dict[key][MEMC_iter][4]) != 2):
                                        bad_input_variables_values_List.append(key)
                                    elif (len(self.input_variables_dict[key][MEMC_iter][2]) == 2 \
                                            and len(self.input_variables_dict[key][MEMC_iter][4]) == 2):
                                        if isinstance(self.input_variables_dict[key][MEMC_iter][0], int) != True \
                                                or str(self.input_variables_dict[key][MEMC_iter][0]) == str(True) \
                                                or str(self.input_variables_dict[key][MEMC_iter][0]) == str(False) \
                                                or isinstance(self.input_variables_dict[key][MEMC_iter][1], str) == False \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][2][0], str) == False \
                                                    and self.input_variables_dict[key][MEMC_iter][2][0] != None ) \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][2][1], str) == False \
                                                and self.input_variables_dict[key][MEMC_iter][2][1] != None) \
                                                or isinstance(self.input_variables_dict[key][MEMC_iter][3], str) == False \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][4][0], str) == False \
                                                and self.input_variables_dict[key][MEMC_iter][4][0] != None) \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][4][1], str) == False \
                                                and self.input_variables_dict[key][MEMC_iter][4][1] != None):
                                            bad_input_variables_values_List.append(key)
                                    else:
                                        bad_input_variables_values_List.append(key)

                                    # check that the atom names match the residues that exist
                                    if self.input_variables_dict[key][MEMC_iter][1] not in \
                                            self.all_atom_names_and_res_pairs_List:
                                        bad_input_variables_values_List.append(key)

                                    elif self.input_variables_dict[key][MEMC_iter][1] in \
                                        self.all_atom_names_and_res_pairs_List:

                                        if self.input_variables_dict[key][MEMC_iter][2][0] not in \
                                                self.all_atom_names_and_res_pairs_dict[
                                                    self.input_variables_dict[key][MEMC_iter][1]]:

                                            bad_input_variables_values_List.append(key)

                                        if self.input_variables_dict[key][MEMC_iter][2][1] not in \
                                                self.all_atom_names_and_res_pairs_dict[
                                                    self.input_variables_dict[key][MEMC_iter][1]]:
                                            bad_input_variables_values_List.append(key)

                                    if self.input_variables_dict[key][MEMC_iter][3] not in \
                                            self.all_atom_names_and_res_pairs_List:
                                        bad_input_variables_values_List.append(key)

                                    elif self.input_variables_dict[key][MEMC_iter][3] in \
                                            self.all_atom_names_and_res_pairs_List:
                                        if self.input_variables_dict[key][MEMC_iter][4][0] not in \
                                                self.all_atom_names_and_res_pairs_dict[
                                                    self.input_variables_dict[key][MEMC_iter][3]]:
                                            bad_input_variables_values_List.append(key)

                                        if self.input_variables_dict[key][MEMC_iter][4][1] not in \
                                                self.all_atom_names_and_res_pairs_dict[
                                                    self.input_variables_dict[key][MEMC_iter][3]]:
                                            bad_input_variables_values_List.append(key)

                                    if isinstance(self.input_variables_dict[key][MEMC_iter][0], int) == True:
                                        if self.input_variables_dict[key][MEMC_iter][0] <= 0:
                                            bad_input_variables_values_List.append(key)

                                else:
                                    bad_input_variables_values_List.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_List:
                    self.MEMC_DataInput = self.input_variables_dict[key]

        #Error out and print the bad input values
        if len(bad_input_variables_values_List) > 0:
            self.all_inputs_pass = False
            # create unique list
            bad_input_variables_values_set = set(bad_input_variables_values_List)
            bad_unique_input_variables_values_List = list(bad_input_variables_values_set)
            print_error_message = 'ERROR: The following input variables have ' \
                                  'bad values: {}'.format(bad_unique_input_variables_values_List)
            raise ValueError(print_error_message)

        else:
            print("INFO: All the input variable passed the intial error checking")


        # check to make sure the VDW FF (ParaTypeCHARMM) is not true for multiple ones
        # (i.e., ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI)
        if (self.ParaTypeCHARMM == True and self.ParaTypeMie == True) \
                or (self.ParaTypeCHARMM == True and self.ParaTypeMARTINI == True) \
                or self.ParaTypeMie == True and self.ParaTypeMARTINI == True :
            self.all_inputs_pass = False
            print_error_message = 'ERROR: there can only be 1 VDW type set to true.  Please set only one of the '\
                                  'ParaTypeCHARMM = {}, ParaTypeMie = {}, ParaTypeMARTINI = {} types ' \
                                  'to True'.format(self.ParaTypeCHARMM, self.ParaTypeMie, self.ParaTypeMARTINI)
            raise ValueError(print_error_message)
        elif self.ParaTypeCHARMM == True or self.ParaTypeMie == True \
                or self.ParaTypeMARTINI == True:
            if self.ParaTypeCHARMM == True:
                self.VDW_type = "ParaTypeCHARMM"
            elif self.ParaTypeMie == True:
                self.VDW_type = "ParaTypeMie"
            elif self.ParaTypeMARTINI == True:
                self.VDW_type = "ParaTypeMARTINI"
        else:
            self.all_inputs_pass = False
            print_error_message = 'ERROR: There no VDW types that are set as True.  Please set only one of the '\
                                  'ParaTypeCHARMM = {}, ParaTypeMie = {}, ParaTypeMARTINI = {} types ' \
                                  'to True'.format(self.ParaTypeCHARMM, self.ParaTypeMie, self.ParaTypeMARTINI)
            raise ValueError(print_error_message)



        # check to see if the moves sum up to 1
        if ensemble_type in ['NVT', 'GCMC']:
            if self.VolFreq != 0:
                self.all_inputs_pass = False
                print_error_message = 'ERROR: The input variable VolFreq is non-zero (0). '\
                                      'VolFreq must be zero (0) for the "NVT", "GEMC_NVT", and "GCMC" ensembles.'
                raise ValueError(print_error_message)

        if ensemble_type in ['NVT', 'NPT']:
            if self.SwapFreq != 0 or self.MEMC_1Freq != 0 or self.MEMC_2Freq !=0 or  self.MEMC_3Freq !=0:
                self.all_inputs_pass = False
                print_error_message = 'ERROR: All the MC move input variables must be non-zero (0) for the '\
                                      'SwapFreq, MEMC_1Freq, MEMC_2Freq, and MEMC_3Freq. ' \
                                      'The SwapFreq, MEMC_1Freq, MEMC_2Freq, '\
                                      'and MEMC_3Freq need to be zero (0) for the "NVT" and "NPT" ensembles.'
                raise ValueError(print_error_message)

        moves_list = [self.DisFreq, self.RotFreq, self.IntraSwapFreq,
                      self.SwapFreq, self.RegrowthFreq, self.CrankShaftFreq,
                      self.VolFreq, self.MultiParticleFreq,
                      self.IntraMEMC_1Freq, self.MEMC_1Freq,
                      self.IntraMEMC_2Freq, self.MEMC_2Freq,
                      self.IntraMEMC_3Freq, self.MEMC_3Freq
                      ]

        if sum(moves_list) <= 1 + 10**(-13) and  sum(moves_list) >= 1 - 10**(-13):
            print('INFO: The sum of the Monte Carlo move ratios = ' +  str("{:.12f}".format(sum(moves_list))))

        else:
            print('INFO: sum(moves_list) = ' + str("{:.12f}".format(sum(moves_list))))
            print('\t DisFreq = ' + str(self.DisFreq))
            print('\t RotFreq = ' + str(self.RotFreq))
            print('\t IntraSwapFreq = ' + str(self.IntraSwapFreq))
            print('\t SwapFreq = ' + str(self.SwapFreq))
            print('\t RegrowthFreq = ' + str(self.RegrowthFreq))
            print('\t CrankShaftFreq = ' + str(self.CrankShaftFreq))
            print('\t VolFreq = ' + str(self.VolFreq))
            print('\t MultiParticleFreq = ' + str(self.MultiParticleFreq))
            print('\t IntraMEMC_1Freq = ' + str(self.IntraMEMC_1Freq))
            print('\t MEMC_1Freq = ' + str(self.MEMC_1Freq))
            print('\t IntraMEMC_2Freq = ' + str(self.IntraMEMC_2Freq))
            print('\t MEMC_2Freq = ' + str(self.MEMC_2Freq))
            print('\t IntraMEMC_3Freq = ' + str(self.IntraMEMC_3Freq))
            print('\t MEMC_3Freq = ' + str(self.MEMC_3Freq))
            self.all_inputs_pass = False
            print_error_message = 'ERROR: The sum of the Monte Carlo move ratios does not equal 1. '\
                                  'Note: The sum that was manually entered may equal 1, but some ' \
                                  'moves may not be valid for the provided ensemble. The moves ' \
                                  'that are invalid for a given ensemble are set to zero. If ' \
                                  'the default moves are not being used, all the move frequencies ' \
                                  'which do not have default values of zero will need to be set manually '\
                                  'so the sum equals (DisFreq, RotFreq, IntraSwapFreq, SwapFreq, '\
                                  'RegrowthFreq, CrankShaftFreq, and VolFreq).'
            raise ValueError(print_error_message)

        # Check that RunSteps >= EqSteps >= AdjSteps
        if self.RunSteps < self.EqSteps or self.RunSteps < self.AdjSteps \
                or self.EqSteps < self.AdjSteps:
            self.all_inputs_pass = False
            print_error_message = 'ERROR: The values must be in this order RunSteps >= EqSteps >= AdjSteps ' \
                                  ' ({} >= {} >= {} )'.format(self.RunSteps, self.EqSteps, self.AdjSteps)
            raise ValueError(print_error_message)

        # check if both the ChemPot and Fugacity are not set to None.  Only one can be used
        if self.Fugacity != None and self.ChemPot != None \
                and self.ensemble_type == 'GCMC':
            self.all_inputs_pass = False
            print_error_message = 'ERROR:  In the GCMC ensemble, both Fugacity and ChemPot are provided. '\
                                  'Add a dictionary for either the Fugacity or ChemPot and set the other ' \
                                  'variable to None. Note: Both the Fugacity or ChemPot and set to None by default'
            raise ValueError(print_error_message)

        # check that both the ChemPot and Fugacity are set to None.  Only one can be used
        if self.Fugacity is None and self.ChemPot is None \
                and self.ensemble_type == 'GCMC':
            warn('ERROR: In the GCMC ensemble, neither Fugacity and ChemPot are provided (i.e., both are None). '
                       'Add a dictionary for either the Fugacity or ChemPot and set the other variable to None. ' 
                       'Note: Both the Fugacity or ChemPot and set to None by default')
            self.all_inputs_pass = False
            print_error_message = 'ERROR: In the GCMC ensemble, neither Fugacity and ChemPot are provided ' \
                                  '(i.e., both are None). Add a dictionary for either the Fugacity or ' \
                                  'ChemPot and set the other variable to None. ' \
                                  'Note: Both the Fugacity or ChemPot and set to None by default'
            raise ValueError(print_error_message)


        # check that MEMC moves rations are > 0 if MEMC_DataInput is used
        if self.MEMC_DataInput != None and (    self.MEMC_1Freq == 0 \
                                                           and self.MEMC_2Freq == 0 \
                                                           and self.MEMC_3Freq == 0 \
                                                           and self.IntraMEMC_1Freq == 0 \
                                                           and self.IntraMEMC_2Freq == 0 \
                                                           and self.IntraMEMC_3Freq == 0):
            self.all_inputs_pass = False
            print_error_message = 'ERROR: The MEMC_DataInput variable is not equal to None, ' \
                                  'but all the MEMC move ratios are '\
                                  'zero (IntraMEMC_1Freq, MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, '\
                                  'IntraMEMC_3Freq, and MEMC_3Freq).'
            raise ValueError(print_error_message)

        if self.MEMC_DataInput is None and (    self.MEMC_1Freq != 0 \
                                                           or self.MEMC_2Freq != 0 \
                                                           or self.MEMC_3Freq != 0 \
                                                           or self.IntraMEMC_1Freq != 0 \
                                                           or self.IntraMEMC_2Freq != 0 \
                                                           or self.IntraMEMC_3Freq != 0):
            self.all_inputs_pass = False
            print_error_message = 'ERROR: The MEMC_DataInput variable is equal to None, ' \
                                  'but at least one of the MEMC move ratios are '\
                                  'all non-zero (IntraMEMC_1Freq, MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, '\
                                  'IntraMEMC_3Freq, and MEMC_3Freq).'
            raise ValueError(print_error_message)

        # ensure the LargeKindBackBone and SmallKindBackBones are provided as appropriate for MEMC-1, MEMC-2, MEMC-3
        if self.MEMC_DataInput != None and (self.MEMC_2Freq > 0
                                                       or self.IntraMEMC_2Freq > 0):
            for MEMC_2_i in range(0, len(self.MEMC_DataInput)):
                if self.MEMC_DataInput[MEMC_2_i][2][0] is None \
                    or self.MEMC_DataInput[MEMC_2_i][2][1] is None \
                        or self.MEMC_DataInput[MEMC_2_i][4][0] is None \
                        or self.MEMC_DataInput[MEMC_2_i][4][1] is None :

                    self.all_inputs_pass = False
                    print_error_message = 'ERROR:  The  LargeKindBackBone and SmallKindBackBones unique ' \
                                          'atom names, strings, both must be provided when using the ' \
                                          'IntraMEMC-2Freq or MEMC-2Freq moves ' \
                                          '(i.e., the LargeKindBackBone and SmallKindBackBones can not be None). '
                    raise ValueError(print_error_message)

        if self.MEMC_DataInput != None and (self.MEMC_3Freq > 0
                                                       or self.IntraMEMC_3Freq > 0):
            for MEMC_3_i in range(0, len(self.MEMC_DataInput)):
                if self.MEMC_DataInput[MEMC_3_i][2][0] is None \
                    or self.MEMC_DataInput[MEMC_3_i][2][1] is None:
                    self.all_inputs_pass = False
                    print_error_message = 'ERROR:  The LargeKindBackBone unique atom names, strings, ' \
                                          'both must be provided when using the IntraMEMC-3Freq or MEMC-3Freq moves ' \
                                          '(i.e., the LargeKindBackBone can not be None).'
                    raise ValueError(print_error_message)

        # check that all required free energy values are provided
        if (self.FreeEnergyCalc != None or self.MoleculeType != None \
                or self.InitialState != None or self.LambdaVDW != None ) \
                and (self.FreeEnergyCalc is None or self.MoleculeType is None \
                or self.InitialState is None or self.LambdaVDW is None ):
            self.all_inputs_pass = False
            print_error_message = 'ERROR: To utilize the free energy calculations all the following ' \
                                  'variables need to be set, and not equal to None: ' \
                                  'FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW.'
            raise ValueError(print_error_message)


        if self.LambdaVDW != None and self.LambdaCoulomb != None \
                and isinstance(self.LambdaVDW, list) == True \
                     and (isinstance(self.LambdaCoulomb, list)) == True:
                if len(self.LambdaVDW) == len(self.LambdaCoulomb) :
                    if self.InitialState + 1 <= len(self.LambdaVDW):
                        for lam_i in range(1, len(self.LambdaVDW)):
                            if self.LambdaVDW[lam_i] < self.LambdaVDW[lam_i -1 ]:
                                self.all_inputs_pass = False
                                print_error_message = 'ERROR: The LambdaVDW list is not in accending order.'
                                raise ValueError(print_error_message)
                            if self.LambdaCoulomb[lam_i] < self.LambdaCoulomb[lam_i -1 ]:
                                self.all_inputs_pass = False
                                print_error_message = 'ERROR:  The LambdaCoulomb list is not in accending order.'
                                raise ValueError(print_error_message)
                    else:
                        self.all_inputs_pass = False
                        print_error_message = 'ERROR: The InitialState integer is greater than the LambdaVDW and '\
                                              'LambdaCoulomb list length.  Note: the InitialState integer starts at 0.'
                        raise ValueError(print_error_message)
                else:
                    self.all_inputs_pass = False
                    print_error_message = 'ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length.'
                    raise ValueError(print_error_message)




    # write the control file
    def write_conf_file(self, conf_filename):
        """
        Writes the GOMC control file

        Parameters
        ----------
        conf_filename; str,
            The path and file name for the control file name, with
            the .conf extension, or no extension.  If no extension is provided, the
            code will add the .conf extension to the provided file name.

        Outputs
        ---------
        Writes the GOMC control file with the name provided via conf_filename

        Returns
        ---------
        If completed without errors; str; "PASSED
        If completed with errors;  None
        """

        # check to see if it is OK to proceed writing the control file
        if self.all_inputs_pass == False:
            print_error_message = 'ERROR: The control file was not written as at least 1 input to the '\
                                  'control file writer was bad.'
            raise ValueError(print_error_message)

        date_time = datetime.datetime.today()

        self.conf_filename = conf_filename

        if isinstance(self.conf_filename, str) == False \
                or isinstance(self.conf_filename, str) is None:
            self.all_inputs_pass = False
            print_error_message = 'ERROR: The control file name (conf_filename) is not provided as a string. '
            raise ValueError(print_error_message)

        if os.path.splitext(self.conf_filename)[1] == '.conf':
            self.conf_filename = conf_filename
            print('INFO: the correct extension for the control file was provided in the file name, .conf '
                  'with control file name = ' + str(self.conf_filename))
        elif os.path.splitext(self.conf_filename)[1] == '':
            self.conf_filename = self.conf_filename + '.conf'
            print('INFO: No extension name was provided for the control file. Therefore, the proper '
                  'extension, .conf, was added.  The new total control file name = ' \
                  + str(self.conf_filename))
        else:
            self.all_inputs_pass = False
            print_error_message = 'ERROR: No extension name or the wrong extension name was provided. ' \
                                  'Please enter a proper extension name, .conf or no extension in the conf_filename '\
                                  'The control file as provided name = ' + str(self.conf_filename)
            raise ValueError(print_error_message)


        # get and open the control file file
        data_control_file = open(self.conf_filename, 'w')

        data_control_file.write('#######################################################'
                                '#########################################\n')
        data_control_file.write("##  This file (" + self.conf_filename \
                                + ') - was created by mBuild using the on ' + str(date_time) + '\n')
        data_control_file.write('#######################################################'
                                '#########################################\n')
        data_control_file.write('\n')
        data_control_file.write('############################################################################\n')
        data_control_file.write('#  ========-------------------- INPUT --------------------------=========== \n')
        data_control_file.write('############################################################################\n')
        data_control_file.write(' \n')
        data_control_file.write('#########################\n')
        data_control_file.write('# enable, step\n')
        data_control_file.write('#########################\n')
        data_control_file.write('Restart \t ' +str(self.Restart)  + '\n')
        data_control_file.write('\n')
        data_control_file.write('####################################\n')
        data_control_file.write('# kind {RESTART, RANDOM, INTSEED}\n')
        data_control_file.write('####################################\n')
        if self.PRNG == "RANDOM":
            data_control_file.write('PRNG \t\t ' + str(self.PRNG ) + '\n')
        elif isinstance(self.PRNG,int):
            data_control_file.write('PRNG \t\t ' + 'INTSEED \n')
            data_control_file.write('Random_Seed \t ' + str(self.PRNG) + '\n')
        data_control_file.write(' \n')
        data_control_file.write('####################################\n')
        data_control_file.write('# FORCE FIELD\n')
        data_control_file.write('####################################\n')
        data_control_file.write(str(self.VDW_type) +'\t\t ' + str(True) + '\n')
        data_control_file.write(' \n')
        data_control_file.write('Parameters \t\t ' + str(self.FF_filename) + '\n')
        data_control_file.write('####################################\n')
        data_control_file.write('# INPUT PDB FILES\n')
        data_control_file.write('####################################\n')
        if self.ensemble_type in ['NVT', 'NPT', 'GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Coordinates 0 \t\t ' + str(self.Coordinates_box_0) + '\n')
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Coordinates 1 \t\t ' + str(self.Coordinates_box_1) + '\n')
        data_control_file.write(' \n')
        data_control_file.write('####################################\n')
        data_control_file.write('# INPUT PSF FILES\n')
        data_control_file.write('####################################\n')
        if self.ensemble_type in ['NVT', 'NPT', 'GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Structure 0 \t\t ' + str(self.Structures_box_0) + '\n')
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Structure 1 \t\t ' + str(self.Structures_box_1) + '\n')
        data_control_file.write(' \n')
        data_control_file.write('############################################################################\n')
        data_control_file.write('#  =======--------------------- SYSTEM --------------------------===========\n')
        data_control_file.write('############################################################################ \n')
        data_control_file.write(' \n')
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT']:
            data_control_file.write('##################################\n')
            data_control_file.write('# GEMC TYPE \n')
            data_control_file.write('##################################\n')
            if self.ensemble_type in 'GEMC_NPT':
                data_control_file.write('GEMC \t\t NPT \n')
            elif self.ensemble_type in 'GEMC_NVT':
                data_control_file.write('GEMC \t\t NVT \n')
        data_control_file.write(' \n')
        data_control_file.write('#############################\n')
        data_control_file.write('# SIMULATION CONDITION\n')
        data_control_file.write('#############################\n')
        if self.ensemble_type in ['GEMC_NPT', 'NPT']:
            data_control_file.write('Pressure \t ' + str(self.Pressure) + '\n')
        data_control_file.write('Temperature \t ' + str(self.Temperature) + '\n')

        if self.ensemble_type in ['GCMC'] and self.ChemPot != None \
                and self.Fugacity is None:
            chem_pot_dict_key_List = dict_keys_to_list(self.ChemPot)
            for chem_pot_iter in range(0, len(chem_pot_dict_key_List)):
                chem_pot_residue_iter = chem_pot_dict_key_List[chem_pot_iter]
                data_control_file.write('ChemPot \t ' + str(chem_pot_residue_iter) + ' \t\t ' \
                               + str(self.ChemPot[chem_pot_residue_iter]) + '\n')

        if self.ensemble_type in ['GCMC'] and self.Fugacity != None \
                and self.ChemPot is None :
            fugacity_iter_dict_key_List = dict_keys_to_list(self.Fugacity)
            for fugacity_iter in range(0, len(fugacity_iter_dict_key_List)):
                fugacity_residue_iter = fugacity_iter_dict_key_List[fugacity_iter]
                data_control_file.write('Fugacity \t ' + str(fugacity_residue_iter) + ' \t\t ' \
                               + str(self.Fugacity[fugacity_residue_iter]) + '\n')

        data_control_file.write(' \n')
        data_control_file.write('Potential \t ' + str(self.Potential) + '\n')
        data_control_file.write('LRC \t\t ' + str(self.LRC) + '\n')
        data_control_file.write('Rcut \t\t ' + str(self.Rcut) + '\n')
        data_control_file.write('RcutLow \t ' + str(self.RcutLow) + '\n')
        if self.Potential == 'SWITCH':
            data_control_file.write('Rswitch \t ' + str(self.Rswitch) + '\n')
        data_control_file.write('Exclude \t ' + str(self.Exclude) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('#############################\n')
        data_control_file.write('# ELECTROSTATIC   \n')
        data_control_file.write('#############################\n')
        data_control_file.write('Ewald \t\t ' + str(self.Ewald) + '\n')
        data_control_file.write('ElectroStatic \t ' + str(self.ElectroStatic) + '\n')
        data_control_file.write('CachedFourier \t ' + str(self.CachedFourier) + '\n')
        data_control_file.write('Tolerance \t ' + str(format(self.Tolerance, '.12f')) + '\n')
        if self.VDW_type in ["ParaTypeMARTINI"]:
            data_control_file.write('Dielectric \t ' + str(format(self.Dielectric)) + '\n')
        data_control_file.write('1-4scaling \t ' + str(self.coul_1_4_scaling) + '\n')
        data_control_file.write(' \n')
        if self.RcutCoulomb_box_0 != None:
            data_control_file.write('RcutCoulomb 0 \t ' + str(self.RcutCoulomb_box_0) + '\n')
        if self.RcutCoulomb_box_0 != None \
                and self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('RcutCoulomb 1 \t ' + str(self.RcutCoulomb_box_1) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# PRESSURE CALCULATION\n')
        data_control_file.write('################################\n')
        data_control_file.write('PressureCalc \t ' + str(self.PressureCalc[0]) + ' \t\t ' \
                               + str(self.PressureCalc[1]) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# STEPS \n')
        data_control_file.write('################################\n')
        data_control_file.write('RunSteps \t ' + str(self.RunSteps) + '\n')
        data_control_file.write('EqSteps \t ' + str(self.EqSteps) + '\n')
        data_control_file.write('AdjSteps \t ' + str(self.AdjSteps) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# MOVE FREQUENCY \n')
        data_control_file.write('################################\n')
        data_control_file.write('DisFreq \t\t ' \
                                + str(self.DisFreq) + '\n')
        data_control_file.write('RotFreq \t\t ' \
                                + str(self.RotFreq) + '\n')
        data_control_file.write('IntraSwapFreq \t\t ' \
                                + str(self.IntraSwapFreq) + '\n')
        data_control_file.write('SwapFreq \t\t ' \
                                + str(self.SwapFreq) + '\n')
        data_control_file.write('RegrowthFreq \t\t ' \
                                + str(self.RegrowthFreq) + '\n')
        data_control_file.write('CrankShaftFreq \t\t ' \
                                + str(self.CrankShaftFreq) + '\n')
        data_control_file.write('VolFreq \t\t ' \
                                + str(self.VolFreq) + '\n')
        data_control_file.write('MultiParticleFreq \t ' \
                                + str(self.MultiParticleFreq) + '\n')
        data_control_file.write('IntraMEMC-1Freq \t ' \
                                + str(self.IntraMEMC_1Freq) + '\n')
        data_control_file.write('MEMC-1Freq \t\t ' \
                                + str(self.MEMC_1Freq) + '\n')
        data_control_file.write('IntraMEMC-2Freq \t ' \
                                + str(self.IntraMEMC_2Freq) + '\n')
        data_control_file.write('MEMC-2Freq \t\t ' \
                                + str(self.MEMC_2Freq) + '\n')
        data_control_file.write('IntraMEMC-3Freq \t ' \
                                + str(self.IntraMEMC_3Freq) + '\n')
        data_control_file.write('MEMC-3Freq \t\t ' \
                                + str(self.MEMC_3Freq) + '\n')
        data_control_file.write(' \n')

        # sort and print the MEMC data if MEMC is used for the simulation
        if self.MEMC_DataInput != None and (    self.MEMC_1Freq > 0 \
                                                           or self.MEMC_2Freq > 0 \
                                                           or self.MEMC_3Freq > 0 \
                                                           or self.IntraMEMC_1Freq > 0 \
                                                           or self.IntraMEMC_2Freq > 0 \
                                                           or self.IntraMEMC_3Freq > 0):

            ExchangeRatio_List = []
            ExchangeLargeKind_List = []
            LargeKindBackBone_List = []
            ExchangeSmallKind_List = []
            SmallKindBackBone_List = []

            for memc_i in range(0, len(self.MEMC_DataInput)):
                ExchangeRatio_List.append(str(self.MEMC_DataInput[memc_i][0]) + " \t\t")
                ExchangeLargeKind_List.append(str(self.MEMC_DataInput[memc_i][1]) + " \t\t")
                LargeKindBackBone_List.append(str(self.MEMC_DataInput[memc_i][2][0]) + "   " \
                                              + str(self.MEMC_DataInput[memc_i][2][1]) + " \t")
                ExchangeSmallKind_List.append(str(self.MEMC_DataInput[memc_i][3]) + " \t\t")
                SmallKindBackBone_List.append(str(self.MEMC_DataInput[memc_i][4][0]) + "   " \
                                              + str(self.MEMC_DataInput[memc_i][4][1]) + " \t")

            ExchangeRatio_str = ''.join(ExchangeRatio_List)
            ExchangeLargeKind_str = ''.join(ExchangeLargeKind_List)
            LargeKindBackBone_str = ''.join(LargeKindBackBone_List)
            ExchangeSmallKind_str = ''.join(ExchangeSmallKind_List)
            SmallKindBackBone_str = ''.join(SmallKindBackBone_List)

            data_control_file.write('###############################\n')
            data_control_file.write('# MEMC PARAMETER \n')
            data_control_file.write('###############################\n')
            data_control_file.write('ExchangeVolumeDim \t ' + str(self.ExchangeVolumeDim[0]) + '\t' \
                                    + str(self.ExchangeVolumeDim[1]) + '\t' \
                                    + str(self.ExchangeVolumeDim[2]) + ' \n')
            data_control_file.write('ExchangeRatio \t\t ' + ExchangeRatio_str + '\n')
            data_control_file.write('ExchangeLargeKind \t ' + ExchangeLargeKind_str + '\n')
            data_control_file.write('ExchangeSmallKind \t ' + ExchangeSmallKind_str + '\n')
            if self.MEMC_DataInput != None and (self.MEMC_2Freq > 0
                                                           or self.IntraMEMC_2Freq > 0
                                                           or self.MEMC_3Freq > 0
                                                           or self.IntraMEMC_3Freq > 0 ):
                data_control_file.write('LargeKindBackBone \t ' + LargeKindBackBone_str + '\n')
            if self.MEMC_DataInput != None and (self.MEMC_2Freq > 0
                                                           or self.IntraMEMC_2Freq > 0 ):
                data_control_file.write('SmallKindBackBone \t ' + SmallKindBackBone_str + '\n')



        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# BOX DIMENSION #, X, Y, Z    (only orthoganal boxes are currently ' \
                                + 'available in this control file writer)\n')
        data_control_file.write('################################\n')
        data_control_file.write('CellBasisVector1 0 \t ' + str(self.x_dim_box_0) \
                                + '\t\t'+ '0.00' + '\t\t'+ '0.00' + '\n')
        data_control_file.write('CellBasisVector2 0 \t ' + '0.00' \
                                + '\t\t' + str(self.y_dim_box_0) + '\t\t' + '0.00' + '\n')
        data_control_file.write('CellBasisVector3 0 \t ' + '0.00' \
                                + '\t\t' + '0.00' + '\t\t' + str(self.z_dim_box_0) + '\n')
        data_control_file.write(' \n')
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('CellBasisVector1 1 \t ' + str(self.x_dim_box_1) \
                                    + '\t\t' + '0.00' + '\t\t' + '0.00' + '\n')
            data_control_file.write('CellBasisVector2 1 \t ' + '0.00' \
                                    + '\t\t' + str(self.y_dim_box_1) + '\t\t' + '0.00' + '\n')
            data_control_file.write('CellBasisVector3 1 \t ' + '0.00' \
                                    + '\t\t' + '0.00' + '\t\t' + str(self.z_dim_box_1) + '\n')
            data_control_file.write(' \n')

        if (self.ensemble_type in ['NPT', 'NVT']) \
                and self.FreeEnergyCalc != None and self.MoleculeType != None \
                and self.InitialState != None and self.LambdaVDW != None \
                and self.LambdaCoulomb != None :

            #make list for number of states, LambdaVDW, and Lambda_Coul, and convert the list to string for printing
            Lambda_states_List = []
            Lambda_VDW_List = []
            Lambda_Coul_List = []
            for lamda_i in range(0, len(self.LambdaVDW )):
                Lambda_states_List.append(str(lamda_i))
                Lambda_VDW_List.append(str(self.LambdaVDW[lamda_i]))
                Lambda_Coul_List.append(str(self.LambdaCoulomb[lamda_i]))

            Lambda_states_str = '\t'.join(Lambda_states_List)
            Lambda_VDW_str = '\t'.join(Lambda_VDW_List)
            Lambda_Coul_str = '\t'.join(Lambda_Coul_List)


            data_control_file.write('##############################\n')
            data_control_file.write('# FREE ENERGY PARAMETERS (only available in NPT and NVT ensembles) \n')
            data_control_file.write('##############################\n')
            data_control_file.write('FreeEnergyCalc \t ' + str(self.FreeEnergyCalc[0]) \
                                    + '\t\t' + str(self.FreeEnergyCalc[1]) + '\n')
            data_control_file.write('MoleculeType \t ' + str(self.MoleculeType[0]) \
                                    + '\t\t' + str(self.MoleculeType[1]) + '\n')
            data_control_file.write('InitialState \t ' + str(self.InitialState) + '\n')
            data_control_file.write('ScalePower \t ' + str(self.ScalePower) + '\n')
            data_control_file.write('ScaleAlpha \t ' + str(self.ScaleAlpha) + '\n')
            data_control_file.write('MinSigma \t ' + str(self.MinSigma) + '\n')
            data_control_file.write('ScaleCoulomb \t ' + str(self.ScaleCoulomb) + '\n')
            data_control_file.write('# States \t ' + str( Lambda_states_str) + '\n')
            data_control_file.write('LambdaVDW \t ' + str(Lambda_VDW_str) + '\n')
            data_control_file.write('LambdaCoulomb \t ' + str(Lambda_Coul_str) + '\n')
            data_control_file.write(' \n')

        if (self.ensemble_type in ['NPT', 'NVT']) \
                and self.FreeEnergyCalc != None and self.MoleculeType != None \
                and self.InitialState != None and self.LambdaVDW != None :

            #make list for number of states, LambdaVDW, and Lambda_Coul, and convert the list to string for printing
            Lambda_states_List = []
            Lambda_VDW_List = []
            Lambda_Coul_List = []
            for lamda_i in range(0, len(self.LambdaVDW )):
                Lambda_states_List.append(str(lamda_i))
                Lambda_VDW_List.append(str(self.LambdaVDW[lamda_i]))
                if self.LambdaCoulomb != None :
                    Lambda_Coul_List.append(str(self.LambdaCoulomb[lamda_i]))

            Lambda_states_str = '\t'.join(Lambda_states_List)
            Lambda_VDW_str = '\t'.join(Lambda_VDW_List)
            if self.LambdaCoulomb != None:
                Lambda_Coul_str = '\t'.join(Lambda_Coul_List)


            data_control_file.write('##############################\n')
            data_control_file.write('# FREE ENERGY PARAMETERS (only available in NPT and NVT ensembles) \n')
            data_control_file.write('##############################\n')
            data_control_file.write('FreeEnergyCalc \t ' + str(self.FreeEnergyCalc[0]) \
                                    + '\t\t' + str(self.FreeEnergyCalc[1]) + '\n')
            data_control_file.write('MoleculeType \t ' + str(self.MoleculeType[0]) \
                                    + '\t\t' + str(self.MoleculeType[1]) + '\n')
            data_control_file.write('InitialState \t ' + str(self.InitialState) + '\n')
            data_control_file.write('ScalePower \t ' + str(self.ScalePower) + '\n')
            data_control_file.write('ScaleAlpha \t ' + str(self.ScaleAlpha) + '\n')
            data_control_file.write('MinSigma \t ' + str(self.MinSigma) + '\n')
            data_control_file.write('ScaleCoulomb \t ' + str(self.ScaleCoulomb) + '\n')
            data_control_file.write('# States \t ' + str( Lambda_states_str) + '\n')
            data_control_file.write('LambdaVDW \t ' + str(Lambda_VDW_str) + '\n')
            if self.LambdaCoulomb != None:
                data_control_file.write('LambdaCoulomb \t ' + str(Lambda_Coul_str) + '\n')
            data_control_file.write(' \n')

        data_control_file.write('##############################\n')
        data_control_file.write('# CBMC TRIALS \n')
        data_control_file.write('##############################\n')
        data_control_file.write('CBMC_First \t ' + str(self.CBMC_First) + '\n')
        data_control_file.write('CBMC_Nth \t ' + str(self.CBMC_Nth) + '\n')
        data_control_file.write('CBMC_Ang \t ' + str(self.CBMC_Ang) + '\n')
        data_control_file.write('CBMC_Dih \t ' + str(self.CBMC_Dih) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('############################################################################\n')
        data_control_file.write('#  =======-------------------- OUTPUT --------------------------=========== \n')
        data_control_file.write('############################################################################\n')
        data_control_file.write(' \n')

        data_control_file.write('##########################\n')
        data_control_file.write('# statistics filename add\n')
        data_control_file.write('##########################\n')
        data_control_file.write('OutputName \t ' + str(self.OutputName) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('#####################################\n')
        data_control_file.write('# enable, frequency \n')
        data_control_file.write('#####################################\n')
        data_control_file.write('RestartFreq  \t\t ' + str(self.RestartFreq[0]) \
                                + '\t\t' + str(self.RestartFreq[1]) + '\n')
        data_control_file.write('CheckpointFreq  \t ' + str(self.CheckpointFreq[0]) \
                                + '\t\t' + str(self.CheckpointFreq[1]) + '\n')
        data_control_file.write('CoordinatesFreq  \t ' + str(self.CoordinatesFreq[0]) \
                                + '\t\t' + str(self.CoordinatesFreq[1]) + '\n')
        data_control_file.write('ConsoleFreq  \t\t ' + str(self.ConsoleFreq[0]) \
                                + '\t\t' + str(self.ConsoleFreq[1]) + '\n')
        data_control_file.write('BlockAverageFreq  \t ' + str(self.BlockAverageFreq[0]) \
                                + '\t\t' + str(self.BlockAverageFreq[1]) + '\n')
        data_control_file.write('HistogramFreq  \t\t ' + str(self.HistogramFreq[0]) \
                                + '\t\t' + str(self.HistogramFreq[1]) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# OutHistSettings \n')
        data_control_file.write('################################\n')
        data_control_file.write('DistName \t ' + str(self.DistName) + '\n')
        data_control_file.write('HistName \t ' + str(self.HistName) + '\n')
        data_control_file.write('RunNumber \t ' + str(self.RunNumber) + '\n')
        data_control_file.write('RunLetter \t ' + str(self.RunLetter) + '\n')
        data_control_file.write('SampleFreq \t ' + str(self.SampleFreq) + '\n')
        data_control_file.write(' \n')

        data_control_file.write('################################## \n')
        data_control_file.write('# enable: blk avg., fluct. \n')
        data_control_file.write('################################## \n')
        data_control_file.write('OutEnergy \t\t ' + str(self.OutEnergy[0]) \
                                + '\t\t' + str(self.OutEnergy[1]) + '\n')
        data_control_file.write('OutPressure \t\t ' + str(self.OutPressure[0]) \
                                + '\t\t' + str(self.OutPressure[1]) + '\n')
        data_control_file.write('OutMolNumber \t\t ' + str(self.OutMolNumber[0]) \
                                + '\t\t' + str(self.OutMolNumber[1]) + '\n')
        data_control_file.write('OutDensity \t\t ' + str(self.OutDensity[0]) \
                                + '\t\t' + str(self.OutDensity[1]) + '\n')
        data_control_file.write('OutVolume \t\t ' + str(self.OutVolume[0]) \
                                + '\t\t' + str(self.OutVolume[1]) + '\n')
        data_control_file.write('OutSurfaceTension \t ' + str(self.OutSurfaceTension[0]) \
                                + '\t\t' + str(self.OutSurfaceTension[1]) + '\n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')
        data_control_file.write(' \n')

        return "GOMC_CONTROL_FILE_WRITTEN"

    def ck_input_variable_True_or_False(self,
                                        input_variables_dict,
                                        key,
                                        bad_user_variable_List):
        """
        Checks if the input variable is either True for false.
        If not either True or False, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != True \
                and input_variables_dict[key] != False \
                or (str(input_variables_dict[key]) != str(True) \
                    and str(input_variables_dict[key]) != str(False)):
            bad_user_variable_List.append(key)

    def ck_input_variable_int_or_float_zero_or_greater(self,
                                                       input_variables_dict,
                                                       key,
                                                       bad_input_variables_values_List):
        """
        Checks if the input variable is an integer or float is zero or greater ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key.

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None and isinstance(
                input_variables_dict[key], int) != True \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_int_zero_or_greater(self,
                                              input_variables_dict,
                                              key,
                                              bad_input_variables_values_List):
        """
        Checks if the input variable is an integer is zero or greater   ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None and isinstance(
                input_variables_dict[key], int) != True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_float_zero_or_greater(self,
                                                input_variables_dict,
                                                key,
                                                bad_input_variables_values_List):
        """
        Checks if the input variable is a float is zero or greater  ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_int_or_float_greater_zero(self,
                                                    input_variables_dict,
                                                    key,
                                                    bad_input_variables_values_List):
        """
        Checks if the input variable is an integer or float is greater than zero.
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None and isinstance(
                input_variables_dict[key], int) != True \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_int_greater_zero(self,
                                           input_variables_dict,
                                           key,
                                           bad_input_variables_values_List):
        """
        Checks if the input variable is an integer greater than zero.
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None and isinstance(
                input_variables_dict[key], int) != True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_float_greater_zero(self,
                                             input_variables_dict,
                                             key,
                                             bad_input_variables_values_List):
        """
        Checks if the input variable is a float greater than zero.
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_float_greater_zero_less_1(self,
                                                    input_variables_dict,
                                                    key,
                                                    bad_input_variables_values_List):
        """
        Checks if the input variable is a float greater than zero and less than 1
        ( 0 < value < 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] <= 0 \
                or input_variables_dict[key] >= 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_int_or_float_zero_to_1(self,
                                                 input_variables_dict,
                                                 key,
                                                 bad_input_variables_values_List):
        """
        Checks if the input variable is an integer or float from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], int) != True \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_float_zero_to_1(self,
                                          input_variables_dict,
                                          key,
                                          bad_input_variables_values_List):
        """
        Checks if the input variable is a float from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], float) != True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_int_zero_to_1(self,
                                        input_variables_dict,
                                        key,
                                        bad_input_variables_values_List):
        """
        Checks if the input variable is an integer from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """
        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], int) != True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_List.append(key)

    def ck_input_variable_list_bool_int_zero_or_greater(self,
                                                        input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_List):
        """
        Checks if the input variable is a list with a bool and integer 0 or greater
        ([bool, int >= 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """
        if isinstance(input_variables_dict[key], list) == False:
            bad_input_variables_values_List.append(key)
        elif isinstance(input_variables_dict[key], list) == True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] != True \
                        and input_variables_dict[key][0] != False) \
                    or isinstance(input_variables_dict[key][1], int) != True \
                    or input_variables_dict[key][1] < 0 \
                    or (str(input_variables_dict[key][0]) != str(True) \
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or str(input_variables_dict[key][1]) == str(True) \
                    or str(input_variables_dict[key][1]) == str(False):
                bad_input_variables_values_List.append(key)

    def ck_input_variable_list_bool_int_greater_zero(self,
                                                     input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_List):
        """
        Checks if the input variable is a list with a bool and integer greater than zero
        ([bool, int > 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """
        if isinstance(input_variables_dict[key], list) == False:
            bad_input_variables_values_List.append(key)
        elif isinstance(input_variables_dict[key], list) == True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] != True \
                        and input_variables_dict[key][0] != False) \
                    or isinstance(input_variables_dict[key][1], int) != True \
                    or input_variables_dict[key][1] <= 0 \
                    or (str(input_variables_dict[key][0]) != str(True) \
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or str(input_variables_dict[key][1]) == str(True) \
                    or str(input_variables_dict[key][1]) == str(False):
                bad_input_variables_values_List.append(key)

    def ck_input_variable_list_residue_str_int_greater_zero(self,
                                                            input_variables_dict,
                                                            key,
                                                            bad_input_variables_values_List):
        """
        Checks if the input variable is a list with a bool and integer greater than zero
        ([bool, int > 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if isinstance(input_variables_dict[key], list) == False:
            bad_input_variables_values_List.append(key)
        elif isinstance(input_variables_dict[key], list) == True:
            if len(input_variables_dict[key]) == 2:
                if isinstance(input_variables_dict[key], list) == True:
                    if isinstance(input_variables_dict[key][0], str) == False \
                            or isinstance(input_variables_dict[key][1], int) == False:
                        bad_input_variables_values_List.append(key)
                    elif isinstance(input_variables_dict[key][0], str) == True:
                        if len(input_variables_dict[key][0]) > 4 \
                                or input_variables_dict[key][0] not in self.residues_List:
                            bad_input_variables_values_List.append(key)

                    if isinstance(input_variables_dict[key][1], int) == True \
                            and input_variables_dict[key][1] <= 0:
                        bad_input_variables_values_List.append(key)

    def ck_input_variable_list_bool_bool(self,
                                         input_variables_dict,
                                         key,
                                         bad_input_variables_values_List):
        """
        Checks if the input variable is a list with a 2 booleans  ([bool, bool]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if isinstance(input_variables_dict[key], list) == False:
            bad_input_variables_values_List.append(key)
        elif isinstance(input_variables_dict[key], list) == True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] != True \
                        and input_variables_dict[key][0] != False) \
                    or (input_variables_dict[key][1] != True \
                        and input_variables_dict[key][1] != False) \
                    or (str(input_variables_dict[key][0]) != str(True) \
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or (str(input_variables_dict[key][1]) != str(True) \
                        and str(input_variables_dict[key][1]) != str(False)):
                bad_input_variables_values_List.append(key)

    def ck_input_variable_list_of_floats_zero_to_1(self,
                                                   input_variables_dict,
                                                   key,
                                                   bad_input_variables_values_List):
        """
        Checks if the input variable is a list of floats between 0 and 1 ([0,0, 0.1, ..., 1.0]).
        Note: the list can be of any length with 0.0 <= float <= 1.0
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None:
            if isinstance(input_variables_dict[key], list) == False:
                bad_input_variables_values_List.append(key)
            elif isinstance(input_variables_dict[key], list) == True:
                if len(input_variables_dict[key]) >= 1:
                    for Lambda_i in range(0, len(input_variables_dict[key])):
                        if isinstance(input_variables_dict[key][Lambda_i], float) == False:
                            bad_input_variables_values_List.append(key)
                        elif isinstance(input_variables_dict[key][Lambda_i], float) == True:
                            if input_variables_dict[key][Lambda_i] < 0.0 \
                                    or input_variables_dict[key][Lambda_i] > 1.0:
                                bad_input_variables_values_List.append(key)
                else:
                    bad_input_variables_values_List.append(key)

    def ck_input_variable_GCMC_dict_str_int_or_float(self,
                                                     input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_List):
        """
        Checks if the input variable is a dictionary with a key = string and
        value = integer or float  ({'str_1' : integer_1 or float_1, ....,
        'str_x' : integer_x or float_x }).
        Note: the dictionary can be of any length
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], dict) != True:
            bad_input_variables_values_List.append(key)

        elif isinstance(input_variables_dict[key], dict) == True:
            keys_List = dict_keys_to_list(input_variables_dict[key])
            for keys_iter_No in range(0, len(keys_List)):
                key_iter = keys_List[keys_iter_No]
                value_iter = input_variables_dict[key][key_iter]

                if key_iter not in self.residues_List:
                    bad_input_variables_values_List.append(key)

                if isinstance(key_iter, str) == False:
                    bad_input_variables_values_List.append(key)
                elif isinstance(key_iter, str) == True:
                    if len(key_iter) > 4:
                        bad_input_variables_values_List.append(key)

                if (isinstance(value_iter, int) == False \
                    and isinstance(value_iter, float) == False) \
                        or str(value_iter) == str(True) \
                        or str(value_iter) == str(False):
                    bad_input_variables_values_List.append(key)

    def ck_input_variable_GCMC_dict_str_int_or_float_zero_or_greater(self,
                                                                     input_variables_dict,
                                                                     key,
                                                                     bad_input_variables_values_List):
        """
        Checks if the input variable is a dictionary with a key = string and
        value = integer or float zero or greater  ({'str_1' : integer_1 or float_1 (>= 0), ....,
        'str_x' : integer_x or float_x (>= 0)}).
        Note: the dictionary can be of any length
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] != None \
                and isinstance(input_variables_dict[key], dict) != True:
            bad_input_variables_values_List.append(key)

        elif isinstance(input_variables_dict[key], dict) == True:
            keys_List = dict_keys_to_list(input_variables_dict[key])
            for keys_iter_No in range(0, len(keys_List)):
                key_iter = keys_List[keys_iter_No]
                value_iter = input_variables_dict[key][key_iter]

                if key_iter not in self.residues_List:
                    bad_input_variables_values_List.append(key)

                if isinstance(key_iter, str) == False:
                    bad_input_variables_values_List.append(key)
                elif isinstance(key_iter, str) == True:
                    if len(key_iter) > 4:
                        bad_input_variables_values_List.append(key)

                if (isinstance(value_iter, int) == False \
                    and isinstance(value_iter, float) == False) \
                        or str(value_iter) == str(True) \
                        or str(value_iter) == str(False) \
                        or value_iter < 0:
                    bad_input_variables_values_List.append(key)

    def ck_input_variable_str_with_no_spaces(self,
                                             input_variables_dict,
                                             key,
                                             bad_input_variables_values_List):
        """
        Checks if the input variable is a string with no spaces.
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict; dict; The user input variable dictionary
        key = str; dictionary key for the user provided input variable list
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Outputs
        ---------
        bad_user_variable_List =  list; A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """
        if isinstance(input_variables_dict[key], str) == True:
            No_spaces_in_OutputName_string = " " in input_variables_dict[key]
        if isinstance(input_variables_dict[key], str) != True \
                or No_spaces_in_OutputName_string == True:
            bad_input_variables_values_List.append(key)



# user callable function
def write_gomc_control_file(charmm_object, conf_filename,  ensemble_type,
                            RunSteps, Temperature, input_variables_dict = None):
    """
    Construct the GOMC configuration input file with the defaults,
    or adding additional data in the input_variable section
    Default setting for the GOMC configuraion files are based upon the
    a educated guess which should result in reasonable sampling for a
    given ensemble/simulation type. However, there is no guarantee that
    the default setting will provide the best or adequate sampling for
    the selected system. The user has the option to modify the
    configuration/contorl files based on the simulation specifics or in to
    optimize the system beyond the standard settings.  These override
    options are available via the keyword arguments in input_variable_dict.
    Parameters
    ----------
    charmm_object : CHARMM object (i.e., Charmm()); the Charmm object from the charmm_writer file
    ensemble_type : str, only accepts 'NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', 'GCMC'
    RunSteps : int; must be an integer greater than zero.
        Sets the total number of simulation steps.
    Temperature : Temperature of system in Kelvin (K)
    input_variables_dict: dict, default = None
        These input variables are optional and override the default settings.
        Changing these variables likely required for more advanced systems.

        The details of the acceptable input variables for the selected
        ensembles can be found by running this python workbook,

            print_valid_ensemble_input_variables('GCMC', description = True)

        which prints the input_variables with their subsection description
        for the selected 'GCMC' ensemble (other ensembles can be set as well).

        Example : input_variables_dict = {'Restart' : False, 'PRNG' : 123,
                                          'ParaTypeCHARMM' : True }

    Returns
    -------
    If completed without errors; str; "PASSED
        and the GOMC input control file is writen
    If completed with errors;  None

    Notes
    -------
    The details of the required inputs for the selected
    ensembles can be found by running this python workbook,

         print_valid_required_input_variables('NVT', description = True)

    which prints the required inputs with their subsection description
    for the selected 'NVT' ensemble (other ensembles can be set as well).
    """

    gomc_control = GOMCControl(charmm_object, ensemble_type,
                               RunSteps, Temperature, input_variables_dict = input_variables_dict )
    test_gomc_control_write_conf_file = gomc_control.write_conf_file(conf_filename)

    if gomc_control.all_inputs_pass == True  and test_gomc_control_write_conf_file == "GOMC_CONTROL_FILE_WRITTEN":

        return "GOMC_CONTROL_FILE_WRITTEN"
    else:
        return None







"""
#testing the functions
print('************************')
print('Printing valid ensemble input variables (start)')
print('************************')
print_valid_ensemble_input_variables('GEMC_NPT', description = True)
print('************************')
print('Printing valid ensemble input variables (end)')
print('************************')
print('************************')
print('Printing required ensemble file (start)')
print('************************')
print_required_ensemble_files('GEMC_NPT', description = True)
print('************************')
print('Printing required ensemble file (end)')
print('************************')

print('************************')
print('Writing a control file for a NVT or NPT example (start)')
print('************************')
test_0_box = GOMCControl( 'NVT', 10, 300, 'FF.inp', 0, \
                              'system_box_0.pdb', 'system_box_0.psf',\
                              45.0, 55, 65.0, \
                              Coordinates_box_1 = None, Structures_box_1 = None, \
                              x_dim_box_1 = None, y_dim_box_1 = None, z_dim_box_1 = None, \
                              input_variables_dict= {'Restart' : False, 'PRNG' : 123,
                                              'ParaTypeCHARMM' : True,
                                              'ParaTypeMARTINI' : False,
                                              'ParaTypeMie' : False,
                                              "RcutCoulomb_box_0" : 2,
                                              "PressureCalc" : [False, 4],
                                              "Tolerance" : .1,
                                              "DisFreq": 0.2,
                                              "RotFreq": 0.2,
                                              "IntraSwapFreq": 0.1,
                                              "RegrowthFreq": 0.1,
                                              "CrankShaftFreq": 0.2,
                                              "MultiParticleFreq": 0.2,
                                                     "CBMC_First": 2,
                                                     "CBMC_Nth": 343,
                                                     "CBMC_Dih": 343,
                                                     "OutputName": 'test_out',
                                                     "RestartFreq": [True, 100],
                                                     "CheckpointFreq": [True, 100],
                                                     "CoordinatesFreq": [True, 100],
                                                     "ConsoleFreq": [True, 100],
                                                     "BlockAverageFreq": [True, 100],
                                                     "HistogramFreq": [True, 100],
                                                     "DistName": 'dis',
                                                     "HistName": 'his',
                                                     "RunNumber": 1,
                                                     "RunLetter": 'a',
                                                     "SampleFreq": 500,
                                                     "FreeEnergyCalc": [True, 10000],
                                                     "MoleculeType": ['WAT', 1],
                                                     "InitialState": 3,
                                                     "LambdaVDW": [0.1, 0.2, 0.4, 0.9],
                                                     "LambdaCoulomb": [0.1, 0.3, 0.8, 0.8],

        } )



test_0_box.write_conf_file('Test_box_0')

print('************************')
print('Writing a control file for a NVT or NPT example (end)')
print('************************')


print('************************')
print('Writing a control file for a GEMC_NVT, GEMC_NPT, or GCMC example (start)')
print('************************')

test_0_and_1_box = GOMCControl( 'GEMC_NPT', 10, 300, '../Generate_liq_vap_boxes/GOMC_water_fake_water_FF.inp', 1 ,\
                                    '../equilb_liq_box_final_files/NPT_350K/restart_60M_to_70M_steps/Output_data_BOX_0_restart.pdb', \
                                    '../equilb_liq_box_final_files/NPT_350K/restart_60M_to_70M_steps/Output_data_merged.psf',\
                                    29.911, 29.911, 29.911, \
                                    Coordinates_box_1 = '../Generate_liq_vap_boxes/Box_1_water_box_vap_300A_L_cubed_box_350K.pdb', \
                                    Structures_box_1 = '../Generate_liq_vap_boxes/Box_1_water_box_vap_300A_L_cubed_box_350K.psf',\
                                    x_dim_box_1 = 300, y_dim_box_1 = 300, z_dim_box_1 = 300, \
                                    input_variables_dict= {'Restart' : False, 'PRNG' : 123,
                                                           'OutputName' : 'test_name',
                                              'ParaTypeCHARMM' : True,
                                              'ParaTypeMARTINI' : False,
                                              'ParaTypeMie' : False,
                                              'Potential' : 'VDW',
                                              'LRC' : False,
                                              'Rcut' : 10,
                                              'RcutLow' : 1,
                                              'Rswitch' : 9,
                                              "RcutCoulomb_box_0" : 0,
                                              "RcutCoulomb_box_1" : 9,
                                              "PressureCalc" : [True, 1000],
                                                           "CBMC_First" : 1 ,
                                              "CBMC_Nth" : 43,
                                              "CBMC_Dih" : 4,
                                              #"ChemPot" : {'WAT': 0},

                                                           "DistName": 'dis',
                                                           "HistName": 'his',
                                                           "RunNumber": 1,
                                                           "RunLetter": 'a',
                                                           "SampleFreq": 500,
                                                           "Dielectric": 10,

                                                           "OutEnergy": [True, True],
                                                           "OutPressure": [True, True],
                                                           "OutMolNumber": [True, True],
                                                           "OutDensity": [True, True],
                                                           "OutVolume": [True, True],
                                                           "OutSurfaceTension": [False, False],
                                                           "DisFreq": 0.09,
                                                           "RotFreq": 0.1,
                                                           "IntraSwapFreq": 0.1,
                                                           "SwapFreq": 0.2,
                                                           "RegrowthFreq": 0.1,
                                                           "CrankShaftFreq": 0.1,
                                                           "MultiParticleFreq": 0.0,
                                                           "MEMC-1Freq" : 0.0,
                                                           "MEMC-2Freq" : 0.00,
                                                           "MEMC-3Freq" : 0.0,
                                                           "IntraMEMC-1Freq": 0.00,
                                                           "IntraMEMC-2Freq" : 0.0,
                                                           "IntraMEMC-3Freq": 0.3,
                                                           "ExchangeVolumeDim": [1, 1, 1],
                                                           "MEMC_DataInput": [
                                                               [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1']],
                                                               [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1']]
                                                           ],


                                              } )

test_0_and_1_box.write_conf_file('Test_box_0_and_1')
print('************************')
print('Writing a control file for a GEMC_NVT, GEMC_NPT, or GCMC example (end)')
print('************************')

"[ExchangeRatio_int, ExchangeLargeKind_str, ExchangeSmallKind_str, [LargeKindBackBone_atom_1_str, LargeKindBackBone_atom_2_str ], [SmallKindBackBone_atom_1_str, SmallKindBackBone_atom_2_str ]"
"""
