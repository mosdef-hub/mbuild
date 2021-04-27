import datetime
import os
import mbuild.formats.charmm_writer as mf_charmm
from warnings import warn


def dict_keys_to_list(dict):
    """
    Converts the dictionary keys into a list

    Parameters
    ----------
    dict : dict
        A provided dictionary

    Returns
    ---------
    list : list
        list of keys from the provided dictionary
    """
    return [key for key in dict.keys()]

def print_valid_required_input_variables(description=False):
    """
    Prints the valid required input, which is necessary to write the GOMC control file.

    Parameters
    ----------
    description : bool, default = False
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Returns
    ---------
    Prints out the valid input variables (user optional) on the screen,
        which can be entered in the GOMC writer. These are the valid input
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
    description :  bool, default = False.
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Returns
    ---------
    required_data : dict or list, default = list.
        If the description = True then a dict is provided with the key and value.
        if the description = False then a list of the dict keys is provided.

    Note:
    Variables and text extracted with permission from the GOMC manual version 2.60.
    Some of the text was modified from its original version.
    Cite: Potoff, Jeffrey; Schwiebert, Loren; et al. GOMC Documentation.
    https://raw.githubusercontent.com/GOMC-WSU/GOMC/master/GOMC_Manual.pdf, 2021.
    """

    required_data = {
        "charmm_object": 'Charmm object, ' 
                         'A Charmm object, which by definition has been parameterized ' 
                         'from the selected force field.',
        "ensemble_type": "Required files or System Info (all ensembles): str, " 
                         "(valid strings are 'NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', or 'GCMC'), " 
                         'the ensemble type for the simulation.',
        "RunSteps": "Required files or System Info (all ensembles): int (> 0), " 
                     "The number or run steps for the simulation.",
        "Temperature":  "Required files or System Info (all ensembles): float or integer (> 0), " 
                        "Temperature of system in Kelvin (K)",
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
    description :  bool, default = False.
        If True, it prints the descriptions of the input_variables (i.e. dict),
        If False, it only prints the input_variables without the descriptions (i.e. list)

    Returns
    ---------
    valid_input_variables : dict or list, default = list.
        If the description = True then a dict is provided with the key and value.
        if the description = False then a list of the dict keys is provided.

    Note:
    Variables and text extracted with permission from the GOMC manual version 2.60.
    Some of the text was modified from its original version.
    Cite: Potoff, Jeffrey; Schwiebert, Loren; et. al. GOMC Documentation.
    https://raw.githubusercontent.com/GOMC-WSU/GOMC/master/GOMC_Manual.pdf, 2021.
    """

    valid_input_variables = {

        # ******************************************************************************************************
        # Definitions in this function are copied to a large extent from the GOMC manual release version 2.60 (start)
        # insert citation here:
        # ******************************************************************************************************
        "Restart": 'Simulation info (all ensembles): boolean, default = {}. ' 
                   'Determines whether to restart the simulation ' 
                   'from restart file (*_restart.pdb and *_restart.psf) or not.'
                   ''.format(_get_default_variables_dict()["Restart"]),
        "RestartCheckpoint": 'Simulation info (all ensembles): boolean, default = {}. ' 
                             'Determines whether to restart the ' 
                             'simulation with the checkpoint file (checkpoint.dat) or not. Restarting the ' 
                             'simulation with checkpoint.dat would result in an identical outcome, as if ' 
                             'previous simulation was continued.'
                             ''.format(_get_default_variables_dict()["RestartCheckpoint"]),
        "PRNG": 'Simulation info (all ensembles): string or int (>= 0) ("RANDOM" or integer), default = {}. ' 
                'PRNG = Pseudo-Random Number Generator (PRNG). ' 
                'There are two (2) options, entering the string, "RANDOM", or a integer.  \n' 
                '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "RANDOM", which selects a random seed number. ' 
                'This will enter the line "PRNG RANDOM" in the gomc configuration file. \n'
                '\t\t\t\t\t\t\t\t\t\t\t\t\t --- integer, which defines the integer seed number ' 
                'for the simulation. ' 
                'This is equivalent to entering the following two lines in the configuration file: ' 
                'line 1 = PRNG INTSEED, ' 
                'line 2 = Random_Seed user_selected_integer. '
                'Example 1: for a random seed enter the string "RANDOM. ' 
                'Example 2: for a specific seed number enter a integer of your choosing. '
                ''.format(_get_default_variables_dict()["PRNG"]),
        "ParaTypeCHARMM": 'Simulation info (all ensembles): boolean, default = {}. ' 
                          'True if a CHARMM forcefield, False otherwise.'
                          ''.format(_get_default_variables_dict()["ParaTypeCHARMM"]),
        "ParaTypeMie": 'Simulation info (all ensembles): boolean, default = {}. ' 
                       'True if a Mie forcefield type, False otherwise.'
                       ''.format(_get_default_variables_dict()["ParaTypeMie"]),
        "ParaTypeMARTINI": 'Simulation info (all ensembles): boolean, default = {}. ' 
                           'True if a MARTINI forcefield, False otherwise.'
                           ''.format(_get_default_variables_dict()["ParaTypeMARTINI"]),
        "RcutCoulomb_box_0": 'Simulation info (all ensembles): int or float (>= 0), default = {}.'
                             'Sets a specific radius in box 0 where the short-range ' 
                             'electrostatic energy will be calculated (i.e., The distance to truncate the ' 
                             'short-range electrostatic energy in box 0.)'
                             'Note: if None, GOMC will default to the Rcut value'
                             ''.format(_get_default_variables_dict()["RcutCoulomb_box_0"]),
        "RcutCoulomb_box_1": 'Simulation info (all ensembles): int or float (>= 0), default = {}.'
                             'Sets a specific radius in box 1 where the short-range  '
                             'electrostatic energy will be calculated. (i.e., The distance to truncate the ' 
                             'short-range electrostatic energy in box 1.)'
                             'Note: if None, GOMC will default to the Rcut value'
                             ''.format(_get_default_variables_dict()["RcutCoulomb_box_1"]),
        "Pressure": 'Simulation info (only GEMC_NPT and NPT): int or float (>= 0), default = {}. ' 
                    'The pressure in bar utilized for the NPT ' 
                    'and GEMC_NPT simulations.'
                    ''.format(_get_default_variables_dict()["Pressure"]),
        "Rcut": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut), default = {}. '
                'Sets a specific radius in Angstroms that non-bonded interaction ' 
                'energy and force will be considered and calculated using defined potential function. ' 
                'The distance in Angstoms to truncate the LJ, Mie, or other VDW type potential at. '
                'Note: Rswitch is only used when the "Potential" = SWITCH. '
                ''.format(_get_default_variables_dict()["Rcut"]),
        "RcutLow": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut), default = {}. '
                   'Sets a specific minimum possible distance in Angstroms that reject ' 
                   'any move that places any atom closer than specified distance. The minimum possible '
                   'distance between any atoms. '  
                   'Sets a specific radius in Angstroms that non-bonded interaction '
                   'Note: Rswitch is only used when the "Potential" = SWITCH. '
                   ''.format(_get_default_variables_dict()["RcutLow"]),
        "LRC": 'Simulation info (all ensembles): boolean, default = {}. ' 
               'If True, the simulation considers the long range tail corrections for the non-bonded VDW or '
               'dispersion interactions. '
               'Note: In case of using SHIFT or SWITCH potential functions, LRC will be ignored.'
               ''.format(_get_default_variables_dict()["LRC"]),
        "Exclude": 'Simulation info (all ensembles): str ' 
                   '(The string inputs are "1-2", "1-3", or "1-4"), default = {}. ' 
                   'Note: In CHARMM force field, the 1-4 interaction needs to be considered. ' 
                   'Choosing "Excude 1-3", will modify 1-4 interaction based on 1-4 parameters ' 
                   'in parameter file. If a kind force field is used, where ' 
                   '1-4 interaction needs to be ignored, such as TraPPE, either Exlcude "1-4" needs to be ' 
                   'chosen or 1-4 parameter needs to be assigned to zero in the parameter file. \n'
                   '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "1-2": All interaction pairs of bonded atoms, ' 
                   'except the ones that separated with one bond, ' 
                   'will be considered and modified using 1-4 parameters defined in parameter file. \n'
                   '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "1-3": All interaction pairs of bonded atoms, ' 
                   'except the ones that separated with one or two ' 
                   'bonds, will be considered and modified using 1-4 parameters defined in parameter file. \n' 
                   '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "1-4": All interaction pairs of bonded atoms, ' 
                   'except the ones that separated with one, ' 
                   'two or three bonds, will be considered using non-bonded parameters defined in parameter file.'
                   ''.format(_get_default_variables_dict()["Exclude"]),
        "Potential": 'Simulation info (all ensembles): str, ["VDW", "EXP6", "SHIFT" or "SWITCH"], default = {}. ' 
                     'Defines the potential function type to calculate non-bonded dispersion interaction ' 
                     'energy and force between atoms. \n' 
                     '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "VDW":    Non-bonded dispersion interaction energy and force ' 
                     'calculated based on n-6 (Lennard - Jones) equation. This function will be discussed ' 
                     'further in the Intermolecular energy and ' 
                     'Virial calculation section. \n'  
                     '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "EXP6":   Non-bonded dispersion interaction energy and force ' 
                     'calculated based on exp-6 (Buckingham potential) equation. \n' 
                     '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "SHIFT":  This option forces the potential energy to be ' 
                     'zero at Rcut distance.  \n' 
                     '\t\t\t\t\t\t\t\t\t\t\t\t\t --- "SWITCH": This option smoothly forces the potential ' 
                     'energy to be zero at Rcut distance and starts modifying the potential at Rswitch ' 
                     'distance. Depending on force field type, specific potential function will be applied. '
                     ''.format(_get_default_variables_dict()["Potential"]),
        "Rswitch": 'Simulation info (all ensembles): int or float (>= 0 and RcutLow < Rswitch < Rcut), default = {}. '
                   'Note: Rswitch is only used when the SWITCH function is used (i.e., "Potential" = SWITCH). '
                   'The Rswitch distance is in Angstrom. If the “SWITCH” function is chosen, ' 
                   'Rswitch needs to be defined, otherwise, the program will be terminated. When using ' 
                   'choosing "SWITCH" as potential function, the Rswitch distance defines where the' 
                   'non-bonded interaction energy modification is started, which is eventually truncated ' 
                   'smoothly at Rcut distance.'
                   ''.format(_get_default_variables_dict()["Rswitch"]),
        "ElectroStatic": 'Simulation info (all ensembles): boolean, default = {}. '
                         'Considers the coulomb interactions or not. '
                         'If True, coulomb interactions are considered and false if not. '
                         'Note: To simulate the polar molecule in MARTINI force field, ElectroStatic needs to be ' 
                         'turn on. The MARTINI force field uses short-range coulomb interaction with constant '
                         'Dielectric of 15.0.' 
                         ''.format(_get_default_variables_dict()["ElectroStatic"]),
        "Ewald":  'Simulation info (all ensembles): boolean, default = {}. '
                  'Considers the standard Ewald summation method for electrostatic calculations. ' 
                  'If True, Ewald summation calculation needs to be considered and false if not. '
                  'Note: By default, GOMC will set ElectroStatic to True if Ewald summation  ' 
                  'method was used to calculate coulomb interaction.'
                  ''.format(_get_default_variables_dict()["Ewald"]),
        "CachedFourier": 'Simulation info (all ensembles): boolean, default = {}. ' 
                         'Considers storing the reciprocal terms for Ewald summation ' 
                         'calculation in order to improve the code performance. This option would increase the code ' 
                         'performance with the cost of memory usage. If True, to store reciprocal terms of Ewald ' 
                         'summation calculation and False if not. ' 
                         'Warning: Monte Carlo moves, such as MEMC-1, MEMC-2, MEMC-3, ' 
                         'IntraMEMC-1, IntraMEMC-2, and IntraMEMC-3 are not support with CachedFourier.'
                         ''.format(_get_default_variables_dict()["CachedFourier"]),
        "Tolerance": 'Simulation info (all ensembles): float (0.0 < float < 1.0), default = {}. ' 
                     'Sets the accuracy in Ewald summation calculation. Ewald separation parameter and number ' 
                     'of reciprocal vectors for the Ewald summation are determined based on the accuracy parameter.'
                     ''.format(_get_default_variables_dict()["Tolerance"]),
        "Dielectric": 'Simulation info (all ensembles): int or float (>= 0.0), default = {}. ' 
                      'Sets dielectric value used in coulomb interaction when the Martini ' 
                      'force field is used. Note: In MARTINI force field, Dielectric needs to be set to 15.0.'
                      ''.format(_get_default_variables_dict()["Dielectric"]),
        "PressureCalc": 'Simulation info (all ensembles): list [bool , int (> 0)] or [bool , step_frequency], ' 
                        'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                        'Calculate the system pressure or not. bool = True, enables the pressure calculation ' 
                        'during the simulation, false disables the calculation. The int/step frequency sets the ' 
                        'frequency of calculating the pressure.'
                        ''.format(_get_default_variables_dict()["PressureCalc"],
                                  _get_default_variables_dict()["PressureCalc"][0],
                                  _get_default_variables_dict()["PressureCalc"][1]),
        "EqSteps": 'Simulation info (all ensembles): int (> 0), '
                   'default = set via formula based on the number of RunSteps or {} max. ' 
                   'Sets the number of steps necessary to equilibrate the system. ' 
                   'Averaging will begin at this step. ' 
                   'Note: In GCMC simulations, the Histogram files will be outputed at EqSteps.'
                   ''.format(_get_default_variables_dict()["EqSteps"]),
        "AdjSteps": 'Simulation info (all ensembles): int (> 0), '
                    'default = set via formula based on the number of RunSteps or {} max. ' 
                    'Sets the number of steps per adjustment of the parameter associated with each move ' 
                    '(e.g. maximum translate distance, maximum rotation, maximum volume exchange, etc.).'
                    ''.format(_get_default_variables_dict()["AdjSteps"]),
        "VDWGeometricSigma": 'Simulation info (all ensembles): boolean, default = {}. '
                             'Use geometric mean, as required by OPLS force field, '
                             'to combining Lennard-Jones sigma parameters for different atom types. '
                             'If set to True, GOMC uses geometric mean to combine Lennard-Jones or VDW sigmas. '
                             'Note: The default setting of VDWGeometricSigma is false, which uses the arithmetic '
                             'mean when combining Lennard-Jones or VDW sigma parameters for different atom types.'
                             ''.format(_get_default_variables_dict()["VDWGeometricSigma"]),
        "useConstantArea": 'Simulation info (only GEMC_NPT and NPT): boolean: default = {}. ' 
                           'Changes the volume of the simulation box by fixing the cross-sectional ' 
                           'area (x-y plane). If true, the volume will change only in z axis, If False, '
                           'the volume of the box will change in a way to maintain the constant axis ratio. '
                           ''.format(_get_default_variables_dict()["useConstantArea"]),
        "FixVolBox0": 'Simulation info (only GEMC_NPT): boolean, default = {}. ' 
                      'Changing the volume of fluid phase (Box 1) to maintain the constant imposed pressure and ' 
                      'Temperature, while keeping the volume of adsorbed phase (Box 0) fixed. Note: By default, ' 
                      'GOMC will set useConstantArea to False if no value was set. It means, the volume of the ' 
                      'box will change in a way to maintain the constant axis ratio.'
                      ''.format(_get_default_variables_dict()["FixVolBox0"]),
        # GCMC only properties
        "ChemPot": 'Simulation info (only GCMC): dict {str (4 dig limit) , int or float}, ' +
                   'default = {} (i.e., the user must set this variable as there is no working default).'''
                   ''.format(_get_default_variables_dict()["ChemPot"]) +
                   'The chemical potentials in GOMC units of energy, K. ' 
                   'There is a 4 character limit for the string/residue name since the PDB/PSF '
                   'files have a 4 character limitation and require and exact match in the conf file. '
                   'Note: These strings must match the residue in the psf and psb files or it will fail. ' 
                   'The name of the residues and their corresponding chemical potential must specified ' 
                   'for every residue in the system (i.e., {"residue_name" : chemical_potential}). ' 
                   'Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY ' 
                   'OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER. '
                   'Example 1 (system with only water):  {"H2O" : -4000} . ' 
                   'Example 2 (system with water and ethanol):  {"H2O" : -4000, "ETH" : -8000} ',
        "Fugacity": 'Simulation info (only GCMC): dict {str , int or float (>= 0)}, ' +
                    'default = {} (i.e., the user must set this variable as there is no working default). ' 
                    ''.format(_get_default_variables_dict()["Fugacity"]) +
                    'The fugacity in GOMC units of pressure, bar. '
                    'There is a 4 character limit for the string/residue name since the PDB/PSF ' 
                    'files have a 4 character limitation and require and exact match in the conf file. ' 
                    'Note: These strings must match the residue in the psf and psb files or it will fail. ' 
                    'The name of the residues and their corresponding fugacity must specified ' 
                    'for every residue in the system (i.e., {"residue_name" : fugacity}). ' 
                    'Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY '
                    'OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER. ' 
                    'Example 1 (system with only water):  {"H2O" : 1} . ' 
                    'Example 2 (system with water and ethanol):  {"H2O" : 0.5, "ETH" : 10} ',

        # CBMC inputs
        "CBMC_First": 'CBMC inputs (all ensembles): int (>= 0), default = {}, ' 
                      'The number of CD-CBMC trials to choose the first atom position' 
                      '(Lennard-Jones trials for first seed growth).'
                      ''.format(_get_default_variables_dict()["CBMC_First"]),
        "CBMC_Nth": 'CBMC inputs (all ensembles): int (>= 0), default = {},  ' 
                    'The number of CD-CBMC trials to choose the later atom positions ' 
                    '(Lennard-Jones trials for first seed growth).'
                    ''.format(_get_default_variables_dict()["CBMC_Nth"]),
        "CBMC_Ang": 'CBMC inputs (all ensembles): int (>= 0), default = {}, ' 
                    'The number of CD-CBMC bending angle trials to perform for geometry ' 
                    '(per the coupled-decoupled CBMC scheme).'
                    ''.format(_get_default_variables_dict()["CBMC_Ang"]),
        "CBMC_Dih": 'CBMC inputs (all ensembles): int (>= 0), default = {}, '
                    'The number of CD-CBMC dihedral angle trials to perform for geometry '
                    '(per the coupled-decoupled CBMC scheme).'
                    ''.format(_get_default_variables_dict()["CBMC_Dih"]),

        # Control file (.conf file ) output controls/parameters
        "OutputName": 'Output Frequency (all ensembles): str (NO SPACES), default = {}. '
                      'The UNIQUE STRING NAME, WITH NO SPACES, which is used for the '
                      'output block average, PDB, and PSF file names.'
                      ''.format(_get_default_variables_dict()["OutputName"]),
        "CoordinatesFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' 
                           '[Generate_data_bool , steps_per_data_output_int], ' 
                           'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                           'Controls output of PDB file (coordinates). ' 
                           'If bool is True, this enables outputting the coordinate files at the ' 
                           'integer frequency (set steps_per_data_output_int), '
                           'while "False" disables outputting the coordinates.'
                           ''.format(_get_default_variables_dict()["CoordinatesFreq"],
                                     _get_default_variables_dict()["CoordinatesFreq"][0],
                                     _get_default_variables_dict()["CoordinatesFreq"][1]),
        "RestartFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' 
                       '[Generate_data_bool , steps_per_data_output_int], ' 
                       'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                       'This creates the PDB and PSF (coordinate and topology) files for restarting the system '
                       'at the set steps_per_data_output_int (frequency) '
                       'If bool is True, this enables outputting the PDB/PSF restart files at the ' 
                       'integer frequency (set steps_per_data_output_int), ' 
                       'while “false” disables outputting the PDB/PSF restart files. '
                       ''.format(_get_default_variables_dict()["RestartFreq"],
                                 _get_default_variables_dict()["RestartFreq"][0],
                                 _get_default_variables_dict()["RestartFreq"][1]),
        "CheckpointFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' 
                          '[Generate_data_bool , steps_per_data_output_int], '
                          'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                          'Controls the output of the last state of simulation at a specified step, in a '
                          'binary file format (checkpoint.dat). Checkpoint file contains the following '
                          'information in full precision: '
                          '(1) Last simulation step that saved into checkpoint file '
                          '(2) Simulation cell dimensions and angles ' 
                          '(3) Maximum amount of displacement (Å), rotation (δ), and volume (Å^3) that is used in ' 
                          'the Displacement, Rotation, MultiParticle, and Volume moves ' 
                          '(4) Number of Monte Carlo move trial and acceptance ' 
                          '(5) All molecule’s coordinates ' 
                          '(6) Random number sequence '
                          'If bool is True, this enables outputting the checkpoint file at the ' 
                          'integer frequency (set steps_per_data_output_int), ' 
                          'while "False" disables outputting the checkpoint file.'
                          ''.format(_get_default_variables_dict()["CheckpointFreq"],
                                    _get_default_variables_dict()["CheckpointFreq"][0],
                                    _get_default_variables_dict()["CheckpointFreq"][1]),
        "ConsoleFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or '
                       '[Generate_data_bool , steps_per_data_output_int], '
                       'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                       'Controls the output to the "console” or log file, which prints the '
                       'acceptance statistics, and run timing info. In addition, instantaneously-selected'
                       'thermodynamic properties will be output to this file.  If bool is True, '
                       'this enables outputting the console data at the integer frequency '
                       '(set steps_per_data_output_int), while "False" disables outputting the console '
                       'data file. '
                       ''.format(_get_default_variables_dict()["ConsoleFreq"],
                                 _get_default_variables_dict()["ConsoleFreq"][0],
                                 _get_default_variables_dict()["ConsoleFreq"][1]),
        "BlockAverageFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' 
                            '[Generate_data_bool , steps_per_data_output_int], '
                            'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                            'Controls the block averages output of selected thermodynamic properties. '
                            'Block averages are averages of thermodynamic values of interest for chunks of the '
                            'simulation (for post-processing of averages or std. dev. in those values).'
                            'If bool is True, this enables outputting the block averaging data/file at the '
                            'integer frequency (set steps_per_data_output_int), '
                            'while "False" disables outputting the block averaging data/file.'
                            ''.format(_get_default_variables_dict()["BlockAverageFreq"],
                                      _get_default_variables_dict()["BlockAverageFreq"][0],
                                      _get_default_variables_dict()["BlockAverageFreq"][1]),
        "HistogramFreq": 'Output Frequency (all ensembles): list [bool , int (> 0)] or ' 
                         '[Generate_data_bool , steps_per_data_output_int], ' 
                         'default = {} or [{} , set via formula based on the number of RunSteps or {} max]. ' 
                         'Controls the histograms. Histograms are a binned listing of observation frequency ' 
                         'for a specific thermodynamic variable. In the GOMC code, they also control the output '
                         'of a file containing energy/molecule samples, ' 
                         'which is only used for the "GCMC" ensemble simulations for histogram reweighting purposes.' 
                         'If bool is True, this enables outputting the data to the histogram data at the '
                         'integer frequency (set steps_per_data_output_int), ' 
                         'while "False" disables outputting the histogram data.'
                         ''.format(_get_default_variables_dict()["HistogramFreq"],
                                   _get_default_variables_dict()["HistogramFreq"][0],
                                   _get_default_variables_dict()["HistogramFreq"][1]),

        # Histogram data
        "DistName": 'Histogram Output (all ensembles): str (NO SPACES), default = {}. ' 
                    'Short phrase which will be combined with RunNumber and RunLetter '
                    'to use in the name of the binned histogram for molecule distribution.'
                    ''.format(_get_default_variables_dict()["DistName"]),
        "HistName": 'Histogram Output (all ensembles): str (NO SPACES), default = {}. ' 
                    'Short phrase, which will be combined with RunNumber and RunLetter, '
                    'to use in the name of the energy/molecule count sample file.' 
                    ''.format(_get_default_variables_dict()["HistName"]),
        "RunNumber": 'Histogram Output (all ensembles): int  ( > 0 ), default = {}. ' 
                     'Sets a number, which is a part of DistName and HistName file name.'
                     ''.format(_get_default_variables_dict()["RunNumber"]),
        "RunLetter": 'Histogram Output (all ensembles): str (1 alphabetic character only), default = {}. '
                     'Sets a letter, which is a part of DistName and HistName file name.'
                     ''.format(_get_default_variables_dict()["RunLetter"]),
        "SampleFreq": 'Histogram Output (all ensembles): int ( > 0 ), default = {}. ' 
                      'The number of steps per histogram sample or frequency.'
                      ''.format(_get_default_variables_dict()["SampleFreq"]),

        # Data output for the console and bulk properties calculations
        "OutEnergy": 'Output Data (all ensembles): [bool, bool], default = {}.   '
                     'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                     'This outputs the energy data into the block averages and console output/log files.'
                     ''.format(_get_default_variables_dict()["OutEnergy"]),
        "OutPressure": 'Output Data (all ensembles): [bool, bool], default = {}.   '
                       'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                       'This outputs the pressure data into the block averages and console output/log files.'
                       ''.format(_get_default_variables_dict()["OutPressure"]),
        "OutMolNumber": 'Output Data (all ensembles): [bool, bool], default = {}.   '
                        'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                        'This outputs the number of molecules data into the block averages and console output/log files.'
                        ''.format(_get_default_variables_dict()["OutMolNumber"]),
        "OutDensity": 'Output Data (all ensembles): [bool, bool], default = {}.   '
                      'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                      'This outputs the density data into the block averages and console output/log files.'
                      ''.format(_get_default_variables_dict()["OutDensity"]),
        "OutVolume": 'Output Data (all ensembles): [bool, bool], default = {}.   '
                     'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                     'This outputs the volume data into the block averages and console output/log files.'
                     ''.format(_get_default_variables_dict()["OutVolume"]),
        "OutSurfaceTension": 'Output Data (all ensembles): [bool, bool], default = {}. ' 
                             'The list provides the booleans to [block_averages_bool, console_output_bool]. '
                             'This outputs the surface tension data into the block averages and console '
                             'output/log files.'
                             ''.format(_get_default_variables_dict()["OutSurfaceTension"]),

        # free energy calculation in NVT and NPT ensembles.
        "FreeEnergyCalc": 'Free Energy Calcs (NVT and NPT only): list [bool , int (> 0)] or '
                          '[Generate_data_bool , steps_per_data_output_int], default = {}. ' 
                          'bool = True enabling free energy calculation during the simulation, false disables '
                          'the calculation. The int/step frequency sets the frequency of calculating the free energy.'
                          ''.format(_get_default_variables_dict()["FreeEnergyCalc"]),
        "MoleculeType": 'Free Energy Calcs (NVT and NPT only): list [str , int (> 0)] or '  
                        '["residue_name" , residue_ID], ' 
                        'The user must set this variable as there is no working default (default = {}). ' 
                        'Note: ONLY 4 characters can be used for the string (i.e., "residue_name"). ' 
                        'Sets the solute molecule kind (residue name) and molecule number (residue ID), '  
                        'which absolute solvation free will be calculated for.'
                        ''.format(_get_default_variables_dict()["MoleculeType"]),
        "InitialState": 'Free Energy Calcs (NVT and NPT only): int (>= 0), '
                        'The user must set this variable as there is no working default (default = {}). ' 
                        'The index of LambdaCoulomb and LambdaVDW vectors. Sets the index of the' 
                        'LambdaCoulomb and LambdaVDW vectors, to determine the simulation lambda value for'
                        'VDW and Coulomb interactions. ' 
                        'WARNING : This must an integer within the vector count of the LambdaVDW and LambdaCoulomb, ' 
                        'in which the counting starts at 0.  '
                        ''.format(_get_default_variables_dict()["InitialState"]),
        "LambdaVDW": 'Free Energy Calcs (NVT and NPT only): list of floats (0 <= floats <= 1), ' 
                     'The user must set this variable as there is no working default (default = {}). ' 
                     'Lambda values for VDW interaction in ascending order. Sets the intermediate '
                     'lambda states to which solute-solvent VDW interactions are scaled. '
                     'WARNING : This list must be the same length as the "LambdaCoulomb" list length. '
                     'WARNING : All lambda values must be stated in the ascending order, otherwise '
                     'the program will terminate.  '
                     'Example of ascending order 1: [0.1, 1.0,]  '
                     'Example of ascending order 2: [0.1, 0.2, 0.4, 0.9] '
                     ''.format(_get_default_variables_dict()["LambdaVDW"]),
        "LambdaCoulomb": 'Free Energy Calcs (NVT and NPT only):  list of floats (0 <= floats <= 1), '
                         'The user must set this variable as there is no working default (default = {}). ' 
                         'Lambda values for Coulombic interaction in ascending order. Sets the intermediate '
                         'lambda states to which solute-solvent Coulombic interactions are scaled. ' 
                         'GOMC defauts to the "LambdaVDW" values for the Coulombic interaction '
                         'if no "LambdaCoulomb" variable is set. '
                         'WARNING : This list must be the same length as the "LambdaVDW" list length. '
                         'WARNING : All lambda values must be stated in the ascending order, otherwise '
                         'the program will terminate.  '
                         'Example of ascending order 1: [0.1, 1.0,]  '
                         'Example of ascending order 2: [0.1, 0.2, 0.4, 0.9] '
                         ''.format(_get_default_variables_dict()["LambdaCoulomb"]),
        "ScaleCoulomb": 'Free Energy Calcs (NVT and NPT only): bool, default = {}, '
                        'Determines to scale the Coulombic interaction non-linearly (soft-core scheme) or not. '
                        'True if the Coulombic interaction needs to be scaled non-linearly. '
                        'False if the Coulombic interaction needs to be scaled linearly. '
                        ''.format(_get_default_variables_dict()["ScaleCoulomb"]),
        "ScalePower": 'Free Energy Calcs (NVT and NPT only): int (>= 0), default = {}, '
                      'The p value in the soft-core scaling scheme, where the distance between '
                      'solute and solvent is scaled non-linearly.'
                      ''.format(_get_default_variables_dict()["ScalePower"]),
        "ScaleAlpha": 'Free Energy Calcs (NVT and NPT only): int or float (>= 0), default = {}, '
                      'The alpha value in the soft-core scaling scheme, where the distance '
                      'between solute and solvent is scaled non-linearly.'
                      ''.format(_get_default_variables_dict()["ScaleAlpha"]),
        "MinSigma": 'Free Energy Calcs (NVT and NPT only): int or float (>= 0), default = {}, '
                    'The minimum sigma value in the soft-core scaling scheme, where the '
                    'distance between solute and solvent is scaled non-linearly.' 
                    ''.format(_get_default_variables_dict()["MinSigma"]),


        # moves without MEMC
        "DisFreq": 'Std. MC moves (all ensembles)                     : ' 
                   'int or float (0 <= value <= 1), default are specific for each ' 
                   'ensemble {}. '
                   'Fractional percentage at which the displacement move will occur ' 
                   '(i.e., fraction of displacement moves). Note: all of the move types'  
                   'are not available in for every ensemble. Note: all of the move fractions' 
                   'must sum to 1, or the control file writer will fail.  '
                   ''.format(_get_default_variables_dict()["DisFreq"]),
        "RotFreq": 'Std. MC moves (all ensembles)                     : '
                   'int or float (0 <= value <= 1), default are specific for each ' 
                   'ensemble {}. ' 
                   'Fractional percentage at which the rotation move will occur ' 
                   '(i.e., fraction of rotation moves). Note: all of the move types' 
                   'are not available in for every ensemble. Note: all of the move fractions' 
                   'must sum to 1, or the control file writer will fail.  '
                   ''.format(_get_default_variables_dict()["RotFreq"]),
        "IntraSwapFreq": 'Std. MC moves (all ensembles)                     : ' 
                         'int or float (0 <= value <= 1), default are specific for each ' 
                         'ensemble {}. ' 
                         'Fractional percentage at which the molecule will be removed from a ' 
                         'box and inserted into the same box using coupled-decoupled configurational-bias'
                         'algorithm. (i.e., fraction of intra-molecule swap moves). Note: all of the move types' 
                         'are not available in for every ensemble. Note: all of the move fractions' 
                         'must sum to 1, or the control file writer will fail.  '
                          ''.format(_get_default_variables_dict()["IntraSwapFreq"]),
        "SwapFreq": 'Std. MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' 
                    'int or float (0 <= value <= 1), default are specific for each '
                    'ensemble {}. ' 
                    'For Gibbs and Grand Canonical (GC) ensemble runs only: Fractional ' 
                    'percentage at which molecule swap move will occur using coupled-decoupled '
                    'configurational-bias. (i.e., fraction of molecule swaps moves). Note: all of the move types'  
                    'are not available in for every ensemble. Note: all of the move fractions'  
                    'must sum to 1, or the control file writer will fail.  '
                     ''.format(_get_default_variables_dict()["SwapFreq"]),
        "RegrowthFreq": 'Std. MC moves (all ensembles)                     : '
                        'int or float (0 <= value <= 1), default are specific for each ' 
                        'ensemble {}. ' 
                        'Fractional percentage at which part of the molecule will be ' 
                        'deleted and then regrown using coupled- decoupled configurational-bias algorithm ' 
                        '(i.e., fraction of molecular growth moves). Note: all of the move types'  
                        'are not available in for every ensemble. Note: all of the move fractions' 
                        'must sum to 1, or the control file writer will fail.  '
                         ''.format(_get_default_variables_dict()["RegrowthFreq"]),
        "CrankShaftFreq": 'Std. MC moves (all ensembles)                     : ' 
                          'int or float (0 <= value <= 1), default are specific for each '
                          'ensemble {}. '
                          'Fractional percentage at which crankshaft move will occur. ' 
                          'In this move, two atoms that are forming angle or dihedral are selected randomly and ' 
                          'form a shaft. Then any atoms or group that are within these two selected atoms, will ' 
                          'rotate around the shaft to sample intra-molecular degree of freedom ' 
                          '(i.e., fraction of crankshaft moves). Note: all of the move types' 
                          'are not available in for every ensemble. Note: all of the move fractions' 
                          'must sum to 1, or the control file writer will fail.  '
                           ''.format(_get_default_variables_dict()["CrankShaftFreq"]),
        "VolFreq": 'Std. MC moves (only  GEMC_NPT  and  NPT )         : ' 
                   'int or float (0 <= value <= 1), default are specific for each ' 
                   'ensemble {}. Fractional percentage at  which a volume move will occur '
                   '(i.e., fraction of Volume moves). ' 
                   'Note: all of the move types are not available in for every ensemble. '
                   'Note: all of the move fractions must sum to 1, or the control file writer will fail. '
                    ''.format(_get_default_variables_dict()["VolFreq"]),
        "MultiParticleFreq": 'Std. MC moves (all ensembles)                     : ' 
                             'int or float (0 <= value <= 1), default are specific for each '
                             'ensemble {}. ' 
                             'Fractional percentage at which multi-particle move will ' 
                             'occur. In this move, all molecules in the selected simulation box will be rigidly ' 
                             'rotated or displaced simultaneously, along the calculated torque or force '
                             'respectively (i.e., fraction of multi-particle moves). Note: all of the move types' 
                             'are not available in for every ensemble. Note: all of the move fractions'  
                             'must sum to 1, or the control file writer will fail.  '
                              ''.format(_get_default_variables_dict()["MultiParticleFreq"]),

        # MEMC moves
        "IntraMEMC-1Freq": 'MEMC MC moves (all ensembles)                     : ' 
                           'int or float (0 <= value <= 1), default are specific for each ' 
                           'ensemble {}. ' 
                           'Fractional percentage at which specified number of small molecule kind will be ' 
                           'exchanged with a specified large molecule kind in defined sub-volume within ' 
                           'same simulation box.  This move need additional information such as ' 
                           'ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.' 
                           'Note: all of the move types are not available in for every ensemble.' 
                           'Note: all of the move fractions must sum to 1, or the control file writer will fail.  '
                            ''.format(_get_default_variables_dict()["IntraMEMC-1Freq"]),
        "MEMC-1Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : '
                      'int or float (0 <= value <= 1), default are specific for each '
                      'ensemble {}. '
                      'Fractional percentage at which specified number of small molecule kind will be exchanged '
                      'with a specified large molecule kind in defined sub-volume, between simulation boxes. '
                      'This move needs additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, ' 
                      'and ExchangeLargeKind.'
                      'Note: all of the move types are not available in for every ensemble.' 
                      'Note: all of the move fractions must sum to 1, or the control file writer will fail.  '
                       ''.format(_get_default_variables_dict()["MEMC-1Freq"]),
        "IntraMEMC-2Freq": 'MEMC MC moves (all ensembles)                     : ' 
                           'int or float (0 <= value <= 1), default are specific for each ' 
                           'ensemble {}. ' 
                           'Fractional percentage at which specified number of small molecule kind '
                           'will be exchanged with a specified large molecule kind in defined sub-volume '
                           'within same simulation box. Backbone of small and large molecule kind will be '
                           'used to insert the large molecule more efficiently. This move need additional '
                           'information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, '
                           'ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
                           'Note: all of the move types are not available in for every ensemble.' 
                           'Note: all of the move fractions must sum to 1, or the control file writer will fail.  '
                            ''.format(_get_default_variables_dict()["IntraMEMC-2Freq"]),
        "MEMC-2Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : '
                      'int or float (0 <= value <= 1), default are specific for each ' 
                      'ensemble {}. ' 
                      'Fractional percentage at which specified number of small molecule kind will be '
                      'exchanged with a specified large molecule kind in defined sub-volume, between'
                      'simulation boxes. Backbone of small and large molecule kind will be used to insert '
                      'the large molecule more efficiently. ' 
                      'This move needs additional information such as ExchangeVolumeDim, ExchangeRatio, ' 
                      'ExchangeSmallKind, ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
                      'Note: all of the move types are not available in for every ensemble.' 
                      'Note: all of the move fractions must sum to 1, or the control file writer will fail.  '
                       ''.format(_get_default_variables_dict()["MEMC-2Freq"]),
        "IntraMEMC-3Freq": 'MEMC MC moves (all ensembles)                     : ' 
                           'int or float (0 <= value <= 1), default are specific for each ' 
                           'ensemble {}. ' 
                           'Fractional percentage at which specified number of small molecule kind will be '
                           'exchanged with a specified large molecule kind in defined sub-volume within same '
                           'simulation box. Specified atom of the large molecule kind will be used to insert '
                           'the large molecule using coupled-decoupled configurational-bias. This move needs '
                           'additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, '
                           'ExchangeLargeKind, and LargeKindBackBone. '
                           'Note: all of the move types are not available in for every ensemble.' 
                           'Note: all of the move fractions must sum to 1, or the control file writer will fail. '
                            ''.format(_get_default_variables_dict()["IntraMEMC-3Freq"]),
        "MEMC-3Freq": 'MEMC MC moves (only GEMC_NPT, GEMC_NVT, and GCMC) : ' 
                      'int or float (0 <= value <= 1), default are specific for each ' 
                      'ensemble {}. '
                      'Fractional percentage at which specified number of small molecule kind will be exchanged '
                      'with a specified large molecule kind in defined sub-volume, between simulation boxes. '
                      'Specified atom of the large molecule kind will be used to insert the large molecule '
                      'using coupled-decoupled configurational-bias. This move need additional information '
                      'such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, ExchangeLargeKind, '
                      'and LargeKindBackBone. '
                      'Note: all of the move types are not available in for every ensemble.' 
                      'Note: all of the move fractions must sum to 1, or the control file writer will fail.  '
                      ''.format(_get_default_variables_dict()["MEMC-3Freq"]),

        # MEMC move parameters
        "ExchangeVolumeDim": 'MEMC parameters (all ensembles)                   : ' 
                             'list of 3 floats or integers ' 
                             '[int or float (> 0), int or float (> 0), int or float (> 0)]'
                             ' or [X-dimension, Y-dimension, Z-dimension)], '
                             'default = {}. '
                             'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, the exchange '
                             'subvolume must be defined. The exchange sub-volume is defined as an orthogonal box ' 
                             'with x, y, and z-dimensions, where small molecule/molecules kind will be selected ' 
                             'from to be exchanged with a large molecule kind. ' 
                             'Note: Currently, the X and Y dimension cannot be set independently (X = Y = max(X, Y)). '
                             'Note: A heuristic for setting good values of the x, y, and z-dimensions is to use'
                             'the geometric size of the large molecule plus 1-2 Å in each dimension. '
                             'Note: In case of exchanging 1 small molecule kind with 1 large molecule kind in '
                             'IntraMEMC-2, IntraMEMC-3, MEMC-2, MEMC-3 Monte Carlo moves, the sub-volume '
                             'dimension has no effect on acceptance rate. '
                             ''.format(_get_default_variables_dict()["ExchangeVolumeDim"]),
        "MEMC_DataInput": 'MEMC parameters (availablity based on selelection): nested lists, ' +
                          'default = {}.  '.format(_get_default_variables_dict()["MEMC_DataInput"]) +
                          'Enter data as a list with some sub-lists as follows: '
                          '[[ExchangeRatio_int (> 0), ExchangeLargeKind_str, ' 
                          '[LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ], '
                          'ExchangeSmallKind_str, '
                          '[SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]], ..., ' 
                          '[ExchangeRatio_int (> 0), ExchangeLargeKind_str, '
                          '[LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ], ' 
                          'ExchangeSmallKind_str, '
                          '[SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]. ' 
                          'NOTE: CURRENTLY ALL THESE INPUTS NEED TO BE SPECIFIED, REGARDLESS OF THE MEMC TYPE ' 
                          'SELECTION. IF THE SmallKindBackBone or LargeKindBackBone IS NOT REQUIRED FOR THE '
                          'MEMC TYPE, '
                          'None CAN BE USED IN PLACE OF A STRING. ' 
                          'Note: These strings must match the residue in the psf and psb files or it will fail. ' 
                          'It is recommended that the user print the Charmm object psf and pdb files '
                          'and review the residue names that match the atom name before using the in '
                          'the  MEMC_DataInput variable input.'
                          'Note: see the below data explanations for the ExchangeRatio, ExchangeSmallKind, '
                          'ExchangeLargeKind, LargeKindBackBone, SmallKindBackBone. ' 
                          "Example 1 (MEMC-1) : [ [1, 'WAT', [None, None], 'wat', [None, None]] , "
                          "[1, 'WAT', [None, None], 'wat', [None, None]] . "
                          "Example 2 (MEMC-2): [ [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1' ]] , "
                          " [1, 'WAT', ['H1', 'H2'], 'wat', ['H1', 'H2' ]] . "
                          "Example 3 (MEMC-3) : [ [2, 'WAT', 'O1', 'H1'], 'wat', [None, None]] , "
                          "[2, 'WAT', ['H1', 'H2'], 'wat', [None, None]] .\n"
                          '\t\t\t\t\t\t\t\t\t\t\t\t\t --- ExchangeRatio     = MEMC parameters (all ensembles): ' 
                          'int (> 0), default = None. The Ratio of exchanging ' 
                          'small molecule/molecules with 1 large molecule. '
                          'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, ' 
                          'the exchange ratio must be defined. ' 
                          'The exchange ratio defines how many small molecule will be ' 
                          'exchanged with 1 large molecule. For each large-small molecule pairs, ' 
                          'one exchange ratio must be defined. \n' 
                          '\t\t\t\t\t\t\t\t\t\t\t\t\t --- ExchangeSmallKind = MEMC parameters (all ensembles): ' 
                          'str, default = None. The small molecule ' 
                          'kind (resname) to be exchanged. ' 
                          'Note: ONLY 4 characters can be used for the strings. ' 
                          'To use all variation of MEMC and Intra-MEMC Monte Carlo moves, '
                          'the small molecule kind to be exchanged with a large molecule '
                          'kind must be defined. Multiple small molecule kind can be specified.  \n'
                          '\t\t\t\t\t\t\t\t\t\t\t\t\t --- ExchangeLargeKind = MEMC parameters (all ensembles): ' 
                          'str, default = None. The large molecule ' 
                          'kind (resname) to be exchanged. '
                          'Note: ONLY 4 characters can be used for the strings. '
                          'To use all variations of MEMC and Intra-MEMC Monte Carlo moves, ' 
                          'the large molecule kind to be exchanged with small molecule ' 
                          'kind must be defined. Multiple large molecule kind can be specified. \n' 
                          '\t\t\t\t\t\t\t\t\t\t\t\t\t --- LargeKindBackBone = MEMC parameters (all ensembles): '
                          'list [str, str] or [None, None], default = None '
                          'Note: ONLY 4 characters can be used for the strings. ' 
                          'The [None, None] values can only be used if that MEMC type does not require them. '
                          'The strings for the the atom name 1 and atom name 2 that belong to the large '
                          'molecule’s backbone (i.e., [str_for_atom_name_1, str_for_atom_name_2]) '
                          'To use MEMC-2, MEMC-3, IntraMEMC-2, and IntraMEMC-3 Monte Carlo moves, the ' 
                          'large molecule backbone must be defined. The backbone of the molecule is defined ' 
                          'as a vector that connects two atoms belong to the large molecule. The large ' 
                          'molecule backbone will be used to align the sub-volume in MEMC-2 and IntraMEMC-2 '
                          'moves, while in MEMC-3 and IntraMEMC-3 moves, it uses the atom name to start ' 
                          'growing the large molecule using coupled-decoupled configurational-bias. For ' 
                          'each large-small molecule pairs, two atom names must be defined. ' 
                          'Note: all atom names in the molecule must be unique. ' 
                          'Note: In MEMC-3 and IntraMEMC-3 Monte Carlo moves, both atom names must be same, ' 
                          'otherwise program will be terminated. ' 
                          'Note: If the large molecule has only one atom (mono atomic molecules), '
                          'same atom name must be used for str_for_atom_name_1 and str_for_atom_name_2 ' 
                          'of the LargeKindBackBone.  \n' 
                          '\t\t\t\t\t\t\t\t\t\t\t\t\t --- SmallKindBackBone = MEMC parameters (all ensembles): '
                          ' list [str, str] or [None, None], default = None '
                          'Note: ONLY 4 characters can be used for the strings. ' 
                          'The [None, None] values can only be used if that MEMC type does not require them.' 
                          'The strings for the the atom name 1 and atom name 2 that belong to the small ' 
                          'molecule’s backbone (i.e., [str_for_atom_name_1, str_for_atom_name_2]) ' 
                          'To use MEMC-2, and IntraMEMC-2 Monte Carlo moves, the small molecule backbone ' 
                          'must be defined. The backbone of the molecule is defined as a vector that '
                          'connects two atoms belong to the small molecule and will be used to align the ' 
                          'sub-volume. For each large-small molecule pairs, two atom names must be defined. '
                          'Note: all atom names in the molecule must be unique. ' 
                          'Note: If the small molecule has only one atom (mono atomic molecules), same atom '
                          'name must be used str_for_atom_name_1 and str_for_atom_name_2 ' 
                          'of the SmallKindBackBone. ',

        # ******************************************************************************************************
        # Definitions in this function are copied to a large extent from the GOMC manual release version 2.60 (end)
        # insert citation here:
        # ******************************************************************************************************
    }
    if description:
        return valid_input_variables
    else:
        return list(valid_input_variables.keys())


def _get_default_variables_dict():
    """
    Provides a dictionary of the default variables inputs and their default settings (user optional).

    Returns
    ---------
    default_input_variables_dict : dict
        Provides a dict of the default variables inputs (user optional)

    """

    default_input_variables_dict = {
        "Restart": False,
        "RestartCheckpoint": False,
        "PRNG": "RANDOM",
        "ParaTypeCHARMM": True,
        "ParaTypeMie": False,
        "ParaTypeMARTINI": False,
        "RcutCoulomb_box_0": None,
        "RcutCoulomb_box_1": None,
        "Pressure": 1.01325,
        "Rcut": 10,
        "RcutLow":  1,
        "LRC": True,
        "Exclude": '1-3',
        "coul_1_4_scaling": None,
        "Potential": 'VDW',
        "Rswitch": 9,
        "ElectroStatic": True,
        "Ewald":  True,
        "CachedFourier": False,
        "Tolerance": 0.00001,
        "Dielectric": 15,
        "PressureCalc": [True, 10000],
        "EqSteps": 1000000,
        "AdjSteps": 1000,
        "VDWGeometricSigma": False,
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
        "CoordinatesFreq": [True, 1000000],
        "RestartFreq": [True, 1000000],
        "CheckpointFreq": [True, 1000000],
        "ConsoleFreq": [True, 10000],
        "BlockAverageFreq": [True, 10000],
        "HistogramFreq": [True, 10000],

        # Histogram data
        "DistName": "dis",
        "HistName": "his",
        "RunNumber": 1,
        "RunLetter": "a",
        "SampleFreq": 500,

        # Data output for the console and bulk properties calculations
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
        "DisFreq":           {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.20, 'GEMC_NPT': 0.19, 'GCMC': 0.15},
        "RotFreq":           {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.20, 'GEMC_NPT': 0.20, 'GCMC': 0.15},
        "IntraSwapFreq":     {'NVT': 0.30, 'NPT': 0.29, 'GEMC_NVT': 0.10, 'GEMC_NPT': 0.10, 'GCMC': 0.10},
        "SwapFreq":          {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.20, 'GEMC_NPT': 0.20, 'GCMC': 0.35},
        "RegrowthFreq":      {'NVT': 0.30, 'NPT': 0.30, 'GEMC_NVT': 0.20, 'GEMC_NPT': 0.20, 'GCMC': 0.15},
        "CrankShaftFreq":    {'NVT': 0.10, 'NPT': 0.10, 'GEMC_NVT': 0.10, 'GEMC_NPT': 0.10, 'GCMC': 0.10},
        "VolFreq":           {'NVT': 0.00, 'NPT': 0.01, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.01, 'GCMC': 0.00},
        "MultiParticleFreq": {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        # MEMC moves
        "IntraMEMC-1Freq":   {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        "MEMC-1Freq":        {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        "IntraMEMC-2Freq":   {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        "MEMC-2Freq":        {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        "IntraMEMC-3Freq":   {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
        "MEMC-3Freq":        {'NVT': 0.00, 'NPT': 0.00, 'GEMC_NVT': 0.00, 'GEMC_NPT': 0.00, 'GCMC': 0.00},
    }

    return default_input_variables_dict


def check_valid_ensemble_files(ensemble_type, testing_ensemble_files_list):
    """
    Checks if all the required ensemble inputs are provided,
        and provides a list of the bad variables in the printed output.

    Parameters
    ----------
    ensemble_type : str, valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
        The ensemble type of the simulation.
    testing_ensemble_files_list  list
        A list containing the required ensemble
        files variables, which will be tested for to see if they are valid.

    Returns
    ---------
    bool
        True is all variables are valid, False otherwise
    """

    bad_key_inputs_List = []

    req_ensemble_files_set = set(_get_required_data(description=False))
    testing_ensemble_files_set = set(testing_ensemble_files_List)

    extra = testing_ensemble_files_set - req_ensemble_files_set
    missing = req_ensemble_files_set - testing_ensemble_files_set

    if len(extra) != 0:
        bad_key_inputs_List.extend(extra)
    if len(missing) != 0:
        bad_key_inputs_List.extend(missing)

    if not bool(missing or extra):
        return True
    else:
        return False


def print_required_input(description=False):
    """
    Prints the required ensemble arguments with an optional description based on the ensemble type

    Parameters
    ----------
    description :  bool, default = False.
        If True, it prints the descriptions of the required ensemble inputs (i.e. dict),
        If False, it only prints the required ensemble inputs without the descriptions (i.e. list)

    Returns
    ---------
    Prints the required ensemble arguments with an optional description based on the ensemble type
    """

    required_data_dict = _get_required_data(description=True)
    required_data_list = _get_required_data(description=False)
    ensemble_has_all_valid_required_data = True
    required_data = _get_required_data()

    for iter in range(0, len(required_data_list)):
        if required_data_list[iter] not in required_data:
            ensemble_has_all_valid_required_data = False

    if ensemble_has_all_valid_required_data:
        for iter_2 in range(0, len(required_data_list)):
            required_data_iter = required_data_list[iter_2]
            if description is False:
                print("{:10s}:    {}".format(str(iter_2), str(required_data_iter)))
            elif description is True:
                print("{:10s}:    {:30s}    {}".format(str(iter_2), str(required_data_iter),
                                                       str(required_data_dict[required_data_iter])))
    else:
        print("ERROR: Some files in this ensemble are not in the required file list")


def check_valid_ensemble_input_variables(ensemble_type, testing_input_variables_list):
    """
    Checks if all the input variables (user optional) inputs are valid for the given
        ensemble, and provides a list of the bad variables in the printed output.

    Parameters
    ----------
    ensemble_type : str, valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
        The ensemble type of the simulation.
    testing_input_variables_list : list
        List containing the optional ensemble input variables which will be
        tested for to see if they are valid.

    Returns
    ---------
    bool:
        Returns a bool (True or False) depending on if all variables
        are valid or not.
    """

    bad_key_inputs_list = []

    valid_input_variables_list = _get_possible_ensemble_input_variables(ensemble_type)
    ensemble_has_valid_input_variables_list = True
    for iter in range(0, len(testing_input_variables_list)):
        if testing_input_variables_list[iter] not in valid_input_variables_list:
            bad_key_inputs_list.append(testing_input_variables_list[iter])
            ensemble_has_valid_input_variables_list = False

    if ensemble_has_valid_input_variables_list:
        return [True, bad_key_inputs_list]

    else:
        return [False, bad_key_inputs_list]


def print_valid_ensemble_input_variables(ensemble_type, description=False):
    """
    Prints the arguments for optional variables brief description based on the ensemble type

    Parameters
    ----------
    ensemble_type = str, valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
        The ensemble type of the simulation.
    description =  bool, default = False.
        If True, it prints the descriptions of the optional variable ensemble inputs (i.e. dict),
        If False, it only prints the  optional variable ensemble inputs without the
        descriptions (i.e. list)

    Returns
    ---------
    Prints the arguments for optional variables brief description based on the ensemble type
    """

    valid_input_variables_dict = _get_all_possible_input_variables(description=True)
    ensemble_has_all_valid_input_variables = True
    all_valid_input_variables = _get_all_possible_input_variables()

    valid_input_variables_list = _get_possible_ensemble_input_variables(ensemble_type)

    for iter in range(0, len(valid_input_variables_list)):
        if valid_input_variables_list[iter] not in all_valid_input_variables:
            ensemble_has_all_valid_input_variables = False

    if ensemble_has_all_valid_input_variables:
        for iter_2 in range(0, len(valid_input_variables_list)):
            ensemble_kwarg_iter = valid_input_variables_list[iter_2]
            if description is False:
                print("{:10s}:    {}".format(str(iter_2), str(ensemble_kwarg_iter)))
            elif description is True:
                print("{:10s}:    {:30s}    {}".format(str(iter_2), str(ensemble_kwarg_iter),
                                                       str(valid_input_variables_dict[ensemble_kwarg_iter])))
    else:
        print("ERROR: Some input_variables in the ensemble are not in the main input_variables list")


def _get_possible_ensemble_input_variables(ensemble_type):
    """
    Provides list of the possible optional input variables based on the ensemble type

    Parameters
    ----------
    ensemble_type : str, valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
        The ensemble type of the simulation.

    Returns
    ---------
    valid_input_variables_list : list
        A list possible optional input variables for the provided ensemble type.
    """
    histogram_output_variables_list = ["DistName", "HistName", "RunNumber", "RunLetter", "SampleFreq"]

    cbmc_variables_list = ["CBMC_First", "CBMC_Nth", "CBMC_Ang", "CBMC_Dih"]

    output_freq_variables_list = ["OutputName", "CoordinatesFreq", "RestartFreq", "CheckpointFreq",
                                  "ConsoleFreq", "BlockAverageFreq", "HistogramFreq"]

    output_data_variables_list = ["OutEnergy", "OutPressure", "OutMolNumber", "OutDensity",
                                  "OutVolume", "OutSurfaceTension"]

    std_mc_moves_variables_list = ["DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                                   "CrankShaftFreq", "VolFreq", "MultiParticleFreq"]

    memc_mc_moves_variables_list = ["IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                                    "IntraMEMC-3Freq", "MEMC-3Freq",
                                    "ExchangeVolumeDim", "MEMC_DataInput"]

    basic_sim_info_variables_list = ["Restart", "RestartCheckpoint", "PRNG",
                                     "ParaTypeCHARMM", "ParaTypeMie", "ParaTypeMARTINI",
                                     "RcutCoulomb_box_0",
                                     "Pressure",
                                     "Rcut", "RcutLow", "LRC", "Exclude", "Potential", "Rswitch",
                                     "ElectroStatic", "Ewald", "CachedFourier", "Tolerance",
                                     "Dielectric", "PressureCalc", "EqSteps", "AdjSteps", "VDWGeometricSigma",
                                     "useConstantArea"]

    if ensemble_type in ['NPT', 'NVT']:
        extra_sim_info_variables_list = []

        free_energy_variables_list = ["FreeEnergyCalc", "MoleculeType", "InitialState", "LambdaVDW",
                                      "LambdaCoulomb", "ScaleCoulomb", "ScalePower", "ScaleAlpha", "MinSigma"]

    elif ensemble_type in ['GEMC_NPT', 'GEMC_NVT']:
        extra_sim_info_variables_list = ["RcutCoulomb_box_1", "FixVolBox0"]

        free_energy_variables_list = []  # always empty for GEMC

    elif ensemble_type == 'GCMC':
        extra_sim_info_variables_list = ["ChemPot", "Fugacity"]

        free_energy_variables_list = []  # always empty for GCMC

    if ensemble_type in ['NPT', 'NVT', 'GEMC_NPT', 'GEMC_NVT', 'GCMC']:
        valid_input_variables_list = basic_sim_info_variables_list + extra_sim_info_variables_list \
                                     + cbmc_variables_list \
                                     + output_freq_variables_list + histogram_output_variables_list \
                                     + output_data_variables_list + free_energy_variables_list \
                                     + std_mc_moves_variables_list + memc_mc_moves_variables_list
    else:
        warn('WARNINR: The ensemble_type selected for the _get_possible_ensemble_input_variables '
              'function is not valid.')
        valid_input_variables_list = None

    return valid_input_variables_list

class GOMCControl():
    def __init__(self, charmm_object, ensemble_type, RunSteps, Temperature, input_variables_dict=None
                 ):
        """
        Constructs the GOMC control file with user selected variables.
        The selected variables and Class attributes are mostly or nearly identical to the
        GOMC command names. For many ensembles, the user may use the default input_variables_dict
        variables, unless they are required for the specific ensemble (Example:
        the GCMC ensemble requires the user to input the Chempot or Fugacity variables,
        or the build will fail.)

        Default settings for the GOMC configuration files are based upon
        an educated guess, which should result in appropriate sampling for a
        given ensemble/simulation type. However, there is no guarantee that
        the default setting will provide the best or adequate sampling for
        the selected system. The user can modify the configuration/control files
        based on the simulation specifics or optimize the system beyond the standard
        settings.  These override options are available via the keyword arguments
        in input_variable_dict.

        Parameters
        ----------
        charmm_object :  Charmm object
            Charmm object is has been parameterized from the selected force field.,
        ensemble_typ : str, ['NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', 'GCMC']
            The ensemble type of the simulation.
        RunSteps : int (>0), must be an integer greater than zero.
            Sets the total number of simulation steps.
        Temperature : float or int (>0), must be an integer greater than zero.
            Temperature of system in Kelvin (K)
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

            # *******************************************************************
            # input_variables_dict options (keys and values) - (start)
            # Note: the input_variables_dict keys are also attributes
            # *******************************************************************
            Restart : boolean, default = False
                Determines whether to restart the simulation from restart file
                (*_restart.pdb and *_restart.psf) or not.
            RestartCheckpoint : boolean, default = False, default = "RANDOM"
                Determines whether to restart the simulation with the checkpoint
                file (checkpoint.dat) or not. Restarting the simulation with checkpoint.dat
                would result in an identical outcome, as if previous simulation was continued.
            PRNG : string or int (>= 0) ("RANDOM" or int)
                PRNG = Pseudo-Random Number Generator (PRNG). There are two (2) options, entering
                the string, "RANDOM", or a integer.
                --- "RANDOM", which selects a random seed number.  This will enter the line
                    "PRNG RANDOM" in the gomc configuration file.
                --- integer, which defines the integer seed number for the simulation. This is
                    equivalent to entering the following two lines in the configuration file:
                    line 1 = PRNG INTSEED
                    line 2 = Random_Seed user_selected_integer.
                Example 1: for a random seed enter the string "RANDOM.
                Example 2: for a specific seed number enter a integer of your choosing.
            ParaTypeCHARMM : boolean, default = True
                True if a CHARMM forcefield, False otherwise.
            ParaTypeMie : boolean, default = False
                True if a Mie forcefield type, False otherwise.
            ParaTypeMARTINI : boolean, default = False
                True if a MARTINI forcefield, False otherwise.
            RcutCoulomb_box_0 : int or float (>= 0), default = None
                Sets a specific radius in box 0 where the short-range electrostatic
                energy will be calculated (i.e., The distance to truncate the
                short-range electrostatic energy in box 0.)
                Note: if None, GOMC will default to the Rcut value
            RcutCoulomb_box_1 : int or float (>= 0), default = None
                Sets a specific radius in box 1 where the short-range electrostatic
                energy will be calculated (i.e., The distance to truncate the
                short-range electrostatic energy in box 0.)
                Note: if None, GOMC will default to the Rcut value
            Pressure : int or float (>= 0), default = 1.01325
                The pressure in bar utilized for the NPT and GEMC_NPT simulations.'
            Rcut : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 10
                Sets a specific radius in Angstroms that non-bonded interaction
                energy and force will be considered and calculated using defined potential function.
                The distance in Angstoms to truncate the LJ, Mie, or other VDW type potential at.
                Note: Rswitch is only used when the "Potential" = SWITCH.
            RcutLow : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 1
                Sets a specific minimum possible distance in Angstroms that reject
                any move that places any atom closer than specified distance.
                The minimum possible distance between any atoms.
                Sets a specific radius in Angstroms that non-bonded interaction
                Note: Rswitch is only used when the "Potential" = SWITCH.
            LRC : boolean, default = True
                If True, the simulation considers the long range tail corrections for the
                non-bonded VDW or dispersion interactions.
                Note: In case of using SHIFT or SWITCH potential functions, LRC will be ignored.
            Exclude : str ["1-2", "1-3", or "1-4"], default = 1-3"
                Note: In CHARMM force field, the 1-4 interaction needs to be considered.
                Choosing "Excude 1-3", will modify 1-4 interaction based on 1-4 parameters
                in parameter file. If a kind force field is used, where 1-4 interaction
                needs to be ignored, such as TraPPE, either Exclude "1-4" needs to be
                chosen or 1-4 parameter needs to be assigned to zero in the parameter file.
                --- "1-2": All interaction pairs of bonded atoms, except the ones that
                    separated with one bond, will be considered and modified using 1-4
                    parameters defined in parameter file.
                --- "1-3": All interaction pairs of bonded atoms, except the ones that
                    separated with one or two bonds, will be considered and modified using
                    1-4 parameters defined in parameter file.
                --- "1-4": All interaction pairs of bonded atoms, except the ones that
                    separated with one, two or three bonds, will be considered using
                    non-bonded parameters defined in parameter file.
            Potential : str, ["VDW", "EXP6", "SHIFT" or "SWITCH"], default = "VDW"
                Defines the potential function type to calculate non-bonded dispersion
                interaction energy and force between atoms.
                ---    "VDW":   Non-bonded dispersion interaction energy and force
                                calculated based on n-6 (Lennard - Jones) equation. This
                                function will be discussed further in the Intermolecular energy
                                and Virial calculation section.
                ---   "EXP6":   Non-bonded dispersion interaction energy and force calculated
                                based on exp-6 (Buckingham potential) equation.
                ---  "SHIFT":   This option forces the potential energy to be zero at Rcut distance.
                --- "SWITCH":   This option smoothly forces the potential energy to be zero at
                                Rcut distance and starts modifying the potential at Rswitch
                                distance. Depending on force field type, specific potential
                                function will be applied.
            Rswitch : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 9
                Note: Rswitch is only used when the SWITCH function is used
                (i.e., "Potential" = SWITCH). The Rswitch distance is in Angstrom. If the
                “SWITCH” function is chosen, Rswitch needs to be defined, otherwise, the
                program will be terminated. When using choosing "SWITCH" as potential function,
                the Rswitch distance defines where the non-bonded interaction energy
                modification is started, which is eventually truncated smoothly at Rcut
                distance.
            ElectroStatic : boolean, default = True
                Considers the coulomb interactions or not. If True, coulomb interactions are
                considered and false if not. Note: To simulate the polar molecule in MARTINI
                force field, ElectroStatic needs to be turn on (i.e., True). The MARTINI force
                field uses short-range coulomb interaction with constant Dielectric of 15.0.
            Ewald : boolean, default = True
                Considers the standard Ewald summation method for electrostatic calculations.
                If True, Ewald summation calculation needs to be considered and false if not.
                Note: By default, GOMC will set ElectroStatic to True if Ewald summation
                method was used to calculate coulomb interaction.
            CachedFourier : boolean, default = False
                Considers storing the reciprocal terms for Ewald summation calculation in
                order to improve the code performance. This option would increase the code
                performance with the cost of memory usage. If True, to store reciprocal
                terms of Ewald summation calculation and False if not.
                Warning: Monte Carlo moves, such as MEMC-1, MEMC-2, MEMC-3,
                IntraMEMC-1, IntraMEMC-2, and IntraMEMC-3 are not support with CachedFourier.
            Tolerance : float (0.0 < float < 1.0), default = 1e-05
                Sets the accuracy in Ewald summation calculation. Ewald separation parameter
                and number of reciprocal vectors for the Ewald summation are determined
                based on the accuracy parameter.
            Dielectric : int or float (>= 0.0), default = 15
                Sets dielectric value used in coulomb interaction when the Martini
                force field is used. Note: In MARTINI force field, Dielectric needs to
                be set to 15.0.
            PressureCalc : list [bool , int (> 0)] or [bool , step_frequency],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Calculate the system pressure or not. bool = True, enables the pressure calculation
                during the simulation, false disables the calculation. The int/step frequency sets the
                frequency of calculating the pressure.
            EqSteps : int (> 0), default = set via formula based on the number of RunSteps or 1M max
                Sets the number of steps necessary to equilibrate the system.
                Averaging will begin at this step.
                Note: In GCMC simulations, the Histogram files will be outputed at EqSteps.
            AdjSteps : int (> 0), default = set via formula based on the number of RunSteps or 1k max
                Sets the number of steps per adjustment of the parameter associated with each move
                (e.g. maximum translate distance, maximum rotation, maximum volume exchange, etc.).
            VDWGeometricSigma: boolean, default = False
                Use geometric mean, as required by OPLS force field, to combining
                Lennard-Jones sigma parameters for different atom types.
                If set to True, GOMC uses geometric mean to combine Lennard-Jones or VDW sigmas.
                Note: The default setting of VDWGeometricSigma is false, which uses the arithmetic
                mean when combining Lennard-Jones or VDW sigma parameters for different atom types.
            useConstantArea : boolean,  default = False
                Changes the volume of the simulation box by fixing the cross-sectional
                area (x-y plane). If True, the volume will change only in z axis,
                If False, the volume of the box will change in a way to maintain the constant
                axis ratio.
            FixVolBox0 : boolean, default = False
                Changing the volume of fluid phase (Box 1) to maintain the constant imposed
                pressure and Temperature, while keeping the volume of adsorbed phase (Box 0) fixed.
                Note: By default, GOMC will set useConstantArea to False if no value was set.
                It means, the volume of the box will change in a way to maintain the constant
                axis ratio.
            ChemPot : dict {str (4 dig limit) , int or float}, default = None
                The chemical potentials in GOMC units of energy, K.
                There is a 4 character limit for the string/residue name since the PDB/PSF
                files have a 4 character limitation and require and exact match in the conf file.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                The name of the residues and their corresponding chemical potential must specified
                for every residue in the system (i.e., {"residue_name" : chemical_potential}).
                Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY
                OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER.
                Example 1 (system with only water):  {"H2O" : -4000} .
                Example 2 (system with water and ethanol):  {"H2O" : -4000, "ETH" : -8000}
            Fugacity : dict {str , int or float (>= 0)}, default = None
                The fugacity in GOMC units of pressure, bar.
                There is a 4 character limit for the string/residue name since the PDB/PSF
                files have a 4 character limitation and require and exact match in the conf file.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                The name of the residues and their corresponding fugacity must specified
                for every residue in the system (i.e., {"residue_name" : fugacity}).
                Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY
                OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER.
                Example 1 (system with only water):  {"H2O" : 1} .
                Example 2 (system with water and ethanol):  {"H2O" : 0.5, "ETH" : 10},
            CBMC_First : int (>= 0), default = 12
                The number of CD-CBMC trials to choose the first atom position
                (Lennard-Jones trials for first seed growth).
            CBMC_Nth : int (>= 0), default = 10
                The Number of CD-CBMC trials to choose the later atom positions
                (Lennard-Jones trials for first seed growth).
            CBMC_Ang : int (>= 0), default = 50
                The Number of CD-CBMC bending angle trials to perform for geometry
                (per the coupled-decoupled CBMC scheme).
            CBMC_Dih : int (>= 0), default = 50
                The Number of CD-CBMC dihedral angle trials to perform for geometry
                (per the coupled-decoupled CBMC scheme).
            OutputName : str (NO SPACES), , default = "Output_data", default = [True, 1M] or
                [True , set via formula based on the number of RunSteps or 1M max]
                The UNIQUE STRING NAME, WITH NO SPACES, which is used for the
                output block average, PDB, and PSF file names.
            CoordinatesFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or M max]
                Controls output of PDB file (coordinates). If bool is True, this
                enables outputting the coordinate files at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting
                the coordinates.
            RestartFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or 1M max]
                This creates the PDB and PSF (coordinate and topology) files for
                restarting the system at the set steps_per_data_output_int (frequency)
                If bool is True, this enables outputting the PDB/PSF restart files at the
                integer frequency (set steps_per_data_output_int), while “false”
                disables outputting the PDB/PSF restart files.
            CheckpointFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or 1M max]
                Controls the output of the last state of simulation at a specified step,
                in a binary file format (checkpoint.dat). Checkpoint file contains the
                following information in full precision:
                    (1) Last simulation step that saved into checkpoint file
                    (2) Simulation cell dimensions and angles
                    (3) Maximum amount of displacement (Å), rotation (δ), and volume (Å^3)
                        that is used in the Displacement, Rotation, MultiParticle, and Volume moves
                    (4) Number of Monte Carlo move trial and acceptance
                    (5) All molecule’s coordinates
                    (6) Random number sequence
                If bool is True, this enables outputing the checkpoint file at the
                integer frequency (set steps_per_data_ouput_int),
                while "False" disables outputting the checkpoint file.'
            ConsoleFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the output to the "console” or log file, which prints the
                acceptance statistics, and run timing info. In addition, instantaneously-selected
                thermodynamic properties will be output to this file.  If bool is True,
                this enables outputting the console data at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting the console
                data file.
            BlockAverageFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the block averages output of selected thermodynamic properties.
                Block averages are averages of thermodynamic values of interest for chunks of the
                simulation (for post-processing of averages or std. dev. in those values).
                If bool is True, this enables outputting the block averaging data/file at the
                integer frequency (set steps_per_data_output_int),  while "False"
                disables outputting the block averaging data/file.
            HistogramFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the histograms. Histograms are a binned listing of observation frequency
                for a specific thermodynamic variable. In the GOMC code, they also control the output
                of a file containing energy/molecule samples, which is only used for the "GCMC"
                ensemble simulations for histogram reweighting purposes. If bool is True, this
                enables outputting the data to the histogram data at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting the histogram
                data.
            DistName : str (NO SPACES), default = "dis"
                Short phrase which will be combined with RunNumber and RunLetter
                to use in the name of the binned histogram for molecule distribution.
            HistName : str (NO SPACES), default = "his"
                Short phrase, which will be combined with RunNumber and RunLetter,
                to use in the name of the energy/molecule count sample file.
            RunNumber : int  ( > 0 ), default = 1
                 Sets a number, which is a part of DistName and HistName file name.
            RunLetter : str (1 alphabetic character only), default = "a"
                Sets a letter, which is a part of DistName and HistName file name.
            SampleFreq : int ( > 0 ), default = 500
                The number of steps per histogram sample or frequency.
            OutEnergy : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the energy data into the block averages and console output/log
            OutPressure : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the pressure data into the block averages and console output/log files.
            OutMolNumber : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the number of molecules data into the block averages and console
                output/log files.
            OutDensity : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the density data into the block averages and console output/log files.
            OutVolume : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the volume data into the block averages and console output/log files.
            OutSurfaceTension : [bool, bool], default = [False, False]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the surface tension data into the block averages and console
                output/log files.
            FreeEnergyCalc : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = None
                bool = True enabling free energy calculation during the simulation, false disables
                the calculation. The int/step frequency sets the frequency of calculating the free energy.
            MoleculeType : list [str , int (> 0)] or ["residue_name" , residue_ID], default = None
                The user must set this variable as there is no working default.
                Note: ONLY 4 characters can be used for the string (i.e., "residue_name").
                Sets the solute molecule kind (residue name) and molecule number (residue ID),
                which absolute solvation free will be calculated for.'
            InitialState : int (>= 0), default = None
                The user must set this variable as there is no working default.
                The index of LambdaCoulomb and LambdaVDW vectors. Sets the index of the
                LambdaCoulomb and LambdaVDW vectors, to determine the simulation lambda value for
                VDW and Coulomb interactions.
                WARNING : This must an integer within the vector count of the LambdaVDW and LambdaCoulomb,
                in which the counting starts at 0.  '
            LambdaVDW : list of floats (0 <= floats <= 1), default = None
                The user must set this variable as there is no working default (default = {}).
                Lambda values for VDW interaction in ascending order. Sets the intermediate
                lambda states to which solute-solvent VDW interactions are scaled.
                WARNING : This list must be the same length as the "LambdaCoulomb" list length.
                WARNING : All lambda values must be stated in the ascending order, otherwise the
                program will terminate.
                Example of ascending order 1: [0.1, 1.0,]
                Example of ascending orde 2: [0.1, 0.2, 0.4, 0.9]
            LambdaCoulomb : list of floats (0 <= floats <= 1), default = None
                Lambda values for Coulombic interaction in ascending order. Sets the intermediate
                lambda states to which solute-solvent Coulombic interactions are scaled.
                GOMC defauts to the "LambdaVDW" values for the Coulombic interaction
                if no "LambdaCoulomb" variable is set.
                WARNING : This list must be the same length as the "LambdaVDW" list length.
                WARNING : All lambda values must be stated in the ascending order, otherwise
                the program will terminate.
                Example of ascending order 1: [0.1, 1.0,]
                Example of ascending order 2: [0.1, 0.2, 0.4, 0.9] '
            ScaleCoulomb : bool, default = False
                Determines to scale the Coulombic interaction non-linearly
                (soft-core scheme) or not.
                True if the Coulombic interaction needs to be scaled non-linearly.
                False if the Coulombic interaction needs to be scaled linearly.
            ScalePower : int (>= 0), default = 2
                The p value in the soft-core scaling scheme, where the distance
                between solute and solvent is scaled non-linearly.
            ScaleAlpha : int or float (>= 0), default = 0.5
                The alpha value in the soft-core scaling scheme, where the distance
                between solute and solvent is scaled non-linearly.
            MinSigma : int or float (>= 0), default = 3
                The minimum sigma value in the soft-core scaling scheme, where the
                distance between solute and solvent is scaled non-linearly.
            DisFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.19, 'GCMC': 0.15}
                Fractional percentage at which the displacement move will occur
                (i.e., fraction of displacement moves).
            RotFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.15}
                Fractional percentage at which the rotation move will occur.
                (i.e., fraction of rotation moves).
            IntraSwapFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.3, 'NPT': 0.29, 'GEMC_NVT': 0.1, 'GEMC_NPT': 0.1, 'GCMC': 0.1}
                Fractional percentage at which the molecule will be removed from a
                box and inserted into the same box using coupled-decoupled configurational-bias
                algorithm. (i.e., fraction of intra-molecule swap moves).
            SwapFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.35}
                For Gibbs and Grand Canonical (GC) ensemble runs only: Fractional
                percentage at which molecule swap move will occur using coupled-decoupled
                configurational-bias. (i.e., fraction of molecule swaps moves).
            RegrowthFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.3, 'NPT': 0.3, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.15}
                Fractional percentage at which part of the molecule will be deleted and
                then regrown using coupled- decoupled configurational-bias algorithm
                (i.e., fraction of molecular growth moves).
            CrankShaftFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.1, 'NPT': 0.1, 'GEMC_NVT': 0.1, 'GEMC_NPT': 0.1, 'GCMC': 0.1}
                Fractional percentage at which crankshaft move will occur.
                In this move, two atoms that are forming angle or dihedral are selected
                randomly and form a shaft. Then any atoms or group that are within these
                two selected atoms, will rotate around the shaft to sample intra-molecular
                degree of freedom (i.e., fraction of crankshaft moves).
            VolFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.01, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.01, 'GCMC': 0.0}
                For isobaric-isothermal (NPT) ensemble and Gibbs ensemble
                (GEMC_NPT and GEMC_NVT) runs only: Fractional percentage at
                which a volume move will occur (i.e., fraction of Volume moves).
            MultiParticleFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which multi-particle move will occur.
                In this move, all molecules in the selected simulation box will be rigidly
                rotated or displaced simultaneously, along the calculated torque or force
                respectively (i.e., fraction of multi-particle moves).
            IntraMEMC_1Freq : int or float (0 <= value <= 1), default are specific for each ensemble
            {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume within
                same simulation box.  This move need additional information such as
                ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.
            MEMC_1Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will
                be exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes.  This move need additional information such as
                ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.
            IntraMEMC_2Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind
                will be exchanged with a specified large molecule kind in defined sub-volume
                within same simulation box. Backbone of small and large molecule kind will be
                used to insert the large molecule more efficiently. This move need additional
                information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            MEMC_2Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes. Backbone of small and large molecule kind will be
                used to insert the large molecule more efficiently. This move need additional
                information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            IntraMEMC_3Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume within same
                simulation box. Specified atom of the large molecule kind will be used to insert
                the large molecule using coupled-decoupled configurational-bias. This move need
                additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, and LargeKindBackBone. '
            MEMC_3Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes.  Specified atom of the large molecule kind will be
                used to insert the large molecule using coupled-decoupled configurational-bias.
                This move need additional information such as ExchangeVolumeDim,
                ExchangeRatio, ExchangeSmallKind, ExchangeLargeKind, and LargeKindBackBone.
            ExchangeVolumeDim : list of 3 floats or integers or [X-dimension, Y-dimension, Z-dimension)],
                default = [1.0, 1.0, 1.0]
                To use all variations of MEMC and Intra-MEMC Monte Carlo moves, the exchange
                subvolume must be defined. The exchange sub-volume is defined as an orthogonal box
                with x, y, and z-dimensions, where small molecule/molecules kind will be selected
                from to be exchanged with a large molecule kind.
                Note: Currently, the X and Y dimension cannot be set independently (X = Y = max(X, Y)).
                Note: A heuristic for setting good values of the x, y, and z-dimensions is to use
                the geometric size of the large molecule plus 1-2 Å in each dimension.
                Note: In case of exchanging 1 small molecule kind with 1 large molecule kind in
                IntraMEMC-2, IntraMEMC-3, MEMC-2, MEMC-3 Monte Carlo moves, the sub-volume
                dimension has no effect on acceptance rate. '
            MEMC_DataInput : nested lists, default = None
                Enter data as a list with some sub-lists as follows:
                [[ExchangeRatio_int (> 0), ExchangeLargeKind_str,
                [LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ],
                ExchangeSmallKind_str, [SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]],
                ...,
                [ExchangeRatio_int (> 0), ExchangeLargeKind_str,
                [LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ],
                ExchangeSmallKind_str, [SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ].
                NOTE: CURRENTLY ALL THESE INPUTS NEED TO BE SPECIFIED, REGARDLESS OF THE MEMC TYPE
                SELECTION. IF THE SmallKindBackBone or LargeKindBackBone IS NOT REQUIRED FOR THE
                MEMC TYPE, None CAN BE USED IN PLACE OF A STRING.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                It is recommended that the user print the Charmm object psf and pdb files
                and review the residue names that match the atom name before using the in
                the MEMC_DataInput variable input.
                Note: see the below data explanations for the ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, LargeKindBackBone, SmallKindBackBone.
                Example 1 (MEMC-1) : [ [1, 'WAT', [None, None], 'wat', [None, None]] ,
                [1, 'WAT', [None, None], 'wat', [None, None]] .
                Example 2 (MEMC-2): [ [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1' ]] ,
                [1, 'WAT', ['H1', 'H2'], 'wat', ['H1', 'H2' ]] .
                Example 3 (MEMC-3) : [ [2, 'WAT', 'O1', 'H1'], 'wat', [None, None]] ,
                [2, 'WAT', ['H1', 'H2'], 'wat', [None, None]] .\n"
                --- ExchangeRatio     = MEMC parameters (all ensembles): int (> 0), default = None
                                        The Ratio of exchanging small molecule/molecules with 1 large molecule.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the exchange ratio must be defined. The exchange ratio defines how
                                        many small molecule will be exchanged with 1 large molecule. For each
                                        large-small molecule pairs, one exchange ratio must be defined.
                --- ExchangeSmallKind = MEMC parameters (all ensembles):  str, default = None
                                        The small molecule kind (resname) to be exchanged.
                                        Note: ONLY 4 characters can be used for the strings.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the small molecule kind to be exchanged with a large molecule
                                        'kind must be defined. Multiple small molecule kind can be specified.
                --- ExchangeLargeKind = MEMC parameters (all ensembles):  str, default = None
                                        The large molecule kind (resname) to be exchanged.
                                        Note: ONLY 4 characters can be used for the strings.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the large molecule kind to be exchanged with small molecule '
                                        kind must be defined. Multiple large molecule kind can be specified.
                --- LargeKindBackBone = MEMC parameters (all ensembles): list [str, str] or [None, None], default = None
                                        Note: ONLY 4 characters can be used for the strings.
                                        The [None, None] values can only be used if that MEMC type does
                                        not require them.  The strings for the the atom name 1 and atom name 2
                                        that belong to the large molecule’s backbone
                                        (i.e., [str_for_atom_name_1, str_for_atom_name_2])
                                        To use MEMC-2, MEMC-3, IntraMEMC-2, and IntraMEMC-3 Monte Carlo moves, the
                                        large molecule backbone must be defined. The backbone of the molecule is defined
                                        as a vector that connects two atoms belong to the large molecule. The large
                                        molecule backbone will be used to align the sub-volume in MEMC-2 and IntraMEMC-2
                                        moves, while in MEMC-3 and IntraMEMC-3 moves, it uses the atom name to start
                                        growing the large molecule using coupled-decoupled configurational-bias. For
                                        each large-small molecule pairs, two atom names must be defined.
                                        Note: all atom names in the molecule must be unique.
                                        Note: In MEMC-3 and IntraMEMC-3 Monte Carlo moves, both atom names must be same,
                                        otherwise program will be terminated.
                                        Note: If the large molecule has only one atom (mono atomic molecules),
                                        same atom name must be used for str_for_atom_name_1 and str_for_atom_name_2
                                        of the LargeKindBackBone.
                --- SmallKindBackBone = MEMC parameters (all ensembles): list [str, str] or [None, None], default = None
                                        Note: ONLY 4 characters can be used for the strings.
                                        The [None, None] values can only be used if that MEMC type does not
                                        require them. The strings for the the atom name 1 and atom name 2 that
                                        belong to the small molecule’s backbone
                                        (i.e., [str_for_atom_name_1, str_for_atom_name_2]) '
                                        To use MEMC-2, and IntraMEMC-2 Monte Carlo moves, the small molecule backbone
                                        must be defined. The backbone of the molecule is defined as a vector that
                                        connects two atoms belong to the small molecule and will be used to align the
                                        sub-volume. For each large-small molecule pairs, two atom names must be defined.
                                        Note: all atom names in the molecule must be unique.
                                        Note: If the small molecule has only one atom (mono atomic molecules), same atom
                                        name must be used str_for_atom_name_1 and str_for_atom_name_2
                                        of the SmallKindBackBone.
            # *******************************************************************
            # input_variables_dict options (keys and values) - (end)
            # Note: the input_variables_dict keys are also attributes
            # *******************************************************************

        Attributes
        ----------
        input_error : bool
            This error is typically incurred from an error in the user's input values.
            However, it could also be due to a bug, provided the user is inputting
            the data as this Class intends.
        all_failed_input_List
        ensemble_typ : str, ['NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', 'GCMC']
            The ensemble type of the simulation.
        RunSteps : int (>0), must be an integer greater than zero.
            Sets the total number of simulation steps.
        Temperature : float or int (>0), must be an integer greater than zero.
            Temperature of system in Kelvin (K)
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
        conf_filename : str
            The name of the GOMC contol file, which will be created.  The extension
            of the GOMC control file can be .conf, or no extension can be provided.
            If no extension is provided, this writer will automatically add the
            .conf extension to the provided string.
        Coordinates_box_0 : str
            The coordinate or PDB file for box 0 in the simulation.
        Coordinates_box_1 : str or None
            The coordinate or PDB file for box 1 in the simulation.  This is only for
            GCMC, GEMC_NVT, and GEMC_NVT simulations. If running a NVT or NPT
            simulation, the value will be None.
        Structures_box_0 : str
            The structure file or PSF file for box 0 in the simulation.
            The coordinate or PDB file for box 1 in the simulation.  This is only for
            GCMC, GEMC_NVT, and GEMC_NVT simulations. If running a NVT or NPT
            simulation, the value will be None.
        Structures_box_1 : str or None
        The structure file or PSF file for box 1 in the simulation.  This is only for
            GCMC, GEMC_NVT, and GEMC_NVT simulations. If running a NVT or NPT
            simulation, the value will be None.
        x_dim_box_0 : float or int
            The x-dimension of box 0.  Currently, only orthogonal boxes are supported.
        y_dim_box_0 : float or int
            The y-dimension of box 0.  Currently, only orthogonal boxes are supported.
        z_dim_box_0 : float or int
            The z-dimension of box 0.  Currently, only orthogonal boxes are supported.
        x_dim_box_1 : float or int
            The x-dimension of box 1.  Currently, only orthogonal boxes are supported.
        y_dim_box_1 : float or int
            The y-dimension of box 1.  Currently, only orthogonal boxes are supported.
        z_dim_box_1 : float or int
            The z-dimension of box 1.  Currently, only orthogonal boxes are supported.
        coul_1_4 : float or int
            The non-bonded 1-4 coulombic scaling factor, which is the
            same for all the residues/molecules, regardless if
            differenct force fields are utilized.
        residues : list, [str, ..., str]
            Labels of unique residues in the Compound. Residues are assigned by
            checking against Compound.name.  Only supply residue names as 4 character
            strings, as the residue names are truncated to 4 characters to fit in the
            psf and pdb file.
        all_res_unique_atom_name_dict : dict, {str : [str, ..., str]}
            A dictionary that provides the residue names (keys) and a list
            of the unique atom names in the residue (value), for the
            combined structures (box 0 and box 1 (if supplied)).
        any input_variables_dict key : varies (see each input_variables_dict key and value)
            Any of the input variables keys is also an Attribute and can be called
            the same way.  Please see the input_variables_dict keys in the
            Parameters section above for all the available attributes.








        Notes
        -------
        The attribute's default values and the specific ensembles they are
        also available with can be accessed by the running
        print_valid_ensemble_input_variables('NPT', description = True)
        command, as the information is dynamically contained here.

        The details of the required inputs for the selected
        ensembles can be found by the following function,
        >>> print_valid_required_input_variables('NVT', description = True)
        which prints the required inputs with their subsection description
        for the selected 'NVT' ensemble (other ensembles can be set as well).
        The box units imported are in nm (standard MoSDeF units).
        The units for this writer are auto-scaled to Angstroms, so they
        can be directly used in the GOMC or NAMD engines.

        Note: all of the move types are not available in for every ensemble.
        Note: all of the move fractions must sum to 1, or the control file
        writer will fail.

        The attribute variables and text extracted with permission from the GOMC
        manual version 2.60. Some of the text was modified from its original version.
        Cite: Potoff, Jeffrey; Schwiebert, Loren; et. al. GOMC Documentation.
        https://raw.githubusercontent.com/GOMC-WSU/GOMC/master/GOMC_Manual.pdf, 2021.
        """

        # set this to check and see if all the input pass
        self.input_error = False

        # set this to check and see if all the input pass
        self.all_failed_input_List = []

        # Check if charmm_object is really a Charmm() object
        if not isinstance(charmm_object, mf_charmm.Charmm):
            self.input_error = True
            print_error_message = 'The variable supplied as a charmm_object ({}}) is not a ' \
                                  'charmm_object ({}})'.format(type(mf_charmm.Charmm), type(mf_charmm.Charmm))
            raise TypeError(print_error_message)

        # check ensemble is a correct type
        print('INFO: ensemble_type = ' + str(ensemble_type))
        if ensemble_type in ['NPT', 'NVT', 'GCMC', 'GEMC_NVT', 'GEMC_NPT']:
            self.ensemble_type = ensemble_type
            print("INFO: All the ensemble (ensemble_type) input passed the initial error checking")
        else:
            self.input_error = True
            print_error_message =  "ERROR: The ensemble type selection of {}  is not a valid ensemble option. " \
                                  "Please choose the 'NPT', 'NVT', 'GEMC_NVT','GEMC_NPT', or 'GCMC' " \
                                  "ensembles".format(ensemble_type)
            raise ValueError(print_error_message)

        self.RunSteps = RunSteps
        self.Temperature = Temperature
        if charmm_object.ff_filename is not None and isinstance(charmm_object.ff_filename, str) is True:
            self.ff_filename = charmm_object.ff_filename
        elif charmm_object.ff_filename is None or isinstance(charmm_object.ff_filename, str) is False:
            self.input_error = True
            print_error_message = "The force field file name was not specified and in the Charmm object ({}})." \
                                  "Therefore, the force field file (.inp) can not be written, and thus, the " \
                                  "GOMC control file (.conf) can not be created. Please use the force field file " \
                                  "name when building the Charmm object ({}})".format(type(mf_charmm.Charmm),
                                                                                      type(mf_charmm.Charmm))
            raise ValueError(print_error_message)

        if charmm_object.filename_box_0 is not None and isinstance(charmm_object.filename_box_0, str) is True:
            self.Coordinates_box_0 = str(charmm_object.filename_box_0) + '.pdb'
            self.Structures_box_0 = str(charmm_object.filename_box_0) + '.psf'
        if charmm_object.filename_box_1 is not None and isinstance(charmm_object.filename_box_1, str) is True:
            self.Coordinates_box_1 = str(charmm_object.filename_box_1) + '.pdb'
            self.Structures_box_1 = str(charmm_object.filename_box_1) + '.psf'
        else:
            self.Coordinates_box_1 = None
            self.Structures_box_1 = None

        self.coul_1_4 = charmm_object.coul_1_4
        self.input_variables_dict = input_variables_dict
        self.residues = charmm_object.residues
        self.all_residues_unique_atom_name_dict = charmm_object.all_res_unique_atom_name_dict

        self.x_dim_box_0 = charmm_object.box_0.maxs[0] * 10   # times 10 to convert from nm to Angstroms
        self.y_dim_box_0 = charmm_object.box_0.maxs[1] * 10   # times 10 to convert from nm to Angstroms
        self.z_dim_box_0 = charmm_object.box_0.maxs[2] * 10   # times 10 to convert from nm to Angstroms

        print('charmm_object.filename_box_1 = ',format(str(charmm_object.filename_box_1)))
        print('type(charmm_object.filename_box_1) = ', format(str(type(charmm_object.filename_box_1))))

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            if charmm_object.filename_box_1 is not None and isinstance(charmm_object.filename_box_1, str) is True:
                self.x_dim_box_1 = charmm_object.box_1.maxs[0] * 10   # times 10 to convert from nm to Angstroms
                self.y_dim_box_1 = charmm_object.box_1.maxs[1] * 10   # times 10 to convert from nm to Angstroms
                self.z_dim_box_1 = charmm_object.box_1.maxs[2] * 10   # times 10 to convert from nm to Angstroms
            else:
                self.x_dim_box_1 = None
                self.y_dim_box_1 = None
                self.z_dim_box_1 = None

        # check if the ensembles have the correct number of boxes in the charmm object
        if self.ensemble_type in ['NVT', 'NPT'] and \
                self.Coordinates_box_1 is not None and self.Structures_box_1 is not None:
            self.input_error = True
            print_error_message = "ERROR: The ensemble type selection of {} is using a Charmm " \
                                  "object with two simulation boxes, and the {} ensemble only accepts " \
                                  "one box (box 0)." \
                                  "".format(ensemble_type, ensemble_type)
            raise ValueError(print_error_message)

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC'] and \
                self.Coordinates_box_1 is None and self.Structures_box_1 is None:
            self.input_error = True
            print_error_message = "ERROR: The ensemble type selection of {} is using a Charmm " \
                                  "object with one simulation boxes, and the {} ensemble only accepts " \
                                  "two boxes (box 0 and box 1)." \
                                  "".format(ensemble_type, ensemble_type)
            raise ValueError(print_error_message)


        # check the box dimensions
        ck_box_dim_is_float_or_int_greater_0(self.x_dim_box_0, 'x', 0, self.ensemble_type)
        ck_box_dim_is_float_or_int_greater_0(self.y_dim_box_0, 'y', 0, self.ensemble_type)
        ck_box_dim_is_float_or_int_greater_0(self.z_dim_box_0, 'z', 0, self.ensemble_type)

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            ck_box_dim_is_float_or_int_greater_0(self.x_dim_box_1, 'x', 1, self.ensemble_type)
            ck_box_dim_is_float_or_int_greater_0(self.y_dim_box_1, 'y', 1, self.ensemble_type)
            ck_box_dim_is_float_or_int_greater_0(self.z_dim_box_1, 'z', 1, self.ensemble_type)

        # the future control file name is entered now as None
        self.conf_filename = None

        # list of bad variable inputs
        bad_input_variables_values_list = []

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
        self.VDWGeometricSigma = default_input_variables_dict['VDWGeometricSigma']
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

        if input_variables_dict is None:
            self.input_variables_dict = {}
        elif isinstance(input_variables_dict, dict) is True:
            self.input_variables_dict = input_variables_dict
        else:
            self.input_error = True
            print_error_message = "ERROR: The input_variables_dict variable is not None or a dictionary. "
            raise ValueError(print_error_message)

        # Create all lower case spelled keywords, and return case specific keywords
        # Also, creates a dict to convert the lower case keys to case sensitive keys
        all_input_var_case_spec_list = _get_all_possible_input_variables(description=False)
        all_input_var_case_unspec_list = []
        all_input_var_case_unspec_to_spec_dict = {}
        for var_i in all_input_var_case_spec_list:
            all_input_var_case_unspec_list.append(var_i.lower())
            all_input_var_case_unspec_to_spec_dict.update({var_i.lower(): var_i})

        # create/fix user case insensitive input variables (input_variables_dict) keys to case sensitive keys
        input_var_dict_orig_keys_list = dict_keys_to_list(self.input_variables_dict)
        for z_j in range(0, len(input_var_dict_orig_keys_list)):
            key_lower = input_var_dict_orig_keys_list[z_j].lower()
            if key_lower in all_input_var_case_unspec_list:
                input_variables_dict[all_input_var_case_unspec_to_spec_dict[key_lower]] = \
                    input_variables_dict.pop(input_var_dict_orig_keys_list[z_j])

        # check that the coulombic 1-4 scalar is : 0 =< 1-4 scalar <=1
        if (isinstance(self.coul_1_4, int) is False
                and isinstance(self.coul_1_4, float) is False) \
                or self.coul_1_4 < 0 or self.coul_1_4 > 1:
            self.input_error = True
            print_error_message = "ERROR: The selected 1-4 Coulombic scalar ({{) is not correct. "\
                                  "The 1-4 Coulombic scalar need to be an integer or float from " \
                                  "0 to 1.".format(self.coul_1_4)
            raise ValueError(print_error_message)

        # check that the Temperature is valid
        if self.Temperature <= 1:
            self.input_error = True
            print_error_message = "ERROR: The selected Temperature ({}) is equal to or less than 1 Kelvin. "\
                                  "Please select a valid Temperature".format(self.Temperature)
            raise ValueError(print_error_message)
        else:
            print("INFO: All the Temperature  (Temperature) input passed the initial error checking")

        # RunSteps
        if not isinstance(self.RunSteps, int) or self.RunSteps <= 0:
            self.input_error = True
            print_error_message = "ERROR: The selected run steps (RunSteps variable = {}) is not "\
                                  "an integer or is less than or equal to 0.".format(self.RunSteps)
            raise ValueError(print_error_message)

        # create a list of the possible required files and check them based on the ensemble
        required_data_list = [self.ff_filename, self.Coordinates_box_0,  self.Structures_box_0]

        if self.Coordinates_box_1 is not None:
            required_data_list.append(self.Coordinates_box_1)
        if self.Structures_box_1 is not None:
            required_data_list.append(self.Structures_box_1)

        if self.ensemble_type in ['NVT', 'NPT']:
            if len(required_data_list) != 3 \
                    or os.path.splitext(self.ff_filename)[1] != '.inp' \
                    or os.path.splitext(self.Coordinates_box_0)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_0)[1] != '.psf':
                self.input_error = True
                print_error_message = 'ERROR: The proper force field, PDB, and psf files were not provided, '\
                                      'or at least their extentions are not correct '\
                                      '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were '\
                                      'provided for the NVT or NPT simulations, which is not allowed'
                raise ValueError(print_error_message)

        else:
            print("INFO: All the required force field, pdb, and psf files for box 0 (.inp, .pdb, and .psf) all " 
                  "passed the intial error checking. Note: the file names and their existance is not confirmed.")

        if self.ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']:
            if len(required_data_list) != 5 \
                    or os.path.splitext(self.ff_filename)[1] != '.inp' \
                    or os.path.splitext(self.Coordinates_box_0)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_0)[1] != '.psf' \
                    or os.path.splitext(self.Coordinates_box_1)[1] != '.pdb' \
                    or os.path.splitext(self.Structures_box_1)[1] != '.psf':
                warn('ERROR: The proper force field, PDB, and psf files were not provided, '
                     'or at least their extentions are not correct '
                     '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were not provided '
                     'for the GEMC_NVT, GEMC_NPT or GCMC simulations, which is not allowed')
                self.input_error = True
                print_error_message = 'ERROR: The proper force field, PDB, and psf files were not provided, '\
                                      'or at least their extentions are not correct '\
                                      '(i.e., not .inp, .pdb, or .psf). Or box 1 PSF and PDB files were not provided '\
                                      'for the GEMC_NVT, GEMC_NPT or GCMC simulations, which is not allowed'
                raise ValueError(print_error_message)
        else:
            print("INFO: All the required force field, pdb, and psf files for box 0 and 1 (.inp, .pdb, and .psf) all "
                  "passed the intial error checking. Note: the file names and their existance is not confirmed.")


        # verify all input variables keys are valid
        input_variables_dict_keys_list = dict_keys_to_list(self.input_variables_dict)
        if check_valid_ensemble_input_variables(self.ensemble_type, input_variables_dict_keys_list)[0] is False:
            returned_ck_bad_inputs_list = check_valid_ensemble_input_variables(self.ensemble_type,
                                                                               input_variables_dict_keys_list)[1]
            self.input_error = True
            print_error_message = "ERROR: All the correct input variables where not provided for the {} " \
                                  "ensemble. Please be sure to check that the keys in the " \
                                  "input variables dictionary (input_variables_dict) is correct, and be aware "\
                                  "that added spaces before or after the variable in any keys " \
                                  "will also give this warning. The bad variable inputs " \
                                  "ensemble inputs = {}".format(self.ensemble_type, returned_ck_bad_inputs_list)
            raise ValueError(print_error_message)

        # verify all input variable values are valid, for their keys
        input_var_keys_list = dict_keys_to_list(self.input_variables_dict)

        possible_ensemble_variables_list = _get_possible_ensemble_input_variables(self.ensemble_type)

        # check to make sure the VDW FF (ParaTypeCHARMM) is set true  for multiple ones by the user
        # (i.e., ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI)
        vdw_ck_list = []
        for vdw_ck_inter in range(0, len(input_var_keys_list)):
            if input_var_keys_list[vdw_ck_inter] in ['ParaTypeCHARMM', 'ParaTypeMie', 'ParaTypeMARTINI'] \
                    and self.input_variables_dict[input_var_keys_list[vdw_ck_inter]] is True:
                vdw_ck_list.append(True)

        if sum(vdw_ck_list) > 1:
            self.input_error = True
            print_error_message = 'ERROR: there can only be 1 VDW set to true.  Please set only one of the '\
                                  'ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI types to True in the ' \
                                  'user variable input'
            raise ValueError(print_error_message)

        # check for MC move ratios and zero all of them if any are in the input variables
        for var_iter in range(0, len(input_var_keys_list)):
            # standard MC moves
            key_move_list = ["DisFreq", "RotFreq", "IntraSwapFreq", "SwapFreq", "RegrowthFreq",
                             "CrankShaftFreq", "VolFreq", "MultiParticleFreq",
                             "IntraMEMC-1Freq", "MEMC-1Freq", "IntraMEMC-2Freq", "MEMC-2Freq",
                             "IntraMEMC-3Freq", "MEMC-3Freq"]

            if input_var_keys_list[var_iter] in key_move_list:
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

            if input_var_keys_list[var_iter] == "Potential":
                self.Potential = self.input_variables_dict["Potential"]

        # check for bad input variables and list the bad ones
        for var_iter in range(0, len(input_var_keys_list)):
            key = "Restart"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Restart = self.input_variables_dict[key]

            key = "RestartCheckpoint"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RestartCheckpoint = self.input_variables_dict[key]

            key = "PRNG"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != "RANDOM" \
                        and isinstance(self.input_variables_dict[key], int) is not True:
                    bad_input_variables_values_list.append(key)
                if isinstance(self.input_variables_dict[key], int) is True:
                    if self.input_variables_dict[key] < 0:
                        bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.PRNG = self.input_variables_dict[key]

            key = "ParaTypeCHARMM"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ParaTypeCHARMM = self.input_variables_dict[key]
                    if self.input_variables_dict[key] is True:
                        self.ParaTypeMie = False
                        self.ParaTypeMARTINI = False

            key = "ParaTypeMie"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ParaTypeMie = self.input_variables_dict[key]
                    if self.input_variables_dict[key] is True:
                        self.ParaTypeCHARMM = False

            key = "ParaTypeMARTINI"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ParaTypeMARTINI = self.input_variables_dict[key]
                    if self.input_variables_dict[key] is True:
                        self.ParaTypeCHARMM = False

            key = "RcutCoulomb_box_0"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RcutCoulomb_box_0 = self.input_variables_dict[key]

            key = "RcutCoulomb_box_1"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RcutCoulomb_box_1 = self.input_variables_dict[key]

            key = "Pressure"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Pressure = self.input_variables_dict[key]

            key = "Rcut"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Rcut = self.input_variables_dict[key]

            key = "RcutLow"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)
                if isinstance(self.input_variables_dict[key], float) \
                        or isinstance(self.input_variables_dict[key], int):
                    if self.input_variables_dict[key] > self.Rcut:
                        bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RcutLow = self.input_variables_dict[key]

            key = "LRC"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.LRC = self.input_variables_dict[key]

            key = "Exclude"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != '1-2' and self.input_variables_dict[
                    key] != '1-3' \
                        and self.input_variables_dict[key] != '1-4':
                    bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Exclude = self.input_variables_dict[key]

            key = "Potential"
            if input_var_keys_list[var_iter] == key:
                if self.input_variables_dict[key] != 'VDW' and self.input_variables_dict[
                    key] != 'EXP6' \
                        and self.input_variables_dict[key] != 'SHIFT' and \
                        self.input_variables_dict[key] != 'SWITCH':
                    bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Potential = self.input_variables_dict[key]

            key = "Rswitch"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if (isinstance(self.input_variables_dict[key], float)
                        or isinstance(self.input_variables_dict[key], int)) \
                        and (isinstance(self.RcutLow, float)
                             or isinstance(self.RcutLow, int)) \
                        and (isinstance(self.Rcut, float)
                             or isinstance(self.Rcut, int)) \
                        and self.Potential == "SWITCH":
                    if self.input_variables_dict[key] <= self.RcutLow \
                            or self.input_variables_dict[key] >= self.Rcut:
                        bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Rswitch = self.input_variables_dict[key]

            key = "ElectroStatic"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ElectroStatic = self.input_variables_dict[key]

            key = "Ewald"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Ewald = self.input_variables_dict[key]

            key = "CachedFourier"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CachedFourier = self.input_variables_dict[key]

            key = "Tolerance"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_float_greater_zero_less_1(self.input_variables_dict,
                                                                 key,
                                                                 bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Tolerance = self.input_variables_dict[key]

            key = "Dielectric"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Dielectric = self.input_variables_dict[key]

            key = "PressureCalc"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_zero_or_greater(self.input_variables_dict,
                                                                     key,
                                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.PressureCalc = self.input_variables_dict[key]

            key = "EqSteps"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.EqSteps = self.input_variables_dict[key]

            key = "AdjSteps"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.AdjSteps = self.input_variables_dict[key]

            key = "VDWGeometricSigma"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.VDWGeometricSigma = self.input_variables_dict[key]

            key = "useConstantArea"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.useConstantArea = self.input_variables_dict[key]

            key = "FixVolBox0"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.FixVolBox0 = self.input_variables_dict[key]

            # ChemPot and Fugacity are only for GCMC
            key = "ChemPot"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_GCMC_dict_str_int_or_float(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ChemPot = self.input_variables_dict[key]

            key = "Fugacity"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_GCMC_dict_str_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                                  key,
                                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.Fugacity = self.input_variables_dict[key]

            key = "CBMC_First"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CBMC_First = self.input_variables_dict[key]

            key = "CBMC_Nth"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CBMC_Nth = self.input_variables_dict[key]

            key = "CBMC_Ang"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CBMC_Ang = self.input_variables_dict[key]

            key = "CBMC_Dih"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CBMC_Dih = self.input_variables_dict[key]

            key = "OutputName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                          key,
                                                          bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutputName = self.input_variables_dict[key]

            key = "CoordinatesFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CoordinatesFreq = self.input_variables_dict[key]

            key = "RestartFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RestartFreq = self.input_variables_dict[key]

            key = "CheckpointFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CheckpointFreq = self.input_variables_dict[key]

            key = "ConsoleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ConsoleFreq = self.input_variables_dict[key]

            key = "BlockAverageFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.BlockAverageFreq = self.input_variables_dict[key]

            key = "HistogramFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.HistogramFreq = self.input_variables_dict[key]

            key = "DistName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                          key,
                                                          bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.DistName = self.input_variables_dict[key]

            key = "HistName"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_str_with_no_spaces(self.input_variables_dict,
                                                          key,
                                                          bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.HistName = self.input_variables_dict[key]

            key = "RunNumber"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RunNumber = self.input_variables_dict[key]

            key = "RunLetter"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], str) is not True:
                    bad_input_variables_values_list.append(key)
                if isinstance(self.input_variables_dict[key], str) is True:
                    if len(self.input_variables_dict[key]) != 1:
                        bad_input_variables_values_list.append(key)
                    elif len(self.input_variables_dict[key]) == 1:
                        is_run_letter_alphabet_char = self.input_variables_dict[key].isalpha()
                        if is_run_letter_alphabet_char is False:
                            bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RunLetter = self.input_variables_dict[key]

            key = "SampleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_greater_zero(self.input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.SampleFreq = self.input_variables_dict[key]

            key = "OutEnergy"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutEnergy = self.input_variables_dict[key]

            key = "OutPressure"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutPressure = self.input_variables_dict[key]

            key = "OutMolNumber"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)
                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutMolNumber = self.input_variables_dict[key]

            key = "OutDensity"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutDensity = self.input_variables_dict[key]

            key = "OutVolume"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutVolume = self.input_variables_dict[key]

            key = "OutSurfaceTension"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_bool(self.input_variables_dict,
                                                      key,
                                                      bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.OutSurfaceTension = self.input_variables_dict[key]

            key = "FreeEnergyCalc"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_bool_int_greater_zero(self.input_variables_dict,
                                                                  key,
                                                                  bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.FreeEnergyCalc = self.input_variables_dict[key]

            key = "MoleculeType"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_residue_str_int_greater_zero(self.input_variables_dict,
                                                                         key,
                                                                         bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MoleculeType = self.input_variables_dict[key]

            key = "InitialState"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.InitialState = self.input_variables_dict[key]

            key = "LambdaVDW"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_of_floats_zero_to_1(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.LambdaVDW = self.input_variables_dict[key]

            key = "LambdaCoulomb"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_list_of_floats_zero_to_1(self.input_variables_dict,
                                                                key,
                                                                bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.LambdaCoulomb = self.input_variables_dict[key]

            key = "ScaleCoulomb"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_true_or_false(self.input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ScaleCoulomb = self.input_variables_dict[key]

            key = "ScalePower"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_zero_or_greater(self.input_variables_dict,
                                                           key,
                                                           bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ScalePower = self.input_variables_dict[key]

            key = "ScaleAlpha"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ScaleAlpha = self.input_variables_dict[key]

            key = "MinSigma"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_or_greater(self.input_variables_dict,
                                                                    key,
                                                                    bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MinSigma = self.input_variables_dict[key]

            # standard MC moves
            key = "DisFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.DisFreq = self.input_variables_dict[key]

            key = "RotFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RotFreq = self.input_variables_dict[key]

            key = "IntraSwapFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.IntraSwapFreq = self.input_variables_dict[key]

            key = "SwapFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.SwapFreq = self.input_variables_dict[key]
                else:
                    self.VolFreq = 0.00

            key = "RegrowthFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.RegrowthFreq = self.input_variables_dict[key]

            key = "CrankShaftFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.CrankShaftFreq = self.input_variables_dict[key]

            key = "VolFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.VolFreq = self.input_variables_dict[key]
                else:
                    self.VolFreq = 0.00

            key = "MultiParticleFreq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MultiParticleFreq = self.input_variables_dict[key]

            # MEMC moves freqencies
            key = "IntraMEMC-1Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.IntraMEMC_1Freq = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_1Freq = 0.00

            key = "MEMC-1Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MEMC_1Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_1Freq = 0.00

            key = "IntraMEMC-2Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.IntraMEMC_2Freq = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_2Freq = 0.00

            key = "MEMC-2Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MEMC_2Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_2Freq = 0.00

            key = "IntraMEMC-3Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.IntraMEMC_3Freq = self.input_variables_dict[key]
                else:
                    self.IntraMEMC_3Freq = 0.00

            key = "MEMC-3Freq"
            if input_var_keys_list[var_iter] == key:
                self.ck_input_variable_int_or_float_zero_to_1(self.input_variables_dict,
                                                              key,
                                                              bad_input_variables_values_list)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MEMC_3Freq = self.input_variables_dict[key]
                else:
                    self.MEMC_3Freq = 0.00

            key = "ExchangeVolumeDim"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], list) is False:
                    bad_input_variables_values_list.append(key)
                elif isinstance(self.input_variables_dict[key], list) is True:
                    if len(self.input_variables_dict[key]) != 3 \
                            or (isinstance(self.input_variables_dict[key][0], float) is not True
                                and isinstance(self.input_variables_dict[key][0], int) is not True) \
                            or (isinstance(self.input_variables_dict[key][1], float) is not True
                                and isinstance(self.input_variables_dict[key][1], int) is not True) \
                            or (isinstance(self.input_variables_dict[key][2], float) is not True
                                and isinstance(self.input_variables_dict[key][2], int) is not True) \
                            or str(self.input_variables_dict[key][0]) == str(True) \
                            or str(self.input_variables_dict[key][0]) == str(False) \
                            or str(self.input_variables_dict[key][1]) == str(True) \
                            or str(self.input_variables_dict[key][1]) == str(False) \
                            or str(self.input_variables_dict[key][2]) == str(True) \
                            or str(self.input_variables_dict[key][2]) == str(False):
                        bad_input_variables_values_list.append(key)
                    elif len(self.input_variables_dict[key]) == 3:
                        if (isinstance(self.input_variables_dict[key][0], float) is True
                            or isinstance(self.input_variables_dict[key][0], int) is True) \
                                and (isinstance(self.input_variables_dict[key][1], float) is True
                                     or isinstance(self.input_variables_dict[key][1], int) is True) \
                                and (isinstance(self.input_variables_dict[key][2], float) is True
                                     or isinstance(self.input_variables_dict[key][2], int) is True):
                            if self.input_variables_dict[key][0] <= 0 \
                                or self.input_variables_dict[key][1] <= 0  \
                                    or self.input_variables_dict[key][2] <= 0:
                                bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.ExchangeVolumeDim = self.input_variables_dict[key]

            key = "MEMC_DataInput"
            if input_var_keys_list[var_iter] == key:
                if isinstance(self.input_variables_dict[key], list) is False:
                    bad_input_variables_values_list.append(key)
                elif isinstance(self.input_variables_dict[key], list) is True:
                    if len(self.input_variables_dict[key]) == 0:
                        bad_input_variables_values_list.append(key)
                    elif len(self.input_variables_dict[key]) > 0:
                        no_memc_combos = len(self.input_variables_dict[key])
                        for MEMC_iter in range(0, no_memc_combos):
                            if isinstance(self.input_variables_dict[key][MEMC_iter], list) is False:
                                bad_input_variables_values_list.append(key)
                            elif isinstance(self.input_variables_dict[key][MEMC_iter], list) is True \
                                    and len(self.input_variables_dict[key][MEMC_iter]) == 5:
                                if isinstance(self.input_variables_dict[key][MEMC_iter][2], list) is False \
                                        or isinstance(self.input_variables_dict[key][MEMC_iter][4], list) is False:
                                    bad_input_variables_values_list.append(key)
                                elif isinstance(self.input_variables_dict[key][MEMC_iter][2], list) is True \
                                        and isinstance(self.input_variables_dict[key][MEMC_iter][4], list) is True:
                                    if (len(self.input_variables_dict[key][MEMC_iter][2]) != 2
                                            or len(self.input_variables_dict[key][MEMC_iter][4]) != 2):
                                        bad_input_variables_values_list.append(key)
                                    elif (len(self.input_variables_dict[key][MEMC_iter][2]) == 2
                                            and len(self.input_variables_dict[key][MEMC_iter][4]) == 2):
                                        if isinstance(self.input_variables_dict[key][MEMC_iter][0], int) is not True \
                                                or str(self.input_variables_dict[key][MEMC_iter][0]) == str(True) \
                                                or str(self.input_variables_dict[key][MEMC_iter][0]) == str(False) \
                                                or isinstance(self.input_variables_dict[key][MEMC_iter][1], str
                                                              ) is False \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][2][0], str
                                                               ) is False
                                                    and self.input_variables_dict[key][MEMC_iter][2][0] is not None) \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][2][1], str
                                                               ) is False
                                                    and self.input_variables_dict[key][MEMC_iter][2][1] is not None) \
                                                or isinstance(self.input_variables_dict[key][MEMC_iter][3], str
                                                              ) is False \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][4][0], str
                                                               ) is False
                                                    and self.input_variables_dict[key][MEMC_iter][4][0] is not None) \
                                                or (isinstance(self.input_variables_dict[key][MEMC_iter][4][1], str
                                                               ) is False
                                                    and self.input_variables_dict[key][MEMC_iter][4][1] is not None):
                                            bad_input_variables_values_list.append(key)
                                    else:
                                        bad_input_variables_values_list.append(key)

                                    all_atom_names_and_res_pairs_keys_list = list(
                                        self.all_residues_unique_atom_name_dict.keys())
                                    # check that the atom names match the residues that exist
                                    if self.input_variables_dict[key][MEMC_iter][1] not in \
                                            all_atom_names_and_res_pairs_keys_list:
                                        bad_input_variables_values_list.append(key)

                                    elif self.input_variables_dict[key][MEMC_iter][1] in \
                                            all_atom_names_and_res_pairs_keys_list:

                                        if self.input_variables_dict[key][MEMC_iter][2][0] not in \
                                                self.all_residues_unique_atom_name_dict[
                                                    self.input_variables_dict[key][MEMC_iter][1]]:

                                            bad_input_variables_values_list.append(key)

                                        if self.input_variables_dict[key][MEMC_iter][2][1] not in \
                                                self.all_residues_unique_atom_name_dict[
                                                    self.input_variables_dict[key][MEMC_iter][1]]:
                                            bad_input_variables_values_list.append(key)

                                    if self.input_variables_dict[key][MEMC_iter][3] not in \
                                            all_atom_names_and_res_pairs_keys_list:
                                        bad_input_variables_values_list.append(key)

                                    elif self.input_variables_dict[key][MEMC_iter][3] in \
                                            all_atom_names_and_res_pairs_keys_list:
                                        if self.input_variables_dict[key][MEMC_iter][4][0] not in \
                                                self.all_residues_unique_atom_name_dict[
                                                    self.input_variables_dict[key][MEMC_iter][3]]:
                                            bad_input_variables_values_list.append(key)

                                        if self.input_variables_dict[key][MEMC_iter][4][1] not in \
                                                self.all_residues_unique_atom_name_dict[
                                                    self.input_variables_dict[key][MEMC_iter][3]]:
                                            bad_input_variables_values_list.append(key)

                                    if isinstance(self.input_variables_dict[key][MEMC_iter][0], int) is True:
                                        if self.input_variables_dict[key][MEMC_iter][0] <= 0:
                                            bad_input_variables_values_list.append(key)

                                else:
                                    bad_input_variables_values_list.append(key)

                if input_var_keys_list[var_iter] == key and key in possible_ensemble_variables_list:
                    self.MEMC_DataInput = self.input_variables_dict[key]

        # Error out and print the bad input values
        if len(bad_input_variables_values_list) > 0:
            self.input_error = True
            # create unique list
            bad_input_variables_values_set = set(bad_input_variables_values_list)
            bad_unique_input_variables_values_list = list(bad_input_variables_values_set)
            print_error_message = 'ERROR: The following input variables have ' \
                                  'bad values (check spelling and for empty spaces in the keys or that ' \
                                  'the values are in the correct form with the acceptable values): ' \
                                  '{}'.format(bad_unique_input_variables_values_list)
            raise ValueError(print_error_message)

        else:
            print("INFO: All the input variable passed the initial error checking")

        # auto calculate the best EqSteps (number of Equilbrium Steps) and Adj_Steps (number of AdjSteps Steps)
        self.EqSteps = scale_gen_freq_for_run_steps_int(self.EqSteps, self.RunSteps)

        self.AdjSteps = scale_gen_freq_for_run_steps_int(self.AdjSteps, self.RunSteps)

        # auto calculate the best RestartFreq  for the number of self.RunSteps
        self.RestartFreq = scale_gen_freq_for_run_steps_list_bool_int(self.RestartFreq, self.RunSteps)

        # auto calculate the best CheckpointFreq  for the number of self.RunSteps
        self.CheckpointFreq = scale_gen_freq_for_run_steps_list_bool_int(self.CheckpointFreq, self.RunSteps)

        # auto calculate the best CoordinatesFreq  for the number of self.RunSteps
        self.CoordinatesFreq = scale_gen_freq_for_run_steps_list_bool_int(self.CoordinatesFreq, self.RunSteps)

        # auto calculate the best ConsoleFreq  for the number of self.RunSteps
        self.ConsoleFreq = scale_gen_freq_for_run_steps_list_bool_int(self.ConsoleFreq, self.RunSteps)

        # auto calculate the best PressureCalc  for the number of self.RunSteps
        self.PressureCalc = scale_gen_freq_for_run_steps_list_bool_int(self.PressureCalc, self.RunSteps)

        # auto calculate the best BlockAverageFreq  for the number of self.RunSteps
        self.BlockAverageFreq = scale_gen_freq_for_run_steps_list_bool_int(self.BlockAverageFreq, self.RunSteps)

        # auto calculate the best HistogramFreq  for the number of self.RunSteps
        self.HistogramFreq = scale_gen_freq_for_run_steps_list_bool_int(self.HistogramFreq, self.RunSteps)

        # auto calculate the best SampleFreq  for the number of self.RunSteps
        self.SampleFreq = scale_gen_freq_for_run_steps_int(self.SampleFreq, self.RunSteps)

        # check to make sure the VDW FF (ParaTypeCHARMM) is not true for multiple ones
        # (i.e., ParaTypeCHARMM, ParaTypeMie, ParaTypeMARTINI)
        if (self.ParaTypeCHARMM is True and self.ParaTypeMie is True) \
                or (self.ParaTypeCHARMM is True and self.ParaTypeMARTINI is True) \
                or self.ParaTypeMie is True and self.ParaTypeMARTINI is True:
            self.input_error = True
            print_error_message = 'ERROR: there can only be 1 VDW type set to true.  Please set only one of the '\
                                  'ParaTypeCHARMM = {}, ParaTypeMie = {}, ParaTypeMARTINI = {} types ' \
                                  'to True'.format(self.ParaTypeCHARMM, self.ParaTypeMie, self.ParaTypeMARTINI)
            raise ValueError(print_error_message)
        elif self.ParaTypeCHARMM is True or self.ParaTypeMie is True \
                or self.ParaTypeMARTINI is True:
            if self.ParaTypeCHARMM is True:
                self.VDW_type = "ParaTypeCHARMM"
            elif self.ParaTypeMie is True:
                self.VDW_type = "ParaTypeMie"
            elif self.ParaTypeMARTINI is True:
                self.VDW_type = "ParaTypeMARTINI"
        else:
            self.input_error = True
            print_error_message = 'ERROR: There no VDW types that are set as True.  Please set only one of the '\
                                  'ParaTypeCHARMM = {}, ParaTypeMie = {}, ParaTypeMARTINI = {} types '\
                                  'to True'.format(self.ParaTypeCHARMM, self.ParaTypeMie, self.ParaTypeMARTINI)
            raise ValueError(print_error_message)

        # check to see if the moves sum up to 1
        if ensemble_type in ['NVT', 'GCMC']:
            if self.VolFreq != 0:
                self.input_error = True
                print_error_message = 'ERROR: The input variable VolFreq is non-zero (0). '\
                                      'VolFreq must be zero (0) for the "NVT", "GEMC_NVT", and "GCMC" ensembles.'
                raise ValueError(print_error_message)

        if ensemble_type in ['NVT', 'NPT']:
            if self.SwapFreq != 0 or self.MEMC_1Freq != 0 or self.MEMC_2Freq != 0 or self.MEMC_3Freq != 0:
                self.input_error = True
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

        if sum(moves_list) <= 1 + 10**(-13) and sum(moves_list) >= 1 - 10**(-13):
            print('INFO: The sum of the Monte Carlo move ratios = ' + str("{:.12f}".format(sum(moves_list))))

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
            self.input_error = True
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
            self.input_error = True
            print_error_message = 'ERROR: The values must be in this order RunSteps >= EqSteps >= AdjSteps ' \
                                  ' ({} >= {} >= {} )'.format(self.RunSteps, self.EqSteps, self.AdjSteps)
            raise ValueError(print_error_message)

        # check if both the ChemPot and Fugacity are not set to None.  Only one can be used
        if self.Fugacity is not None and self.ChemPot is not None \
                and self.ensemble_type == 'GCMC':
            self.input_error = True
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
            self.input_error = True
            print_error_message = 'ERROR: In the GCMC ensemble, neither Fugacity and ChemPot are provided ' \
                                  '(i.e., both are None). Add a dictionary for either the Fugacity or ' \
                                  'ChemPot and set the other variable to None. ' \
                                  'Note: Both the Fugacity or ChemPot and set to None by default'
            raise ValueError(print_error_message)

        # check that MEMC moves rations are > 0 if MEMC_DataInput is used
        if self.MEMC_DataInput is not None and (self.MEMC_1Freq == 0
                                                and self.MEMC_2Freq == 0
                                                and self.MEMC_3Freq == 0
                                                and self.IntraMEMC_1Freq == 0
                                                and self.IntraMEMC_2Freq == 0
                                                and self.IntraMEMC_3Freq == 0):
            self.input_error = True
            print_error_message = 'ERROR: The MEMC_DataInput variable is not equal to None, ' \
                                  'but all the MEMC move ratios are '\
                                  'zero (IntraMEMC_1Freq, MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, '\
                                  'IntraMEMC_3Freq, and MEMC_3Freq).'
            raise ValueError(print_error_message)

        if self.MEMC_DataInput is None and (self.MEMC_1Freq != 0
                                            or self.MEMC_2Freq != 0
                                            or self.MEMC_3Freq != 0
                                            or self.IntraMEMC_1Freq != 0
                                            or self.IntraMEMC_2Freq != 0
                                            or self.IntraMEMC_3Freq != 0):
            self.input_error = True
            print_error_message = 'ERROR: The MEMC_DataInput variable is equal to None, ' \
                                  'but at least one of the MEMC move ratios are '\
                                  'all non-zero (IntraMEMC_1Freq, MEMC_1Freq, IntraMEMC_2Freq, MEMC_2Freq, '\
                                  'IntraMEMC_3Freq, and MEMC_3Freq).'
            raise ValueError(print_error_message)

        # ensure the LargeKindBackBone and SmallKindBackBones are provided as appropriate for MEMC-1, MEMC-2, MEMC-3
        if self.MEMC_DataInput is not None and (self.MEMC_2Freq > 0 or self.IntraMEMC_2Freq > 0):
            for MEMC_2_i in range(0, len(self.MEMC_DataInput)):
                if self.MEMC_DataInput[MEMC_2_i][2][0] is None \
                    or self.MEMC_DataInput[MEMC_2_i][2][1] is None \
                        or self.MEMC_DataInput[MEMC_2_i][4][0] is None \
                        or self.MEMC_DataInput[MEMC_2_i][4][1] is None:

                    self.input_error = True
                    print_error_message = 'ERROR:  The  LargeKindBackBone and SmallKindBackBones unique ' \
                                          'atom names, strings, both must be provided when using the ' \
                                          'IntraMEMC-2Freq or MEMC-2Freq moves ' \
                                          '(i.e., the LargeKindBackBone and SmallKindBackBones can not be None). '
                    raise ValueError(print_error_message)

        if self.MEMC_DataInput is not None and (self.MEMC_3Freq > 0 or self.IntraMEMC_3Freq > 0):
            for MEMC_3_i in range(0, len(self.MEMC_DataInput)):
                if self.MEMC_DataInput[MEMC_3_i][2][0] is None or self.MEMC_DataInput[MEMC_3_i][2][1] is None:
                    self.input_error = True
                    print_error_message = 'ERROR:  The LargeKindBackBone unique atom names, strings, '\
                                          'both must be provided when using the IntraMEMC-3Freq or MEMC-3Freq moves '\
                                          '(i.e., the LargeKindBackBone can not be None).'
                    raise ValueError(print_error_message)

        # check that all required free energy values are provided
        if (self.FreeEnergyCalc is not None or self.MoleculeType is not None
            or self.InitialState is not None or self.LambdaVDW is not None) \
                and (self.FreeEnergyCalc is None or self.MoleculeType is None
                     or self.InitialState is None or self.LambdaVDW is None):
            self.input_error = True
            print_error_message = 'ERROR: To utilize the free energy calculations all the following ' \
                                  'variables need to be set, and not equal to None: ' \
                                  'FreeEnergyCalc, MoleculeType, InitialState, LambdaVDW.'
            raise ValueError(print_error_message)

        if self.LambdaVDW is not None and self.LambdaCoulomb is not None \
                and isinstance(self.LambdaVDW, list) is True and (isinstance(self.LambdaCoulomb, list)) is True:
            if len(self.LambdaVDW) == len(self.LambdaCoulomb):
                if self.InitialState + 1 <= len(self.LambdaVDW):
                    for lam_i in range(1, len(self.LambdaVDW)):
                        if self.LambdaVDW[lam_i] < self.LambdaVDW[lam_i - 1]:
                            self.input_error = True
                            print_error_message = 'ERROR: The LambdaVDW list is not in accending order.'
                            raise ValueError(print_error_message)
                        if self.LambdaCoulomb[lam_i] < self.LambdaCoulomb[lam_i - 1]:
                            self.input_error = True
                            print_error_message = 'ERROR:  The LambdaCoulomb list is not in accending order.'
                            raise ValueError(print_error_message)
                else:
                    self.input_error = True
                    print_error_message = 'ERROR: The InitialState integer is greater than the LambdaVDW and '\
                                          'LambdaCoulomb list length.  Note: the InitialState integer starts at 0.'
                    raise ValueError(print_error_message)
            else:
                self.input_error = True
                print_error_message = 'ERROR: The LambdaVDW and LambdaCoulomb list must be of equal length.'
                raise ValueError(print_error_message)

    # write the control file
    def write_conf_file(self, conf_filename):
        """
        Writes the GOMC control file.

        Parameters
        ----------
        conf_filename: str
            The path and file name for the control file name, with
            the .conf extension, or no extension.  If no extension is provided, the
            code will add the .conf extension to the provided file name.

        Returns
        ---------
        Writes the GOMC control file with the name provided via conf_filename

            If completed without errors: str, "PASSED
            If completed with errors :  None
        """

        # check to see if it is OK to proceed writing the control file
        if self.input_error is True:
            print_error_message = 'ERROR: The control file was not written as at least 1 input to the '\
                                  'control file writer was bad.'
            raise ValueError(print_error_message)

        date_time = datetime.datetime.today()

        self.conf_filename = conf_filename

        if isinstance(self.conf_filename, str) is False \
                or isinstance(self.conf_filename, str) is None:
            self.input_error = True
            print_error_message = 'ERROR: The control file name (conf_filename) is not provided as a string. '
            raise ValueError(print_error_message)

        if os.path.splitext(self.conf_filename)[1] == '.conf':
            self.conf_filename = conf_filename
            print('INFO: the correct extension for the control file was provided in the file name, .conf '
                  'with control file name = ' + str(self.conf_filename))
        elif os.path.splitext(self.conf_filename)[1] == '':
            self.conf_filename = self.conf_filename + '.conf'
            print('INFO: No extension name was provided for the control file. Therefore, the proper '
                  'extension, .conf, was added.  The new total control file name = {}'.format(self.conf_filename))
        else:
            self.input_error = True
            print_error_message = 'ERROR: No extension name or the wrong extension name was provided. ' \
                                  'Please enter a proper extension name, .conf or no extension in the conf_filename '\
                                  'The control file as provided name = {}'.format(self.conf_filename)
            raise ValueError(print_error_message)

        # get and open the control file file
        data_control_file = open(self.conf_filename, 'w')

        data_control_file.write('#######################################################'
                                '#########################################\n')
        data_control_file.write("##  This file ({}) - was created by mBuild "
                                "using the on {}\n".format(self.conf_filename, date_time))
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
        data_control_file.write('Restart \t {}\n'.format(self.Restart))
        data_control_file.write('\n')
        data_control_file.write('RestartCheckpoint \t {}\n'.format(self.RestartCheckpoint))
        data_control_file.write('\n')
        data_control_file.write('####################################\n')
        data_control_file.write('# kind {RESTART, RANDOM, INTSEED}\n')
        data_control_file.write('####################################\n')
        if self.PRNG == "RANDOM":
            data_control_file.write('PRNG \t\t {}\n'.format(self.PRNG))
        elif isinstance(self.PRNG, int):
            data_control_file.write('PRNG \t\t ' + 'INTSEED \n')
            data_control_file.write('Random_Seed \t {}\n'.format(self.PRNG))
        data_control_file.write(' \n')
        data_control_file.write('####################################\n')
        data_control_file.write('# FORCE FIELD\n')
        data_control_file.write('####################################\n')
        data_control_file.write('{}\t\t {}\n'.format(self.VDW_type, True))
        data_control_file.write(' \n')
        data_control_file.write('Parameters \t\t {}\n'.format(self.ff_filename))
        data_control_file.write('####################################\n')
        data_control_file.write('# INPUT PDB FILES\n')
        data_control_file.write('####################################\n')
        if self.ensemble_type in ['NVT', 'NPT', 'GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Coordinates 0 \t\t {}\n'.format(self.Coordinates_box_0))
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Coordinates 1 \t\t {}\n'.format(self.Coordinates_box_1))
        data_control_file.write(' \n')
        data_control_file.write('####################################\n')
        data_control_file.write('# INPUT PSF FILES\n')
        data_control_file.write('####################################\n')
        if self.ensemble_type in ['NVT', 'NPT', 'GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Structure 0 \t\t {}\n'.format(self.Structures_box_0))
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('Structure 1 \t\t {}\n'.format(self.Structures_box_1))
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
        data_control_file.write('Temperature \t\t {}\n'.format(self.Temperature))
        if self.ensemble_type in ['GEMC_NPT', 'NPT']:
            data_control_file.write('Pressure \t\t {}\n'.format(self.Pressure))
            data_control_file.write('useConstantArea \t {}\n'.format(self.useConstantArea))

        if self.ensemble_type in ['GEMC_NPT'] and self.FixVolBox0 is True:
            data_control_file.write('FixVolBox0 \t\t {}\n'.format(self.FixVolBox0))

        if self.ensemble_type in ['GCMC'] and self.ChemPot is not None \
                and self.Fugacity is None:
            chem_pot_dict_key_list = dict_keys_to_list(self.ChemPot)
            for chem_pot_iter in range(0, len(chem_pot_dict_key_list)):
                chem_pot_residue_iter = chem_pot_dict_key_list[chem_pot_iter]
                data_control_file.write('ChemPot \t\t {} \t\t {}\n'.format(chem_pot_residue_iter,
                                                                         self.ChemPot[chem_pot_residue_iter]))

        if self.ensemble_type in ['GCMC'] and self.Fugacity is not None \
                and self.ChemPot is None:
            fugacity_iter_dict_key_list = dict_keys_to_list(self.Fugacity)
            for fugacity_iter in range(0, len(fugacity_iter_dict_key_list)):
                fugacity_residue_iter = fugacity_iter_dict_key_list[fugacity_iter]
                data_control_file.write('Fugacity \t\t {} \t\t {}\n'.format(fugacity_residue_iter,
                                                                          self.Fugacity[fugacity_residue_iter]))

        data_control_file.write(' \n')
        data_control_file.write('Potential \t\t {}\n'.format(self.Potential))
        data_control_file.write('LRC \t\t\t {}\n'.format(self.LRC))
        data_control_file.write('Rcut \t\t\t {}\n'.format(self.Rcut))
        data_control_file.write('RcutLow \t\t {}\n'.format(self.RcutLow))
        if self.Potential == 'SWITCH':
            data_control_file.write('Rswitch \t\t {}\n'.format(self.Rswitch))
        data_control_file.write('Exclude \t\t {}\n'.format(self.Exclude))
        if self.VDWGeometricSigma is True:
            data_control_file.write('VDWGeometricSigma \t {}\n'.format(self.VDWGeometricSigma))
        data_control_file.write(' \n')

        data_control_file.write('#############################\n')
        data_control_file.write('# ELECTROSTATIC   \n')
        data_control_file.write('#############################\n')
        data_control_file.write('Ewald \t\t {}\n'.format(self.Ewald))
        data_control_file.write('ElectroStatic \t {}\n'.format(self.ElectroStatic))
        data_control_file.write('CachedFourier \t {}\n'.format(self.CachedFourier))
        data_control_file.write('Tolerance \t {}\n'.format(format(self.Tolerance, '.12f')))
        if self.VDW_type in ["ParaTypeMARTINI"]:
            data_control_file.write('Dielectric \t {}\n'.format(format(self.Dielectric)))
        data_control_file.write('1-4scaling \t {}\n'.format(self.coul_1_4))
        data_control_file.write(' \n')
        if self.RcutCoulomb_box_0 is not None:
            data_control_file.write('RcutCoulomb 0 \t {}\n'.format(self.RcutCoulomb_box_0))
        if self.RcutCoulomb_box_1 is not None \
                and self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write('RcutCoulomb 1 \t {}\n'.format(self.RcutCoulomb_box_1))
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# PRESSURE CALCULATION\n')
        data_control_file.write('################################\n')
        data_control_file.write('PressureCalc \t {} \t\t {}\n'.format(self.PressureCalc[0], self.PressureCalc[1]))
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# STEPS \n')
        data_control_file.write('################################\n')
        data_control_file.write('RunSteps \t {}\n'.format(self.RunSteps))
        data_control_file.write('EqSteps \t {}\n'.format(self.EqSteps))
        data_control_file.write('AdjSteps \t {}\n'.format(self.AdjSteps))
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# MOVE FREQUENCY \n')
        data_control_file.write('################################\n')
        data_control_file.write('DisFreq \t\t {}\n'.format(self.DisFreq))
        data_control_file.write('RotFreq \t\t {}\n'.format(self.RotFreq))
        data_control_file.write('IntraSwapFreq \t\t {}\n'.format(self.IntraSwapFreq))
        data_control_file.write('SwapFreq \t\t {}\n'.format(self.SwapFreq))
        data_control_file.write('RegrowthFreq \t\t {}\n'.format(self.RegrowthFreq))
        data_control_file.write('CrankShaftFreq \t\t {}\n'.format(self.CrankShaftFreq))
        data_control_file.write('VolFreq \t\t {}\n'.format(self.VolFreq))
        data_control_file.write('MultiParticleFreq \t {}\n'.format(self.MultiParticleFreq))
        data_control_file.write('IntraMEMC-1Freq \t {}\n'.format(self.IntraMEMC_1Freq))
        data_control_file.write('MEMC-1Freq \t\t {}\n'.format(self.MEMC_1Freq))
        data_control_file.write('IntraMEMC-2Freq \t {}\n'.format(self.IntraMEMC_2Freq))
        data_control_file.write('MEMC-2Freq \t\t {}\n'.format(self.MEMC_2Freq))
        data_control_file.write('IntraMEMC-3Freq \t {}\n'.format(self.IntraMEMC_3Freq))
        data_control_file.write('MEMC-3Freq \t\t {}\n'.format(self.MEMC_3Freq))
        data_control_file.write(' \n')

        # sort and print the MEMC data if MEMC is used for the simulation
        if self.MEMC_DataInput is not None and (self.MEMC_1Freq > 0
                                                or self.MEMC_2Freq > 0
                                                or self.MEMC_3Freq > 0
                                                or self.IntraMEMC_1Freq > 0
                                                or self.IntraMEMC_2Freq > 0
                                                or self.IntraMEMC_3Freq > 0):

            ExchangeRatio_list = []
            ExchangeLargeKind_list = []
            LargeKindBackBone_list = []
            ExchangeSmallKind_list = []
            SmallKindBackBone_list = []

            for memc_i in range(0, len(self.MEMC_DataInput)):
                ExchangeRatio_list.append(str("{} \t\t".format(self.MEMC_DataInput[memc_i][0])))
                ExchangeLargeKind_list.append("{} \t\t".format(self.MEMC_DataInput[memc_i][1]))
                LargeKindBackBone_list.append("{}   {} \t".format(self.MEMC_DataInput[memc_i][2][0],
                                                                  self.MEMC_DataInput[memc_i][2][1]))
                ExchangeSmallKind_list.append("{} \t\t".format(self.MEMC_DataInput[memc_i][3]))
                SmallKindBackBone_list.append("{}   {} \t".format(self.MEMC_DataInput[memc_i][4][0],
                                                                  self.MEMC_DataInput[memc_i][4][1]))

            ExchangeRatio_str = ''.join(ExchangeRatio_list)
            ExchangeLargeKind_str = ''.join(ExchangeLargeKind_list)
            LargeKindBackBone_str = ''.join(LargeKindBackBone_list)
            ExchangeSmallKind_str = ''.join(ExchangeSmallKind_list)
            SmallKindBackBone_str = ''.join(SmallKindBackBone_list)

            data_control_file.write('###############################\n')
            data_control_file.write('# MEMC PARAMETER \n')
            data_control_file.write('###############################\n')
            data_control_file.write('ExchangeVolumeDim \t {}\t{}\t{} \n'.format(self.ExchangeVolumeDim[0],
                                                                                self.ExchangeVolumeDim[1],
                                                                                self.ExchangeVolumeDim[2]
                                                                                )
                                    )
            data_control_file.write('ExchangeRatio \t\t {}\n'.format(ExchangeRatio_str))
            data_control_file.write('ExchangeLargeKind \t {}\n'.format(ExchangeLargeKind_str))
            data_control_file.write('ExchangeSmallKind \t {}\n'.format(ExchangeSmallKind_str))
            if self.MEMC_DataInput is not None and (self.MEMC_2Freq > 0
                                                    or self.IntraMEMC_2Freq > 0
                                                    or self.MEMC_3Freq > 0
                                                    or self.IntraMEMC_3Freq > 0):
                data_control_file.write('LargeKindBackBone \t {}\n'.format(LargeKindBackBone_str))
            if self.MEMC_DataInput is not None and (self.MEMC_2Freq > 0
                                                    or self.IntraMEMC_2Freq > 0):
                data_control_file.write('SmallKindBackBone \t {}\n'.format(SmallKindBackBone_str))

        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# BOX DIMENSION #, X, Y, Z    (only orthoganal boxes are currently '
                                'available in this control file writer)\n')
        data_control_file.write('################################\n')
        data_control_file.write(
            'CellBasisVector1 0 \t {}\t\t0.00\t\t0.00\n'.format(self.x_dim_box_0))
        data_control_file.write(
            'CellBasisVector2 0 \t 0.00\t\t{}\t\t0.00\n'.format(self.y_dim_box_0))
        data_control_file.write(
            'CellBasisVector3 0 \t 0.00\t\t0.00\t\t{}\n'.format(self.z_dim_box_0))
        data_control_file.write(' \n')
        if self.ensemble_type in ['GEMC_NPT', 'GEMC_NVT', 'GCMC']:
            data_control_file.write(
                'CellBasisVector1 1 \t {}\t\t0.00\t\t0.00\n'.format(self.x_dim_box_1))
            data_control_file.write(
                'CellBasisVector2 1 \t 0.00\t\t{}\t\t0.00\n'.format(self.y_dim_box_1))
            data_control_file.write(
                'CellBasisVector3 1 \t 0.00\t\t0.00\t\t{}\n'.format(self.z_dim_box_1))
            data_control_file.write(' \n')

        if (self.ensemble_type in ['NPT', 'NVT']) \
                and self.FreeEnergyCalc is not None and self.MoleculeType is not None \
                and self.InitialState is not None and self.LambdaVDW is not None \
                and self.LambdaCoulomb is not None:

            # make list for number of states, LambdaVDW, and Lambda_Coul, and convert the list to string for printing
            Lambda_states_list = []
            Lambda_VDW_list = []
            Lambda_Coul_list = []
            for lamda_i in range(0, len(self.LambdaVDW)):
                Lambda_states_list.append(str(lamda_i))
                Lambda_VDW_list.append(str(self.LambdaVDW[lamda_i]))
                Lambda_Coul_list.append(str(self.LambdaCoulomb[lamda_i]))

            Lambda_states_str = '\t'.join(Lambda_states_list)
            Lambda_VDW_str = '\t'.join(Lambda_VDW_list)
            Lambda_Coul_str = '\t'.join(Lambda_Coul_list)

            data_control_file.write('##############################\n')
            data_control_file.write('# FREE ENERGY PARAMETERS (only available in NPT and NVT ensembles) \n')
            data_control_file.write('##############################\n')
            data_control_file.write('FreeEnergyCalc \t {}\t\t{}\n'.format(self.FreeEnergyCalc[0],
                                                                          self.FreeEnergyCalc[1]))
            data_control_file.write('MoleculeType \t {}\t\t{}\n'.format(self.MoleculeType[0],
                                                                        self.MoleculeType[1]))
            data_control_file.write('InitialState \t {}\n'.format(self.InitialState))
            data_control_file.write('ScalePower \t {}\n'.format(self.ScalePower))
            data_control_file.write('ScaleAlpha \t {}\n'.format(self.ScaleAlpha))
            data_control_file.write('MinSigma \t {}\n'.format(self.MinSigma))
            data_control_file.write('ScaleCoulomb \t {}\n'.format(self.ScaleCoulomb))
            data_control_file.write('# States \t {}\n'.format(Lambda_states_str))
            data_control_file.write('LambdaVDW \t {}\n'.format(Lambda_VDW_str))
            data_control_file.write('LambdaCoulomb \t {}\n'.format(Lambda_Coul_str))
            data_control_file.write(' \n')

        if (self.ensemble_type in ['NPT', 'NVT']) \
                and self.FreeEnergyCalc is not None and self.MoleculeType is not None \
                and self.InitialState is not None and self.LambdaVDW is not None:

            # make list for number of states, LambdaVDW, and Lambda_Coul, and convert the list to string for printing
            Lambda_states_list = []
            Lambda_VDW_list = []
            Lambda_Coul_list = []
            for lamda_i in range(0, len(self.LambdaVDW)):
                Lambda_states_list.append(str(lamda_i))
                Lambda_VDW_list.append(str(self.LambdaVDW[lamda_i]))
                if self.LambdaCoulomb is not None:
                    Lambda_Coul_list.append(str(self.LambdaCoulomb[lamda_i]))

            Lambda_states_str = '\t'.join(Lambda_states_list)
            Lambda_VDW_str = '\t'.join(Lambda_VDW_list)
            if self.LambdaCoulomb is not None:
                Lambda_Coul_str = '\t'.join(Lambda_Coul_list)

            data_control_file.write('##############################\n')
            data_control_file.write('# FREE ENERGY PARAMETERS (only available in NPT and NVT ensembles) \n')
            data_control_file.write('##############################\n')
            data_control_file.write('FreeEnergyCalc \t {}\t\t{}\n'.format(self.FreeEnergyCalc[0],
                                                                          self.FreeEnergyCalc[1]))
            data_control_file.write('MoleculeType \t {}\t\t{}\n'.format(self.MoleculeType[0],
                                                                        self.MoleculeType[1]))
            data_control_file.write('InitialState \t {}\n'.format(self.InitialState))
            data_control_file.write('ScalePower \t {}\n'.format(self.ScalePower))
            data_control_file.write('ScaleAlpha \t {}\n'.format(self.ScaleAlpha))
            data_control_file.write('MinSigma \t {}\n'.format(self.MinSigma))
            data_control_file.write('ScaleCoulomb \t {}\n'.format(self.ScaleCoulomb))
            data_control_file.write('# States \t {}\n'.format(Lambda_states_str))
            data_control_file.write('LambdaVDW \t {}\n'.format(Lambda_VDW_str))
            if self.LambdaCoulomb is not None:
                data_control_file.write('LambdaCoulomb \t {}\n'.format(Lambda_Coul_str))
            data_control_file.write(' \n')

        data_control_file.write('##############################\n')
        data_control_file.write('# CBMC TRIALS \n')
        data_control_file.write('##############################\n')
        data_control_file.write('CBMC_First \t {}\n'.format(self.CBMC_First))
        data_control_file.write('CBMC_Nth \t {}\n'.format(self.CBMC_Nth))
        data_control_file.write('CBMC_Ang \t {}\n'.format(self.CBMC_Ang))
        data_control_file.write('CBMC_Dih \t {}\n'.format(self.CBMC_Dih))
        data_control_file.write(' \n')

        data_control_file.write('############################################################################\n')
        data_control_file.write('#  =======-------------------- OUTPUT --------------------------=========== \n')
        data_control_file.write('############################################################################\n')
        data_control_file.write(' \n')

        data_control_file.write('##########################\n')
        data_control_file.write('# statistics filename add\n')
        data_control_file.write('##########################\n')
        data_control_file.write('OutputName \t {}\n'.format(self.OutputName))
        data_control_file.write(' \n')

        data_control_file.write('#####################################\n')
        data_control_file.write('# enable, frequency \n')
        data_control_file.write('#####################################\n')
        data_control_file.write('RestartFreq  \t\t {}\t\t{}\n'.format(self.RestartFreq[0],
                                                                      self.RestartFreq[1]))
        data_control_file.write('CheckpointFreq  \t {}\t\t{}\n'.format(self.CheckpointFreq[0],
                                                                       self.CheckpointFreq[1]))
        data_control_file.write('CoordinatesFreq  \t {}\t\t{}\n'.format(self.CoordinatesFreq[0],
                                                                        self.CoordinatesFreq[1]))
        data_control_file.write('ConsoleFreq  \t\t {}\t\t{}\n'.format(self.ConsoleFreq[0],
                                                                      self.ConsoleFreq[1]))
        data_control_file.write('BlockAverageFreq  \t {}\t\t{}\n'.format(self.BlockAverageFreq[0],
                                                                         self.BlockAverageFreq[1]))
        data_control_file.write('HistogramFreq  \t\t {}\t\t{}\n'.format(self.HistogramFreq[0],
                                                                        self.HistogramFreq[1]))
        data_control_file.write(' \n')

        data_control_file.write('################################\n')
        data_control_file.write('# OutHistSettings \n')
        data_control_file.write('################################\n')
        data_control_file.write('DistName \t {}\n'.format(self.DistName))
        data_control_file.write('HistName \t {}\n'.format(self.HistName))
        data_control_file.write('RunNumber \t {}\n'.format(self.RunNumber))
        data_control_file.write('RunLetter \t {}\n'.format(self.RunLetter))
        data_control_file.write('SampleFreq \t {}\n'.format(self.SampleFreq))
        data_control_file.write(' \n')

        data_control_file.write('################################## \n')
        data_control_file.write('# enable: blk avg., fluct. \n')
        data_control_file.write('################################## \n')
        data_control_file.write('OutEnergy \t\t {}\t\t{}\n'.format(self.OutEnergy[0], self.OutEnergy[1]))
        data_control_file.write('OutPressure \t\t {}\t\t{}\n'.format(self.OutPressure[0], self.OutPressure[1]))
        data_control_file.write('OutMolNumber \t\t {}\t\t{}\n'.format(self.OutMolNumber[0], self.OutMolNumber[1]))
        data_control_file.write('OutDensity \t\t {}\t\t{}\n'.format(self.OutDensity[0], self.OutDensity[1]))
        data_control_file.write('OutVolume \t\t {}\t\t{}\n'.format(self.OutVolume[0], self.OutVolume[1]))
        data_control_file.write('OutSurfaceTension \t {}\t\t{}\n'.format(self.OutSurfaceTension[0],
                                                                         self.OutSurfaceTension[1]))
        data_control_file.write('\n')
        data_control_file.write('\n')

        return "GOMC_CONTROL_FILE_WRITTEN"


    def ck_input_variable_true_or_false(self,
                                        input_variables_dict,
                                        key,
                                        bad_user_variable_list):
        """
        Checks if the input variable is either True for False.
        If not either True or False, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_user_variable_list: list
            A list to append with the bad dict_key user inputs
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.

        Returns
        ---------
        bad_user_variable_list : list,
            A list to append with the bad dict_key user inputs,
            which is appended upon detecting a bad user input variable.
            Note: This list is intented to be printed with all the bad input variables
            so the user can fix them after upon a failed GOMC conf file writing attempt.
        """

        if input_variables_dict[key] is not True \
                and input_variables_dict[key] is not False \
                or (str(input_variables_dict[key]) == str(True)
                    and str(input_variables_dict[key]) == str(False)):
            bad_user_variable_list.append(key)

    def ck_input_variable_int_or_float_zero_or_greater(self,
                                                       input_variables_dict,
                                                       key,
                                                       bad_input_variables_values_list):
        """
        Checks if the input variable is an integer or float is zero or greater ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key.

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None and isinstance(
                input_variables_dict[key], int) is not True \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_int_zero_or_greater(self,
                                              input_variables_dict,
                                              key,
                                              bad_input_variables_values_list):
        """
        Checks if the input variable is an integer is zero or greater ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None and isinstance(
                input_variables_dict[key], int) is not True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_float_zero_or_greater(self,
                                                input_variables_dict,
                                                key,
                                                bad_input_variables_values_list):
        """
        Checks if the input variable is a float is zero or greater ( value >= 0 ).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] < 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_int_or_float_greater_zero(self,
                                                    input_variables_dict,
                                                    key,
                                                    bad_input_variables_values_list):
        """
        Checks if the input variable is an integer or float is greater than zero ( value > 0).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None and isinstance(
                input_variables_dict[key], int) is not True \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_int_greater_zero(self,
                                           input_variables_dict,
                                           key,
                                           bad_input_variables_values_list):
        """
        Checks if the input variable is an integer greater than zero ( value > 0)..
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None and isinstance(
                input_variables_dict[key], int) is not True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_float_greater_zero(self,
                                             input_variables_dict,
                                             key,
                                             bad_input_variables_values_list):
        """
        Checks if the input variable is a float greater than zero ( value > 0).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] <= 0 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_float_greater_zero_less_1(self,
                                                    input_variables_dict,
                                                    key,
                                                    bad_input_variables_values_list):
        """
        Checks if the input variable is a float greater than zero and less than 1
        ( 0 < value < 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] <= 0 \
                or input_variables_dict[key] >= 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_int_or_float_zero_to_1(self,
                                                 input_variables_dict,
                                                 key,
                                                 bad_input_variables_values_list):
        """
        Checks if the input variable is an integer or float from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], int) is not True \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_float_zero_to_1(self,
                                          input_variables_dict,
                                          key,
                                          bad_input_variables_values_list):
        """
        Checks if the input variable is a float from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], float) is not True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_int_zero_to_1(self,
                                        input_variables_dict,
                                        key,
                                        bad_input_variables_values_list):
        """
        Checks if the input variable is an integer from 0 to 1 ( 0 =< value <= 1).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], int) is not True \
                or input_variables_dict[key] < 0 \
                or input_variables_dict[key] > 1 \
                or str(input_variables_dict[key]) == str(True) \
                or str(input_variables_dict[key]) == str(False):
            bad_input_variables_values_list.append(key)

    def ck_input_variable_list_bool_int_zero_or_greater(self,
                                                        input_variables_dict,
                                                        key,
                                                        bad_input_variables_values_list):
        """
        Checks if the input variable is a list with a bool and integer 0 or greater
        ([bool, int >= 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if isinstance(input_variables_dict[key], list) is False:
            bad_input_variables_values_list.append(key)
        elif isinstance(input_variables_dict[key], list) is True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] is not True
                        and input_variables_dict[key][0] is not False) \
                    or isinstance(input_variables_dict[key][1], int) is not True \
                    or input_variables_dict[key][1] < 0 \
                    or (str(input_variables_dict[key][0]) != str(True)
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or str(input_variables_dict[key][1]) == str(True) \
                    or str(input_variables_dict[key][1]) == str(False):
                bad_input_variables_values_list.append(key)

    def ck_input_variable_list_bool_int_greater_zero(self,
                                                     input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list):
        """
        Checks if the input variable is a list with a bool and integer greater than zero
        ([bool, int > 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if isinstance(input_variables_dict[key], list) is False:
            bad_input_variables_values_list.append(key)
        elif isinstance(input_variables_dict[key], list) is True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] is not True
                        and input_variables_dict[key][0] is not False) \
                    or isinstance(input_variables_dict[key][1], int) is not True \
                    or input_variables_dict[key][1] <= 0 \
                    or (str(input_variables_dict[key][0]) != str(True)
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or str(input_variables_dict[key][1]) == str(True) \
                    or str(input_variables_dict[key][1]) == str(False):
                bad_input_variables_values_list.append(key)

    def ck_input_variable_list_residue_str_int_greater_zero(self,
                                                            input_variables_dict,
                                                            key,
                                                            bad_input_variables_values_list):
        """
        Checks if the input variable is a list with a str and integer greater than zero
        ([str, int > 0 ]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if isinstance(input_variables_dict[key], list) is False:
            bad_input_variables_values_list.append(key)
        elif isinstance(input_variables_dict[key], list) is True:
            if len(input_variables_dict[key]) == 2:
                if isinstance(input_variables_dict[key], list) is True:
                    if isinstance(input_variables_dict[key][0], str) is False \
                            or isinstance(input_variables_dict[key][1], int) is False:
                        bad_input_variables_values_list.append(key)
                    elif isinstance(input_variables_dict[key][0], str) is True:
                        if len(input_variables_dict[key][0]) > 4 \
                                or input_variables_dict[key][0] not in self.residues:
                            bad_input_variables_values_list.append(key)

                    if isinstance(input_variables_dict[key][1], int) is True \
                            and input_variables_dict[key][1] <= 0:
                        bad_input_variables_values_list.append(key)

    def ck_input_variable_list_bool_bool(self,
                                         input_variables_dict,
                                         key,
                                         bad_input_variables_values_list):
        """
        Checks if the input variable is a list with a 2 booleans  ([bool, bool]).
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if isinstance(input_variables_dict[key], list) is False:
            bad_input_variables_values_list.append(key)
        elif isinstance(input_variables_dict[key], list) is True:
            if len(input_variables_dict[key]) != 2 \
                    or (input_variables_dict[key][0] is not True
                        and input_variables_dict[key][0] is not False) \
                    or (input_variables_dict[key][1] is not True
                        and input_variables_dict[key][1] is not False) \
                    or (str(input_variables_dict[key][0]) != str(True)
                        and str(input_variables_dict[key][0]) != str(False)) \
                    or (str(input_variables_dict[key][1]) != str(True)
                        and str(input_variables_dict[key][1]) != str(False)):
                bad_input_variables_values_list.append(key)

    def ck_input_variable_list_of_floats_zero_to_1(self,
                                                   input_variables_dict,
                                                   key,
                                                   bad_input_variables_values_list):
        """
        Checks if the input variable is a list of floats between 0 and 1 ([0,0, 0.1, ..., 1.0]).
        Note: the list can be of any length with 0.0 <= float <= 1.0
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None:
            if isinstance(input_variables_dict[key], list) is False:
                bad_input_variables_values_list.append(key)
            elif isinstance(input_variables_dict[key], list) is True:
                if len(input_variables_dict[key]) >= 1:
                    for Lambda_i in range(0, len(input_variables_dict[key])):
                        if isinstance(input_variables_dict[key][Lambda_i], float) is False:
                            bad_input_variables_values_list.append(key)
                        elif isinstance(input_variables_dict[key][Lambda_i], float) is True:
                            if input_variables_dict[key][Lambda_i] < 0.0 \
                                    or input_variables_dict[key][Lambda_i] > 1.0:
                                bad_input_variables_values_list.append(key)
                else:
                    bad_input_variables_values_list.append(key)

    def ck_input_variable_GCMC_dict_str_int_or_float(self,
                                                     input_variables_dict,
                                                     key,
                                                     bad_input_variables_values_list):
        """
        Checks if the input variable is a dictionary with a key = string and
        value = integer or float  ({'str_1' : integer_1 or float_1, ....,
        'str_x' : integer_x or float_x }).
        Note: the dictionary can be of any length
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], dict) is not True:
            bad_input_variables_values_list.append(key)

        elif isinstance(input_variables_dict[key], dict) is True:
            keys_list = dict_keys_to_list(input_variables_dict[key])
            for keys_iter_No in range(0, len(keys_list)):
                key_iter = keys_list[keys_iter_No]
                value_iter = input_variables_dict[key][key_iter]

                if key_iter not in self.residues:
                    bad_input_variables_values_list.append(key)

                if isinstance(key_iter, str) is False:
                    bad_input_variables_values_list.append(key)
                elif isinstance(key_iter, str) is True:
                    if len(key_iter) > 4:
                        bad_input_variables_values_list.append(key)

                if (isinstance(value_iter, int) is False
                    and isinstance(value_iter, float) is False) \
                        or str(value_iter) == str(True) \
                        or str(value_iter) == str(False):
                    bad_input_variables_values_list.append(key)

    def ck_input_variable_GCMC_dict_str_int_or_float_zero_or_greater(self,
                                                                     input_variables_dict,
                                                                     key,
                                                                     bad_input_variables_values_list):
        """
        Checks if the input variable is a dictionary with a key = string and
        value = integer or float zero or greater  ({'str_1' : integer_1 or float_1 (>= 0), ....,
        'str_x' : integer_x or float_x (>= 0)}).
        Note: the dictionary can be of any length
        If not, the provided list is appended with the bad with the dict_key

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if input_variables_dict[key] is not None \
                and isinstance(input_variables_dict[key], dict) is not True:
            bad_input_variables_values_list.append(key)

        elif isinstance(input_variables_dict[key], dict) is True:
            keys_list = dict_keys_to_list(input_variables_dict[key])
            for keys_iter_No in range(0, len(keys_list)):
                key_iter = keys_list[keys_iter_No]
                value_iter = input_variables_dict[key][key_iter]

                if key_iter not in self.residues:
                    bad_input_variables_values_list.append(key)

                if isinstance(key_iter, str) is False:
                    bad_input_variables_values_list.append(key)
                elif isinstance(key_iter, str) is True:
                    if len(key_iter) > 4:
                        bad_input_variables_values_list.append(key)

                if (isinstance(value_iter, int) is False
                    and isinstance(value_iter, float) is False) \
                        or str(value_iter) == str(True) \
                        or str(value_iter) == str(False) \
                        or value_iter < 0:
                    bad_input_variables_values_list.append(key)

    def ck_input_variable_str_with_no_spaces(self,
                                             input_variables_dict,
                                             key,
                                             bad_input_variables_values_list):
        """
        Checks if the input variable is a string with no spaces.
        If not, the provided list is appended with the bad with the dict_key.

        Parameters
        ----------
        input_variables_dict: dict
            The user input variable dictionary
        key: str
            Dictionary key for the user provided input variable list
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.

        Returns
        ---------
        bad_input_variables_values_list: list
            A list to append with the bad variable user inputs
            so the user can see which variable input values are bad.
        """

        if isinstance(input_variables_dict[key], str) is True:
            no_spaces_in_OutputName_string = " " in input_variables_dict[key]
        if isinstance(input_variables_dict[key], str) is not True \
                or no_spaces_in_OutputName_string is True:
            bad_input_variables_values_list.append(key)

def scale_gen_freq_for_run_steps_list_bool_int(charmm_variable,
                                               run_steps
                                               ):
    """
    Scales the frequency of the output to a a more realistic value,
    if the output frequency does not make sense based on the
    total number of simulation run steps.

    Parameters
    ----------
    charmm_variable : GOMCControl object variable list, [bool, int]
        This only for the frequency output variables of the GOMCControl object
    run_steps : int (>0), must be an integer greater than zero.
        This is the GOMCControl object varaible which sets the total number of simulation steps.
        This should be the RunSteps variable in the GOMCControl object.

    Returns
    ---------
    charmm_variable : GOMCControl object variable list [bool, int]
        A rescaled and appropriate value for the frequency output variables of the
        GOMCControl object, based on the RunSteps in the simulation.
    """

    set_max_steps_charmm_variable = charmm_variable[1]

    if run_steps / 10 >= set_max_steps_charmm_variable and run_steps / 10 >= 1:
        charmm_variable[1] = int(set_max_steps_charmm_variable)
    elif run_steps / 10 >= 1:
        charmm_variable[1] = int(run_steps / 10)
    else:
        charmm_variable[1] = int(1)

    return charmm_variable

def scale_gen_freq_for_run_steps_int(charmm_variable,
                                     run_steps
                                     ):
    """
    Scales the frequency of the output to a a more realistic value,
    if the output frequency does not make sense based on the
    total number of simulation run steps.

    Parameters
    ----------
    charmm_variable : GOMCControl object variable, int
        This only for the frequency output variables of the GOMCControl object
    run_steps : int (>0), must be an integer greater than zero.
        This is the GOMCControl object varaible which sets the total number of simulation steps.
        This should be the RunSteps variable in the GOMCControl object.

    Returns
    ---------
    charmm_variable : GOMCControl object variable list, int
        A rescaled and appropriate value for the frequency output variables of the
        GOMCControl object, based on the RunSteps in the simulation.
    """

    set_max_steps_charmm_variable = charmm_variable

    if run_steps / 10 >= set_max_steps_charmm_variable and run_steps / 10 >= 1:
        charmm_variable = int(set_max_steps_charmm_variable)
    elif run_steps / 10 >= 1:
        charmm_variable = int(run_steps / 10)
    else:
        charmm_variable = int(1)

    return charmm_variable

def ck_box_dim_is_float_or_int_greater_0(charmm_box_dimension,
                                         dimension,
                                         box_no,
                                         ensemble_type
                                         ):
    """
    Scales the frequency of the output to a a more realistic value,
    if the output frequency does not make sense based on the
    total number of simulation run steps.

    Parameters
    ----------
    charmm_box_dimension : GOMCControl object box dimension variable
        This is the variable that contains the box input x, y, or z dimensions for box 0 or 1.
    dimension : str (Only enter 'x', 'y', or 'z')
        This is the dimension of the charmm_box_dimension variable.
        Only enter 'x', 'y', or 'z', but it will not error out if you do not.
    box_no : int (only 0 or 1)
        This is the box number which is defined in the mbuild.Charmm object
        Only enter only 0 or 1, but it will not error out if you do not.
    ensemble_type : str, valid options are 'NVT', 'NPT', 'GEMC_NVT', 'GEMC_NPT', 'GCMC'
        The ensemble type of the simulation.

    Returns
    ---------
    If charmm_box_dimension is an int or float (>0) : None
    If charmm_box_dimension is not an int or float (>0) : raise ValueError

    """

    if (isinstance(charmm_box_dimension, int) is False
        and isinstance(charmm_box_dimension, float) is False) \
            or charmm_box_dimension <= 0:
        print_error_message = "ERROR: The {}-dimension for box {} is, {} , not an integer, float, "
        "or is <= 0.".format(dimension, box_no, charmm_box_dimension)
        raise ValueError(print_error_message)

    if (ensemble_type in ['GEMC_NVT', 'GEMC_NPT', 'GCMC']) and charmm_box_dimension is None and box_no == 1:
        print_error_message = "ERROR: The {}-dimension for box {} was not provided.  " \
                              "The {}-dimension for box {}, " \
                              " is required to be an an integer, float, and be > 0.".format(dimension, box_no,
                                                                                            dimension, box_no)
        raise ValueError(print_error_message)

    return None

# user callable function to write the GOMC control file
def write_gomc_control_file(charmm_object, conf_filename,  ensemble_type,
                            RunSteps, Temperature, input_variables_dict=None):
    """
    The usable command that creates the GOMCControl object and write
    the GOMC control file via the GOMCControl.write_conf_file function

    Constructs the GOMC GOMCControl object with the defaults,
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
        charmm_object :  Charmm object
            Charmm object is has been parameterized from the selected force field.,
        ensemble_typ : str, ['NVT', 'NPT', 'GEMC_NPT', 'GCMC-NVT', 'GCMC']
            The ensemble type of the simulation.
        RunSteps : int (>0), must be an integer greater than zero.
            Sets the total number of simulation steps.
        Temperature : float or int (>0), must be an integer greater than zero.
            Temperature of system in Kelvin (K)
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

            # *******************************************************************
            # input_variables_dict options (keys and values) - (start)
            # Note: the input_variables_dict keys are also attributes
            # *******************************************************************
            Restart : boolean, default = False
                Determines whether to restart the simulation from restart file
                (*_restart.pdb and *_restart.psf) or not.
            RestartCheckpoint : boolean, default = False, default = "RANDOM"
                Determines whether to restart the simulation with the checkpoint
                file (checkpoint.dat) or not. Restarting the simulation with checkpoint.dat
                would result in an identical outcome, as if previous simulation was continued.
            PRNG : string or int (>= 0) ("RANDOM" or int)
                PRNG = Pseudo-Random Number Generator (PRNG). There are two (2) options, entering
                the string, "RANDOM", or a integer.
                --- "RANDOM", which selects a random seed number.  This will enter the line
                    "PRNG RANDOM" in the gomc configuration file.
                --- integer, which defines the integer seed number for the simulation. This is
                    equivalent to entering the following two lines in the configuration file:
                    line 1 = PRNG INTSEED
                    line 2 = Random_Seed user_selected_integer.
                Example 1: for a random seed enter the string "RANDOM.
                Example 2: for a specific seed number enter a integer of your choosing.
            ParaTypeCHARMM : boolean, default = True
                True if a CHARMM forcefield, False otherwise.
            ParaTypeMie : boolean, default = False
                True if a Mie forcefield type, False otherwise.
            ParaTypeMARTINI : boolean, default = False
                True if a MARTINI forcefield, False otherwise.
            RcutCoulomb_box_0 : int or float (>= 0), default = None
                Sets a specific radius in box 0 where the short-range electrostatic
                energy will be calculated (i.e., The distance to truncate the
                short-range electrostatic energy in box 0.)
                Note: if None, GOMC will default to the Rcut value
            RcutCoulomb_box_1 : int or float (>= 0), default = None
                Sets a specific radius in box 1 where the short-range electrostatic
                energy will be calculated (i.e., The distance to truncate the
                short-range electrostatic energy in box 0.)
                Note: if None, GOMC will default to the Rcut value
            Pressure : int or float (>= 0), default = 1.01325
                The pressure in bar utilized for the NPT and GEMC_NPT simulations.'
            Rcut : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 10
                Sets a specific radius in Angstroms that non-bonded interaction
                energy and force will be considered and calculated using defined potential function.
                The distance in Angstoms to truncate the LJ, Mie, or other VDW type potential at.
                Note: Rswitch is only used when the "Potential" = SWITCH.
            RcutLow : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 1
                Sets a specific minimum possible distance in Angstroms that reject
                any move that places any atom closer than specified distance.
                The minimum possible distance between any atoms.
                Sets a specific radius in Angstroms that non-bonded interaction
                Note: Rswitch is only used when the "Potential" = SWITCH.
            LRC : boolean, default = True
                If True, the simulation considers the long range tail corrections for the
                non-bonded VDW or dispersion interactions.
                Note: In case of using SHIFT or SWITCH potential functions, LRC will be ignored.
            Exclude : str ["1-2", "1-3", or "1-4"], default = 1-3"
                Note: In CHARMM force field, the 1-4 interaction needs to be considered.
                Choosing "Excude 1-3", will modify 1-4 interaction based on 1-4 parameters
                in parameter file. If a kind force field is used, where 1-4 interaction
                needs to be ignored, such as TraPPE, either Exclude "1-4" needs to be
                chosen or 1-4 parameter needs to be assigned to zero in the parameter file.
                --- "1-2": All interaction pairs of bonded atoms, except the ones that
                    separated with one bond, will be considered and modified using 1-4
                    parameters defined in parameter file.
                --- "1-3": All interaction pairs of bonded atoms, except the ones that
                    separated with one or two bonds, will be considered and modified using
                    1-4 parameters defined in parameter file.
                --- "1-4": All interaction pairs of bonded atoms, except the ones that
                    separated with one, two or three bonds, will be considered using
                    non-bonded parameters defined in parameter file.
            Potential : str, ["VDW", "EXP6", "SHIFT" or "SWITCH"], default = "VDW"
                Defines the potential function type to calculate non-bonded dispersion
                interaction energy and force between atoms.
                ---    "VDW":   Non-bonded dispersion interaction energy and force
                                calculated based on n-6 (Lennard - Jones) equation. This
                                function will be discussed further in the Intermolecular energy
                                and Virial calculation section.
                ---   "EXP6":   Non-bonded dispersion interaction energy and force calculated
                                based on exp-6 (Buckingham potential) equation.
                ---  "SHIFT":   This option forces the potential energy to be zero at Rcut distance.
                --- "SWITCH":   This option smoothly forces the potential energy to be zero at
                                Rcut distance and starts modifying the potential at Rswitch
                                distance. Depending on force field type, specific potential
                                function will be applied.
            Rswitch : int or float (>= 0 and RcutLow < Rswitch < Rcut), default = 9
                Note: Rswitch is only used when the SWITCH function is used
                (i.e., "Potential" = SWITCH). The Rswitch distance is in Angstrom. If the
                “SWITCH” function is chosen, Rswitch needs to be defined, otherwise, the
                program will be terminated. When using choosing "SWITCH" as potential function,
                the Rswitch distance defines where the non-bonded interaction energy
                modification is started, which is eventually truncated smoothly at Rcut
                distance.
            ElectroStatic : boolean, default = True
                Considers the coulomb interactions or not. If True, coulomb interactions are
                considered and false if not. Note: To simulate the polar molecule in MARTINI
                force field, ElectroStatic needs to be turn on (i.e., True). The MARTINI force
                field uses short-range coulomb interaction with constant Dielectric of 15.0.
            Ewald : boolean, default = True
                Considers the standard Ewald summation method for electrostatic calculations.
                If True, Ewald summation calculation needs to be considered and false if not.
                Note: By default, GOMC will set ElectroStatic to True if Ewald summation
                method was used to calculate coulomb interaction.
            CachedFourier : boolean, default = False
                Considers storing the reciprocal terms for Ewald summation calculation in
                order to improve the code performance. This option would increase the code
                performance with the cost of memory usage. If True, to store reciprocal
                terms of Ewald summation calculation and False if not.
                Warning: Monte Carlo moves, such as MEMC-1, MEMC-2, MEMC-3,
                IntraMEMC-1, IntraMEMC-2, and IntraMEMC-3 are not support with CachedFourier.
            Tolerance : float (0.0 < float < 1.0), default = 1e-05
                Sets the accuracy in Ewald summation calculation. Ewald separation parameter
                and number of reciprocal vectors for the Ewald summation are determined
                based on the accuracy parameter.
            Dielectric : int or float (>= 0.0), default = 15
                Sets dielectric value used in coulomb interaction when the Martini
                force field is used. Note: In MARTINI force field, Dielectric needs to
                be set to 15.0.
            PressureCalc : list [bool , int (> 0)] or [bool , step_frequency],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Calculate the system pressure or not. bool = True, enables the pressure calculation
                during the simulation, false disables the calculation. The int/step frequency sets the
                frequency of calculating the pressure.
            EqSteps : int (> 0), default = set via formula based on the number of RunSteps or 1M max
                Sets the number of steps necessary to equilibrate the system.
                Averaging will begin at this step.
                Note: In GCMC simulations, the Histogram files will be outputed at EqSteps.
            AdjSteps : int (> 0), default = set via formula based on the number of RunSteps or 1k max
                Sets the number of steps per adjustment of the parameter associated with each move
                (e.g. maximum translate distance, maximum rotation, maximum volume exchange, etc.).
            VDWGeometricSigma: boolean, default = False
                Use geometric mean, as required by OPLS force field, to combining
                Lennard-Jones sigma parameters for different atom types.
                If set to True, GOMC uses geometric mean to combine Lennard-Jones or VDW sigmas.
                Note: The default setting of VDWGeometricSigma is false, which uses the arithmetic
                mean when combining Lennard-Jones or VDW sigma parameters for different atom types.
            useConstantArea : boolean,  default = False
                Changes the volume of the simulation box by fixing the cross-sectional
                area (x-y plane). If True, the volume will change only in z axis,
                If False, the volume of the box will change in a way to maintain the constant
                axis ratio.
            FixVolBox0 : boolean, default = False
                Changing the volume of fluid phase (Box 1) to maintain the constant imposed
                pressure and Temperature, while keeping the volume of adsorbed phase (Box 0) fixed.
                Note: By default, GOMC will set useConstantArea to False if no value was set.
                It means, the volume of the box will change in a way to maintain the constant
                axis ratio.
            ChemPot : dict {str (4 dig limit) , int or float}, default = None
                The chemical potentials in GOMC units of energy, K.
                There is a 4 character limit for the string/residue name since the PDB/PSF
                files have a 4 character limitation and require and exact match in the conf file.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                The name of the residues and their corresponding chemical potential must specified
                for every residue in the system (i.e., {"residue_name" : chemical_potential}).
                Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY
                OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER.
                Example 1 (system with only water):  {"H2O" : -4000} .
                Example 2 (system with water and ethanol):  {"H2O" : -4000, "ETH" : -8000}
            Fugacity : dict {str , int or float (>= 0)}, default = None
                The fugacity in GOMC units of pressure, bar.
                There is a 4 character limit for the string/residue name since the PDB/PSF
                files have a 4 character limitation and require and exact match in the conf file.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                The name of the residues and their corresponding fugacity must specified
                for every residue in the system (i.e., {"residue_name" : fugacity}).
                Note: IF 2 KEYS WITH THE SAME STRING/RESIDUE ARE PROVIDED, ONE WILL BE AUTOMATICALLY
                OVERWRITTEN AND NO ERROR WILL BE THROWN IN THIS CONTROL FILE WRITER.
                Example 1 (system with only water):  {"H2O" : 1} .
                Example 2 (system with water and ethanol):  {"H2O" : 0.5, "ETH" : 10},
            CBMC_First : int (>= 0), default = 12
                The number of CD-CBMC trials to choose the first atom position
                (Lennard-Jones trials for first seed growth).
            CBMC_Nth : int (>= 0), default = 10
                The Number of CD-CBMC trials to choose the later atom positions
                (Lennard-Jones trials for first seed growth).
            CBMC_Ang : int (>= 0), default = 50
                The Number of CD-CBMC bending angle trials to perform for geometry
                (per the coupled-decoupled CBMC scheme).
            CBMC_Dih : int (>= 0), default = 50
                The Number of CD-CBMC dihedral angle trials to perform for geometry
                (per the coupled-decoupled CBMC scheme).
            OutputName : str (NO SPACES), , default = "Output_data", default = [True, 1M] or
                [True , set via formula based on the number of RunSteps or 1M max]
                The UNIQUE STRING NAME, WITH NO SPACES, which is used for the
                output block average, PDB, and PSF file names.
            CoordinatesFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or M max]
                Controls output of PDB file (coordinates). If bool is True, this
                enables outputting the coordinate files at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting
                the coordinates.
            RestartFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or 1M max]
                This creates the PDB and PSF (coordinate and topology) files for
                restarting the system at the set steps_per_data_output_int (frequency)
                If bool is True, this enables outputting the PDB/PSF restart files at the
                integer frequency (set steps_per_data_output_int), while “false”
                disables outputting the PDB/PSF restart files.
            CheckpointFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 1M] or [True , set via formula based on the number of RunSteps or 1M max]
                Controls the output of the last state of simulation at a specified step,
                in a binary file format (checkpoint.dat). Checkpoint file contains the
                following information in full precision:
                    (1) Last simulation step that saved into checkpoint file
                    (2) Simulation cell dimensions and angles
                    (3) Maximum amount of displacement (Å), rotation (δ), and volume (Å^3)
                        that is used in the Displacement, Rotation, MultiParticle, and Volume moves
                    (4) Number of Monte Carlo move trial and acceptance
                    (5) All molecule’s coordinates
                    (6) Random number sequence
                If bool is True, this enables outputing the checkpoint file at the
                integer frequency (set steps_per_data_ouput_int),
                while "False" disables outputting the checkpoint file.'
            ConsoleFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the output to the "console” or log file, which prints the
                acceptance statistics, and run timing info. In addition, instantaneously-selected
                thermodynamic properties will be output to this file.  If bool is True,
                this enables outputting the console data at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting the console
                data file.
            BlockAverageFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the block averages output of selected thermodynamic properties.
                Block averages are averages of thermodynamic values of interest for chunks of the
                simulation (for post-processing of averages or std. dev. in those values).
                If bool is True, this enables outputting the block averaging data/file at the
                integer frequency (set steps_per_data_output_int),  while "False"
                disables outputting the block averaging data/file.
            HistogramFreq : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = [True, 10k] or [True , set via formula based on the number of RunSteps or 10k max]
                Controls the histograms. Histograms are a binned listing of observation frequency
                for a specific thermodynamic variable. In the GOMC code, they also control the output
                of a file containing energy/molecule samples, which is only used for the "GCMC"
                ensemble simulations for histogram reweighting purposes. If bool is True, this
                enables outputting the data to the histogram data at the integer frequency
                (set steps_per_data_output_int), while "False" disables outputting the histogram
                data.
            DistName : str (NO SPACES), default = "dis"
                Short phrase which will be combined with RunNumber and RunLetter
                to use in the name of the binned histogram for molecule distribution.
            HistName : str (NO SPACES), default = "his"
                Short phrase, which will be combined with RunNumber and RunLetter,
                to use in the name of the energy/molecule count sample file.
            RunNumber : int  ( > 0 ), default = 1
                 Sets a number, which is a part of DistName and HistName file name.
            RunLetter : str (1 alphabetic character only), default = "a"
                Sets a letter, which is a part of DistName and HistName file name.
            SampleFreq : int ( > 0 ), default = 500
                The number of steps per histogram sample or frequency.
            OutEnergy : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the energy data into the block averages and console output/log
            OutPressure : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the pressure data into the block averages and console output/log files.
            OutMolNumber : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the number of molecules data into the block averages and console
                output/log files.
            OutDensity : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the density data into the block averages and console output/log files.
            OutVolume : [bool, bool], default = [True, True]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the volume data into the block averages and console output/log files.
            OutSurfaceTension : [bool, bool], default = [False, False]
                The list provides the booleans to [block_averages_bool, console_output_bool].
                This outputs the surface tension data into the block averages and console
                output/log files.
            FreeEnergyCalc : list [bool , int (> 0)] or [Generate_data_bool , steps_per_data_output_int],
                default = None
                bool = True enabling free energy calculation during the simulation, false disables
                the calculation. The int/step frequency sets the frequency of calculating the free energy.
            MoleculeType : list [str , int (> 0)] or ["residue_name" , residue_ID], default = None
                The user must set this variable as there is no working default.
                Note: ONLY 4 characters can be used for the string (i.e., "residue_name").
                Sets the solute molecule kind (residue name) and molecule number (residue ID),
                which absolute solvation free will be calculated for.'
            InitialState : int (>= 0), default = None
                The user must set this variable as there is no working default.
                The index of LambdaCoulomb and LambdaVDW vectors. Sets the index of the
                LambdaCoulomb and LambdaVDW vectors, to determine the simulation lambda value for
                VDW and Coulomb interactions.
                WARNING : This must an integer within the vector count of the LambdaVDW and LambdaCoulomb,
                in which the counting starts at 0.  '
            LambdaVDW : list of floats (0 <= floats <= 1), default = None
                The user must set this variable as there is no working default (default = {}).
                Lambda values for VDW interaction in ascending order. Sets the intermediate
                lambda states to which solute-solvent VDW interactions are scaled.
                WARNING : This list must be the same length as the "LambdaCoulomb" list length.
                WARNING : All lambda values must be stated in the ascending order, otherwise the
                program will terminate.
                Example of ascending order 1: [0.1, 1.0,]
                Example of ascending orde 2: [0.1, 0.2, 0.4, 0.9]
            LambdaCoulomb : list of floats (0 <= floats <= 1), default = None
                Lambda values for Coulombic interaction in ascending order. Sets the intermediate
                lambda states to which solute-solvent Coulombic interactions are scaled.
                GOMC defauts to the "LambdaVDW" values for the Coulombic interaction
                if no "LambdaCoulomb" variable is set.
                WARNING : This list must be the same length as the "LambdaVDW" list length.
                WARNING : All lambda values must be stated in the ascending order, otherwise
                the program will terminate.
                Example of ascending order 1: [0.1, 1.0,]
                Example of ascending order 2: [0.1, 0.2, 0.4, 0.9] '
            ScaleCoulomb : bool, default = False
                Determines to scale the Coulombic interaction non-linearly
                (soft-core scheme) or not.
                True if the Coulombic interaction needs to be scaled non-linearly.
                False if the Coulombic interaction needs to be scaled linearly.
            ScalePower : int (>= 0), default = 2
                The p value in the soft-core scaling scheme, where the distance
                between solute and solvent is scaled non-linearly.
            ScaleAlpha : int or float (>= 0), default = 0.5
                The alpha value in the soft-core scaling scheme, where the distance
                between solute and solvent is scaled non-linearly.
            MinSigma : int or float (>= 0), default = 3
                The minimum sigma value in the soft-core scaling scheme, where the
                distance between solute and solvent is scaled non-linearly.
            DisFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.19, 'GCMC': 0.15}
                Fractional percentage at which the displacement move will occur
                (i.e., fraction of displacement moves).
            RotFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.15, 'NPT': 0.15, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.15}
                Fractional percentage at which the rotation move will occur.
                (i.e., fraction of rotation moves).
            IntraSwapFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.3, 'NPT': 0.29, 'GEMC_NVT': 0.1, 'GEMC_NPT': 0.1, 'GCMC': 0.1}
                Fractional percentage at which the molecule will be removed from a
                box and inserted into the same box using coupled-decoupled configurational-bias
                algorithm. (i.e., fraction of intra-molecule swap moves).
            SwapFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.35}
                For Gibbs and Grand Canonical (GC) ensemble runs only: Fractional
                percentage at which molecule swap move will occur using coupled-decoupled
                configurational-bias. (i.e., fraction of molecule swaps moves).
            RegrowthFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.3, 'NPT': 0.3, 'GEMC_NVT': 0.2, 'GEMC_NPT': 0.2, 'GCMC': 0.15}
                Fractional percentage at which part of the molecule will be deleted and
                then regrown using coupled- decoupled configurational-bias algorithm
                (i.e., fraction of molecular growth moves).
            CrankShaftFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.1, 'NPT': 0.1, 'GEMC_NVT': 0.1, 'GEMC_NPT': 0.1, 'GCMC': 0.1}
                Fractional percentage at which crankshaft move will occur.
                In this move, two atoms that are forming angle or dihedral are selected
                randomly and form a shaft. Then any atoms or group that are within these
                two selected atoms, will rotate around the shaft to sample intra-molecular
                degree of freedom (i.e., fraction of crankshaft moves).
            VolFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.01, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.01, 'GCMC': 0.0}
                For isobaric-isothermal (NPT) ensemble and Gibbs ensemble
                (GEMC_NPT and GEMC_NVT) runs only: Fractional percentage at
                which a volume move will occur (i.e., fraction of Volume moves).
            MultiParticleFreq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which multi-particle move will occur.
                In this move, all molecules in the selected simulation box will be rigidly
                rotated or displaced simultaneously, along the calculated torque or force
                respectively (i.e., fraction of multi-particle moves).
            IntraMEMC_1Freq : int or float (0 <= value <= 1), default are specific for each ensemble
            {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume within
                same simulation box.  This move need additional information such as
                ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.
            MEMC_1Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will
                be exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes.  This move need additional information such as
                ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind, and ExchangeLargeKind.
            IntraMEMC_2Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind
                will be exchanged with a specified large molecule kind in defined sub-volume
                within same simulation box. Backbone of small and large molecule kind will be
                used to insert the large molecule more efficiently. This move need additional
                information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            MEMC_2Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes. Backbone of small and large molecule kind will be
                used to insert the large molecule more efficiently. This move need additional
                information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, SmallKindBackBone, and LargeKindBackBone. '
            IntraMEMC_3Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume within same
                simulation box. Specified atom of the large molecule kind will be used to insert
                the large molecule using coupled-decoupled configurational-bias. This move need
                additional information such as ExchangeVolumeDim, ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, and LargeKindBackBone. '
            MEMC_3Freq : int or float (0 <= value <= 1), default are specific for each ensemble
                {'NVT': 0.0, 'NPT': 0.0, 'GEMC_NVT': 0.0, 'GEMC_NPT': 0.0, 'GCMC': 0.0}
                Fractional percentage at which specified number of small molecule kind will be
                exchanged with a specified large molecule kind in defined sub-volume,
                between simulation boxes.  Specified atom of the large molecule kind will be
                used to insert the large molecule using coupled-decoupled configurational-bias.
                This move need additional information such as ExchangeVolumeDim,
                ExchangeRatio, ExchangeSmallKind, ExchangeLargeKind, and LargeKindBackBone.
            ExchangeVolumeDim : list of 3 floats or integers or [X-dimension, Y-dimension, Z-dimension)],
                default = [1.0, 1.0, 1.0]
                To use all variations of MEMC and Intra-MEMC Monte Carlo moves, the exchange
                subvolume must be defined. The exchange sub-volume is defined as an orthogonal box
                with x, y, and z-dimensions, where small molecule/molecules kind will be selected
                from to be exchanged with a large molecule kind.
                Note: Currently, the X and Y dimension cannot be set independently (X = Y = max(X, Y)).
                Note: A heuristic for setting good values of the x, y, and z-dimensions is to use
                the geometric size of the large molecule plus 1-2 Å in each dimension.
                Note: In case of exchanging 1 small molecule kind with 1 large molecule kind in
                IntraMEMC-2, IntraMEMC-3, MEMC-2, MEMC-3 Monte Carlo moves, the sub-volume
                dimension has no effect on acceptance rate. '
            MEMC_DataInput : nested lists, default = None
                Enter data as a list with some sub-lists as follows:
                [[ExchangeRatio_int (> 0), ExchangeLargeKind_str,
                [LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ],
                ExchangeSmallKind_str, [SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ]],
                ...,
                [ExchangeRatio_int (> 0), ExchangeLargeKind_str,
                [LargeKindBackBone_atom_1_str_or_NONE, LargeKindBackBone_atom_2_str_or_NONE ],
                ExchangeSmallKind_str, [SmallKindBackBone_atom_1_str_or_NONE, SmallKindBackBone_atom_2_str_or_NONE ].
                NOTE: CURRENTLY ALL THESE INPUTS NEED TO BE SPECIFIED, REGARDLESS OF THE MEMC TYPE
                SELECTION. IF THE SmallKindBackBone or LargeKindBackBone IS NOT REQUIRED FOR THE
                MEMC TYPE, None CAN BE USED IN PLACE OF A STRING.
                Note: These strings must match the residue in the psf and psb files or it will fail.
                It is recommended that the user print the Charmm object psf and pdb files
                and review the residue names that match the atom name before using the in
                the MEMC_DataInput variable input.
                Note: see the below data explanations for the ExchangeRatio, ExchangeSmallKind,
                ExchangeLargeKind, LargeKindBackBone, SmallKindBackBone.
                Example 1 (MEMC-1) : [ [1, 'WAT', [None, None], 'wat', [None, None]] ,
                [1, 'WAT', [None, None], 'wat', [None, None]] .
                Example 2 (MEMC-2): [ [1, 'WAT', ['O1', 'H1'], 'wat', ['O1', 'H1' ]] ,
                [1, 'WAT', ['H1', 'H2'], 'wat', ['H1', 'H2' ]] .
                Example 3 (MEMC-3) : [ [2, 'WAT', 'O1', 'H1'], 'wat', [None, None]] ,
                [2, 'WAT', ['H1', 'H2'], 'wat', [None, None]] .\n"
                --- ExchangeRatio     = MEMC parameters (all ensembles): int (> 0), default = None
                                        The Ratio of exchanging small molecule/molecules with 1 large molecule.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the exchange ratio must be defined. The exchange ratio defines how
                                        many small molecule will be exchanged with 1 large molecule. For each
                                        large-small molecule pairs, one exchange ratio must be defined.
                --- ExchangeSmallKind = MEMC parameters (all ensembles):  str, default = None
                                        The small molecule kind (resname) to be exchanged.
                                        Note: ONLY 4 characters can be used for the strings.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the small molecule kind to be exchanged with a large molecule
                                        'kind must be defined. Multiple small molecule kind can be specified.
                --- ExchangeLargeKind = MEMC parameters (all ensembles):  str, default = None
                                        The large molecule kind (resname) to be exchanged.
                                        Note: ONLY 4 characters can be used for the strings.
                                        To use all variation of MEMC and Intra-MEMC Monte Carlo moves,
                                        the large molecule kind to be exchanged with small molecule '
                                        kind must be defined. Multiple large molecule kind can be specified.
                --- LargeKindBackBone = MEMC parameters (all ensembles): list [str, str] or [None, None], default = None
                                        Note: ONLY 4 characters can be used for the strings.
                                        The [None, None] values can only be used if that MEMC type does
                                        not require them.  The strings for the the atom name 1 and atom name 2
                                        that belong to the large molecule’s backbone
                                        (i.e., [str_for_atom_name_1, str_for_atom_name_2])
                                        To use MEMC-2, MEMC-3, IntraMEMC-2, and IntraMEMC-3 Monte Carlo moves, the
                                        large molecule backbone must be defined. The backbone of the molecule is defined
                                        as a vector that connects two atoms belong to the large molecule. The large
                                        molecule backbone will be used to align the sub-volume in MEMC-2 and IntraMEMC-2
                                        moves, while in MEMC-3 and IntraMEMC-3 moves, it uses the atom name to start
                                        growing the large molecule using coupled-decoupled configurational-bias. For
                                        each large-small molecule pairs, two atom names must be defined.
                                        Note: all atom names in the molecule must be unique.
                                        Note: In MEMC-3 and IntraMEMC-3 Monte Carlo moves, both atom names must be same,
                                        otherwise program will be terminated.
                                        Note: If the large molecule has only one atom (mono atomic molecules),
                                        same atom name must be used for str_for_atom_name_1 and str_for_atom_name_2
                                        of the LargeKindBackBone.
                --- SmallKindBackBone = MEMC parameters (all ensembles): list [str, str] or [None, None], default = None
                                        Note: ONLY 4 characters can be used for the strings.
                                        The [None, None] values can only be used if that MEMC type does not
                                        require them. The strings for the the atom name 1 and atom name 2 that
                                        belong to the small molecule’s backbone
                                        (i.e., [str_for_atom_name_1, str_for_atom_name_2]) '
                                        To use MEMC-2, and IntraMEMC-2 Monte Carlo moves, the small molecule backbone
                                        must be defined. The backbone of the molecule is defined as a vector that
                                        connects two atoms belong to the small molecule and will be used to align the
                                        sub-volume. For each large-small molecule pairs, two atom names must be defined.
                                        Note: all atom names in the molecule must be unique.
                                        Note: If the small molecule has only one atom (mono atomic molecules), same atom
                                        name must be used str_for_atom_name_1 and str_for_atom_name_2
                                        of the SmallKindBackBone.
            # *******************************************************************
            # input_variables_dict options (keys and values) - (end)
            # Note: the input_variables_dict keys are also attributes
            # *******************************************************************

    Notes
    -------
    The user input variables (input_variables_dict) and the specific
    ensembles they are also available with can be accessed by the running
    print_valid_ensemble_input_variables('NPT', description = True)
    command, as the information is dynamically contained here.

    The details of the required inputs for the selected
    ensembles can be found by running this python workbook,

    >>> print_valid_required_input_variables('NVT', description = True)

    which prints the required inputs with their subsection description
    for the selected 'NVT' ensemble (other ensembles can be set as well).

    The details of the required inputs for the selected
    ensembles can be found by the following function,
    >>> print_valid_required_input_variables('NVT', description = True)
    which prints the required inputs with their subsection description
    for the selected 'NVT' ensemble (other ensembles can be set as well).
    The box units imported are in nm (standard MoSDeF units).
    The units for this writer are auto-scaled to Angstroms, so they
    can be directly used in the GOMC or NAMD engines.

    Note: all of the move types are not available in for every ensemble.
    Note: all of the move fractions must sum to 1, or the control file
    writer will fail.

    The input variables (input_variables_dict) and text extracted with permission from
    the GOMC manual version 2.60. Some of the text was modified from its original version.
    Cite: Potoff, Jeffrey; Schwiebert, Loren; et. al. GOMC Documentation.
    https://raw.githubusercontent.com/GOMC-WSU/GOMC/master/GOMC_Manual.pdf, 2021.

    Returns
    -------
    If completed without errors: str
        returns "GOMC_CONTROL_FILE_WRITTEN" when the GOMC input control file is writen
    If completed with errors:  None
    """

    #write the control file and return a testable value
    gomc_control = GOMCControl(charmm_object, ensemble_type,
                               RunSteps, Temperature, input_variables_dict=input_variables_dict)
    test_gomc_control_write_conf_file = gomc_control.write_conf_file(conf_filename)

    if gomc_control.input_error is False and test_gomc_control_write_conf_file == "GOMC_CONTROL_FILE_WRITTEN":

        return "GOMC_CONTROL_FILE_WRITTEN"
    else:
        return None
