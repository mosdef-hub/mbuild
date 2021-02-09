import pytest
import mbuild as mb

from mbuild.tests.base_test import BaseTest
from mbuild.formats import charmm_writer
from mbuild.formats.charmm_writer import Charmm
from mbuild.utils.io import has_foyer
from mbuild.utils.conversion import base10_to_base16_alph_num
from mbuild.utils.conversion import base10_to_base26_alph_num
from mbuild.utils.conversion import base10_to_base52_alph_num
from mbuild.utils.conversion import base10_to_base62_alph_num
from mbuild.utils.specific_ff_to_residue  import specific_ff_to_residue
from foyer.forcefields import forcefields
from collections import OrderedDict


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestCharmmWriterData(BaseTest):

    def test_save(self, EthaneGOMC):
        Charmm(EthaneGOMC, 'ethane', FF_filename='ethane',
                                   residues = [EthaneGOMC.name], forcefield_selection = 'oplsaa')

    def test_save_charmm_GOMC_FF(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'charmm_data', FF_filename='charmm_data',
                          residues = [EthaneGOMC.name], forcefield_selection = 'oplsaa')
        charmm.write_inp()

        out_GOMC = open('charmm_data.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName' \
               ' (i.e., atoms_type_per_utilized_FF)' in line:
                assert len(out_GOMC[i + 1].split('!')[0].split()) == 3
                assert out_GOMC[i + 1].split('!')[0].split()[0:3] == ['*', 'A', '12.010780']
                assert len(out_GOMC[i + 2].split('!')[0].split()) == 3
                assert out_GOMC[i + 2].split('!')[0].split()[0:3] == ['*', 'B', '1.007947']
                assert out_GOMC[i + 1].split()[4:5] == ['opls_135_ETH']
                assert out_GOMC[i + 2].split()[4:5] == ['opls_140_ETH']

            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                bond_types = [['A', 'B', '340.0', '1.09'], ['A', 'A', '268.0', '1.529']]
                assert len(out_GOMC[i + 1].split('!')[0].split()) == 4
                assert len(out_GOMC[i + 2].split('!')[0].split()) == 4
                if out_GOMC[i + 1].split('!')[0].split()[0:4] == bond_types[0]:
                    assert out_GOMC[i + 1].split('!')[0].split()[0:4] == bond_types[0]
                    assert out_GOMC[i + 2].split('!')[0].split()[0:4] == bond_types[1]
                else:
                    out_GOMC[i + 1].split('!')[0].split()[0:4] == bond_types[1]
                    assert out_GOMC[i + 1].split('!')[0].split()[0:4] == bond_types[1]
                    assert out_GOMC[i + 2].split('!')[0].split()[0:4] == bond_types[0]

            elif '!atom_types 		Ktheta	Theta0			  atoms_types_per_utilized_FF' in line:
                angle_types = [['A', 'A', 'B', '37.5', '110.70000'], ['B', 'A', 'B', '33.0', '107.80000']]
                assert len(out_GOMC[i + 1].split('!')[0].split()) == 5
                assert len(out_GOMC[i + 2].split('!')[0].split()) == 5
                if out_GOMC[i + 1].split('!')[0].split()[0:5] == angle_types[0]:
                    assert out_GOMC[i + 1].split('!')[0].split()[0:5] == angle_types[0]
                    assert out_GOMC[i + 2].split('!')[0].split()[0:5] == angle_types[1]
                else:
                    out_GOMC[i + 1].split('!')[0].split()[0:4] == angle_types[1]
                    assert out_GOMC[i + 1].split('!')[0].split()[0:5] == angle_types[1]
                    assert out_GOMC[i + 2].split('!')[0].split()[0:5] == angle_types[0]

            elif '!atom_types 			Kchi		n	delta		  atoms_types_per_utilized_FF' in line:
                dihed_types = [['B', 'A', 'A', 'B', '0.300000', '0', '90.0'],
                               ['B', 'A', 'A', 'B', '0.000000', '1', '180.0'],
                               ['B', 'A', 'A', 'B', '0.000000', '2', '0.0'],
                               ['B', 'A', 'A', 'B', '-0.150000', '3', '180.0'],
                               ['B', 'A', 'A', 'B', '0.000000', '4', '0.0'],
                               ['B', 'A', 'A', 'B', '0.000000', '5', '180.0']
                               ]
                for j in range(0, len(dihed_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] ==dihed_types[j]

            elif '!atype 	ignored	epsilon 	Rmin/2 		ignored	eps,1-4		Rmin/2,1-4' \
                 '		  atom_type_per_utilized_FF' in line:
                NB_types = [['A', '0.00', '-0.066000000', '1.96430858454', '0.00', '-0.033000000',	'1.96430858454'],
                            ['B', '0.00', '-0.030000000', '1.40307756039', '0.00', '-0.015000000',	'1.40307756039']]

                for j in range(0, len(NB_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] == NB_types[j]

            else:
                pass

    def test_save_charmm_psf(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'charmm_data', FF_filename='charmm_data',
                          residues = [EthaneGOMC.name], forcefield_selection = 'oplsaa')
        charmm.write_psf()

        out_GOMC = open('charmm_data.psf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '8 !NATOM' in line:
                Atom_type_charge_etc_list = [['1', 'SYS', '1', 'ETH', 'C1', 'A', '-0.180000', '12.0108'],
                                             ['2', 'SYS', '1', 'ETH', 'C2', 'A', '-0.180000', '12.0108'],
                                             ['3', 'SYS', '1', 'ETH', 'H1', 'B', '0.060000', '1.0079'],
                                             ['4', 'SYS', '1', 'ETH', 'H2', 'B', '0.060000', '1.0079'],
                                             ['5', 'SYS', '1', 'ETH', 'H3', 'B', '0.060000', '1.0079'],
                                             ['6', 'SYS', '1', 'ETH', 'H4', 'B', '0.060000', '1.0079'],
                                             ['7', 'SYS', '1', 'ETH', 'H5', 'B', '0.060000', '1.0079'],
                                             ['8', 'SYS', '1', 'ETH', 'H6', 'B', '0.060000', '1.0079']
                                             ]
                for j in range(0, len(Atom_type_charge_etc_list)):
                    assert out_GOMC[i + 1 + j ].split()[0:8] == Atom_type_charge_etc_list[j]

            else:
                pass


    def test_save_charmm_pdb(self, EthaneGOMC):
        charmm = Charmm(EthaneGOMC, 'charmm_data', FF_filename='charmm_data',
                          residues = [EthaneGOMC.name], forcefield_selection = 'oplsaa')
        charmm.write_pdb()

        out_GOMC = open('charmm_data.pdb', 'r').readlines()

        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', 'A', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', 'A', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', 'A', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', 'A', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', 'A', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', 'A', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', 'A', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', 'A', '1']
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00','C'], ['1.00', '0.00','C'], ['1.00', '0.00','H'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','H'], ['1.00', '0.00','H'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','H'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] == Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass


    def test_save_charmm_UA_GOMC_FF(self, TwoPropanolUA):
        charmm = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                          residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                          bead_to_atom_name_dict= {'_CH3' : 'C'})
        charmm.write_inp()

        out_GOMC = open('charmm_data_UA.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName ' \
               '(i.e., atoms_type_per_utilized_FF)' in line:
                atom_types_1 = [['*', 'A', '15.035000'], ['*', 'B', '13.019000'],
                              ['*', 'D', '15.999430'], ['*', 'C', '1.007947']]
                atom_types_2= [['CH3_sp3_POL'], ['CH_O_POL'], ['O_POL'],['H_POL']]

                for j in range(0, len(atom_types_1)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 3
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:3] == atom_types_1[j]
                    assert out_GOMC[i + 1 + j].split()[4:5] == atom_types_2[j]

            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                bond_types = [['C', 'D', '600.40152964', '0.945'], ['B', 'D', '600.40152964', '1.43'],
                              ['A', 'B', '600.40152964', '1.54']]
                total_bonds_evaluated = []
                total_bonds_evaluated_reorg = []
                for j in range(0, len(bond_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 4
                    if out_GOMC[i + 1 + j].split('!')[0].split()[0:4] == bond_types[0] or bond_types[1] or \
                            bond_types[2]:
                        total_bonds_evaluated.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                for k in range(0, len(bond_types)):
                    if bond_types[k] in total_bonds_evaluated:
                        total_bonds_evaluated_reorg.append(bond_types[k])
                assert total_bonds_evaluated_reorg == bond_types

            elif '!atom_types 		Ktheta	Theta0			  atoms_types_per_utilized_FF' in line:
                angle_types = [['A', 'B', 'A', '62.10013026', '112.00007'], ['A', 'B', 'D', '50.0775', '109.46989'],
                               ['B', 'D', 'C', '55.04555449', '108.49987']]
                total_angles_evaluated = []
                total_angles_evaluated_reorg = []
                for j in range(0, len(angle_types )):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 5
                    if out_GOMC[i + 1 + j].split('!')[0].split()[0:5]==angle_types[0] or angle_types[1] \
                            or angle_types[2]:
                        total_angles_evaluated.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:5])
                for k in range(0, len(angle_types)):
                    if angle_types[k] in  total_angles_evaluated:
                        total_angles_evaluated_reorg.append(angle_types[k])
                assert total_angles_evaluated_reorg == angle_types

            elif '!atom_types 			Kchi		n	delta		  atoms_types_per_utilized_FF' in line:
                dihedral_types = [['A', 'B', 'D', 'C', '0.647232', '0', '90.0'],
                                  ['A', 'B', 'D', 'C', '-0.392135', '1', '180.0'],
                                  ['A', 'B', 'D', 'C', '-0.062518', '2', '0.0'],
                                  ['A', 'B', 'D', 'C', '0.345615', '3', '180.0'],
                                  ['A', 'B', 'D', 'C', '0.000000', '4', '0.0'],
                                  ['A', 'B', 'D', 'C', '0.000000', '5', '180.0']
                                  ]
                for j in range(0, len(dihedral_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] == dihedral_types[j]

            elif '!atype 	ignored	epsilon 	Rmin/2 		ignored	eps,1-4		Rmin/2,1-4		  ' \
                 'atom_type_per_utilized_FF' in line:
                NB_types =    [['A', '0.00', '-0.194745937', '2.10461634058', '0.00', '-0.000000000', '2.10461634058'],
                               ['B', '0.00', '-0.019872012', '2.43013033459', '0.00', '-0.000000000', '2.43013033459'],
                               ['D', '0.00', '-0.184809990', '1.69491769295', '0.00', '-0.000000000', '1.69491769295'],
                               ['C', '0.00', '-0.000000000', '5.61231024155', '0.00', '-0.000000000', '5.61231024155']
                               ]

                for j in range(0, len(NB_types)):
                    assert len(out_GOMC[i + 1+ j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1+ j].split('!')[0].split()[0:7] == NB_types[j]

            else:
                pass


    def test_save_charmm_UA_psf(self, TwoPropanolUA):
        charmm = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                          residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                          bead_to_atom_name_dict= {'_CH3' : 'C'})
        charmm.write_psf()

        out_GOMC = open('charmm_data_UA.psf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '5 !NATOM' in line:
                Atom_type_charge_etc_list = [['1', 'SYS', '1', 'POL', 'C1', 'A', '0.000000', '15.0350'],
                                             ['2', 'SYS', '1', 'POL', 'BD1', 'B', '0.265000', '13.0190'],
                                             ['3', 'SYS', '1', 'POL', 'O1', 'D', '-0.700000', '15.9994'],
                                             ['4', 'SYS', '1', 'POL', 'H1', 'C', '0.435000', '1.0079'],
                                             ['5', 'SYS', '1', 'POL', 'C2', 'A', '0.000000', '15.0350'],
                                             ]

                for j in range(0, len(Atom_type_charge_etc_list)):
                    assert out_GOMC[i + 1 + j].split()[0:8] ==  Atom_type_charge_etc_list[j]

            else:
                pass


    def test_save_charmm_UA_pdb(self, TwoPropanolUA):
        charmm = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                          residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                          bead_to_atom_name_dict= {'_CH3' : 'C'})
        charmm.write_pdb()

        out_GOMC = open('charmm_data_UA.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'POL', 'A', '1'],
                                             ['ATOM', '2', 'BD1', 'POL', 'A', '1'],
                                             ['ATOM', '3', 'O1', 'POL', 'A', '1'],
                                             ['ATOM', '4', 'H1', 'POL', 'A', '1'],
                                             ['ATOM', '5', 'C2', 'POL', 'A', '1'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00', 'EP'], ['1.00', '0.00','EP'], ['1.00', '0.00','O'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','EP'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass


    def test_charmm_pdb_fix_angle_bond_fix_atoms(self, EthaneGOMC,EthanolGOMC):
        test_box_ethane_propane = mb.fill_box(compound=[EthaneGOMC,EthanolGOMC],
                                  n_compounds= [1,1] ,
                                  box=[2.0, 2.0, 2.0])
        charmm = Charmm(test_box_ethane_propane, 'Test_fixes_angle_bond_atoms',
                        FF_filename='Test_fixes_angle_bond_atoms',
                        residues=[EthanolGOMC.name, EthaneGOMC.name],
                        forcefield_selection='oplsaa',
                        fix_residue=[EthaneGOMC.name],
                        fix_residue_in_box=[EthanolGOMC.name],
                        fix_res_bonds_angles=[EthaneGOMC.name]
                        )
        charmm.write_inp()
        charmm.write_pdb()

        out_GOMC = open('Test_fixes_angle_bond_atoms.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName ' \
               '(i.e., atoms_type_per_utilized_FF)' in line:
                mass_type_1 = [ ['*', 'A', '12.010780'], ['*', 'C', '1.007947'], ['*', 'B', '12.010780'],
                                ['*', 'G', '12.010780'], ['*', 'E', '15.999430'], ['*', 'D', '1.007947'],
                                ['*', 'F', '1.007947']
                                ]
                mass_type_2 = [ ['opls_135_ETH'], ['opls_140_ETH'], ['opls_135_ETO'], ['opls_157_ETO'],
                                ['opls_154_ETO'], ['opls_140_ETO'], ['opls_155_ETO'] ]

                for j in range(0, len(mass_type_1)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 3
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:3] ==  mass_type_1[j]
                    assert out_GOMC[i + 1+ j].split()[4:5] == mass_type_2[j]


            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                bond_types = [['D', 'G', '340.0', '1.09'], ['E', 'G', '320.0', '1.41'],
                                    ['E', 'F', '553.0', '0.945'], ['A', 'C', '999999999999', '1.09'],
                                    ['B', 'D', '340.0', '1.09'], ['A', 'A', '999999999999', '1.529'],
                                    ['B', 'G', '268.0', '1.529']]
                total_bonds_evaluated = []
                total_fixed_bonds = []
                for j in range(0,  7):
                    total_bonds_evaluated.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                    if out_GOMC[i + 1 + j].split('!')[0].split()[2:3] == ['999999999999']:
                        total_fixed_bonds.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                assert total_bonds_evaluated.sort() == bond_types.sort()
                assert len(total_fixed_bonds) == 2

            elif '!atom_types 		Ktheta	Theta0			  atoms_types_per_utilized_FF' in line:
                fixed_angle_types = [['A',	'A',	'C',	'999999999999',	'110.70000'],
                                     ['C',	'A',	'C',	'999999999999',	'107.80000']
                                     ]
                total_angles_evaluated = []
                total_fixed_angles = []
                for j in range(0, 9):
                    if out_GOMC[i + 1 + j].split('!')[0].split()[0:4] == (fixed_angle_types[0] or fixed_angle_types[1]):
                        total_angles_evaluated.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                    if out_GOMC[i + 1 + j].split('!')[0].split()[3:4] == ['999999999999']:
                        total_fixed_angles.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                assert total_angles_evaluated.sort() == total_angles_evaluated.sort()
                assert len(total_fixed_angles) == len(fixed_angle_types)

            else:
                pass

        out_GOMC = open('Test_fixes_angle_bond_atoms.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] == ['CRYST1', '20.000', '20.000', '20.000',
                                                    '90.00', '90.00', '90.00']

        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', 'A', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', 'A', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', 'A', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', 'A', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', 'A', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', 'A', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', 'A', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', 'A', '1'],
                                             ['ATOM', '9', 'C1', 'ETO', 'A', '2'],
                                             ['ATOM', '10', 'C2', 'ETO', 'A', '2'],
                                             ['ATOM', '11', 'O1', 'ETO', 'A', '2'],
                                             ['ATOM', '12', 'H1', 'ETO', 'A', '2'],
                                             ['ATOM', '13', 'H2', 'ETO', 'A', '2'],
                                             ['ATOM', '14', 'H3', 'ETO', 'A', '2'],
                                             ['ATOM', '15', 'H4', 'ETO', 'A', '2'],
                                             ['ATOM', '16', 'H5', 'ETO', 'A', '2'],
                                             ['ATOM', '17', 'H6', 'ETO', 'A', '2'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '1.00', 'C'], ['1.00', '1.00', 'C'], ['1.00', '1.00', 'H'],
                                             ['1.00', '1.00', 'H'], ['1.00', '1.00', 'H'], ['1.00', '1.00', 'H'],
                                             ['1.00', '1.00', 'H'], ['1.00', '1.00', 'H'],
                                             ['1.00', '2.00', 'C'], ['1.00', '2.00', 'C'], ['1.00', '2.00', 'O'],
                                             ['1.00', '2.00', 'H'], ['1.00', '2.00', 'H'], ['1.00', '2.00', 'H'],
                                             ['1.00', '2.00', 'H'], ['1.00', '2.00', 'H'], ['1.00', '2.00', 'H']]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass

    def test_charmm_pdb_residue_reorder_box_sizing(self, EthanolGOMC, EthaneGOMC):
        test_box_ethane_EthanolGOMC = mb.fill_box(compound=[EthanolGOMC, EthaneGOMC],
                                  n_compounds= [1,1] ,
                                  box=[2.0, 2.0, 2.0])
        charmm = Charmm(test_box_ethane_EthanolGOMC, 'residue_reorder_box_sizing_box_0',
                          structure_box_1 = EthaneGOMC,
                          filename_box_1 = 'residue_reorder_box_sizing_box_1',
                          FF_filename='residue_reorder_box',
                          residues = [EthaneGOMC.name, EthanolGOMC.name],
                          forcefield_selection = { EthanolGOMC.name : 'oplsaa', EthaneGOMC.name : 'oplsaa'},
                          fix_residue=None,
                          fix_residue_in_box = None,
                          fix_res_bonds_angles=None,
                          reorder_res_in_pdb_psf=True,
                          box_0 = [3, 3, 3],
                          box_1 =[4, 4, 4],
                          bead_to_atom_name_dict={'_CH3': 'C'}
                          )

        charmm.write_pdb()


    def test_charmm_pdb_no_differenc_1_4_coul_scalars(self, TwoPropanolUA, EthaneGOMC):
        test_box_ethane_TwoPropanolUA = mb.fill_box(compound=[TwoPropanolUA, EthaneGOMC],
                                                    n_compounds=[1, 1],
                                                    box=[2.0, 2.0, 2.0])
        try:
            Test_value = Charmm(test_box_ethane_TwoPropanolUA, 'residue_reorder_box_sizing_box_0',
                            structure_box_1=EthaneGOMC,
                            filename_box_1='residue_reorder_box_sizing_box_1',
                            FF_filename='residue_reorder_box',
                            residues=[TwoPropanolUA.name, EthaneGOMC.name],
                            forcefield_selection= {TwoPropanolUA.name : 'trappe-ua', EthaneGOMC.name : 'oplsaa'} ,
                            fix_residue=None,
                            fix_residue_in_box=None,
                            fix_res_bonds_angles=None,
                            reorder_res_in_pdb_psf=False,
                            box_0=[3, 3, 3],
                            box_1=[4, 4, 4],
                            bead_to_atom_name_dict={'_CH3': 'C'}
                            )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_pdb_residue_reorder_and_FF_filename_box_sizing(self, EthanolGOMC, EthaneGOMC):
        test_box_ethane_EthanolGOMC = mb.fill_box(compound=[EthanolGOMC, EthaneGOMC],
                                                    n_compounds=[1, 1],
                                                    box=[2.0, 2.0, 2.0])
        charmm = Charmm(test_box_ethane_EthanolGOMC, 'residue_reorder_box_sizing_box_0',
                        structure_box_1=EthaneGOMC,
                        filename_box_1='residue_reorder_box_sizing_box_1',
                        FF_filename=None,
                        residues=[EthaneGOMC.name, EthanolGOMC.name],
                        forcefield_selection= str(forcefields.get_ff_path()[0]) +  '/xml/' + 'oplsaa.xml' ,
                        fix_residue=None,
                        fix_residue_in_box=None,
                        fix_res_bonds_angles=None,
                        reorder_res_in_pdb_psf=True,
                        box_0=[3, 3, 3],
                        box_1=[4, 4, 4],
                        bead_to_atom_name_dict={'_CH3': 'C'}
                        )
        charmm.write_pdb()



        out_GOMC = open('residue_reorder_box_sizing_box_0.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] ==['CRYST1', '30.000', '30.000',  '30.000',
                                                       '90.00', '90.00', '90.00']
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', 'A', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', 'A', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', 'A', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', 'A', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', 'A', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', 'A', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', 'A', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', 'A', '1'],
                                             ['ATOM', '9', 'C1', 'ETO', 'A', '2'],
                                             ['ATOM', '10', 'C2', 'ETO', 'A', '2'],
                                             ['ATOM', '11', 'O1', 'ETO', 'A', '2'],
                                             ['ATOM', '12', 'H1', 'ETO', 'A', '2'],
                                             ['ATOM', '13', 'H2', 'ETO', 'A', '2'],
                                             ['ATOM', '14', 'H3', 'ETO', 'A', '2'],
                                             ['ATOM', '15', 'H4', 'ETO', 'A', '2'],
                                             ['ATOM', '16', 'H5', 'ETO', 'A', '2'],
                                             ['ATOM', '17', 'H6', 'ETO', 'A', '2']
                                             ]

                Atom_type_res_part_2_list = [['1.00', '0.00', 'C'], ['1.00', '0.00', 'C'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'C'], ['1.00', '0.00', 'C'], ['1.00', '0.00', 'O'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H']
                                             ]
                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass

        out_GOMC = open('residue_reorder_box_sizing_box_1.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] == ['CRYST1', '40.000', '40.000', '40.000',
                                                    '90.00', '90.00', '90.00']
            else:
                pass




    # test utils base 10 to base 16 converter
    def test_base_10_to_base_16(self):
        List_Base_10_and_16 = [[15, 'f'], [16, '10'], [17, '11'], [200, 'c8'], [1000, '3e8'], [5000, '1388'],
                              [int(16**3 - 1), 'fff'], [int(16**3), '1000']
                               ]

        for test_base_16_iter in  range(0, len(List_Base_10_and_16)):
            test_10_iter = List_Base_10_and_16[ test_base_16_iter][0]
            test_16_iter = List_Base_10_and_16[test_base_16_iter][1]
            assert str(base10_to_base16_alph_num(test_10_iter)) == str(test_16_iter)

        unique_entries_base_16_List = []
        for test_unique_base_16 in range(0, 16**2):
            unique_entries_base_16_List.append(base10_to_base16_alph_num(test_unique_base_16))

        verified_unique_entries_base_16_List = np.unique(unique_entries_base_16_List)
        assert len(verified_unique_entries_base_16_List) == len(unique_entries_base_16_List)

        add_same_values_List = ['1', 'a']
        for add_same_base_16 in range(0, len(add_same_values_List)):
            verified_unique_entries_base_16_List = np.append(verified_unique_entries_base_16_List,
                                                             add_same_values_List[add_same_base_16])
        assert len(verified_unique_entries_base_16_List) - len(add_same_values_List) == len(unique_entries_base_16_List)


    # test utils base 10 to base 26 converter
    def test_base_10_to_base_26(self):
        List_Base_10_and_26 = [[0, 'A'], [5, 'F'], [25, 'Z'], [26, 'BA'],
                               [200, 'HS'], [1000, 'BMM'], [5000, 'HKI'],
                               [int(26**3 - 1), 'ZZZ'], [int(26**3 ), 'BAAA']]

        for test_base_26_iter in range(0, len(List_Base_10_and_26)):
            test_10_iter = List_Base_10_and_26[test_base_26_iter][0]
            test_26_iter = List_Base_10_and_26[test_base_26_iter][1]
            assert str(base10_to_base26_alph(test_10_iter)) == str(test_26_iter)

        unique_entries_base_26_List = []
        for test_unique_base_26 in range(0, 26**2):
            unique_entries_base_26_List.append(base10_to_base26_alph(test_unique_base_26))

        verified_unique_entries_base_26_List = np.unique(unique_entries_base_26_List)
        assert len(verified_unique_entries_base_26_List) == len(unique_entries_base_26_List)

        add_same_values_List = ['1', 'a']
        for add_same_base_26 in range(0, len(add_same_values_List)):
            verified_unique_entries_base_26_List = np.append(verified_unique_entries_base_26_List,
                                                             add_same_values_List[add_same_base_26])
        assert len(verified_unique_entries_base_26_List) - len(add_same_values_List) == len(
            unique_entries_base_26_List)

    # test utils base 10 to base 52 converter
    def test_base_10_to_base_52(self):
        List_Base_10_and_52 = [[17, 'R'], [51, 'z'], [52, 'BA'], [53, 'BB'],
                               [200, 'Ds'], [1000, 'TM'], [5000, 'BsI'],
                              [int(52**3 - 1), 'zzz'], [int(52**3), 'BAAA']
                               ]

        for test_base_52_iter in range(0, len(List_Base_10_and_52)):
            test_10_iter = List_Base_10_and_52[test_base_52_iter][0]
            test_52_iter = List_Base_10_and_52[test_base_52_iter][1]
            assert str(base10_to_base52_alph(test_10_iter)) == str(test_52_iter)

        unique_entries_base_52_List = []
        for test_unique_base_52 in range(0, 52**2):
            unique_entries_base_52_List.append(base10_to_base52_alph(test_unique_base_52))

        verified_unique_entries_base_52_List = np.unique(unique_entries_base_52_List)
        assert len(verified_unique_entries_base_52_List) == len(unique_entries_base_52_List)

        add_same_values_List = ['1', 'a']
        for add_same_base_52 in range(0, len(add_same_values_List)):
            verified_unique_entries_base_52_List = np.append(verified_unique_entries_base_52_List,
                                                             add_same_values_List[add_same_base_52])
        assert len(verified_unique_entries_base_52_List) - len(add_same_values_List) == len(
            unique_entries_base_52_List)

    # test utils base 10 to base 62 converter
    def test_base_10_to_base_62(self):
        List_Base_10_and_62 = [[17, 'H'], [61, 'z'], [62, '10'], [63, '11'], [200, '3E'], [1000, 'G8'], [5000, '1Ie'],
                              [int(62**3 - 1), 'zzz'], [int(62**3), '1000']
                               ]

        for test_base_62_iter in  range(0, len(List_Base_10_and_62)):
            test_10_iter = List_Base_10_and_62[ test_base_62_iter][0]
            test_62_iter = List_Base_10_and_62[test_base_62_iter][1]
            assert str(base10_to_base62_alph_num(test_10_iter)) == str(test_62_iter)

        unique_entries_base_62_List = []
        for test_unique_base_62 in range(0,62**2):
            unique_entries_base_62_List.append(base10_to_base62_alph_num(test_unique_base_62))

        verified_unique_entries_base_62_List = np.unique(unique_entries_base_62_List )
        assert len(verified_unique_entries_base_62_List) == len(unique_entries_base_62_List)

        add_same_values_List = ['1', 'a']
        for add_same_base_62 in range(0, len(add_same_values_List)):
            verified_unique_entries_base_62_List = np.append(verified_unique_entries_base_62_List,
                                                             add_same_values_List[add_same_base_62])
        assert len(verified_unique_entries_base_62_List)-len(add_same_values_List) == len(unique_entries_base_62_List)


    # Tests for the mbuild.utils.specific_FF_to_residue.Specific_FF_to_residue() function

    def test_Specific_FF_to_box_value_negative(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection='oplsaa',
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[1,-2,3],
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_box_value_str(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection='oplsaa',
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[1,'2',3],
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None


    def test_Specific_FF_FF_is_None(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection=None,
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_wrong_FF_extention(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection='oplsaa.pdb',
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_all_residue_not_input(self, EthaneGOMC, EthanolGOMC):
        box = mb.fill_box(compound=[EthaneGOMC, EthanolGOMC],
                                      box=[1, 1, 1], n_compounds=[1,1])

        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(box ,
                                                            forcefield_selection='oplsaa',
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_residue_FFselection_not_dict(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection='oplsaa',
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )
        print('here 999999999999999999999999999')
        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_residue_is_None(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name : 'oplsaa'},
                                                            residues=None,
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_residue_reorder_not_True_or_False(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name : 'oplsaa'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=None,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_box_one_dim_is_negative(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[-2,3,4,5],
                                                            boxes_for_simulation=1
                                                            )
        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_box_one_dim_is_string(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=["string", 3, 4, 5],
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_simulation_boxes_not_1_or_2(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[2, 3, 4, 5],
                                                            boxes_for_simulation=3
                                                            )

        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None

    def test_Specific_FF_to_residue_FFselection_wrong_path(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name: 'oplsaa.xml'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[4, 5, 6],
                                                            boxes_for_simulation=1
                                                            )
        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None


    def test_Specific_FF_to_residue_FFselection_run(self, EthaneGOMC):
        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(EthaneGOMC,
                                                            forcefield_selection={EthaneGOMC.name: forcefields.get_ff_path()[0]
                                                                                               +'/xml/'+'oplsaa.xml'},
                                                            residues=[EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[4,5,6],
                                                            boxes_for_simulation=1
                                                            )
        assert str(Test_value_0) == "<Structure 8 atoms; 1 residues; 7 bonds; PBC (orthogonal); parametrized>"
        assert Test_value_1 == {'ETH': 0.5}
        assert Test_value_2 == {'ETH': 0.5}
        assert Test_value_3 ==  ['ETH']


    def test_Specific_FF_to_no_atoms_in_residue(self):
        Empty_compound = mb.Compound()

        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(Empty_compound,
                                                            forcefield_selection={'Empty_compound': 'oplsaa'},
                                                            residues=[],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=[5, 6, 7],
                                                            boxes_for_simulation=1
                                                            )
        print('Test_value_0 = ' +str(Test_value_0))
        assert Test_value_0 is None
        assert Test_value_1 is None
        assert Test_value_2 is None
        assert Test_value_3 is None




    def test_charmm_correct_residue_format(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename=None,
                                residues = [EthaneGOMC.name],
                                forcefield_selection = {EthaneGOMC.name : 'oplsaa'},
                                )

        except:
            Test_value = "TEST_FAILED"

        assert Test_value.input_error is False

    def test_charmm_residue_not_None(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1=None,
                                filename_box_1=None,
                                FF_filename=None,
                                residues=EthaneGOMC.name,
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'}
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"

    def test_charmm_residue_string(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1=None,
                                filename_box_1=None,
                                FF_filename=None,
                                residues='EthaneGOMC.name',
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'}
                                )

        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"

    def test_charmm_residue_is_None(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1=None,
                                filename_box_1=None,
                                FF_filename=None,
                                residues=None,
                                forcefield_selections={EthaneGOMC.name: 'oplsaa'}
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"

    def test_charmm_filename_0_is_string(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 0,
                                structure_box_1=None,
                                filename_box_1=None,
                                FF_filename=None,
                                residues= [EthaneGOMC.name],
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'}
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_filename_box_1_is_string(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = EthaneGOMC,
                                filename_box_1 = ['box_0'],
                                FF_filename=None,
                                residues= EthaneGOMC.name,
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_box_1_not_None_no_structure_box_1(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename=None,
                                residues=[EthaneGOMC.name],
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                box_1=[4, 4, 4],
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_GOMC_filename_not_string(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename=0,
                                residues=[EthaneGOMC.name],
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_GOMC_filename_ext_not_dot_inp(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename='box.test',
                                residues=[EthaneGOMC.name],
                                forcefield_selection={EthaneGOMC.name: 'oplsaa'},
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_FFselection_not_dict(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename='box.test',
                                residues=[EthaneGOMC.name],
                                forcefield_selection=['oplsaa', 'oplsaa'],
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_FFselection_string(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                structure_box_1 = None,
                                filename_box_1 = None,
                                FF_filename='box.test',
                                residues=[EthaneGOMC.name],
                                forcefield_selection='oplsaa',
                                )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_Residue_name_not_in_residues(self, EthaneGOMC):
        try:
            Test_value = Charmm(EthaneGOMC, 'box_0',
                                  structure_box_1 = None,
                                  filename_box_1 = None,
                                  FF_filename='box.test',
                                  residues=["XXX"],
                                  forcefield_selection='oplsaa',
                                  )
        except:
            Test_value = "TEST_FAILED"

        assert Test_value == "TEST_FAILED"


    def test_charmm_Methane_test_no_children(self, MethaneUAGOMC):

        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(MethaneUAGOMC,
                                                            forcefield_selection={MethaneUAGOMC.name: 'trappe-ua'},
                                                            residues=[MethaneUAGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )
        print('Test_value_3 = ' +str(Test_value_3))
        assert str(Test_value_0) == "<Structure 1 atoms; 1 residues; 0 bonds; PBC (orthogonal); NOT parametrized>"
        assert Test_value_1 == {'_CH4': 0.0}
        assert Test_value_2 == {'_CH4': 0.0}
        assert Test_value_3 == ['_CH4']


    def test_charmm_a_few_mbuild_layers(self, EthaneGOMC, EthanolGOMC):
        box_reservior_1 = mb.fill_box(compound=[EthaneGOMC],
                                    box=[1, 1, 1], n_compounds=[1])
        box_reservior_1.periodicity[0] = 2
        box_reservior_1.periodicity[1] = 2
        box_reservior_1.periodicity[1] = 2
        box_reservior_2 = mb.fill_box(compound=[EthanolGOMC],
                                    box=[1, 1, 1], n_compounds=[1])
        box_reservior_2.translate([0, 0, 1])
        box_reservior_1.add(box_reservior_2, inherit_periodicity=False)


        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(box_reservior_1,
                                                            forcefield_selection={EthanolGOMC.name: 'oplsaa',
                                                                                  EthaneGOMC.name:  'oplsaa',},
                                                            residues=[EthanolGOMC.name, EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert str(Test_value_0) == '<Structure 17 atoms; 2 residues; 15 bonds; PBC (orthogonal); parametrized>'
        assert Test_value_1 == {'ETO': 0.5, 'ETH': 0.5}
        assert Test_value_2 == {'ETO': 0.5, 'ETH': 0.5}
        assert Test_value_3 == ['ETH', 'ETO']

    def test_charmm_all_residues_not_in_dict(self, EthaneGOMC, EthanolGOMC):
        box_reservior_1 = mb.fill_box(compound=[EthaneGOMC],
                                      box=[1, 1, 1], n_compounds=[1])
        box_reservior_1.periodicity[0] = 2
        box_reservior_1.periodicity[1] = 2
        box_reservior_1.periodicity[1] = 2
        box_reservior_2 = mb.fill_box(compound=[EthanolGOMC],
                                      box=[1, 1, 1], n_compounds=[1])
        box_reservior_2.translate([0, 0, 1])
        box_reservior_1.add(box_reservior_2, inherit_periodicity=False)


        Test_value_0, Test_value_1, \
        Test_value_2, Test_value_3 = specific_ff_to_residue(box_reservior_1,
                                                            forcefield_selection={EthanolGOMC.name: 'oplsaa' },
                                                            residues=[EthanolGOMC.name, EthaneGOMC.name],
                                                            reorder_res_in_pdb_psf=False,
                                                            box=None,
                                                            boxes_for_simulation=1
                                                            )

        assert Test_value_0 == None
        assert Test_value_1 == None
        assert Test_value_2 == None
        assert Test_value_3 == None


    def test_FFselection_string(self, TwoPropanolUA):
        charmm = Charmm(TwoPropanolUA, 'FFselection_string', FF_filename='FFselection_string',
                        residues = [TwoPropanolUA.name],
                        forcefield_selection = forcefields.get_ff_path()[0]+'/xml/'+'trappe-ua.xml',
                        bead_to_atom_name_dict= {'_CH3' : 'C'})
        charmm.write_pdb()

        out_GOMC = open('FFselection_string.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'POL', 'A', '1'],
                                             ['ATOM', '2', 'BD1', 'POL', 'A', '1'],
                                             ['ATOM', '3', 'O1', 'POL', 'A', '1'],
                                             ['ATOM', '4', 'H1', 'POL', 'A', '1'],
                                             ['ATOM', '5', 'C2', 'POL', 'A', '1'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00', 'EP'], ['1.00', '0.00','EP'], ['1.00', '0.00','O'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','EP'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass



    def test_FFselection_list(self, TwoPropanolUA):
        try:
            Test_value_0 = Charmm(TwoPropanolUA, 'S', FF_filename='S',
                                  residues = [TwoPropanolUA.name],
                                  forcefield_selection = [str(forcefields.get_ff_path()[0])+'/xml/'+'trappe-ua.xml'],
                                  bead_to_atom_name_dict= {'_CH3' : 'C'},
                                  )
        except:
            Test_value_0 = "TEST_FAILED"

        assert Test_value_0 == "TEST_FAILED"



    def test_residuals_not_a_string(self, TwoPropanolUA):
        try:
            Test_value_0 = Charmm(TwoPropanolUA, 'box_0', FF_filename='box_0',
                                  residues = TwoPropanolUA.name,
                                  forcefield_selection = {TwoPropanolUA.name: 'trappe-ua' },
                                  bead_to_atom_name_dict= {'_CH3' : 'C'},
                                  )

        except:
            Test_value_0 = "TEST_FAILED"

        assert Test_value_0 == "TEST_FAILED"



    #Charmm writer sub-function testing
    def test_charmm_bond_reorder_angle_UreyBradleys(self, TwoPropanolGOMC, EthanolGOMC):
        box_reservior_0 = mb.fill_box(compound=[TwoPropanolGOMC, EthanolGOMC],
                                      box=[2, 2, 2], n_compounds=[2,2])

        try:
            structure_FF, \
            coulomb14scaler_dict, \
            LJ14scaler_dict, \
            residues_applied_list  =  specific_ff_to_residue(box_reservior_0,
                                                             forcefield_selection={TwoPropanolGOMC.name: 'oplsaa',
                                                                               EthanolGOMC.name: 'oplsaa' },
                                                             residues=[EthanolGOMC.name, TwoPropanolGOMC.name],
                                                             reorder_res_in_pdb_psf=False,
                                                             box=None,
                                                             boxes_for_simulation=1
                                                             )

            sigma_conversion_factor = 1
            epsilon_conversion_factor = 1
            #reversed the bond order so it fixes itself
            bonds_1 = [[bond.atom1.idx + 1, bond.atom2.idx + 1] for bond in structure_FF.bonds]
            bond_types_1, unique_bond_types_1 = charmm_writer._get_bond_types(structure_FF,
                                                                              sigma_conversion_factor,
                                                                              epsilon_conversion_factor)

            bonds_2 = [[bond.atom2.idx + 1, bond.atom1.idx + 1] for bond in structure_FF.bonds]
            bond_types_2, unique_bond_types_2 = charmm_writer._get_bond_types(structure_FF,
                                                                              sigma_conversion_factor,
                                                                              epsilon_conversion_factor)



            assert bonds_1 != bonds_2
            assert bond_types_1 == bond_types_2
            assert unique_bond_types_1 ==  unique_bond_types_2

            # test for error if trying to use urey_bradleys in th angles
            use_urey_bradleys = True
            angle_types_1, unique_angle_types_1 = charmm_writer._get_angle_types(
                structure_FF,
                sigma_conversion_factor,
                epsilon_conversion_factor,
                use_urey_bradleys=use_urey_bradleys
            )

            assert angle_types_1 is None
            assert unique_angle_types_1 is None

            Test_value_0 = "TEST_PASSED"

        except:
            Test_value_0 = "TEST_FAILED"

        assert Test_value_0 == "TEST_PASSED"

    # test for error if trying to use  use_dihedrals and impropers in the dihedrals (i.e. only RB torsion allowed)
    def test_charmm_dihedral_reorder(self, EthylEtherGOMC, MethlyEtherGOMC):
        try:
            box_reservior_0 = mb.fill_box(compound=[EthylEtherGOMC, MethlyEtherGOMC],
                                          box=[10, 10, 10], n_compounds=[10, 10])

            structure_FF, \
            coulomb14scaler_dict, \
            LJ14scaler_dict, \
            residues_applied_list = specific_ff_to_residue(box_reservior_0,
                                                           forcefield_selection={EthylEtherGOMC.name: 'oplsaa',
                                                                                 MethlyEtherGOMC.name: 'oplsaa'},
                                                           residues=[EthylEtherGOMC.name, MethlyEtherGOMC.name],
                                                           reorder_res_in_pdb_psf=False,
                                                           box=None,
                                                           boxes_for_simulation=1
                                                           )

            use_rb_torsions_1 = False
            use_dihedrals_1 = True
            epsilon_conversion_factor = 1
            lj_unit = 1 / epsilon_conversion_factor
            dihedral_types_1, unique_dihedral_types_1 = charmm_writer._get_dihedral_types(structure_FF,
                                                                                          use_rb_torsions_1,
                                                                                          use_dihedrals_1,
                                                                                          epsilon_conversion_factor)
            assert dihedral_types_1 is None
            assert unique_dihedral_types_1  is None

            use_rb_torsions_2 = True
            use_dihedrals_2 = False

            dihedral_types_2, unique_dihedral_types_2 = charmm_writer._get_dihedral_types(structure_FF,
                                                                                          use_rb_torsions_2,
                                                                                          use_dihedrals_2,
                                                                                          epsilon_conversion_factor)

            unique_dih_typ_unsorted_2 = dict(enumerate(set([(round(dihedral.type.c0 * lj_unit, 3),
                                                             round(dihedral.type.c1 * lj_unit, 3),
                                                             round(dihedral.type.c2 * lj_unit, 3),
                                                             round(dihedral.type.c3 * lj_unit, 3),
                                                             round(dihedral.type.c4 * lj_unit, 3),
                                                             round(dihedral.type.c5 * lj_unit, 3),
                                                             round(dihedral.type.scee, 1),
                                                             round(dihedral.type.scnb, 1),
                                                             dihedral.atom1.type, dihedral.atom2.type,
                                                             dihedral.atom3.type, dihedral.atom4.type,
                                                             dihedral.atom1.residue.name, dihedral.atom2.residue.name,
                                                             dihedral.atom3.residue.name, dihedral.atom4.residue.name
                                                             ) for dihedral in structure_FF.rb_torsions])))

            unique_dih_typ_unsorted_2 = OrderedDict([(y, x + 1) for x, y in unique_dih_typ_unsorted_2.items()])

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
            improper_types_1, unique_improper_types_1 = charmm_writer._get_impropers(structure_FF,
                                                                                     epsilon_conversion_factor)

            assert str(improper_types_1) == '[]'
            assert str(unique_improper_types_1) == '{}'

            Test_value_0 = "TEST_PASSED"

        except:
            Test_value_0 = "TEST_FAILED"

        assert Test_value_0 == "TEST_PASSED"

    def test_charmm_angle_reorder(self, EthylEtherGOMC, MethlyEtherGOMC):
        try:
            box_reservior_0 = mb.fill_box(compound=[EthylEtherGOMC, MethlyEtherGOMC],
                                          box=[10, 10, 10], n_compounds=[10, 10])

            structure_FF, \
            coulomb14scaler_dict, \
            LJ14scaler_dict, \
            residues_applied_list = specific_ff_to_residue(box_reservior_0,
                                                           forcefield_selection={EthylEtherGOMC.name: 'oplsaa',
                                                                                 MethlyEtherGOMC.name: 'oplsaa'},
                                                           residues=[EthylEtherGOMC.name, MethlyEtherGOMC.name],
                                                           reorder_res_in_pdb_psf=False,
                                                           box=None,
                                                           boxes_for_simulation=1
                                                           )
            sigma_conversion_factor = 1
            epsilon_conversion_factor = 1
            use_urey_bradleys = False
            angle_types_1, unique_angle_types_1 = charmm_writer._get_angle_types(
                structure_FF,
                sigma_conversion_factor,
                epsilon_conversion_factor,
                use_urey_bradleys
            )

            #note this sorts all the posssible combinations, so this should be the same as the double check (i.e, both 10)
            unique_angle_types_1_unsorted = dict(enumerate(set([(round(angle.type.k * (
                    sigma_conversion_factor ** 2 / epsilon_conversion_factor), 3),
                                                      round(angle.type.theteq, 3),
                                                      angle.atom2.type,
                                                      tuple(sorted((angle.atom1.type, angle.atom3.type))),
                                                      angle.atom1.residue.name, angle.atom2.residue.name,
                                                      angle.atom3.residue.name
                                                      ) for angle in structure_FF.angles])))
            unique_angle_types_1_unsorted = OrderedDict([(y, x + 1) for x, y in unique_angle_types_1_unsorted.items()])

            assert len(unique_angle_types_1_unsorted) == 10
            assert len(unique_angle_types_1) == 10

            Test_value_0 = "TEST_PASSED"

        except:
            Test_value_0 = "TEST_FAILED"

        assert Test_value_0 == "TEST_PASSED"


    def test_Bead_AtomName_equal_3(self, TwoPropanolUA):
        # testing def unique_atom_naming in charmm_writer, expecting when failing
        try:
            box_reservior_0 = mb.fill_box(compound=[TwoPropanolUA],
                                          box=[10, 10, 10], n_compounds=[10])

            value_0 = Charmm(box_reservior_0, 'test_Bead_AtomName_equal_3',
                             FF_filename='test_Bead_AtomName_equal_3',
                             residues = [TwoPropanolUA.name],
                             forcefield_selection = 'trappe-ua',
                             bead_to_atom_name_dict= {'_CH3' : 'Cx', '_HC' : 'Cxx'}
                             )
            value_0.write_inp()
            value_0.write_pdb()
            value_0.write_psf()


        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"


    def test_residue_string(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                              residues = 'TwoPropanolUA.name', forcefield_selection = 'trappe-ua',
                              bead_to_atom_name_dict= {'_CH3' : 'C'}
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_res_bonds_angles_string(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                              residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                              bead_to_atom_name_dict= {'_CH3' : 'C'}, fix_res_bonds_angles='TwoPropanolUA.name'
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_res_bonds_angles_residue_not_in_system(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                              residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                              bead_to_atom_name_dict= {'_CH3' : 'C'}, fix_res_bonds_angles=['WNG']
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_residue_string(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                             bead_to_atom_name_dict= {'_CH3' : 'C'}, fix_residue = 'TwoPropanolUA.name'
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_residue_string_residue_not_in_system(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues = [TwoPropanolUA.name], forcefield_selection = 'trappe-ua',
                             bead_to_atom_name_dict= {'_CH3' : 'C'}, fix_residue = ['WNG']
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_residue_in_box_string(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'},
                             fix_residue_in_box='TwoPropanolUA.name'
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_fix_residue_in_box_string_residue_not_in_system(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'}, fix_residue_in_box=['WNG']
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_bead_to_atom_name_dict_list(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict=['_CH3', 'C']
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_0_4_dims(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA_box_0',
                             FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'},
                             box_0=[4, 5, 6, 6], box_1=[3, 4, 5]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_1_4_dims(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA_box_0',
                             structure_box_1=TwoPropanolUA, filename_box_1='charmm_data_UA_box_1',
                             FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'},
                             box_0=[4, 5, 6], box_1=[3, 4, 5, 6]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_0_negative_dims(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'}, box_0=[-3, 4, 5]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_1_negative_dims(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA_box_0',
                             structure_box_1=TwoPropanolUA, filename_box_1='charmm_data_UA_box_1',
                             FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'},
                             box_0=[4, 5, 6], box_1=[-3, 4, 5]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_0_string_dims(self, TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA', FF_filename='charmm_data_UA',
                             residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                             bead_to_atom_name_dict={'_CH3': 'C'}, box_0=['string', 5, 6]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_box_1_string_dims(self,  TwoPropanolUA):
        try:
            value_0 = Charmm(TwoPropanolUA, 'charmm_data_UA_box_0', FF_filename='charmm_data_UA',
                                        structure_box_1= TwoPropanolUA, filename_box_1='charmm_data_UA_box_1',
                                        residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                                        bead_to_atom_name_dict={'_CH3': 'C'},
                                        box_0=[ 4, 5, 6], box_1=['string', 5, 6]
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_2_box_residues_not_all_listed(self, EthaneGOMC, EthanolGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthanolGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthanolGOMC.name, EthanolGOMC.name], forcefield_selection='oplsaa',
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_all_residues_are_listed(self, EthaneGOMC, EthanolGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                                        structure_box_1=EthanolGOMC, filename_box_1='charmm_data_box_1',
                                        FF_filename='charmm_data',
                                        residues=[EthanolGOMC.name], forcefield_selection='oplsaa',
                                        )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"


    def test_box_1_empty_test_1(self, TwoPropanolUA):
        Empty_compound = mb.Compound()

        charmm = Charmm(TwoPropanolUA, 'charmm_filled_box',
                         structure_box_1= Empty_compound, filename_box_1='charmm_empty_box',
                         FF_filename='charmm_empty_box.inp',
                         residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                         bead_to_atom_name_dict={'_CH3': 'C'},
                         box_0=[4, 5, 6], box_1=[3, 4, 5]
                         )
        charmm.write_pdb()
        charmm.write_psf()
        out_GOMC = open('charmm_empty_box.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] == ['CRYST1', '30.000', '40.000', '50.000',
                                                    '90.00', '90.00', '90.00']
                assert  out_GOMC[i + 1 ].split() == ['END']

            else:
                pass

        out_GOMC = open('charmm_filled_box.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'POL', 'A', '1'],
                                             ['ATOM', '2', 'BD1', 'POL', 'A', '1'],
                                             ['ATOM', '3', 'O1', 'POL', 'A', '1'],
                                             ['ATOM', '4', 'H1', 'POL', 'A', '1'],
                                             ['ATOM', '5', 'C2', 'POL', 'A', '1'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00', 'EP'], ['1.00', '0.00','EP'], ['1.00', '0.00','O'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','EP'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass

    def test_box_1_empty_test_2(self, TwoPropanolUA):
        Empty_compound = mb.Compound()

        charmm = Charmm(Empty_compound, 'charmm_empty_box',
                         structure_box_1= TwoPropanolUA, filename_box_1='charmm_filled_box',
                         FF_filename='charmm_empty_box',
                         residues=[TwoPropanolUA.name], forcefield_selection='trappe-ua',
                         bead_to_atom_name_dict={'_CH3': 'C'},
                         box_0=[4, 5, 6], box_1=[3, 4, 5]
                         )
        charmm.write_pdb()
        charmm.write_psf()
        out_GOMC = open('charmm_empty_box.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] == ['CRYST1', '40.000', '50.000', '60.000',
                                                    '90.00', '90.00', '90.00']
                assert  out_GOMC[i + 1 ].split() == ['END']

            else:
                pass

        out_GOMC = open('charmm_filled_box.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'POL', 'A', '1'],
                                             ['ATOM', '2', 'BD1', 'POL', 'A', '1'],
                                             ['ATOM', '3', 'O1', 'POL', 'A', '1'],
                                             ['ATOM', '4', 'H1', 'POL', 'A', '1'],
                                             ['ATOM', '5', 'C2', 'POL', 'A', '1'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00', 'EP'], ['1.00', '0.00','EP'], ['1.00', '0.00','O'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','EP'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass

    def test_structure_box_0_not_mb_Compound(self, EthaneGOMC):
        try:
            value_0 = Charmm('EthaneGOMC', 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"


    def test_structure_box_1_not_mb_Compound(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1='EthaneGOMC', filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_residue_list_not_entered(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1='EthaneGOMC', filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=None, forcefield_selection='oplsaa',
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_FF_dict_not_entered(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1='EthaneGOMC', filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name], forcefield_selection=None,
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"


    def test_residues_not_None_not_not_list(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1='EthaneGOMC', filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues='EthaneGOMC.name', forcefield_selection='oplsaa',
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_Mie_non_bonded_type(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                             non_bonded_type = 'Mie'
                             )
            charmm.write_inp()
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_other_non_bonded_type(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name], forcefield_selection='oplsaa',
                             non_bonded_type = 'OTH'
                             )
            charmm.write_inp()
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"


    def test_diff_1_4_coul_scalers(self, EthaneGOMC, TwoPropanolUA):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=TwoPropanolUA, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name, TwoPropanolUA.name],
                             forcefield_selection = {EthaneGOMC.name : 'oplsaa', TwoPropanolUA.name : 'trappe-ua'},
                             )
        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_write_inp_wo_FFfilename(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename=None,
                             residues=[EthaneGOMC.name],
                             forcefield_selection = 'oplsaa',
                             )
            value_0.write_inp()

        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_write_inp_with_2_boxes(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name],
                             forcefield_selection='oplsaa',
                             )
            value_0.write_inp()

        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

    def test_extra_residue_specified(self, EthaneGOMC):
        try:
            value_0 = Charmm(EthaneGOMC, 'charmm_data_box_0',
                             structure_box_1=EthaneGOMC, filename_box_1='charmm_data_box_1',
                             FF_filename='charmm_data',
                             residues=[EthaneGOMC.name, 'XXX'],
                             forcefield_selection='oplsaa',
                             )

        except:
            value_0 = "TEST_FAILED"

        assert value_0 == "TEST_FAILED"

