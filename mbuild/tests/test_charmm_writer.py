import pytest
import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.formats.charmm_writer import charmm_psf_psb_FF
from mbuild.utils.io import has_foyer



ethane = mb.load('CC', smiles=True)
ethane.name = "ETH"

ethanol = mb.load('CCO', smiles=True)
ethanol.name = "ETO"



class TWO_PROPANOL(mb.Compound):
    def __init__(self, name="POL"):
        super(TWO_PROPANOL, self).__init__()
        self.name = name

        CH3_1_1 = mb.Particle(pos=[0.2, 0.0, 0.0], name='_CH3')
        HC_1_1 = mb.Particle(pos=[0.4, 0.0, 0.0], name='_HC')
        O_1_1 = mb.Particle(pos=[0.8, 0.0, 0.0], name='O')
        H_1_1 = mb.Particle(pos=[1.0, 0.0, 0.0], name='H')
        CH3_1_2 = mb.Particle(pos=[0.6, 0.0, 0.0], name='_CH3')
        self.add([CH3_1_1, HC_1_1, O_1_1, H_1_1, CH3_1_2 ])

        port_R_CH3_1_1 = mb.Port(anchor=CH3_1_1, orientation=[0.1, 0, 0], separation=0.05)
        port_L_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[-0.1, 0, 0], separation=0.05)
        port_R_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[0.1, 0, 0], separation=0.05)
        port_D_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[0, -0.1, 0], separation=0.05)
        port_L_CH3_1_2 = mb.Port(anchor=CH3_1_2, orientation=[-0.1, 0, 0], separation=0.05)
        port_L_O_1_1 = mb.Port(anchor=O_1_1, orientation=[-0.1, 0, 0], separation=0.05)
        port_R_O_1_1 = mb.Port(anchor=O_1_1, orientation=[0.1, 0, 0], separation=0.05)
        port_L_H_1_1 = mb.Port(anchor=H_1_1, orientation=[-0.1, 0, 0], separation=0.05)

        self.add(port_R_CH3_1_1, label='port_R_CH3_1_1')
        self.add(port_L_HC_1_1, label = 'port_L_HC_1_1')
        self.add(port_R_HC_1_1, label='port_R_HC_1_1')
        self.add(port_L_CH3_1_2, label = 'port_L_CH3_1_2')
        self.add(port_D_HC_1_1, label='port_D_HC_1_1')
        self.add(port_L_O_1_1, label = 'port_L_O_1_1')
        self.add(port_R_O_1_1, label = 'port_R_O_1_1')
        self.add(port_L_H_1_1, label = 'port_L_H_1_1')


        mb.force_overlap(move_this=HC_1_1,
                         from_positions=self['port_L_HC_1_1'],
                         to_positions=self['port_R_CH3_1_1'])
        mb.force_overlap(move_this=CH3_1_2,
                         from_positions=self['port_L_CH3_1_2'],
                         to_positions=self['port_R_HC_1_1'])
        mb.force_overlap(move_this=O_1_1,
                         from_positions=self['port_L_O_1_1'],
                         to_positions=self['port_D_HC_1_1'])
        mb.force_overlap(move_this=H_1_1,
                         from_positions=self['port_L_H_1_1'],
                         to_positions=self['port_R_O_1_1'])

two_propanol_UA = TWO_PROPANOL()
two_propanol_UA.name = "POL"
two_propanol_UA.energy_minimize(forcefield = 'trappe-ua'  , steps=10**9)


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestCharmmWriterData(BaseTest):

    def test_save(self, ethane):
        charmm_psf_psb_FF(ethane, 'ethane', GOMC_FF_filename='ethane',
                          residues = [ethane.name], forcefield_names = 'oplsaa')

    def test_save_charmm_GOMC_FF(self):
        charmm_psf_psb_FF(ethane, 'charmm_data', GOMC_FF_filename='charmm_data',
                          residues = [ethane.name], forcefield_names = 'oplsaa')

        out_GOMC = open('charmm_data.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName' \
               ' (i.e., atoms_type_per_utilized_FF)' in line:
                assert len(out_GOMC[i + 1].split('!')[0].split()) == 3
                assert out_GOMC[i + 1].split('!')[0].split()[0:3] == ['*', '1', '12.010780']
                assert len(out_GOMC[i + 2].split('!')[0].split()) == 3
                assert out_GOMC[i + 2].split('!')[0].split()[0:3] == ['*', '2', '1.007947']
                assert out_GOMC[i + 1].split()[4:5] == ['opls_135_ETH']
                assert out_GOMC[i + 2].split()[4:5] == ['opls_140_ETH']

            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                bond_types = [['1', '2', '340.0', '1.09'], ['1', '1', '268.0', '1.529']]
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
                angle_types = [['1', '1', '2', '37.5', '110.70000'], ['2', '1', '2', '33.0', '107.80000']]
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
                dihed_types = [['2', '1', '1', '2', '0.300000', '0', '90.0'],
                               ['2', '1', '1', '2', '0.000000', '1', '180.0'],
                               ['2', '1', '1', '2', '0.000000', '2', '0.0'],
                               ['2', '1', '1', '2', '-0.150000', '3', '180.0'],
                               ['2', '1', '1', '2', '0.000000', '4', '0.0'],
                               ['2', '1', '1', '2', '0.000000', '5', '180.0']
                               ]
                for j in range(0, len(dihed_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] ==dihed_types[j]

            elif '!atype 	ignored	epsilon 	Rmin/2 		ignored	eps,1-4		Rmin/2,1-4' \
                 '		  atom_type_per_utilized_FF' in line:
                NB_types = [['1', '0.00', '-0.066000000', '1.96430858454', '0.00', '-0.033000000',	'0.98215429227'],
                               ['2', '0.00', '-0.030000000', '1.40307756039', '0.00', '-0.015000000',	'0.70153878019']]

                for j in range(0, len(NB_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] == NB_types[j]

            else:
                pass

    def test_save_charmm_psf(self):
        charmm_psf_psb_FF(ethane, 'charmm_data', GOMC_FF_filename='charmm_data',
                          residues = [ethane.name], forcefield_names = 'oplsaa')

        out_GOMC = open('charmm_data.psf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '8 !NATOM' in line:
                Atom_type_charge_etc_list = [['1', 'SYS', '1', 'ETH', 'C1', '1', '-0.180000', '12.0108'],
                                             ['2', 'SYS', '1', 'ETH', 'C2', '1', '-0.180000', '12.0108'],
                                             ['3', 'SYS', '1', 'ETH', 'H1', '2', '0.060000', '1.0079'],
                                             ['4', 'SYS', '1', 'ETH', 'H2', '2', '0.060000', '1.0079'],
                                             ['5', 'SYS', '1', 'ETH', 'H3', '2', '0.060000', '1.0079'],
                                             ['6', 'SYS', '1', 'ETH', 'H4', '2', '0.060000', '1.0079'],
                                             ['7', 'SYS', '1', 'ETH', 'H5', '2', '0.060000', '1.0079'],
                                             ['8', 'SYS', '1', 'ETH', 'H6', '2', '0.060000', '1.0079']
                                             ]
                for j in range(0, len(Atom_type_charge_etc_list)):
                    assert out_GOMC[i + 1 + j ].split()[0:8] == Atom_type_charge_etc_list[j]

            else:
                pass


    def test_save_charmm_pdb(self):
        charmm_psf_psb_FF(ethane, 'charmm_data', GOMC_FF_filename='charmm_data',
                          residues = [ethane.name], forcefield_names = 'oplsaa')


        out_GOMC = open('charmm_data.pdb', 'r').readlines()

        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', '1', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', '1', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', '1', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', '1', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', '1', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', '1', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', '1', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', '1', '1']
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00','C'], ['1.00', '0.00','C'], ['1.00', '0.00','H'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','H'], ['1.00', '0.00','H'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','H'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] == Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass


    def test_save_charmm_UA_GOMC_FF(self, two_propanol_UA):
        charmm_psf_psb_FF(two_propanol_UA, 'charmm_data_UA', GOMC_FF_filename='charmm_data_UA',
                          residues = [two_propanol_UA.name], forcefield_names = 'trappe-ua',
                          Bead_to_atom_name_dict= {'_CH3' : 'C'})

        out_GOMC = open('charmm_data_UA.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName ' \
               '(i.e., atoms_type_per_utilized_FF)' in line:
                atom_types_1 = [['*', '1', '15.035000'], ['*', '2', '13.019000'],
                              ['*', '4', '15.999430'], ['*', '3', '1.007947']]
                atom_types_2= [['CH3_sp3_POL'], ['CH_O_POL'], ['O_POL'],['H_POL']]

                for j in range(0, len(atom_types_1)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 3
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:3] == atom_types_1[j]
                    assert out_GOMC[i + 1 + j].split()[4:5] == atom_types_2[j]

            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                bond_types = [['3', '4', '600.402', '0.945'], ['2', '4', '600.402', '1.43'],
                              ['1', '2', '600.402', '1.54']]
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
                angle_types = [['1', '2', '1', '62.1', '112.00000'], ['1', '2', '4', '50.078', '109.47000'],
                               ['2', '4', '3', '55.046', '108.50000']]
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
                dihedral_types = [['1', '2', '4', '3', '0.648000', '0', '90.0'],
                                  ['1', '2', '4', '3', '-0.392500', '1', '180.0'],
                                  ['1', '2', '4', '3', '-0.062500', '2', '0.0'],
                                  ['1', '2', '4', '3', '0.345500', '3', '180.0'],
                                  ['1', '2', '4', '3', '0.000000', '4', '0.0'],
                                  ['1', '2', '4', '3', '0.000000', '5', '180.0']
                                  ]
                for j in range(0, len(dihedral_types)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:7] == dihedral_types[j]

            elif '!atype 	ignored	epsilon 	Rmin/2 		ignored	eps,1-4		Rmin/2,1-4		  ' \
                 'atom_type_per_utilized_FF' in line:
                angle_types = [['1', '0.00', '-0.194745937', '2.10461634058', '0.00', '-0.000000000', '0.00000000000'],
                               ['2', '0.00', '-0.019872012', '2.43013033459', '0.00', '-0.000000000', '0.00000000000'],
                               ['4', '0.00', '-0.184809990', '1.69491769295', '0.00', '-0.000000000', '0.00000000000'],
                               ['3', '0.00', '-0.000000000', '5.61231024155', '0.00', '-0.000000000', '0.00000000000']
                               ]

                for j in range(0, len(angle_types)):
                    assert len(out_GOMC[i + 1+ j].split('!')[0].split()) == 7
                    assert out_GOMC[i + 1+ j].split('!')[0].split()[0:7] == angle_types[j]

            else:
                pass


    def test_save_charmm_UA_psf(self):
        charmm_psf_psb_FF(two_propanol_UA, 'charmm_data_UA', GOMC_FF_filename='charmm_data_UA',
                          residues = [two_propanol_UA.name], forcefield_names = 'trappe-ua',
                          Bead_to_atom_name_dict= {'_CH3' : 'C'})


        out_GOMC = open('charmm_data_UA.psf', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '5 !NATOM' in line:
                Atom_type_charge_etc_list = [['1', 'SYS', '1', 'POL', 'C1', '1', '0.000000', '15.0350'],
                                             ['2', 'SYS', '1', 'POL', 'BD1', '2', '0.265000', '13.0190'],
                                             ['3', 'SYS', '1', 'POL', 'O1', '4', '-0.700000', '15.9994'],
                                             ['4', 'SYS', '1', 'POL', 'H1', '3', '0.435000', '1.0079'],
                                             ['5', 'SYS', '1', 'POL', 'C2', '1', '0.000000', '15.0350'],
                                             ]

                for j in range(0, len(Atom_type_charge_etc_list)):
                    assert out_GOMC[i + 1 + j].split()[0:8] ==  Atom_type_charge_etc_list[j]

            else:
                pass


    def test_save_charmm_UA_pdb(self):
        charmm_psf_psb_FF(two_propanol_UA, 'charmm_data_UA', GOMC_FF_filename='charmm_data_UA',
                          residues = [two_propanol_UA.name], forcefield_names = 'trappe-ua',
                          Bead_to_atom_name_dict= {'_CH3' : 'C'})

        out_GOMC = open('charmm_data_UA.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'POL', '1', '1'],
                                             ['ATOM', '2', 'BD1', 'POL', '1', '1'],
                                             ['ATOM', '3', 'O1', 'POL', '1', '1'],
                                             ['ATOM', '4', 'H1', 'POL', '1', '1'],
                                             ['ATOM', '5', 'C2', 'POL', '1', '1'],
                                             ]
                Atom_type_res_part_2_list = [['1.00', '0.00', 'EP'], ['1.00', '0.00','EP'], ['1.00', '0.00','O'],
                                             ['1.00', '0.00','H'], ['1.00', '0.00','EP'] ]

                for j in range(0, len(Atom_type_res_part_1_list)):
                    assert out_GOMC[i + 1 + j].split()[0:6] ==  Atom_type_res_part_1_list[j]
                    assert out_GOMC[i + 1 + j].split()[9:12] == Atom_type_res_part_2_list[j]

            else:
                pass


    def test_charmm_pdb_fix_angle_bond_fix_atoms(self):
        test_box_ethane_propane = mb.fill_box(compound=[ethane,ethanol],
                                  n_compounds= [1,1] ,
                                  box=[2.0, 2.0, 2.0])
        charmm_psf_psb_FF(test_box_ethane_propane, 'Test_fixes_angle_bond_atoms',
                          GOMC_FF_filename='Test_fixes_angle_bond_atoms',
                          residues = [ethane.name, ethanol.name], forcefield_names = 'oplsaa',
                          fix_residue=[ethane.name],
                          fix_residue_in_box = [ethanol.name],
                          fix_res_bonds_angles=[ethane.name]
                          )

        out_GOMC = open('Test_fixes_angle_bond_atoms.inp', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if '! atom_types 	mass 		  atomTypeForceFieldName_ResidueName ' \
               '(i.e., atoms_type_per_utilized_FF)' in line:
                mass_type_1 = [ ['*', '1', '12.010780'], ['*', '3', '1.007947'], ['*', '2', '12.010780'],
                                ['*', '7', '12.010780'], ['*', '5', '15.999430'], ['*', '4', '1.007947'],
                                ['*', '6', '1.007947']
                                ]
                mass_type_2 = [ ['opls_135_ETH'], ['opls_140_ETH'], ['opls_135_ETO'], ['opls_157_ETO'],
                                ['opls_154_ETO'], ['opls_140_ETO'], ['opls_155_ETO'] ]

                for j in range(0, len(mass_type_1)):
                    assert len(out_GOMC[i + 1 + j].split('!')[0].split()) == 3
                    assert out_GOMC[i + 1 + j].split('!')[0].split()[0:3] ==  mass_type_1[j]
                    assert out_GOMC[i + 1+ j].split()[4:5] == mass_type_2[j]


            elif '!atom_types 	 Kb	b0 		  atoms_types_per_utilized_FF' in line:
                fixed_bond_types = [['1', '1', '999999999999', '1.529'], ['1', '3', '999999999999', '1.09'] ]
                total_bonds_evaluated = []
                total_fixed_bonds = []
                for j in range(0,  7):
                    if out_GOMC[i + 1 + j].split('!')[0].split()[0:4] == (fixed_bond_types[0] or fixed_bond_types[1]):
                        total_bonds_evaluated.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                    if out_GOMC[i + 1 + j].split('!')[0].split()[2:3] == ['999999999999']:
                        total_fixed_bonds.append(out_GOMC[i + 1 + j].split('!')[0].split()[0:4])
                assert total_bonds_evaluated.sort() == total_bonds_evaluated.sort()
                assert len(total_fixed_bonds) == len(fixed_bond_types)

            elif '!atom_types 		Ktheta	Theta0			  atoms_types_per_utilized_FF' in line:
                fixed_angle_types = [['1',	'1',	'3',	'999999999999',	'110.70000'],
                                     ['3',	'1',	'3',	'999999999999',	'107.80000']
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
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', '1', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', '1', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', '1', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', '1', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', '1', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', '1', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', '1', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', '1', '1'],
                                             ['ATOM', '9', 'C1', 'ETO', '1', '2'],
                                             ['ATOM', '10', 'C2', 'ETO', '1', '2'],
                                             ['ATOM', '11', 'O1', 'ETO', '1', '2'],
                                             ['ATOM', '12', 'H1', 'ETO', '1', '2'],
                                             ['ATOM', '13', 'H2', 'ETO', '1', '2'],
                                             ['ATOM', '14', 'H3', 'ETO', '1', '2'],
                                             ['ATOM', '15', 'H4', 'ETO', '1', '2'],
                                             ['ATOM', '16', 'H5', 'ETO', '1', '2'],
                                             ['ATOM', '17', 'H6', 'ETO', '1', '2'],
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

    def test_charmm_pdb_residue_reorder_box_sizing(self):
        test_box_ethane_two_propanol_UA = mb.fill_box(compound=[two_propanol_UA, ethane],
                                  n_compounds= [1,1] ,
                                  box=[2.0, 2.0, 2.0])
        charmm_psf_psb_FF(test_box_ethane_two_propanol_UA, 'residue_reorder_box_sizing_box_0',
                          structure_1 = ethane,
                          filename_1 = 'residue_reorder_box_sizing_box_1',
                          GOMC_FF_filename=None,
                          residues = [ethane.name, two_propanol_UA.name],
                          forcefield_names = {ethane.name : 'oplsaa', two_propanol_UA.name : 'trappe-ua'},
                          fix_residue=None,
                          fix_residue_in_box = None,
                          fix_res_bonds_angles=None,
                          reorder_res_in_pdb_psf=True,
                          box_0 = [3, 3, 3],
                          box_1 =[4, 4, 4],
                          Bead_to_atom_name_dict={'_CH3': 'C'}
                          )

        out_GOMC = open('residue_reorder_box_sizing_box_0.pdb', 'r').readlines()
        for i, line in enumerate(out_GOMC):
            if 'CRYST1' in line:
                assert out_GOMC[i].split()[0:7] ==['CRYST1', '30.000', '30.000',  '30.000',
                                                       '90.00', '90.00', '90.00']
            if 'CRYST1' in line:
                Atom_type_res_part_1_list = [['ATOM', '1', 'C1', 'ETH', '1', '1'],
                                             ['ATOM', '2', 'C2', 'ETH', '1', '1'],
                                             ['ATOM', '3', 'H1', 'ETH', '1', '1'],
                                             ['ATOM', '4', 'H2', 'ETH', '1', '1'],
                                             ['ATOM', '5', 'H3', 'ETH', '1', '1'],
                                             ['ATOM', '6', 'H4', 'ETH', '1', '1'],
                                             ['ATOM', '7', 'H5', 'ETH', '1', '1'],
                                             ['ATOM', '8', 'H6', 'ETH', '1', '1'],
                                             ['ATOM', '9', 'C1', 'POL', '1', '2'],
                                             ['ATOM', '10', 'BD1', 'POL', '1', '2'],
                                             ['ATOM', '11', 'O1', 'POL', '1', '2'],
                                             ['ATOM', '12', 'H1', 'POL', '1', '2'],
                                             ['ATOM', '13', 'C2', 'POL', '1', '2']
                                             ]

                Atom_type_res_part_2_list = [['1.00', '0.00', 'C'], ['1.00', '0.00', 'C'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'H'],
                                             ['1.00', '0.00', 'EP'], ['1.00', '0.00', 'EP'], ['1.00', '0.00', 'O'],
                                             ['1.00', '0.00', 'H'], ['1.00', '0.00', 'EP']
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