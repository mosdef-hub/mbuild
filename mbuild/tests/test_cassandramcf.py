import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import has_foyer


class TestCassandraMCF(BaseTest):

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.mcf', forcefield_name='oplsaa', angle_style='harmonic', dihedral_style='opls')

        mcf_data = []
        with open('ethane-opls.mcf') as f:
            for line in f:
                mcf_data.append(line.strip().split()) 

        for idx,line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == 'Atom_Info':
                    atom_section_start = idx
                elif line[1] == 'Bond_Info':
                    bond_section_start = idx
                elif line[1] == 'Angle_Info':
                    angle_section_start = idx
                elif line[1] == 'Dihedral_Info':
                    dihedral_section_start = idx
                elif line[1] == 'Improper_Info':
                    improper_section_start = idx
                elif line[1] == 'Fragment_Info':
                    fragment_section_start = idx
                elif line[1] == 'Fragment_Connectivity':
                    fragment_conn_start = idx
        
        # Check a some atom info
        assert mcf_data[atom_section_start+1][0] == '8'
        assert mcf_data[atom_section_start+2][1] == 'opls_135'
        assert mcf_data[atom_section_start+2][3] == '12.011'
        assert mcf_data[atom_section_start+2][4] == '-0.180'
        assert mcf_data[atom_section_start+2][5] == 'LJ'
        assert mcf_data[atom_section_start+2][6] == '33.212'
        assert mcf_data[atom_section_start+2][7] == '3.500'

        # Bond info
        assert mcf_data[bond_section_start+1][0] == '7'
        assert mcf_data[bond_section_start+2][3] == 'fixed'
        assert mcf_data[bond_section_start+2][4] == '1.090'

        # Angle info
        assert mcf_data[angle_section_start+1][0] == '12'
        assert mcf_data[angle_section_start+2][4] == 'harmonic'
        assert mcf_data[angle_section_start+2][5] == '18870.7'
        assert mcf_data[angle_section_start+2][6] == '110.70'

        assert mcf_data[dihedral_section_start+1][0] == '9'
        assert mcf_data[dihedral_section_start+2][5] == 'OPLS'
        assert mcf_data[dihedral_section_start+2][6] == '0.000'
        assert mcf_data[dihedral_section_start+2][7] == '0.000'
        assert mcf_data[dihedral_section_start+2][8] == '-0.000'
        assert mcf_data[dihedral_section_start+2][9] == '0.628'

        assert mcf_data[improper_section_start+1][0] == '0'
        assert mcf_data[fragment_section_start+1][0] == '2'

        # Check fragment connectivity
        assert mcf_data[fragment_conn_start+1][0] == '1'
        assert mcf_data[fragment_conn_start+2][0] == '1'
        assert mcf_data[fragment_conn_start+2][1] == '1'
        assert mcf_data[fragment_conn_start+2][2] == '2'


    def test_save_ring_forcefield(self,benzene):
        benzene.save(filename='benzene-opls.mcf', forcefield_name='oplsaa', angle_style='fixed', dihedral_style='opls')

        mcf_data = []
        with open('benzene-opls.mcf') as f:
            for line in f:
                mcf_data.append(line.strip().split()) 

        for idx,line in enumerate(mcf_data):
            if len(line) > 1:
                if line[1] == 'Atom_Info':
                    atom_section_start = idx
                elif line[1] == 'Bond_Info':
                    bond_section_start = idx
                elif line[1] == 'Angle_Info':
                    angle_section_start = idx
                elif line[1] == 'Dihedral_Info':
                    dihedral_section_start = idx
                elif line[1] == 'Improper_Info':
                    improper_section_start = idx
                elif line[1] == 'Fragment_Info':
                    fragment_section_start = idx
                elif line[1] == 'Fragment_Connectivity':
                    fragment_conn_start = idx
        
        # Check a some atom info
        assert mcf_data[atom_section_start+1][0] == '12'
        assert mcf_data[atom_section_start+2][1] == 'opls_145'
        assert mcf_data[atom_section_start+2][3] == '12.011'
        assert mcf_data[atom_section_start+2][4] == '-0.115'
        assert mcf_data[atom_section_start+2][5] == 'LJ'
        assert mcf_data[atom_section_start+2][6] == '35.225'
        assert mcf_data[atom_section_start+2][7] == '3.550'
        assert mcf_data[atom_section_start+2][8] == 'ring'

        # Bond info
        assert mcf_data[bond_section_start+1][0] == '12'
        assert mcf_data[bond_section_start+2][3] == 'fixed'
        assert mcf_data[bond_section_start+2][4] == '1.400'

        # Angle info
        assert mcf_data[angle_section_start+1][0] == '18'
        assert mcf_data[angle_section_start+2][4] == 'fixed'
        assert mcf_data[angle_section_start+2][5] == '120.00'

        assert mcf_data[dihedral_section_start+1][0] == '24'
        assert mcf_data[dihedral_section_start+2][5] == 'OPLS'
        assert mcf_data[dihedral_section_start+2][6] == '0.000'
        assert mcf_data[dihedral_section_start+2][7] == '-0.000'
        assert mcf_data[dihedral_section_start+2][8] == '15.167'
        assert mcf_data[dihedral_section_start+2][9] == '-0.000'

        assert mcf_data[improper_section_start+1][0] == '0'
        assert mcf_data[fragment_section_start+1][0] == '1'

        # Check fragment connectivity
        assert mcf_data[fragment_conn_start+1][0] == '0'


