import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.utils.io import has_foyer


@pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
class TestLammpsData(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.lammps')

    @pytest.mark.parametrize('unit_style',['real', 'lj'])
    def test_save_forcefield(self, ethane, unit_style):
        ethane.save(filename='ethane-opls.lammps',
                forcefield_name='oplsaa', unit_style=unit_style)

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    def test_save_charmm(self):
        cmpd = mb.load(get_fn('charmm_dihedral.mol2'))
        for i in cmpd.particles():
            i.name = "_{}".format(i.name)
        structure = cmpd.to_parmed(box=cmpd.boundingbox, 
                                    residues=set([p.parent.name for \
                                                 p in cmpd.particles()]))

        from foyer import Forcefield
        ff = Forcefield(forcefield_files=[get_fn('charmm_truncated.xml')])
        structure = ff.apply(structure, assert_dihedral_params=False)

        from mbuild.formats.lammpsdata import write_lammpsdata
        write_lammpsdata(structure, 'charmm_dihedral.lammps')
        out_lammps = open('charmm_dihedral.lammps', 'r').readlines()
        for i, line in enumerate(out_lammps):
            if 'Angle Coeffs' in line:
                assert '# charmm' in line
                assert '#\tk(kcal/mol/rad^2)\t\ttheteq(deg)\tk(kcal/mol/angstrom^2)\treq(angstrom)\n' in out_lammps[i+1]
                assert len(out_lammps[i+2].split('#')[0].split()) == 5
            elif 'Dihedral Coeffs' in line:
                assert '# charmm' in line
                assert '#k, n, phi, weight' in out_lammps[i+1]
                assert len(out_lammps[i+2].split('#')[0].split()) == 5
            else:
                pass

    @pytest.mark.skipif(not has_foyer, reason="Foyer package not installed")
    @pytest.mark.parametrize('unit_style', ['real', 'lj'])
    def test_save_box(self, ethane, unit_style):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.lammps',
                forcefield_name='oplsaa', box=box,
                unit_style=unit_style)

    def test_nbfix(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        # Add nbfixes
        types = list(set([a.atom_type for a in structure.atoms]))
        types[0].add_nbfix(types[1].name, 1.2, 2.1)
        types[1].add_nbfix(types[0].name, 1.2, 2.1)
        write_lammpsdata(filename='nbfix.lammps', structure=structure)

        checked_section = False
        with open('nbfix.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'PairIJ Coeffs' in line:
                    fi.readline()
                    line = fi.readline().partition('#')[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float),
                        [1, 1, 0.066, 3.5])
                    line = fi.readline().partition('#')[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float),
                        [1, 2, 2.1, 1.06907846])
                    line = fi.readline().partition('#')[0]
                    assert np.allclose(
                        np.asarray(line.split(), dtype=float),
                        [2, 2, 0.03, 2.5])
                    line = fi.readline()
                    checked_section = True
                # Break if PairIJ Coeffs is not found
                if 'Atoms' in line:
                    break

    def test_save_triclinic_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]), angles=[60, 70, 80])
        ethane.save(filename='triclinic-box.lammps', forcefield_name='oplsaa', box=box)

    @pytest.mark.parametrize(
        'atom_style, n_columns',
        [('full', 7), ('atomic', 5), ('molecular', 6), ('charge', 6)]
    ) 
    def test_writing_atom_styles(self, ethane, atom_style, n_columns):
        ethane.save(filename='ethane.lammps', atom_style=atom_style)
        with open('ethane.lammps', 'r') as f:
            for line in f:
                if "Atoms" not in line:
                    continue
                    atoms_header = next(f)
                    first_atom_line = next(f)
                    columns = first_atom_line.split("\t")
                    assert len(columns) == n_columns

    def test_resid(self, ethane, methane):
        structure = ethane.to_parmed() + methane.to_parmed()
        n_atoms = len(structure.atoms)
        write_lammpsdata(structure, 'compound.lammps')
        res_list = list()
        with open('compound.lammps', 'r') as f:
            for i,line in enumerate(f):
                if 'Atoms' in line:
                    break
        atom_lines = open('compound.lammps', 'r').readlines()[i+2:i+n_atoms+2]
        for line in atom_lines:
            res_list.append(line.rstrip().split()[1])

        assert set(res_list) == set(['1', '0'])

    def test_lj_box(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section: 
                line = fi.readline()
                if 'dihedral types' in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 2.04)
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 2.268)
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 1.898857)
                    checked_section = True

    def test_lj_masses(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'Masses' in line:
                    fi.readline()
                    line = float(fi.readline().split()[1])
                    assert np.isclose(line, 1.00)
                    checked_section = True

    def test_lj_pairs(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'Pair Coeffs' in line:
                    fi.readline()
                    line = fi.readline().split()
                    epsilon = float(line[1])
                    sigma = float(line[2])
                    assert np.isclose(epsilon, 1.00)
                    assert np.isclose(sigma, 1.00)
                    checked_section = True

    def test_lj_bonds(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'Bond Coeffs' in line:
                    fi.readline()
                    bonds = list()
                    bonds.append(float(fi.readline().split()[1]))
                    bonds.append(float(fi.readline().split()[1]))
                    assert np.allclose(sorted(bonds), [49742.424, 63106.06])
                    checked_section = True

    def test_lj_angles(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'Angle Coeffs' in line:
                    fi.readline()
                    angles = list()
                    angles.append(float(fi.readline().split()[1]))
                    angles.append(float(fi.readline().split()[1]))

                    assert np.allclose(sorted(angles), [6125.0, 6960.227])
                    checked_section = True

    def test_lj_dihedrals(self, ethane):
        from foyer import Forcefield

        OPLSAA = Forcefield(name='oplsaa')
        structure = OPLSAA.apply(ethane)
        write_lammpsdata(filename='lj.lammps', structure=structure,
                unit_style='lj')

        checked_section = False
        with open('lj.lammps', 'r') as fi:
            while not checked_section:
                line = fi.readline()
                if 'Dihedral Coeffs' in line:
                    fi.readline()
                    dihedrals = fi.readline().split()[1:5]
                    dihedrals = [float(i) for i in dihedrals]
                    assert np.allclose(dihedrals, [0.0005, 0.0, 4.5455, -0.0])
                    checked_section=True
