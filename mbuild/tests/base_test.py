import numpy as np
import pytest

import mbuild as mb
from mbuild.utils.geometry import calc_dihedral
from mbuild.utils.io import get_fn


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def ethane(self):
        from mbuild.examples import Ethane
        return Ethane()

    @pytest.fixture
    def methane(self):
        from mbuild.examples import Methane
        return Methane()

    @pytest.fixture
    def h2o(self):
        from mbuild.lib.moieties import H2O
        return H2O()

    @pytest.fixture
    def ch2(self):
        from mbuild.lib.moieties import CH2
        return CH2()

    @pytest.fixture
    def ester(self):
        from mbuild.lib.moieties import Ester
        return Ester()

    @pytest.fixture
    def ch3(self):
        from mbuild.lib.moieties import CH3
        return CH3()

    @pytest.fixture
    def c3(self):
        from mbuild.lib.atoms import C3
        return C3()

    @pytest.fixture
    def n4(self):
        from mbuild.lib.atoms import N4
        return N4()

    @pytest.fixture
    def hydrogen(self):
        from mbuild.lib.atoms import H
        return H()

    @pytest.fixture
    def betacristobalite(self):
        from mbuild.lib.surfaces import Betacristobalite
        return Betacristobalite()

    @pytest.fixture
    def propyl(self):
        from mbuild.examples import Alkane
        return Alkane(3, cap_front=True, cap_end=False)

    @pytest.fixture
    def hexane(self, propyl):
        class Hexane(mb.Compound):
            def __init__(self):
                super(Hexane, self).__init__()

                self.add(propyl, 'propyl1')
                self.add(mb.clone(propyl), 'propyl2')

                mb.force_overlap(self['propyl1'],
                                 self['propyl1']['down'],
                                 self['propyl2']['down'])
        return Hexane()

    @pytest.fixture
    def octane(self):
        from mbuild.examples import Alkane
        return Alkane(8, cap_front=True, cap_end=True)

    @pytest.fixture
    def sixpoints(self):
        molecule = mb.Compound()
        molecule.add(mb.Particle(name='C', pos=[5, 5, 5]), label='middle')
        molecule.add(mb.Particle(name='C', pos=[6, 5, 5]), label='right')
        molecule.add(mb.Particle(name='C', pos=[4, 5, 5]), label='left')
        molecule.add(mb.Port(anchor=molecule[0]), label='up')
        molecule['up'].translate([0, 1, 0])
        molecule.add(mb.Port(anchor=molecule[0]), label='down')
        molecule['down'].translate([0, -1, 0])
        molecule.add(mb.Particle(name='C', pos=[5, 5, 6]), label='front')
        molecule.add(mb.Particle(name='C', pos=[5, 5, 4]), label='back')
        molecule.generate_bonds('C', 'C', 0.9, 1.1)
        return molecule

    @pytest.fixture
    def benzene(self):
        compound = mb.load(get_fn('benzene.mol2'))
        compound.name = 'Benzene'
        return compound

    @pytest.fixture
    def rigid_benzene(self):
        compound = mb.load(get_fn('benzene.mol2'))
        compound.name = 'Benzene'
        compound.label_rigid_bodies()
        return compound

    @pytest.fixture
    def benzene_from_parts(self):
        ch = mb.load(get_fn('ch.mol2'))
        ch.name = 'CH'
        mb.translate(ch, -ch[0].pos)       
        ch.add(mb.Port(anchor=ch[0], separation=0.07), 'a')
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0], separation=0.07), 'b')
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))
        ch_copy = mb.clone(ch)

        benzene = mb.Compound(name='Benzene')
        benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch_copy)
            mb.force_overlap(move_this=ch_new,
                             from_positions=ch_new['a'],
                             to_positions=current['b'])
            current = ch_new
            benzene.add(ch_new)

        carbons = [p for p in benzene.particles_by_name('C')]
        benzene.add_bond((carbons[0],carbons[-1]))

        return benzene

    @pytest.fixture
    def box_of_benzenes(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4]) 
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        return filled

    @pytest.fixture
    def rigid_ch(self):
        ch = mb.load(get_fn('ch.mol2'))
        ch.name = 'CH'
        ch.label_rigid_bodies()
        mb.translate(ch, -ch[0].pos)    
        ch.add(mb.Port(anchor=ch[0]), 'a')
        mb.translate(ch['a'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['a'], 120.0 * (np.pi/180.0))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        mb.translate(ch['b'], [0, 0.07, 0]) 
        mb.rotate_around_z(ch['b'], -120.0 * (np.pi/180.0))
        return ch

    @pytest.fixture
    def silane(self):
        from mbuild.lib.moieties import Silane
        return Silane()

    @pytest.fixture
    def chf(self):
        class CHF(mb.Compound):
            def __init__(self):
                super(CHF, self).__init__()
                carbon = mb.Particle(name='C', pos=[0.0, 0.0, 0.0])
                hydrogen = mb.Particle(name='H', pos=[0.0, -0.15, 0.0])
                fluorine = mb.Particle(name='F', pos=[0.0, 0.15, 0.0])
                self.add([carbon, hydrogen, fluorine])
                self.add_bond((carbon, hydrogen))
                self.add_bond((carbon, fluorine))
        return CHF()

    @pytest.fixture
    def connect_and_reconnect(self, chf):
        def _connect_and_reconnect(chf, bond_vector):
            first = mb.clone(chf)
            second = mb.clone(chf)
            first.add(mb.Port(anchor=first[0], orientation=bond_vector,
                separation=0.075), label='up')
            second.add(mb.Port(anchor=second[0], orientation=-bond_vector,
                separation=0.075), label='down')
            c2h2f2 = mb.Compound(subcompounds=(first, second))
            mb.force_overlap(first, first['up'], second['down'])
            fccf_dihedral_init = calc_dihedral(first[2].pos, first[0].pos,
                second[0].pos, second[2].pos)
            c2h2f2.remove_bond((first[0], second[0]))
            mb.force_overlap(first, first['port[0]'], second['port[0]'])
            fccf_dihedral_final = calc_dihedral(first[2].pos, first[0].pos,
                second[0].pos, second[2].pos)
            return fccf_dihedral_init, fccf_dihedral_final
        return _connect_and_reconnect

    @pytest.fixture
    def alkane_monolayer(self):
        from mbuild.examples import AlkaneMonolayer
        pattern = mb.Random2DPattern(50, seed=1)
        monolayer = AlkaneMonolayer(pattern=pattern)
        return monolayer
