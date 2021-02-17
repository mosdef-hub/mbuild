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
        from mbuild.lib.molecules import Ethane
        return Ethane()

    @pytest.fixture
    def methane(self):
        from mbuild.lib.molecules import Methane
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
        from mbuild.lib.recipes import Alkane
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
        from mbuild.lib.recipes import Alkane
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
        ch.translate(-ch[0].pos)       
        ch.add(mb.Port(anchor=ch[0], separation=0.07), 'a')
        ch['a'].rotate(120.0 * (np.pi/180.0), around=np.asarray([0, 0, 1]))
        ch.add(mb.Port(anchor=ch[0], separation=0.07), 'b')
        ch['b'].rotate(-120.0 * (np.pi/180.0), around=np.asarray([0, 0, 1]))
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
        ch.translate(-ch[0].pos)    
        ch.add(mb.Port(anchor=ch[0]), 'a')
        ch['a'].translate([0, 0.07, 0]) 
        ch['a'].rotate(120.0 * (np.pi/180.0), around=np.asarray([0, 0, 1]))

        ch.add(mb.Port(anchor=ch[0]), 'b')
        ch['b'].translate([0, 0.07, 0]) 
        ch['b'].rotate(-120.0 * (np.pi/180.0), around=np.asarray([0, 0, 1]))

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
    def copper_cell(self):
        copper = mb.Compound(name='Cu')
        lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        spacing = [.36149, .36149, .36149]
        copper_locations = [[0., 0., 0.], [.5, .5, 0.],
                [.5, 0., .5], [0., .5, .5]]
        basis =  {'Cu' : copper_locations}
        copper_lattice = mb.Lattice(lattice_spacing = spacing,
                lattice_vectors=lattice_vector,
                lattice_points=basis)
        copper_dict = {'Cu': copper}
        copper_pillar = copper_lattice.populate(x=3, y=3, z=1,
                compound_dict=copper_dict)
        return copper_pillar

    @pytest.fixture
    def graphene(self):
        carbon = mb.Compound(name='C')
        angles = [90, 90, 120]
        carbon_locations = [[0, 0, 0], [2/3, 1/3, 0]]
        basis = {'C' : carbon_locations}
        graphene = mb.Lattice(lattice_spacing=[.2456, .2456, 0],
                               angles=angles, lattice_points=basis)
        carbon_dict = {'C' : carbon}
        graphene_cell = graphene.populate(compound_dict=carbon_dict,
                                          x=3, y=3, z=1)
        return graphene_cell

    @pytest.fixture
    def cscl_crystal(self):
        cesium = mb.Compound(name='Cs')
        chlorine = mb.Compound(name='Cl')
        spacing = [.4123, .4123, .4123]
        basis = {'Cs' : [[0.5, 0.5, 0.5]], 'Cl' : [[0, 0, 0]]}
        cscl_lattice = mb.Lattice(spacing, lattice_points=basis)
        
        
        cscl_dict = {'Cs' : cesium, 'Cl' : chlorine}
        cscl_compound = cscl_lattice.populate(x=3, y=3, z=1,
                                              compound_dict=cscl_dict)
        return cscl_compound



    @pytest.fixture
    def ethane_gomc(self):
        ethane_gomc = mb.load('CC', smiles=True)
        ethane_gomc.name = "ETH"

        return ethane_gomc

    @pytest.fixture
    def ethanol_gomc(self):
        ethanol_gomc = mb.load('CCO', smiles=True)
        ethanol_gomc.name = "ETO"

        return ethanol_gomc

    @pytest.fixture
    def methane_ua_gomc(self):
        methane_ua_gomc = mb.Compound(name="_CH4")

        return methane_ua_gomc

    @pytest.fixture
    def two_propanol_gomc(self):
        two_propanol_gomc = mb.load('CC(C)O', smiles=True)
        two_propanol_gomc.name = "TPR"
        return two_propanol_gomc

    @pytest.fixture
    def ethyl_ether_gomc(self):
        ethyl_ether_gomc = mb.load('CCOCC', smiles=True)
        ethyl_ether_gomc.name = "ETE"
        return ethyl_ether_gomc

    @pytest.fixture
    def methyl_ether_gomc(self):
        methyl_ether_gomc = mb.load('COC', smiles=True)
        methyl_ether_gomc.name = "MTE"
        return methyl_ether_gomc

    @pytest.fixture
    def two_propanol_ua(self):
        class TwoPropanolUA(mb.Compound):
            def __init__(self):
                super(TwoPropanolUA, self).__init__()
                self.name = "POL"

                CH3_1_1 = mb.Particle(pos=[0.2, 0.0, 0.0], name='_CH3')
                HC_1_1 = mb.Particle(pos=[0.4, 0.0, 0.0], name='_HC')
                O_1_1 = mb.Particle(pos=[0.8, 0.0, 0.0], name='O')
                H_1_1 = mb.Particle(pos=[1.0, 0.0, 0.0], name='H')
                CH3_1_2 = mb.Particle(pos=[0.6, 0.0, 0.0], name='_CH3')
                self.add([CH3_1_1, HC_1_1, O_1_1, H_1_1, CH3_1_2])

                port_R_CH3_1_1 = mb.Port(anchor=CH3_1_1, orientation=[0.1, 0, 0], separation=0.05)
                port_L_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[-0.1, 0, 0], separation=0.05)
                port_R_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[0.1, 0, 0], separation=0.05)
                port_D_HC_1_1 = mb.Port(anchor=HC_1_1, orientation=[0, -0.1, 0], separation=0.05)
                port_L_CH3_1_2 = mb.Port(anchor=CH3_1_2, orientation=[-0.1, 0, 0], separation=0.05)
                port_L_O_1_1 = mb.Port(anchor=O_1_1, orientation=[-0.1, 0, 0], separation=0.05)
                port_R_O_1_1 = mb.Port(anchor=O_1_1, orientation=[0.1, 0, 0], separation=0.05)
                port_L_H_1_1 = mb.Port(anchor=H_1_1, orientation=[-0.1, 0, 0], separation=0.05)

                self.add(port_R_CH3_1_1, label='port_R_CH3_1_1')
                self.add(port_L_HC_1_1, label='port_L_HC_1_1')
                self.add(port_R_HC_1_1, label='port_R_HC_1_1')
                self.add(port_L_CH3_1_2, label='port_L_CH3_1_2')
                self.add(port_D_HC_1_1, label='port_D_HC_1_1')
                self.add(port_L_O_1_1, label='port_L_O_1_1')
                self.add(port_R_O_1_1, label='port_R_O_1_1')
                self.add(port_L_H_1_1, label='port_L_H_1_1')

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

                self.energy_minimize(forcefield='trappe-ua', steps=10 ** 9)

        return TwoPropanolUA()
