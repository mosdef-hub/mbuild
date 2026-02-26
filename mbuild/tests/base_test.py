import numpy as np
import pytest

import mbuild as mb
from mbuild import Polymer
from mbuild.utils.geometry import calc_dihedral
from mbuild.utils.io import get_fn


def radius_of_gyration(coordinates):
    """Calculate the square radius of gyration for a set of coordinates using the geometric center."""
    coordinates = np.array(coordinates)
    geometric_center = np.mean(coordinates, axis=0)
    squared_distances = np.sum((coordinates - geometric_center) ** 2, axis=1)
    rg2 = np.mean(squared_distances)
    return rg2


class BaseTest:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def ethane(self):
        from mbuild.lib.molecules import Ethane

        return Ethane()

    @pytest.fixture
    def ethane_monomer(self):
        ethane = mb.load("C{<}C{>}", smiles=True)
        return ethane

    @pytest.fixture
    def ethane_chain(self):
        chain = Polymer()
        ethane = mb.load("C{<}C{>}", smiles=True)
        chain.add_monomer(ethane, head_tag=">", tail_tag="<", separation=0.145)
        return chain

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

                self.add(propyl, "propyl1")
                self.add(mb.clone(propyl), "propyl2")

                mb.force_overlap(
                    self["propyl1"],
                    self["propyl1"]["down"],
                    self["propyl2"]["down"],
                )

        return Hexane()

    @pytest.fixture
    def octane(self):
        from mbuild.lib.recipes import Alkane

        return Alkane(8, cap_front=True, cap_end=True)

    @pytest.fixture
    def sixpoints(self):
        molecule = mb.Compound()
        molecule.add(mb.Particle(name="C", pos=[5, 5, 5]), label="middle")
        molecule.add(mb.Particle(name="C", pos=[6, 5, 5]), label="right")
        molecule.add(mb.Particle(name="C", pos=[4, 5, 5]), label="left")
        molecule.add(mb.Port(anchor=molecule[0]), label="up")
        molecule["up"].translate([0, 1, 0])
        molecule.add(mb.Port(anchor=molecule[0]), label="down")
        molecule["down"].translate([0, -1, 0])
        molecule.add(mb.Particle(name="C", pos=[5, 5, 6]), label="front")
        molecule.add(mb.Particle(name="C", pos=[5, 5, 4]), label="back")
        molecule.generate_bonds("C", "C", 0.9, 1.1)
        return molecule

    @pytest.fixture
    def benzene(self):
        compound = mb.load(get_fn("benzene.mol2"))
        compound.name = "Benzene"
        return compound

    @pytest.fixture
    def benzene_from_SMILES(self):
        compound = mb.load("c1ccccc1", smiles=True)
        compound.name = "Benzene"
        return compound

    @pytest.fixture
    def benzene_from_parts(self):
        ch = mb.load(get_fn("ch.mol2"))
        ch.name = "CH"
        ch.translate(-ch[0].pos)
        ch.add(mb.Port(anchor=ch[0], separation=0.07), "a")
        ch["a"].rotate(120.0 * (np.pi / 180.0), around=np.asarray([0, 0, 1]))
        ch.add(mb.Port(anchor=ch[0], separation=0.07), "b")
        ch["b"].rotate(-120.0 * (np.pi / 180.0), around=np.asarray([0, 0, 1]))
        ch_copy = mb.clone(ch)

        benzene = mb.Compound(name="Benzene")
        benzene.add(ch)
        current = ch

        for _ in range(5):
            ch_new = mb.clone(ch_copy)
            mb.force_overlap(
                move_this=ch_new,
                from_positions=ch_new["a"],
                to_positions=current["b"],
            )
            current = ch_new
            benzene.add(ch_new)

        carbons = [p for p in benzene.particles_by_name("C")]
        benzene.add_bond((carbons[0], carbons[-1]))

        return benzene

    @pytest.fixture
    def box_of_benzenes(self, benzene):
        n_benzenes = 10
        benzene.name = "Benzene"
        filled = mb.fill_box(benzene, n_compounds=n_benzenes, box=[0, 0, 0, 4, 4, 4])
        return filled

    @pytest.fixture
    def silane(self):
        from mbuild.lib.moieties import Silane

        return Silane()

    @pytest.fixture
    def chf(self):
        class CHF(mb.Compound):
            def __init__(self):
                super(CHF, self).__init__()
                carbon = mb.Particle(name="C", pos=[0.0, 0.0, 0.0])
                hydrogen = mb.Particle(name="H", pos=[0.0, -0.15, 0.0])
                fluorine = mb.Particle(name="F", pos=[0.0, 0.15, 0.0])
                self.add([carbon, hydrogen, fluorine])
                self.add_bond((carbon, hydrogen))
                self.add_bond((carbon, fluorine))

        return CHF()

    @pytest.fixture
    def connect_and_reconnect(self, chf):
        def _connect_and_reconnect(chf, bond_vector):
            first = mb.clone(chf)
            second = mb.clone(chf)
            first.add(
                mb.Port(anchor=first[0], orientation=bond_vector, separation=0.075),
                label="up",
            )
            second.add(
                mb.Port(anchor=second[0], orientation=-bond_vector, separation=0.075),
                label="down",
            )
            c2h2f2 = mb.Compound(subcompounds=(first, second))
            mb.force_overlap(first, first["up"], second["down"])
            fccf_dihedral_init = calc_dihedral(
                first[2].pos, first[0].pos, second[0].pos, second[2].pos
            )
            c2h2f2.remove_bond((first[0], second[0]))
            mb.force_overlap(first, first["port[0]"], second["port[0]"])
            fccf_dihedral_final = calc_dihedral(
                first[2].pos, first[0].pos, second[0].pos, second[2].pos
            )
            return fccf_dihedral_init, fccf_dihedral_final

        return _connect_and_reconnect

    @pytest.fixture
    def copper_cell(self):
        copper = mb.Compound(name="Cu")
        lattice_vector = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        spacing = [0.36149, 0.36149, 0.36149]
        copper_locations = [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ]
        basis = {"Cu": copper_locations}
        copper_lattice = mb.Lattice(
            lattice_spacing=spacing,
            lattice_vectors=lattice_vector,
            lattice_points=basis,
        )
        copper_dict = {"Cu": copper}
        copper_pillar = copper_lattice.populate(
            x=3, y=3, z=1, compound_dict=copper_dict
        )
        return copper_pillar

    @pytest.fixture
    def graphene(self):
        carbon = mb.Compound(name="C")
        angles = [90, 90, 120]
        carbon_locations = [[0, 0, 0], [2 / 3, 1 / 3, 0]]
        basis = {"C": carbon_locations}
        graphene = mb.Lattice(
            lattice_spacing=[0.2456, 0.2456, 0],
            angles=angles,
            lattice_points=basis,
        )
        carbon_dict = {"C": carbon}
        graphene_cell = graphene.populate(compound_dict=carbon_dict, x=3, y=3, z=1)
        return graphene_cell

    @pytest.fixture
    def cscl_crystal(self):
        cesium = mb.Compound(name="Cs")
        chlorine = mb.Compound(name="Cl")
        spacing = [0.4123, 0.4123, 0.4123]
        basis = {"Cs": [[0.5, 0.5, 0.5]], "Cl": [[0, 0, 0]]}
        cscl_lattice = mb.Lattice(spacing, lattice_points=basis)

        cscl_dict = {"Cs": cesium, "Cl": chlorine}
        cscl_compound = cscl_lattice.populate(x=3, y=3, z=1, compound_dict=cscl_dict)
        return cscl_compound

    @pytest.fixture
    def gilmerite(self):
        """Structure taken from:

        Sarp H., Cerny R., European Journal of Mineralogy, 11 (1999) p.549-555,
        Gilmarite, Cu3(AsO4)(OH)3, a new mineral:, its description and crystal
        structure.
        """
        gilmerite = mb.Compound()
        gilmerite.box = mb.Box(
            lengths=[5.44500017, 5.87300015, 5.10400009],
            angles=[114.94999695, 93.05000305, 91.91999817],
        )
        gilmerite.add(
            mb.Particle(name="As", pos=[5.43783569e-01, 1.54457900e-04, 4.61488000e-05])
        )
        gilmerite.add(mb.Particle(name="Cu", pos=[0.00615697, 0.28454988, 0.14878373]))
        gilmerite.add(mb.Particle(name="Cu", pos=[0.28204174, 0.13594167, 0.16465892]))
        gilmerite.add(mb.Particle(name="Cu", pos=[0.28377297, 0.43348074, 0.17854972]))
        gilmerite.add(mb.Particle(name="O", pos=[0.08094841, 0.14370937, 0.03276565]))
        gilmerite.add(mb.Particle(name="O", pos=[0.08451646, 0.45650253, 0.03415011]))
        gilmerite.add(mb.Particle(name="O", pos=[0.36793642, 0.28482277, 0.07106915]))
        gilmerite.add(mb.Particle(name="O", pos=[0.39833167, -0.00501952, 0.08583677]))
        gilmerite.add(mb.Particle(name="O", pos=[0.17415644, 0.27864442, 0.24828055]))
        gilmerite.add(mb.Particle(name="O", pos=[0.24149659, -0.00539473, 0.29073744]))
        gilmerite.add(mb.Particle(name="O", pos=[0.45795937, 0.36822546, 0.30135167]))
        return gilmerite
