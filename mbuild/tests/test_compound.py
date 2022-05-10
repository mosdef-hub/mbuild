import os
import time

import numpy as np
import parmed as pmd
import pytest

import mbuild as mb
from mbuild import Box, Compound, Particle, Port
from mbuild.exceptions import MBuildError
from mbuild.tests.base_test import BaseTest
from mbuild.utils.exceptions import RemovedFuncError
from mbuild.utils.geometry import calc_dihedral
from mbuild.utils.io import (
    get_fn,
    has_foyer,
    has_freud,
    has_intermol,
    has_mdtraj,
    has_networkx,
    has_openbabel,
    has_py3Dmol,
    has_rdkit,
    import_,
)

try:
    import nglview

    has_nglview = True
    del nglview
except ImportError:
    has_nglview = False


class TestCompound(BaseTest):
    def test_load_and_create(self):
        mb.load(get_fn("methyl.pdb"))

    def test_load_conversion(self, ethane, h2o):
        compound = Compound([ethane, h2o])
        parm = compound.to_parmed()
        traj = compound.to_trajectory()
        belmol = compound.to_pybel()

        for topo in [compound, parm, traj, belmol]:
            topo_converted = mb.load(topo)
            assert isinstance(topo_converted, Compound)
            assert topo_converted.n_particles == 11
            assert (
                len([at for at in topo_converted.particles() if at.name == "C"])
                == 2
            )
            assert (
                len([at for at in topo_converted.particles() if at.name == "H"])
                == 8
            )
            assert (
                len([at for at in topo_converted.particles() if at.name == "O"])
                == 1
            )

        for topo in [parm, traj]:
            new_topo = mb.load(compound)
            new_topo.xyz = np.random.random(topo_converted.xyz.shape)
            new_topo = mb.load(topo, compound=new_topo, coords_only=True)
            assert np.allclose(mb.load(topo).xyz, new_topo.xyz)

        # Extra test
        test = pmd.load_file(get_fn("styrene.mol2"), structure=True)
        assert isinstance(test, pmd.Structure)
        test_converted1 = mb.load(test)
        test_converted2 = Compound()
        test_converted2.from_parmed(test)

        assert isinstance(test_converted1, Compound)
        assert test_converted1.n_particles == len(test.atoms)
        assert test_converted2.n_particles == test_converted1.n_particles
        assert test_converted1.n_bonds == len(test.bonds)
        assert test_converted2.n_bonds == test_converted2.n_bonds

        test_converted1.xyz = np.random.random(test_converted1.xyz.shape)
        test_converted1 = mb.load(
            test, compound=test_converted1, coords_only=True
        )
        test_converted2.xyz = np.random.random(test_converted2.xyz.shape)
        test_converted2.from_parmed(test, coords_only=True)
        assert np.allclose(test_converted1.xyz, test_converted2.xyz)

    def test_load_xyz(self):
        class MyCompound(Compound):
            def __init__(self):
                super(MyCompound, self).__init__()

                mb.load(get_fn("ethane.xyz"), compound=self)

        myethane = MyCompound()
        assert myethane.n_particles == 8

    def test_update_from_file(self, ch3):
        ch3.update_coordinates(get_fn("methyl.pdb"))

    def test_load_protein(self):
        # Testing the loading function with complicated protein,
        # The protein file is taken from RCSB protein data bank
        # https://www.rcsb.org/structure/6M03
        protein = mb.load(get_fn("6m03.pdb"))

        # Asserting the protein having correct number of chains
        if has_mdtraj:
            # Loading using mdtraj
            assert len(protein.children) == 2

            chain1 = protein.children[0]
            chain2 = protein.children[1]  # Guess what this is: CG Water

            # Asserting the main chains having the correct number of residues
            # as well number of atoms and bonds
            # Main protein chain
            assert len(chain1.children) == 306
            assert chain1.n_particles == 2367
            assert chain1.n_bonds == 2420

            # CG Water
            assert len(chain2.children) == 87
            assert chain2.n_particles == 87
            assert chain2.n_bonds == 0
        else:
            # Loading using parmed
            # Chains info is lossed
            assert len(protein.children) == 393

    def test_save_simple(self, ch3):
        extensions = [".xyz", ".pdb", ".mol2", ".json", ".sdf"]
        for ext in extensions:
            outfile = "methyl_out" + ext
            ch3.save(filename=outfile)
            assert os.path.exists(outfile)

    def test_save_json_loop(self, ethane):
        ethane.save("ethane.json", show_ports=True)
        ethane_copy = mb.load("ethane.json")
        assert ethane.n_particles == ethane_copy.n_particles
        assert ethane.n_bonds == ethane_copy.n_bonds
        assert len(ethane.children) == len(ethane_copy.children)

    def test_save_box(self, ch3):
        extensions = [".mol2", ".pdb", ".hoomdxml", ".gro", ".sdf"]
        box_attributes = ["lengths"]
        custom_box = Box(lengths=[0.8, 0.8, 0.8], angles=[90, 90, 90])
        for ext in extensions:
            outfile_padded = "padded_methyl" + ext
            outfile_custom = "custom_methyl" + ext
            ch3.save(filename=outfile_padded, box=None, overwrite=True)
            ch3.save(filename=outfile_custom, box=custom_box, overwrite=True)
            padded_ch3 = mb.load(outfile_padded)
            custom_ch3 = mb.load(outfile_custom)
            for attr in box_attributes:
                pad_attr = getattr(padded_ch3.get_boundingbox(), attr)
                custom_attr = getattr(custom_ch3.get_boundingbox(), attr)
                assert np.array_equal(pad_attr, custom_attr)

    def test_save_overwrite(self, ch3):
        extensions = [".gsd", ".hoomdxml", ".lammps", ".lmp", ".top", ".gro"]
        for ext in extensions:
            outfile = "lyhtem" + ext
            ch3.save(filename=outfile)
            ch3.save(filename=outfile, overwrite=True)
            with pytest.raises(IOError):
                ch3.save(filename=outfile, overwrite=False)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield(self, methane):
        exts = [
            ".gsd",
            ".hoomdxml",
            ".lammps",
            ".lmp",
            ".top",
            ".gro",
            ".mol2",
            ".pdb",
            ".xyz",
            ".sdf",
        ]
        for ext in exts:
            methane.save(
                "lythem" + ext, forcefield_name="oplsaa", overwrite=True
            )

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield_with_file(self, methane):
        exts = [
            ".gsd",
            ".hoomdxml",
            ".lammps",
            ".lmp",
            ".top",
            ".gro",
            ".mol2",
            ".pdb",
            ".xyz",
            ".sdf",
        ]
        for ext in exts:
            methane.save(
                "lythem" + ext,
                forcefield_files=get_fn("methane_oplssaa.xml"),
                overwrite=True,
            )

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    @pytest.mark.parametrize(
        "ff_filename,foyer_kwargs",
        [
            ("ethane-angle-typo.xml", {"assert_angle_params": False}),
            ("ethane-dihedral-typo.xml", {"assert_dihedral_params": False}),
        ],
    )
    def test_save_missing_topo_params(self, ff_filename, foyer_kwargs):
        """Test that the user is notified if not all topology parameters are found."""
        from foyer.tests.utils import get_fn

        ethane = mb.load(get_fn("ethane.mol2"))
        with pytest.raises(Exception):
            ethane.save("ethane.mol2", forcefield_files=get_fn(ff_filename))
        with pytest.warns(UserWarning):
            ethane.save(
                "ethane.mol2",
                forcefield_files=get_fn(ff_filename),
                overwrite=True,
                foyer_kwargs=foyer_kwargs,
            )

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_forcefield_with_file_foyer_kwargs(self, methane):
        foyer_kwargs = {"assert_improper_params": True}
        with pytest.raises(Exception):
            methane.save(
                "lythem.hoomdxml",
                forcefield_files=get_fn("methane_oplssaa.xml"),
                overwrite=True,
                foyer_kwargs=foyer_kwargs,
            )
        methane.save(
            "lythem.hoomdxml",
            forcefield_files=get_fn("methane_oplssaa.xml"),
            overwrite=True,
            foyer_kwargs={},
        )

    def test_save_resnames(self, ch3, h2o):
        system = Compound([ch3, h2o])
        system.save("resnames.gro", residues=["CH3", "H2O"])
        struct = pmd.load_file("resnames.gro")

        assert struct.residues[0].name == "CH3"
        assert struct.residues[1].name == "H2O"

    def test_save_resnames_single(self, c3, n4):
        system = Compound([c3, n4])
        system.save("resnames_single.gro", residues=["C3", "N4"])
        struct = pmd.load_file("resnames_single.gro")
        assert struct.residues[0].number == 1
        assert struct.residues[1].number == 2

    def test_save_residue_map(self, methane):
        filled = mb.fill_box(methane, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        t0 = time.time()
        filled.save("filled.mol2", forcefield_name="oplsaa", residues="Methane")
        t1 = time.time()
        foyer_kwargs = {"use_residue_map": False}
        filled.save(
            "filled.mol2",
            forcefield_name="oplsaa",
            overwrite=True,
            residues="Methane",
            foyer_kwargs=foyer_kwargs,
        )
        t2 = time.time()
        assert (t2 - t1) > (t1 - t0)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_references(self, methane):
        foyer_kwargs = {"references_file": "methane.bib"}
        methane.save(
            "methyl.mol2", forcefield_name="oplsaa", foyer_kwargs=foyer_kwargs
        )
        assert os.path.isfile("methane.bib")

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_save_combining_rule(self, methane):
        combining_rules = ["lorentz", "geometric"]
        gmx_rules = {"lorentz": 2, "geometric": 3}
        for combining_rule in combining_rules:
            if combining_rule == "geometric":
                methane.save(
                    "methane.top",
                    forcefield_name="oplsaa",
                    combining_rule=combining_rule,
                    overwrite=True,
                )
            else:
                with pytest.warns(UserWarning):
                    methane.save(
                        "methane.top",
                        forcefield_name="oplsaa",
                        combining_rule=combining_rule,
                        overwrite=True,
                    )
            with open("methane.top") as fp:
                for i, line in enumerate(fp):
                    if i == 18:
                        gmx_rule = int(line.split()[1])
                        assert gmx_rule == gmx_rules[combining_rule]

    def test_clone_with_box(self, ethane):
        ethane.box = ethane.get_boundingbox()
        ethane.periodicity = (True, True, False)
        ethane_clone = mb.clone(ethane)
        assert np.all(ethane.xyz == ethane_clone.xyz)
        assert np.all(
            [p.name for p in ethane.particles()]
            == [p.name for p in ethane_clone.particles()]
        )
        assert len(ethane.children) == len(ethane_clone.children)
        assert ethane_clone.mass == ethane.mass
        assert ethane.periodicity == ethane_clone.periodicity

    def test_batch_add(self, ethane, h2o):
        compound = Compound()
        compound.add([ethane, h2o])
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds1(self, ethane):
        compound = Compound(ethane)
        assert compound.n_particles == 8
        assert compound.n_bonds == 7

    def test_init_with_subcompounds2(self, ethane, h2o):
        compound = Compound([ethane, h2o])
        assert compound.n_particles == 8 + 3
        assert compound.n_bonds == 7 + 2

    def test_init_with_subcompounds3(self, ethane, h2o):
        compound = Compound([ethane, [h2o, mb.clone(h2o)]])
        assert compound.n_particles == 8 + 2 * 3
        assert compound.n_bonds == 7 + 2 * 2

    def test_init_mass(self):
        element = mb.Compound(name="carbon", element="C")
        assert np.allclose(12.011, element.mass, atol=1e-5)

        bead = mb.Compound(name="A", mass=1.0)
        assert bead.mass == 1.0

        bead_overwrite = mb.Compound(name="A", element="C", mass=1.0)
        assert bead_overwrite.mass == 1.0

        bead_no_mass = mb.Compound(name="A")
        assert bead_no_mass.mass == 0.0

    def test_init_with_bad_mass(self):
        with pytest.raises(MBuildError):
            a = mb.Compound(name="A")
            b = mb.Compound(name="B")
            compound = mb.Compound(subcompounds=[a, b], mass=2.0)

        with pytest.raises(ValueError):
            mb.Compound(name="A", mass=-1.0)

    def test_init_with_bad_name(self):
        with pytest.raises(ValueError):
            Compound(name=1)

    def test_add_wrong_input(self, ethane):
        with pytest.raises(ValueError):
            ethane.add("water")

    def test_mass_property(self, h2o):
        methane = mb.load("C", smiles=True)
        assert np.allclose(methane.mass, 16.043, atol=1e-5)

        assert np.allclose(h2o.mass, 18.015, atol=1e-5)

        system = mb.fill_box(compound=h2o, n_compounds=5, box=[0.5, 0.5, 0.5])
        assert np.allclose(system.mass, 5 * h2o.mass, atol=1e-5)

    def test_mass_setter(self, ethane):
        comp = mb.Compound(name="A", mass=1.0)
        comp.mass = 2.0
        assert comp.mass == 2.0

        carbon = mb.Compound(name="carbon", element="C")
        carbon.mass = 1.0
        assert carbon.mass == 1.0

        carbon.mass = np.array(12.011)
        assert carbon.mass == 12.011

        carbon.mass = "1.0"
        assert carbon.mass == 1.0

        for p in ethane.particles():
            p.mass = 1.0
        assert ethane.mass == 8.0

    def test_mass_setter_wrong_input(self, methane):
        with pytest.raises(MBuildError):
            methane.mass = 1.0

        with pytest.raises(ValueError):
            methane[0].mass = -1.0

    def test_reset_parent_mass(self, methane):
        A = mb.Compound()
        A.mass = 1.0
        A.add(methane)
        assert A._mass == 0.0
        assert A.mass == methane.mass

    def test_mass_add_port(self):
        A = mb.Compound(mass=2.0)
        A.add(mb.Port())
        assert A.mass == 2.0

    def test_add_existing_parent(self, ethane, h2o):
        water_in_water = mb.clone(h2o)
        h2o.add(water_in_water)
        with pytest.raises(MBuildError):
            ethane.add(water_in_water)

    def test_add_label_exists(self, ethane, h2o):
        ethane.add(h2o, label="water")
        with pytest.raises(MBuildError):
            ethane.add(mb.clone(h2o), label="water")

    def test_set_pos(self, ethane):
        with pytest.raises(MBuildError):
            ethane.pos = [0, 0, 0]

    def test_xyz(self, ch3):
        xyz = ch3.xyz
        assert xyz.shape == (4, 3)

        xyz = ch3.xyz_with_ports
        assert xyz.shape == (12, 3)

    def test_xyz_setter_bad_shape(self):
        single_compound = Compound()
        with pytest.raises(ValueError):
            single_compound.xyz = np.zeros(shape=(4, 10))
        with pytest.raises(ValueError):
            single_compound.xyz_with_ports = np.zeros(shape=(4, 10))

    def test_particles_by_name(self, ethane):
        assert sum(1 for _ in ethane.particles()) == 8

        only_H = ethane.particles_by_name("H")
        assert sum(1 for _ in only_H) == 6

        only_C = ethane.particles_by_name("C")
        assert sum(1 for _ in only_C) == 2

    def test_particles_in_range(self, ethane):
        group = ethane.particles_in_range(ethane[0], 0.141)
        assert sum([1 for x in group if x.name == "H"]) == 3
        assert sum([1 for x in group if x.name == "C"]) == 2

        group = ethane.particles_in_range(ethane[0], 0.141, max_particles=4)
        assert sum([1 for x in group if x.name == "H"]) == 3
        assert sum([1 for x in group if x.name == "C"]) == 1

    def test_generate_bonds(self, ch3):
        ch3.generate_bonds("H", "H", dmin=0.01, dmax=2.0)
        assert ch3.n_bonds == 3 + 3

    @pytest.mark.skipif(not has_freud, reason="Freud not installed.")
    def test_freud_generated_bonds_periodicity(self, ch3):
        bounding_box = ch3.get_boundingbox()

        ch3_clone = mb.clone(ch3)
        ch3_clone.box = mb.Box(lengths=[max(bounding_box.lengths) + 1] * 3)
        ch3_clone.periodicity = (True, True, True)
        ch3_clone.freud_generate_bonds(
            "H", "H", dmin=0.01, dmax=0.2, exclude_ii=True
        )
        assert ch3_clone.n_bonds == 3 + 3

        ch3_clone2 = mb.clone(ch3)
        ch3_clone2.box = mb.Box(lengths=[max(bounding_box.lengths) + 1] * 3)
        ch3_clone2.periodicity = (True, True, False)
        ch3_clone2.freud_generate_bonds(
            "H", "H", dmin=0.01, dmax=0.2, exclude_ii=True
        )
        assert ch3_clone2.n_bonds == 3 + 3

    @pytest.mark.skipif(not has_freud, reason="Freud not installed.")
    def test_freud_generate_bonds(self, ch3):
        bounding_box = ch3.get_boundingbox()
        ch3.box = mb.Box(lengths=[max(bounding_box.lengths) + 1] * 3)
        ch3.freud_generate_bonds("H", "H", dmin=0.01, dmax=0.2, exclude_ii=True)
        assert ch3.n_bonds == 3 + 3

    @pytest.mark.skipif(not has_freud, reason="Freud not installed.")
    def test_freud_generate_bonds_expected(self, ch3):
        bounding_box = ch3.get_boundingbox()
        ch3.box = mb.Box(lengths=[max(bounding_box.lengths) + 1] * 3)
        ch3.freud_generate_bonds("H", "H", dmin=0.01, dmax=0.1, exclude_ii=True)
        assert ch3.n_bonds == 3

    def test_remove_from_box(self, ethane):
        n_ethanes = 5
        box = mb.fill_box(ethane, n_ethanes, [3, 3, 3])
        box.remove(box.children[3])

        n_ethanes -= 1
        assert box.n_particles == n_ethanes * ethane.n_particles
        assert len(box.children) == n_ethanes
        assert box.n_bonds == n_ethanes * ethane.n_bonds
        assert (
            len(
                [
                    meth.referenced_ports()
                    for eth in box.children
                    for meth in eth.children
                ]
            )
            == 2 * n_ethanes
        )

    def test_remove(self, ethane):
        # create and remove a subcompound

        ethane1 = mb.clone(ethane)
        hydrogens = ethane1.particles_by_name("H")
        ethane1.remove(hydrogens)

        assert ethane1.n_particles == 2
        assert ethane1.n_bonds == 1
        for part in ethane1.children:
            assert part.n_bonds == 0
            assert part.n_particles == 1
            assert len(part.children) == 4
        assert len(ethane1.children) == 2

        carbons = ethane1.particles_by_name("C")
        ethane1.remove(carbons)
        assert ethane1.n_particles == 1  # left with the highest Compound
        assert ethane1.n_bonds == 0
        assert len(ethane1.children) == 0  # left with highest Compound

        # Test remove all particles belong to a single child of an Ethane
        ethane2 = mb.clone(ethane)
        CH3_particles = list(ethane2.children[0].particles())
        ethane2.remove(CH3_particles)
        assert len(ethane2.children) == 1
        assert len(ethane2.children[0].children) == 5  # 4 particles + 1 port

        # Test remove a subcompound
        ethane3 = mb.clone(ethane)
        ethane3.remove(ethane3.children[0])
        assert len(ethane3.children) == 1
        assert len(ethane3.children[0].children) == 5  # 4 particles + 1 port

        # Test remove an entire compound
        ethane4 = mb.clone(ethane)
        ethane4.remove(ethane4)
        assert ethane4.n_particles == 1  # left with the highest Compound
        assert ethane4.n_bonds == 0
        assert len(ethane4.children) == 0  # left with highest Compound

        # Test remove one subcompound and part of another
        ethane5 = mb.clone(ethane)
        ethane5.remove(
            [particle for particle in ethane5.children[0].particles()]
            + [ethane5.children[1].children[0]]
        )
        assert ethane5.n_particles == 3  # three hydrogens
        assert ethane5.n_bonds == 0
        assert len(ethane5.children[0].children) == 6  # 3 hydrogens + 3 ports
        assert len(ethane5.children) == 1

    def test_remove_many(self, ethane):
        ethane.remove([ethane.children[0], ethane.children[1]])

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        for part in ethane.children:
            assert isinstance(part, Port)

    def test_remove_subcompound(self, ethane):
        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 4
        assert ethane.n_bonds == 3
        assert len(ethane.children) == 1
        # Still contains a port
        assert len(ethane.children[0].children) == 5

        methyl = ethane.children[0]
        ethane.remove(methyl)

        assert ethane.n_particles == 1
        assert ethane._n_particles() == 0
        assert ethane.n_bonds == 0
        assert len(ethane.children) == 0

    def test_remove_no_bond_graph(self):
        compound = Compound()
        particle = Compound(name="C", pos=[0, 0, 0])
        compound.add(particle, "test-particle")
        compound.remove(particle)
        assert particle not in compound.particles()

    def test_remove_bond(self, ch3):
        ch_bond = list(ch3.bonds())[0]
        ch3.remove_bond(ch_bond)
        assert ch3.n_bonds == 2

        with pytest.warns(UserWarning):
            ch3.remove_bond(ch_bond)

    def test_port_does_not_exist(self, ethane):
        with pytest.raises(MBuildError):
            ethane["not_port"]

    def test_center(self, methane):
        assert np.array_equal(methane.center, np.array([0, 0, 0]))
        for orientation in np.identity(3):
            separation = 0.2
            port = Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(
                port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15
            )
            port = Port(
                anchor=methane[0],
                orientation=orientation,
                separation=separation,
            )
            assert np.allclose(
                port.center, separation * orientation, atol=1e-15
            )
        np.random.seed(0)
        for orientation in np.random.rand(5, 3):
            port = Port(anchor=methane[0], orientation=orientation)
            assert np.allclose(
                port.center, np.array([0.0, 0.0, 0.0]), atol=1e-15
            )
            port = Port(
                anchor=methane[0],
                orientation=orientation,
                separation=separation,
            )
            assert np.allclose(
                port.center,
                separation * orientation / np.linalg.norm(orientation),
                atol=1e-15,
            )

    def test_single_particle(self):
        part = Particle(name="A")
        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == part
        assert len(list(part.ancestors())) == 0
        assert next(part.particles_by_name("A")) == part

    def test_name(self):
        with pytest.raises(ValueError):
            Compound(name=1)

    def test_particle_in_particle(self):
        part = Particle(name="A")
        parent = Compound(part)

        assert part.n_particles == 1
        assert len(list(part.particles())) == 1
        assert part.xyz.shape == (1, 3)
        assert part.root == parent
        assert len(list(part.ancestors())) == 1
        assert next(part.particles_by_name("A")) == part

        assert parent.n_particles == 1
        assert len(list(parent.particles())) == 1
        assert parent.xyz.shape == (1, 3)
        assert parent.root == parent
        assert len(list(parent.ancestors())) == 0
        assert next(parent.particles_by_name("A")) == part

    def test_flatten_eth(self, ethane):
        # Before flattening
        assert len(ethane.children) == 2
        assert ethane.n_particles == 8
        assert ethane.n_bonds == 7

        # Flatten with inplace = False
        copy = ethane.flatten(inplace=False)
        assert ethane.n_particles == copy.n_particles == len(copy.children)
        assert ethane.n_bonds == copy.n_bonds

        # After flattening
        ethane.flatten()
        assert len(ethane.children) == ethane.n_particles == 8
        assert ethane.n_bonds == 7

    def test_flatten_box_of_eth(self, ethane):
        box_of_eth = mb.fill_box(compound=ethane, n_compounds=2, box=[1, 1, 1])
        # Before flattening
        assert len(box_of_eth.children) == 2
        assert box_of_eth.n_bonds == 7 * 2
        assert box_of_eth.n_particles == 8 * 2

        # After flattening
        box_of_eth.flatten()
        assert len(box_of_eth.children) == box_of_eth.n_particles == 8 * 2
        assert box_of_eth.n_bonds == 7 * 2

    def test_flatten_with_port(self, ethane):
        ethane.remove(ethane[2])
        original_ports = ethane.all_ports()
        ethane.flatten()
        assert len(ethane.all_ports()) == len(original_ports)
        assert ethane.all_ports()[0] == original_ports[0]

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_reload(self):
        # Create a compound and write it to file.
        p3ht1 = mb.load("CCCCCCC1=C(SC(=C1)C)C", smiles=True, backend="pybel")
        p3ht1.save("p3ht1.pdb")

        # Create another compound, rotate it and write it to file.
        p3ht2 = mb.load("CCCCCCC1=C(SC(=C1)C)C", smiles=True, backend="pybel")
        p3ht2.rotate(np.pi / 2, [0, 0, 1])
        p3ht2.save("p3ht2.pdb")

        # Load p3ht2.pdb into p3ht1, modifying the atom positions of p3ht1.
        p3ht1.update_coordinates("p3ht2.pdb")
        p3ht1.save("modified_p3ht1.pdb")

        assert p3ht1.n_particles == 33
        assert p3ht1.n_bonds == 33

    @pytest.mark.parametrize(
        "extension", [(".xyz"), (".pdb"), (".mol2"), (".gro")]
    )
    def test_update_coordinates(self, ethane, extension):
        ethane_clone = mb.clone(ethane)
        ethane_clone.xyz += [1, 1, 1]

        fn = "ethane_clone" + extension
        ethane_clone.save(fn)
        ethane.update_coordinates(fn)

        new_file = mb.load(fn)
        assert np.allclose(ethane.xyz, ethane_clone.xyz, atol=1e-3)
        assert np.allclose(ethane.xyz, new_file.xyz)

    def test_update_coordinates_no_hierarchy(self):
        mycomp = Compound()
        myclone = mb.clone(mycomp)
        myclone.xyz += 1

        myclone.save("myclone.pdb", overwrite=True)

        assert np.allclose(mycomp.xyz, np.array([0, 0, 0]))
        mycomp.update_coordinates("myclone.pdb")
        assert np.allclose(mycomp.xyz, np.array([1, 1, 1]))
        ref = mb.load("myclone.pdb")
        assert np.allclose(mycomp.xyz, ref.xyz)

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_to_trajectory(self, ethane, c3, n4):
        traj = ethane.to_trajectory()
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 1

        traj = ethane.to_trajectory(residues="CH3")
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 1
        assert traj.n_residues == 2
        assert "CH3" in [res.name for res in traj.top.residues]
        assert all(res.n_atoms == 4 for res in traj.top.residues)

        traj = ethane.to_trajectory(chains="CH3")
        assert traj.n_atoms == 8
        assert traj.top.n_bonds == 7
        assert traj.n_chains == 2
        assert traj.n_residues == 2
        assert all(chain.n_atoms == 4 for chain in traj.top.chains)
        assert all(chain.n_residues == 1 for chain in traj.top.chains)

        system = Compound([c3, n4])
        traj = system.to_trajectory(residues=["C", "N"])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 1
        assert traj.n_residues == 2

        traj = system.to_trajectory(chains=["C", "N"])
        assert traj.n_atoms == 2
        assert traj.top.n_bonds == 0
        assert traj.n_chains == 2
        assert traj.n_residues == 2

        methyl = next(iter(ethane.children))
        traj = methyl.to_trajectory()
        assert traj.n_atoms == 4
        assert traj.top.n_bonds == 3
        assert traj.n_chains == 1
        assert traj.n_residues == 1

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_box_mdtraj(self, ethane):
        box = Box(lengths=[4.0, 5.0, 6.0], angles=[90, 90, 90])
        assert ethane.box is None
        # (1) Specify box
        traj = ethane.to_trajectory(box=box)
        assert np.allclose(traj.unitcell_lengths, [4.0, 5.0, 6.0])
        # (2) Extract from bounding box
        assert ethane.box is None
        traj = ethane.to_trajectory()
        assert np.allclose(
            # Comparison in nm
            traj.unitcell_lengths,
            np.array(ethane.get_boundingbox().lengths) + 0.5,
        )
        # (3) Box as part of compound
        ethane.box = Box.from_mins_maxs_angles(
            mins=np.zeros(3), maxs=8.0 * np.ones(3), angles=[90.0, 90.0, 90.0]
        )
        traj = ethane.to_trajectory()
        assert np.allclose(traj.unitcell_lengths, np.asarray([8, 8, 8]))

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_resnames_mdtraj(self, h2o, ethane):
        system = Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(residues=["Ethane", "H2O"])
        residues = list(traj.top.residues)
        assert traj.n_residues == 3
        assert residues[0].name == "H2O"
        assert residues[1].name == "H2O"
        assert residues[2].name == "Ethane"

        traj = system.to_trajectory(residues="Ethane")
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == "RES"
        assert residues[1].name == "Ethane"

        traj = system.to_trajectory(residues=["Ethane"])
        residues = list(traj.top.residues)
        assert traj.n_residues == 2
        assert residues[0].name == "RES"
        assert residues[1].name == "Ethane"

        traj = system.to_trajectory()
        residues = list(traj.top.residues)
        assert traj.n_residues == 1
        assert residues[0].name == "RES"

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_chainnames_mdtraj(self, h2o, ethane):
        system = Compound([h2o, mb.clone(h2o), ethane])
        traj = system.to_trajectory(chains=["Ethane", "H2O"])
        assert traj.n_chains == 3

        traj = system.to_trajectory(chains="Ethane")
        assert traj.n_chains == 2

        traj = system.to_trajectory(chains=["Ethane"])
        assert traj.n_chains == 2

        traj = system.to_trajectory()
        assert traj.n_chains == 1

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_mdtraj_box(self, h2o):
        compound = Compound()
        compound.add(h2o)
        tilted_box = Box(lengths=[2.0, 2.0, 2.0], angles=[60.0, 80.0, 100.0])
        trajectory = compound.to_trajectory(box=tilted_box)
        assert np.all(np.isclose(trajectory.unitcell_lengths, [2.0, 2.0, 2.0]))
        assert np.all(
            np.isclose(trajectory.unitcell_angles, [60.0, 80.0, 100.0])
        )

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion1(self, ethane, h2o):
        compound = Compound([ethane, h2o])

        intermol_system = compound.to_intermol()
        assert len(intermol_system.molecule_types) == 1
        assert "Compound" in intermol_system.molecule_types
        assert len(intermol_system.molecule_types["Compound"].bond_forces) == 9

        assert len(intermol_system.molecule_types["Compound"].molecules) == 1
        molecules = list(intermol_system.molecule_types["Compound"].molecules)
        assert len(molecules[0].atoms) == 11

    @pytest.mark.skipif(not has_intermol, reason="InterMol is not installed")
    def test_intermol_conversion2(self, ethane, h2o):
        # 2 distinct Ethane objects.
        compound = Compound([ethane, mb.clone(ethane), h2o])

        molecule_types = [ethane.name, h2o.name]
        intermol_system = compound.to_intermol(molecule_types=molecule_types)
        assert len(intermol_system.molecule_types) == 2
        assert "Ethane" in intermol_system.molecule_types
        assert "H2O" in intermol_system.molecule_types

        assert len(intermol_system.molecule_types["Ethane"].bond_forces) == 7
        assert len(intermol_system.molecule_types["H2O"].bond_forces) == 2

        assert len(intermol_system.molecule_types["Ethane"].molecules) == 2
        ethanes = list(intermol_system.molecule_types["Ethane"].molecules)
        assert len(ethanes[0].atoms) == len(ethanes[1].atoms) == 8

        assert len(intermol_system.molecule_types["H2O"].molecules) == 1
        h2os = list(intermol_system.molecule_types["H2O"].molecules)
        assert len(h2os[0].atoms) == 3

    def test_parmed_conversion(self, ethane, h2o):
        compound = Compound([ethane, h2o])

        structure = compound.to_parmed()
        assert structure.title == "Compound"

        structure = compound.to_parmed(title="eth_h2o")
        assert structure.title == "eth_h2o"

        assert len(structure.atoms) == 11
        assert len([at for at in structure.atoms if at.element == 6]) == 2
        assert len([at for at in structure.atoms if at.element == 1]) == 8
        assert len([at for at in structure.atoms if at.element == 8]) == 1

        assert len(structure.bonds) == 9

        assert sum(len(res.atoms) for res in structure.residues) == len(
            structure.atoms
        )

        compound2 = Compound()
        compound2.from_parmed(structure)

        assert compound2.n_particles == 11
        assert len([at for at in compound2.particles() if at.name == "C"]) == 2
        assert len([at for at in compound2.particles() if at.name == "H"]) == 8
        assert len([at for at in compound2.particles() if at.name == "O"]) == 1

        assert compound2.n_bonds == 9

        compound3 = mb.clone(compound2)
        compound3.xyz = np.random.random(compound3.xyz.shape)
        compound3.from_parmed(structure, coords_only=True)

        assert np.allclose(compound2.xyz, compound3.xyz)

    def test_box_parmed(self, ethane):
        box = Box(lengths=[4.0, 5.0, 6.0], angles=[90, 90, 90])
        assert ethane.box is None
        # (1) Specify box
        pmd = ethane.to_parmed(box=box)
        assert np.allclose(
            # Comparison in nm
            0.1 * pmd.box[0:3],
            box.lengths,
        )

        # (2) Extract from bounding box
        assert ethane.box is None
        pmd = ethane.to_parmed()
        assert np.allclose(
            # Comparison in nm
            0.1 * pmd.box[0:3],
            np.array(ethane.get_boundingbox().lengths) + 0.5,
        )

        # (3) Extract from self.box
        ethane.box = box
        pmd = ethane.to_parmed()
        assert ethane.box is not None
        assert np.allclose(
            # Comparison in nm
            0.1 * pmd.box[0:3],
            ethane.box.lengths,
        )

    def test_fillbox_then_parmed(self):
        # This test would fail with the old to_parmed code (pre PR #699)

        bead = Compound(name="Bead")
        box = Box.from_mins_maxs_angles(
            mins=(2, 2, 2), maxs=(3, 3, 3), angles=(90.0, 90.0, 90.0)
        )
        bead_box = mb.fill_box(bead, 100, box=[2, 2, 2, 3, 3, 3])
        bead_box_in_pmd = bead_box.to_parmed(box=box)

        assert isinstance(bead_box_in_pmd, pmd.Structure)
        assert len(bead_box_in_pmd.atoms) == 100
        assert (
            bead_box_in_pmd.box
            == np.array([10.0, 10.0, 10.0, 90.0, 90.0, 90.0])
        ).all()

    def test_resnames_parmed(self, h2o, ethane):
        system = Compound([h2o, mb.clone(h2o), ethane])
        struct = system.to_parmed(residues=["Ethane", "H2O"])
        assert len(struct.residues) == 3
        assert struct.residues[0].name == "H2O"
        assert struct.residues[1].name == "H2O"
        assert struct.residues[2].name == "Ethane"
        assert sum(len(res.atoms) for res in struct.residues) == len(
            struct.atoms
        )

        struct = system.to_parmed(residues="Ethane")
        assert len(struct.residues) == 2
        assert struct.residues[0].name == "RES"
        assert struct.residues[1].name == "Ethane"
        assert sum(len(res.atoms) for res in struct.residues) == len(
            struct.atoms
        )

        struct = system.to_parmed()
        assert len(struct.residues) == 1
        assert struct.residues[0].name == "RES"
        assert sum(len(res.atoms) for res in struct.residues) == len(
            struct.atoms
        )

        struct = system.to_parmed(infer_residues=True)
        assert len(struct.residues) == 3
        assert struct.residues[0].name == "H2O"
        assert struct.residues[1].name == "H2O"
        assert struct.residues[2].name == "Ethane"
        assert sum(len(res.atoms) for res in struct.residues) == len(
            struct.atoms
        )

    def test_parmed_element_guess(self):
        compound = Particle(name="foobar")
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

        compound = Particle(name="XXXXXX")
        with pytest.warns(UserWarning):
            _ = compound.to_parmed()

    def test_parmed_box(self, h2o):
        compound = Compound()
        compound.add(h2o)
        tilted_box = Box(lengths=[2.0, 2.0, 2.0], angles=[90.0, 90.0, 120.0])
        structure = compound.to_parmed(box=tilted_box)
        assert np.all(
            np.isclose(structure.box, [20, 20, 20, 90.0, 90.0, 120.0])
        )

    def test_parmed_box_with_periodicity(self, h2o):
        compound = Compound()
        compound.add(h2o)
        compound.box = Box(lengths=[2.0, 2.0, 2.0], angles=[90.0, 90.0, 120.0])
        structure = compound.to_parmed()
        assert np.all(
            np.isclose(structure.box, [20, 20, 20, 90.0, 90.0, 120.0])
        )

    def test_min_periodic_dist(self, ethane):
        compound = Compound(ethane)
        C_pos = np.array(
            [atom.pos for atom in list(compound.particles_by_name("C"))]
        )
        assert (
            round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.14
        )
        compound.box = Box(lengths=[0.2, 0.2, 0.2], angles=[90, 90, 90])
        assert (
            round(compound.min_periodic_distance(C_pos[0], C_pos[1]), 2) == 0.06
        )

    def test_bond_graph(self, ch3):
        compound = Compound()
        compound.add(ch3)
        assert compound.n_bonds == 3
        assert all(
            compound.bond_graph.has_node(particle)
            for particle in ch3.particles()
        )

        ch3_nobonds = mb.clone(ch3)
        for bond in ch3_nobonds.bonds():
            ch3_nobonds.remove_bond(bond)
        compound.add(ch3_nobonds)
        assert compound.n_bonds == 3
        assert not any(
            compound.bond_graph.has_node(particle)
            for particle in ch3_nobonds.particles()
        )

        carbons = list(compound.particles_by_name("C"))
        compound.add_bond((carbons[0], carbons[1]))
        assert compound.n_bonds == 4
        assert all(
            compound.bond_graph.has_node(particle) for particle in carbons
        )
        assert any(
            compound.bond_graph.has_node(particle)
            for particle in ch3_nobonds.particles()
        )

        compound.remove_bond((carbons[0], carbons[1]))
        assert not any(
            compound.bond_graph.has_node(particle)
            for particle in ch3_nobonds.particles()
        )

    def test_is_independent(self):
        from mbuild.lib.molecules import WaterSPC

        # Test the normal case
        box_of_water = mb.fill_box(WaterSPC(), n_compounds=100, box=[3, 3, 3])
        assert box_of_water.is_independent()

        for particle in box_of_water.particles():
            assert not particle.is_independent()

        for child in box_of_water.children:
            assert child.is_independent()

        # Test the case where there is no bond
        top_comp = mb.Compound(name="top")
        mid_comp = mb.Compound(name="mid")
        bot_comp = mb.Compound(name="bot")

        top_comp.add(mid_comp)
        mid_comp.add(bot_comp)

        assert top_comp.is_independent()
        assert mid_comp.is_independent()
        assert bot_comp.is_independent()

    def test_is_independent_progression(self):
        from mbuild.lib.moieties import CH3

        ch3_1 = CH3()
        ch3_2 = CH3()
        assert ch3_1.is_independent()
        assert ch3_2.is_independent()

        eth = mb.Compound()
        eth.add(ch3_1, "methyl1")
        eth.add(ch3_2, "methyl2")
        assert ch3_1.is_independent()
        assert ch3_2.is_independent()

        mb.force_overlap(
            eth["methyl1"], eth["methyl1"]["up"], eth["methyl2"]["up"]
        )
        assert not ch3_1.is_independent()
        assert not ch3_2.is_independent()

        eth.remove_bond(list(eth.particles_by_name("C")))
        assert ch3_1.is_independent()
        assert ch3_2.is_independent()

    def test_update_coords_update_ports(self, ch2):
        distances = np.round(
            [
                ch2.min_periodic_distance(port.pos, ch2[0].pos)
                for port in ch2.referenced_ports()
            ],
            5,
        )
        orientations = np.round(
            [port.pos - port.anchor.pos for port in ch2.referenced_ports()], 5
        )

        ch2_clone = mb.clone(ch2)
        ch2_clone[0].pos += [1, 1, 1]
        ch2_clone.save("ch2-shift.pdb")

        ch2.update_coordinates("ch2-shift.pdb")
        updated_distances = np.round(
            [
                ch2.min_periodic_distance(port.pos, ch2[0].pos)
                for port in ch2.referenced_ports()
            ],
            5,
        )
        updated_orientations = np.round(
            [port.pos - port.anchor.pos for port in ch2.referenced_ports()], 5
        )

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    def test_charge(self, ch2, ch3):
        compound = Compound(charge=2.0)
        assert compound.charge == 2.0
        compound2 = Compound()
        assert compound2.charge == 0.0

        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        ch3[0].charge = 1.0
        compound.add([ch2, ch3])
        assert compound.charge == 1.25
        assert ch2.charge == 0.25
        assert compound[0].charge == 0.5

        with pytest.raises(AttributeError):
            compound.charge = 2.0

    def test_charge_subcompounds(self, ch2, ch3):
        ch2[0].charge = 0.5
        ch2[1].charge = -0.25
        compound = Compound(subcompounds=ch2)
        assert compound.charge == 0.25

        with pytest.raises(MBuildError):
            compound = Compound(subcompounds=ch3, charge=1.0)

    def test_charge_neutrality_warn(self, benzene):
        benzene[0].charge = 0.25
        with pytest.warns(UserWarning):
            benzene.save("charge-test.mol2")

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_energy_minimize(self, octane):
        octane.energy_minimize()

    @pytest.mark.skipif(has_openbabel, reason="Open Babel package is installed")
    def test_energy_minimize_openbabel_warn(self, octane):
        with pytest.raises(MBuildError):
            octane.energy_minimize()

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_energy_minimize_ff(self, octane):
        for ff in ["UFF", "GAFF", "MMFF94", "MMFF94s", "Ghemical"]:
            octane.energy_minimize(forcefield=ff)
        with pytest.raises(IOError):
            octane.energy_minimize(forcefield="fakeFF")

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_energy_minimize_algorithm(self, octane):
        for algorithm in ["cg", "steep", "md"]:
            octane.energy_minimize(algorithm=algorithm)
        with pytest.raises(MBuildError):
            octane.energy_minimize(algorithm="fakeAlg")

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_energy_minimize_non_element(self, octane):
        for particle in octane.particles():
            particle.element = None
        # Pass with element inference from names
        octane.energy_minimize()
        for particle in octane.particles():
            particle.name = "Q"
        # Fail once names don't match elements
        with pytest.raises(MBuildError):
            octane.energy_minimize()

    @pytest.mark.skipif(
        not has_openbabel, reason="Open Babel package not installed"
    )
    def test_energy_minimize_ports(self, octane):
        distances = np.round(
            [
                octane.min_periodic_distance(port.pos, port.anchor.pos)
                for port in octane.all_ports()
            ],
            5,
        )
        orientations = np.round(
            [port.pos - port.anchor.pos for port in octane.all_ports()], 5
        )

        octane.energy_minimize()

        updated_distances = np.round(
            [
                octane.min_periodic_distance(port.pos, port.anchor.pos)
                for port in octane.all_ports()
            ],
            5,
        )
        updated_orientations = np.round(
            [port.pos - port.anchor.pos for port in octane.all_ports()], 5
        )

        assert np.array_equal(distances, updated_distances)
        assert np.array_equal(orientations, updated_orientations)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_energy_minimize_openmm(self, octane):
        octane.energy_minimize(forcefield="oplsaa")

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    @pytest.mark.parametrize(
        "constraints", ["AllBonds", "HBonds", "HAngles", None]
    )
    def test_energy_minimize_openmm_constraints(self, octane, constraints):
        octane.energy_minimize(forcefield="oplsaa", constraints=constraints)

    def test_energy_minimize_openmm_invalid_constraints(self, octane):
        with pytest.raises(ValueError):
            octane.energy_minimize(forcefield="oplsaa", constraints="boo")

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_energy_minimize_openmm_xml(self, octane):
        octane.energy_minimize(forcefield=get_fn("small_oplsaa.xml"))

    def test_clone_outside_containment(self, ch2, ch3):
        compound = Compound()
        compound.add(ch2)
        mb.force_overlap(ch3, ch3["up"], ch2["up"])
        with pytest.raises(MBuildError):
            ch3_clone = mb.clone(ch3)

    def test_load_nonelement_mol2(self):
        mb.load(get_fn("benzene-nonelement.mol2"), backend="mdtraj")
        mb.load(get_fn("benzene-nonelement.mol2"), backend="parmed")

    def test_load_nonatom_mdtraj_mol2(self):
        # First atom name and element are incorrect
        # Loading with MDTraj should raise an error
        with pytest.raises(KeyError):
            mb.load(get_fn("benzene-nonatom-nonelement.mol2"), backend="mdtraj")

    def test_siliane_bond_number(self, silane):
        assert silane.n_bonds == 4

    def test_add_bond_remove_ports(self, hydrogen):
        h_clone = mb.clone(hydrogen)
        h2 = Compound(subcompounds=(hydrogen, h_clone))
        assert len(h2.all_ports()) == 2
        assert len(hydrogen.all_ports()) == 1
        assert len(h_clone.all_ports()) == 1

        mb.force_overlap(h_clone, h_clone["up"], hydrogen["up"])
        assert len(h2.all_ports()) == 0
        assert len(hydrogen.all_ports()) == 0
        assert len(h_clone.all_ports()) == 0

    def test_pruning_ghost_ports(self, ethane):
        eth1 = ethane
        h1 = eth1.particles_by_name("H")
        eth1.remove(h1)
        assert len(eth1.all_ports()) == 6
        for port in eth1.all_ports():
            assert port.anchor in list(eth1.particles())

        eth2 = mb.load("CC", smiles=True)
        h2 = eth2[-1]
        eth2.remove(h2)
        assert len(eth2.all_ports()) == len(eth2.available_ports()) == 1
        assert eth2.all_ports()[0] is eth2.available_ports()[0]

    def test_remove_bond_add_ports(self, hydrogen):
        h_clone = mb.clone(hydrogen)
        h2 = Compound(subcompounds=(hydrogen, h_clone))
        mb.force_overlap(h_clone, h_clone["up"], hydrogen["up"])
        h2.remove_bond((h2[0], h2[1]))
        assert len(h2.all_ports()) == 2
        assert len(hydrogen.all_ports()) == 1
        assert len(h_clone.all_ports()) == 1

    def test_reconnect_keeps_structure_x(self, chf, connect_and_reconnect):
        bond_vector = np.array([1, 0, 0])
        angle1, angle2 = connect_and_reconnect(chf, bond_vector)
        assert np.isclose(angle1, angle2, atol=1e-6)

    def test_reconnect_keeps_structure_y(self, chf, connect_and_reconnect):
        chf.spin(np.pi / 2, [1, 0, 0])
        bond_vector = np.array([0, 1, 0])
        angle1, angle2 = connect_and_reconnect(chf, bond_vector)
        assert np.isclose(angle1, angle2, atol=1e-6)

    def test_reconnect_keeps_structure_z(self, chf, connect_and_reconnect):
        bond_vector = np.array([0, 0, 1])
        angle1, angle2 = connect_and_reconnect(chf, bond_vector)
        assert np.isclose(angle1, angle2, atol=1e-6)

    def test_reconnect_keeps_structure_random(self, chf, connect_and_reconnect):
        np.random.seed(92)
        for _ in range(5):
            bond_vector = np.random.random(3) - 0.5
            angle1, angle2 = connect_and_reconnect(chf, bond_vector)
            assert np.isclose(angle1, angle2, atol=1e-6)

    def test_smarts_from_string(self):
        p3ht = mb.load("CCCCCCC1=C(SC(=C1)C)C", smiles=True, backend="pybel")
        assert p3ht.n_bonds == 33
        assert p3ht.n_particles == 33

    def test_smarts_from_file(self):
        p3ht = mb.load(get_fn("p3ht.smi"), smiles=True, backend="pybel")
        assert p3ht.n_bonds == 33
        assert p3ht.n_particles == 33

    @pytest.mark.skipif(not has_networkx, reason="NetworkX is not installed")
    def test_to_networkx(self):
        comp = Compound()
        comp.name = "Parent"

        for n in range(2):
            child = Compound()
            child.name = "c_{}".format(n)
            comp.add(child)
            for m in range(3):
                child_child = Compound()
                child_child.name = "c_{0}_{1}".format(m, n)
                child.add(child_child)

        graph = comp.to_networkx()

        assert graph.number_of_edges() == 8
        assert graph.number_of_nodes() == 9

        assert all([isinstance(n, Compound) for n in graph.nodes()])

    @pytest.mark.skipif(not has_networkx, reason="NetworkX is not installed")
    def test_to_networkx_no_hierarchy(self):
        comp = Compound()
        comp.name = "Parent"

        graph = comp.to_networkx()

        assert graph.number_of_edges() == 0
        assert graph.number_of_nodes() == 1

        assert all([isinstance(n, Compound) for n in graph.nodes()])

    @pytest.mark.skipif(not has_networkx, reason="NetworkX is not installed")
    def test_to_networkx_names_only(self):
        comp = Compound()
        comp.name = "Parent"

        for n in range(2):
            child = Compound()
            child.name = "c_{}".format(n)
            comp.add(child)
            for m in range(3):
                child_child = Compound()
                child_child.name = "c_{0}_{1}".format(m, n)
                child.add(child_child)

        graph = comp.to_networkx(names_only=True)

        assert graph.number_of_edges() == 8
        assert graph.number_of_nodes() == 9

        assert all([isinstance(n, str) for n in graph.nodes()])

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_from_trajectory(self):
        if has_mdtraj:
            mdtraj = import_("mdtraj")
        comp = Compound()
        traj = mdtraj.load(get_fn("spc.pdb"))
        comp.from_trajectory(traj)
        assert comp.children[0].name == "SPC"

    def test_from_parmed(self):
        comp = Compound()
        struc = pmd.load_file(get_fn("spc.pdb"))
        comp.from_parmed(struc)
        assert comp.children[0].name == "SPC"

    @pytest.mark.skipif(not has_mdtraj, reason="MDTraj not installed")
    def test_complex_from_trajectory(self):
        if has_mdtraj:
            mdtraj = import_("mdtraj")
        comp = Compound()
        traj = mdtraj.load(get_fn("pro_but.pdb"))
        comp.from_trajectory(traj)
        assert comp.children[0].children[0].name == "pro"
        assert comp.children[1].children[0].name == "but"

    def test_complex_from_parmed(self):
        comp = Compound()
        struc = pmd.load_file(get_fn("pro_but.pdb"))
        comp.from_parmed(struc)
        assert comp.children[0].name == "pro"
        assert comp.children[1].name == "but"

    @pytest.mark.skipif(not has_networkx, reason="NetworkX is not installed")
    def test_to_networkx_names_only_with_same_names(self):
        comp = Compound()
        comp.name = "compound"

        for n in range(2):
            child = Compound()
            child.name = "sub_compound"
            comp.add(child)
            for m in range(3):
                child_child = Compound()
                child_child.name = "sub_sub_compound"
                child.add(child_child)

        graph = comp.to_networkx(names_only=True)

        assert graph.number_of_edges() == 8
        assert graph.number_of_nodes() == 9

        assert all([isinstance(n, str) for n in graph.nodes()])

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_to_pybel(self, ethane):
        pybel_mol = ethane.to_pybel(box=None)
        assert pybel_mol.OBMol.NumAtoms() == 8
        assert pybel_mol.OBMol.NumBonds() == 7
        assert np.allclose(
            [
                pybel_mol.unitcell.GetA(),
                pybel_mol.unitcell.GetB(),
                pybel_mol.unitcell.GetC(),
            ],
            [2.139999, 2.9380001, 1.646],
        )

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_from_pybel(self):
        pybel = import_("pybel")
        benzene = list(pybel.readfile("mol2", get_fn("benzene.mol2")))[0]
        cmpd = Compound()
        cmpd.from_pybel(benzene)
        assert benzene.OBMol.NumAtoms() == cmpd.n_particles
        assert benzene.OBMol.NumBonds() == cmpd.n_bonds

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_to_pybel_residues(self, ethane):
        pybel_mol = ethane.to_pybel(box=None, residues="Ethane")
        assert "Ethane" in pybel_mol.residues[0].name

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_to_more_pybel_residues(self, methane, ethane):
        box = mb.fill_box(
            [methane, ethane],
            n_compounds=[3, 3],
            box=Box(lengths=[10, 10, 10], angles=[90.0, 90.0, 90.0]),
        )
        pybel_mol = box.to_pybel(box=None, residues=["Ethane", "Methane"])
        pybel_mol_resnames = {a.name for a in pybel_mol.residues}
        assert "Ethane" in pybel_mol_resnames
        assert "Methane" in pybel_mol_resnames

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_from_pybel_residues(self):
        pybel = import_("pybel")
        pybel_mol = list(pybel.readfile("mol2", get_fn("methyl.mol2")))[0]
        cmpd = Compound()
        cmpd.from_pybel(pybel_mol)
        assert "LIG1" in cmpd.children[0].name

    @pytest.mark.parametrize("extension", ["pdb", "sdf"])
    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_from_pybel_molecule(self, extension):
        pybel = import_("pybel")
        chol = list(
            pybel.readfile(
                extension, get_fn("cholesterol.{}".format(extension))
            )
        )[0]
        cmpd = mb.Compound()
        cmpd.from_pybel(chol)
        assert chol.OBMol.NumAtoms() == cmpd.n_particles
        assert chol.OBMol.NumBonds() == cmpd.n_bonds
        first_atom = chol.OBMol.GetAtom(1)
        assert np.allclose(
            cmpd[0].pos,
            [
                first_atom.GetX() / 10,
                first_atom.GetY() / 10,
                first_atom.GetZ() / 10,
            ],
        )

    @pytest.mark.parametrize(
        "test_smiles", ["CCO", "CCCCCCCC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
    )
    @pytest.mark.skipif(not has_rdkit, reason="RDKit is not installed")
    def test_from_rdkit_smiles(self, test_smiles):
        pos = list()
        for _ in range(3):
            cmpd = mb.load(test_smiles, smiles=True, backend="rdkit", seed=29)
            pos.append(cmpd.xyz)
        assert (
            np.diff(np.vstack(pos).reshape(len(pos), -1), axis=0) == 0
        ).all()

    @pytest.mark.parametrize("bad_smiles", ["F[P-](F)(F)(F)(F)F"])
    @pytest.mark.skipif(not has_rdkit, reason="RDKit is not installed")
    def test_incorrect_rdkit_smiles(self, bad_smiles):
        with pytest.raises(
            MBuildError,
            match=r"RDKit was unable to generate " r"3D coordinates",
        ):
            mb.load(bad_smiles, smiles=True, backend="rdkit", seed=29)

    @pytest.mark.skipif(not has_openbabel, reason="Pybel is not installed")
    def test_get_smiles(self):
        test_strings = ["CCO", "CCCCCCCC", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]
        for test_string in test_strings:
            my_cmp = mb.load(test_string, smiles=True, backend="pybel")
            assert my_cmp.get_smiles() == test_string

    def test_sdf(self, methane):
        methane.save("methane.sdf")
        sdf_string = mb.load("methane.sdf")
        assert np.allclose(methane.xyz, sdf_string.xyz, atol=1e-5)

    def test_load_multiple_sdf(self, methane):
        filled = mb.fill_box(methane, n_compounds=10, box=Box([4, 4, 4]))
        filled.save("methane.sdf")
        sdf_string = mb.load("methane.sdf")

    def test_save_multiple_sdf(self, methane):
        filled = mb.fill_box(methane, n_compounds=10, box=[0, 0, 0, 4, 4, 4])
        filled.save("methane.sdf")
        sdf_string = mb.load("methane.sdf")
        assert np.allclose(filled.xyz, sdf_string.xyz, atol=1e-5)

    def test_box(self):
        angles = [90.0, 90.0, 90.0]
        lengths = [3.0, 3.0, 3.0]
        compound = Compound()
        assert compound.box == None
        compound.box = Box([3.0, 3.0, 3.0])
        assert np.allclose(compound.box.lengths, lengths)
        assert np.allclose(compound.box.angles, angles)
        with pytest.raises(TypeError, match=r"specified as an mbuild.Box"):
            compound.box = "Hello, world"
        with pytest.raises(TypeError, match=r"specified as an mbuild.Box"):
            compound.box = [3.0, 3.0, 3.0]
        port = Port()
        assert port.box == None
        with pytest.raises(ValueError, match=r"cannot have"):
            port.box = Box(lengths=lengths, angles=angles)

        compound = Compound()
        subcomp = Compound(box=Box(lengths=lengths, angles=angles))
        compound.add(subcomp)
        assert np.allclose(compound.box.lengths, lengths)
        assert np.allclose(compound.box.angles, angles)
        compound = Compound(box=Box([3.0, 3.0, 3.0]))
        subcomp = Compound(box=Box([6.0, 6.0, 6.0], angles=[90.0, 90.0, 120.0]))
        with pytest.warns(UserWarning):
            compound.add(subcomp)
        assert np.allclose(compound.box.lengths, [3.0, 3.0, 3.0])
        assert np.allclose(compound.box.angles, [90.0, 90.0, 90.0])
        compound = Compound(box=Box(lengths=lengths, angles=angles))
        subcomp = Compound(
            box=Box(lengths=[6.0, 6.0, 6.0], angles=[90.0, 90.0, 120.0])
        )
        compound.add(subcomp, inherit_box=True)
        assert np.allclose(compound.box.lengths, [6.0, 6.0, 6.0])
        assert np.allclose(compound.box.angles, [90.0, 90.0, 120.0])
        compound = Compound(box=Box(lengths=lengths, angles=angles))
        subcomp = Compound()
        with pytest.warns(UserWarning):
            compound.add(subcomp, inherit_box=True)
        assert np.allclose(compound.box.lengths, [3.0, 3.0, 3.0])
        assert np.allclose(compound.box.angles, [90.0, 90.0, 90.0])
        compound = Compound()
        carbon = Compound(name="C")
        compound.add(carbon)
        compound.box = Box(lengths=lengths, angles=angles)
        nitrogen = Compound(name="N", pos=[4, 3, 3])
        with pytest.warns(UserWarning):
            compound.add(nitrogen)
        compound.box = Box(lengths=[5.0, 4.0, 4.0], angles=angles)
        with pytest.warns(UserWarning):
            compound.box = Box(lengths=[1.0, 1.0, 1.0], angles=angles)

    def test_get_boundingbox_extrema(self):
        h1 = mb.Compound()
        h2 = mb.Compound()
        h1.pos = [-0.07590747, 0.00182889, 0.00211742]
        h2.pos = [0.07590747, -0.00182889, -0.00211742]
        container = mb.Compound([h1, h2])
        distances = container.maxs - container.mins
        with pytest.raises(
            MBuildError, match=r"The vectors to define the box are co\-linear\,"
        ):
            container.get_boundingbox()
        distance_list = [val for val in distances]
        distance_list = [val + 1.0 for val in distance_list]
        np.testing.assert_almost_equal(
            container.get_boundingbox(pad_box=1.0).lengths,
            distance_list,
            decimal=6,
        )

        distance_list = [val for val in distances]
        distance_list[0] = distance_list[0] + 1.0
        distance_list[1] = distance_list[1] + 2.0
        distance_list[2] = distance_list[2] + 3.0
        np.testing.assert_almost_equal(
            container.get_boundingbox(pad_box=[1.0, 2.0, 3.0]).lengths,
            distance_list,
            decimal=6,
        )

    @pytest.mark.parametrize(
        "bad_value", [[1.0, 2.0], set([1, 2, 3]), {"x": 1.0}]
    )
    def test_get_boundingbox_error(self, bad_value):
        with pytest.raises(TypeError):
            meth = mb.load(get_fn("methyl.pdb"))
            meth.get_boundingbox(pad_box=bad_value)

    @pytest.mark.skipif(not has_py3Dmol, reason="Py3Dmol is not installed")
    def test_visualize_py3dmol(self, ethane):
        py3Dmol = import_("py3Dmol")
        vis_object = ethane._visualize_py3dmol()
        assert isinstance(vis_object, py3Dmol.view)

    @pytest.mark.skipif(not has_nglview, reason="NGLView is not installed")
    def test_visualize_nglview(self, ethane):
        nglview = import_("nglview")
        vis_object = ethane._visualize_nglview()
        assert isinstance(
            vis_object.component_0, nglview.component.ComponentViewer
        )

    def test_element(self):
        from ele import Elements

        na_compound = Compound(element="Na")
        assert na_compound.element == Elements.Na
        na_compound = Compound(element="NA")
        assert na_compound.element == Elements.Na
        na_compound = Compound(element="na")
        assert na_compound.element == Elements.Na
        co_compound = Compound(element="Co")
        assert co_compound.element != Elements.Na

        na_compound_clone = mb.clone(na_compound)
        assert na_compound_clone.element == Elements.Na
        container = Compound()
        container.add(na_compound)
        container.add(na_compound_clone)
        for child in container.children:
            assert child.element == Elements.Na

        na_compound = Compound()
        na_compound.element = "Na"
        assert na_compound.element == Elements.Na

    def test_invalid_element(self):
        from ele.exceptions import ElementError

        with pytest.raises(ElementError, match=r"No element with symbol"):
            na_compound = Compound(element="sodium")
        with pytest.raises(ElementError, match=r"No element with symbol"):
            na_compound = Compound(element="")

    def test_get_by_element(self):
        from ele import Elements
        from ele.exceptions import ElementError

        container = Compound()
        na = Compound(element="Na")
        na2 = Compound(element="Na")
        co = Compound(element="Co")
        container.add([na, na2, co])
        element_list = [c.element for c in container.particles_by_element("Na")]
        assert len(element_list) == 2
        for item in element_list:
            assert item == Elements.Na

        with pytest.raises(ElementError, match=r"No element with symbol"):
            element_list = [
                c.element for c in container.particles_by_element("sod")
            ]

    @pytest.mark.parametrize("backend", ["pybel", "rdkit"])
    def test_elements_from_smiles(self, backend):
        mol = mb.load("COC", smiles=True, backend=backend)
        for particle in mol.particles():
            assert particle.element is not None

    def test_mins_maxs(self, benzene):
        assert np.allclose(benzene.mins, [-0.2267, -0.15422, 0.0])
        assert np.allclose(benzene.maxs, [0.20318, 0.34207, 0.0])

    def test_periodicity_raises(self):
        with pytest.raises(ValueError):
            Compound(periodicity=[True, True])
        with pytest.raises(TypeError):
            Compound(periodicity=[1, 2, 3])

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_xyz_setter(self):
        # Test taken from issue 829
        import foyer

        # Build a CG hexane with mBuild
        ethane = mb.Compound()
        C1 = mb.Compound(name="_CH3")
        C2 = mb.Compound(name="_CH3")

        ethane.add([C1, C2])
        ethane.add_bond([C1, C2])

        # This causes a failure on the ff.apply line
        ethane.xyz = [[0.1, 0, 0], [0.2, 0, 0]]
        ff = foyer.forcefields.load_TRAPPE_UA()
        # This fails prior to applying PR # 892
        ff.apply(ethane)

    def test_ordered_bonds(self):
        ethane = mb.load("CC", smiles=True)
        ethane2 = mb.load("CC", smiles=True)
        for bond2, bond in zip(ethane2.bonds(), ethane.bonds()):
            assert bond2[0].name == bond[0].name
            assert all(bond2[0].pos == bond[0].pos)
