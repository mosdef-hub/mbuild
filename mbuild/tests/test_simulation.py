import hoomd
import numpy as np
import pytest
from gmso import ForceField

import mbuild as mb
from mbuild.simulation import (
    ForcesHandler,
    HoomdSimulation,
    OpenMMSimulation,
    PropertyLogger,
)
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer, has_hoomd


def _make_simple_compound(n=5, spacing=0.15):
    """Create a simple linear compound with bonds and proper elements."""

    cpd = mb.Compound()
    for i in range(n):
        p = mb.Compound(name="C", pos=[i * spacing, 0.0, 0.0], element="C")
        cpd.add(p)
    particles = list(cpd.particles())
    for i in range(len(particles) - 1):
        cpd.add_bond((particles[i], particles[i + 1]))
    return cpd


def _make_simple_path(n=8, spacing=0.15):
    """Create a simple Path object with name attributes on nodes."""
    import networkx as nx

    from mbuild.path.build import Path

    coords = np.array([[i * spacing, 0.0, 0.0] for i in range(n)])
    bond_graph = nx.Graph()
    for i in range(n):
        bond_graph.add_node(i, name="_A")
    for i in range(n - 1):
        bond_graph.add_edge(i, i + 1)
    return Path(coordinates=coords, bond_graph=bond_graph)


class TestForcesHandler(BaseTest):
    """Tests for the ForcesHandler class."""

    @pytest.fixture
    def oplsaa(self):
        return ForceField("oplsaa")

    @pytest.fixture
    def sim(self, octane, oplsaa):
        return HoomdSimulation(octane, oplsaa, r_cut=0.3, run_on_gpu=False)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_default_init(self):
        """Test default ForcesHandler initialization."""
        fh = ForcesHandler()
        assert fh.dpd == 0.0
        assert fh.scale_forces["bond"] == 1
        assert fh.scale_forces["angle"] == 1
        assert fh.scale_forces["lj"] == 1
        assert fh.scale_forces["charge"] == 0
        assert fh.scale_forces["opls"] == 0
        assert fh.scale_forces["periodic"] == 0
        assert fh.scale_forces["improper"] == 0

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_custom_init(self):
        """Test ForcesHandler with custom scaling."""
        fh = ForcesHandler(dpd=5.0, scale_bond=0.5, scale_angle=0.2, scale_lj=0.8)
        assert fh.dpd == 5.0
        assert fh.scale_forces["bond"] == 0.5
        assert fh.scale_forces["angle"] == 0.2
        assert fh.scale_forces["lj"] == 0.8

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_sim_dpd(self, sim):
        """Test that DPD replaces LJ when dpd > 0."""
        fh = ForcesHandler(dpd=1.0, scale_bond=0, scale_angle=0)
        fh.scale_sim(sim)
        lj_force = sim._get_force(hoomd.md.pair.LJ)
        dpd_force = sim._get_force(hoomd.md.pair.DPDConservative, active_only=True)
        for param in dpd_force.params:
            assert dpd_force.params[param]["A"] == 1.0
            assert dpd_force.r_cut[param] == lj_force.params[param]["sigma"]

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_sim_excludes_zero_forces(self, sim):
        """Test that forces with scale=0 are excluded from active_forces."""
        fh = ForcesHandler(scale_lj=1, scale_bond=0, scale_angle=0)
        fh.scale_sim(sim)
        force_types = [type(f) for f in sim.active_forces]
        assert hoomd.md.bond.Harmonic not in force_types
        assert hoomd.md.angle.Harmonic not in force_types

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_sim_includes_nonzero_forces(self, sim):
        """Test that forces with scale>0 are included in active_forces."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        fh.scale_sim(sim)
        force_types = [type(f) for f in sim.active_forces]
        assert hoomd.md.pair.LJ in force_types
        assert hoomd.md.bond.Harmonic in force_types
        assert hoomd.md.angle.Harmonic in force_types

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_bond_values(self, sim):
        """Test that bond k values are properly scaled."""
        bond = sim._get_force(hoomd.md.bond.Harmonic)
        orig_params = {p: dict(bond.params[p]) for p in bond.params}

        ForcesHandler(scale_bond=0.5, scale_angle=0).scale_sim(sim)
        for p in bond.params:
            assert np.isclose(bond.params[p]["k"], orig_params[p]["k"] * 0.5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_angle_values(self, sim):
        """Test that angle k values are properly scaled."""
        angle = sim._get_force(hoomd.md.angle.Harmonic)
        orig_params = {p: dict(angle.params[p]) for p in angle.params}

        ForcesHandler(scale_bond=0, scale_angle=0.3).scale_sim(sim)
        for p in angle.params:
            assert np.isclose(angle.params[p]["k"], orig_params[p]["k"] * 0.3)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_opls(self, sim):
        """Test that OPLS dihedral forces can be activated."""
        fh = ForcesHandler(scale_opls=0.5)
        fh.scale_sim(sim)
        force_types = [type(f) for f in sim.active_forces]
        assert hoomd.md.dihedral.OPLS in force_types

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_restores_on_second_call(self, sim):
        """Test that calling scale_sim again restores original params before scaling."""
        bond = sim._get_force(hoomd.md.bond.Harmonic)
        angle = sim._get_force(hoomd.md.angle.Harmonic)
        orig_bond_params = {p: dict(bond.params[p]) for p in bond.params}
        orig_angle_params = {p: dict(angle.params[p]) for p in angle.params}

        ForcesHandler(scale_bond=0.5, scale_angle=0.5).scale_sim(sim)
        for p in bond.params:
            assert np.isclose(bond.params[p]["k"], orig_bond_params[p]["k"] * 0.5)
        for p in angle.params:
            assert np.isclose(angle.params[p]["k"], orig_angle_params[p]["k"] * 0.5)

        ForcesHandler(scale_bond=1, scale_angle=1).scale_sim(sim)
        for p in bond.params:
            assert np.isclose(bond.params[p]["k"], orig_bond_params[p]["k"])
        for p in angle.params:
            assert np.isclose(angle.params[p]["k"], orig_angle_params[p]["k"])

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_forces_dict_populated(self, sim):
        """Test that forcesDict is populated after scale_sim."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        fh.scale_sim(sim)
        assert "lj" in fh.forcesDict
        assert "bond" in fh.forcesDict
        assert "angle" in fh.forcesDict

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_missing_force_skipped(self):
        """Test that missing forces (e.g., no dihedrals) are skipped gracefully."""
        cpd = mb.Compound()
        for i in range(3):
            p = mb.Compound(name="C", pos=[i * 0.15, 0.0, 0.0], element="C")
            cpd.add(p)
        particles = list(cpd.particles())
        for i in range(len(particles) - 1):
            cpd.add_bond((particles[i], particles[i + 1]))

        sim = HoomdSimulation(cpd, forcefield=None, r_cut=0.5, run_on_gpu=False)
        # No OPLS dihedrals exist in the default system
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1, scale_opls=1)
        fh.scale_sim(sim)  # Should not raise


class TestHoomdSimulation(BaseTest):
    """Tests for the HoomdSimulation class."""

    @pytest.fixture
    def oplsaa(self):
        return ForceField("oplsaa")

    @pytest.fixture
    def sim(self, octane, oplsaa):
        return HoomdSimulation(octane, oplsaa, r_cut=0.3, run_on_gpu=False)

    @pytest.fixture
    def sim_no_ff(self, octane):
        return HoomdSimulation(octane, forcefield=None, r_cut=0.3, run_on_gpu=False)

    # ── Initialization tests ──────────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_with_compound(self, octane, oplsaa):
        """Test initialization with an mb.Compound."""
        sim = HoomdSimulation(octane, oplsaa, r_cut=0.3, run_on_gpu=False)
        assert sim._is_path is False
        assert sim.compound is octane
        assert len(sim.forces) > 0

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_with_path(self):
        """Test initialization with a Path object."""
        path = _make_simple_path(n=6)
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        assert sim._is_path is True
        assert sim.compound.n_particles == 6

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_invalid_input(self):
        """Test that invalid input raises TypeError."""
        with pytest.raises(TypeError):
            HoomdSimulation("invalid", forcefield=None, r_cut=0.3)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_no_forcefield_creates_defaults(self, sim_no_ff):
        """Test that default forces are created when no FF is given."""
        force_types = [type(f) for f in sim_no_ff.forces]
        assert hoomd.md.pair.LJ in force_types
        assert hoomd.md.bond.Harmonic in force_types

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_with_forcefield(self, sim):
        """Test that forces are created from the provided forcefield."""
        force_types = [type(f) for f in sim.forces]
        assert hoomd.md.pair.LJ in force_types
        assert hoomd.md.bond.Harmonic in force_types
        assert hoomd.md.angle.Harmonic in force_types

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_init_with_logger(self, octane, oplsaa):
        """Test initialization with a PropertyLogger."""
        logger = PropertyLogger(properties=["potential_energy"])
        sim = HoomdSimulation(
            octane, oplsaa, r_cut=0.3, run_on_gpu=False, logger=logger
        )
        assert sim.logger is logger

    # ── get_force tests ───────────────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_get_force_exists(self, sim):
        """Test retrieving existing forces."""
        force = sim._get_force(hoomd.md.pair.LJ)
        assert force is not None
        assert isinstance(force, hoomd.md.pair.LJ)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_get_force_not_found(self, sim_no_ff):
        """Test that ValueError is raised for non-existent force type."""
        with pytest.raises(ValueError):
            sim_no_ff._get_force(hoomd.md.dihedral.OPLS)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_get_force_active_only(self, sim):
        """Test get_force with active_only flag."""
        fh = ForcesHandler(scale_lj=1, scale_bond=0, scale_angle=0)
        fh.scale_sim(sim)
        force = sim._get_force(hoomd.md.pair.LJ, active_only=True)
        assert force is not None
        with pytest.raises(ValueError):
            sim._get_force(hoomd.md.bond.Harmonic, active_only=True)

    # ── Simulation method tests ───────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_nvt(self, sim):
        """Test NVT simulation method."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        old_coords = sim.compound.xyz.copy()
        sim.nvt(forces_handler=fh, n_steps=20, kT=2.0, dt=1e-4, tau=1e-2)
        assert not np.allclose(old_coords, sim.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_nvt_no_forces_handler(self, sim_no_ff):
        """Test NVT with forces_handler=None uses all forces."""
        old_coords = sim_no_ff.compound.xyz.copy()
        sim_no_ff.nvt(forces_handler=None, n_steps=20, kT=2.0, dt=1e-4, tau=1e-2)
        assert not np.allclose(old_coords, sim_no_ff.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fire(self, sim):
        """Test FIRE minimization method."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        old_coords = sim.compound.xyz.copy()
        sim.fire(forces_handler=fh, n_steps=20)
        assert not np.allclose(old_coords, sim.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fire_no_forces_handler(self, sim_no_ff):
        """Test FIRE with no forces handler."""
        old_coords = sim_no_ff.compound.xyz.copy()
        sim_no_ff.fire(forces_handler=None, n_steps=20)
        assert not np.allclose(old_coords, sim_no_ff.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fire_multiple_iterations(self, sim):
        """Test FIRE with multiple iterations."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        old_coords = sim.compound.xyz.copy()
        sim.fire(forces_handler=fh, n_steps=200, n_iterations=3)
        new_coords = sim.compound.xyz
        assert not np.allclose(old_coords, new_coords, atol=1e-4)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_cap_displacement(self, sim):
        """Test capped displacement method."""
        fh = ForcesHandler(dpd=1, scale_bond=0.1, scale_angle=0)
        old_coords = sim.compound.xyz.copy()
        sim.cap_displacement(forces_handler=fh, dt=1, max_displacement=1, n_steps=10)
        assert not np.allclose(old_coords, sim.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_cap_displacement_no_forces_handler(self, sim_no_ff):
        """Test capped displacement with forces_handler=None."""
        old_coords = sim_no_ff.compound.xyz.copy()
        sim_no_ff.cap_displacement(
            forces_handler=None, dt=1, max_displacement=1, n_steps=10
        )
        assert not np.allclose(old_coords, sim_no_ff.compound.xyz, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_cap_displacement_force_types(self, sim):
        """Test that only expected forces are active after cap_displacement."""
        fh = ForcesHandler(dpd=1, scale_bond=0.1, scale_angle=0)
        sim.cap_displacement(forces_handler=fh, dt=1, max_displacement=1, n_steps=10)
        force_types = set(type(f) for f in sim.active_forces)
        assert force_types == {hoomd.md.pair.DPDConservative, hoomd.md.bond.Harmonic}

    # ── Energy tracking tests ─────────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_get_energy_empty(self, sim):
        """Test get_energy returns None before any run."""
        assert sim.get_energy() is None

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_get_energy_after_run(self, sim):
        """Test get_energy returns data after running."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        sim.cap_displacement(forces_handler=fh, dt=1, max_displacement=1, n_steps=10)
        energy = sim.get_energy()
        assert energy is not None
        assert len(energy) > 0
        for key in energy:
            assert len(energy[key]) == 2

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_energies_accumulate(self, sim):
        """Test that energies accumulate over multiple calls."""
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        sim.cap_displacement(forces_handler=fh, dt=1, max_displacement=1, n_steps=5)
        sim.cap_displacement(forces_handler=fh, dt=1, max_displacement=1, n_steps=5)
        energy = sim.get_energy()
        for key in energy:
            assert len(energy[key]) == 3

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_property_logger(self, octane, oplsaa):
        """Test that PropertyLogger records data during runs."""
        logger = PropertyLogger(properties=["potential_energy"])
        sim = HoomdSimulation(
            octane, oplsaa, r_cut=0.3, run_on_gpu=False, logger=logger
        )
        fh = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        sim.nvt(forces_handler=fh, n_steps=20, kT=2.0, dt=1e-4, tau=1e-2)
        assert len(logger.data["potential_energy"]) > 0

    # ── Path-specific tests ───────────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_fire(self):
        """Test FIRE minimization on a Path."""
        path = _make_simple_path(n=10)
        old_coords = path.coordinates.copy()
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        sim.fire(forces_handler=None, n_steps=50)
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_nvt(self):
        """Test NVT simulation on a Path."""
        path = _make_simple_path(n=10)
        old_coords = path.coordinates.copy()
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        sim.nvt(forces_handler=None, n_steps=20, kT=1.0, dt=1e-4, tau=1e-2)
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_syncs_back(self):
        """Test that Path coordinates are synced back after simulation."""
        path = _make_simple_path(n=5)
        old_coords = path.coordinates.copy()
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        sim.fire(forces_handler=None, n_steps=50)
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)
        assert np.allclose(sim.compound.xyz, path.coordinates)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_cap_displacement(self):
        """Test capped displacement on a Path."""
        path = _make_simple_path(n=8)
        old_coords = path.coordinates.copy()
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        sim.cap_displacement(
            forces_handler=None, dt=1, max_displacement=0.5, n_steps=10
        )
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)

    # ── Integration group tests ───────────────────────────────────────────

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_integrate_group_all(self, octane, oplsaa):
        """Test that default integration group includes all particles."""
        sim = HoomdSimulation(octane, oplsaa, r_cut=0.3, run_on_gpu=False)
        group = sim._get_integrate_group()
        assert isinstance(group, hoomd.filter.All)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_integrate_group_both_raises(self, octane, oplsaa):
        """Test that passing both integrate and fixed compounds raises."""
        children = list(octane.children)
        if len(children) >= 2:
            sim = HoomdSimulation(
                octane,
                oplsaa,
                r_cut=0.3,
                run_on_gpu=False,
                integrate_compounds=[children[0]],
                fixed_compounds=[children[1]],
            )
            with pytest.raises(RuntimeError):
                sim._get_integrate_group()


class TestOpenMMSimulation(BaseTest):
    """Tests for the OpenMMSimulation class."""

    @pytest.fixture
    def simple_compound(self):
        return _make_simple_compound(n=5, spacing=0.15)

    @pytest.fixture
    def simple_path(self):
        return _make_simple_path(n=8)

    @pytest.fixture
    def compound_with_box(self):
        """Compound with box, particles centered inside."""
        cpd = mb.Compound()
        for i in range(5):
            p = mb.Compound(
                name="C",
                pos=[4.0 + i * 0.5, 5.0, 5.0],  # centered in box
                element="C",
            )
            cpd.add(p)
        particles = list(cpd.particles())
        for i in range(len(particles) - 1):
            cpd.add_bond((particles[i], particles[i + 1]))
        cpd.box = mb.Box([10, 10, 10])
        return cpd

    @pytest.fixture
    def compound_no_bonds(self):
        """Compound with no bonds (isolated particles)."""

        cpd = mb.Compound()
        for i in range(3):
            p = mb.Compound(name="C", pos=[i * 0.4, 0.0, 0.0], element="C")
            cpd.add(p)
        return cpd

    # ── Initialization tests ──────────────────────────────────────────────

    def test_init_with_compound(self, simple_compound):
        """Test initialization with a Compound."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        assert sim._is_path is False
        assert sim.compound is simple_compound
        assert sim.system is not None
        assert sim.topology is not None

    def test_init_with_path(self, simple_path):
        """Test initialization with a Path."""
        sim = OpenMMSimulation(simple_path, forcefield=None)
        assert sim._is_path is True
        assert sim.compound.n_particles == 8

    def test_init_invalid_input(self):
        """Test that invalid input raises TypeError."""
        with pytest.raises(TypeError):
            OpenMMSimulation("not_valid", forcefield=None)

    def test_init_invalid_constraints(self, simple_compound):
        """Test that invalid constraints raise ValueError."""
        with pytest.raises(ValueError):
            OpenMMSimulation(simple_compound, forcefield=None, constraints="BadValue")

    @pytest.mark.parametrize("constraints", [None, "AllBonds", "HBonds", "HAngles"])
    def test_init_valid_constraints(self, simple_compound, constraints):
        """Test that all valid constraint options are accepted."""
        sim = OpenMMSimulation(
            simple_compound, forcefield=None, constraints=constraints
        )
        assert sim is not None

    def test_init_with_logger(self, simple_compound):
        """Test initialization with a PropertyLogger."""
        logger = PropertyLogger(properties=["potential_energy"])
        sim = OpenMMSimulation(simple_compound, forcefield=None, logger=logger)
        assert sim.logger is logger

    def test_init_default_forces_with_bonds(self, simple_compound):
        """Test that default system has LJ + bond forces."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        n_forces = sim.system.getNumForces()
        assert n_forces >= 2  # CustomNonbondedForce + HarmonicBondForce

    def test_init_default_forces_no_bonds(self, compound_no_bonds):
        """Test default system with only LJ (no bonds)."""
        sim = OpenMMSimulation(compound_no_bonds, forcefield=None)
        assert sim.system.getNumForces() >= 1

    def test_init_with_box_uses_periodic(self, compound_with_box):
        """Test that a compound with a box uses periodic boundary conditions."""
        import openmm

        sim = OpenMMSimulation(compound_with_box, forcefield=None)
        for i in range(sim.system.getNumForces()):
            force = sim.system.getForce(i)
            if isinstance(force, openmm.CustomNonbondedForce):
                assert (
                    force.getNonbondedMethod()
                    == openmm.CustomNonbondedForce.CutoffPeriodic
                )
                break

    def test_init_without_box_uses_nocutoff(self, simple_compound):
        """Test that a compound without a box uses NoCutoff."""
        import openmm

        simple_compound.box = None
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        for i in range(sim.system.getNumForces()):
            force = sim.system.getForce(i)
            if isinstance(force, openmm.CustomNonbondedForce):
                assert (
                    force.getNonbondedMethod()
                    == openmm.CustomNonbondedForce.NoCutoff
                )
                break

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_init_with_forcefield(self, octane):
        """Test initialization with a foyer forcefield."""
        sim = OpenMMSimulation(octane, forcefield="oplsaa")
        assert sim.system is not None

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_init_with_xml_forcefield(self, octane):
        """Test initialization with an XML forcefield file."""
        sim = OpenMMSimulation(octane, forcefield=get_fn("small_oplsaa.xml"))
        assert sim.system is not None

    def test_init_platform_cpu(self, simple_compound):
        """Test that CPU platform is set correctly."""
        sim = OpenMMSimulation(simple_compound, forcefield=None, platform="CPU")
        assert sim._platform.getName() == "CPU"

    def test_init_nthreads(self, simple_compound):
        """Test that nthreads is set."""
        sim = OpenMMSimulation(simple_compound, forcefield=None, nthreads=2)
        assert sim._platform.getPropertyDefaultValue("Threads") == "2"

    # ── minimize tests ────────────────────────────────────────────────────

    def test_minimize_compound(self, simple_compound):
        """Test energy minimization on a Compound."""
        old_coords = simple_compound.xyz.copy()
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=100)
        assert not np.allclose(old_coords, simple_compound.xyz, atol=1e-6)

    def test_minimize_path(self, simple_path):
        """Test energy minimization on a Path."""
        old_coords = simple_path.coordinates.copy()
        sim = OpenMMSimulation(simple_path, forcefield=None)
        sim.minimize(n_steps=100)
        assert not np.allclose(old_coords, simple_path.coordinates, atol=1e-6)

    def test_minimize_records_energy(self, simple_compound):
        """Test that minimize records energy."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        energy = sim.get_energy()
        assert energy is not None
        assert "potential_energy" in energy
        assert len(energy["potential_energy"]) == 1

    def test_minimize_custom_tolerance(self, simple_compound):
        """Test minimize with custom tolerance."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=100, tolerance=1.0)
        assert sim.get_energy() is not None

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_minimize_with_forcefield(self, octane):
        """Test minimize with a real forcefield."""
        # Distort the molecule away from its minimum
        for i, p in enumerate(octane.particles()):
            p.pos = p.pos + np.array([0.01 * (i % 3), -0.01 * (i % 2), 0.005 * i])
        old_coords = octane.xyz.copy()
        sim = OpenMMSimulation(octane, forcefield="oplsaa")
        sim.minimize(n_steps=100)
        assert not np.allclose(old_coords, octane.xyz, atol=1e-4)

    # ── nvt tests ─────────────────────────────────────────────────────────

    def test_nvt_compound(self, simple_compound):
        """Test NVT simulation on a Compound."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=1000)
        old_coords = simple_compound.xyz.copy()
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        assert not np.allclose(old_coords, simple_compound.xyz, atol=1e-6)

    def test_nvt_path(self, simple_path):
        """Test NVT simulation on a Path."""
        sim = OpenMMSimulation(simple_path, forcefield=None)
        sim.minimize(n_steps=100)
        old_coords = simple_path.coordinates.copy()
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        assert not np.allclose(old_coords, simple_path.coordinates, atol=1e-6)

    def test_nvt_records_energy(self, simple_compound):
        """Test that NVT records energy."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        energy = sim.get_energy()
        assert energy is not None
        # minimize + nvt = 2 entries
        assert len(energy["potential_energy"]) == 2

    def test_nvt_with_report_interval(self, simple_compound):
        """Test NVT with periodic reporting."""
        logger = PropertyLogger(properties=["potential_energy", "kinetic_energy"])
        sim = OpenMMSimulation(simple_compound, forcefield=None, logger=logger)
        sim.nvt(n_steps=100, kT=300, dt=0.0001, report_interval=25)
        assert len(logger.data["potential_energy"]) == 4
        assert len(logger.data["kinetic_energy"]) == 4

    def test_nvt_custom_friction(self, simple_compound):
        """Test NVT with custom friction coefficient."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        sim.nvt(n_steps=20, kT=300, dt=0.0001, friction=5.0)
        assert sim.get_energy() is not None

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_nvt_with_forcefield(self, octane):
        """Test NVT with a real forcefield."""
        octane.box = mb.Box([5, 5, 5])
        sim = OpenMMSimulation(octane, forcefield="oplsaa")
        sim.minimize(n_steps=200)
        old_coords = octane.xyz.copy()
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        assert not np.allclose(old_coords, octane.xyz, atol=1e-6)

    # ── get_energy tests ──────────────────────────────────────────────────

    def test_get_energy_empty(self, simple_compound):
        """Test get_energy returns None before any run."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        assert sim.get_energy() is None

    def test_get_energy_accumulates(self, simple_compound):
        """Test that energies accumulate across multiple calls."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        sim.nvt(n_steps=20, kT=300, dt=0.0001)
        energy = sim.get_energy()
        assert len(energy["potential_energy"]) == 2

    # ── PropertyLogger integration ────────────────────────────────────────

    def test_logger_records_pe(self, simple_compound):
        """Test logger records potential energy."""
        logger = PropertyLogger(properties=["potential_energy"])
        sim = OpenMMSimulation(simple_compound, forcefield=None, logger=logger)
        sim.minimize(n_steps=50)
        assert len(logger.data["potential_energy"]) == 1
        assert isinstance(logger.data["potential_energy"][0], float)

    def test_logger_records_ke(self, simple_compound):
        """Test logger records kinetic energy during dynamics."""
        logger = PropertyLogger(properties=["kinetic_energy"])
        sim = OpenMMSimulation(simple_compound, forcefield=None, logger=logger)
        sim.nvt(n_steps=50, kT=300, dt=0.0001, report_interval=25)
        assert len(logger.data["kinetic_energy"]) == 2

    def test_logger_as_dict(self, simple_compound):
        """Test logger.as_dict() returns numpy arrays."""
        logger = PropertyLogger(properties=["potential_energy"])
        sim = OpenMMSimulation(simple_compound, forcefield=None, logger=logger)
        sim.minimize(n_steps=50)
        data = logger.as_dict()
        assert isinstance(data["potential_energy"], np.ndarray)

    def test_logger_reset(self):
        """Test logger reset clears all data."""
        logger = PropertyLogger(properties=["potential_energy", "kinetic_energy"])
        logger.data["potential_energy"].append(1.0)
        logger.data["kinetic_energy"].append(2.0)
        logger.reset()
        assert len(logger.data["potential_energy"]) == 0
        assert len(logger.data["kinetic_energy"]) == 0

    # ── Path sync-back tests ──────────────────────────────────────────────

    def test_path_syncs_back_minimize(self, simple_path):
        """Test Path coordinates are updated after minimize."""
        old_coords = simple_path.coordinates.copy()
        sim = OpenMMSimulation(simple_path, forcefield=None)
        sim.minimize(n_steps=100)
        assert not np.allclose(old_coords, simple_path.coordinates, atol=1e-6)
        assert np.allclose(sim.compound.xyz, simple_path.coordinates)

    def test_path_syncs_back_nvt(self, simple_path):
        """Test Path coordinates are updated after NVT."""
        sim = OpenMMSimulation(simple_path, forcefield=None)
        sim.minimize(n_steps=100)
        old_coords = simple_path.coordinates.copy()
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        assert not np.allclose(old_coords, simple_path.coordinates, atol=1e-6)
        assert np.allclose(sim.compound.xyz, simple_path.coordinates)

    # ── Multiple calls / reuse tests ──────────────────────────────────────

    def test_minimize_then_nvt(self, simple_compound):
        """Test calling minimize then nvt sequentially."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        coords_after_min = simple_compound.xyz.copy()
        sim.nvt(n_steps=50, kT=300, dt=0.0001)
        assert not np.allclose(coords_after_min, simple_compound.xyz, atol=1e-6)
        energy = sim.get_energy()
        assert len(energy["potential_energy"]) == 2
