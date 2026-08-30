import hoomd
import numpy as np
import pytest
from gmso import ForceField

import mbuild as mb
from mbuild.simulation import (
    ForcesHandler,
    HoomdSimulation,
    OpenMMSimulation,
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


def _make_two_type_path(n=6, spacing=0.15):
    """Create a straight-line Path with alternating A/B bead types."""
    from mbuild.path.build import straight_line
    from mbuild.path.namers import CyclicNamer

    return straight_line(spacing=spacing, N=n, bead_name=CyclicNamer(["A", "B"]))


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
    def test_get_energy_nan_for_inactive_forces(self, sim):
        """Forces inactive in a frame read NaN, not zero."""
        fh = ForcesHandler(dpd=1, scale_bond=0.1, scale_angle=0)
        sim.cap_displacement(
            forces_handler=fh, dt=1e-3, max_displacement=1e-3, n_steps=5
        )
        sim.fire(n_steps=5)
        energy = sim.get_energy()

        def key(cls):
            return cls.__module__ + "." + cls.__name__

        # DPD replaces LJ only while the handler is applied.
        assert np.isnan(energy[key(hoomd.md.pair.DPDConservative)][-1])
        assert np.all(np.isnan(energy[key(hoomd.md.pair.LJ)][:-1]))
        assert not np.any(np.isnan(energy["potential_energy"]))

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
            forces_handler=None, dt=1, max_displacement=0.001, n_steps=10
        )
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_explicit_box(self):
        """An explicit box=[Lx, Ly, Lz] is used verbatim and stays fixed."""
        path = _make_simple_path(n=6)
        box = [5.0, 6.0, 7.0]
        for _ in range(2):
            sim = HoomdSimulation(path, forcefield=None, box=box, run_on_gpu=False)
            assert np.allclose(sim.state.box.L, box)
            sim.fire(forces_handler=None, n_steps=10)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_per_type_radius(self):
        """Per-bead-type radii set unlike-pair WCA sigma to the arithmetic mean."""
        from mbuild.utils.simulation.path_forces import PathForcefield

        path = _make_two_type_path(n=6)
        pff = PathForcefield(radius={"A": 0.4, "B": 0.8})
        sim = HoomdSimulation(path, forcefield=pff, r_cut=0.5, run_on_gpu=False)
        wca = sim.forces[0]
        assert wca.params[("A", "A")]["sigma"] == pytest.approx(0.4)
        assert wca.params[("B", "B")]["sigma"] == pytest.approx(0.8)
        assert wca.params[("A", "B")]["sigma"] == pytest.approx(0.6)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_per_type_bond_length(self):
        """Bond length keyed by bond-type name sizes the FENE bond."""
        from mbuild.utils.simulation.path_forces import PathForcefield

        path = _make_two_type_path(n=6)
        pff = PathForcefield(radius=0.4, bond_length={"A-B": 0.25})
        sim = HoomdSimulation(path, forcefield=pff, r_cut=0.5, run_on_gpu=False)
        fene = sim.forces[1]
        assert fene.params["A-B"]["sigma"] == pytest.approx(0.25)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_auto_per_type(self):
        """forcefield=None derives per-type radii/lengths and relaxes a 2-type path."""
        path = _make_two_type_path(n=8)
        old_coords = path.coordinates.copy()
        sim = HoomdSimulation(path, forcefield=None, r_cut=0.5, run_on_gpu=False)
        wca_pairs = set(sim.forces[0].params.keys())
        assert {("A", "A"), ("A", "B"), ("B", "B")} <= wca_pairs
        sim.fire(forces_handler=None, n_steps=50)
        assert not np.allclose(old_coords, path.coordinates, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_missing_type_raises(self):
        """A per-type radius dict missing a bead type raises a clear error."""
        from mbuild.utils.simulation.path_forces import PathForcefield

        path = _make_two_type_path(n=6)
        with pytest.raises(KeyError, match="type"):
            HoomdSimulation(
                path, forcefield=PathForcefield(radius={"A": 0.4}), run_on_gpu=False
            )

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_path_per_type_angle_partial(self):
        """A per-type angle dict covers only its types and still runs."""
        from mbuild.path.points import AnglesSampler
        from mbuild.utils.simulation.path_forces import PathForcefield

        path = _make_two_type_path(n=6)
        sampler = AnglesSampler("normal", {"loc": 2.5, "scale": 0.2})
        pff = PathForcefield(radius=0.15, angles={"A-B-A": sampler})
        sim = HoomdSimulation(path, forcefield=pff, r_cut=0.5, run_on_gpu=False)
        angle = sim.forces[-1]
        assert isinstance(angle, hoomd.md.angle.Table)
        # Covered type carries a nonzero potential; uncovered type is flat.
        assert np.any(np.asarray(angle.params["A-B-A"]["U"]) != 0)
        assert np.all(np.asarray(angle.params["B-A-B"]["U"]) == 0)
        sim.cap_displacement(dt=1, max_displacement=0.003, n_steps=50)
        sim.fire(forces_handler=None, n_steps=50)

    def test_auto_bond_length_median_ignores_wrapping(self):
        """Median sampling ignores periodic-wrap bonds seen as long outliers."""
        import gsd.hoomd

        from mbuild.utils.simulation.path_forces import (
            _auto_bond_lengths,
            _auto_radii,
        )

        frame = gsd.hoomd.Frame()
        frame.particles.N = 5
        # Four in a row 0.25 apart, plus one placed far away.
        frame.particles.position = [
            [0, 0, 0],
            [0.25, 0, 0],
            [0.5, 0, 0],
            [0.75, 0, 0],
            [4.5, 0, 0],
        ]
        frame.particles.types = ["A"]
        frame.particles.typeid = [0] * 5
        frame.configuration.box = [10, 10, 10, 0, 0, 0]
        # Three real 0.25 bonds and one long outlier (0,4), like a wrap artifact.
        frame.bonds.N = 4
        frame.bonds.types = ["A-A"]
        frame.bonds.typeid = [0, 0, 0, 0]
        frame.bonds.group = [[0, 1], [1, 2], [2, 3], [0, 4]]
        assert _auto_bond_lengths(frame)["A-A"] == pytest.approx(0.25)
        assert _auto_radii(frame)["A"] == pytest.approx(0.25 / 0.97)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_stress_uff(self):
        """UFF on a 60-atom molecule spanning all 10 elements and every hybridization."""
        smiles = "N#CC(=O)OC(O)C(Br)(Cl)c1ccc(o1)C=Nc1ncc(F)cc1C#CCSCP(C)CCN(C)CI"
        comp = mb.load(smiles, smiles=True)
        old_coords = comp.xyz.copy()
        sim = HoomdSimulation(
            comp, forcefield=None, r_cut=0.5, box_buffer=6, run_on_gpu=False
        )
        # LJ + harmonic bonds + harmonic angles + periodic dihedrals.
        assert len(sim.forces) == 4
        sim.fire(forces_handler=None, n_steps=100)
        assert np.isfinite(comp.xyz).all()
        assert not np.allclose(old_coords, comp.xyz, atol=1e-5)

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

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_integrate_compounds_only_selected_moves(self):
        """Test that only integrate_compounds move while the rest stay fixed."""
        methane = mb.load("C", smiles=True)
        methane.name = "Methane"
        ethanol = mb.load("CCO", smiles=True)
        ethanol.name = "Ethanol"
        box = mb.fill_box(
            compound=[methane, ethanol], n_compounds=[10, 1], box=[3, 3, 3]
        )
        ethanol_compound = box.children[-1]
        sim = HoomdSimulation(
            system=box,
            integrate_compounds=[ethanol_compound],
            r_cut=1.0,
            box_buffer=2,
            run_on_gpu=False,
        )

        ethanol_ids = {id(p) for p in ethanol_compound.particles()}
        methane_particles = [p for p in box.particles() if id(p) not in ethanol_ids]
        methane_before = np.array([np.copy(p.pos) for p in methane_particles])
        ethanol_before = np.copy(ethanol_compound.xyz)

        sim.nvt(n_steps=50000, kT=2, dt=1e-4, tau=1e-2)

        methane_after = np.array([p.pos for p in methane_particles])
        ethanol_after = ethanol_compound.xyz

        # Non-integrated methane coords stay fixed
        assert np.allclose(methane_before, methane_after, atol=1e-2)
        # Integrated ethanol coords change
        assert not np.allclose(ethanol_before, ethanol_after, atol=1e-2)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fixed_compounds_anchor_stays_put(self):
        """Test that a fixed_compounds particle does not move during a run."""
        chain = mb.load("CCCCCCC", smiles=True)
        chain.translate_to((0, 0, 0))
        sim = HoomdSimulation(
            system=chain, fixed_compounds=[3], r_cut=1.0, box_buffer=4, run_on_gpu=False
        )

        anchor_pos_before = np.copy(chain.xyz[3])
        other_pos_before = np.copy(chain.xyz[0])

        sim.nvt(n_steps=50000, kT=2, dt=1e-4, tau=1e-2)

        anchor_pos_after = chain.xyz[3]
        other_pos_after = chain.xyz[0]

        # The fixed anchor particle stays put
        assert np.allclose(anchor_pos_before, anchor_pos_after, atol=1e-2)
        # A non-fixed particle moves
        assert not np.allclose(other_pos_before, other_pos_after, atol=1e-2)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_positions_returned_in_original_frame(self):
        """Regression: sync-back must undo HOOMD's box-centering translation."""
        methane = mb.load("C", smiles=True)
        box = mb.fill_box(compound=[methane], n_compounds=[5], box=[3, 3, 3])
        frozen_indices = list(range(box.n_particles - 1))
        sim = HoomdSimulation(
            system=box,
            fixed_compounds=frozen_indices,
            r_cut=1.0,
            box_buffer=2,
            run_on_gpu=False,
        )
        frozen_before = np.copy(sim.compound.xyz[frozen_indices])
        sim.nvt(n_steps=100, kT=2, dt=1e-4, tau=1e-2)
        frozen_after = sim.compound.xyz[frozen_indices]
        # The fixed particles do not move, so their absolute coordinates must be
        assert np.allclose(frozen_before, frozen_after, atol=1e-6)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_box_update_bounding_box_within_target(self):
        """box_update should leave the compound fitting inside the target box."""
        methane = mb.load("C", smiles=True)
        methane.name = "Methane"
        box = mb.fill_box(compound=[methane], n_compounds=[8], box=[3, 3, 3])
        sim = HoomdSimulation(system=box, r_cut=1.0, box_buffer=2, run_on_gpu=False)

        target = hoomd.Box(5, 5, 5)
        sim.box_update(
            n_steps=3000,
            kT=1.0,
            dt=1e-4,
            tau=1e-2,
            target_box=target,
            update_period=50,
        )

        bbox_lengths = np.asarray(sim.compound.get_boundingbox().lengths)
        assert np.all(bbox_lengths <= np.asarray(target.L))

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_box_update_no_half_box_shift(self):
        """Resizing must adjust the frame offset by dL/2, not leave it stale."""
        methane = mb.load("C", smiles=True)
        methane.name = "Methane"
        box = mb.fill_box(compound=[methane], n_compounds=[8], box=[3, 3, 3])
        sim = HoomdSimulation(system=box, r_cut=1.0, box_buffer=2, run_on_gpu=False)

        sim.box_update(
            n_steps=3000,
            kT=1.0,
            dt=1e-4,
            tau=1e-2,
            target_box=hoomd.Box(5, 5, 5),
            update_period=50,
        )

        L_final = np.asarray(sim.state.box.L)
        xyz = sim.compound.xyz
        # Box corner is at the origin: all coords lie within [0, L_final].
        # A stale offset would shift these by ~dL/2 (1 nm here) and go negative.
        assert np.all(xyz >= -1e-4)
        assert np.all(xyz <= L_final + 1e-4)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_recover_undoes_only_the_last_run(self):
        """recover() rewinds the last run's coordinates and drops its energies."""
        octane = mb.load("CCCCCCCC", smiles=True)
        sim = HoomdSimulation(system=octane, box=[3, 3, 3], run_on_gpu=False)

        sim.fire(n_steps=10)
        after_fire = np.copy(sim.compound.xyz)
        n_energies = len(sim.energies)

        sim.cap_displacement(n_steps=10, dt=1e-3)
        assert not np.allclose(sim.compound.xyz, after_fire)

        sim.recover()
        assert np.allclose(sim.compound.xyz, after_fire, atol=1e-5)
        assert len(sim.energies) == n_energies

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_recover_rewinds_cached_snapshot(self):
        """A Simulation built after recover() starts from the restored coords."""
        octane = mb.load("CCCCCCCC", smiles=True)
        sim = HoomdSimulation(system=octane, box=[3, 3, 3], run_on_gpu=False)
        before = np.copy(sim.compound.xyz)

        sim.fire(n_steps=10)
        sim.recover()

        reused = HoomdSimulation(system=octane, box=[3, 3, 3], kick=False)
        assert np.allclose(reused.compound.xyz, before, atol=1e-5)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fixed_compounds_accepts_numpy_indices(self):
        """Numpy integer indices work the same as built-in ints."""
        chain = mb.load("CCCCCCCC", smiles=True)
        sim = HoomdSimulation(
            system=chain,
            fixed_compounds=np.arange(3),
            box_buffer=4,
            run_on_gpu=False,
        )
        sim.nvt(n_steps=10, kT=1.0, dt=1e-4, tau=1e-2)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_integrate_compounds_mixed_types_raises(self, octane):
        """A mix of Compounds and indices is rejected with a clear error."""
        with pytest.raises(TypeError, match="integrate_compounds"):
            HoomdSimulation(
                system=octane,
                integrate_compounds=[list(octane.particles())[0], 1],
                box_buffer=4,
                run_on_gpu=False,
            )

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_box_update_syncs_compound_box(self):
        """The compound's box must track the resized simulation box."""
        methane = mb.load("C", smiles=True)
        methane.name = "Methane"
        box = mb.fill_box(compound=[methane], n_compounds=[8], box=[3, 3, 3])
        sim = HoomdSimulation(system=box, r_cut=1.0, box_buffer=2, run_on_gpu=False)

        sim.box_update(
            n_steps=3000,
            kT=1.0,
            dt=1e-4,
            tau=1e-2,
            target_box=hoomd.Box(5, 5, 5),
            update_period=50,
        )

        assert np.allclose(sim.compound.box.lengths, sim.state.box.L, atol=1e-4)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_nvt_leaves_compound_box_unchanged(self):
        """Runs that do not resize must leave the compound's box alone."""
        methane = mb.load("C", smiles=True)
        methane.name = "Methane"
        box = mb.fill_box(compound=[methane], n_compounds=[8], box=[3, 3, 3])
        sim = HoomdSimulation(system=box, r_cut=1.0, box_buffer=2, run_on_gpu=False)
        lengths_before = np.asarray(sim.compound.box.lengths)

        sim.nvt(n_steps=100, kT=1.0, dt=1e-4, tau=1e-2)

        assert np.allclose(sim.compound.box.lengths, lengths_before, atol=1e-4)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_second_sim_from_same_compound(self, oplsaa):
        """A compound can be re-simulated after its sim data has been cached."""
        octane = mb.load("CCCCCCCC", smiles=True)
        sim = HoomdSimulation(system=octane, box=[3, 3, 3], run_on_gpu=False)
        sim.fire(n_steps=10)
        reused = HoomdSimulation(system=octane, box=[3, 3, 3], run_on_gpu=False)
        assert reused.forces is sim.forces

        rebuilt = HoomdSimulation(
            system=octane, forcefield=oplsaa, box=[3, 3, 3], run_on_gpu=False
        )
        assert rebuilt.forces is not sim.forces

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_cached_sim_data_respects_box_and_r_cut(self):
        """Changing the box or cutoff must not silently reuse cached sim data."""
        octane = mb.load("CCCCCCCC", smiles=True)
        sim = HoomdSimulation(system=octane, box=[3, 3, 3], run_on_gpu=False)

        bigger = HoomdSimulation(system=octane, box=[9, 9, 9], run_on_gpu=False)
        assert np.allclose(bigger.state.box.L, [9, 9, 9])

        wider_cut = HoomdSimulation(
            system=octane, box=[9, 9, 9], r_cut=2.5, run_on_gpu=False
        )
        assert wider_cut.forces is not bigger.forces
        assert wider_cut.forces is not sim.forces

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_relax_path(self):
        from mbuild.path.build import Path

        bond_length = 1  # target, starting at lengths of (2,1)
        bead_radius = 0.5

        xcoords = np.concat((np.arange(0, 10, 2), np.arange(1, 10, 4)))
        xcoords.sort()
        coords = np.zeros((len(xcoords), 3))
        coords[:, 0] = xcoords
        path = Path(coords)
        path.form_linear_bond_graph()
        with pytest.raises(RuntimeError):
            path.relax(bead_radius=bead_radius, bond_length=bond_length, btype="fene")
        path.relax(bead_radius=bead_radius, bond_length=bond_length)
        bond_lengths = [
            np.linalg.norm(path.coordinates[i] - path.coordinates[i + 1])
            for i in range(len(path) - 1)
        ]
        assert path
        assert max(bond_lengths) < bond_length + 0.05  # converges to target


class TestOpenMMSimulation(BaseTest):
    """Tests for the OpenMMSimulation class."""

    @pytest.fixture
    def simple_compound(self):
        return _make_simple_compound(n=5, spacing=0.15)

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

    def test_init_with_compound(self, simple_compound):
        """Test initialization with a Compound."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        assert sim._is_path is False
        assert sim.compound is simple_compound
        assert sim.system is not None
        assert sim.topology is not None

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
                    force.getNonbondedMethod() == openmm.CustomNonbondedForce.NoCutoff
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

    def test_minimize_compound(self, simple_compound):
        """Test energy minimization on a Compound."""
        old_coords = simple_compound.xyz.copy()
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=100)
        assert not np.allclose(old_coords, simple_compound.xyz, atol=1e-6)

    def test_minimize_records_energy(self, simple_compound):
        """Test that minimize records energy."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        energy = sim.get_energy()
        assert energy is not None
        assert "potential_energy" in energy
        assert len(energy["potential_energy"]) == 1

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

    def test_nvt_with_report_interval(self, simple_compound):
        """Test NVT records an energy frame per reporting chunk."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.nvt(n_steps=100, T=300, dt=0.0001, report_interval=25)
        assert len(sim.get_energy()["potential_energy"]) == 4

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_nvt_with_forcefield(self, octane):
        """Test NVT with a real forcefield."""
        octane.box = mb.Box([5, 5, 5])
        sim = OpenMMSimulation(octane, forcefield="oplsaa")
        sim.minimize(n_steps=200)
        old_coords = octane.xyz.copy()
        sim.nvt(n_steps=50, T=300, dt=0.0001)
        assert not np.allclose(old_coords, octane.xyz, atol=1e-6)

    def test_get_energy_empty(self, simple_compound):
        """Test get_energy returns None before any run."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        assert sim.get_energy() is None

    def test_minimize_then_nvt(self, simple_compound):
        """Test calling minimize then nvt sequentially."""
        sim = OpenMMSimulation(simple_compound, forcefield=None)
        sim.minimize(n_steps=50)
        coords_after_min = simple_compound.xyz.copy()
        sim.nvt(n_steps=50, T=300, dt=0.0001)
        assert not np.allclose(coords_after_min, simple_compound.xyz, atol=1e-6)
        energy = sim.get_energy()
        assert len(energy["potential_energy"]) == 2
