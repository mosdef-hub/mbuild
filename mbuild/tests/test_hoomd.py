import xml.etree.ElementTree

import numpy as np
import packaging.version
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer, has_gsd, has_hoomd, import_

if has_hoomd:
    import hoomd

    if "version" in dir(hoomd):
        hoomd_version = packaging.version.parse(hoomd.version.version)
    else:
        hoomd_version = packaging.version.parse(hoomd.__version__)


@pytest.mark.skipif(not has_hoomd, reason="HOOMD is not installed")
class TestHoomdAny(BaseTest):
    def test_empty_initial_snapshot(self):
        import hoomd

        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        part = mb.Compound(name="Ar")
        box = mb.Box(lengths=[5, 5, 5], angles=[90, 90, 90])
        system = mb.fill_box(part, n_compounds=10, box=box)

        if hoomd_version.major == 2:
            hoomd.context.initialize("")
            init_snap = hoomd.data.make_snapshot(
                N=0, box=hoomd.data.boxdim(L=10)
            )
        else:
            init_snap = hoomd.Snapshot()
            init_snap.configuration.box = hoomd.Box.cube(L=10)

        with pytest.raises(RuntimeError):
            snap, _ = to_hoomdsnapshot(system, hoomd_snapshot=init_snap)

    def test_compound_from_snapshot(self, ethane):
        from mbuild.formats.hoomd_snapshot import (
            from_snapshot,
            to_hoomdsnapshot,
        )

        lengths = [5, 5, 5]
        filled = mb.fill_box(ethane, n_compounds=5, box=mb.Box(lengths))
        snap, _ = to_hoomdsnapshot(filled)
        new_filled = from_snapshot(snap, scale=0.1)

        assert filled.n_bonds == new_filled.n_bonds
        assert filled.n_particles == new_filled.n_particles

        assert np.array_equal(filled.box.angles, new_filled.box.angles)
        assert np.array_equal(filled.box.lengths, new_filled.box.lengths)

        for i in range(filled.n_particles):
            assert np.allclose(filled[i].pos, new_filled[i].pos)

    @pytest.mark.skipif(not has_gsd, reason="gsd is not installed")
    def test_compound_from_gsdsnapshot(self, ethane):
        import gsd.hoomd

        from mbuild.formats.hoomd_snapshot import (
            from_snapshot,
            to_hoomdsnapshot,
        )

        lengths = [5, 5, 5]
        filled = mb.fill_box(ethane, n_compounds=5, box=mb.Box(lengths))
        snap, _ = to_hoomdsnapshot(filled)

        # copy attributes from the snapshot to a gsd snapshot
        gsd_snap = gsd.hoomd.Frame()
        gsd_snap.particles.N = snap.particles.N
        gsd_snap.particles.types = snap.particles.types
        gsd_snap.particles.typeid = snap.particles.typeid
        gsd_snap.particles.position = snap.particles.position
        if hoomd_version.major == 2:
            gsd_snap.configuration.box = np.array(
                [
                    snap.box.Lx,
                    snap.box.Ly,
                    snap.box.Lz,
                    snap.box.xy,
                    snap.box.xy,
                    snap.box.yz,
                ]
            )
        else:
            gsd_snap.configuration.box = snap.configuration.box

        gsd_snap.bonds.N = snap.bonds.N
        gsd_snap.bonds.group = snap.bonds.group
        gsd_snap.particles.charge = snap.particles.charge
        gsd_snap.validate()

        new_filled = from_snapshot(gsd_snap, scale=0.1)

        assert filled.n_bonds == new_filled.n_bonds
        assert filled.n_particles == new_filled.n_particles

        assert np.array_equal(filled.box.angles, new_filled.box.angles)
        assert np.array_equal(filled.box.lengths, new_filled.box.lengths)

        for i in range(filled.n_particles):
            assert np.allclose(filled[i].pos, new_filled[i].pos)

    def test_compound_to_snapshot(self, ethane):
        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        snap, _ = to_hoomdsnapshot(ethane)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 0

    def test_particles_to_snapshot(self):
        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        part = mb.Compound(name="Ar")
        box = mb.Box(lengths=[5, 5, 5], angles=[90, 90, 90])
        system = mb.fill_box(part, n_compounds=10, box=box)
        snap, _ = to_hoomdsnapshot(system)

        assert snap.particles.N == 10
        assert snap.bonds.N == 0
        assert snap.angles.N == 0

    def test_snapshot_from_initial(self):
        import hoomd

        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        part = mb.Compound(name="Ar")
        box = mb.Box(lengths=[5, 5, 5], angles=[90, 90, 90])
        system = mb.fill_box(part, n_compounds=10, box=box)
        if hoomd_version.major == 2:
            hoomd.context.initialize("")
            init_snap = hoomd.data.make_snapshot(
                N=10, box=hoomd.data.boxdim(L=10)
            )
        else:
            init_snap = hoomd.Snapshot()
            init_snap.particles.N = 10
            init_snap.configuration.box = hoomd.Box.cube(L=10)

        snap, _ = to_hoomdsnapshot(system, hoomd_snapshot=init_snap)

        assert snap.particles.N == 20
        assert snap.bonds.N == 0
        assert snap.angles.N == 0
        if hoomd_version.major == 2:
            assert (snap.box.Lx, snap.box.Ly, snap.box.Lz) == (50, 50, 50)
            assert (snap.box.xy, snap.box.xz, snap.box.yz) == (0, 0, 0)
        else:
            np.testing.assert_allclose(
                snap.configuration.box, [50, 50, 50, 0, 0, 0]
            )

    def test_bad_input_to_snapshot(self):
        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        with pytest.raises(ValueError):
            to_hoomdsnapshot("fake_object")

    def test_non_param_struc_to_snapshot(self, ethane):
        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        structure = ethane.to_parmed()
        snap, _ = to_hoomdsnapshot(structure)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 0

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_param_structure_to_snapshot(self, ethane):
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_snapshot import to_hoomdsnapshot

        ff = Forcefield(name="oplsaa")
        structure = ff.apply(ethane)
        snap, _ = to_hoomdsnapshot(structure)

        assert snap.particles.N == 8
        assert snap.bonds.N == 7
        assert snap.angles.N == 12
        assert snap.dihedrals.N == 9
        assert snap.pairs.N == 9


@pytest.mark.skipif(
    not has_hoomd or hoomd_version.major != 2,
    reason="HOOMD is not installed or is the wrong version",
)
class TestHoomdSimulation(BaseTest):
    def test_bad_input_to_hoomdsimulation(self):
        from mbuild.formats.hoomd_simulation import create_hoomd_simulation

        with pytest.raises(ValueError):
            create_hoomd_simulation("fake_object", 2)

    def test_compound_to_hoomdsimulation(self, ethane):
        from mbuild.formats.hoomd_simulation import create_hoomd_simulation

        with pytest.raises(ValueError):
            create_hoomd_simulation(ethane, 2.5)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_structure_to_hoomdsimulation(self, ethane):
        import hoomd
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_simulation import create_hoomd_simulation

        ff = Forcefield(name="oplsaa")
        structure = ff.apply(ethane)
        sim = hoomd.context.SimulationContext()
        with sim:
            create_hoomd_simulation(structure, 2.5)

            sim_forces = hoomd.context.current.forces
            pair_force = import_("hoomd.md.pair")
            charge_force = import_("hoomd.md.charge")
            special_pair_force = import_("hoomd.md.special_pair")
            bond_force = import_("hoomd.md.bond")
            angle_force = import_("hoomd.md.angle")
            dihedral_force = import_("hoomd.md.dihedral")

            assert isinstance(sim_forces[0], pair_force.lj)
            assert isinstance(sim_forces[1], charge_force.pppm)
            assert isinstance(sim_forces[2], pair_force.ewald)
            assert isinstance(sim_forces[3], special_pair_force.lj)
            assert isinstance(sim_forces[4], special_pair_force.coulomb)
            assert isinstance(sim_forces[5], bond_force.harmonic)
            assert isinstance(sim_forces[6], angle_force.harmonic)
            assert isinstance(sim_forces[7], dihedral_force.opls)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_lj_to_hoomdsimulation(self):
        import hoomd
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_simulation import create_hoomd_simulation

        box = mb.Compound()
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        ff = Forcefield(forcefield_files=get_fn("lj.xml"))
        structure = ff.apply(box)
        structure.box = [10, 10, 10, 90, 90, 90]
        sim = hoomd.context.SimulationContext()
        with sim:
            create_hoomd_simulation(structure, 2.5)
            sim_forces = hoomd.context.current.forces
            pair_force = import_("hoomd.md.pair")

            assert isinstance(sim_forces[0], pair_force.lj)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    @pytest.mark.skipif(not has_gsd, reason="GSD is not installed")
    def test_hoomdsimulation_restart(self):
        import gsd.hoomd
        import hoomd
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_simulation import create_hoomd_simulation

        box = mb.Compound()
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        ff = Forcefield(forcefield_files=get_fn("lj.xml"))
        structure = ff.apply(box)
        structure.box = [10, 10, 10, 90, 90, 90]
        sim = hoomd.context.SimulationContext()
        with sim:
            hoomd_obj, ref_vals = create_hoomd_simulation(
                structure, 2.5, restart=get_fn("restart.gsd")
            )
            sim_forces = hoomd.context.current.forces
            pair_force = import_("hoomd.md.pair")

            assert isinstance(sim_forces[0], pair_force.lj)

        snap = hoomd_obj[0]
        with gsd.hoomd.open(get_fn("restart.gsd")) as f:
            rsnap = f[0]
        assert np.array_equal(snap.particles.position, rsnap.particles.position)

    def test_hoomdsimulation_nlist(self, ethane):
        hoomd_simulation = import_("mbuild.formats.hoomd_simulation")
        hoomd = import_("hoomd")
        hoomd.md = import_("hoomd.md")

        with pytest.raises(ValueError):
            hoomd_simulation.create_hoomd_simulation(
                ethane, 2.5, nlist=hoomd.md.nlist.tree
            )


@pytest.mark.skipif(
    not has_hoomd or hoomd_version.major < 3,
    reason="HOOMD is not installed or is the wrong version",
)
class TestHoomdForcefield(BaseTest):
    def test_bad_input_to_hoomd_forcefield(self):
        from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

        with pytest.raises(ValueError):
            create_hoomd_forcefield("fake_object", 2.5)

    def test_compound_to_hoomd_forcefield(self, ethane):
        from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

        with pytest.raises(ValueError):
            create_hoomd_forcefield(ethane, 2.5)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_structure_to_hoomd_forcefield(self, ethane):
        import hoomd
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

        ff = Forcefield(name="oplsaa")
        structure = ff.apply(ethane)

        snapshot, forces, ref_values = create_hoomd_forcefield(structure, 2.5)

        assert isinstance(forces[0], hoomd.md.pair.LJ)
        assert isinstance(forces[1], hoomd.md.pair.Ewald)
        assert isinstance(forces[2], hoomd.md.long_range.pppm.Coulomb)
        assert isinstance(forces[3], hoomd.md.special_pair.LJ)
        assert isinstance(forces[4], hoomd.md.special_pair.Coulomb)
        assert isinstance(forces[5], hoomd.md.bond.Harmonic)
        assert isinstance(forces[6], hoomd.md.angle.Harmonic)
        assert isinstance(forces[7], hoomd.md.dihedral.OPLS)

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_lj_to_hoomd_forcefield(self):
        import hoomd
        from foyer.forcefield import Forcefield

        from mbuild.formats.hoomd_forcefield import create_hoomd_forcefield

        box = mb.Compound()
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        box.add(mb.Compound(name="Ar", pos=[1, 1, 1]))
        ff = Forcefield(forcefield_files=get_fn("lj.xml"))
        structure = ff.apply(box)
        structure.box = [10, 10, 10, 90, 90, 90]

        snapshot, forces, ref_values = create_hoomd_forcefield(structure, 2.5)

        assert isinstance(forces[0], hoomd.md.pair.LJ)
