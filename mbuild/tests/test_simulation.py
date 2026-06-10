import sys

import hoomd
import numpy as np
import pytest
from gmso import ForceField

import mbuild as mb
from mbuild.exceptions import MBuildError
from mbuild.simulation import (
    ForcesHandler,
    HoomdSimulation,
    energy_minimize,
    hoomd_cap_displacement,
    hoomd_fire,
)
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, has_foyer, has_hoomd, has_openbabel


class TestSimulation(BaseTest):
    def test_energy_minimize(self, octane):
        energy_minimize(compound=octane)

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_shift_com(self, octane):
        com_old = octane.pos
        energy_minimize(compound=octane)
        # check to see if COM of energy minimized Compound
        # has been shifted back to the original COM
        assert np.allclose(com_old, octane.pos)

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_shift_anchor(self, octane):
        anchor_compound = octane.labels["chain"].labels["CH3[0]"]
        pos_old = anchor_compound.pos
        energy_minimize(compound=octane, anchor=anchor_compound)
        # check to see if COM of the anchor Compound
        # has been shifted back to the original COM
        assert np.allclose(pos_old, anchor_compound.pos)

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_fix_compounds(self, octane):
        methyl_end0 = octane.labels["chain"].labels["CH3[0]"]
        methyl_end1 = octane.labels["chain"].labels["CH3[0]"]
        carbon_end = octane.labels["chain"].labels["CH3[0]"].labels["C[0]"]
        not_in_compound = mb.Compound(name="H")

        # fix the whole molecule and make sure positions are close
        # given stochastic nature and use of restraining springs
        # we need to have a pretty loose tolerance for checking
        old_com = octane.pos
        energy_minimize(
            compound=octane,
            fixed_compounds=octane,
            shift_com=False,
            constraint_factor=1e6,
        )
        assert np.allclose(octane.pos, old_com, rtol=1e-2, atol=1e-2)

        # primarily focus on checking inputs are parsed correctly
        energy_minimize(compound=octane, fixed_compounds=[octane])
        energy_minimize(compound=octane, fixed_compounds=carbon_end)
        energy_minimize(compound=octane, fixed_compounds=methyl_end0)
        energy_minimize(compound=octane, fixed_compounds=[methyl_end0])
        energy_minimize(
            compound=octane, fixed_compounds=[methyl_end0, (True, True, True)]
        )

        energy_minimize(
            compound=octane, fixed_compounds=[methyl_end0, (True, True, False)]
        )
        energy_minimize(
            compound=octane, fixed_compounds=[methyl_end0, [True, True, False]]
        )
        energy_minimize(
            compound=octane, fixed_compounds=[methyl_end0, (True, False, False)]
        )
        energy_minimize(
            compound=octane, fixed_compounds=[methyl_end0, (False, False, False)]
        )

        energy_minimize(compound=octane, fixed_compounds=[methyl_end0, methyl_end1])
        energy_minimize(compound=octane, fixed_compounds=[[methyl_end0], [methyl_end1]])
        energy_minimize(
            compound=octane,
            fixed_compounds=[
                [methyl_end0, (True, True, True)],
                [methyl_end1, (True, True, True)],
            ],
        )

        with pytest.raises(MBuildError):
            energy_minimize(compound=octane, fixed_compounds=not_in_compound)
        with pytest.raises(MBuildError):
            energy_minimize(compound=octane, fixed_compounds=[not_in_compound])
        with pytest.raises(MBuildError):
            energy_minimize(
                compound=octane, fixed_compounds=[12323.3, (True, False, False)]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane,
                fixed_compounds=[methyl_end0, (True, False, False, False)],
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, True, False, False]
            )
        with pytest.raises(Exception):
            energy_minimize(compound=octane, fixed_compounds=[methyl_end0, True])
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane,
                fixed_compounds=[methyl_end0, [True, False, False, False]],
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, (True, False)]
            )

        with pytest.raises(Exception):
            energy_minimize(compound=octane, fixed_compounds=[methyl_end0, (True)])

        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, ("True", True, True)]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, (True, "True", True)]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, (True, True, "True")]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, ("True", True, "True")]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, (True, "True", "True")]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, ("True", "True", True)]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, ("True", "True", "True")]
            )
        with pytest.raises(Exception):
            energy_minimize(
                compound=octane, fixed_compounds=[methyl_end0, (123.0, 231, "True")]
            )

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_ignore_compounds(self, octane):
        methyl_end0 = octane.labels["chain"].labels["CH3[0]"]
        methyl_end1 = octane.labels["chain"].labels["CH3[1]"]
        carbon_end = octane.labels["chain"].labels["CH3[0]"].labels["C[0]"]
        not_in_compound = mb.Compound(name="H")

        # fix the whole molecule and make sure positions are close
        # given stochastic nature and use of restraining springs
        # we need to have a pretty loose tolerance for checking
        old_com = octane.pos
        energy_minimize(
            compound=octane,
            ignore_compounds=octane,
            shift_com=False,
            constraint_factor=1e6,
        )
        assert np.allclose(octane.pos, old_com, rtol=1e-2, atol=1e-2)

        # primarily focus on checking inputs are parsed correctly
        energy_minimize(compound=octane, ignore_compounds=[octane])
        energy_minimize(compound=octane, ignore_compounds=carbon_end)
        energy_minimize(compound=octane, ignore_compounds=methyl_end0)
        energy_minimize(compound=octane, ignore_compounds=[methyl_end0])
        energy_minimize(compound=octane, ignore_compounds=[methyl_end0, methyl_end1])
        energy_minimize(
            compound=octane, ignore_compounds=[[methyl_end0], [methyl_end1]]
        )

        with pytest.raises(MBuildError):
            energy_minimize(compound=octane, ignore_compounds=not_in_compound)
        with pytest.raises(MBuildError):
            energy_minimize(compound=octane, ignore_compounds=[1231, 123124])

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_distance_constraints(self, octane):
        methyl_end0 = octane.labels["chain"].labels["CH3[0]"]
        methyl_end1 = octane.labels["chain"].labels["CH3[1]"]

        carbon_end0 = octane.labels["chain"].labels["CH3[0]"].labels["C[0]"]
        carbon_end1 = octane.labels["chain"].labels["CH3[1]"].labels["C[0]"]
        h_end0 = octane.labels["chain"].labels["CH3[0]"].labels["H[0]"]

        not_in_compound = mb.Compound(name="H")

        # given stochastic nature and use of restraining springs
        # we need to have a pretty loose tolerance for checking
        energy_minimize(
            compound=octane,
            distance_constraints=[(carbon_end0, carbon_end1), 0.7],
            constraint_factor=1e20,
        )
        assert np.allclose(
            np.linalg.norm(carbon_end0.pos - carbon_end1.pos),
            0.7,
            rtol=1e-2,
            atol=1e-2,
        )

        energy_minimize(
            compound=octane, distance_constraints=[[(carbon_end0, carbon_end1), 0.7]]
        )
        energy_minimize(
            compound=octane,
            distance_constraints=[
                [(carbon_end0, carbon_end1), 0.7],
                [(carbon_end0, h_end0), 0.1],
            ],
        )

        with pytest.raises(MBuildError):
            energy_minimize(
                compound=octane,
                distance_constraints=[(carbon_end0, not_in_compound), 0.7],
            )
        with pytest.raises(MBuildError):
            energy_minimize(
                compound=octane, distance_constraints=[(carbon_end0, carbon_end0), 0.7]
            )
        with pytest.raises(MBuildError):
            energy_minimize(
                compound=octane, distance_constraints=[(methyl_end0, carbon_end1), 0.7]
            )
        with pytest.raises(MBuildError):
            energy_minimize(
                compound=octane, distance_constraints=[(methyl_end0, methyl_end1), 0.7]
            )

    @pytest.mark.skipif(has_openbabel, reason="Open Babel package is installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_openbabel_warn(self, octane):
        with pytest.raises(MBuildError):
            energy_minimize(compound=octane)

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_ff(self, octane):
        for ff in ["UFF", "GAFF", "MMFF94", "MMFF94s", "Ghemical"]:
            energy_minimize(compound=octane, forcefield=ff)
        with pytest.raises(IOError):
            energy_minimize(compound=octane, forcefield="fakeFF")

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_algorithm(self, octane):
        for algorithm in ["cg", "steep", "md"]:
            energy_minimize(compound=octane, algorithm=algorithm)
        with pytest.raises(MBuildError):
            energy_minimize(compound=octane, algorithm="fakeAlg")

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
    )
    def test_energy_minimize_non_element(self, octane):
        for particle in octane.particles():
            particle.element = None
        # Pass with element inference from names
        energy_minimize(compound=octane)
        for particle in octane.particles():
            particle.name = "Q"
            particle.element = None

        # Fail once names cannot be set as elements
        with pytest.raises(MBuildError):
            energy_minimize(compound=octane)

    @pytest.mark.skipif(not has_openbabel, reason="Open Babel not installed")
    @pytest.mark.skipif(
        "win" in sys.platform, reason="Unknown issue with Window's Open Babel "
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

        energy_minimize(compound=octane)

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
        energy_minimize(compound=octane, forcefield="oplsaa")

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    @pytest.mark.parametrize("constraints", ["AllBonds", "HBonds", "HAngles", None])
    def test_energy_minimize_openmm_constraints(self, octane, constraints):
        energy_minimize(compound=octane, forcefield="oplsaa", constraints=constraints)

    def test_energy_minimize_openmm_invalid_constraints(self, octane):
        with pytest.raises(ValueError):
            energy_minimize(compound=octane, forcefield="oplsaa", constraints="boo")

    @pytest.mark.skipif(not has_foyer, reason="Foyer is not installed")
    def test_energy_minimize_openmm_xml(self, octane):
        energy_minimize(compound=octane, forcefield=get_fn("small_oplsaa.xml"))


class TestSimulationHoomd(BaseTest):
    """Test minimization methods using hoomd-blue."""

    @pytest.fixture
    def oplsaa(self):
        return ForceField("oplsaa")

    @pytest.fixture
    def sim(self, octane, oplsaa):
        return HoomdSimulation(octane, oplsaa, r_cut=0.3, run_on_gpu=False)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_forces_handler(self, sim):
        force = sim.get_force(hoomd.md.pair.LJ)
        assert len(force.params.keys()) == 6
        force = sim.get_force(hoomd.md.pair.Ewald)
        assert force
        force = sim.get_force(hoomd.md.long_range.pppm.Coulomb)
        assert force
        force = sim.get_force(hoomd.md.special_pair.LJ)
        assert len(force.params.keys()) == 5
        force = sim.get_force(hoomd.md.special_pair.Coulomb)
        assert force
        force = sim.get_force(hoomd.md.bond.Harmonic)
        assert len(force.params.keys()) == 2
        force = sim.get_force(hoomd.md.angle.Harmonic)
        assert len(force.params.keys()) == 3
        force = sim.get_force(hoomd.md.dihedral.OPLS)
        assert len(force.params.keys()) == 3

    def test_scale_forces(self, sim):
        A = 1
        ffhandler = ForcesHandler(dpd=A, scale_bond=0, scale_angle=0)
        ffhandler.scale_sim(sim)
        lj_force = sim.get_force(hoomd.md.pair.LJ)
        force = sim.get_force(hoomd.md.pair.DPDConservative, active_only=True)
        for param in force.params:
            assert force.params[param]["A"] == A
            assert force.r_cut[param] == lj_force.params[param]["sigma"]
        forceTypes = [type(force) for force in sim.active_forces]
        allforceTypes = [type(force) for force in sim.forces]
        assert hoomd.md.bond.Harmonic not in forceTypes
        assert hoomd.md.bond.Harmonic in allforceTypes
        assert hoomd.md.angle.Harmonic not in forceTypes
        assert hoomd.md.angle.Harmonic in allforceTypes
        assert hoomd.md.dihedral.OPLS not in forceTypes
        assert hoomd.md.dihedral.OPLS in allforceTypes

    def test_scale_opls(self, sim):
        ffhandler = ForcesHandler(scale_opls=0.5)
        ffhandler.scale_sim(sim)
        forceTypes = [type(force) for force in sim.active_forces]
        assert hoomd.md.dihedral.OPLS in forceTypes

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_fire(self, sim):
        ffhandler = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        cpd = sim.compound
        old_coords = cpd.xyz.copy()
        hoomd_fire(cpd, sim, ffhandler, n_steps=20)
        new_coords = cpd.xyz
        assert not np.allclose(old_coords, new_coords, atol=1e-5)
        # TODO: Get and validate energies here

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_capped_displacement(self, sim):
        bond = sim.get_force(hoomd.md.bond.Harmonic)
        orig_bond_params = {p: dict(bond.params[p]) for p in bond.params}

        ffhandler = ForcesHandler(dpd=1, scale_bond=0.1, scale_angle=0)
        cpd = sim.compound
        old_coords = cpd.xyz.copy()
        hoomd_cap_displacement(
            cpd, sim, ffhandler, dt=1, max_displacement=1, n_steps=10
        )
        new_coords = cpd.xyz
        assert not np.allclose(old_coords, new_coords, atol=1e-5)
        forcesTypesSet = set([type(force) for force in sim.active_forces])
        assert forcesTypesSet == set(
            (
                hoomd.md.pair.DPDConservative,
                hoomd.md.bond.Harmonic,
            )
        )
        for p in bond.params:
            assert np.isclose(bond.params[p]["k"], orig_bond_params[p]["k"] * 0.1)

    @pytest.mark.skipif(not has_hoomd, reason="hoomd is not installed")
    def test_scale_forces_restores_on_second_call(self, sim):
        bond = sim.get_force(hoomd.md.bond.Harmonic)
        angle = sim.get_force(hoomd.md.angle.Harmonic)
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
    def test_get_energy(self, sim):
        cpd = sim.compound
        ffhandler = ForcesHandler(scale_lj=1, scale_bond=1, scale_angle=1)
        assert sim.get_energy() is None
        hoomd_cap_displacement(
            cpd, sim, ffhandler, dt=1, max_displacement=1, n_steps=10
        )
        energyDict = sim.get_energy()
        assert np.allclose(
            energyDict["hoomd.md.pair.pair.LJ"], np.array(([-1.25498152, 0.0]))
        )
        assert np.allclose(
            energyDict["hoomd.md.bond.Harmonic"],
            np.array(([1.35651551e02, 3.05080224e06])),
            rtol=1e-2
        )
        assert np.allclose(
            energyDict["hoomd.md.angle.Harmonic"],
            np.array(([3942.05658645, 19129.72619001])),
            rtol=1e-2
        )
