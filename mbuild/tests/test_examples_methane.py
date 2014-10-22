import pytest

__author__ = 'sallai'


class TestExamplesMethane:

    @pytest.fixture
    def methane(self):
        from mbuild.examples.methane.methane import Methane
        return Methane()

    def test_numbers(self, methane):
        assert methane.n_atoms == 5
        assert methane.n_bonds == 4

    def test_atom_types(self, methane):
        assert methane.atoms[0].kind == "C"
        assert methane.hc[0].kind == "H"
        assert methane.hc[1].kind == "H"
        assert methane.hc[2].kind == "H"
        assert methane.hc[3].kind == "H"

