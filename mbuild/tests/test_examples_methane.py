import pytest

__author__ = 'sallai'


class TestExamplesMethane:

    @pytest.fixture
    def methane(self):
        from mbuild.examples.methane.methane import Methane
        return Methane()

    def test_n_atoms(self, methane):
        assert methane.n_atoms() == 5
        assert methane.n_bonds() == 4

    def test_atom_types(self, methane):
        assert methane.c.kind == "C"
        assert methane.h1.kind == "H"
        assert methane.h2.kind == "H"
        assert methane.h3.kind == "H"
        assert methane.h4.kind == "H"

