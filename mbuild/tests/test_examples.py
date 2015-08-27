from mbuild.tests.base_test import BaseTest


class TestExamples(BaseTest):

    def test_alkane(self):
        import mbuild.examples.alkane.alkane as example
        example.main()

    def test_alkane_monolayer(self):
        import mbuild.examples.alkane_monolayer.alkane_monolayer as example
        example.main()

    def test_bilayer(self):
        import mbuild.examples.bilayer.bilayer as example
        example.main()

    def test_ethane(self):
        import mbuild.examples.ethane.ethane as example
        ethane = example.main()
        assert ethane.n_atoms == 8
        assert ethane.n_bonds == 7

    def test_methane(self):
        import mbuild.examples.methane.methane as example
        methane = example.main()

        assert methane.n_atoms == 5
        assert methane.n_bonds == 4

        assert methane.atoms[0].name == "C"
        assert methane.HC[0].name == "H"
        assert methane.HC[1].name == "H"
        assert methane.HC[2].name == "H"
        assert methane.HC[3].name == "H"

    def test_pmpc(self):
        import mbuild.examples.pmpc.pmpc_brush_layer as example
        example.main()

    def test_reload(self):
        import mbuild.examples.reload.reload as example
        example.main()

    def test_solvate(self):
        import mbuild.examples.solvate.solvate as example
        example.main()

    def test_tnp(self):
        import mbuild.examples.tnp.tnp as example
        example.main()

    def test_tnp_box(self):
        import mbuild.examples.tnp.tnp_box as example
        example.main()
