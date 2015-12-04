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
        assert ethane.n_particles == 8
        assert ethane.n_bonds == 7

    def test_methane(self):
        import mbuild.examples.methane.methane as example
        methane = example.main()

        assert len(methane.children) == 5
        assert methane.n_particles == 5

        assert len(methane.bond_graph) == 5
        assert methane.n_bonds == 4

        assert next(methane.particles()).name == 'C'
        assert methane[0].name == 'C'
        assert (methane[1].name == methane[2].name ==
                methane[3].name == methane[4].name == "H")


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
