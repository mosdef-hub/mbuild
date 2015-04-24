import pytest


class BaseTest:

    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()

    @pytest.fixture
    def ethane(self):
        from mbuild.examples.ethane.ethane import Ethane
        return Ethane()

    @pytest.fixture
    def methane(self):
        from mbuild.examples.methane.methane import Methane
        return Methane()

    @pytest.fixture
    def h2o(self):
        from mbuild.components.small_groups.h2o import H2O
        return H2O()

    @pytest.fixture
    def ch2(self):
        from mbuild.components.small_groups.ch2 import CH2
        return CH2()

    @pytest.fixture
    def betacristobalite(self):
        from mbuild.components.surfaces.betacristobalite import Betacristobalite
        return Betacristobalite()
