from mbuild.tests.base_test import BaseTest


class TestGromacs(BaseTest):

    def test_write(self, ethane):
        ethane.save('ethane.top')

    def test_update_from_file(self, ethane):
        pass

    def test_read_write_compare(self, ethane):
        pass

if __name__ == "__main__":
    from mbuild.examples.ethane.ethane import Ethane
    et = Ethane()
    TestGromacs().test_write(et)
