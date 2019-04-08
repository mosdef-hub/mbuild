import difflib
import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn
from mbuild.utils.validation import assert_port_exists


class TestUtils(BaseTest):

    def test_assert_port_exists(self, ch2):
        assert_port_exists('up', ch2)
        with pytest.raises(ValueError):
            assert_port_exists('dog', ch2)

    def test_structure_reproducibility(self, alkane_monolayer):
        filename = 'monolayer-tmp.pdb'
        alkane_monolayer.save(filename)
        with open(get_fn('monolayer.pdb')) as file1:
            with open('monolayer-tmp.pdb') as file2:
                diff = difflib.ndiff(file1.readlines(), file2.readlines())
        changes = [l for l in diff if l.startswith('+ ') or l.startswith('- ')]
        assert not changes
