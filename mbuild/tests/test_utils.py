import pytest
from mbuild.tests.base_test import BaseTest
from mbuild.utils.validation import assert_port_exists


class TestUtils(BaseTest):

    def test_assert_port_exists(self, ch2):
        assert_port_exists('up', ch2)
        with pytest.raises(ValueError):
            assert_port_exists('dog', ch2)
