import pytest
import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestPattern(BaseTest):

    def test_apply_to_compound(self, betacristobalite, alkyl, ch3):
        pattern = mb.Random2DPattern(90)
        chains, backfills = pattern.apply_to_compound(
            guest=alkyl, host=betacristobalite, backfill=ch3)
        assert len(chains) == 90
        assert len(backfills) == 10

        with pytest.raises(AssertionError):
            pattern = mb.Random2DPattern(101)
            chains, backfills = pattern.apply_to_compound(
                guest=alkyl, host=betacristobalite, backfill=ch3)
