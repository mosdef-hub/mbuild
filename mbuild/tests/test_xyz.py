import numpy as np
import pytest

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestXYZ(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.xyz')
