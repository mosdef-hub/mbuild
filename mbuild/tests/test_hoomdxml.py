import mbuild as mb
import numpy as np
from mbuild.tests.base_test import BaseTest


class TestHoomdXML(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.hoomdxml')

    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.hoomdxml',forcefield='opls')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0,2.0,2.0]))
        ethane.save(filename='ethane-box.hoomdxml',forcefield='opls',box=box)
