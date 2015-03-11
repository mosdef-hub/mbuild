import pytest
import numpy as np

import mbuild as mb
from mbuild.examples.ethane.ethane import Ethane
from mbuild.testing.tools import get_fn

from mbuild.tests.base_test import BaseTest


class TestHoomdXml(BaseTest):

    @pytest.fixture
    def molecule(self):
        lj_units = {'mass': 72.0,
                    'distance': 0.6, 
                    'energy': 0.4}
        compound = mb.load(get_fn('ecer2.hoomdxml'), lj_units=lj_units)
        return compound

    def test_write(self, molecule):
        molecule.save('ecer2-saved.hoomdxml')

    def test_update_from_file(self):
        ethane = Ethane()
        ethane.update_coordinates(get_fn("ethane.hoomdxml"))

    def test_units(self, molecule):
        positions = np.array([
             [-0.855480134487, -3.11021924019, -0.716033160686],
             [-0.907484650612, -3.15630435944, -0.192944049835],
             [-0.875647962093, -3.17938971519, 0.370299994946],
             [-0.89426112175, -3.14589977264, 0.94112187624],
             [-0.907657623291, -3.27201843262, 1.52672827244],
             [-0.828304767609, -3.54459166527, 2.13118267059],
             [-0.493888318539, -4.01740121841, 2.17989301682],
             [-0.169051453471, -3.64267301559, 2.06289696693],
             [-0.0267458446324, -3.55486989021, 1.46613144875],
             [-0.126119762659, -3.60231733322, 0.928843021393],
             [0.0277806278318, -3.6437664032, 0.310265809298],
             [0.00889551546425, -3.54814314842, -0.252536147833],
             [-0.121015898883, -3.50847649574, -0.732142329216]])
        loaded = molecule.to_trajectory().xyz / 0.6
        assert loaded.all() == positions.all()
