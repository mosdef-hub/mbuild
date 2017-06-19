import numpy as np
import pytest
import xml.etree.ElementTree

import mbuild as mb
from mbuild.tests.base_test import BaseTest


class TestHoomdXML(BaseTest):

    def test_save(self, ethane):
        ethane.save(filename='ethane.hoomdxml')

    def test_save_forcefield(self, ethane):
        ethane.save(filename='ethane-opls.hoomdxml', forcefield_name='oplsaa')

    def test_save_box(self, ethane):
        box = mb.Box(lengths=np.array([2.0, 2.0, 2.0]))
        ethane.save(filename='ethane-box.hoomdxml', box=box)

    def test_rigid(self, benzene):
        n_benzenes = 10
        benzene.name = 'Benzene'
        filled = mb.fill_box(benzene,
                             n_compounds=n_benzenes,
                             box=[0, 0, 0, 4, 4, 4]) 
        filled.label_rigid_bodies(discrete_bodies='Benzene', rigid_particles='C')
        filled.save(filename='benzene.hoomdxml')

        xml_file = xml.etree.ElementTree.parse('benzene.hoomdxml').getroot()
        body_text = xml_file[0].find('body').text
        rigid_bodies = [int(body) for body in body_text.split('\n') if body]
        for body_id in range(10):
            assert rigid_bodies.count(body_id) == 6
        assert rigid_bodies.count(-1) == n_benzenes * 6

    def test_number_in_each_section(self, box_of_benzenes):
        box_of_benzenes.save(filename='benzene.hoomdxml', forcefield_name='oplsaa')
        xml_file = xml.etree.ElementTree.parse('benzene.hoomdxml').getroot()
        for attribute in ['position', 'type', 'mass', 'charge']: 
            body_text = xml_file[0].find(attribute).text
            list_of_things = [x for x in body_text.split('\n') if x]
            assert len(list_of_things) == 12*10
        for attribute, number in [('bond', 12), ('angle', 18), ('dihedral', 24)]:
            body_text = xml_file[0].find(attribute).text
            list_of_things = [x for x in body_text.split('\n') if x]
            assert len(list_of_things) == number*10
