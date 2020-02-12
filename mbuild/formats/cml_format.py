from lxml import etree as ET

import mbuild as mb
from mbuild.exceptions import MBuildError 
import numpy as np

def compound_from_cml(cml_file):
    """
    Convert the given cml file into a mb.Compound
    
    Given an input cml file, this method scans for the particles, bond information
    as well as other hierachical information regarding the compound and returns a
    mb.Compound

    Parameters
    ----------
    cml_file: (path, str) Path of the cml file

    Returns
    -------
    Compound: (mb.Compound) corresponding compound
    
    """

    # cml schema 3 has certain level of flexibiliy, so different software may have
    # slight variation in their schema. This initial code works with the cml save
    # out from Avogadro2, but I will work on making this work with other schema of cml
    # either do some type scanning or do partial string matching. 
    compound_tree = ET.parse(cml_file)
    compound_root = compound_tree.getroot()
    moleculeArray = compound_tree.find("{http://www.xml-cml.org/schema}atomArray")
    bondArray = compound_tree.find("{http://www.xml-cml.org/schema}bondArray")

    compound = mb.Compound()
    particles_dict = {}
    for molecule in moleculeArray:
        xyz = [molecule.attrib['x3'], molecule.attrib['y3'], molecule.attrib['z3']]
        particle = mb.Particle(name=molecule.attrib['elementType'], pos=xyz)
        particles_dict[molecule.attrib['id']] = particle
        compound.add(particle)

    for bond in bondArray:
        atom1, atom2 = bond.attrib['atomRefs2'].split()
        compound.add_bond((particles_dict[atom1], particles_dict[atom2]))
       
    return compound

 
def compound_to_cml(cmpd, file_path):
    """Convert the mb.Compound into equivelent json representation

    This method takes in the mb.Compound and tries to save the
    information of the mb.Compound into a cml file.
    Parameters
    ----------
    cmpd: mb.Compound
    file_path: str, path to save the cml file

    """
    root = ET.Element("molecules", nsmap={None:"http://www.xml-cml.org/schema",
                                          "cml":"http://www.xml-cml.org/dict/cml",
                                          "units":"http://www.xml-cml.org/units/units",
                                          "xsd":"http://www.w3c.org/2001/XMLSchema",
                                          "iupac":"http://www.iupac.org"})
    _write_atoms(cmpd, root, cmpd.particles())
    _write_bonds(cmpd, root, cmpd.bonds())

    tree = ET.ElementTree(root)
    tree.write(file_path, pretty_print=True, xml_declaration=True, encoding="utf-8")


def _write_atoms(self, root, atoms):
    atomArray = ET.SubElement(root, 'atomArray')
    for atom in atoms:
        atom_cml = ET.SubElement(atomArray, 'atom')
        atom_cml.set('id', str(id(atom)))
        atom_cml.set('elementType', atom.name)
        atom_cml.set('x3', str(round(atom.xyz[0][0],5)))
        atom_cml.set('y3', str(round(atom.xyz[0][1],5)))
        atom_cml.set('z3', str(round(atom.xyz[0][2],5)))


def _write_bonds(self, root, bonds):
    bondArray = ET.SubElement(root, 'bondArray')
    for bond in bonds:
        bond_cml = ET.SubElement(bondArray, 'bond')
        bond_cml.set('atomRefs2', (str(id(bond[0]))+' '+str(id(bond[1]))))
        bond_cml.set('order', 'None')
