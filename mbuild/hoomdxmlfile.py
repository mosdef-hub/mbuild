from warnings import warn
import pdb

import numpy as np

from atom import Atom
from bond import Bond
from compound import Compound

def load_hoomdxml(filename, component=None):
    """Load a HOOMD XML file into a Compound.

    If no Compound is specified, this function will return a new Compound.

    Args:
        filename (str): Path to the HOOMD XML file to be read.
        component (Compound, optional): Optionally read the HOOMD XML data into an
            existing Compound.

    Returns:
        component (Compound): The Compound containing the HOOMD XML file's data.

    """
    if component is None:
        component = Compound()
    pass


def write_hoomdxml(component, box=None, filename='mbuild.xml'):
    """Output a Compound as a HOOMD XML file.

    Args:
        component (Compound): The Compound to be output.
        filename (str, optional): Path of the output file.

    """
    n_atoms = len([atom for atom in component.atoms() if atom.kind != "G"])
    n_bonds = len(list(component.bonds()))

    with open(filename, 'w') as xml_file:
        xml_file.write("""<?xml version="1.3" encoding="UTF-8"?>\n""")
        xml_file.write("""<hoomd_xml>\n""")
        xml_file.write("""<configuration time_step="0">\n""")

        box_lengths = np.array([0.0, 0.0, 0.0])
        if box is not None:
            warn("Using box dimensions to overwrite periodicity of Compound. "
                    "You should have a very good reason for doing this!.")
            box_lengths = box.lengths
        elif (box is None) and np.all(component.periodicity):
            # Logic is redundant but saves time calculating boundingbox for
            # larger systems.
            box_lengths = component.periodicty
        else:
            bounding_box = component.boundingbox()
            for i, dim in enumerate(component.periodicity):
                if float(dim) != 0.0:
                    box_lengths[i] = dim
                else:
                    box_lengths[i] = bounding_box.lengths[i]
        lx, ly, lz = box_lengths
        xml_file.write("""<box lx="{0}" ly="{1}" lz="{2}"/>\n""".format(lx, ly, lz))

        atom_types = list()
        id_to_idx = dict()
        atom_idx = 0
        xml_file.write("""<position num="{0}">\n""".format(n_atoms))
        for atom in component.atoms():
            if atom.kind != "G":
                atom_types.append(atom.kind)
                x, y, z = atom.pos
                xml_file.write("{0:8.5f} {1:8.5f} {2:8.5f}\n".format(x, y, z))

                id_to_idx[id(atom)] = atom_idx
                atom_idx += 1
        xml_file.write("</position>\n")

        xml_file.write("""<type num="{0}">\n""".format(n_atoms))
        for atom_type in atom_types:
            xml_file.write("{0}\n".format(atom_type))
        xml_file.write("</type>\n")

        xml_file.write("<bond>\n")
        for bond in component.bonds():
            atom1_id = id_to_idx[id(bond.atom1)]
            atom2_id = id_to_idx[id(bond.atom2)]
            xml_file.write("{0} {1} {2}\n".format(bond.kind, atom1_id, atom2_id))
        xml_file.write("</bond>\n")

        xml_file.write("</configuration>\n")
        xml_file.write("</hoomd_xml>\n")
