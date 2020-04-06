import os
import sys
import tempfile
from warnings import warning
import numpy as np

import mdtraj as md
from mdtraj.core.element import get_by_symbol
import parmed as pmd
from parmed.periodic_table import AtomicNum, element_by mass, Mass, Element

import mbuild as mb
from mbuild.compound import Compound
from mbuild.formats.xyz import read_xyz
from mbuild.formats.json_formats import compound_from_json

def load(filename_or_object, relative_to_module=None, compound=None,
         coords_only=False, rigid=False, smiles=False, infer_hierarchy=True,
         backend=None, **kwargs):
    """Load a file or an existing topology into an mbuild compound.

    Files are read using the MDTraj package unless the `use_parmed` argument is
    specified as True. Please refer to http://mdtraj.org/1.8.0/load_functions.html
    for formats supported by MDTraj and https://parmed.github.io/ParmEd/html/
    readwrite.html for formats supported by ParmEd.

    Parameters
    ----------
    filename_or_object : str, mdtraj.Trajectory, parmed.Structure, mbuild.Compound,
            pybel.Molecule
        Name of the file or topology from which to load atom and bond information.
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory, look for the file
        where this module is defined. This is typically used in Compound
        classes that will be instantiated from a different directory
        (such as the Compounds located in mbuild.lib).
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information into.
        New structure will be added to the existing compound as a sub compound.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    backend : str, optional, default=None
        Backend used to load structure from file. If not specified, a default
        backend (extension specific) will be used.
    smiles: bool, optional, default=False
        Use Open Babel to parse filename as a SMILES string
        or file containing a SMILES string.
    infer_hierarchy : bool, optional, default=True
        If True, infer hierarchy from chains and residues
    **kwargs : keyword arguments
        Key word arguments passed to mdTraj for loading.

    Returns
    -------
    compound : mb.Compound

    """
    # First check if we are loading from an object)
    if not isinstance(filename_or_object, str):
        return load_object()
    # Second check if we are loading SMILES strings
    elif smiles:
        return load_smiles()
    # Last, if none of the above, load from file
    else:
        return load_file()


def load_object():
    """ Helper function to load an object into a mb.Compound

    Functions to load on-disk object to mb.Compound, supporting conversion
    from mb.Compound, pmd.Structure, md.Trajectory, and pybel.Molecule. If the
    compound flag is supplied, this will add the object as a subcompound of that
    compound.

    Parameters
    ----------

    Return
    ------

    """

    return None

def load_smiles(smiles_or_filename, compound):
    """ Helper function to load a SMILES string

    Loading SMILES string from a string, a list, or a file using pybel.
    Must have pybel packages installed.

    Parameters
    ----------

    Return
    ------

    """

    return None

def load_file():
    """ Helper function to load from files

    Loading and converting a topology to mb.Compound from file. User can specify
    prefered backend, or else it will be handle by default prefered backend.

    Parameters
    ----------

    Return
    ------


def from_mbuild():
    """ Backend-specific loading function - mbuild

    """

def from_parmed():
    """ Backend-specific loading function - parmed

    """

def from_trajectory():
    """ Backend-specific loading function - mdtraj

    """

def from_pybel():
    """ Backend-specific loading function - pybel

    """
