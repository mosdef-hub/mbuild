from __future__ import division

import sys
import tempfile
from distutils.spawn import find_executable
from subprocess import Popen, PIPE

from mbuild.compound import Compound
from mbuild.exceptions import MBuildError
from mbuild.box import Box
from mbuild import clone

__all__ = ['fill_box', 'fill_region', 'solvate']

PACKMOL = find_executable('packmol')
PACKMOL_HEADER = """
tolerance {0:.1f}
filetype pdb
output {1}
seed {2}

"""
PACKMOL_SOLUTE = """
structure {0}
    number 1
    center
    fixed {1:.3f} {2:.3f} {3:.3f} 0. 0. 0.
end structure
"""
PACKMOL_BOX = """
structure {0}
    number {1:d}
    inside box {2:.3f} {3:.3f} {4:.3f} {5:.3f} {6:.3f} {7:.3f}
end structure
"""


def fill_box(compound, n_compounds, box, overlap=0.2, seed=12345):
    """Fill a box with a compound using packmol.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
    n_compounds : int or list of int
    box : mb.Box
    overlap : float

    Returns
    -------
    filled : mb.Compound

    """
    if not PACKMOL:
        msg = "Packmol not found."
        if sys.platform.startswith("win"):
            msg = (msg + " If packmol is already installed, make sure that the "
                         "packmol.exe is on the path.")
        raise IOError(msg)

    box = _validate_box(box)
    if not isinstance(compound, (list, set)):
        compound = [compound]
    if not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]

    # In angstroms for packmol.
    box_mins = box.mins * 10
    box_maxs = box.maxs * 10
    overlap *= 10

    # Build the input file for each compound and call packmol.
    filled_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    input_text = PACKMOL_HEADER.format(overlap, filled_pdb, seed)

    for comp, m_compounds in zip(compound, n_compounds):
        m_compounds = int(m_compounds)
        compound_pdb = tempfile.mkstemp(suffix='.pdb')[1]
        comp.save(compound_pdb, overwrite=True)
        input_text += PACKMOL_BOX.format(compound_pdb, m_compounds,
                           box_mins[0], box_mins[1], box_mins[2],
                           box_maxs[0], box_maxs[1], box_maxs[2])

    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)
    if err:
        _packmol_error(out, err)

    # Create the topology and update the coordinates.
    filled = Compound()
    for comp, m_compounds in zip(compound, n_compounds):
        for _ in range(m_compounds):
            filled.add(clone(comp))
    filled.update_coordinates(filled_pdb)
    return filled


def fill_region(compound, n_compounds, region, overlap=0.2, seed=12345):
    """Fill a region of a box with a compound using packmol.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
    n_compounds : int or list of int
    region : mb.Box or list of mb.Box
    overlap : float

    Returns
    -------
    filled : mb.Compound

    If using mulitple regions and compounds, the nth value in each list are used in order.
    For example, if the third compound will be put in the third region using the third value in n_compounds.
    """
    if not PACKMOL:
        msg = "Packmol not found."
        if sys.platform.startswith("win"):
            msg = (msg + " If packmol is already installed, make sure that the "
                         "packmol.exe is on the path.")
        raise IOError(msg)

    if not isinstance(compound, (list, set)):
        compound = [compound]
    if not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]
    if not any(isinstance(reg, list) for reg in region): # See if region is a single region or list
        region = [region]
    region = [_validate_box(reg) for reg in region]

    # In angstroms for packmol.
    overlap *= 10

    # Build the input file and call packmol.
    filled_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    input_text = PACKMOL_HEADER.format(overlap, filled_pdb, seed)

    for comp, m_compounds, reg in zip(compound, n_compounds, region):
        m_compounds = int(m_compounds)
        compound_pdb = tempfile.mkstemp(suffix='.pdb')[1]
        comp.save(compound_pdb, overwrite=True)
        reg_mins = reg.mins * 10
        reg_maxs = reg.maxs * 10
        input_text += PACKMOL_BOX.format(compound_pdb, m_compounds,
                                        reg_mins[0], reg_mins[1], reg_mins[2],
                                        reg_maxs[0], reg_maxs[1], reg_maxs[2])

    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)
    if err:
        _packmol_error(out, err)

    # Create the topology and update the coordinates.
    filled = Compound()
    for comp, m_compounds in zip(compound, n_compounds):
        for _ in range(m_compounds):
            filled.add(clone(comp))
    filled.update_coordinates(filled_pdb)
    return filled


def solvate(solute, solvent, n_solvent, box, overlap=0.2, seed=12345):
    """Solvate a compound in a box of solvent using packmol.

    Parameters
    ----------
    solute : mb.Compound
    solvent : mb.Compound
    n_solvent : int
    box : mb.Box
    overlap : float

    Returns
    -------
    solvated : mb.Compound

    """
    if not PACKMOL:
        raise IOError("Packmol not found")

    box = _validate_box(box)
    if not isinstance(solvent, (list, set)):
        solvent = [solvent]
    if not isinstance(n_solvent, (list, set)):
        n_solvent = [n_solvent]

    # In angstroms for packmol.
    box_mins = box.mins * 10
    box_maxs = box.maxs * 10
    overlap *= 10
    center_solute = (box_maxs + box_mins) / 2

    # Build the input file for each compound and call packmol.
    solvated_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    solute_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    solute.save(solute_pdb, overwrite=True)
    input_text = (PACKMOL_HEADER.format(overlap, solvated_pdb, seed) +
                  PACKMOL_SOLUTE.format(solute_pdb, *center_solute))

    for solv, m_solvent in zip(solvent, n_solvent):
        m_solvent = int(m_solvent)
        solvent_pdb = tempfile.mkstemp(suffix='.pdb')[1]
        solv.save(solvent_pdb, overwrite=True)
        input_text += PACKMOL_BOX.format(solvent_pdb, m_solvent,
                           box_mins[0], box_mins[1], box_mins[2],
                           box_maxs[0], box_maxs[1], box_maxs[2])

    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)
    if err:
        _packmol_error(out, err)

    # Create the topology and update the coordinates.
    solvated = Compound()
    solvated.add(solute)
    for solv, m_solvent in zip(solvent, n_solvent):
        for _ in range(m_solvent):
            solvated.add(clone(solv))
    solvated.update_coordinates(solvated_pdb)
    return solvated


def _validate_box(box):
    if isinstance(box, (list, tuple)):
        if len(box) == 3:
            box = Box(lengths=box)
        elif len(box) == 6:
            box = Box(mins=box[:3], maxs=box[3:])

    if not isinstance(box, Box):
        raise MBuildError('Unknown format for `box` parameter. Must pass a'
                          ' list/tuple of length 3 (box lengths) or length'
                          ' 6 (box mins and maxes) or an mbuild.Box object.')
    return box


def _packmol_error(out, err):
    """Log packmol output to files. """
    with open('log.txt', 'w') as log_file, open('err.txt', 'w') as err_file:
        log_file.write(out)
        err_file.write(err)
    raise RuntimeError("PACKMOL failed. See 'err.txt' and 'log.txt'")
