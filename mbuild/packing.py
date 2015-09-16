from __future__ import division

import sys
import tempfile
from distutils.spawn import find_executable
from subprocess import Popen, PIPE

from mbuild.compound import Compound
from mbuild.box import Box
from mbuild import clone

__all__ = ['fill_box', 'solvate']

PACKMOL = find_executable('packmol')
PACKMOL_HEADER = """
tolerance {0:.1f}
filetype pdb
output {1}

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
    inside box 0. 0. 0. {2:.3f} {3:.3f} {4:.3f}
end structure
"""


def fill_box(compound, n_compounds, box, overlap=0.2):
    """Fill a box with a compound using packmol.

    Parameters
    ----------
    compound : mb.Compound
    n_compounds : int
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

    if isinstance(box, (list, tuple)):
        box = Box(lengths=box)

    compound_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    compound.save(compound_pdb)
    filled_pdb = tempfile.mkstemp(suffix='.pdb')[1]

    # In angstroms for packmol.
    box_lengths = box.lengths * 10
    overlap *= 10

    # Build the input file and call packmol.
    input_text = (PACKMOL_HEADER.format(overlap, filled_pdb) +
                  PACKMOL_BOX.format(compound_pdb, n_compounds, *box_lengths))

    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)
    if err:
        _packmol_error(out, err)

    # Create the topology and update the coordinates.
    filled = Compound()
    for _ in range(n_compounds):
        filled.add(clone(compound))
    filled.update_coordinates(filled_pdb)
    return filled


def solvate(solute, solvent, n_solvent, box, overlap=0.2):
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

    if isinstance(box, (list, tuple)):
        box = Box(lengths=box)

    solute_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    solute.save(solute_pdb)
    solvent_pdb = tempfile.mkstemp(suffix='.pdb')[1]
    solvent.save(solvent_pdb)
    solvated_pdb = tempfile.mkstemp(suffix='.pdb')[1]

    # In angstroms for packmol.
    box_lengths = box.lengths * 10
    overlap *= 10
    # center_solute = (-solute.center) * 10
    center_solute = box_lengths/2

    # Build the input file and call packmol.
    input_text = (PACKMOL_HEADER.format(overlap, solvated_pdb) +
                  PACKMOL_SOLUTE.format(solute_pdb, *center_solute) +
                  PACKMOL_BOX.format(solvent_pdb, n_solvent, *box_lengths))

    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)
    if err:
        _packmol_error(out, err)

    # Create the topology and update the coordinates.
    solvated = Compound()
    solvated.add(solute)
    for _ in range(n_solvent):
        solvated.add(clone(solvent))
    solvated.update_coordinates(solvated_pdb)
    return solvated


def _packmol_error(out, err):
    """Log packmol output to files. """
    with open('log.txt', 'w') as log_file, open('err.txt', 'w') as err_file:
        log_file.write(out)
        err_file.write(err)
    raise RuntimeError("PACKMOL failed. See 'err.txt' and 'log.txt'")
