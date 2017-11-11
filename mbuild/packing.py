from __future__ import division

import sys
import os
import tempfile
import warnings
from distutils.spawn import find_executable
from subprocess import Popen, PIPE

import numpy as np

from mbuild.compound import Compound
from mbuild.exceptions import MBuildError
from mbuild.box import Box
from mbuild import clone

__all__ = ['fill_box', 'fill_region', 'solvate']

PACKMOL = find_executable('packmol')
PACKMOL_HEADER = """
tolerance {0:.16f}
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


def fill_box(compound, n_compounds=None, box=None, density=None, overlap=0.2,
             seed=12345, edge=0.2, compound_ratio=None,
             aspect_ratio=None, temp_file=None):
    """Fill a box with a compound using packmol.

    Two arguments of `n_compounds, box, and density` must be specified.

    If `n_compounds` and `box` are not None, the specified number of
    n_compounds will be inserted into a box of the specified size.

    If `n_compounds` and `density` are not None, the corresponding box
    size will be calculated internally. In this case, `n_compounds`
    must be an int and not a list of int.

    If `box` and `density` are not None, the corresponding number of
    compounds will be calculated internally.

    For the cases in which `box` is not specified but generated internally,
    the default behavior is to calculate a cubic box. Optionally,
    `aspect_ratio` can be passed to generate a non-cubic box.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
        Compound or list of compounds to be put in box.
    n_compounds : int or list of int
        Number of compounds to be put in box.
    box : mb.Box
        Box to be filled by compounds.
    density : float, units kg/m^3, default=None
        Target density for the system in macroscale units. If not None, one of
        `n_compounds` or `box`, but not both, must be specified.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    edge : float, units nm, default=0.2
        Buffer at the edge of the box to not place molecules. This is necessary
        in some systems because PACKMOL does not account for periodic boundary
        conditions in its optimization.
    compound_ratio : list, default=None
        Ratio of number of each compound to be put in box. Only used in the
        case of `density` and `box` having been specified, `n_compounds` not
        specified, and more than one `compound`.
    aspect_ratio : list of float
        If a non-cubic box is desired, the ratio of box lengths in the x, y,
        and z directions.
    temp_file : str, default=None
        File name to write PACKMOL's raw output to.

    Returns
    -------
    filled : mb.Compound
    
    """
    _check_packmol(PACKMOL)

    arg_count = 3 - [n_compounds, box, density].count(None)
    if arg_count != 2:
        msg = ("Exactly 2 of `n_compounds`, `box`, and `density` "
            "must be specified. {} were given.".format(arg_count))
        raise ValueError(msg)

    if box is not None:
        box = _validate_box(box)
    if not isinstance(compound, (list, set)):
        compound = [compound]
    if n_compounds is not None and not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]

    if compound is not None and n_compounds is not None:
        if len(compound) != len(n_compounds):
            msg = ("`compound` and `n_compounds` must be of equal length.")
            raise ValueError(msg)

    if density is not None:
        if box is None and n_compounds is not None:
            total_mass = np.sum([n*np.sum([a.mass for a in c.to_parmed().atoms])
                for c,n in zip(compound, n_compounds)])
            # Conversion from (amu/(kg/m^3))**(1/3) to nm
            L = (total_mass/density)**(1/3)*1.1841763
            if aspect_ratio is None:
                box = _validate_box(Box(3*[L]))
            else:
                L *= np.prod(aspect_ratio) ** (-1/3)
                box = _validate_box(Box([val*L for val in aspect_ratio]))
        if n_compounds is None and box is not None:
            if len(compound) == 1:
                compound_mass = np.sum([a.mass for a in compound[0].to_parmed().atoms])
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_compounds = [int(density/compound_mass*np.prod(box.lengths)*.60224)]
            else:
                if compound_ratio is None:
                    msg = ("Determing `n_compounds` from `density` and `box` "
                           "for systems with more than one compound type requires"
                           "`compound_ratio`")
                    raise ValueError(msg)
                if len(compound) != len(compound_ratio):
                    msg = ("Length of `compound_ratio` must equal length of "
                           "`compound`")
                    raise ValueError(msg)
                prototype_mass = 0
                for c, r in zip(compound, compound_ratio):
                    prototype_mass += r * np.sum([a.mass for a in c.to_parmed().atoms])
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_prototypes = int(density/prototype_mass*np.prod(box.lengths)*.60224)
                n_compounds = list()
                for c in compound_ratio:
                    n_compounds.append(int(n_prototypes * c))

    # In angstroms for packmol.
    box_mins = box.mins * 10
    box_maxs = box.maxs * 10
    overlap *= 10
    
    # Apply edge buffer 
    box_maxs -= edge * 10

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

    _run_packmol(input_text, filled_pdb, temp_file)

    # Create the topology and update the coordinates.
    filled = Compound()
    for comp, m_compounds in zip(compound, n_compounds):
        for _ in range(m_compounds):
            filled.add(clone(comp))
    filled.update_coordinates(filled_pdb)
    filled.periodicity = np.asarray(box.lengths, dtype=np.float32)
    return filled


def fill_region(compound, n_compounds, region, overlap=0.2,
                seed=12345, edge=0.2, temp_file=None):
    """Fill a region of a box with a compound using packmol.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
        Compound or list of compounds to be put in region.
    n_compounds : int or list of int
        Number of compounds to be put in region.
    region : mb.Box or list of mb.Box
        Region to be filled by compounds.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    edge : float, units nm, default=0.2
        Buffer at the edge of the region to not place molecules. This is
        necessary in some systems because PACKMOL does not account for
        periodic boundary conditions in its optimization.
    temp_file : str, default=None
        File name to write PACKMOL's raw output to.

    Returns
    -------
    filled : mb.Compound

    If using mulitple regions and compounds, the nth value in each list are used in order.
    For example, if the third compound will be put in the third region using the third value in n_compounds.
    """
    _check_packmol(PACKMOL)

    if not isinstance(compound, (list, set)):
        compound = [compound]
    if not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]

    if compound is not None and n_compounds is not None:
        if len(compound) != len(n_compounds):
            msg = ("`compound` and `n_compounds` must be of equal length.")
            raise ValueError(msg)

    # See if region is a single region or list
    if isinstance(region, Box): # Cannot iterate over boxes
        region = [region]
    elif not any(isinstance(reg, (list, set, Box)) for reg in region):
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
        reg_maxs -= edge * 10 # Apply edge buffer
        input_text += PACKMOL_BOX.format(compound_pdb, m_compounds,
                                        reg_mins[0], reg_mins[1], reg_mins[2],
                                        reg_maxs[0], reg_maxs[1], reg_maxs[2])

    _run_packmol(input_text, filled_pdb, temp_file)

    # Create the topology and update the coordinates.
    filled = Compound()
    for comp, m_compounds in zip(compound, n_compounds):
        for _ in range(m_compounds):
            filled.add(clone(comp))
    filled.update_coordinates(filled_pdb)
    return filled


def solvate(solute, solvent, n_solvent, box, overlap=0.2,
            seed=12345, edge=0.2, temp_file=None):
    """Solvate a compound in a box of solvent using packmol.

    Parameters
    ----------
    solute : mb.Compound
        Compound to be placed in a box and solvated.
    solvent : mb.Compound
        Compound to solvate the box.
    n_solvent : int
        Number of solvents to be put in box.
    box : mb.Box
        Box to be filled by compounds.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    edge : float, units nm, default=0.2
        Buffer at the edge of the box to not place molecules. This is necessary
        in some systems because PACKMOL does not account for periodic boundary
        conditions in its optimization.
    temp_file : str, default=None
        File name to write PACKMOL's raw output to.

    Returns
    -------
    solvated : mb.Compound

    """
    _check_packmol(PACKMOL)

    box = _validate_box(box)
    if not isinstance(solvent, (list, set)):
        solvent = [solvent]
    if not isinstance(n_solvent, (list, set)):
        n_solvent = [n_solvent]

    if len(solvent) != len(n_solvent):
        msg = ("`n_solvent` and `n_solvent` must be of equal length.")
        raise ValueError(msg)

    # In angstroms for packmol.
    box_mins = box.mins * 10
    box_maxs = box.maxs * 10
    overlap *= 10
    center_solute = (box_maxs + box_mins) / 2
    
    # Apply edge buffer 
    box_maxs -= edge * 10

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

    _run_packmol(input_text, solvated_pdb, temp_file)

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

def _run_packmol(input_text, filled_pdb, temp_file):
    proc = Popen(PACKMOL, stdin=PIPE, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate(input=input_text)

    if 'WITHOUT PERFECT PACKING' in out:
        msg = ("Packmol finished with imperfect packing. Using "
               "the .pdb_FORCED file instead. This may not be a "
               "sufficient packing result.")
        warnings.warn(msg)
        os.system('cp {0}_FORCED {0}'.format(filled_pdb))
    if 'ERROR' in out:
        _packmol_error(out, err)

    if temp_file is not None:
        os.system('cp {0} {1}'.format(filled_pdb, os.path.join(temp_file)))

def _check_packmol(PACKMOL):
    if not PACKMOL:
        msg = "Packmol not found."
        if sys.platform.startswith("win"):
            msg = (msg + " If packmol is already installed, make sure that the "
                         "packmol.exe is on the path.")
        raise IOError(msg)
