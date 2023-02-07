"""mBuild packing module: a wrapper for PACKMOL.

http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml
"""
import os
import sys
import tempfile
import warnings
from distutils.spawn import find_executable
from itertools import zip_longest
from subprocess import PIPE, Popen

import numpy as np

from mbuild import clone
from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.exceptions import MBuildError

__all__ = ["fill_box", "fill_region", "fill_sphere", "solvate"]

PACKMOL = find_executable("packmol")
PACKMOL_HEADER = """
tolerance {0:.16f}
filetype xyz
output {1}
seed {2}
sidemax {3}
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
    {8}
end structure
"""
PACKMOL_SPHERE = """
structure {0}
    number {1:d}
    inside sphere {2:.3f} {3:.3f} {4:.3f} {5:.3f}
    {6}
end structure
"""

PACKMOL_CONSTRAIN = """
constrain_rotation x 0. 0.
constrain_rotation y 0. 0.
constrain_rotation z 0. 0.
"""


def fill_box(
    compound,
    n_compounds=None,
    box=None,
    density=None,
    overlap=0.2,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    compound_ratio=None,
    aspect_ratio=None,
    fix_orientation=False,
    temp_file=None,
    update_port_locations=False,
):
    """Fill a box with an `mbuild.compound` or `Compound` s using PACKMOL.

    `fill_box` takes a single `Compound` or a list of `Compound` s and
    returns a `Compound` that has been filled to specification by PACKMOL.

    When filling a system, two arguments of `n_compounds` , `box` , and
    `density` must be specified.

    If `n_compounds` and `box` are not None, the specified number of
    compounds will be inserted into a box of the specified size.

    If `n_compounds` and `density` are not None, the corresponding box size
    will be calculated internally. In this case, `n_compounds` must be an int
    and not a list of int.

    If `box` and `density` are not None, the corresponding number of
    compounds will be calculated internally.

    For the cases in which `box` is not specified but generated internally,
    the default behavior is to calculate a cubic box. Optionally,
    `aspect_ratio` can be passed to generate a non-cubic box.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
        Compound or list of compounds to fill in box.
    n_compounds : int or list of int
        Number of compounds to be filled in box.
    box : mb.Box
        Box to be filled by compounds.
    density : float, units :math:`kg/m^3` , default=None
        Target density for the system in macroscale units. If not None, one of
        `n_compounds` or `box` , but not both, must be specified.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    sidemax : float, optional, default=100.0
        Needed to build an initial approximation of the molecule distribution in
        PACKMOL. All system coordinates must fit with in +/- sidemax, so
        increase sidemax accordingly to your final box size.
    edge : float, units nm, default=0.2
        Buffer at the edge of the box to not place molecules. This is necessary
        in some systems because PACKMOL does not account for periodic boundary
        conditions in its optimization.
    compound_ratio : list, default=None
        Ratio of number of each compound to be put in box. Only used in the case
        of `density` and `box` having been specified, `n_compounds` not
        specified, and more than one `compound` .
    aspect_ratio : list of float
        If a non-cubic box is desired, the ratio of box lengths in the x, y, and
        z directions.
    fix_orientation : bool or list of bools
        Specify that compounds should not be rotated when filling the box,
        default=False.
    temp_file : str, default=None
        File name to write PACKMOL raw output to.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds can be
        rotated, port orientation may be incorrect.

    Returns
    -------
    filled : mb.Compound
    """
    # check that the user has the PACKMOL binary on their PATH
    _check_packmol(PACKMOL)

    arg_count = 3 - [n_compounds, box, density].count(None)
    if arg_count != 2:
        raise ValueError(
            "Exactly 2 of `n_compounds`, `box`, and `density` must be "
            "specified. {} were given.".format(arg_count)
        )

    if box is not None:
        (box, my_mins, my_maxs) = _validate_box(box)
    if not isinstance(compound, (list, set)):
        compound = [compound]
    if n_compounds is not None and not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]
    if not isinstance(fix_orientation, (list, set)):
        fix_orientation = [fix_orientation] * len(compound)

    if compound is not None and n_compounds is not None:
        if len(compound) != len(n_compounds):
            raise ValueError(
                "`compound` and `n_compounds` must be of equal length."
            )

    if compound is not None:
        if len(compound) != len(fix_orientation):
            raise ValueError(
                "`compound`, `n_compounds`, and `fix_orientation` must be of "
                "equal length."
            )

    if density is not None:
        total_mass = _validate_mass(compound, n_compounds)
        if box is None and n_compounds is not None:
            # Conversion from (amu/(kg/m^3))**(1/3) to nm
            L = (total_mass / density) ** (1 / 3) * 1.1841763
            if aspect_ratio is None:
                (box, my_mins, my_maxs) = _validate_box(
                    Box(lengths=[L, L, L], angles=[90.0, 90.0, 90.0])
                )
            else:
                L *= np.prod(aspect_ratio) ** (-1 / 3)
                (box, my_mins, my_maxs) = _validate_box(
                    [val * L for val in aspect_ratio]
                )
        if n_compounds is None and box is not None:
            if len(compound) == 1:
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_compounds = [
                    int(
                        density
                        / total_mass
                        * np.prod(np.asarray(box.lengths))
                        * 0.60224
                    )
                ]
            else:
                if compound_ratio is None:
                    raise ValueError(
                        "Determing `n_compounds` from `density` and `box` for "
                        "systems with more than one compound type requires "
                        "`compound_ratio`"
                    )
                if len(compound) != len(compound_ratio):
                    raise ValueError(
                        "Length of `compound_ratio` must equal length of "
                        "`compound`"
                    )
                prototype_mass = 0
                for c, r in zip(compound, compound_ratio):
                    prototype_mass += r * c.mass
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_prototypes = int(
                    density
                    / prototype_mass
                    * np.prod(np.asarray(box.lengths))
                    * 0.60224
                )
                n_compounds = list()
                for c in compound_ratio:
                    n_compounds.append(int(n_prototypes * c))

    # Convert nm to angstroms for PACKMOL.
    # TODO how to handle box starts and ends
    # not really a thing a box would know?
    box_mins = np.asarray(my_mins) * 10
    box_maxs = np.asarray(my_maxs) * 10
    overlap *= 10

    # Apply 0.2nm edge buffer
    box_maxs = [a_max - (edge * 10) for a_max in box_maxs]
    box_mins = [a_min + (edge * 10) for a_min in box_mins]

    # Build the input file for each compound and call packmol.
    filled_xyz = _new_xyz_file()

    # create a list to contain the file handles for the compound temp files
    compound_xyz_list = list()
    try:
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10
        )
        for comp, m_compounds, rotate in zip(
            compound, n_compounds, fix_orientation
        ):
            m_compounds = int(m_compounds)

            compound_xyz = _new_xyz_file()
            compound_xyz_list.append(compound_xyz)

            comp.save(compound_xyz.name, overwrite=True)
            input_text += PACKMOL_BOX.format(
                compound_xyz.name,
                m_compounds,
                box_mins[0],
                box_mins[1],
                box_mins[2],
                box_maxs[0],
                box_maxs[1],
                box_maxs[2],
                PACKMOL_CONSTRAIN if rotate else "",
            )
        _run_packmol(input_text, filled_xyz, temp_file)
        # Create the topology and update the coordinates.
        filled = Compound()
        filled = _create_topology(filled, compound, n_compounds)
        filled.update_coordinates(
            filled_xyz.name, update_port_locations=update_port_locations
        )
        filled.box = box

    # ensure that the temporary files are removed from the machine after filling
    finally:
        for file_handle in compound_xyz_list:
            file_handle.close()
            os.unlink(file_handle.name)
        filled_xyz.close()
        os.unlink(filled_xyz.name)
    return filled


def fill_region(
    compound,
    n_compounds,
    region,
    overlap=0.2,
    bounds=None,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    fix_orientation=False,
    temp_file=None,
    update_port_locations=False,
):
    """Fill a region of a box with `mbuild.Compound` (s) using PACKMOL.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
        Compound or list of compounds to fill in region.
    n_compounds : int or list of ints
        Number of compounds to be put in region.
    region : mb.Box or list of mb.Box
        Region to be filled by compounds.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    sidemax : float, optional, default=100.0
        Needed to build an initial approximation of the molecule distribution in
        PACKMOL. All system coordinates must fit with in +/- sidemax, so
        increase sidemax accordingly to your final box size.
    edge : float, units nm, default=0.2
        Buffer at the edge of the region to not place molecules. This is
        necessary in some systems because PACKMOL does not account for periodic
        boundary conditions in its optimization.
    fix_orientation : bool or list of bools
        Specify that compounds should not be rotated when filling the box,
        default=False.
    bounds : list-like of floats [minx, miny, minz, maxx, maxy, maxz], units nm, default=None
        Bounding within box to pack compounds, if you want to pack within a bounding
        area that is not the full extent of the region, bounds are required.
    temp_file : str, default=None
        File name to write PACKMOL raw output to.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds can be
        rotated, port orientation may be incorrect.

    Returns
    -------
    filled : mb.Compound

    If using mulitple regions and compounds, the nth value in each list are used
    in order.
    For example, the third compound will be put in the third region using the
    third value in n_compounds.
    """
    # check that the user has the PACKMOL binary on their PATH
    _check_packmol(PACKMOL)
    if not isinstance(compound, (list, set)):
        compound = [compound]
    if not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]
    if not isinstance(fix_orientation, (list, set)):
        fix_orientation = [fix_orientation] * len(compound)

    if compound is not None and n_compounds is not None:
        if len(compound) != len(n_compounds):
            raise ValueError(
                "`compound` and `n_compounds` must be of equal length."
            )
    if compound is not None:
        if len(compound) != len(fix_orientation):
            raise ValueError(
                "`compound`, `n_compounds`, and `fix_orientation` must be of "
                "equal length."
            )

    # See if region is a single region or list
    my_regions = []
    if isinstance(region, Box):  # Cannot iterate over boxes
        my_regions.append(region)
    # if region is a list of boxes or a list of lists of floats append to my_regions, otherwise the list is expected to be a list of floats
    elif isinstance(region, list):
        for reg in region:
            if isinstance(reg, (list, Box)):
                my_regions.append(reg)
            else:
                raise ValueError(
                    f"list contents expected to be mbuild.Box or list of floats, provided: {type(reg)}"
                )
    else:
        raise ValueError(
            f"expected a list of type: list or mbuild.Box, was provided {region} of type: {type(region)}"
        )
    container = []
    if not bounds:
        bounds = []
    for bound, reg in zip_longest(bounds, my_regions, fillvalue=None):
        if bound is None:
            container.append(reg)
        else:
            container.append(bound)
    container = [_validate_box(bounding) for bounding in container]

    # In angstroms for packmol.
    overlap *= 10

    # Build the input file and call packmol.
    filled_xyz = _new_xyz_file()

    # List to hold file handles for the temporary compounds
    compound_xyz_list = list()
    try:
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10
        )
        for comp, m_compounds, rotate, items_n in zip(
            compound, n_compounds, fix_orientation, container
        ):
            m_compounds = int(m_compounds)

            compound_xyz = _new_xyz_file()
            compound_xyz_list.append(compound_xyz)

            comp.save(compound_xyz.name, overwrite=True)
            # TODO how to handle these mins and maxs of this system
            # box should not have any idea of mins and maxs
            my_min = items_n[1]
            my_max = items_n[2]
            reg_mins = np.asarray(my_min) * 10.0
            reg_maxs = np.asarray(my_max) * 10.0

            reg_maxs -= edge * 10  # Apply edge buffer
            input_text += PACKMOL_BOX.format(
                compound_xyz.name,
                m_compounds,
                reg_mins[0],
                reg_mins[1],
                reg_mins[2],
                reg_maxs[0],
                reg_maxs[1],
                reg_maxs[2],
                PACKMOL_CONSTRAIN if rotate else "",
            )

        _run_packmol(input_text, filled_xyz, temp_file)

        # Create the topology and update the coordinates.
        filled = Compound()
        filled = _create_topology(filled, compound, n_compounds)
        filled.update_coordinates(
            filled_xyz.name, update_port_locations=update_port_locations
        )
    finally:
        for file_handle in compound_xyz_list:
            file_handle.close()
            os.unlink(file_handle.name)
        filled_xyz.close()
        os.unlink(filled_xyz.name)
    return filled


def fill_sphere(
    compound,
    sphere,
    n_compounds=None,
    density=None,
    overlap=0.2,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    compound_ratio=None,
    fix_orientation=False,
    temp_file=None,
    update_port_locations=False,
):
    """Fill a sphere with a compound using PACKMOL.

    One argument of `n_compounds and density` must be specified.

    If `n_compounds` is not None, the specified number of n_compounds will be
    inserted into a sphere of the specified size.

    If `density` is not None, the corresponding number of compounds will be
    calculated internally.

    Parameters
    ----------
    compound : mb.Compound or list of mb.Compound
        Compound or list of compounds to be put in box.
    sphere : list, units nm
        Sphere coordinates in the form [x_center, y_center, z_center, radius]
    n_compounds : int or list of int
        Number of compounds to be put in box.
    density : float, units :math:`kg/m^3`, default=None
        Target density for the sphere in macroscale units.
    overlap : float, units nm, default=0.2
        Minimum separation between atoms of different molecules.
    seed : int, default=12345
        Random seed to be passed to PACKMOL.
    sidemax : float, optional, default=100.0
        Needed to build an initial approximation of the molecule distribution in
        PACKMOL. All system coordinates must fit with in +/- sidemax, so
        increase sidemax accordingly to your final sphere size
    edge : float, units nm, default=0.2
        Buffer at the edge of the sphere to not place molecules. This is
        necessary in some systems because PACKMOL does not account for periodic
        boundary conditions in its optimization.
    compound_ratio : list, default=None
        Ratio of number of each compound to be put in sphere. Only used in the
        case of `density` having been specified, `n_compounds` not specified,
        and more than one `compound`.
    fix_orientation : bool or list of bools
        Specify that compounds should not be rotated when filling the sphere,
        default=False.
    temp_file : str, default=None
        File name to write PACKMOL raw output to.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds can be
        rotated, port orientation may be incorrect.

    Returns
    -------
    filled : mb.Compound
    """
    _check_packmol(PACKMOL)

    arg_count = 2 - [n_compounds, density].count(None)
    if arg_count != 1:
        raise ValueError(
            "Exactly 1 of `n_compounds` and `density` must be specified. "
            "{} were given.".format(arg_count)
        )

    if isinstance(sphere, (list, set, tuple)):
        if len(sphere) != 4:
            raise ValueError("`sphere` must be a list of len 4")
    else:
        raise TypeError("`sphere` must be a list")

    if not isinstance(compound, (list, set)):
        compound = [compound]
    if n_compounds is not None and not isinstance(n_compounds, (list, set)):
        n_compounds = [n_compounds]
    if not isinstance(fix_orientation, (list, set)):
        fix_orientation = [fix_orientation] * len(compound)

    if compound is not None and n_compounds is not None:
        if len(compound) != len(n_compounds):
            raise ValueError(
                "`compound` and `n_compounds` must be of equal length."
            )

    if compound is not None:
        if len(compound) != len(fix_orientation):
            raise ValueError(
                "`compound`, `n_compounds`, and `fix_orientation` must be of "
                "equal length."
            )

    for coord in sphere[:3]:
        if coord < sphere[3]:
            raise ValueError(
                "`sphere` center coordinates must be greater than radius."
            )

    # Apply edge buffer
    radius = sphere[3] - edge

    if density is not None:
        total_mass = _validate_mass(compound, n_compounds)
        if n_compounds is None:
            if len(compound) == 1:
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_compounds = [
                    int(
                        density
                        / total_mass
                        * (4 / 3 * np.pi * radius**3)
                        * 0.60224
                    )
                ]
            else:
                if compound_ratio is None:
                    raise ValueError(
                        "Determing `n_compounds` from `density` for systems "
                        "with more than one compound type requires"
                        "`compound_ratio`"
                    )
                if len(compound) != len(compound_ratio):
                    raise ValueError(
                        "Length of `compound_ratio` must equal length of "
                        "`compound`"
                    )
                prototype_mass = 0
                for c, r in zip(compound, compound_ratio):
                    prototype_mass += r * c.mass
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_prototypes = int(
                    density
                    / prototype_mass
                    * (4 / 3 * np.pi * radius**3)
                    * 0.60224
                )
                n_compounds = list()
                for c in compound_ratio:
                    n_compounds.append(int(n_prototypes * c))

    # In angstroms for packmol.
    sphere = np.multiply(sphere, 10)
    radius *= 10
    overlap *= 10

    # Build the input file for each compound and call packmol.
    filled_xyz = _new_xyz_file()

    # List to hold file handles for the temporary compounds
    compound_xyz_list = list()
    try:
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10
        )
        for comp, m_compounds, rotate in zip(
            compound, n_compounds, fix_orientation
        ):
            m_compounds = int(m_compounds)

            compound_xyz = _new_xyz_file()
            compound_xyz_list.append(compound_xyz)

            comp.save(compound_xyz.name, overwrite=True)
            input_text += PACKMOL_SPHERE.format(
                compound_xyz.name,
                m_compounds,
                sphere[0],
                sphere[1],
                sphere[2],
                radius,
                PACKMOL_CONSTRAIN if rotate else "",
            )
        _run_packmol(input_text, filled_xyz, temp_file)

        # Create the topology and update the coordinates.
        filled = Compound()
        filled = _create_topology(filled, compound, n_compounds)
        filled.update_coordinates(
            filled_xyz.name, update_port_locations=update_port_locations
        )
    finally:
        for file_handle in compound_xyz_list:
            file_handle.close()
            os.unlink(file_handle.name)
        filled_xyz.close()
        os.unlink(filled_xyz.name)
    return filled


def solvate(
    solute,
    solvent,
    n_solvent,
    box,
    overlap=0.2,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    fix_orientation=False,
    temp_file=None,
    update_port_locations=False,
):
    """Solvate a compound in a box of solvent using PACKMOL.

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
    sidemax : float, optional, default=100.0
        Needed to build an initial approximation of the molecule distribution in
        PACKMOL. All system coordinates must fit with in +/- sidemax, so
        increase sidemax accordingly to your final box size
    edge : float, units nm, default=0.2
        Buffer at the edge of the box to not place molecules. This is necessary
        in some systems because PACKMOL does not account for periodic boundary
        conditions in its optimization.
    fix_orientation : bool
        Specify if solvent should not be rotated when filling box,
        default=False.
    temp_file : str, default=None
        File name to write PACKMOL raw output to.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds can be
        rotated, port orientation may be incorrect.

    Returns
    -------
    solvated : mb.Compound
    """
    # check that the user has the PACKMOL binary on their PATH
    _check_packmol(PACKMOL)

    (box, min_tmp, max_tmp) = _validate_box(box)
    if not isinstance(solvent, (list, set)):
        solvent = [solvent]
    if not isinstance(n_solvent, (list, set)):
        n_solvent = [n_solvent]
    if not isinstance(fix_orientation, (list, set)):
        fix_orientation = [fix_orientation] * len(solvent)

    if len(solvent) != len(n_solvent):
        raise ValueError("`n_solvent` and `n_solvent` must be of equal length.")

    # In angstroms for packmol.
    box_mins = np.asarray(min_tmp) * 10
    box_maxs = np.asarray(max_tmp) * 10
    overlap *= 10
    center_solute = (box_maxs + box_mins) / 2

    # Apply edge buffer
    box_maxs = np.subtract(box_maxs, edge * 10)
    box_mins = np.add(box_mins, edge * 10)
    # Build the input file for each compound and call packmol.
    solvated_xyz = _new_xyz_file()
    solute_xyz = _new_xyz_file()

    # generate list of temp files for the solvents
    solvent_xyz_list = list()
    try:
        solute.save(solute_xyz.name, overwrite=True)
        input_text = PACKMOL_HEADER.format(
            overlap, solvated_xyz.name, seed, sidemax * 10
        ) + PACKMOL_SOLUTE.format(solute_xyz.name, *center_solute.tolist())

        for solv, m_solvent, rotate in zip(solvent, n_solvent, fix_orientation):
            m_solvent = int(m_solvent)
            solvent_xyz = _new_xyz_file()
            solvent_xyz_list.append(solvent_xyz)

            solv.save(solvent_xyz.name, overwrite=True)
            input_text += PACKMOL_BOX.format(
                solvent_xyz.name,
                m_solvent,
                box_mins[0],
                box_mins[1],
                box_mins[2],
                box_maxs[0],
                box_maxs[1],
                box_maxs[2],
                PACKMOL_CONSTRAIN if rotate else "",
            )
        _run_packmol(input_text, solvated_xyz, temp_file)

        # Create the topology and update the coordinates.
        solvated = Compound()
        solvated.add(clone(solute))
        solvated = _create_topology(solvated, solvent, n_solvent)
        solvated.update_coordinates(
            solvated_xyz.name, update_port_locations=update_port_locations
        )

    finally:
        for file_handle in solvent_xyz_list:
            file_handle.close()
            os.unlink(file_handle.name)
        solvated_xyz.close()
        solute_xyz.close()
        os.unlink(solvated_xyz.name)
        os.unlink(solute_xyz.name)
    return solvated


def _validate_mass(compound, n_compounds):
    """Check the mass of the compounds passed into the packing functions.

    Raises an error if the total mass is zero, and density cannot be used to
    find box size or number of compounds.
    Returns a warning of any subcompound in compound has a mass of zero.
    """
    if n_compounds is None:
        n_compounds = [1] * len(compound)
    found_zero_mass = False
    total_mass = 0
    for c, n in zip(compound, n_compounds):
        comp_masses = np.array([c._particle_mass(p) for p in c.particles()])
        if 0.0 in comp_masses or None in comp_masses:
            found_zero_mass = True
        comp_masses[comp_masses == None] = 0.0
        total_mass += np.sum(comp_masses) * n

    if total_mass == 0:
        raise MBuildError(
            "The total mass of your compound(s) is zero "
            "In order to use density when packing a box, the mass of "
            "the compounds must be set. See the doc strings of the "
            "Compound() class in compound.py for more information "
            "on how mass is handled."
        )
    if found_zero_mass:
        warnings.warn(
            "Some of the compounds or subcompounds in `compound` "
            "have a mass of zero/None. This may have an effect on "
            "density calculations"
        )
    return total_mass


def _validate_box(box):
    """Ensure that the box passed by the user can be formatted as an mbuild.Box.

    Parameters
    ----------
    box : mbuild.Box or a tuple or list thereof
        Box or inputs to `mbuild.Box` to generate a `mbuild.Box`.

    Returns
    -------
    box : mbuild.Box
    mins : list-like
    maxs : list-like
    """
    if isinstance(box, (list, tuple)):
        if len(box) == 3:
            mins = [0.0, 0.0, 0.0]
            maxs = box
            box = Box.from_mins_maxs_angles(
                mins=mins, maxs=maxs, angles=(90.0, 90.0, 90.0)
            )
        elif len(box) == 6:
            mins = box[:3]
            maxs = box[3:]
            box = Box.from_mins_maxs_angles(
                mins=mins, maxs=maxs, angles=(90.0, 90.0, 90.0)
            )
        else:
            raise MBuildError(
                "Unknown format for `box` parameter. Must pass a"
                " list/tuple of length 3 (box lengths) or length"
                " 6 (box mins and maxes) or an mbuild.Box object."
            )

    elif isinstance(box, Box):
        mins = [0.0, 0.0, 0.0]
        maxs = box.lengths
    else:
        raise MBuildError(
            "Unknown format for `box` parameter. Must pass a list/tuple of "
            "length 3 (box lengths) or length 6 (box mins and maxes) or an "
            "mbuild.Box object."
        )
    return (box, mins, maxs)


def _new_xyz_file():
    """Generate PDB file using tempfile.NamedTemporaryFile.

    Return
    ------
    _ : file-object
        Temporary PDB file.
    """
    return tempfile.NamedTemporaryFile(suffix=".xyz", delete=False)


def _create_topology(container, comp_to_add, n_compounds):
    """Return updated mBuild compound with new coordinates.

    Parameters
    ----------
    container : mb.Compound, required
        Compound containing the updated system generated by PACKMOL.
    comp_to_add : mb.Compound or list of mb.Compounds, required
        Compound(s) to add to the container.
    container : int or list of int, required
        Amount of comp_to_add to container.

    Return
    ------
    container : mb.Compound
        Compound with added compounds from PACKMOL.
    """
    for comp, m_compound in zip(comp_to_add, n_compounds):
        for _ in range(m_compound):
            container.add(clone(comp))
    return container


def _packmol_error(out, err):
    """Log packmol output to files."""
    with open("log.txt", "w") as log_file:
        log_file.write(out)
    raise RuntimeError("PACKMOL failed. See 'log.txt'")


def _run_packmol(input_text, filled_xyz, temp_file):
    """Call PACKMOL to pack system based on the input text.

    Parameters
    ----------
    input_text : str, required
        String formatted in the input file syntax for PACKMOL.
    filled_xyz : `tempfile` object, required
        Tempfile that will store the results of PACKMOL packing.
    temp_file : `tempfile` object, required
        Where to copy the filled tempfile.
    """
    # Create input file
    packmol_inp = tempfile.NamedTemporaryFile(
        mode="w", delete=False, prefix="packmol-", suffix=".inp"
    )
    packmol_inp.write(input_text)
    packmol_inp.close()

    proc = Popen(
        "{} < {}".format(PACKMOL, packmol_inp.name),
        stdin=PIPE,
        stdout=PIPE,
        stderr=PIPE,
        universal_newlines=True,
        shell=True,
    )
    out, err = proc.communicate()

    if "WITHOUT PERFECT PACKING" in out:
        warnings.warn(
            "Packmol finished with imperfect packing. Using the .xyz_FORCED "
            "file instead. This may not be a sufficient packing result."
        )
        os.system("cp {0}_FORCED {0}".format(filled_xyz.name))

    if "ERROR" in out or proc.returncode != 0:
        _packmol_error(out, err)
    else:
        # Delete input file if success
        os.remove(packmol_inp.name)

    if temp_file is not None:
        os.system("cp {0} {1}".format(filled_xyz.name, os.path.join(temp_file)))


def _check_packmol(PACKMOL):  # pragma: no cover
    if not PACKMOL:
        msg = "Packmol not found."
        if sys.platform.startswith("win"):
            msg = (
                msg + " If packmol is already installed, make sure that the "
                "packmol.exe is on the path."
            )
        raise IOError(msg)
