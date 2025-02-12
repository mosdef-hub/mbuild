"""mBuild packing module: a wrapper for PACKMOL.

http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml
"""

import os
import shutil
import sys
import tempfile
import warnings
from itertools import zip_longest
from subprocess import PIPE, Popen

import numpy as np

from mbuild import clone
from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.exceptions import MBuildError

__all__ = ["fill_box", "fill_region", "fill_sphere", "solvate"]

PACKMOL = shutil.which("packmol")
PACKMOL_HEADER = """
tolerance {0:.16f}
filetype xyz
output {1}
seed {2}
sidemax {3}
{4}
{5}
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
    {2}
    {3}
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


def check_packmol_args(custom_args):
    # List of all available packmol_inputs.
    # Only file-level arguments can be passed.
    allowed_args = [
        "maxit",  # int
        "nloop",  # int
        "fbins",  # float
        "discale",  # float
        "movefrac",  # float
        "avoid_overlap",  # On/Off (empty string "" is on)
        "precision",  # float
        "movebadrandom",  # On/off (empty string "" is on)
        "use_short_tol",  # On/off (empty string "" is on)
        "short_tol_dist",  # float
        "short_tol_scale",  # float
        "tolerance",
        "seed",
        "sidemax",
    ]
    default_args = ["tolerance", "seed", "sidemax"]
    for key in custom_args:
        if key not in allowed_args:
            raise ValueError(
                f"PACKMOL argument {key} is not usable in `packmol_args`. "
                f"Availble arguments that can be set are {allowed_args}."
                "Only file-level arguments can be set with `packmol_args`."
                "See https://m3g.github.io/packmol/userguide.shtml#run"
            )
        if key in default_args:
            warnings.warn(
                f"The PACKMOL argument {key} was passed to `packmol_args`, "
                "but should be set using the corresponding function parameters. "
                "The value passed to the function will be used. "
                "See the function's doc strings for more information."
            )


def fill_box(
    compound,
    n_compounds=None,
    box=None,
    use_pbc=False,
    density=None,
    overlap=0.2,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    compound_ratio=None,
    aspect_ratio=None,
    fix_orientation=False,
    temp_file=None,
    save_packmol_input=False,
    update_port_locations=False,
    packmol_args=None,
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
    use_pbc : bool, default=False
        If True, applies periodic boundary conditions
        when placing molecules.
        `edge` must be zero when use_pbc is True.
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
    save_packmol_input : bool, default=False,
        If true, saves the file `packmol.inp` to the
        current working directory.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds
        can be rotated, port orientation may be incorrect.
    packmol_args : dict
        Dictionary where the key, value pairs are options and their
        corresponding keyword arguments for PACKMOL. Some PACKMOL options
        do not require a specified keyword. In this case, the value in
        the dictionary should be an empty string e.g. {'movebadrandom':""}
        These commands are placed at the header of the PACKMOL input file
        and therefore applied to all structures. NOTE: The PACKMOL options
        for seed and tolerance are specified by the function parameters
        seed and overlap.
        Other command options can be found in the PACKMOL userguide:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

    Notes
    -----
    The packmol_args parameter is designed to only accept file-level
    PACKMOL arguments, as opposed to molecule level arguments.

    The allowed arguments include the following:

        "maxit",  # int
        "nloop",  # int
        "fbins",  # float
        "discale",  # float
        "movefrac",  # float
        "avoid_overlap",  # On/Off (empty string "" is on)
        "precision",  # float
        "movebadrandom",  # On/off (empty string "" is on)
        "use_short_tol",  # On/off (empty string "" is on)
        "short_tol_dist",  # float
        "short_tol_scale",  # float

    The following PACKMOL arguments should be specified
    by the mbuild.packing method's parameters:

    "tolerance", "seed", "sidemax"

    See the PACKMOL documentation for detailed explanations of these arguments:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

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
            f"specified. {arg_count} were given."
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
            raise ValueError("`compound` and `n_compounds` must be of equal length.")

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
                        "Length of `compound_ratio` must equal length of `compound`"
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
    if use_pbc and edge != 0:
        raise ValueError("edge must be 0 if use_pbc is set to True")
    # Apply edge buffer
    box_maxs = [a_max - (edge * 10) for a_max in box_maxs]
    box_mins = [a_min + (edge * 10) for a_min in box_mins]
    box_arg = box_mins + box_maxs

    # generate string of addl. packmol inputs given in packmol_args
    packmol_commands = ""
    if packmol_args:
        check_packmol_args(packmol_args)
        for arg, val in packmol_args.items():
            packmol_commands += f"{arg} {val} \n"

    # Build the input file for each compound and call packmol.
    filled_xyz = _new_xyz_file()

    # create a list to contain the file handles for the compound temp files
    compound_xyz_list = list()
    try:
        if use_pbc:
            pbc_arg = "pbc {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}".format(
                *box_arg
            )
            fill_arg = ""
            periodicity = (True, True, True)
        else:
            fill_arg = (
                "inside box {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}".format(
                    *box_arg
                )
            )
            pbc_arg = ""
            periodicity = (False, False, False)
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10, packmol_commands, pbc_arg
        )
        for comp, m_compounds, rotate in zip(compound, n_compounds, fix_orientation):
            m_compounds = int(m_compounds)

            compound_xyz = _new_xyz_file()
            compound_xyz_list.append(compound_xyz)

            comp.save(compound_xyz.name, overwrite=True)
            input_text += PACKMOL_BOX.format(
                compound_xyz.name,
                m_compounds,
                fill_arg,
                PACKMOL_CONSTRAIN if rotate else "",
            )
        _run_packmol(input_text, filled_xyz, temp_file, save_packmol_input)
        # Create the topology and update the coordinates.
        filled = Compound(periodicity=periodicity)
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
    save_packmol_input=False,
    update_port_locations=False,
    packmol_args=None,
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
    bounds : list-like of list-likes of floats [[min_x, min_y, min_z, max_x, max_y, max_z], ...], units nm, default=None
        Bounding(s) within box to pack compound(s). To pack within a bounding
        area that is not the full extent of the region, bounds are required.
        Each item of `compound` must have its own bound specified. Use `None`
        to indicate a given compound is not bounded, e.g.
        `[ [0., 0., 1., 2., 2., 2.], None]` to bound only the first element
        of `compound` and not the second.
    temp_file : str, default=None
        File name to write PACKMOL raw output to.
    save_packmol_input : bool, default=False,
        If true, saves the file `packmol.inp` to the
        current working directory.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds
        can be rotated, port orientation may be incorrect.
    packmol_args : dict
        Dictionary where the key, value pairs are options and their
        corresponding keyword arguments for PACKMOL. Some PACKMOL options
        do not require a specified keyword. In this case, the value in
        the dictionary should be an empty string e.g. {'movebadrandom':""}
        These commands are placed at the header of the PACKMOL input file
        and therefore applied to all structures. NOTE: The PACKMOL options
        for seed and tolerance are specified by the function parameters
        seed and overlap.
        Other command options can be found in the PACKMOL userguide:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

    Notes
    -----
    The packmol_args parameter is designed to only accept file-level
    PACKMOL arguments, as opposed to molecule level arguments.

    The allowed arguments include the following:

        "maxit",  # int
        "nloop",  # int
        "fbins",  # float
        "discale",  # float
        "movefrac",  # float
        "avoid_overlap",  # On/Off (empty string "" is on)
        "precision",  # float
        "movebadrandom",  # On/off (empty string "" is on)
        "use_short_tol",  # On/off (empty string "" is on)
        "short_tol_dist",  # float
        "short_tol_scale",  # float

    The following PACKMOL arguments should be specified
    by the mbuild.packing method's parameters:

    "tolerance", "seed", "sidemax"

    See the PACKMOL documentation for detailed explanations of these arguments:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

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
            raise ValueError("`compound` and `n_compounds` must be of equal length.")
    if compound is not None:
        if len(compound) != len(fix_orientation):
            raise ValueError(
                "`compound`, `n_compounds`, and `fix_orientation` must be of "
                "equal length."
            )
    if bounds is not None:
        if not isinstance(bounds, (list)):
            raise TypeError(
                "`bounds` must be a list of one or more bounding boxes "
                "and/or `None` for each item in `compound`."
            )
        if len(bounds) != len(n_compounds):
            raise ValueError(
                "`compound` and `bounds` must be of equal length. Use `None` "
                "for non-bounded items in `compound`."
            )
        for bound in bounds:
            if not isinstance(bound, (Box, list)):
                raise ValueError(
                    "Each bound in `bounds` must be `None`, `Box`, or a "
                    "list of [min_x, min_y, min_z, max_x, max_y, max_z]."
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
    if bounds is None:
        bounds = []
    for bound, reg in zip_longest(bounds, my_regions, fillvalue=None):
        if bound is None:
            container.append(reg)
        else:
            container.append(bound)
    container = [_validate_box(bounding) for bounding in container]

    # In angstroms for packmol.
    overlap *= 10

    # generate string of addl. packmol inputs given in packmol_args
    packmol_commands = ""
    if packmol_args:
        check_packmol_args(packmol_args)
        for arg, val in packmol_args.items():
            packmol_commands += f"{arg} {val} \n"

    # Build the input file and call packmol.
    filled_xyz = _new_xyz_file()

    # List to hold file handles for the temporary compounds
    compound_xyz_list = list()
    try:
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10, packmol_commands, ""
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
            box_arg = list(reg_mins) + list(reg_maxs)
            fill_arg = (
                "inside box {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}".format(
                    *box_arg
                )
            )
            input_text += PACKMOL_BOX.format(
                compound_xyz.name,
                m_compounds,
                fill_arg,
                PACKMOL_CONSTRAIN if rotate else "",
            )

        _run_packmol(input_text, filled_xyz, temp_file, save_packmol_input)

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
    save_packmol_input=False,
    update_port_locations=False,
    packmol_args=None,
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
    save_packmol_input : bool, default=False,
        If true, saves the file `packmol.inp` to the
        current working directory.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds
        can be rotated, port orientation may be incorrect.
    packmol_args : dict
        Dictionary where the key, value pairs are options and their
        corresponding keyword arguments for PACKMOL. Some PACKMOL options
        do not require a specified keyword. In this case, the value in
        the dictionary should be an empty string e.g. {'movebadrandom':""}
        These commands are placed at the header of the PACKMOL input file
        and therefore applied to all structures. NOTE: The PACKMOL options
        for seed and tolerance are specified by the function parameters
        seed and overlap.
        Other command options can be found in the PACKMOL userguide:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

    Notes
    -----
    The packmol_args parameter is designed to only accept file-level
    PACKMOL arguments, as opposed to molecule level arguments.

    The allowed arguments include the following:

        "maxit",  # int
        "nloop",  # int
        "fbins",  # float
        "discale",  # float
        "movefrac",  # float
        "avoid_overlap",  # On/Off (empty string "" is on)
        "precision",  # float
        "movebadrandom",  # On/off (empty string "" is on)
        "use_short_tol",  # On/off (empty string "" is on)
        "short_tol_dist",  # float
        "short_tol_scale",  # float

    The following PACKMOL arguments should be specified
    by the mbuild.packing method's parameters:

    "tolerance", "seed", "sidemax"

    See the PACKMOL documentation for detailed explanations of these arguments:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

    Returns
    -------
    filled : mb.Compound
    """
    _check_packmol(PACKMOL)

    arg_count = 2 - [n_compounds, density].count(None)
    if arg_count != 1:
        raise ValueError(
            "Exactly 1 of `n_compounds` and `density` must be specified. "
            f"{arg_count} were given."
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
            raise ValueError("`compound` and `n_compounds` must be of equal length.")

    if compound is not None:
        if len(compound) != len(fix_orientation):
            raise ValueError(
                "`compound`, `n_compounds`, and `fix_orientation` must be of "
                "equal length."
            )

    for coord in sphere[:3]:
        if coord < sphere[3]:
            raise ValueError("`sphere` center coordinates must be greater than radius.")

    # Apply edge buffer
    radius = sphere[3] - edge

    if density is not None:
        total_mass = _validate_mass(compound, n_compounds)
        if n_compounds is None:
            if len(compound) == 1:
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_compounds = [
                    int(density / total_mass * (4 / 3 * np.pi * radius**3) * 0.60224)
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
                        "Length of `compound_ratio` must equal length of `compound`"
                    )
                prototype_mass = 0
                for c, r in zip(compound, compound_ratio):
                    prototype_mass += r * c.mass
                # Conversion from kg/m^3 / amu * nm^3 to dimensionless units
                n_prototypes = int(
                    density / prototype_mass * (4 / 3 * np.pi * radius**3) * 0.60224
                )
                n_compounds = list()
                for c in compound_ratio:
                    n_compounds.append(int(n_prototypes * c))

    # In angstroms for packmol.
    sphere = np.multiply(sphere, 10)
    radius *= 10
    overlap *= 10

    # generate string of addl. packmol inputs given in packmol_args
    packmol_commands = ""
    if packmol_args:
        check_packmol_args(packmol_args)
        for arg, val in packmol_args.items():
            packmol_commands += f"{arg} {val} \n"

    # Build the input file for each compound and call packmol.
    filled_xyz = _new_xyz_file()

    # List to hold file handles for the temporary compounds
    compound_xyz_list = list()
    try:
        input_text = PACKMOL_HEADER.format(
            overlap, filled_xyz.name, seed, sidemax * 10, packmol_commands, ""
        )
        for comp, m_compounds, rotate in zip(compound, n_compounds, fix_orientation):
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
        _run_packmol(input_text, filled_xyz, temp_file, save_packmol_input)

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
    use_pbc=False,
    overlap=0.2,
    seed=12345,
    sidemax=100.0,
    edge=0.2,
    fix_orientation=False,
    temp_file=None,
    save_packmol_input=False,
    update_port_locations=False,
    center_solute=True,
    packmol_args=None,
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
    use_pbc : bool, default=False
        If True, applies periodic boundary conditions
        when placing molecules.
        `edge` must be zero when use_pbc is True.
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
    save_packmol_input : bool, default=False,
        If true, saves the file `packmol.inp` to the
        current working directory.
    update_port_locations : bool, default=False
        After packing, port locations can be updated, but since compounds
        can be rotated, port orientation may be incorrect.
    center_solute : bool, optional, default=True
        Move solute center of mass to the center of the `mb.Box` used.
    packmol_args : dict
        Dictionary where the key, value pairs are options and their
        corresponding keyword arguments for PACKMOL. Some PACKMOL options
        do not require a specified keyword. In this case, the value in
        the dictionary should be an empty string e.g. {'movebadrandom':""}
        These commands are placed at the header of the PACKMOL input file
        and therefore applied to all structures. NOTE: The PACKMOL options
        for seed and tolerance are specified by the function parameters
        seed and overlap.
        Other command options can be found in the PACKMOL userguide:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

    Notes
    -----
    The packmol_args parameter is designed to only accept file-level
    PACKMOL arguments, as opposed to molecule level arguments.

    The allowed arguments include the following:

        "maxit",  # int
        "nloop",  # int
        "fbins",  # float
        "discale",  # float
        "movefrac",  # float
        "avoid_overlap",  # On/Off (empty string "" is on)
        "precision",  # float
        "movebadrandom",  # On/off (empty string "" is on)
        "use_short_tol",  # On/off (empty string "" is on)
        "short_tol_dist",  # float
        "short_tol_scale",  # float

    The following PACKMOL arguments should be specified
    by the mbuild.packing method's parameters:

    "tolerance", "seed", "sidemax"

    See the PACKMOL documentation for detailed explanations of these arguments:
        http://www.ime.unicamp.br/~martinez/packmol/userguide.shtml

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
    if use_pbc and edge != 0:
        raise ValueError("edge must be 0 if use_pbc is set to True")
    if center_solute:
        center_solute = (box_maxs + box_mins) / 2
    else:
        center_solute = solute.pos * 10
    # Apply edge buffer
    box_maxs = np.subtract(box_maxs, edge * 10)
    box_mins = np.add(box_mins, edge * 10)
    box_arg = list(box_mins) + list(box_maxs)
    # generate string of addl. packmol inputs given in packmol_args
    packmol_commands = ""
    if packmol_args:
        check_packmol_args(packmol_args)
        for arg, val in packmol_args.items():
            packmol_commands += f"{arg} {val} \n"

    # Build the input file for each compound and call packmol.
    solvated_xyz = _new_xyz_file()
    solute_xyz = _new_xyz_file()

    # generate list of temp files for the solvents
    solvent_xyz_list = list()
    try:
        if use_pbc:
            pbc_arg = "pbc {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}".format(
                *box_arg
            )
            fill_arg = ""
            periodicity = (True, True, True)
        else:
            fill_arg = (
                "inside box {0:.3f} {1:.3f} {2:.3f} {3:.3f} {4:.3f} {5:.3f}".format(
                    *box_arg
                )
            )
            pbc_arg = ""
            periodicity = (False, False, False)
        solute.save(solute_xyz.name, overwrite=True)
        input_text = PACKMOL_HEADER.format(
            overlap, solvated_xyz.name, seed, sidemax * 10, packmol_commands, pbc_arg
        ) + PACKMOL_SOLUTE.format(solute_xyz.name, *center_solute.tolist())

        for solv, m_solvent, rotate in zip(solvent, n_solvent, fix_orientation):
            m_solvent = int(m_solvent)
            solvent_xyz = _new_xyz_file()
            solvent_xyz_list.append(solvent_xyz)

            solv.save(solvent_xyz.name, overwrite=True)
            input_text += PACKMOL_BOX.format(
                solvent_xyz.name,
                m_solvent,
                fill_arg,
                PACKMOL_CONSTRAIN if rotate else "",
            )
        _run_packmol(input_text, solvated_xyz, temp_file, save_packmol_input)

        # Create the topology and update the coordinates.
        solvated = Compound(periodicity=periodicity)
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
        comp_masses[comp_masses is None] = 0.0
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
    container_list = []
    for comp, m_compound in zip(comp_to_add, n_compounds):
        for _ in range(m_compound):
            container_list.append(clone(comp))

    container.add(container_list)
    return container


def _packmol_error(out, err):
    """Log packmol output to files."""
    with open("log.txt", "w") as log_file:
        log_file.write(out)
    raise RuntimeError("PACKMOL failed. See 'log.txt'")


def _run_packmol(input_text, filled_xyz, temp_file, save_packmol_input):
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
    # Save PACKMOL file to cwd
    if save_packmol_input:
        new_file_path = "packmol.inp"
        with (
            open(packmol_inp.name, "r") as inp_file,
            open(new_file_path, "w") as new_file,
        ):
            shutil.copyfileobj(inp_file, new_file)

    proc = Popen(
        f"{PACKMOL} < {packmol_inp.name}",
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
        os.system(f"cp {filled_xyz.name}_FORCED {filled_xyz.name}")

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
