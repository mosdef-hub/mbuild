"""GSD format.

https://gsd.readthedocs.io/en/stable/
"""

import gsd.hoomd

import mbuild as mb
from mbuild.conversion import to_hoomdsnapshot

__all__ = ["write_gsd"]


def write_gsd(
    compound,
    filename,
    identify_connections=True,
    ref_distance=1.0,
    ref_mass=1.0,
    rigid_bodies=None,
    shift_coords=True,
    write_special_pairs=True,
    **kwargs,
):
    """Output a GSD file (HOOMD-Blue topology format).

    Parameters
    ----------
    compound : mb.Compound
        mBuild compound to save to the GSD format.
    filename : str
        Path of the output file.
    identify_connections : bool
        If `True`, then infer angles and dihedrals in the topology.
    ref_distance : float, optional, default=1.0
        Reference distance for conversion to reduced units
    ref_mass : float, optional, default=1.0
        Reference mass for conversion to reduced units
    rigid_bodies : list of int, optional, default=None
        List of rigid body information. An integer value is required for each
        atom corresponding to the index of the rigid body the particle is to be
        associated with. A value of None indicates the atom is not part of a
        rigid body.
    shift_coords : bool, optional, default=True
        Shift coordinates from (0, L) to (-L/2, L/2) if necessary.
    write_special_pairs : bool, optional, default=True
        Writes out special pair information necessary to correctly use the OPLS
        fudged 1,4 interactions in HOOMD-Blue.

    Notes
    -----
    See gmso.external.convert_hoomd on how to make sure the GSD file contains
    forcefield information (e.g. atom types, angle types, etc.).

    """
    snapshot = to_hoomdsnapshot(
        compound=compound,
        identify_connections=identify_connections,
        ref_distance=ref_distance,
        ref_mass=ref_mass,
        rigid_bodies=rigid_bodies,
        shift_coords=shift_coords,
        write_special_pairs=write_special_pairs,
    )

    with gsd.hoomd.open(filename, "w") as traj:
        traj.append(snapshot)
