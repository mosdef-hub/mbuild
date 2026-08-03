"""Utility to build forces for relaxing coarse bead paths without a chemical forcefield.

Paths are relaxed with a coarse Kremer-Grest model built from the geometry
the path was generated with: a WCA pair potential for excluded volume, a
FENE+WCA bond, and, only when the caller supplies a sampler, a tabulated angle
and/or dihedral potential.
"""

import logging
import math

import numpy as np

logger = logging.getLogger(__name__)

TWO_POW_1_6 = 2.0 ** (1.0 / 6.0)

# Kremer-Grest equilibrium bond length as a fraction of the WCA sigma.
KG_BOND_SIGMA_RATIO = 0.97


class PathForcefield:
    """Coarse Kremer-Grest forcefield parameters for a bead Path.

    Passed to a simulation in place of a chemical forcefield to override the
    geometry-derived defaults. Any value left None falls back to its geometric
    default: the mean bond length, radius = bond_length / 0.97, and no angle or
    dihedral term.

    Parameters
    ----------
    radius : float or None, default None
        WCA sigma (bead exclusion diameter, center-to-center).
    bond_length : float or None, default None
        Target bond length for the FENE+WCA bond.
    angles : mbuild.path.points.AnglesSampler or None, default None
        Optional bond-angle distribution for a tabulated angle potential.
    dihedrals : object or None, default None
        Optional dihedral distribution for a tabulated dihedral potential.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.
    per_type : dict or None, default None
        Reserved for future per-bead-type parameters, e.g.
        {bead_type: {"radius": ..., "bond_length": ...}}. Not yet used.
    """

    def __init__(
        self,
        radius=None,
        bond_length=None,
        angles=None,
        dihedrals=None,
        epsilon=1.0,
        per_type=None,
    ):
        self.radius = radius
        self.bond_length = bond_length
        self.angles = angles
        self.dihedrals = dihedrals
        self.epsilon = epsilon
        self.per_type = per_type


def generate_ff_from_path(
    path,
    snap,
    radius=None,
    bond_length=None,
    angles_sampler=None,
    dihedrals_sampler=None,
    epsilon=1.0,
):
    """Build a coarse forcefield for a path from its geometry.

    When not given explicitly, the bond length is the mean over the path's
    bonds and the bead radius follows the Kremer-Grest relation l0 = sigma * 0.97
    Angle and dihedral terms are added only when a sampler is supplied.

    Parameters
    ----------
    path : mbuild.path.build.Path
        Path providing the geometry the bond length is derived from.
    snap : gsd.hoomd.Frame
        Snapshot providing particle, bond, angle, and dihedral type names.
    radius : float or None, default None
        WCA sigma. If None, derived from bond_length via the Kremer-Grest
        relation.
    bond_length : float or None, default None
        Target bond length. If None, the mean over the path's bonds is used.
    angles_sampler : mbuild.path.points.AnglesSampler or None, default None
        Optional bond-angle distribution for a tabulated angle potential.
    dihedrals_sampler : object or None, default None
        Optional dihedral distribution for a tabulated dihedral potential.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.

    Returns
    -------
    list
        HOOMD force objects, as returned by build_path_ff.
    """
    if bond_length is None:
        bond_length = _mean_bond_length(path)
    if radius is None:
        radius = bond_length / KG_BOND_SIGMA_RATIO
    return build_path_ff(
        snap,
        radius=radius,
        bond_length=bond_length,
        angle_sampler=angles_sampler,
        dihedral_sampler=dihedrals_sampler,
        epsilon=epsilon,
    )


def build_path_ff(
    snap,
    radius,
    bond_length,
    angle_sampler=None,
    dihedral_sampler=None,
    epsilon=1.0,
    table_width=100,
):
    """Build coarse Kremer-Grest HOOMD forces for a bead path.

    Every bead, bond, angle and dihedral type in the
    snapshot receives the same corresponding parameters.

    Parameters
    ----------
    snap : gsd.hoomd.Frame
        Snapshot providing particle, bond, angle, and dihedral type names.
    radius : float
        WCA sigma (bead exclusion diameter, center-to-center).
    bond_length : float
        Target bond length, used to size the FENE+WCA bond.
    angle_sampler : mbuild.path.points.AnglesSampler or None, default None
        If given, a tabulated angle potential is Boltzmann-inverted from it.
    dihedral_sampler : object or None, default None
        If given, a tabulated dihedral potential is Boltzmann-inverted from it.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.
    table_width : int, default 100
        Resolution of the tabulated angle and dihedral potentials.

    Returns
    -------
    list
        HOOMD force objects: [WCA pairs, FENE+WCA bonds, (angle), (dihedral)].
    """
    import hoomd

    forces = []

    # WCA excluded volume between non-bonded beads. 1-2 pairs are excluded and
    # handled by the bond's own WCA term.
    r_cut = TWO_POW_1_6 * radius
    nlist = hoomd.md.nlist.Cell(buffer=0.2, exclusions=("bond",))
    wca = hoomd.md.pair.LJ(nlist=nlist, default_r_cut=r_cut, mode="shift")
    ptypes = list(set(snap.particles.types))
    for i, t1 in enumerate(ptypes):
        for t2 in ptypes[i:]:
            wca.params[(t1, t2)] = dict(sigma=radius, epsilon=epsilon)
    forces.append(wca)

    if snap.bonds.N > 0:
        sigma_bond = bond_length
        if bond_length > 1.33 * radius:
            logger.warning(
                "build_path_ff: bond_length (%.3g) exceeds ~1.33*radius (%.3g); "
                "FENE+WCA no longer guarantees non-crossing bonds.",
                bond_length,
                radius,
            )
        fene = hoomd.md.bond.FENEWCA()
        for bt in set(snap.bonds.types):
            fene.params[bt] = dict(
                k=30.0 * epsilon / sigma_bond**2,
                r0=1.5 * sigma_bond,
                epsilon=epsilon,
                sigma=sigma_bond,
                delta=0.0,
            )
        forces.append(fene)

    # Tabulated angle from a sampled bond-angle distribution.
    if angle_sampler is not None and snap.angles.N > 0:
        _, U, tau = angle_table_from_sampler(angle_sampler, table_width)
        angle = hoomd.md.angle.Table(width=len(U))
        for angle_type in set(snap.angles.types):
            angle.params[angle_type] = dict(U=U, tau=tau)
        forces.append(angle)

    # Tabulated dihedral from a sampled dihedral distribution.
    if dihedral_sampler is not None and snap.dihedrals.N > 0:
        _, U, tau = dihedral_table_from_sampler(dihedral_sampler, table_width)
        dihedral = hoomd.md.dihedral.Table(width=len(U))
        for dih_type in set(snap.dihedrals.types):
            dihedral.params[dih_type] = dict(U=U, tau=tau)
        forces.append(dihedral)

    return forces


def angle_table_from_sampler(sampler, n_bins=100, jacobian=True, n_samples=500_000):
    """Boltzmann-invert an AnglesSampler into a tabulated angle potential.
    The Boltzmann-inversion assumes kT = 1.0

    Parameters
    ----------
    sampler : mbuild.path.points.AnglesSampler
        Provides sampled bond angles (radians) via sampler.sample(size).
    n_bins : int, default 100
        Number of evenly spaced theta points from 0 to pi.
    jacobian : bool, default True
        Divide the sampled density by sin(theta) before inverting, so that an
        MD angle distribution reproduces the sampled one.
    n_samples : int, default 500000
        Number of angles drawn to estimate the distribution.

    Returns
    -------
    theta, U, tau : numpy.ndarray
        Grid (radians), potential energy, and torque (-dU/dtheta), each of
        length n_bins, ready for hoomd.md.angle.Table.
    """
    theta = np.linspace(0.0, math.pi, n_bins)
    samples = np.clip(
        np.asarray(sampler.sample(size=n_samples), dtype=float), 0.0, math.pi
    )
    density, edges = np.histogram(
        samples, bins=n_bins, range=(0.0, math.pi), density=True
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    prob = np.interp(theta, centers, density)

    prob = np.maximum(prob, prob.max() * 1e-8)
    if jacobian:
        prob = prob / np.clip(np.sin(theta), 1e-6, None)

    # Assuming kT=1.
    # TODO: Possibly add kT parameter for paramaterized Boltzmann sampling.
    U = -1 * np.log(prob)

    # Beyond the sampled range the density is noise, and dividing it by a small
    # sin(theta) can produce wells at theta = 0 or pi. Force U to rise
    # monotonically toward each pole outside the sample quantiles.
    q_lo, q_hi = np.quantile(samples, [0.005, 0.995])
    lo = int(np.searchsorted(theta, q_lo))
    hi = int(np.searchsorted(theta, q_hi))
    for i in range(lo - 1, -1, -1):
        U[i] = max(U[i], U[i + 1])
    for i in range(hi, n_bins):
        U[i] = max(U[i], U[i - 1])

    U -= U.min()
    tau = -np.gradient(U, theta)
    return theta, U, tau


def dihedral_table_from_sampler(sampler, n_bins=100, n_samples=500_000):
    """Boltzmann-invert a dihedral-angle sampler into a tabulated potential.
    The Boltzmann-inversion assumes kT = 1.0

    Parameters
    ----------
    sampler : object with a sample(size) method
        Provides sampled dihedral angles (radians) in [-pi, pi].
    n_bins : int, default 100
        Number of evenly spaced phi points from -pi to pi.
    n_samples : int, default 500000
        Number of angles drawn to estimate the distribution.

    Returns
    -------
    phi, U, tau : numpy.ndarray
        Grid (radians), potential energy, and torque (-dU/dphi), each of
        length n_bins, ready for hoomd.md.dihedral.Table.
    """
    phi = np.linspace(-math.pi, math.pi, n_bins)
    samples = np.clip(
        np.asarray(sampler.sample(size=n_samples), dtype=float), -math.pi, math.pi
    )
    density, edges = np.histogram(
        samples, bins=n_bins, range=(-math.pi, math.pi), density=True
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    prob = np.interp(phi, centers, density)

    prob = np.maximum(prob, prob.max() * 1e-8)
    # Assuming kT=1. TODO: Possibly add kT parameter for Boltzmann sampling.
    U = -1 * np.log(prob)
    U -= U.min()
    tau = -np.gradient(U, phi)
    return phi, U, tau


def _mean_bond_length(path):
    """Return the mean center-to-center distance over the path's bonds.

    Parameters
    ----------
    path : mbuild.path.build.Path
        Path whose coordinates and bond_graph define the bonds.

    Returns
    -------
    float
        Mean bond length over every edge in the bond graph.

    Raises
    ------
    ValueError
        If the path has no bonds, so no length can be derived. Pass radius
        and bond_length explicitly through build_path_ff instead.
    """
    coords = np.asarray(path.coordinates)
    edges = list(path.bond_graph.edges())
    if not edges:
        raise ValueError(
            "Cannot derive a bond length from a path with no bonds. "
            "Use build_path_ff with explicit radius and bond_length."
        )
    lengths = [float(np.linalg.norm(coords[u] - coords[v])) for u, v in edges]
    return float(np.mean(lengths))
