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

# Kremer-Grest equilibrium bond length as a fraction of the WCA sigma.
KG_BOND_SIGMA_RATIO = 0.97

# Used as default r_cut to give purely repulsive (WCA)
TWO_POW_1_6 = 2.0 ** (1.0 / 6.0)


class PathForcefield:
    """Coarse Kremer-Grest forcefield parameters for a bead Path.

    Passed to a simulation in place of a chemical forcefield to override the
    geometry-derived defaults. Any value left None falls back to its geometric
    default: per-type bond lengths sampled from the system, per-bead-type radii
    from those lengths, and no angle or dihedral term.

    Each parameter may be a single value applied to every type, or a dict for
    per-type parameterization keyed by the type names of the built system (its
    gmso/HOOMD snapshot types), for example radius={"A": 0.5, "B": 0.8} and
    bond_length={"A-B": 0.3}. radius keys are bead-type names; bond_length,
    angles and dihedrals keys are the bond/angle/dihedral-type names gmso
    assigns (sorted bead names joined by "-", e.g. "A-B", "A-B-A"). radius and
    bond_length dicts must cover every type; an angles/dihedrals dict need only
    cover the types it should act on.

    Parameters
    ----------
    radius : float or dict or None, default None
        WCA sigma (bead exclusion diameter, center-to-center), keyed by bead
        type. Unlike-pair sigma is the arithmetic mean of the two bead radii.
    bond_length : float or dict or None, default None
        Target bond length for the FENE+WCA bond, keyed by bond-type name.
    angles : mbuild.path.points.AnglesSampler or dict or None, default None
        Optional bond-angle distribution(s) for a tabulated angle potential,
        keyed by angle-type name.
    dihedrals : object or dict or None, default None
        Optional dihedral distribution(s) for a tabulated dihedral potential,
        keyed by dihedral-type name.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.
    r_cut : float, default 2**(1/6)
        Pair cutoff in units of sigma. The default gives purely repulsive
        WCA/KG pairs; larger values include the attractive portion of the LJ
        well, which helps close excluded-volume voids during relaxation.
    """

    def __init__(
        self,
        radius=None,
        bond_length=None,
        btype=None,
        angles=None,
        dihedrals=None,
        epsilon=1.0,
        r_cut=TWO_POW_1_6,
    ):
        self.radius = radius
        self.bond_length = bond_length
        self.btype = btype
        self.angles = angles
        self.dihedrals = dihedrals
        self.epsilon = epsilon
        self.r_cut = r_cut


def generate_ff_from_path(
    snap,
    radius=None,
    bond_length=None,
    btype="harmonic",
    angles_sampler=None,
    dihedrals_sampler=None,
    epsilon=1.0,
    r_cut=TWO_POW_1_6,
):
    """Build a coarse Kremer-Grest forcefield from a system snapshot.

    Resolves each parameter against the snapshot's type names and passes fully
    keyed dicts to build_path_ff. A scalar is broadcast to every type; None is
    auto-derived from the snapshot geometry (per-bond-type median lengths, and
    per-bead-type radii from the incident bond lengths via bond_l = sigma * 0.97);
    a dict is used as given, keyed by type name. Angle and dihedral terms are
    added only when a sampler is supplied.

    Parameters
    ----------
    snap : gsd.hoomd.Frame
        Snapshot providing the particle, bond, angle, and dihedral types and
        the geometry the defaults are derived from.
    radius : float or dict or None, default None
        WCA sigma, keyed by bead type. If None, derived from incident bonds.
    bond_length : float or dict or None, default None
        Target bond length, keyed by bond-type name. If None, per-bond-type
        medians over the snapshot's bonds are used.
    btype : str, optional, default='harmonic'
        Which bond type to use by default. Harmonic bonds are useful when the current
        bond_length distribution is far from the bond_length target.
    angles_sampler : mbuild.path.points.AnglesSampler or dict or None, default None
        Optional bond-angle distribution(s), keyed by angle-type name.
    dihedrals_sampler : object or dict or None, default None
        Optional dihedral distribution(s), keyed by dihedral-type name.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.
    r_cut : float, default 2**(1/6)
        Pair cutoff in units of sigma. The default gives purely repulsive
        WCA/KG pairs; larger values include the attractive portion of the well.

    Notes
    -----
    Using the default values gives a Kremer-Grest bead-spring model with a purely
    repulsive and shifted pair force (WCA), FENE bonds, and no angle or dihedral forces.
    To include attractive portions of the LJ force, use larger values of r_cut.

    Returns
    -------
    list
        HOOMD force objects, as returned by build_path_ff.
    """
    radius = _resolve_param(
        radius, list(snap.particles.types or []), "radius", lambda: _auto_radii(snap)
    )
    bond_length = _resolve_param(
        bond_length,
        list(snap.bonds.types or []),
        "bond_length",
        lambda: _auto_bond_lengths(snap),
    )
    angle_sampler = _resolve_sampler(angles_sampler, list(snap.angles.types or []))
    dihedral_sampler = _resolve_sampler(
        dihedrals_sampler, list(snap.dihedrals.types or [])
    )
    return build_path_ff(
        snap,
        radius=radius,
        bond_length=bond_length,
        btype=btype,
        angle_sampler=angle_sampler,
        dihedral_sampler=dihedral_sampler,
        epsilon=epsilon,
        r_cut=r_cut,
    )


def build_path_ff(
    snap,
    radius,
    bond_length,
    btype="harmonic",
    angle_sampler=None,
    dihedral_sampler=None,
    epsilon=1.0,
    table_width=100,
    r_cut=TWO_POW_1_6,
):
    """Build coarse Kremer-Grest HOOMD forces for a bead path.

    All parameters are dicts keyed by the snapshot's own type names. radius is
    keyed by bead (particle) type and bond_length by bond-type name, and both
    must cover every type present. angle_sampler/dihedral_sampler are None or
    dicts keyed by angle/dihedral-type name and may cover a subset (uncovered
    types get a flat, zero table). LJ cross-pairs sigma is the arithmetic
    mean of the two bead radii.

    Parameters
    ----------
    snap : gsd.hoomd.Frame
        Snapshot providing particle, bond, angle, and dihedral type names.
    radius : dict
        {bead_type: WCA sigma}, the center-to-center bead exclusion diameter.
    bond_length : dict
        {bond_type: target bond length}, used to size the Harmonic or FENE+WCA bond.
    btype : str, default="harmonic"
        The approach for sampling bonds. Can be either fene or harmonic
    angle_sampler : dict or None, default None
        {angle_type: sampler}; each is Boltzmann-inverted into a tabulated term.
    dihedral_sampler : dict or None, default None
        {dihedral_type: sampler}, Boltzmann-inverted into a tabulated term.
    epsilon : float, default 1.0
        Energy scale for the WCA and FENE potentials.
    table_width : int, default 100
        Resolution of the tabulated angle and dihedral potentials.
    r_cut : float, default 2**(1/6)
        Pair cutoff in units of sigma (per-pair cutoff is r_cut * sigma_ij).
        The default truncates at the LJ minimum for purely repulsive pairs
        (WCA/KG); larger values include the attractive portion of the LJ well.

    Returns
    -------
    list
        HOOMD force objects: [LJ pairs, FENE+WCA bonds, (angle), (dihedral)].
    """
    import hoomd

    forces = []
    ptypes = list(snap.particles.types)

    nlist = hoomd.md.nlist.Cell(buffer=0.2, exclusions=("bond",))
    lj = hoomd.md.pair.LJ(
        nlist=nlist, default_r_cut=r_cut * max(radius.values()), mode="shift"
    )
    for i, t1 in enumerate(ptypes):
        for t2 in ptypes[i:]:
            sigma_ij = 0.5 * (radius[t1] + radius[t2])
            lj.params[(t1, t2)] = dict(sigma=sigma_ij, epsilon=epsilon)
            lj.r_cut[(t1, t2)] = r_cut * sigma_ij
    forces.append(lj)

    if snap.bonds.N > 0:
        ptypeid = np.asarray(snap.particles.typeid)
        group = np.asarray(snap.bonds.group)
        btypeid = np.asarray(snap.bonds.typeid)
        if btype is None or btype.lower() == "harmonic":
            bond_template = hoomd.md.bond.Harmonic()
            for i, bt in enumerate(snap.bonds.types):
                sigma_bond = bond_length[bt]
                a, b = group[np.argmax(btypeid == i)]
                energy_units = epsilon / sigma_bond**2
                bond_template.params[bt] = dict(
                    k=100.0 * energy_units,
                    r0=sigma_bond,
                )
        elif btype.lower() == "fene":
            bond_template = hoomd.md.bond.FENEWCA()
            for i, bt in enumerate(snap.bonds.types):
                sigma_bond = bond_length[bt]
                a, b = group[np.argmax(btypeid == i)]
                sigma_pair = 0.5 * (
                    radius[ptypes[ptypeid[a]]] + radius[ptypes[ptypeid[b]]]
                )
                if sigma_bond > 1.33 * sigma_pair:
                    logger.warning(
                        "build_path_ff: bond_length (%.3g) for %s exceeds "
                        "~1.33*radius (%.3g); FENE+WCA no longer guarantees "
                        "non-crossing bonds.",
                        sigma_bond,
                        bt,
                        sigma_pair,
                    )
                bond_template.params[bt] = dict(
                    k=30.0 * epsilon / sigma_bond**2,
                    r0=1.5 * sigma_bond,
                    epsilon=epsilon,
                    sigma=sigma_bond,
                    delta=0.0,
                )
        else:
            raise ValueError(
                f"Incorrect argument {btype=}. Please use one of 'fene' or 'harmonic'."
            )

        forces.append(bond_template)

    # Tabulated angle(s). Types the sampler dict does not cover get a flat
    # (zero) table so HOOMD sees every angle type specified.
    if angle_sampler and snap.angles.N > 0:
        angle = hoomd.md.angle.Table(width=table_width)
        for at in snap.angles.types:
            if at in angle_sampler:
                _, U, tau = angle_table_from_sampler(angle_sampler[at], table_width)
            else:
                U = tau = np.zeros(table_width)
            angle.params[at] = dict(U=U, tau=tau)
        forces.append(angle)

    if dihedral_sampler and snap.dihedrals.N > 0:
        dihedral = hoomd.md.dihedral.Table(width=table_width)
        for dt in snap.dihedrals.types:
            if dt in dihedral_sampler:
                _, U, tau = dihedral_table_from_sampler(
                    dihedral_sampler[dt], table_width
                )
            else:
                U = tau = np.zeros(table_width)
            dihedral.params[dt] = dict(U=U, tau=tau)
        forces.append(dihedral)

    return forces


def _resolve_param(value, type_names, label, auto_fn):
    """Resolve radius/bond_length to a {type_name: float} dict over every type.

    None auto-derives from the snapshot; a scalar broadcasts to all types; a
    dict is used as given and must cover every type name.
    """
    if value is None:
        return auto_fn()
    if isinstance(value, dict):
        missing = [t for t in type_names if t not in value]
        if missing:
            raise KeyError(f"{label} has no entry for type(s) {missing}.")
        return {t: float(value[t]) for t in type_names}
    return {t: float(value) for t in type_names}


def _resolve_sampler(sampler, type_names):
    """Resolve a sampler (or dict of them) to a {type_name: sampler} dict, or None.

    A single sampler broadcasts to every type; a dict is used as given (keyed
    by type name) and may cover only the types it should act on.
    """
    if sampler is None:
        return None
    if isinstance(sampler, dict):
        return {t: sampler[t] for t in type_names if t in sampler} or None
    return {t: sampler for t in type_names}


def _bond_lengths(snap):
    """Per-bond center-to-center lengths from the snapshot, minimum-imaged.

    Assumes an orthorhombic box: only the Lx/Ly/Lz box lengths are used and any
    triclinic tilt factors are ignored (paths produce orthorhombic boxes). The
    minimum image also assumes bonds shorter than half the box.
    """
    pos = np.asarray(snap.particles.position)
    group = np.asarray(snap.bonds.group)
    d = pos[group[:, 0]] - pos[group[:, 1]]
    box = np.asarray(snap.configuration.box[:3], dtype=float)
    periodic = box > 0
    d[:, periodic] -= box[periodic] * np.round(d[:, periodic] / box[periodic])
    return np.linalg.norm(d, axis=1)


def _auto_bond_lengths(snap):
    """{bond_type: median length} from the snapshot geometry.

    The median (not the mean) is robust to bonds that wrap across a periodic
    boundary whose true period is not the snapshot box, which appear as a few
    large outliers and would otherwise inflate the sampled length.
    """
    if snap.bonds.N == 0:
        return {}
    lengths = _bond_lengths(snap)
    btypeid = np.asarray(snap.bonds.typeid)
    return {
        bt: float(np.median(lengths[btypeid == i]))
        for i, bt in enumerate(snap.bonds.types)
    }


def _auto_radii(snap):
    """{bead_type: radius} from median incident bond length via l0 = sigma * 0.97.

    Uses the median for robustness to periodic-wrap bonds (see
    _auto_bond_lengths). Bead types with no incident bond fall back to the
    overall median length.
    """
    if snap.bonds.N == 0:
        raise ValueError(
            "Cannot derive radii from a system with no bonds. "
            "Pass radius (and bond_length) explicitly."
        )
    lengths = _bond_lengths(snap)
    group = np.asarray(snap.bonds.group)
    ptypeid = np.asarray(snap.particles.typeid)
    ptypes = list(snap.particles.types)
    incident = {t: [] for t in ptypes}
    for (a, b), length in zip(group, lengths):
        incident[ptypes[ptypeid[a]]].append(length)
        incident[ptypes[ptypeid[b]]].append(length)
    fallback = float(np.median(lengths))
    return {
        t: (float(np.median(v)) if v else fallback) / KG_BOND_SIGMA_RATIO
        for t, v in incident.items()
    }


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

    # U(theta) = -kTlog(P(theta)) --> Assume kT = 1
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
