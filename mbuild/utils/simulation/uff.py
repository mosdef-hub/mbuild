"""Rule-based Universal Force Field (UFF) parameter generation for mBuild.

Implements UFF (Rappe et al., *J. Am. Chem. Soc.* 1992, 114, 10024) as an
on-the-fly generator: bonded and nonbonded parameters are computed directly
from element + bond order. Hybridization and aromaticity are read from the
bond-graph bond orders, not perceived from geometry or ring topology.

References
----------
Rappe, Casewit, Colwell, Goddard, Skiff. *J. Am. Chem. Soc.* 1992, 114, 10024.
Parameter values follow the canonical UFF table as distributed with Open Babel.
"""

import logging
import math

logger = logging.getLogger(__name__)

KCAL_TO_KJ = 4.184
ANG_TO_NM = 0.1
# kcal/mol/Ang^2 -> kJ/mol/nm^2  == 4.184 / 0.01
BOND_K_TO_KJ_NM2 = KCAL_TO_KJ / (ANG_TO_NM**2)

# UFF universal constants
LAMBDA = 0.1332  # bond-order correction coefficient (Rappe eq. 3)
G_BOND = 664.12  # bond force-constant prefactor, kcal/mol (== 2 * 332.06)
G_ANGLE = 664.12  # angle force-constant prefactor, kcal/mol

# UFF atom-type parameters (organic / main-group subset)
# Columns: r1 (valence bond radius, Ang), theta0 (equilibrium angle, deg),
#          x1 (vdW minimum distance r_min, Ang), D1 (vdW well depth, kcal/mol),
#          Zstar (effective charge), Vi (sp3 torsional barrier, kcal/mol),
#          Uj (sp2 torsional constant), chi (GMP electronegativity).
_COLS = ("r1", "theta0", "x1", "D1", "Zstar", "Vi", "Uj", "chi")
_TABLE = {
    #  label      r1     theta0   x1      D1     Zstar   Vi     Uj    chi
    "H_":     (0.354, 180.00, 2.886, 0.044, 0.712, 0.000, 0.0, 4.528),
    "C_3":    (0.757, 109.47, 3.851, 0.105, 1.912, 2.119, 2.0, 5.343),
    "C_R":    (0.729, 120.00, 3.851, 0.105, 1.912, 0.000, 2.0, 5.343),
    "C_2":    (0.732, 120.00, 3.851, 0.105, 1.912, 0.000, 2.0, 5.343),
    "C_1":    (0.706, 180.00, 3.851, 0.105, 1.912, 0.000, 2.0, 5.343),
    "N_3":    (0.700, 106.70, 3.660, 0.069, 2.544, 0.450, 2.0, 6.899),
    "N_R":    (0.699, 120.00, 3.660, 0.069, 2.544, 0.000, 2.0, 6.899),
    "N_2":    (0.685, 111.20, 3.660, 0.069, 2.544, 0.000, 2.0, 6.899),
    "N_1":    (0.656, 180.00, 3.660, 0.069, 2.544, 0.000, 2.0, 6.899),
    "O_3":    (0.658, 104.51, 3.500, 0.060, 2.300, 0.018, 2.0, 8.741),
    "O_R":    (0.680, 110.00, 3.500, 0.060, 2.300, 0.000, 2.0, 8.741),
    "O_2":    (0.634, 120.00, 3.500, 0.060, 2.300, 0.000, 2.0, 8.741),
    "O_1":    (0.639, 180.00, 3.500, 0.060, 2.300, 0.000, 2.0, 8.741),
    "F_":     (0.668, 180.00, 3.364, 0.050, 1.735, 0.000, 2.0, 10.874),
    "Cl":     (1.044, 180.00, 3.947, 0.227, 2.348, 0.000, 2.0, 8.564),
    "Br":     (1.192, 180.00, 4.189, 0.251, 2.519, 0.000, 2.0, 7.790),
    "I_":     (1.382, 180.00, 4.500, 0.339, 2.650, 0.000, 2.0, 6.822),
    "S_3+2":  (1.064, 92.10, 4.035, 0.274, 2.703, 0.484, 2.0, 6.928),
    "P_3+3":  (1.101, 93.80, 4.147, 0.305, 2.863, 2.400, 2.0, 5.463),
}
UFF_PARAMS = {label: dict(zip(_COLS, vals)) for label, vals in _TABLE.items()}

# Which UFF label to use per element, keyed by hybridization category derived
# from bond order.  Elements with a single supported type ignore the category.
_HYB_TYPES = {
    "H": {"*": "H_"},
    "C": {"sp3": "C_3", "sp2": "C_2", "ar": "C_R", "sp": "C_1"},
    "N": {"sp3": "N_3", "sp2": "N_2", "ar": "N_R", "sp": "N_1"},
    "O": {"sp3": "O_3", "sp2": "O_2", "ar": "O_R", "sp": "O_1"},
    "F": {"*": "F_"},
    "Cl": {"*": "Cl"},
    "Br": {"*": "Br"},
    "I": {"*": "I_"},
    "S": {"*": "S_3+2"},
    "P": {"*": "P_3+3"},
}


def uff_forces(top, snap, order_map, r_cut=1.0, backend="hoomd"):
    """Generate UFF forces as backend-native force objects.

    Requires the topology to already be UFF-typed and the snapshot built from
    that typed topology, so their type strings stay consistent with the forces.

    Parameters
    ----------
    top : gmso.Topology
        UFF-typed topology with connections identified.
    snap : gsd.hoomd.Frame
        Snapshot built from the UFF-typed topology (types are UFF labels).
    order_map : dict
        ``(i, j) -> bond order`` map, as returned by :func:`assign_uff_types`.
    r_cut : float, default 1.0
        Nonbonded cutoff (nm).
    backend : {"hoomd", "openmm"}, default "hoomd"
        Force representation to emit.

    Returns
    -------
    list
        Backend-native force objects.
    """
    if backend == "hoomd":
        return _hoomd_forces(top, snap, order_map, r_cut)

    if backend == "openmm":
        # snap is HOOMD-specific and unused here; OpenMM forces reference atom
        # indices (site order), which the caller's System particles must match.
        return _openmm_forces(top, order_map, r_cut)

    raise ValueError(f"Unknown backend '{backend}'. Options: 'hoomd', 'openmm'.")


class UFFTypingError(Exception):
    """Raised when a compound cannot be typed with the lean UFF subset."""


def assign_uff_types(top, compound):
    """Type every site in a GMSO topology with a UFF atom-type label.

    Sets ``site.name`` to the UFF label (e.g. ``"C_3"``) so that a GSD snapshot
    built afterward carries UFF-consistent particle / bond / angle type strings.
    Returns the per-index label list and an ``(i, j) -> order`` bond-order map.

    Hybridization and aromaticity come from the bond-graph bond orders, so
    aromatic bonds must be marked as such in the input to be typed as resonant.

    Parameters
    ----------
    top : gmso.Topology
        Topology whose sites are relabeled in place. Site order must match the
        order of ``compound.particles()``.
    compound : mb.Compound
        Source of elements and bond orders (via ``compound.bond_graph``).
    """
    particles = list(compound.particles())
    index_of = {p: i for i, p in enumerate(particles)}
    n = len(particles)

    # Numeric bond orders keyed by sorted index pair, plus per-atom order lists.
    order_map = {}
    atom_orders = {i: [] for i in range(n)}
    assumed_single = False
    for u, v, order in compound.bond_graph.edges.data("bond_order"):
        num, assumed = _num_order(order)
        assumed_single = assumed_single or assumed
        i, j = index_of[u], index_of[v]
        order_map[(min(i, j), max(i, j))] = num
        atom_orders[i].append(num)
        atom_orders[j].append(num)

    if assumed_single:
        logger.warning(
            "UFF: one or more bonds had missing or unrecognized bond order; "
            "assumed single bonds. Set bond_order for correct sp/sp2/sp3 and "
            "resonance typing."
        )

    labels = []
    for i, p in enumerate(particles):
        element = getattr(p, "element", None)
        if element is None or getattr(element, "symbol", None) is None:
            raise UFFTypingError(
                f"Particle {i} ('{p.name}') has no element set; UFF requires "
                "elements. Build the compound with element information."
            )
        labels.append(_label_for(element.symbol, atom_orders[i]))

    for site, label in zip(top.sites, labels):
        site.name = label

    return labels, order_map


_ORDER_STRINGS = {"single": 1.0, "double": 2.0, "triple": 3.0, "aromatic": 1.5}


def _num_order(order):
    """Map an mBuild bond-order value to a numeric UFF bond order.

    Returns (order, assumed), where assumed is True when the value was missing
    or unrecognized and a single bond was substituted.
    """
    if order in (None, "unspecified", "default", 0.0):
        return 1.0, True
    if isinstance(order, str):
        if order in _ORDER_STRINGS:
            return _ORDER_STRINGS[order], False
        return 1.0, True
    return float(order), False


def _hyb_category(orders):
    """Classify hybridization from the set of bond orders on an atom."""
    if any(abs(o - 1.5) < 1e-6 for o in orders):
        return "ar"
    max_order = max(orders) if orders else 1.0
    if max_order >= 2.5:
        return "sp"
    if max_order >= 1.5:
        return "sp2"
    return "sp3"


def _label_for(element_symbol, orders):
    """Return the UFF atom-type label for an element + its bond orders."""
    types = _HYB_TYPES.get(element_symbol)
    if types is None:
        raise UFFTypingError(
            f"Element '{element_symbol}' is not in the lean UFF subset. "
            f"Supported: {sorted(_HYB_TYPES)}."
        )
    if "*" in types:
        return types["*"]
    return types[_hyb_category(orders)]


def _bond_rest_length(ti, tj, order):
    """UFF natural bond length r_ij (Angstrom) for types ti, tj at bond order."""
    pi, pj = UFF_PARAMS[ti], UFF_PARAMS[tj]
    ri, rj = pi["r1"], pj["r1"]
    # Bond-order correction (Rappe eq. 3).
    r_bo = -LAMBDA * (ri + rj) * math.log(order)
    # Electronegativity correction (Rappe eq. 4).
    chi_i, chi_j = pi["chi"], pj["chi"]
    r_en = ri * rj * (math.sqrt(chi_i) - math.sqrt(chi_j)) ** 2 / (
        chi_i * ri + chi_j * rj
    )
    return ri + rj + r_bo - r_en


def _bond_force_constant(ti, tj, r_ij):
    """UFF bond force constant (kcal/mol/Ang^2) for E = 0.5 k (r - r_ij)^2."""
    return G_BOND * UFF_PARAMS[ti]["Zstar"] * UFF_PARAMS[tj]["Zstar"] / r_ij**3


def _angle_force_constant(ti, tj, tk, r_ij, r_jk, theta0):
    """UFF angle force constant (kcal/mol/rad^2), == curvature at theta0.

    ``theta0`` is in radians. Equal to the second derivative of the UFF
    cosine-Fourier angle term at equilibrium, i.e. the harmonic stiffness.
    """
    c = math.cos(theta0)
    r_ik2 = r_ij**2 + r_jk**2 - 2.0 * r_ij * r_jk * c  # law of cosines
    zi, zk = UFF_PARAMS[ti]["Zstar"], UFF_PARAMS[tk]["Zstar"]
    return (
        G_ANGLE
        * (zi * zk / r_ik2**2.5)
        * r_ij
        * r_jk
        * (3.0 * r_ij * r_jk * (1.0 - c * c) - r_ik2 * c)
    )


def _hyb_of(label):
    """Return 'sp3' | 'sp2' | 'sp' | None from a UFF label suffix."""
    if label.endswith("_3") or label.endswith("+2") or label.endswith("+3"):
        return "sp3"
    if label.endswith("_2") or label.endswith("_R"):
        return "sp2"
    if label.endswith("_1"):
        return "sp"
    return None


def _torsion_params(tj, tk, order_jk):
    """UFF torsion parameters for the central bond j-k.

    Returns ``(V, n, cos_n_phi0)`` in kcal/mol, or ``None`` if the case is
    outside the lean subset. ``V`` is the *total* barrier before dividing by
    the number of torsions about the bond.
    """
    hj, hk = _hyb_of(tj), _hyb_of(tk)
    vj, vk = UFF_PARAMS[tj], UFF_PARAMS[tk]

    if hj == "sp3" and hk == "sp3":
        # n = 3, phi0 = 180 deg -> cos(n phi0) = -1.
        return math.sqrt(vj["Vi"] * vk["Vi"]), 3, -1.0
    if {hj, hk} == {"sp2", "sp3"}:
        # Single bond, one sp2 one sp3: n = 6, phi0 = 0 -> cos = +1, V = 1.0.
        return 1.0, 6, 1.0
    if hj == "sp2" and hk == "sp2":
        # Conjugated / double bond: n = 2, phi0 = 180 -> cos = +1.
        v = 5.0 * math.sqrt(vj["Uj"] * vk["Uj"]) * (1.0 + 4.18 * math.log(order_jk))
        return v, 2, 1.0
    return None  # sp centers, group-16 specials, etc.not currently handled.


def _set_both(force, type_tuple, params, seen):
    """Register ``params`` on a HOOMD force under a type string and its reverse.

    The GSD snapshot's connection type-string ordering is not guaranteed to
    match ours; UFF bonded terms are symmetric under member reversal, so setting
    both orientations guarantees a match without changing any parameter value.
    """
    forward = "-".join(type_tuple)
    reverse = "-".join(type_tuple[::-1])
    force.params[forward] = params
    if reverse != forward:
        force.params[reverse] = params
    seen.update({type_tuple, type_tuple[::-1]})


def _hoomd_forces(top, snap, order_map, r_cut):
    """Build UFF HOOMD forces from a UFF-typed topology + snapshot."""
    import hoomd

    sites = list(top.sites)
    site_idx = {s: i for i, s in enumerate(sites)}

    def order_between(i, j):
        return order_map.get((min(i, j), max(i, j)), 1.0)

    particle_types = list(set(snap.particles.types))
    # UFF computes vdW for atoms >= 3 bonds apart: exclude 1-2 (bond) and 1-3
    # neighbors; 1-4 and beyond are included at full strength (UFF does not
    # scale 1-4, unlike OPLS). HOOMD's default nlist excludes only 'bond'.
    nlist = hoomd.md.nlist.Cell(buffer=0.4, exclusions=("bond", "1-3"))
    lj = hoomd.md.pair.LJ(nlist=nlist, default_r_cut=r_cut)
    two_pow = 2.0 ** (1.0 / 6.0)
    for a, t1 in enumerate(particle_types):
        for t2 in particle_types[a:]:
            x1, x2 = UFF_PARAMS[t1]["x1"], UFF_PARAMS[t2]["x1"]
            d1, d2 = UFF_PARAMS[t1]["D1"], UFF_PARAMS[t2]["D1"]
            # Geometric combining (UFF): r_min = sqrt(x_i x_j), D = sqrt(D_i D_j).
            sigma = math.sqrt(x1 * x2) / two_pow * ANG_TO_NM
            epsilon = math.sqrt(d1 * d2) * KCAL_TO_KJ
            lj.params[(t1, t2)] = dict(sigma=sigma, epsilon=epsilon)
    forces = [lj]

    if top.n_bonds > 0:
        bond_force = hoomd.md.bond.Harmonic()
        seen = set()
        for bond in top.bonds:
            m0, m1 = bond.connection_members
            i, j = site_idx[m0], site_idx[m1]
            ti, tj = sites[i].name, sites[j].name
            if (ti, tj) in seen:
                continue
            key = tuple(sorted((ti, tj)))
            order = order_between(i, j)
            r_ij = _bond_rest_length(key[0], key[1], order)
            k = _bond_force_constant(key[0], key[1], r_ij)
            params = dict(k=k * BOND_K_TO_KJ_NM2, r0=r_ij * ANG_TO_NM)
            # Register both orderings: the GSD snapshot's type-string ordering is
            # not guaranteed to match ours, and the potential is symmetric.
            _set_both(bond_force, (ti, tj), params, seen)
        forces.append(bond_force)

    if top.n_angles > 0:
        angle_force = hoomd.md.angle.Harmonic()
        seen = set()
        for angle in top.angles:
            m0, m1, m2 = angle.connection_members
            i, j, k_idx = site_idx[m0], site_idx[m1], site_idx[m2]
            ti, tj, tk = sites[i].name, sites[j].name, sites[k_idx].name
            if (ti, tj, tk) in seen:
                continue
            theta0 = math.radians(UFF_PARAMS[tj]["theta0"])
            # Harmonic approximation is invalid at linear centers (sin -> 0).
            if abs(math.sin(theta0)) < 1e-3:
                logger.warning(
                    "UFF: skipping linear angle type %s-%s-%s (theta0=180); the "
                    "lean harmonic approximation cannot represent it.",
                    ti, tj, tk,
                )
                seen.update({(ti, tj, tk), (tk, tj, ti)})
                continue
            r_ij = _bond_rest_length(ti, tj, order_between(i, j))
            r_jk = _bond_rest_length(tj, tk, order_between(j, k_idx))
            k_theta = _angle_force_constant(ti, tj, tk, r_ij, r_jk, theta0)
            params = dict(k=k_theta * KCAL_TO_KJ, t0=theta0)
            _set_both(angle_force, (ti, tj, tk), params, seen)
        forces.append(angle_force)

    if top.n_dihedrals > 0:
        # Count torsions about each central bond for the multiplicity divisor.
        mult = {}
        for dih in top.dihedrals:
            _, mj, mk, _ = dih.connection_members
            bond_key = tuple(sorted((site_idx[mj], site_idx[mk])))
            mult[bond_key] = mult.get(bond_key, 0) + 1

        dihedral_force = hoomd.md.dihedral.Periodic()
        seen = set()
        for dih in top.dihedrals:
            m0, m1, m2, m3 = dih.connection_members
            i, j, k_idx, ll = (
                site_idx[m0], site_idx[m1], site_idx[m2], site_idx[m3]
            )
            ti, tj, tk, tl = (
                sites[i].name, sites[j].name, sites[k_idx].name, sites[ll].name
            )
            if (ti, tj, tk, tl) in seen:
                continue
            params = _torsion_params(tj, tk, order_between(j, k_idx))
            if params is None:
                seen.update({(ti, tj, tk, tl), (tl, tk, tj, ti)})
                continue
            v_total, n, cos_n_phi0 = params
            bond_key = tuple(sorted((j, k_idx)))
            v = v_total / mult[bond_key]
            if abs(v) < 1e-9:
                seen.update({(ti, tj, tk, tl), (tl, tk, tj, ti)})
                continue
            # UFF: E = 0.5 V (1 - cos(n phi0) cos(n phi)).
            # HOOMD Periodic: E = k (1 + d cos(n phi - phi0)); phi0 = 0,
            # k = 0.5 V, d = -cos(n phi0) (which is +/-1).
            params = dict(
                k=0.5 * v * KCAL_TO_KJ,
                d=int(round(-cos_n_phi0)),
                n=n,
                phi0=0.0,
            )
            _set_both(dihedral_force, (ti, tj, tk, tl), params, seen)
        forces.append(dihedral_force)

    return forces


def _openmm_forces(top, order_map, r_cut):
    """Build UFF forces as a list of openmm.Force objects.

    The forces are added to an openmm.System by the caller. Each force
    references atom indices in site order, so the System's particles must be
    added in that same order.

    Force types:
      bonds use HarmonicBondForce.
      angles use CustomAngleForce with the UFF cosine-Fourier energy, plus a
        1 + cos(theta) form for linear centers.
      torsions use PeriodicTorsionForce.
      vdW uses CustomNonbondedForce with geometric combining of sigma and
        epsilon. Pairs 1-2 and 1-3 are excluded; 1-4 and beyond are included at
        full strength.

    Inversion (improper) terms are not generated.

    Parameters
    ----------
    top : gmso.Topology
        UFF-typed topology with connections identified.
    order_map : dict
        (i, j) -> bond order map, as returned by assign_uff_types.
    r_cut : float
        Nonbonded cutoff (nm), used only when top.box is set.
    """
    import openmm

    sites = list(top.sites)
    site_idx = {s: i for i, s in enumerate(sites)}

    def order_between(i, j):
        return order_map.get((min(i, j), max(i, j)), 1.0)

    forces = []
    two_pow = 2.0 ** (1.0 / 6.0)

    # --- vdW (LJ 12-6) with UFF geometric combining ---------------------------
    vdw = openmm.CustomNonbondedForce(
        "4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
        " sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)"
    )
    vdw.addPerParticleParameter("sigma")
    vdw.addPerParticleParameter("epsilon")
    for s in sites:
        p = UFF_PARAMS[s.name]
        vdw.addParticle(
            [p["x1"] / two_pow * ANG_TO_NM, p["D1"] * KCAL_TO_KJ]
        )
    # UFF computes vdW for atoms >= 3 bonds apart: exclude 1-2 and 1-3.
    vdw.createExclusionsFromBonds([(i, j) for (i, j) in order_map], 2)
    if top.box is not None:
        vdw.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        vdw.setCutoffDistance(r_cut)
    else:
        vdw.setNonbondedMethod(openmm.CustomNonbondedForce.NoCutoff)
    forces.append(vdw)

    # --- Bonds (harmonic) -----------------------------------------------------
    if top.n_bonds > 0:
        bond_force = openmm.HarmonicBondForce()
        for bond in top.bonds:
            m0, m1 = bond.connection_members
            i, j = site_idx[m0], site_idx[m1]
            ti, tj = sites[i].name, sites[j].name
            r_ij = _bond_rest_length(ti, tj, order_between(i, j))
            k = _bond_force_constant(ti, tj, r_ij)
            bond_force.addBond(
                i, j, r_ij * ANG_TO_NM, k * BOND_K_TO_KJ_NM2
            )
        forces.append(bond_force)

    # --- Angles (exact UFF cosine-Fourier, + linear form) ---------------------
    if top.n_angles > 0:
        general = openmm.CustomAngleForce(
            "K*(C0 + C1*cos(theta) + C2*cos(2*theta))"
        )
        for name in ("K", "C0", "C1", "C2"):
            general.addPerAngleParameter(name)
        linear = openmm.CustomAngleForce("K*(1 + cos(theta))")
        linear.addPerAngleParameter("K")
        for angle in top.angles:
            m0, m1, m2 = angle.connection_members
            i, j, k_idx = site_idx[m0], site_idx[m1], site_idx[m2]
            ti, tj, tk = sites[i].name, sites[j].name, sites[k_idx].name
            theta0 = math.radians(UFF_PARAMS[tj]["theta0"])
            r_ij = _bond_rest_length(ti, tj, order_between(i, j))
            r_jk = _bond_rest_length(tj, tk, order_between(j, k_idx))
            k_ang = _angle_force_constant(ti, tj, tk, r_ij, r_jk, theta0)
            k_ang *= KCAL_TO_KJ
            if abs(math.sin(theta0)) < 1e-3:  # linear center
                linear.addAngle(i, j, k_idx, [k_ang])
            else:
                c = math.cos(theta0)
                c2 = 1.0 / (4.0 * math.sin(theta0) ** 2)
                c1 = -4.0 * c2 * c
                c0 = c2 * (2.0 * c * c + 1.0)
                general.addAngle(i, j, k_idx, [k_ang, c0, c1, c2])
        if general.getNumAngles() > 0:
            forces.append(general)
        if linear.getNumAngles() > 0:
            forces.append(linear)

    # --- Torsions (principal cases) -------------------------------------------
    if top.n_dihedrals > 0:
        mult = {}
        for dih in top.dihedrals:
            _, mj, mk, _ = dih.connection_members
            key = tuple(sorted((site_idx[mj], site_idx[mk])))
            mult[key] = mult.get(key, 0) + 1

        torsion = openmm.PeriodicTorsionForce()
        for dih in top.dihedrals:
            m0, m1, m2, m3 = dih.connection_members
            i, j, k_idx, ll = (
                site_idx[m0], site_idx[m1], site_idx[m2], site_idx[m3]
            )
            tj, tk = sites[j].name, sites[k_idx].name
            params = _torsion_params(tj, tk, order_between(j, k_idx))
            if params is None:
                continue
            v_total, n, cos_n_phi0 = params
            v = v_total / mult[tuple(sorted((j, k_idx)))]
            if abs(v) < 1e-9:
                continue
            # UFF: E = 0.5 V (1 - cos(n phi0) cos(n phi)).
            # OpenMM: k (1 + cos(n phi - phase)); k = 0.5 V,
            # phase = 0 when cos(n phi0) = -1 else pi.
            phase = 0.0 if cos_n_phi0 < 0 else math.pi
            torsion.addTorsion(
                i, j, k_idx, ll, n, phase, 0.5 * v * KCAL_TO_KJ
            )
        if torsion.getNumTorsions() > 0:
            forces.append(torsion)

    return forces


