"""mBuild recipe for a single-walled carbon nanotube."""

from math import gcd

import numpy as np

from mbuild import Box, Compound, Particle, Port, clone, force_overlap
from mbuild.coordinate_transform import unit_vector
from mbuild.utils.validation import assert_port_exists

CC_BOND = 0.142
GRAPHENE_A = CC_BOND * np.sqrt(3)

# C-X bond lengths (nm) for the single-atom caps. C-O is a carbonyl.
CAP_BOND_LENGTHS = {
    "H": 0.109,
    "O": 0.122,
    "F": 0.135,
    "N": 0.140,
    "Cl": 0.177,
    "S": 0.182,
}


class CarbonNanotube(Compound):
    """A single-walled carbon nanotube of arbitrary chirality.

    The tube is built by tiling the graphene unit cell over the (n, m)
    nanotube unit cell and rolling the resulting sheet onto a cylinder.
    The tube axis is aligned with the z-axis and centered on x = y = 0.

    The chiral indices (n, m) say how the graphene sheet is wrapped: the sheet
    is rolled so that the atom at n * a1 + m * a2 lands on the atom at the
    origin, where a1 and a2 are the graphene lattice vectors. That vector wraps
    once around the tube, so its length is the circumference and larger indices
    give a wider tube. (n, n) gives an armchair tube, (n, 0) a zigzag tube, and
    any other pair a chiral tube, whose carbon rings spiral along the axis.

    Parameters
    ----------
    n : int, optional, default=5
        First chiral index, the number of a1 steps around the circumference.
        Must be positive.
    m : int, optional, default=5
        Second chiral index, the number of a2 steps around the circumference.
        Must be between 0 and n; (n, m) tubes with m > n are mirror images of
        the (m, n) tube.
    length : float, optional, default=2.0
        Approximate length of the tube in nm. The tube is built from a whole
        number of unit cells, so the actual length is the nearest multiple of
        the unit cell length (and is never shorter than one unit cell).
    radius : float, optional, default=None
        Target radius of the tube in nm. If given, (n, m) are chosen to give
        the closest achievable radius for the requested `chirality`, and the
        `n` and `m` arguments are ignored.
    chirality : str, optional, default="armchair"
        Either "armchair" ((n, n)) or "zigzag" ((n, 0)), named for the shape
        traced by the carbon rings around the circumference. Only used when
        `radius` is given.
    bond_tolerance : float, optional, default=0.02
        Bonds are added between carbons closer than CC_BOND + bond_tolerance.
    periodic : bool, optional, default=False
        If True, make the tube periodic along z and set the compound box.
        A periodic tube has no open ends, so `cap` is ignored.
    cap : str or Compound or None, optional, default="H"
        What to bond to the carbons at the open ends. Either an element symbol
        from CAP_BOND_LENGTHS, a Compound carrying a port labeled "up", or None
        to leave the ports exposed for the caller to fill. Capping a tube end
        terminates the dangling bonds but does not close off the bore.

    Attributes
    ----------
    n, m : int
        Chiral indices of the tube.
    radius : float
        Radius of the tube in nm.
    unit_cell_length : float
        Length of one nanotube unit cell along z in nm.
    cap_ports : list of Port
        Ports on the end carbons, empty once they have been capped.

    Examples
    --------
    A 2 nm long (6, 4) chiral tube, and an infinitely long (5, 5) armchair
    tube periodic along its axis.

    >>> from mbuild.lib.recipes import CarbonNanotube
    >>> chiral = CarbonNanotube(n=6, m=4, length=2.0)
    >>> armchair = CarbonNanotube(n=5, m=5, length=3.0, periodic=True)

    The chiral indices can also be chosen from a target radius.

    >>> zigzag = CarbonNanotube(radius=0.5, chirality="zigzag", length=3.0)
    >>> zigzag.n, zigzag.m
    (13, 0)

    Ends can be capped with another species, or left as ports to attach
    something larger.

    >>> fluorinated = CarbonNanotube(n=6, m=6, length=1.0, cap="F")
    >>> open_ended = CarbonNanotube(n=6, m=6, length=1.0, cap=None)
    >>> len(open_ended.cap_ports)
    24

    Notes
    -----
    Cutting the tube to length can leave carbons bonded to only one neighbor,
    which are not a physical tube edge. These are trimmed before capping, so
    every end carbon has two neighbors and one dangling bond.

    Adapted from the Nanotube-Builder mBuild recipe by M. Whitehead:
    https://github.com/whitehml/Nanotube-Builder

    References
    ----------
    .. [1] Dresselhaus, M. S., Dresselhaus, G., Saito, R. "Physics of carbon
           nanotubes." (1995) Carbon 33, 883-891
    """

    def __init__(
        self,
        n=5,
        m=5,
        length=2.0,
        radius=None,
        chirality="armchair",
        bond_tolerance=0.02,
        periodic=False,
        cap="H",
    ):
        super().__init__()

        if radius is not None:
            n, m = _indices_from_radius(radius, chirality)
        n, m = int(n), int(m)
        if n < 1:
            raise ValueError(f"`n` must be a positive integer, got {n}.")
        if not 0 <= m <= n:
            raise ValueError(f"`m` must satisfy 0 <= m <= n, got m={m}, n={n}.")
        if length <= 0:
            raise ValueError(f"`length` must be positive, got {length}.")
        if isinstance(cap, str) and cap not in CAP_BOND_LENGTHS:
            raise ValueError(
                f"'{cap}' is not one of the supported cap atoms "
                f"({', '.join(CAP_BOND_LENGTHS)}). Pass a Compound with a port "
                "labeled 'up' to attach anything else."
            )

        self.n = n
        self.m = m
        self.radius = float(GRAPHENE_A * np.sqrt(n**2 + n * m + m**2) / (2 * np.pi))

        cell_coords, self.unit_cell_length = _unit_cell(n, m)
        n_cells = max(1, round(length / self.unit_cell_length))

        for i in range(n_cells):
            for u, v in cell_coords:
                theta = 2 * np.pi * u
                pos = (
                    self.radius * np.cos(theta),
                    self.radius * np.sin(theta),
                    (v + i) * self.unit_cell_length,
                )
                self.add(Particle(name="C", element="C", pos=pos))

        if periodic:
            self.periodicity = (False, False, True)
            xy = 2 * (self.radius + CC_BOND)
            self.box = Box([xy, xy, n_cells * self.unit_cell_length])

        self.generate_bonds(
            name_a="C", name_b="C", dmin=0.0, dmax=CC_BOND + bond_tolerance
        )
        self.cap_ports = []
        if not periodic:
            # generate_bonds assigns a bounding box when the compound has none
            self.box = None
            self._trim_ends()
            self._add_cap_ports()
            if cap is not None:
                self._cap_ends(cap)

    def _trim_ends(self):
        """Remove end carbons left with fewer than two neighbors by the cut."""
        while True:
            lonely = [
                p
                for p in self.particles()
                if len(list(self.bond_graph.neighbors(p))) < 2
            ]
            if not lonely:
                break
            self.remove(lonely)
        # Compound.remove leaves a port on each particle it disconnects
        self.remove(self.all_ports())

    def _add_cap_ports(self):
        """Add a port along the missing sp2 direction of each end carbon."""
        for particle in list(self.particles()):
            neighbors = list(self.bond_graph.neighbors(particle))
            if len(neighbors) != 2:
                continue
            bonds = [unit_vector(n.pos - particle.pos) for n in neighbors]
            port = Port(
                anchor=particle,
                orientation=-unit_vector(bonds[0] + bonds[1]),
                separation=CC_BOND / 2,
            )
            self.add(port, label="cap_port[$]")
            self.cap_ports.append(port)

    def _cap_ends(self, cap):
        """Attach a copy of `cap` to every port."""
        for port in list(self.cap_ports):
            if isinstance(cap, str):
                capper = Compound(name=cap)
                capper.add(Particle(name=cap, element=cap))
                separation = CAP_BOND_LENGTHS[cap] - CC_BOND / 2
                capper.add(Port(anchor=capper[0], separation=separation), "up")
            else:
                capper = clone(cap)
                assert_port_exists("up", capper)
            self.add(capper)
            force_overlap(capper, capper.labels["up"], port)
        self.cap_ports = []


def _unit_cell(n, m):
    """Return fractional (circumferential, axial) coordinates of one unit cell."""
    d_r = gcd(2 * n + m, 2 * m + n)
    # Chiral vector C and translation vector T in the graphene basis
    c_hat = np.array([n, m])
    t_hat = np.array([(2 * m + n) / d_r, -(2 * n + m) / d_r])

    a1 = GRAPHENE_A * np.array([1.0, 0.0])
    a2 = GRAPHENE_A * np.array([0.5, np.sqrt(3) / 2])
    basis = np.array([a1, a2])
    c_vec = c_hat @ basis
    t_vec = t_hat @ basis

    # Fractional coordinates along (C, T) for a point r satisfy r = frac @ [C, T]
    to_frac = np.linalg.inv(np.array([c_vec, t_vec]).T)

    n_hex = 2 * (n**2 + n * m + m**2) // d_r
    sublattice = [np.zeros(2), (a1 + a2) / 3]

    # The unit cell is the parallelogram spanned by C and T; its corners bound
    # the range of lattice indices that can fall inside it.
    corners = np.array([[0, 0], c_hat, t_hat, c_hat + t_hat])
    lo = np.floor(corners.min(axis=0)).astype(int) - 1
    hi = np.ceil(corners.max(axis=0)).astype(int) + 1

    tol = 1e-6
    coords = []
    for i in range(lo[0], hi[0] + 1):
        for j in range(lo[1], hi[1] + 1):
            for offset in sublattice:
                frac = to_frac @ (i * a1 + j * a2 + offset)
                if np.all(frac >= -tol) and np.all(frac < 1 - tol):
                    coords.append(frac)

    if len(coords) != 2 * n_hex:
        raise RuntimeError(
            f"Found {len(coords)} atoms in the ({n}, {m}) unit cell, "
            f"expected {2 * n_hex}."
        )
    return np.array(coords), float(np.linalg.norm(t_vec))


def _indices_from_radius(radius, chirality):
    """Return the (n, m) pair whose radius is closest to `radius`."""
    if chirality not in ("armchair", "zigzag"):
        raise ValueError(
            f"`chirality` must be 'armchair' or 'zigzag', got '{chirality}'."
        )
    m_over_n = 1 if chirality == "armchair" else 0
    scale = float(GRAPHENE_A * np.sqrt(1 + m_over_n + m_over_n**2) / (2 * np.pi))
    n = round(radius / scale)
    if n < 1:
        raise ValueError(
            f"A radius of {radius} nm is smaller than the smallest {chirality} "
            f"tube ({scale:.3f} nm)."
        )
    return n, n * m_over_n
