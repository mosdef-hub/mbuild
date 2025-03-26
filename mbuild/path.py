"""Molecular paths and templates"""

from abc import ABC, abstractmethod

import freud
import numpy as np

from mbuild import Compound
from mbuild.utils.geometry import bounding_box


class Path(ABC):
    def __init__(self, N=None, max_attempts=10000):
        self.max_attempts = max_attempts
        self.attempts = 0
        self.compound = None
        self.N = N
        self.nlist = None
        if N:
            self.coordinates = np.zeros((N, 3))
        # Not every path will know N ahead of time (Lamellae)
        else:
            self.coordinates = []
        # Generate dict values now?
        # Entry for each node, initialize with empty lists
        self.bonds = []

    """
    A path is basically a bond graph with coordinates/positions
    assigned to the nodes. This is kind of what Compound is already.

    The interesting part is building up/creating the path.
    This follows an algorithm to generate next coordinates.
    Any useful path generation algorithm will include
    a rejection/acception step. Basically end up with
    monte carlo.

    Is Path basically going to be a simple Monte carlo
    class that others can inherit from?

    Classes that inherit from path will have their own
    verions of next_coordinate, check path, etc..
    We can define abstract methods here.

    """

    @abstractmethod
    def generate(self):
        pass

    def neighbor_list(self, r_max, coordinates=None, box=None):
        if coordinates is None:
            coordinates = self.coordinates
        if box is None:
            box = bounding_box(coordinates)
        freud_box = freud.box.Box(Lx=box[0], Ly=box[1], Lz=box[2])
        aq = freud.locality.AABBQuery(freud_box, coordinates)
        aq_query = aq.query(
            query_points=coordinates,
            query_args=dict(r_min=0.0, r_max=r_max, exclude_ii=True),
        )
        nlist = aq_query.toNeighborList()
        return nlist

    def to_compound(self):
        compound = Compound()
        for xyz in self.coordinates:
            compound.add(Compound(name="Bead", pos=xyz))
        for bond_group in self.bonds:
            compound.add_bond([compound[bond_group[0]], compound[bond_group[1]]])
        return compound

    def apply_mapping(self):
        pass

    def _path_history(self):
        pass

    @abstractmethod
    def _next_coordinate(self):
        pass

    @abstractmethod
    def _check_path(self):
        pass


class HardSphereRandomWalk(Path):
    def __init__(
        self,
        N,
        bond_length,
        radius,
        min_angle,
        max_angle,
        max_attempts,
        seed,
        tolerance=1e-5,
    ):
        self.bond_length = bond_length
        self.radius = radius
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.seed = seed
        self.tolerance = tolerance
        self.count = 0
        super(HardSphereRandomWalk, self).__init__(N=N, max_attempts=max_attempts)

    def generate(self):
        np.random.seed(self.seed)
        # With fixed bond lengths, the first move is always accepted
        self.coordinates[1] = self._next_coordinate(pos1=self.coordinates[0])
        self.bonds.append([0, 1])
        self.count += 1  # We already have 1 accepted move
        while self.count < self.N - 1:
            new_xyz = self._next_coordinate(
                pos1=self.coordinates[self.count],
                pos2=self.coordinates[self.count - 1],
            )
            self.coordinates[self.count + 1] = new_xyz
            if self._check_path():
                self.bonds.append((self.count, self.count + 1))
                self.count += 1
                self.attempts += 1
            else:
                self.coordinates[self.count + 1] = np.zeros(3)
                self.attempts += 1
            if self.attempts == self.max_attempts and self.count < self.N:
                raise RuntimeError(
                    "The maximum number attempts allowed have passed, and only ",
                    f"{self.count} sucsessful attempts were completed.",
                    "Try changing the parameters and running again.",
                )

    def _next_coordinate(self, pos1, pos2=None):
        if pos2 is None:
            phi = np.random.uniform(0, 2 * np.pi)
            theta = np.random.uniform(0, np.pi)
            next_pos = np.array(
                [
                    self.bond_length * np.sin(theta) * np.cos(phi),
                    self.bond_length * np.sin(theta) * np.sin(phi),
                    self.bond_length * np.cos(theta),
                ]
            )
        else:  # Get the last bond vector, use angle range with last 2 coords.
            v1 = pos2 - pos1
            v1_norm = v1 / np.linalg.norm(v1)
            theta = np.random.uniform(self.min_angle, self.max_angle)
            r = np.random.rand(3) - 0.5
            r_perp = r - np.dot(r, v1_norm) * v1_norm
            r_perp_norm = r_perp / np.linalg.norm(r_perp)
            v2 = np.cos(theta) * v1_norm + np.sin(theta) * r_perp_norm
            next_pos = v2 * self.bond_length

        return pos1 + next_pos

    def _check_path(self):
        """Use neighbor_list to check for pairs within a distance smaller than the radius"""
        # Grow box size as number of steps grows
        box_length = self.count * self.radius * 2.01
        # Only need neighbor list for accepted moves + current trial move
        coordinates = self.coordinates[: self.count + 2]
        nlist = self.neighbor_list(
            coordinates=coordinates,
            r_max=self.radius - self.tolerance,
            box=[box_length, box_length, box_length],
        )
        if len(nlist.distances) > 0:  # Particle pairs found within the particle radius
            return False
        else:
            return True


class Lamellae(Path):
    def __init__(self):
        super(Lamellae, self).__init__()
