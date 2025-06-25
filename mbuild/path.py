"""Molecular paths and templates"""

import math
from abc import abstractmethod

import freud
import numpy as np

from mbuild import Compound
from mbuild.utils.geometry import bounding_box

"""
A path is basically a bond graph with coordinates/positions
assigned to the nodes. This is kind of what Compound is already.

The interesting and challenging part is building up/creating the path.
This follows an algorithm to generate next coordinates.
Any random path generation algorithm will include
a rejection/acception step. We basically end up with
Monte Carlo. Some path algorithms won't be random (lamellae)

Is Path essentially going to be a simple Monte carlo-ish
class that others can inherit from then implement their own approach?

Classes that inherit from path will have their own
verions of next_coordinate(), check_path(), etc..
We can define abstract methods for these in Path.
We can put universally useful methods in Path as well (e.g., n_list(), to_compound(), etc..).

Some paths (lamellar structures) would kind of just do
everything in generate() without having to use
next_coordinate() or check_path(). These would still need to be
defined, but just left empty and/or always return True in the case of check_path.
Maybe that means these kinds of "paths" need a different data structure?

Do we just have RandomPath and DeterministicPath?

RandomPath ideas:
- Random walk (tons of possibilities here)
- Branching
- Multiple self-avoiding random walks

DeterministicPath ideas:
- Lamellar layers
- DNA strands

Some combination of these?
- Lamellar + random walk to generate semi-crystalline like structures?
-
"""


class Path:
    def __init__(self, N=None, bond_graph=None, coordinates=None):
        self.bond_graph = bond_graph
        # Only N is defined, make empty coordinates
        if N is not None and coordinates is None:
            self.N = N
            self.coordinates = np.zeros((N, 3))
        # Only coordinates is defined, set N from length
        elif coordinates is not None and N is None:
            self.N = len(coordinates)
            self.coordinates = coordinates
        # Neither is defined, use list for coordinates
        elif N is None and coordinates is None:
            self.N = N
            self.coordinates = []
        else:
            raise ValueError("Specify either one of N and coordinates, or neither")
        self.generate()

    @classmethod
    def from_coordinates(cls, coordinates, bond_graph=None):
        return cls(coordinates=coordinates, bond_graph=bond_graph, N=None)

    @abstractmethod
    def generate(self):
        """Abstract class for running a Path generation algorithm

        This method should:
        -----------------
            - Set initial conditions
            - Implement Path generation steps by calling _next_coordinate() and _check_path()
            - Update bonding info depending on specific path approach
                - Ex) Random walk will always bond consecutive beads together
            - Handle cases of next coordiante acceptance
            - Handle cases of next coordinate rejection
        """
        pass

    @abstractmethod
    def _next_coordinate(self):
        """Algorithm to generate the next coordinate in the path"""
        pass

    @abstractmethod
    def _check_path(self):
        """Algorithm to accept/reject trial move of the current path"""
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

    def to_compound(self, bead_name="Bead", bead_mass=1):
        """Visualize a path as an mBuild Compound"""
        compound = Compound()
        for xyz in self.coordinates:
            compound.add(Compound(name=bead_name, mass=bead_mass, pos=xyz))
        if self.bond_graph:
            compound.set_bond_graph(self.bond_graph)
        return compound

    def apply_mapping(self):
        # TODO: Finish, add logic to align orientation with path site pos and bond graph
        """Mapping other compounds onto a Path's coordinates"""
        pass

    def _path_history(self):
        """Maybe this is a method that can be used optionally.
        We could add a save_history parameter to __init__.
        Depending on the approach, saving histories might add additional computation time and resources.
        """
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
        bond_graph=None,
    ):
        self.bond_length = bond_length
        self.radius = radius
        self.min_angle = min_angle
        self.max_angle = max_angle
        self.seed = seed
        self.tolerance = tolerance
        self.max_attempts = max_attempts
        self.attempts = 0
        self.count = 0
        super(HardSphereRandomWalk, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        np.random.seed(self.seed)
        # With fixed bond lengths, the first move is always accepted
        phi = np.random.uniform(0, 2 * np.pi)
        theta = np.random.uniform(0, np.pi)
        next_pos = np.array(
            [
                self.bond_length * np.sin(theta) * np.cos(phi),
                self.bond_length * np.sin(theta) * np.sin(phi),
                self.bond_length * np.cos(theta),
            ]
        )
        self.coordinates[1] = next_pos
        self.count += 1  # We already have 1 accepted move
        while self.count < self.N - 1:
            new_xyz = self._next_coordinate(
                pos1=self.coordinates[self.count],
                pos2=self.coordinates[self.count - 1],
            )
            self.coordinates[self.count + 1] = new_xyz
            if self._check_path():
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

    def _next_coordinate(self, pos1, pos2):
        # Vector formed by previous 2 coordinates
        v1 = pos2 - pos1
        v1_norm = v1 / np.linalg.norm(v1)
        theta = np.random.uniform(self.min_angle, self.max_angle)
        # Pick random vector and center around origin (0,0,0)
        r = np.random.rand(3) - 0.5
        r_perp = r - np.dot(r, v1_norm) * v1_norm
        r_perp_norm = r_perp / np.linalg.norm(r_perp)
        # New vector, rotated relative to v1
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


class Lamellar(Path):
    def __init__(
        self, num_layers, layer_separation, layer_length, bond_length, bond_graph=None
    ):
        self.num_layers = num_layers
        self.layer_separation = layer_separation
        self.layer_length = layer_length
        self.bond_length = bond_length
        super(Lamellar, self).__init__(N=None, bond_graph=bond_graph)

    def generate(self):
        layer_spacing = np.arange(0, self.layer_length, self.bond_length)
        # Info for generating coords of the curves between layers
        r = self.layer_separation / 2
        arc_length = r * np.pi
        arc_num_points = math.floor(arc_length / self.bond_length)
        arc_angle = np.pi / (arc_num_points + 1)  # incremental angle
        arc_angles = np.linspace(arc_angle, np.pi, arc_num_points, endpoint=False)
        for i in range(self.num_layers):
            if i % 2 == 0:  # Even layer; build from left to right
                layer = [
                    np.array([self.layer_separation * i, y, 0]) for y in layer_spacing
                ]
                # Mid-point between this and next layer; use to get curve coords.
                origin = layer[-1] + np.array([r, 0, 0])
                arc = [
                    origin + np.array([-np.cos(theta), np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            else:  # Odd layer; build from right to left
                layer = [
                    np.array([self.layer_separation * i, y, 0])
                    for y in layer_spacing[::-1]
                ]
                # Mid-point between this and next layer; use to get curve coords.
                origin = layer[-1] + np.array([r, 0, 0])
                arc = [
                    origin + np.array([-np.cos(theta), -np.sin(theta), 0]) * r
                    for theta in arc_angles
                ]
            if i != self.num_layers - 1:
                self.coordinates.extend(layer + arc)
            else:
                self.coordinates.extend(layer)

    def _next_coordinate(self):
        pass

    def _check_path(self):
        pass


class StraightLine(Path):
    def __init__(self, spacing, N, direction=(1, 0, 0), bond_graph=None):
        self.spacing = spacing
        self.direction = np.asarray(direction)
        super(StraightLine, self).__init__(N=N, bond_graph=bond_graph)

    def generate(self):
        self.coordinates = np.array(
            [np.zeros(3) + i * self.spacing * self.direction for i in range(self.N)]
        )

    def _next_coordinate(self):
        pass

    def _check_path(self):
        pass
