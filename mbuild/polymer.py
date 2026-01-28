"""Recipe for an mBuild polymer."""

import itertools as it

import numpy as np

from mbuild import clone
from mbuild.compound import Compound
from mbuild.coordinate_transform import (
    force_overlap,
    x_axis_transform,
    y_axis_transform,
    z_axis_transform,
)
from mbuild.lib.atoms import H
from mbuild.port import Port
from mbuild.simulation import energy_minimize as e_min
from mbuild.utils.validation import assert_port_exists

__all__ = ["Polymer"]


class Polymer(Compound):
    """Connect one or more components in a specified sequence.

    Attributes
    ----------
    monomers : list of mbuild.Compounds
        The compound(s) to replicate. Add to this list using the add_monomers
        method.
    end_groups : list of mbuild.Compounds
        The compound to cap the end of the polymer. Add to this list using the
        add_end_groups method.

    Methods
    -------
    add_monomer(monomer, indices, separation, port_labels, orientation, replace)
        Use to add a monomer compound to Polymer.monomers

    add_end_groups(compound, index, separation, orientation, replace)
        Use to add an end group compound to Polymer.end_groups

    build_from_path(path, sequence, add_hydrogens, bond_head_tail, energy_minimize)
        Used to set the polymer configuration from a pre-determined mbuild.path.Path
        Use this to create polymers that follow a random walk, or form
        lamellar order, ring polymers and more. See the ``mbuild.path`` module.

    build(n, sequence)
        Use to create a single polymer compound. This method uses the compounds
        created by calling the add_monomer and add_end_group methods.

    Notes
    -----
    There are two different approaches to using the Polymer class to create
    polymers

    1) Pass in already created Compound instances to the monomers and end_groups
    parameters when creating a Polymer instance:

        You can then call the Polymer.build() method to create a polymer.
        This approach can be used if the compounds being passed into the Polymer
        instance already have the ports created, and correct atomic structure to
        allow for the monomer-monomer and monomer-end group bonds. These
        compounds are used as-is when creating the polymer chain.

        Example
        -------
        >>> chain = Polymer(
        ...     monomers=[mb.Compound],
        ...     end_groups = [mb.Compound, mb.Compound]
        ... )
        >>> chain.build(n=5)

    2) Use the add_monomer() and add_end_group() methods:

        These functions are there to help with the creation of mb.Ports, which
        are required by the build() method, and they help with the removal of
        any atoms (hydrogens) that are occupying what should be the monomer-
        monomer and monomer-end group bonding sites.  With this approach, create
        a Polymer() instance, then call the add_monomer() and add_end_group()
        methods before calling the build() method.

        Example
        -------
        >>> chain = Polymer()
        >>> chain.add_monomer(mb.Compound)
        >>> chain.add_end_groups(mb.Compound)
        >>> chain.build(n=5)

        Refer to the method specific doc strings to see the correct use.
    """

    def __init__(self, monomers=None, end_groups=None):
        super(Polymer, self).__init__()
        self._monomers = monomers or []
        self._end_groups = end_groups or [None, None]
        if not isinstance(self._end_groups, list):
            raise ValueError(
                "Please provide two end groups in a list; "
                f"you provided {self._end_groups}"
            )
        elif len(self._end_groups) != 2:
            raise ValueError(
                f"Please provide two end groups; you provided {len(self._end_groups)}"
            )
        self._port_labels = ["up", "down"]
        self._headtail = [None, None]
        self.head_port = None
        self.tail_port = None

    @property
    def monomers(self):
        """Get the monomers.

        monomers cannot be set. Use add_monomer method instead.
        """
        return self._monomers

    @property
    def end_groups(self):
        """Get the end groups.

        end_groups cannot be set. Use add_end_group method instead.
        """
        return self._end_groups

    def build_from_path(
        self,
        path,
        sequence="A",
        add_hydrogens=True,
        bond_head_tail=False,
        energy_minimize=False,
    ):
        """Build the polymer chain to a pre-defined configuraiton.

        See ``mbuild.path.Path`` to available polymer configurations
        or use mbuild.path.Path.from_coordinates() to manually set
        a polymer chan configuraiton.

        Parameters
        ----------
        path : mbuild.path.Path, required
            The configuration the polymer will be mapped to after connecting
            all of the components (monomers and end groups).
            The path will determine the number of monomer sites.
        sequence : str, optional, default 'A'
            A string of characters where each unique character represents one
            repetition of a monomer. Characters in `sequence` are assigned to
            monomers in the order they appear in `Polymer.monomers`.
            The characters in `sequence` are assigned to the compounds in the
            in the order that they appear in the Polymer.monomers list.
            For example, 'AB' where 'A'corresponds to the first compound
            added to Polymer.monomers and 'B' to the second compound.
        add_hydrogens : bool, default True
            If True and an ``end_groups`` compound is None, then the head or tail
            of the polymer will be capped off with hydrogen atoms. If end group
            compounds exist, then they will be used.
            If False and an end group compound is None, then the head or tail
            port will be exposed in the polymer.
        bond_head_tail : bool, default False
            If set to ``True``, then a bond between the head and tail groups (monomers or end groups)
            is created. This does not create a periodic bond (see Polymer.create_periodic_bond).
            The ``add_hydrogens`` parameter must be set to ``False`` to create this bond.
            This is useful for creating ring polymers, knot polymers.
            See ``mbuild.path.Cyclic`` and ``mbuild.path.Knot``.
        energy_minimize : bool, default True
            If ``True`` then relax the bonds and angles that may be distorted from mapping
            the atomistic polymer to a path.
            This uses the capped displacement methods in ``mbuild.simulation``.

        Notes
        -----
        Energy Minimization:

            It is required to run energy minimization on chains built from a path to relax
            the topology to a suitable starting point.
            mBuild contains multiple energy minimization approaches in the ``mbuild.simulation`` module.

            When ``energy_minimize`` is set to ``True``, this method uses the OpenBabel based energy
            minimization method with the UFF force field. We have found that once polymers
            reach a size on the order of ~500 atoms, significantly faster energy minimization can be
            achieved using one of the hoomd-based methods in ``mbuild.simulation``.
            In that case, set ``energy_minimize=False`` and pass the resulting polymer
            compound into one of these methods.
            See: ``mbuild.simulation.hoomd_cap_displacement``, and ``mbuild.simulation.hoomd_fire``.

        """
        n = len(path.coordinates) - sum([1 for i in self.end_groups if i is not None])
        self.build(
            n=n // len(sequence),
            sequence=sequence,
            add_hydrogens=add_hydrogens,
            bond_head_tail=bond_head_tail,
            coordinates=path.coordinates,
        )
        if energy_minimize:
            e_min(self)

    def set_monomer_positions(self, coordinates, energy_minimize=False):
        """Shift monomers so that their center of mass matches a set of pre-defined coordinates.

        Parameters
        ----------
        coordinates : np.ndarray, shape=(N,3)
            Set of x,y,z coordinatess

        Notes
        -----
        Energy Minimization:

            It is required to run energy minimization on chains built from a path to relax
            the topology to a suitable starting point.
            mBuild contains multiple energy minimization approaches in the ``mbuild.simulation`` module.

            When ``energy_minimize`` is set to ``True``, this method uses the OpenBabel based energy
            minimization method with the UFF force field. We have found that once polymers
            reach a size on the order of ~500 atoms, significantly faster energy minimization can be
            achieved using one of the hoomd-based methods in ``mbuild.simulation``
            In that case, set ``energy_minimize=False`` and pass the resulting polymer
            compound into one of these methods.
            See: `hoomd_cap_displacement`, and `hoomd_fire`.
        """
        for i, xyz in enumerate(coordinates):
            self.children[i].translate_to(xyz)
        if energy_minimize:
            e_min(self)

    def build(
        self,
        n,
        sequence="A",
        add_hydrogens=True,
        bond_head_tail=False,
        coordinates=None,
    ):
        """Connect one or more components in a specified sequence.

        Uses the compounds that are stored in Polymer.monomers and
        Polymer.end_groups.

        If no capping method is used, i.e., if ``add_hydrogens == False``
        and ``Polymer.end_groups is None``, then the ports are exposed as,
        ``Polymer.head_port`` and ``Polymer.tail_port``.

        Parameters
        ----------
        n : int
            The number of times to replicate the sequence.
        sequence : str, optional, default 'A'
            A string of characters where each unique character represents one
            repetition of a monomer. Characters in `sequence` are assigned to
            monomers in the order they appear in `Polymer.monomers`.
            The characters in `sequence` are assigned to the compounds in the
            in the order that they appear in the Polymer.monomers list.
            For example, 'AB' where 'A'corresponds to the first compound
            added to Polymer.monomers and 'B' to the second compound.
        add_hydrogens : bool, default True
            If True and an ``end_groups`` compound is None, then the head or tail
            of the polymer will be capped off with hydrogen atoms. If end group
            compounds exist, then they will be used.
            If False and an end group compound is None, then the head or tail
            port will be exposed in the polymer.
        bond_head_tail : bool, default False
            If set to ``True``, then a bond between the head and tail groups (monomers or end groups)
            is created. This does not create a periodic bond (see Polymer.create_periodic_bond).
            The ``add_hydrogens`` parameter must be set to ``False`` to create this bond.
            This is useful for creating ring polymers, knot polymers.
            See ``mbuild.path.Cyclic`` and ``mbuild.path.Knot`` and ``Polymer.build_from_path``.
        coordinates : np.ndarray (n, 3)
            Manually pass in the monomer coordinates.
            Each monomer (including head and end groups) will be translated to the
            corresponding coordinate. A rotatation will be performed to attempt to
            align the monomer head-tail vector with the neighboring monomer orientations.
            Further energy minimization may be required to remove overlapping particles and relax bonds.
            See ``mbuild.simulation.hoomd_cap_displacement`` and ``mbuild.simulation.hoomd_fire``.
            See ``Polyer.build_from_path`` to build a polymer chain from an ``mbuild.path.Path``.

        Notes
        -----
        The chain conformations obtained from this method are often difficult to use
        in later steps of creating systems of polymers. See the alternative build
        method ``Polymer.build_from_path``.
        """
        if add_hydrogens and bond_head_tail:
            raise ValueError(
                "In order to bond head and tail ports, the Polymer instance cannot contain "
                "end_groups and add_hydrogens must be set to `False`"
            )
        if n < 1:
            raise ValueError("n must be 1 or more")
        head_and_tail_compounds = []
        repeat_compounds = []

        for monomer in self._monomers:
            for label in self._port_labels:
                assert_port_exists(label, monomer)

        unique_seq_ids = sorted(set(sequence))

        if len(self._monomers) != len(unique_seq_ids):
            raise ValueError(
                "Number of monomers passed to `Polymer` class must match "
                "number of unique entries in the specified sequence."
            )

        # 'A': monomer_1, 'B': monomer_2....
        seq_map = dict(zip(unique_seq_ids, self._monomers))
        last_part = None
        site_count = 0 if not self._end_groups[0] else 1
        for n_added, seq_item in enumerate(it.cycle(sequence)):
            this_part = clone(seq_map[seq_item])
            repeat_compounds.append(this_part)
            if last_part is None:
                first_part = this_part
                if coordinates is not None:
                    this_part.translate_to(coordinates[site_count])
                    this_part_head = this_part.labels[self._port_labels[1]].anchor
                    this_part_tail = this_part.labels[self._port_labels[0]].anchor
                    v1 = this_part_head.pos - this_part_tail.pos
                    v1 /= np.linalg.norm(v1)
                    v2 = coordinates[site_count + 1] - coordinates[site_count]
                    v2 /= np.linalg.norm(v2)
                    normal = np.cross(v1, v2)
                    angle = np.arccos(
                        v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    )
                    if angle > np.pi / 2:
                        angle = np.pi - angle
                    # Center of mass needs to be at origin for rotation
                    this_part.translate_to((0, 0, 0))
                    this_part.rotate(around=normal, theta=angle)
                    this_part.translate_to(coordinates[site_count])
            else:
                # Transform this part, such that its bottom port is rotated
                # and translated to the last parts top port.
                force_overlap(
                    move_this=this_part,
                    from_positions=this_part.labels[self._port_labels[0]],
                    to_positions=last_part.labels[self._port_labels[1]],
                )
                if coordinates is not None:
                    try:
                        this_part.translate_to(coordinates[site_count])
                        this_part_head = this_part.labels[self._port_labels[1]].anchor
                        this_part_tail = this_part.labels[self._port_labels[0]].anchor
                        # Get this parts head-tail vector (head port pos - tail port pos)
                        v1 = this_part_head.pos - this_part_tail.pos
                        v1 /= np.linalg.norm(v1)
                        # v2 is the vector between the previous site pos and next site pos
                        v2 = coordinates[site_count + 1] - coordinates[site_count - 1]
                        v2 /= np.linalg.norm(v2)
                        normal = np.cross(v1, v2)
                        angle = np.arccos(
                            v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                        )
                        if angle > np.pi / 2:
                            angle = np.pi - angle
                        # Center of mass needs to be at origin for rotation
                        this_part.translate_to((0, 0, 0))
                        this_part.rotate(around=normal, theta=angle)
                        this_part.translate_to(coordinates[site_count])
                    except IndexError:
                        pass
            last_part = this_part
            site_count += 1
            if n_added == n * len(sequence) - 1:
                break
        self.head_port = first_part["up"] if not first_part["up"].used else None
        self.tail_port = last_part["down"] if not last_part["down"].used else None

        head_tail = [self.head_port, self.tail_port]
        for i, compound in enumerate(self._end_groups):
            if compound is not None:
                if self._headtail[i] is not None:
                    head_tail[i].update_separation(self._headtail[i])
                anchor_particle = compound.labels["up"].anchor
                force_overlap(
                    move_this=compound,
                    from_positions=compound.labels["up"],
                    to_positions=head_tail[i],
                )
                if coordinates is not None:
                    head = anchor_particle.pos
                    tail = compound.center
                    v1 = head - tail
                    v1 /= np.linalg.norm(v1)
                    if i == 0:  # v2 is between first monomer 0 and end group
                        end_group_index = 0
                        v2 = coordinates[1] - coordinates[0]
                    elif i == 1:  # v2 is between end group and last monomer
                        end_group_index = -1
                        v2 = coordinates[-1] - coordinates[-2]
                    v2 /= np.linalg.norm(v2)
                    normal = np.cross(v1, v2)
                    angle = np.arccos(
                        v1.dot(v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                    )
                    if angle > np.pi / 2:
                        angle = np.pi - angle
                    # Center of mass needs to be at origin for rotation
                    compound.translate_to((0, 0, 0))
                    compound.rotate(around=normal, theta=angle)
                    compound.translate_to(coordinates[end_group_index])

                head_and_tail_compounds.append(compound)
                head_tail[i] = None
            else:
                if add_hydrogens:
                    hydrogen = H()
                    # Defaut to 1/2 H-C bond len
                    head_tail[i].update_separation(0.0547)
                    hydrogen["up"].update_separation(0.0547)
                    head_tail[i].parent.add(hydrogen)
                    force_overlap(hydrogen, hydrogen["up"], head_tail[i])
                    head_tail[i] = None
                else:
                    # if None, hoist port to polymer level
                    self.add(
                        head_tail[i],
                        self._port_labels[i],
                        containment=False,
                    )

            if head_tail[i] is None and i == 0:
                self.head_port = None
            elif head_tail[i] is None and i == 1:
                self.tail_port = None

        port_ids = [
            id(x) for x in self.available_ports()
        ]  # prevent overlooking down port and incorrectly removing
        for port in self.all_ports():
            if id(port) not in port_ids:
                self.remove(port)
        if self.end_groups[0]:
            self.add(self.end_groups[0])
        for compound in repeat_compounds:
            self.add(new_child=compound, label="monomer[$]")
        if self.end_groups[1]:
            self.add(self.end_groups[1])

        if bond_head_tail:
            force_overlap(self, self.head_port, self.tail_port)

    def add_monomer(
        self,
        compound,
        head_tag=None,
        tail_tag=None,
        separation=None,
        head_orientation=None,
        tail_orientation=None,
        bond_order=1,
        remove_hydrogens=True,
    ):
        """Add a Compound to self.monomers.

        The monomers will be used to build the polymer. Call this function for
        each unique monomer to be used in the polymer.

        Notes
        -----
        Using the 'replace' and 'indices' parameters:

        The atoms in an mbuild compound can be identified by their index
        numbers. For example, an ethane compound with the index number next to
        each atom::

                    H(4)    H(7)
                     |      |
             H(3) - C(0) - C(1) - H(6)
                     |      |
                    H(2)   H(5)

        If replace=True, then this function removes the hydrogen atoms that are
        occupying where the C-C bond should occur between monomers.
        It is required that you specify which atoms should be removed which is
        achieved by the `indices` parameter.

        In this example, you would remove H(2) and H(7) by indicating indices
        [2, 7]. The resulting structure of the polymer can vary wildly depending
        on your choice for `indices`, so you will have to test out different
        combinations to find the two that result in the desired structure.

        Parameters
        ----------
        compound : mbuild.Compound
            A compound of the individual monomer
        indices : list of int of length 2
            The particle indices of compound that represent the polymer
            bonding sites. You can specify the indices of particles that will
            be replaced by the polymer bond, or indices of particles that act
            as the bonding sites. See the 'replace' parameter notes.
        separation : float, units nm
            The bond length desired at the monomer-monomer bonding site.
            (separation / 2) is used to set the length of each port
        orientation : list of array-like, shape=(3,) of length 2,
            default=[None, None]
            Vector along which to orient the port
            If replace = True, and orientation = None,
            the orientation of the bond between the particle being
            removed and the anchor particle is used.
            Recommended behavior is to leave orientation set to None
            if you are using replace=True.
        replace : Bool, required, default=True
            If True, then the particles identified by bonding_indices
            will be removed and ports are added to the particles they
            were initially bonded to. Only use replace=True in the case
            that bonding_indices point to hydrogen atoms bonded to the
            desired monomer-monomer bonding site particles.
            If False, then the particles identified by bonding_indices
            will have ports added, and no particles are removed from
            the monomer compound.
        """
        comp = clone(compound)
        remove_hydrogens = []

        head = [p for p in comp.particles() if p.particle_tag == head_tag][0]
        head_hydrogens = [p for p in head.direct_bonds() if p.name == "H"]
        if len(head_hydrogens) != 0 and head_orientation is None:
            bond_vectors = [h.pos - head.pos for h in head_hydrogens[:bond_order]]
            head_orientation = np.sum(bond_vectors, axis=0)
            head_orientation /= np.linalg.norm(head_orientation)

        head_port = Port(
            anchor=head,
            orientation=head_orientation,
            separation=separation / 2,
        )
        comp.add(head_port, label="up")
        for p in head_hydrogens[:bond_order]:
            remove_hydrogens.append(p)

        tail = [p for p in comp.particles() if p.particle_tag == tail_tag][0]
        tail_hydrogens = [p for p in tail.direct_bonds() if p.name == "H"]
        if len(tail_hydrogens) != 0 and tail_orientation is None:
            bond_vectors = [h.pos - tail.pos for h in tail_hydrogens]
            tail_orientation = np.sum(bond_vectors, axis=0)
            tail_orientation /= np.linalg.norm(tail_orientation)

        tail_port = Port(
            anchor=tail,
            orientation=tail_orientation,
            separation=separation / 2,
        )
        comp.add(tail_port, label="down")
        for p in tail_hydrogens[:bond_order]:
            remove_hydrogens.append(p)

        for p in remove_hydrogens:
            comp.remove(p)

        self._monomers.append(comp)

    def add_end_groups(
        self,
        compound,
        bond_tag=None,
        separation=None,
        orientation=None,
        bond_order=1,
        label="head",
        duplicate=True,
    ):
        """Add an mBuild compound to self.end_groups.

        End groups will be used to cap the polymer. Call this function for each
        unique end group compound to be used in the polymer, or call it once
        with duplicate=True if the head and tail end groups are the same.

        Notes
        -----
        Refer to the docstring notes of the add_monomer() function for an
        explanation of the correct way to use the `replace` and `index`
        parameters.

        Parameters
        ----------
        compound : mbuild.Compound
            A compound of the end group structure
        index : int
            The particle index in compound that represents the bonding
            site between the end group and polymer.
            You can specify the index of a particle that will
            be replaced by the polymer bond or that acts as the bonding site.
            See the `replace` parameter notes.
        separation : float
            The bond length (units nm) desired between monomer and end-group.
        orientation : array-like, shape=(3,), default None
            Vector along which to orient the port
            If `replace=True` and `orientation=None`, the orientation of the
            bond between the particle being removed and the anchor particle is
            used.
            Recommended behavior is to leave orientation set to None if you
            are using `replace=True`.
        label : str, default 'head'
            Whether to add the end group to the 'head or 'tail' of the polymer.
            If `duplicate=True`, `label` is ignored.
        duplicate : Bool, default True
            If True, then `compound` is duplicated and added to
            `Polymer.end_groups` twice. Set to True if you want the same end
            group compound at the head and tail of the polymer. If that's the
            case, you only need to call `add_end_groups()` once.
            If False, `compound` is not duplicated and only one instance of the
            end group structure is added to `Polymer.end_groups`. You can call
            `add_end_groups()` a second time to add another end group.
        """
        comp = clone(compound)
        remove_hydrogens = []

        head = [p for p in comp.particles() if p.particle_tag == bond_tag][0]
        head_hydrogens = [p for p in head.direct_bonds() if p.name == "H"]
        if len(head_hydrogens) != 0 and orientation is None:
            bond_vectors = [h.pos - head.pos for h in head_hydrogens[:bond_order]]
            orientation = np.sum(bond_vectors, axis=0)
            orientation /= np.linalg.norm(orientation)

        head_port = Port(
            anchor=head,
            orientation=orientation,
            separation=separation / 2,
        )
        comp.add(head_port, label="up")
        for p in head_hydrogens[:bond_order]:
            remove_hydrogens.append(p)

        if duplicate:
            self._end_groups = [comp, clone(comp)]
            self._headtail = [separation / 2, separation / 2]
        else:
            if label.lower() == "head":
                self._end_groups[0] = comp
                self._headtail[0] = separation / 2
            elif label.lower() == "tail":
                self._end_groups[1] = comp
                self._headtail[1] = separation / 2
            else:
                raise ValueError("Label must be 'head' or 'tail'")

    def create_periodic_bond(self, axis="z"):
        """Align and bond the end points of a polymer along an axis.

        Parameters
        ----------
        axis : str, default="z"
            Axis along which to orient the polymer taken as the line connected the
            free ports of the end group. May be "x", "y", or "z".
        """
        if self.head_port is None or self.tail_port is None:
            raise ValueError(
                "Polymer head_port and/or tail_port could not be found. Be sure to initialize "
                "this object without end_groups and build with add_hydrogens=False."
            )

        if axis.lower() == "x":
            x_axis_transform(
                compound=self,
                new_origin=self.head_port.pos,
                point_on_x_axis=self.tail_port.pos,
            )
        elif axis.lower() == "y":
            y_axis_transform(
                compound=self,
                new_origin=self.head_port.pos,
                point_on_y_axis=self.tail_port.pos,
            )
        elif axis.lower() == "z":
            z_axis_transform(
                compound=self,
                new_origin=self.head_port.pos,
                point_on_z_axis=self.tail_port.pos,
            )
        else:
            raise ValueError("axis must be either: 'x', 'y', or 'z'")

        force_overlap(self, self.head_port, self.tail_port)
