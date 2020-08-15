__all__ = ['clone', 'Compound', 'Particle']

import os
import tempfile
import itertools
import ele
import numpy as np

from collections import OrderedDict, Iterable
from copy import deepcopy
from warnings import warn
from oset import oset as OrderedSet

from ele.element import Element

from mbuild import conversion
from mbuild.bond_graph import BondGraph
from mbuild.box import Box
from mbuild.exceptions import MBuildError
from mbuild.utils.decorators import deprecated
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.utils.io import run_from_ipython, import_
from mbuild.utils.jsutils import overwrite_nglview_default
from mbuild.coordinate_transform import _translate, _rotate


def clone(existing_compound, clone_of=None, root_container=None):
    """A faster alternative to deepcopying.

    Does not resolve circular dependencies. This should be safe provided
    you never try to add the top of a Compound hierarchy to a
    sub-Compound.

    Parameters
    ----------
    existing_compound : mb.Compound
        Existing Compound that will be copied

    Other Parameters
    ----------------
    clone_of : dict, optional
    root_container : mb.Compound, optional

    """
    if clone_of is None:
        clone_of = dict()

    newone = existing_compound._clone(clone_of=clone_of,
                                      root_container=root_container)
    existing_compound._clone_bonds(clone_of=clone_of)
    return newone


class Compound(object):
    """A building block in the mBuild hierarchy.

    Compound is the superclass of all composite building blocks in the mBuild
    hierarchy. That is, all composite building blocks must inherit from
    compound, either directly or indirectly. The design of Compound follows the
    Composite design pattern (Gamma, Erich; Richard Helm; Ralph Johnson; John
    M. Vlissides (1995). Design Patterns: Elements of Reusable Object-Oriented
    Software. Addison-Wesley. p. 395. ISBN 0-201-63361-2.), with Compound being
    the composite, and Particle playing the role of the primitive (leaf) part,
    where Particle is in fact simply an alias to the Compound class.

    Compound maintains a list of children (other Compounds contained within),
    and provides a means to tag the children with labels, so that the compounds
    can be easily looked up later. Labels may also point to objects outside the
    Compound's containment hierarchy. Compound has built-in support for copying
    and deepcopying Compound hierarchies, enumerating particles or bonds in the
    hierarchy, proximity based searches, visualization, I/O operations, and a
    number of other convenience methods.

    Parameters
    ----------
    subcompounds : mb.Compound or list of mb.Compound, optional, default=None
        One or more compounds to be added to self.
    name : str, optional, default=self.__class__.__name__
        The type of Compound.
    pos : np.ndarray, shape=(3,), dtype=float, optional, default=[0, 0, 0]
        The position of the Compound in Cartestian space
    charge : float, optional, default=0.0
        Currently not used. Likely removed in next release.
    periodicity : np.ndarray, shape=(3,), dtype=float, optional, default=[0, 0, 0]
        The periodic lengths of the Compound in the x, y and z directions.
        Defaults to zeros which is treated as non-periodic.
    port_particle : bool, optional, default=False
        Whether or not this Compound is part of a Port
    box : mb.Box, optional
        The simulation box containing the compound. Also accounts for the
        periodicity. Defaults to None which is treated as non-periodic.
    element: str, optional, default=None
        The one or two character element symbol

    Attributes
    ----------
    bond_graph : mb.BondGraph
        Graph-like object that stores bond information for this Compound
    children : OrderedSet
        Contains all children (other Compounds).
    labels : OrderedDict
        Labels to Compound/Atom mappings. These do not necessarily need not be
        in self.children.
    parent : mb.Compound
        The parent Compound that contains this part. Can be None if this
        compound is the root of the containment hierarchy.
    referrers : set
        Other compounds that reference this part with labels.
    rigid_id : int, default=None
        The ID of the rigid body that this Compound belongs to.  Only Particles
        (the bottom of the containment hierarchy) can have integer values for
        `rigid_id`. Compounds containing rigid particles will always have
        `rigid_id == None`. See also `contains_rigid`.
    boundingbox : mb.Box
        The bounds (xmin, xmax, ymin, ymax, zmin, zmax) of particles in Compound
    center
    contains_rigid
    max_rigid_id
    n_particles
    n_bonds
    root
    xyz
    xyz_with_ports

    """

    def __init__(self, subcompounds=None, name=None, pos=None, charge=0.0,
                 periodicity=None, box=None, element=None,
                 port_particle=False):
        super(Compound, self).__init__()

        if name:
            if not isinstance(name, str):
                raise ValueError(
                    'Compound.name should be a string. You passed '
                    '{}'.format(name))
            self.name = name
        else:
            self.name = self.__class__.__name__

        # A periodicity of zero in any direction is treated as non-periodic.
        if periodicity is None:
            self._periodicity = np.array([0.0, 0.0, 0.0])
        else:
            self._periodicity = np.asarray(periodicity)

        if pos is not None:
            self._pos = np.asarray(pos, dtype=float)
        else:
            self._pos = np.zeros(3)

        self.parent = None
        self.children = OrderedSet()
        self.labels = OrderedDict()
        self.referrers = set()

        self.bond_graph = None
        self.port_particle = port_particle

        self._rigid_id = None
        self._contains_rigid = False
        self._check_if_contains_rigid_bodies = False

        self.box = box
        self.element = element

        # self.add() must be called after labels and children are initialized.
        if subcompounds:
            if charge:
                raise MBuildError(
                    'Cannot set the charge of a Compound containing '
                    'subcompounds.')
            self.add(subcompounds)
            self._charge = 0.0
        else:
            self._charge = charge

    def particles(self, include_ports=False):
        """Return all Particles of the Compound.

        Parameters
        ----------
        include_ports : bool, optional, default=False
            Include port particles

        Yields
        -------
        mb.Compound
            The next Particle in the Compound

        """
        if not self.children:
            yield self
        else:
            for particle in self._particles(include_ports):
                yield particle

    def _particles(self, include_ports=False):
        """Return all Particles of the Compound. """
        for child in self.successors():
            if not child.children:
                if include_ports or not child.port_particle:
                    yield child

    def successors(self):
        """Yield Compounds below self in the hierarchy.

        Yields
        -------
        mb.Compound
            The next Particle below self in the hierarchy

        """
        if not self.children:
            return
        for part in self.children:
            # Parts local to the current Compound.
            yield part
            # Parts further down the hierarchy.
            for subpart in part.successors():
                yield subpart

    @property
    def n_particles(self):
        """Return the number of Particles in the Compound.

        Returns
        -------
        int
            The number of Particles in the Compound

        """
        if not self.children:
            return 1
        else:
            return self._n_particles(include_ports=False)

    def _n_particles(self, include_ports=False):
        """Return the number of Particles in the Compound. """
        return sum(1 for _ in self._particles(include_ports))

    def _contains_only_ports(self):
        for part in self.children:
            if not part.port_particle:
                return False
        return True

    def ancestors(self):
        """Generate all ancestors of the Compound recursively.

        Yields
        ------
        mb.Compound
            The next Compound above self in the hierarchy

        """
        if self.parent is not None:
            yield self.parent
            for ancestor in self.parent.ancestors():
                yield ancestor

    @property
    def root(self):
        """The Compound at the top of self's hierarchy.

        Returns
        -------
        mb.Compound
            The Compound at the top of self's hierarchy

        """
        parent = None
        for parent in self.ancestors():
            pass
        if parent is None:
            return self
        return parent

    def particles_by_name(self, name):
        """Return all Particles of the Compound with a specific name

        Parameters
        ----------
        name : str
            Only particles with this name are returned

        Yields
        ------
        mb.Compound
            The next Particle in the Compound with the user-specified name

        """
        for particle in self.particles():
            if particle.name == name:
                yield particle

    def particles_by_element(self, element):
        """Return all Particles of the Compound with a specific name

        Parameters
        ----------
        name : str or ele.Element
            element abbreviation or element 

        Yields
        ------
        mb.Compound
            The next Particle in the Compound with the user-specified element

        """
        if not isinstance(element, Element):
            element = ele.element_from_symbol(element)
        for particle in self.particles():
            if particle.element == element:
                yield particle

    @property
    def charge(self):
        return sum([particle._charge for particle in self.particles()])

    @charge.setter
    def charge(self, value):
        if self._contains_only_ports():
            self._charge = value
        else:
            raise AttributeError(
                "charge is immutable for Compounds that are "
                "not at the bottom of the containment hierarchy.")

    @property
    def rigid_id(self):
        return self._rigid_id

    @rigid_id.setter
    def rigid_id(self, value):
        if self._contains_only_ports():
            self._rigid_id = value
            for ancestor in self.ancestors():
                ancestor._check_if_contains_rigid_bodies = True
        else:
            raise AttributeError(
                "rigid_id is immutable for Compounds that are "
                "not at the bottom of the containment hierarchy.")

    @property
    def contains_rigid(self):
        """Returns True if the Compound contains rigid bodies

        If the Compound contains any particle with a rigid_id != None
        then contains_rigid will return True. If the Compound has no
        children (i.e. the Compound resides at the bottom of the containment
        hierarchy) then contains_rigid will return False.

        Returns
        -------
        bool
            True if the Compound contains any particle with a rigid_id != None

        Notes
        -----
        The private variable '_check_if_contains_rigid_bodies' is used to help
        cache the status of 'contains_rigid'. If '_check_if_contains_rigid_bodies'
        is False, then the rigid body containment of the Compound has not changed,
        and the particle tree is not traversed, boosting performance.

        """
        if self._check_if_contains_rigid_bodies:
            self._check_if_contains_rigid_bodies = False
            if any(particle.rigid_id is not None for particle in self._particles()):
                self._contains_rigid = True
            else:
                self._contains_rigid = False
        return self._contains_rigid

    @property
    def max_rigid_id(self):
        """Returns the maximum rigid body ID contained in the Compound.

        This is usually used by compound.root to determine the maximum
        rigid_id in the containment hierarchy.

        Returns
        -------
        int or None
            The maximum rigid body ID contained in the Compound. If no
            rigid body IDs are found, None is returned

        """
        try:
            return max([particle.rigid_id for particle in self.particles()
                        if particle.rigid_id is not None])
        except ValueError:
            return

    def rigid_particles(self, rigid_id=None):
        """Generate all particles in rigid bodies.

        If a rigid_id is specified, then this function will only yield particles
        with a matching rigid_id.

        Parameters
        ----------
        rigid_id : int, optional
            Include only particles with this rigid body ID

        Yields
        ------
        mb.Compound
            The next particle with a rigid_id that is not None, or the next
            particle with a matching rigid_id if specified

        """
        for particle in self.particles():
            if rigid_id is not None:
                if particle.rigid_id == rigid_id:
                    yield particle
            else:
                if particle.rigid_id is not None:
                    yield particle

    def label_rigid_bodies(self, discrete_bodies=None, rigid_particles=None):
        """Designate which Compounds should be treated as rigid bodies

        If no arguments are provided, this function will treat the compound
        as a single rigid body by providing all particles in `self` with the
        same rigid_id. If `discrete_bodies` is not None, each instance of
        a Compound with a name found in `discrete_bodies` will be treated as a
        unique rigid body. If `rigid_particles` is not None, only Particles
        (Compounds at the bottom of the containment hierarchy) matching this name
        will be considered part of the rigid body.

        Parameters
        ----------
        discrete_bodies : str or list of str, optional, default=None
            Name(s) of Compound instances to be treated as unique rigid bodies.
            Compound instances matching this (these) name(s) will be provided
            with unique rigid_ids
        rigid_particles : str or list of str, optional, default=None
            Name(s) of Compound instances at the bottom of the containment
            hierarchy (Particles) to be included in rigid bodies. Only Particles
            matching this (these) name(s) will have their rigid_ids altered to
            match the rigid body number.

        Examples
        --------
        Creating a rigid benzene

        >>> import mbuild as mb
        >>> from mbuild.utils.io import get_fn
        >>> benzene = mb.load(get_fn('benzene.mol2'))
        >>> benzene.label_rigid_bodies()

        Creating a semi-rigid benzene, where only the carbons are treated as
        a rigid body

        >>> import mbuild as mb
        >>> from mbuild.utils.io import get_fn
        >>> benzene = mb.load(get_fn('benzene.mol2'))
        >>> benzene.label_rigid_bodies(rigid_particles='C')

        Create a box of rigid benzenes, where each benzene has a unique rigid
        body ID.

        >>> import mbuild as mb
        >>> from mbuild.utils.io import get_fn
        >>> benzene = mb.load(get_fn('benzene.mol2'))
        >>> benzene.name = 'Benzene'
        >>> filled = mb.fill_box(benzene,
        ...                      n_compounds=10,
        ...                      box=[0, 0, 0, 4, 4, 4])
        >>> filled.label_rigid_bodies(distinct_bodies='Benzene')

        Create a box of semi-rigid benzenes, where each benzene has a unique
        rigid body ID and only the carbon portion is treated as rigid.

        >>> import mbuild as mb
        >>> from mbuild.utils.io import get_fn
        >>> benzene = mb.load(get_fn('benzene.mol2'))
        >>> benzene.name = 'Benzene'
        >>> filled = mb.fill_box(benzene,
        ...                      n_compounds=10,
        ...                      box=[0, 0, 0, 4, 4, 4])
        >>> filled.label_rigid_bodies(distinct_bodies='Benzene',
        ...                           rigid_particles='C')

        """
        if discrete_bodies is not None:
            if isinstance(discrete_bodies, str):
                discrete_bodies = [discrete_bodies]
        if rigid_particles is not None:
            if isinstance(rigid_particles, str):
                rigid_particles = [rigid_particles]

        if self.root.max_rigid_id is not None:
            rigid_id = self.root.max_rigid_id + 1
            warn("{} rigid bodies already exist.  Incrementing 'rigid_id'"
                 "starting from {}.".format(rigid_id, rigid_id))
        else:
            rigid_id = 0

        for successor in self.successors():
            if discrete_bodies and successor.name not in discrete_bodies:
                continue
            for particle in successor.particles():
                if rigid_particles and particle.name not in rigid_particles:
                    continue
                particle.rigid_id = rigid_id
            if discrete_bodies:
                rigid_id += 1

    def unlabel_rigid_bodies(self):
        """Remove all rigid body labels from the Compound """
        self._check_if_contains_rigid_bodies = True
        for child in self.children:
            child._check_if_contains_rigid_bodies = True
        for particle in self.particles():
            particle.rigid_id = None

    def _increment_rigid_ids(self, increment):
        """Increment the rigid_id of all rigid Particles in a Compound

        Adds `increment` to the rigid_id of all Particles in `self` that
        already have an integer rigid_id.
        """
        for particle in self.particles():
            if particle.rigid_id is not None:
                particle.rigid_id += increment

    def _reorder_rigid_ids(self):
        """Reorder rigid body IDs ensuring consecutiveness.

        Primarily used internally to ensure consecutive rigid_ids following
        removal of a Compound.

        """
        max_rigid = self.max_rigid_id
        unique_rigid_ids = sorted(
            set([p.rigid_id for p in self.rigid_particles()]))
        n_unique_rigid = len(unique_rigid_ids)
        if max_rigid and n_unique_rigid != max_rigid + 1:
            missing_rigid_id = (
                unique_rigid_ids[-1] * (unique_rigid_ids[-1] + 1)) / 2 - sum(unique_rigid_ids)
            for successor in self.successors():
                if successor.rigid_id is not None:
                    if successor.rigid_id > missing_rigid_id:
                        successor.rigid_id -= 1
            if self.rigid_id:
                if self.rigid_id > missing_rigid_id:
                    self.rigid_id -= 1

    def add(self, new_child, label=None, containment=True, replace=False,
            inherit_periodicity=True, inherit_box=False, reset_rigid_ids=True):
        """Add a part to the Compound.

        Note:
            This does not necessarily add the part to self.children but may
            instead be used to add a reference to the part to self.labels. See
            'containment' argument.

        Parameters
        ----------
        new_child : mb.Compound or list-like of mb.Compound
            The object(s) to be added to this Compound.
        label : str, optional
            A descriptive string for the part.
        containment : bool, optional, default=True
            Add the part to self.children.
        replace : bool, optional, default=True
            Replace the label if it already exists.
        inherit_periodicity : bool, optional, default=True
            Replace the periodicity of self with the periodicity of the
            Compound being added
        inherit_box: bool, optional, default=False
            Replace the box of self with the box of the Compound being added
        reset_rigid_ids : bool, optional, default=True
            If the Compound to be added contains rigid bodies, reset the
            rigid_ids such that values remain distinct from rigid_ids
            already present in `self`. Can be set to False if attempting
            to add Compounds to an existing rigid body.

        """
        # Support batch add via lists, tuples and sets.
        if (isinstance(new_child, Iterable) and
                not isinstance(new_child, str)):
            for child in new_child:
                self.add(child, reset_rigid_ids=reset_rigid_ids)
            return

        if not isinstance(new_child, Compound):
            raise ValueError('Only objects that inherit from mbuild.Compound '
                             'can be added to Compounds. You tried to add '
                             '"{}".'.format(new_child))

        if new_child.contains_rigid or new_child.rigid_id is not None:
            if self.contains_rigid and reset_rigid_ids:
                new_child._increment_rigid_ids(increment=self.max_rigid_id + 1)
            self._check_if_contains_rigid_bodies = True
        if self.rigid_id is not None:
            self.rigid_id = None

        # Create children and labels on the first add operation
        if self.children is None:
            self.children = OrderedSet()
        if self.labels is None:
            self.labels = OrderedDict()

        if containment:
            if new_child.parent is not None:
                raise MBuildError('Part {} already has a parent: {}'.format(
                    new_child, new_child.parent))
            self.children.add(new_child)
            new_child.parent = self

            if new_child.bond_graph is not None:
                if self.root.bond_graph is None:
                    self.root.bond_graph = new_child.bond_graph
                else:
                    self.root.bond_graph.compose(new_child.bond_graph)

                new_child.bond_graph = None

        # Add new_part to labels. Does not currently support batch add.
        if label is None:
            label = '{0}[$]'.format(new_child.__class__.__name__)

        if label.endswith('[$]'):
            label = label[:-3]
            if label not in self.labels:
                self.labels[label] = []
            label_pattern = label + '[{}]'

            count = len(self.labels[label])
            self.labels[label].append(new_child)
            label = label_pattern.format(count)

        if not replace and label in self.labels:
            raise MBuildError('Label "{0}" already exists in {1}.'.format(
                label, self))
        else:
            self.labels[label] = new_child
        new_child.referrers.add(self)

        if (inherit_periodicity and isinstance(new_child, Compound) and
                new_child.periodicity.any()):
            self.periodicity = new_child.periodicity

        # If parent has no box --> inherit child box
        # If parent has box --> keep unless inherit_box == True
        # If inherit_box == True, parent box != None, child_box == None,
        # keep parent box anyway and warn
        # If inherit_box == False, child.box != None, warn
        if self.box is None:
            if new_child.box is not None:
                self.box = new_child.box
        else:
            if inherit_box:
                if new_child.box is None:
                    warn(
                        "The Compound you are adding has no box but "
                        "inherit_box=True. The box of the original "
                        "Compound will remain unchanged."
                    )
                else:
                    self.box = new_child.box
            else:
                if new_child.box is not None:
                    warn(
                        "The Compound you are adding has a box. "
                        "The box of the parent compound will be used. Use "
                        "inherit_box = True if you wish to replace the parent "
                        "compound box with that of Compound being added."
                    )

        # Check that bounding box is within box after adding compound
        if self.box:
            if (self.box.lengths < self.boundingbox.lengths).any():
                warn(
                    "After adding new Compound, Compound.box.lengths < "
                    "Compound.boundingbox.lengths. There may be particles "
                    "outside of the defined simulation box"
                )


    def remove(self, objs_to_remove):
        """ Cleanly remove children from the Compound.

        Parameters
        ----------
        objs_to_remove : mb.Compound or list of mb.Compound
            The Compound(s) to be removed from self

        """
        # Preprocessing and validating input type
        from mbuild.port import Port
        if not hasattr(objs_to_remove, '__iter__'):
            objs_to_remove = [objs_to_remove]
        objs_to_remove = set(objs_to_remove)

        # If nothing is to be remove, do nothing
        if len(objs_to_remove) == 0:
            return

        # Remove Port objects separately
        ports_removed = set()
        for obj in objs_to_remove:
            if isinstance(obj, Port):
                ports_removed.add(obj)
                self._remove(obj)
                obj.parent.children.remove(obj)
                self._remove_references(obj)

        objs_to_remove = objs_to_remove - ports_removed

        # Get particles to remove
        particles_to_remove = set([particle for obj in objs_to_remove
                                            for particle in obj.particles()])

        # Recursively get container compounds to remove
        to_remove = list()

        def _check_if_empty(child):
            if child in to_remove:
                return
            if set(child.particles()).issubset(particles_to_remove):
                if child.parent:
                    to_remove.append(child)
                    _check_if_empty(child.parent)
                else:
                    warn("This will remove all particles in "
                            "compound {}".format(self))
            return

        for particle in particles_to_remove:
            _check_if_empty(particle)

        # Fix rigid_ids and remove obj from bondgraph
        for removed_part in to_remove:
            self._remove(removed_part)

        # Remove references to object
        for removed_part in to_remove:
            if removed_part.parent is not None:
                removed_part.parent.children.remove(removed_part)
            self._remove_references(removed_part)

        # Remove ghost ports
        all_ports_list = list(self.all_ports())
        for port in all_ports_list:
            if port.anchor not in [i for i in self.particles()]:
                port.parent.children.remove(port)

        # Check and reorder rigid id
        for _ in particles_to_remove:
            if self.contains_rigid:
                self.root._reorder_rigid_ids()


    def _remove(self, removed_part):
        """Worker for remove(). Fixes rigid IDs and removes bonds"""
        if removed_part.rigid_id is not None:
            for ancestor in removed_part.ancestors():
                ancestor._check_if_contains_rigid_bodies = True
        if self.root.bond_graph and self.root.bond_graph.has_node(
                removed_part):
            for neighbor in self.root.bond_graph.neighbors(
                    removed_part):
                self.root.remove_bond((removed_part, neighbor))
            self.root.bond_graph.remove_node(removed_part)


    def _remove_references(self, removed_part):
        """Remove labels pointing to this part and vice versa. """
        removed_part.parent = None

        # Remove labels in the hierarchy pointing to this part.
        referrers_to_remove = set()
        for referrer in removed_part.referrers:
            if removed_part not in referrer.ancestors():
                for label, referred_part in list(referrer.labels.items()):
                    if referred_part is removed_part:
                        del referrer.labels[label]
                        referrers_to_remove.add(referrer)
        removed_part.referrers -= referrers_to_remove

        # Remove labels in this part pointing into the hierarchy.
        labels_to_delete = []
        if isinstance(removed_part, Compound):
            for label, part in list(removed_part.labels.items()):
                if not isinstance(part, Compound):
                    for p in part:
                        self._remove_references(p)
                elif removed_part not in part.ancestors():
                    try:
                        part.referrers.discard(removed_part)
                    except KeyError:
                        pass
                    else:
                        labels_to_delete.append(label)
        for label in labels_to_delete:
            removed_part.labels.pop(label, None)

    def referenced_ports(self):
        """Return all Ports referenced by this Compound.

        Returns
        -------
        list of mb.Compound
            A list of all ports referenced by the Compound

        """
        from mbuild.port import Port
        return [port for port in self.labels.values()
                if isinstance(port, Port)]

    def all_ports(self):
        """Return all Ports referenced by this Compound and its successors

        Returns
        -------
        list of mb.Compound
            A list of all Ports referenced by this Compound and its successors

        """
        from mbuild.port import Port
        return [successor for successor in self.successors()
                if isinstance(successor, Port)]

    def available_ports(self):
        """Return all unoccupied Ports referenced by this Compound.

        Returns
        -------
        list of mb.Compound
            A list of all unoccupied ports referenced by the Compound

        """
        from mbuild.port import Port
        return [port for port in self.labels.values()
                if isinstance(port, Port) and not port.used]

    def bonds(self):
        """Return all bonds in the Compound and sub-Compounds.

        Yields
        -------
        tuple of mb.Compound
            The next bond in the Compound

        See Also
        --------
        bond_graph.edges_iter : Iterates over all edges in a BondGraph

        """
        if self.root.bond_graph:
            if self.root == self:
                return self.root.bond_graph.edges_iter()
            else:
                return self.root.bond_graph.subgraph(
                    self.particles()).edges_iter()
        else:
            return iter(())

    @property
    def n_bonds(self):
        """Return the number of bonds in the Compound.

        Returns
        -------
        int
            The number of bonds in the Compound

        """
        return sum(1 for _ in self.bonds())

    def add_bond(self, particle_pair):
        """Add a bond between two Particles.

        Parameters
        ----------
        particle_pair : indexable object, length=2, dtype=mb.Compound
            The pair of Particles to add a bond between

        """
        if self.root.bond_graph is None:
            self.root.bond_graph = BondGraph()

        self.root.bond_graph.add_edge(particle_pair[0], particle_pair[1])

    def generate_bonds(self, name_a, name_b, dmin, dmax):
        """Add Bonds between all pairs of types a/b within [dmin, dmax].

        Parameters
        ----------
        name_a : str
            The name of one of the Particles to be in each bond
        name_b : str
            The name of the other Particle to be in each bond
        dmin : float
            The minimum distance between Particles for considering a bond
        dmax : float
            The maximum distance between Particles for considering a bond

        """
        particle_kdtree = PeriodicCKDTree(
            data=self.xyz, bounds=self.periodicity)
        particle_array = np.array(list(self.particles()))
        added_bonds = list()
        for p1 in self.particles_by_name(name_a):
            nearest = self.particles_in_range(p1, dmax, max_particles=20,
                                              particle_kdtree=particle_kdtree,
                                              particle_array=particle_array)
            for p2 in nearest:
                if p2 == p1:
                    continue
                bond_tuple = (p1, p2) if id(p1) < id(p2) else (p2, p1)
                if bond_tuple in added_bonds:
                    continue
                min_dist = self.min_periodic_distance(p2.pos, p1.pos)
                if (p2.name == name_b) and (dmin <= min_dist <= dmax):
                    self.add_bond((p1, p2))
                    added_bonds.append(bond_tuple)

    def remove_bond(self, particle_pair):
        """Deletes a bond between a pair of Particles

        Parameters
        ----------
        particle_pair : indexable object, length=2, dtype=mb.Compound
            The pair of Particles to remove the bond between

        """
        from mbuild.port import Port
        if self.root.bond_graph is None or not self.root.bond_graph.has_edge(
                *particle_pair):
            warn("Bond between {} and {} doesn't exist!".format(*particle_pair))
            return
        self.root.bond_graph.remove_edge(*particle_pair)
        bond_vector = particle_pair[0].pos - particle_pair[1].pos
        if np.allclose(bond_vector, np.zeros(3)):
            warn("Particles {} and {} overlap! Ports will not be added."
                 "".format(*particle_pair))
            return
        distance = np.linalg.norm(bond_vector)
        particle_pair[0].parent.add(Port(anchor=particle_pair[0],
                                         orientation=-bond_vector,
                                         separation=distance / 2), 'port[$]')
        particle_pair[1].parent.add(Port(anchor=particle_pair[1],
                                         orientation=bond_vector,
                                         separation=distance / 2), 'port[$]')

    @property
    def pos(self):
        if not self.children:
            return self._pos
        else:
            return self.center

    @pos.setter
    def pos(self, value):
        if not self.children:
            self._pos = value
        else:
            raise MBuildError('Cannot set position on a Compound that has'
                              ' children.')

    @property
    def periodicity(self):
        return self._periodicity

    @periodicity.setter
    def periodicity(self, periods):
        self._periodicity = np.array(periods)

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self, box):
        if box is not None and type(box) != Box:
            raise TypeError("box must be specified as an mbuild.Box")
        if self.port_particle and box is not None:
            raise ValueError("Ports cannot have a box")
        # TODO: Fix this for non-orthogonal boxes
        # Make sure the box is bigger than the bounding box
        if box is not None:
            if (box.lengths < self.boundingbox.lengths).any():
                warn(
                    "Compound.box.lengths < Compound.boundingbox.lengths. "
                    "There may be particles outside of the defined "
                    "simulation box."
                )
        self._box = box

    @property
    def element(self):
        return self._element
    
    @element.setter
    def element(self, element):
        if element is None:
            self._element = None
        elif isinstance(element, Element):
            self._element = element
        else:
            self._element = ele.element_from_symbol(element)

    @property
    def xyz(self):
        """Return all particle coordinates in this compound.

        Returns
        -------
        pos : np.ndarray, shape=(n, 3), dtype=float
            Array with the positions of all particles.
        """
        if not self.children:
            pos = np.expand_dims(self._pos, axis=0)
        else:
            arr = np.fromiter(itertools.chain.from_iterable(
                particle.pos for particle in self.particles()), dtype=float)
            pos = arr.reshape((-1, 3))
        return pos

    @property
    def xyz_with_ports(self):
        """Return all particle coordinates in this compound including ports.

        Returns
        -------
        pos : np.ndarray, shape=(n, 3), dtype=float
            Array with the positions of all particles and ports.

        """
        if not self.children:
            pos = self._pos
        else:
            arr = np.fromiter(
                itertools.chain.from_iterable(
                    particle.pos for particle in self.particles(
                        include_ports=True)), dtype=float)
            pos = arr.reshape((-1, 3))
        return pos

    @xyz.setter
    def xyz(self, arrnx3):
        """Set the positions of the particles in the Compound, excluding the Ports.

        This function does not set the position of the ports.

        Parameters
        ----------
        arrnx3 : np.ndarray, shape=(n,3), dtype=float
            The new particle positions

        """
        if not self.children:
            if not arrnx3.shape[0] == 1:
                raise ValueError(
                    'Trying to set position of {} with more than one'
                    'coordinate: {}'.format(
                        self, arrnx3))
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(
                self._particles(
                    include_ports=False), arrnx3):
                atom.pos = coords

    @xyz_with_ports.setter
    def xyz_with_ports(self, arrnx3):
        """Set the positions of the particles in the Compound, including the Ports.

        Parameters
        ----------
        arrnx3 : np.ndarray, shape=(n,3), dtype=float
            The new particle positions

        """
        if not self.children:
            if not arrnx3.shape[0] == 1:
                raise ValueError(
                    'Trying to set position of {} with more than one'
                    'coordinate: {}'.format(
                        self, arrnx3))
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(
                self._particles(
                    include_ports=True), arrnx3):
                atom.pos = coords

    @property
    def center(self):
        """The cartesian center of the Compound based on its Particles.

        Returns
        -------
        np.ndarray, shape=(3,), dtype=float
            The cartesian center of the Compound based on its Particles

        """

        if np.all(np.isfinite(self.xyz)):
            return np.mean(self.xyz, axis=0)

    @property
    def boundingbox(self):
        """Compute the bounding box of the compound.

        Returns
        -------
        mb.Box
            The bounding box for this Compound

        """
        xyz = self.xyz
        return Box(mins=xyz.min(axis=0), maxs=xyz.max(axis=0))

    def min_periodic_distance(self, xyz0, xyz1):
        """Vectorized distance calculation considering minimum image.

        Parameters
        ----------
        xyz0 : np.ndarray, shape=(3,), dtype=float
            Coordinates of first point
        xyz1 : np.ndarray, shape=(3,), dtype=float
            Coordinates of second point

        Returns
        -------
        float
            Vectorized distance between the two points following minimum
            image convention

        """
        d = np.abs(xyz0 - xyz1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def particles_in_range(
            self,
            compound,
            dmax,
            max_particles=20,
            particle_kdtree=None,
            particle_array=None):
        """Find particles within a specified range of another particle.

        Parameters
        ----------
        compound : mb.Compound
            Reference particle to find other particles in range of
        dmax : float
            Maximum distance from 'compound' to look for Particles
        max_particles : int, optional, default=20
            Maximum number of Particles to return
        particle_kdtree : mb.PeriodicCKDTree, optional
            KD-tree for looking up nearest neighbors. If not provided, a KD-
            tree will be generated from all Particles in self
        particle_array : np.ndarray, shape=(n,), dtype=mb.Compound, optional
            Array of possible particles to consider for return. If not
            provided, this defaults to all Particles in self

        Returns
        -------
        np.ndarray, shape=(n,), dtype=mb.Compound
            Particles in range of compound according to user-defined limits

        See Also
        --------
        periodic_kdtree.PerioidicCKDTree : mBuild implementation of kd-trees
        scipy.spatial.ckdtree : Further details on kd-trees

        """
        if particle_kdtree is None:
            particle_kdtree = PeriodicCKDTree(
                data=self.xyz, bounds=self.periodicity)
        _, idxs = particle_kdtree.query(
            compound.pos, k=max_particles, distance_upper_bound=dmax)
        idxs = idxs[idxs != self.n_particles]
        if particle_array is None:
            particle_array = np.array(list(self.particles()))
        return particle_array[idxs]

    def visualize(self, show_ports=False,
            backend='py3dmol', color_scheme={}): # pragma: no cover
        """Visualize the Compound using py3dmol (default) or nglview.

        Allows for visualization of a Compound within a Jupyter Notebook.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        backend : str, optional, default='py3dmol'
            Specify the backend package to visualize compounds
            Currently supported: py3dmol, nglview
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}

        """
        viz_pkg = {'nglview': self._visualize_nglview,
                'py3dmol': self._visualize_py3dmol}
        if run_from_ipython():
            if backend.lower() in viz_pkg:
                return viz_pkg[backend.lower()](show_ports=show_ports,
                        color_scheme=color_scheme)
            else:
                raise RuntimeError("Unsupported visualization " +
                        "backend ({}). ".format(backend) +
                        "Currently supported backends include nglview and py3dmol")

        else:
            raise RuntimeError('Visualization is only supported in Jupyter '
                               'Notebooks.')

    def _visualize_py3dmol(self, show_ports=False, color_scheme={}):
        """Visualize the Compound using py3Dmol.

        Allows for visualization of a Compound within a Jupyter Notebook.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        color_scheme : dict, optional
            Specify coloring for non-elemental particles
            keys are strings of the particle names
            values are strings of the colors
            i.e. {'_CGBEAD': 'blue'}


        Returns
        ------
        view : py3Dmol.view

        """
        py3Dmol = import_('py3Dmol')

        cloned = clone(self)

        modified_color_scheme = {}
        for name, color in color_scheme.items():
            # Py3dmol does some element string conversions,
            # first character is as-is, rest of the characters are lowercase
            new_name = name[0] + name[1:].lower()
            modified_color_scheme[new_name] = color
            modified_color_scheme[name] = color

        for particle in cloned.particles():
            if not particle.name:
                particle.name = 'UNK'
        tmp_dir = tempfile.mkdtemp()
        cloned.save(os.path.join(tmp_dir, 'tmp.mol2'),
                  show_ports=show_ports,
                  overwrite=True)

        view = py3Dmol.view()
        with open(os.path.join(tmp_dir, 'tmp.mol2'), 'r') as f:
            view.addModel(f.read(), 'mol2', keepH=True)

        view.setStyle({'stick': {'radius': 0.2,
                                'color':'grey'},
                        'sphere': {'scale': 0.3,
                                    'colorscheme':modified_color_scheme}})
        view.zoomTo()

        return view

    def _visualize_nglview(self, show_ports=False, color_scheme={}):
        """Visualize the Compound using nglview.

        Allows for visualization of a Compound within a Jupyter Notebook.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
            """
        nglview = import_('nglview')
        mdtraj = import_('mdtraj')
        from mdtraj.geometry.sasa import _ATOMIC_RADII
        remove_digits = lambda x: ''.join(i for i in x if not i.isdigit()
                                              or i == '_')
        for particle in self.particles():
            particle.name = remove_digits(particle.name).upper()
            if not particle.name:
                particle.name = 'UNK'
        tmp_dir = tempfile.mkdtemp()
        self.save(os.path.join(tmp_dir, 'tmp.mol2'),
                  show_ports=show_ports,
                  overwrite=True)
        widget = nglview.show_file(os.path.join(tmp_dir, 'tmp.mol2'))
        widget.clear()
        widget.add_ball_and_stick(cylinderOnly=True)
        elements = set([particle.name for particle in self.particles()])
        scale = 50.0
        for element in elements:
            try:
                widget.add_ball_and_stick('_{}'.format(
                    element.upper()), aspect_ratio=_ATOMIC_RADII[element.title()]**1.5 * scale)
            except KeyError:
                ids = [str(i) for i, particle in enumerate(self.particles())
                       if particle.name == element]
                widget.add_ball_and_stick(
                    '@{}'.format(
                        ','.join(ids)),
                    aspect_ratio=0.17**1.5 * scale,
                    color='grey')
        if show_ports:
            widget.add_ball_and_stick('_VS',
                                      aspect_ratio=1.0, color='#991f00')
        overwrite_nglview_default(widget)
        return widget

    def update_coordinates(self, filename, update_port_locations=True):
        """Update the coordinates of this Compound from a file.

        Parameters
        ----------
        filename : str
            Name of file from which to load coordinates. Supported file types
            are the same as those supported by load()
        update_port_locations : bool, optional, default=True
            Update the locations of Ports so that they are shifted along with
            their anchor particles.  Note: This conserves the location of
            Ports with respect to the anchor Particle, but does not conserve
            the orientation of Ports with respect to the molecule as a whole.

        See Also
        --------
        load : Load coordinates from a file

        """
        if update_port_locations:
            xyz_init = self.xyz
            self = conversion.load(filename, compound=self, coords_only=True)
            self._update_port_locations(xyz_init)
        else:
            self = conversion.load(filename, compound=self, coords_only=True)

    def _update_port_locations(self, initial_coordinates):
        """Adjust port locations after particles have moved

        Compares the locations of Particles between 'self' and an array of
        reference coordinates.  Shifts Ports in accordance with how far anchors
        have been moved.  This conserves the location of Ports with respect to
        their anchor Particles, but does not conserve the orientation of Ports
        with respect to the molecule as a whole.

        Parameters
        ----------
        initial_coordinates : np.ndarray, shape=(n, 3), dtype=float
            Reference coordinates to use for comparing how far anchor Particles
            have shifted.

        """
        particles = list(self.particles())
        for port in self.all_ports():
            if port.anchor:
                idx = particles.index(port.anchor)
                shift = particles[idx].pos - initial_coordinates[idx]
                port.translate(shift)

    def _kick(self):
        """Slightly adjust all coordinates in a Compound

        Provides a slight adjustment to coordinates to kick them out of local
        energy minima.
        """
        xyz_init = self.xyz
        for particle in self.particles():
            particle.pos += (np.random.rand(3,) - 0.5) / 100
        self._update_port_locations(xyz_init)

    warning_message = 'Please use Compound.energy_minimize()'

    @deprecated(warning_message)
    def energy_minimization(self, forcefield='UFF', steps=1000, **kwargs):
        self.energy_minimize(forcefield=forcefield, steps=steps, **kwargs)

    def energy_minimize(self, forcefield='UFF', steps=1000, **kwargs):
        """Perform an energy minimization on a Compound

        Default behavior utilizes Open Babel (http://openbabel.org/docs/dev/)
        to perform an energy minimization/geometry optimization on a
        Compound by applying a generic force field

        Can also utilize OpenMM (http://openmm.org/) to energy minimize
        after atomtyping a Compound using
        Foyer (https://github.com/mosdef-hub/foyer) to apply a forcefield
        XML file that contains valid SMARTS strings.

        This function is primarily intended to be used on smaller components,
        with sizes on the order of 10's to 100's of particles, as the energy
        minimization scales poorly with the number of particles.

        Parameters
        ----------
        steps : int, optional, default=1000
            The number of optimization iterations
        forcefield : str, optional, default='UFF'
            The generic force field to apply to the Compound for minimization.
            Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', and 'Ghemical'.
            Please refer to the Open Babel documentation (http://open-babel.
            readthedocs.io/en/latest/Forcefields/Overview.html) when considering
            your choice of force field.
            Utilizing OpenMM for energy minimization requires a forcefield
            XML file with valid SMARTS strings. Please refer to (http://docs.
            openmm.org/7.0.0/userguide/application.html#creating-force-fields)
            for more information.


        Keyword Arguments
        ------------
        algorithm : str, optional, default='cg'
            The energy minimization algorithm.  Valid options are 'steep',
            'cg', and 'md', corresponding to steepest descent, conjugate
            gradient, and equilibrium molecular dynamics respectively.
            For _energy_minimize_openbabel
        scale_bonds : float, optional, default=1
            Scales the bond force constant (1 is completely on).
            For _energy_minimize_openmm
        scale_angles : float, optional, default=1
            Scales the angle force constant (1 is completely on)
            For _energy_minimize_openmm
        scale_torsions : float, optional, default=1
            Scales the torsional force constants (1 is completely on)
            For _energy_minimize_openmm
            Note: Only Ryckaert-Bellemans style torsions are currently supported
        scale_nonbonded : float, optional, default=1
            Scales epsilon (1 is completely on)
            For _energy_minimize_openmm

        References
        ----------
        If using _energy_minimize_openmm(), please cite:
        .. [1] P. Eastman, M. S. Friedrichs, J. D. Chodera, R. J. Radmer,
               C. M. Bruns, J. P. Ku, K. A. Beauchamp, T. J. Lane,
               L.-P. Wang, D. Shukla, T. Tye, M. Houston, T. Stich,
               C. Klein, M. R. Shirts, and V. S. Pande.
               "OpenMM 4: A Reusable, Extensible, Hardware Independent
               Library for High Performance Molecular Simulation."
               J. Chem. Theor. Comput. 9(1): 461-469. (2013).


        If using _energy_minimize_openbabel(), please cite:
        .. [1] O'Boyle, N.M.; Banck, M.; James, C.A.; Morley, C.;
               Vandermeersch, T.; Hutchison, G.R. "Open Babel: An open
               chemical toolbox." (2011) J. Cheminf. 3, 33

        .. [2] Open Babel, version X.X.X http://openbabel.org, (installed
               Month Year)

        If using the 'MMFF94' force field please also cite the following:
        .. [3] T.A. Halgren, "Merck molecular force field. I. Basis, form,
               scope, parameterization, and performance of MMFF94." (1996)
               J. Comput. Chem. 17, 490-519
        .. [4] T.A. Halgren, "Merck molecular force field. II. MMFF94 van der
               Waals and electrostatic parameters for intermolecular
               interactions." (1996) J. Comput. Chem. 17, 520-552
        .. [5] T.A. Halgren, "Merck molecular force field. III. Molecular
               geometries and vibrational frequencies for MMFF94." (1996)
               J. Comput. Chem. 17, 553-586
        .. [6] T.A. Halgren and R.B. Nachbar, "Merck molecular force field.
               IV. Conformational energies and geometries for MMFF94." (1996)
               J. Comput. Chem. 17, 587-615
        .. [7] T.A. Halgren, "Merck molecular force field. V. Extension of
               MMFF94 using experimental data, additional computational data,
               and empirical rules." (1996) J. Comput. Chem. 17, 616-641

        If using the 'MMFF94s' force field please cite the above along with:
        .. [8] T.A. Halgren, "MMFF VI. MMFF94s option for energy minimization
               studies." (1999) J. Comput. Chem. 20, 720-729

        If using the 'UFF' force field please cite the following:
        .. [3] Rappe, A.K., Casewit, C.J., Colwell, K.S., Goddard, W.A. III,
               Skiff, W.M. "UFF, a full periodic table force field for
               molecular mechanics and molecular dynamics simulations." (1992)
               J. Am. Chem. Soc. 114, 10024-10039

        If using the 'GAFF' force field please cite the following:
        .. [3] Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A., Case, D.A.
               "Development and testing of a general AMBER force field" (2004)
               J. Comput. Chem. 25, 1157-1174

        If using the 'Ghemical' force field please cite the following:
        .. [3] T. Hassinen and M. Perakyla, "New energy terms for reduced
               protein models implemented in an off-lattice force field" (2001)
               J. Comput. Chem. 22, 1229-1242



        """
        tmp_dir = tempfile.mkdtemp()
        original = clone(self)
        self._kick()
        self.save(os.path.join(tmp_dir, 'un-minimized.mol2'))
        extension = os.path.splitext(forcefield)[-1]
        openbabel_ffs = ['MMFF94', 'MMFF94s', 'UFF', 'GAFF', 'Ghemical']
        if forcefield in openbabel_ffs:
            self._energy_minimize_openbabel(tmp_dir, forcefield=forcefield,
                                            steps=steps, **kwargs)
        elif extension == '.xml':
            self._energy_minimize_openmm(tmp_dir, forcefield_files=forcefield,
                                         forcefield_name=None,
                                         steps=steps, **kwargs)
        else:
            self._energy_minimize_openmm(tmp_dir, forcefield_files=None,
                                         forcefield_name=forcefield,
                                         steps=steps, **kwargs)

        self.update_coordinates(os.path.join(tmp_dir, 'minimized.pdb'))

    def _energy_minimize_openmm(
            self,
            tmp_dir,
            forcefield_files=None,
            forcefield_name=None,
            steps=1000,
            scale_bonds=1,
            scale_angles=1,
            scale_torsions=1,
            scale_nonbonded=1):
        """ Perform energy minimization using OpenMM

        Converts an mBuild Compound to a ParmEd Structure,
        applies a forcefield using Foyer, and creates an OpenMM System.

        Parameters
        ----------
        forcefield_files : str or list of str, optional, default=None
            Forcefield files to load
        forcefield_name : str, optional, default=None
            Apply a named forcefield to the output file using the `foyer`
            package, e.g. 'oplsaa'. Forcefields listed here:
            https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields
        steps : int, optional, default=1000
            Number of energy minimization iterations
        scale_bonds : float, optional, default=1
            Scales the bond force constant (1 is completely on)
        scale_angles : float, optiona, default=1
            Scales the angle force constant (1 is completely on)
        scale_torsions : float, optional, default=1
            Scales the torsional force constants (1 is completely on)
        scale_nonbonded : float, optional, default=1
            Scales epsilon (1 is completely on)


        Notes
        -----
        Assumes a particular organization for the force groups
        (HarmonicBondForce, HarmonicAngleForce, RBTorsionForce, NonBondedForce)

        References
        ----------

        .. [1] P. Eastman, M. S. Friedrichs, J. D. Chodera, R. J. Radmer,
               C. M. Bruns, J. P. Ku, K. A. Beauchamp, T. J. Lane,
               L.-P. Wang, D. Shukla, T. Tye, M. Houston, T. Stich,
               C. Klein, M. R. Shirts, and V. S. Pande.
               "OpenMM 4: A Reusable, Extensible, Hardware Independent
               Library for High Performance Molecular Simulation."
               J. Chem. Theor. Comput. 9(1): 461-469. (2013).



        """
        foyer = import_('foyer')

        to_parmed = self.to_parmed()
        ff = foyer.Forcefield(forcefield_files=forcefield_files, name=forcefield_name)
        to_parmed = ff.apply(to_parmed)

        from simtk.openmm.app.simulation import Simulation
        from simtk.openmm.app.pdbreporter import PDBReporter
        from simtk.openmm.openmm import LangevinIntegrator
        import simtk.unit as u

        system = to_parmed.createSystem() # Create an OpenMM System
        # Create a Langenvin Integrator in OpenMM
        integrator = LangevinIntegrator(298 * u.kelvin, 1 / u.picosecond,
                                        0.002 * u.picoseconds)
        # Create Simulation object in OpenMM
        simulation = Simulation(to_parmed.topology, system, integrator)

        # Loop through forces in OpenMM System and set parameters
        for force in system.getForces():
            if type(force).__name__ == "HarmonicBondForce":
                for bond_index in range(force.getNumBonds()):
                    atom1, atom2, r0, k = force.getBondParameters(bond_index)
                    force.setBondParameters(bond_index,
                                            atom1, atom2,
                                            r0, k * scale_bonds)
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "HarmonicAngleForce":
                for angle_index in range(force.getNumAngles()):
                    atom1, atom2, atom3, r0, k = force.getAngleParameters(
                        angle_index)
                    force.setAngleParameters(angle_index,
                                             atom1, atom2, atom3,
                                             r0, k * scale_angles)
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "RBTorsionForce":
                for torsion_index in range(force.getNumTorsions()):
                    atom1, atom2, atom3, atom4, c0, c1, c2, c3, c4, c5 = force.getTorsionParameters(
                        torsion_index)
                    force.setTorsionParameters(
                        torsion_index,
                        atom1,
                        atom2,
                        atom3,
                        atom4,
                        c0 * scale_torsions,
                        c1 * scale_torsions,
                        c2 * scale_torsions,
                        c3 * scale_torsions,
                        c4 * scale_torsions,
                        c5 * scale_torsions)
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "NonbondedForce":
                for nb_index in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(
                        nb_index)
                    force.setParticleParameters(nb_index,
                                                charge, sigma,
                                                epsilon * scale_nonbonded)
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "CMMotionRemover":
                pass

            else:
                warn(
                    'OpenMM Force {} is '
                    'not currently supported in _energy_minimize_openmm. '
                    'This Force will not be updated!'.format(
                        type(force).__name__))

        simulation.context.setPositions(to_parmed.positions)
        # Run energy minimization through OpenMM
        simulation.minimizeEnergy(maxIterations=steps)
        reporter = PDBReporter(os.path.join(tmp_dir, 'minimized.pdb'), 1)
        reporter.report(
            simulation,
            simulation.context.getState(
                getPositions=True))

    def _energy_minimize_openbabel(self, tmp_dir, steps=1000, algorithm='cg',
                                   forcefield='UFF'):
        """Perform an energy minimization on a Compound

        Utilizes Open Babel (http://openbabel.org/docs/dev/) to perform an
        energy minimization/geometry optimization on a Compound by applying
        a generic force field.

        This function is primarily intended to be used on smaller components,
        with sizes on the order of 10's to 100's of particles, as the energy
        minimization scales poorly with the number of particles.

        Parameters
        ----------
        steps : int, optionl, default=1000
            The number of optimization iterations
        algorithm : str, optional, default='cg'
            The energy minimization algorithm.  Valid options are 'steep',
            'cg', and 'md', corresponding to steepest descent, conjugate
            gradient, and equilibrium molecular dynamics respectively.
        forcefield : str, optional, default='UFF'
            The generic force field to apply to the Compound for minimization.
            Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', and 'Ghemical'.
            Please refer to the Open Babel documentation (http://open-babel.
            readthedocs.io/en/latest/Forcefields/Overview.html) when considering
            your choice of force field.

        References
        ----------
        .. [1] O'Boyle, N.M.; Banck, M.; James, C.A.; Morley, C.;
               Vandermeersch, T.; Hutchison, G.R. "Open Babel: An open
               chemical toolbox." (2011) J. Cheminf. 3, 33
        .. [2] Open Babel, version X.X.X http://openbabel.org, (installed
               Month Year)

        If using the 'MMFF94' force field please also cite the following:
        .. [3] T.A. Halgren, "Merck molecular force field. I. Basis, form,
               scope, parameterization, and performance of MMFF94." (1996)
               J. Comput. Chem. 17, 490-519
        .. [4] T.A. Halgren, "Merck molecular force field. II. MMFF94 van der
               Waals and electrostatic parameters for intermolecular
               interactions." (1996) J. Comput. Chem. 17, 520-552
        .. [5] T.A. Halgren, "Merck molecular force field. III. Molecular
               geometries and vibrational frequencies for MMFF94." (1996)
               J. Comput. Chem. 17, 553-586
        .. [6] T.A. Halgren and R.B. Nachbar, "Merck molecular force field.
               IV. Conformational energies and geometries for MMFF94." (1996)
               J. Comput. Chem. 17, 587-615
        .. [7] T.A. Halgren, "Merck molecular force field. V. Extension of
               MMFF94 using experimental data, additional computational data,
               and empirical rules." (1996) J. Comput. Chem. 17, 616-641

        If using the 'MMFF94s' force field please cite the above along with:
        .. [8] T.A. Halgren, "MMFF VI. MMFF94s option for energy minimization
               studies." (1999) J. Comput. Chem. 20, 720-729

        If using the 'UFF' force field please cite the following:
        .. [3] Rappe, A.K., Casewit, C.J., Colwell, K.S., Goddard, W.A. III,
               Skiff, W.M. "UFF, a full periodic table force field for
               molecular mechanics and molecular dynamics simulations." (1992)
               J. Am. Chem. Soc. 114, 10024-10039

        If using the 'GAFF' force field please cite the following:
        .. [3] Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A., Case, D.A.
               "Development and testing of a general AMBER force field" (2004)
               J. Comput. Chem. 25, 1157-1174

        If using the 'Ghemical' force field please cite the following:
        .. [3] T. Hassinen and M. Perakyla, "New energy terms for reduced
               protein models implemented in an off-lattice force field" (2001)
               J. Comput. Chem. 22, 1229-1242
        """

        openbabel = import_('openbabel')
        md = import_('mdtraj')
        from mdtraj.core.element import get_by_symbol
        for particle in self.particles():
            try:
                get_by_symbol(particle.name)
            except KeyError:
                raise MBuildError("Element name {} not recognized. Cannot "
                                  "perform minimization."
                                  "".format(particle.name))

        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol2", "pdb")
        mol = openbabel.OBMol()

        obConversion.ReadFile(mol, os.path.join(tmp_dir, "un-minimized.mol2"))

        ff = openbabel.OBForceField.FindForceField(forcefield)
        if ff is None:
            raise MBuildError("Force field '{}' not supported for energy "
                              "minimization. Valid force fields are 'MMFF94', "
                              "'MMFF94s', 'UFF', 'GAFF', and 'Ghemical'."
                              "".format(forcefield))
        warn(
            "Performing energy minimization using the Open Babel package. Please "
            "refer to the documentation to find the appropriate citations for "
            "Open Babel and the {} force field".format(forcefield))
        ff.Setup(mol)
        if algorithm == 'steep':
            ff.SteepestDescent(steps)
        elif algorithm == 'md':
            ff.MolecularDynamicsTakeNSteps(steps, 300)
        elif algorithm == 'cg':
            ff.ConjugateGradients(steps)
        else:
            raise MBuildError("Invalid minimization algorithm. Valid options "
                              "are 'steep', 'cg', and 'md'.")
        ff.UpdateCoordinates(mol)

        obConversion.WriteFile(mol, os.path.join(tmp_dir, 'minimized.pdb'))

    def save(self, filename, show_ports=False, forcefield_name=None,
             forcefield_files=None, forcefield_debug=False, box=None,
             overwrite=False, residues=None, combining_rule='lorentz',
             foyer_kwargs=None, **kwargs):
        """Save the Compound to a file.

        Parameters
        ----------
        filename : str
            Filesystem path in which to save the trajectory. The extension or
            prefix will be parsed and control the format. Supported
            extensions are: 'hoomdxml', 'gsd', 'gro', 'top',
            'lammps', 'lmp', 'mcf'
        show_ports : bool, optional, default=False
            Save ports contained within the compound.
        forcefield_files : str, optional, default=None
            Apply a forcefield to the output file using a forcefield provided
            by the `foyer` package.
        forcefield_name : str, optional, default=None
            Apply a named forcefield to the output file using the `foyer`
            package, e.g. 'oplsaa'. Forcefields listed here:
            https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields
        forcefield_debug : bool, optional, default=False
            Choose level of verbosity when applying a forcefield through `foyer`.
            Specifically, when missing atom types in the forcefield xml file,
            determine if the warning is condensed or verbose.
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be written to the output file. If 'None', a
            bounding box is used with 0.25nm buffers at each face to avoid
            overlapping atoms.
        overwrite : bool, optional, default=False
            Overwrite if the filename already exists
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        combining_rule : str, optional, default='lorentz'
            Specify the combining rule for nonbonded interactions. Only relevant
            when the `foyer` package is used to apply a forcefield. Valid
            options are 'lorentz' and 'geometric', specifying Lorentz-Berthelot
            and geometric combining rules respectively.
        foyer_kwargs : dict, optional, default=None
            Keyword arguments to provide to `foyer.Forcefield.apply`.
        **kwargs
            Depending on the file extension these will be passed to either
            `write_gsd`, `write_hoomdxml`, `write_lammpsdata`,
            `write_mcf`, or `parmed.Structure.save`.
            See https://parmed.github.io/ParmEd/html/structobj/parmed.structure.Structure.html#parmed.structure.Structure.save

        Other Parameters
        ----------------
        ref_distance : float, optional, default=1.0
            Normalization factor used when saving to .gsd and .hoomdxml formats
            for converting distance values to reduced units.
        ref_energy : float, optional, default=1.0
            Normalization factor used when saving to .gsd and .hoomdxml formats
            for converting energy values to reduced units.
        ref_mass : float, optional, default=1.0
            Normalization factor used when saving to .gsd and .hoomdxml formats
            for converting mass values to reduced units.
        atom_style: str, default='full'
            Defines the style of atoms to be saved in a LAMMPS data file. The following atom
            styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
            see http://lammps.sandia.gov/doc/atom_style.html for more
            information on atom styles.
        unit_style: str, default='real'
            Defines to unit style to be save in a LAMMPS data file.  Defaults to 'real' units.
            Current styles are supported: 'real', 'lj'
            see https://lammps.sandia.gov/doc/99/units.html for more information
            on unit styles

        Notes
        ------
        When saving the compound as a json, only the following arguments are used:
            - filename
            - show_ports

        See Also
        --------
        conversion.save : Main saver logic
        formats.gsdwrite.write_gsd : Write to GSD format
        formats.hoomdxml.write_hoomdxml : Write to Hoomd XML format
        formats.xyzwriter.write_xyz : Write to XYZ format
        formats.lammpsdata.write_lammpsdata : Write to LAMMPS data format
        formats.cassandramcf.write_mcf : Write to Cassandra MCF format
        formats.json_formats.compound_to_json : Write to a json file

        """
        conversion.save(self, filename, show_ports, forcefield_name,
             forcefield_files, forcefield_debug, box,
             overwrite, residues, combining_rule, foyer_kwargs, **kwargs)


    def translate(self, by):
        """Translate the Compound by a vector

        Parameters
        ----------
        by : np.ndarray, shape=(3,), dtype=float

        """
        new_positions = _translate(self.xyz_with_ports, by)
        self.xyz_with_ports = new_positions

    def translate_to(self, pos):
        """Translate the Compound to a specific position

        Parameters
        ----------
        pos : np.ndarray, shape=3(,), dtype=float

        """
        self.translate(pos - self.center)

    def rotate(self, theta, around):
        """Rotate Compound around an arbitrary vector.

        Parameters
        ----------
        theta : float
            The angle by which to rotate the Compound, in radians.
        around : np.ndarray, shape=(3,), dtype=float
            The vector about which to rotate the Compound.

        """
        new_positions = _rotate(self.xyz_with_ports, theta, around)
        self.xyz_with_ports = new_positions

    def spin(self, theta, around):
        """Rotate Compound in place around an arbitrary vector.

        Parameters
        ----------
        theta : float
            The angle by which to rotate the Compound, in radians.
        around : np.ndarray, shape=(3,), dtype=float
            The axis about which to spin the Compound.

        """
        around = np.asarray(around).reshape(3)
        center_pos = self.center
        self.translate(-center_pos)
        self.rotate(theta, around)
        self.translate(center_pos)

    # Interface to Trajectory for reading/writing .pdb and .mol2 files.
    # -----------------------------------------------------------------
    def from_trajectory(self, traj, frame=-1, coords_only=False,
            infer_hierarchy=True):
        """Extract atoms and bonds from a md.Trajectory.

        Will create sub-compounds for every chain if there is more than one
        and sub-sub-compounds for every residue.

        Parameters
        ----------
        traj : mdtraj.Trajectory
            The trajectory to load.
        frame : int, optional, default=-1 (last)
            The frame to take coordinates from.
        coords_only : bool, optional, default=False
            Only read coordinate information
        infer_hierarchy : bool, optional, default=True
            If True, infer compound hierarchy from chains and residues

        See also
        --------
        mbuild.conversion.from_trajectory
        """
        conversion.from_trajectory(traj=traj, compound=self, frame=frame,
            coords_only=coords_only, infer_hierarchy=True)

    def to_trajectory(self, show_ports=False, chains=None,
                      residues=None, box=None):
        """Convert to an md.Trajectory and flatten the compound.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Include all port atoms when converting to trajectory.
        chains : mb.Compound or list of mb.Compound
            Chain types to add to the topology
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be used when converting to a `Trajectory`.
            If 'None', a bounding box is used with a 0.5nm buffer in each
            dimension. to avoid overlapping atoms, unless `self.periodicity`
            is not None, in which case those values are used for the
            box lengths.

        Returns
        -------
        trajectory : md.Trajectory

        See also
        --------
        mbuild.conversion.to_trajectory

        """
        return conversion.to_trajectory(compound=self, show_ports=show_ports,
            chains=chains, residues=residues, box=box)

    def from_parmed(self, structure, coords_only=False,
            infer_hierarchy=True):
        """Extract atoms and bonds from a pmd.Structure.

        Will create sub-compounds for every chain if there is more than one
        and sub-sub-compounds for every residue.

        Parameters
        ----------
        structure : pmd.Structure
            The structure to load.
        coords_only : bool
            Set preexisting atoms in compound to coordinates given by structure.
        infer_hierarchy : bool, optional, default=True
            If true, infer compound hierarchy from chains and residues
        """
        conversion.from_parmed(structure=structure,compound=self,
            coords_only=coords_only, infer_hierarchy=infer_hierarchy)

    def to_parmed(self, box=None, title='', residues=None, show_ports=False,
            infer_residues=False):
        """Create a ParmEd Structure from a Compound.

        Parameters
        ----------
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be used when converting to a `Structure`.
            If 'None', a bounding box is used with 0.25nm buffers at
            each face to avoid overlapping atoms, unless `self.periodicity`
            is not None, in which case those values are used for the
            box lengths.
        title : str, optional, default=self.name
            Title/name of the ParmEd Structure
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        show_ports : boolean, optional, default=False
            Include all port atoms when converting to a `Structure`.
        infer_residues : bool, optional, default=False
            Attempt to assign residues based on names of children.

        Returns
        -------
        parmed.structure.Structure
            ParmEd Structure object converted from self

        See Also
        --------
        mbuild.conversion.to_parmed
        parmed.structure.Structure : Details on the ParmEd Structure object

        """
        return conversion.to_parmed(compound=self, box=box, title=title,
            residues=residues, show_ports=show_ports, infer_residues=infer_residues)

    def to_networkx(self, names_only=False):
        """Create a NetworkX graph representing the hierarchy of a Compound.

        Parameters
        ----------
        names_only : bool, optional, default=False
        Store only the names of the
            compounds in the graph, appended with their IDs, for distinction even
            if they have the same name. When set to False, the default behavior,
            the nodes are the compounds themselves.

        Returns
        -------
        G : networkx.DiGraph

        Notes
        -----
        This digraph is not the bondgraph of the compound.

        See Also
        --------
        mbuild.conversion.to_networkx
        mbuild.bond_graph
        """
        return conversion.to_networkx(compound=self, names_only=names_only)

    def to_pybel(self, box=None, title='', residues=None, show_ports=False,
            infer_residues=False):
        """ Create a pybel.Molecule from a Compound

        Parameters
        ---------
        box : mb.Box, def None
        title : str, optional, default=self.name
            Title/name of the ParmEd Structure
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        show_ports : boolean, optional, default=False
            Include all port atoms when converting to a `Structure`.
        infer_residues : bool, optional, default=False
            Attempt to assign residues based on names of children

        Returns
        ------
        pybel.Molecule

        See also
        --------
        mbuild.conversion.to_pybel

        Notes
        -----
        Most of the mb.Compound is first converted to openbabel.OBMol
        And then pybel creates a pybel.Molecule from the OBMol
        Bond orders are assumed to be 1
        OBMol atom indexing starts at 1, with spatial dimension Angstrom

        """
        return conversion.to_pybel(compound=self, box=box, title=title,
            residues=residues, show_ports=show_ports)

    def from_pybel(self, pybel_mol, use_element=True, coords_only=False,
            infer_hierarchy=True, ignore_box_warn=False):
        """Create a Compound from a Pybel.Molecule

        Parameters
        ---------
        pybel_mol: pybel.Molecule
        use_element : bool, default True
            If True, construct mb Particles based on the pybel Atom's element.
            If False, construcs mb Particles based on the pybel Atom's type
        coords_only : bool, default False
            Set preexisting atoms in compound to coordinates given by
            structure.  Note: Not yet implemented, included only for parity
            with other conversion functions
        infer_hierarchy : bool, optional, default=True
            If True, infer hierarchy from residues
        ignore_box_warn : bool, optional, default=False
            If True, ignore warning if no box is present.

        See also
        --------
        mbuild.conversion.from_pybel

        """
        conversion.from_pybel(pybel_mol=pybel_mol, compound=self,
            use_element=use_element, coords_only=coords_only,
            ignore_box_warn=ignore_box_warn)

    def to_intermol(self, molecule_types=None): # pragma: no cover
        """Create an InterMol system from a Compound.

        Parameters
        ----------
        molecule_types : list or tuple of subclasses of Compound

        Returns
        -------
        intermol_system : intermol.system.System

        See also
        --------
        mbuild.conversion.to_intermol

        """
        return conversion.to_intermol(compound=self, molecule_types=None)

    def get_smiles(self):
        """Get SMILES string for compound

        Bond order is guessed with pybel and may lead to incorrect SMILES
        strings.

        Returns
        -------
        smiles_string: str
        """

        pybel_cmp = self.to_pybel()
        pybel_cmp.OBMol.PerceiveBondOrders()
        # we only need the smiles string
        smiles = pybel_cmp.write().split()[0]
        return smiles

    def __getitem__(self, selection):
        if isinstance(selection, int):
            return list(self.particles())[selection]
        if isinstance(selection, str):
            if selection not in self.labels:
                raise MBuildError('{}[\'{}\'] does not exist.'.format(self.name,selection))
            return self.labels.get(selection)

    def __repr__(self):
        descr = list('<')
        descr.append(self.name + ' ')

        if self.children:
            descr.append('{:d} particles, '.format(self.n_particles))
            if any(self.periodicity):
                descr.append('periodicity: {}, '.format(self.periodicity))
            else:
                descr.append('non-periodic, ')
        else:
            descr.append('pos=({: .4f},{: .4f},{: .4f}), '.format(*self.pos))

        descr.append('{:d} bonds, '.format(self.n_bonds))

        descr.append('id: {}>'.format(id(self)))
        return ''.join(descr)

    def _clone(self, clone_of=None, root_container=None):
        """A faster alternative to deepcopying.

        Does not resolve circular dependencies. This should be safe provided
        you never try to add the top of a Compound hierarchy to a
        sub-Compound. Clones compound hierarchy only, not the bonds.
        """
        if root_container is None:
            root_container = self
        if clone_of is None:
            clone_of = dict()

        # If this compound has already been cloned, return that.
        if self in clone_of:
            return clone_of[self]

        # Otherwise we make a new clone.
        cls = self.__class__
        newone = cls.__new__(cls)

        # Remember that we're cloning the new one of self.
        clone_of[self] = newone

        newone.name = deepcopy(self.name)
        newone._element = deepcopy(self.element)
        newone.periodicity = deepcopy(self.periodicity)
        newone._pos = deepcopy(self._pos)
        newone.port_particle = deepcopy(self.port_particle)
        newone._check_if_contains_rigid_bodies = deepcopy(
            self._check_if_contains_rigid_bodies)
        newone._contains_rigid = deepcopy(self._contains_rigid)
        newone._rigid_id = deepcopy(self._rigid_id)
        newone._charge = deepcopy(self._charge)
        if hasattr(self, 'index'):
            newone.index = deepcopy(self.index)

        if self.children is None:
            newone.children = None
        else:
            newone.children = OrderedSet()
        # Parent should be None initially.
        newone.parent = None
        newone.labels = OrderedDict()
        newone.referrers = set()
        newone.bond_graph = None

        # Add children to clone.
        if self.children:
            for child in self.children:
                newchild = child._clone(clone_of, root_container)
                newone.children.add(newchild)
                newchild.parent = newone

        # Copy labels, except bonds with atoms outside the hierarchy.
        if self.labels:
            for label, compound in self.labels.items():
                if not isinstance(compound, list):
                    newone.labels[label] = compound._clone(
                        clone_of, root_container)
                    compound.referrers.add(clone_of[compound])
                else:
                    # compound is a list of compounds, so we create an empty
                    # list, and add the clones of the original list elements.
                    newone.labels[label] = []
                    for subpart in compound:
                        newone.labels[label].append(
                            subpart._clone(clone_of, root_container))
                        # Referrers must have been handled already, or the will
                        # be handled

        newone.box = deepcopy(self.box)
        return newone

    def _clone_bonds(self, clone_of=None):
        """While cloning, clone the bond of the source compound to clone compound"""
        newone = clone_of[self]
        for c1, c2 in self.bonds():
            try:
                newone.add_bond((clone_of[c1], clone_of[c2]))
            except KeyError:
                raise MBuildError(
                    "Cloning failed. Compound contains bonds to "
                    "Particles outside of its containment hierarchy.")


Particle = Compound
