"""Module for working with mBuild Compounds."""

import itertools
import logging
import os
import tempfile
from collections import OrderedDict
from collections.abc import Iterable
from copy import deepcopy
from typing import Sequence

import ele
import networkx as nx
import numpy as np
from boltons.setutils import IndexedSet
from ele.element import Element, element_from_name, element_from_symbol
from ele.exceptions import ElementError
from treelib import Tree

from mbuild import conversion
from mbuild.bond_graph import BondGraph
from mbuild.box import Box
from mbuild.coordinate_transform import _rotate, _translate
from mbuild.exceptions import MBuildError
from mbuild.periodic_kdtree import PeriodicKDTree
from mbuild.utils.io import import_, run_from_ipython
from mbuild.utils.jsutils import overwrite_nglview_default

__all__ = ["clone", "Compound", "Particle"]

logger = logging.getLogger(__name__)


def clone(existing_compound, clone_of=None, root_container=None):
    """Clone Compound.

    A faster alternative to deepcopying. Does not resolve circular dependencies.
    This should be safe provided you never try to add the top of a Compound
    hierarchy to a sub-Compound.

    Parameters
    ----------
    existing_compound : mbuild.Compound
        Existing Compound that will be copied
    clone_of : dict, optional, default None,
    root_container : mb.Compound, optional, default None,
    """
    if clone_of is None:
        clone_of = dict()

    newone = existing_compound._clone(clone_of=clone_of, root_container=root_container)
    existing_compound._clone_bonds(clone_of=clone_of)
    return newone


class Compound(object):
    """A building block in the mBuild hierarchy.

    Compound is the superclass of all composite building blocks in the mBuild
    hierarchy. That is, all composite building blocks must inherit from
    compound, either directly or indirectly. The design of Compound follows the
    Composite design pattern::

        @book{DesignPatterns,
            author = "Gamma, Erich and Helm, Richard and Johnson, Ralph and
            Vlissides, John M.",
            title = "Design Patterns",
            subtitle = "Elements of Reusable Object-Oriented Software",
            year = "1995",
            publisher = "Addison-Wesley",
            note = "p. 395",
            ISBN = "0-201-63361-2",
        }

    with Compound being the composite, and Particle playing the role of the
    primitive (leaf) part, where Particle is in fact simply an alias to the
    Compound class.

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
    mass : float, optional, default=None
        The mass of the compound. If none is set, then will try to
        infer the mass from a compound's element attribute.
        If neither `mass` or `element` are specified, then the
        mass will be None.
    charge : float, optional, default=0.0
        Currently not used. Likely removed in next release.
    periodicity : tuple of bools, length=3, optional, default=None
        Whether the Compound is periodic in the x, y, and z directions.
        If None is provided, the periodicity is set to (False, False, False)
        which is non-periodic in all directions.
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
    children : list
        Contains all children (other Compounds).
    labels : OrderedDict
        Labels to Compound/Atom mappings. These do not necessarily need not be
        in self.children.
    parent : mb.Compound
        The parent Compound that contains this part. Can be None if this
        compound is the root of the containment hierarchy.
    referrers : set
        Other compounds that reference this part with labels.
    boundingbox : mb.Box
        The bounds (xmin, xmax, ymin, ymax, zmin, zmax) of particles in Compound
    center
    mass
    n_particles
    n_bonds
    root
    xyz
    xyz_with_ports
    """

    def __init__(
        self,
        subcompounds=None,
        name=None,
        pos=None,
        mass=None,
        charge=None,
        periodicity=None,
        box=None,
        element=None,
        port_particle=False,
    ):
        super(Compound, self).__init__()

        if name:
            if not isinstance(name, str):
                raise ValueError(
                    f"Compound.name should be a string. You passed {name}."
                )
            self.name = name
        else:
            self.name = self.__class__.__name__

        if pos is not None:
            self._pos = np.asarray(pos, dtype=float)
        else:
            self._pos = np.zeros(3)

        self.parent = None
        self.children = list()
        self.labels = OrderedDict()
        self.referrers = set()

        self.bond_graph = BondGraph()
        self.bond_graph.add_node(self)

        self.port_particle = port_particle

        self.element = element
        if mass and float(mass) < 0.0:
            raise ValueError("Cannot set a Compound mass value less than zero")
        self._box = box
        if periodicity is not None:
            self.periodicity = periodicity
        else:
            self.periodicity = (False, False, False)
        # self.add() must be called after labels and children are initialized.
        if subcompounds:
            if charge:
                raise MBuildError(
                    "Can't set the charge of a Compound with subcompounds."
                )
            if mass:
                raise MBuildError(
                    "Can't set the mass of a Compound with subcompounds. "
                )
            self._charge = None
            self._mass = mass
            self.add(subcompounds)
        else:
            self._charge = charge
            self._mass = mass

    def particles(self, include_ports=False):
        """Return all Particles of the Compound.

        Parameters
        ----------
        include_ports : bool, optional, default=False
            Include port particles

        Yields
        ------
        mb.Compound
            The next Particle in the Compound
        """
        if not self.children:
            yield self
        else:
            for particle in self._particles(include_ports):
                yield particle

    def _particles(self, include_ports=False):
        """Return all Particles of the Compound."""
        for child in self.successors():
            if not child.children:
                if include_ports or not child.port_particle:
                    yield child

    def successors(self):
        """Yield Compounds below self in the hierarchy.

        Yields
        ------
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
        int,
            The number of Particles in the Compound

        """
        if not self.children:
            return 1
        else:
            return self._n_particles(include_ports=False)

    def _n_particles(self, include_ports=False):
        """Return the number of Particles in the Compound."""
        return sum(1 for _ in self._particles(include_ports))

    def _contains_only_ports(self):
        if self.children:
            for part in self.children:
                if not part.port_particle:
                    return False
        return True

    def print_hierarchy(self, print_full=False, index=None, show_tree=True):
        """Print the hierarchy of the Compound.

        Parameters
        ----------
        print_full: bool, optional, default=False
            The full hierarchy will be printed, rather than condensing
            compounds with identical topologies.
            Topologies are considered identical if they have the same name,
            contain the number and names of children,
            contain the same number and names of particles,
            and the same number of bonds.
        index: int, optional, default=None
            Print the branch of the first level of the hiearchy
            corresponding to the value specified by index.
            This only applies when print_full is True.
        show_tree: bool, optional, default=True
            If False, do not print the tree to the screen.

        Returns
        -------
        tree, treelib.tree.Tree, hierarchy of the compound as a tree

        """
        tree = Tree()

        # loop through the hierarchy saving the data to an array hh
        if print_full:
            hh = [h for h in self._get_hierarchy()]
        else:
            hh = [h for h in self._get_hierarchy_nodup()]

        # if our compound does not have any children we need to call n_direct_bonds instead of n_bonds
        if len(self.children) == 0:
            n_bonds = self.n_direct_bonds
        else:
            n_bonds = self.n_bonds

        # add the top level compound to create the top level of the tree
        # note that node identifiers passed as the second argument
        # correspond to the compound id
        tree.create_node(
            f"{self.name}, {self.n_particles} particles, {n_bonds} bonds, {len(self.children)} children",
            f"{id(self)}",
        )

        # if index is specified, ensure we are not selecting an index out of range
        if index is not None:
            if index >= len(self.children):
                raise MBuildError(
                    f"Index {index} out of range. The number of first level nodes in the tree is {len(self.children)}."
                )

        count = -1

        for h in hh:
            if len(h["comp"].children) == 0:
                n_bonds = h["comp"].n_direct_bonds
            else:
                n_bonds = h["comp"].n_bonds
            if h["level"] == 0:
                count = count + 1
            if print_full:
                if index is None:
                    tree.create_node(
                        f"[{h['comp'].name}]: {h['comp'].n_particles} particles, {n_bonds} bonds, {len(h['comp'].children)} children",
                        f"{h['comp_id']}",
                        f"{h['parent_id']}",
                    )
                elif count == index:
                    tree.create_node(
                        f"[{h['comp'].name}]: {h['comp'].n_particles} particles, {n_bonds} bonds, {len(h['comp'].children)} children",
                        f"{h['comp_id']}",
                        f"{h['parent_id']}",
                    )
            else:
                tree.create_node(
                    f"[{h['comp'].name} x {h['n_dup']}], {h['comp'].n_particles} particles, {n_bonds} bonds, {len(h['comp'].children)} children",
                    f"{h['comp_id']}",
                    f"{h['parent_id']}",
                )
        if show_tree:
            print(tree)
        return tree

    def _get_hierarchy(self, level=0):
        """Return an array of dictionaries corresponding to hierarchy of the compound, recursively."""
        if not self.children:
            return
        for child in self.children:
            yield {
                "level": level,
                "parent_id": id(self),
                "comp_id": id(child),
                "comp": child,
            }
            for subchild in child._get_hierarchy(level + 1):
                yield subchild

    def _get_hierarchy_nodup(self, level=0):
        """Return an array of dictionaries corresponding to hierarchy of the compound, recursively.

        This routine will identify any duplicate compounds at a given level, including the number of
        duplicates for each compound. Compounds are considered to be identical if the name,
        number of children, and number of particles are the same at the same level.
        """
        if not self.children:
            return

        duplicates = {}
        for child in self.children:
            part_string = "".join([part.name for part in child.particles()])
            child_string = "".join([child.name for child in child.children])

            if len(child.children) == 0:
                n_bonds = child.n_direct_bonds
            else:
                n_bonds = child.n_bonds

            identifier = f"{child.name}_{len(child.children)}_{child_string}_{child.n_particles}_{part_string}_{n_bonds}"

            if identifier not in duplicates:
                duplicates[identifier] = [1, True]
            else:
                duplicates[identifier][0] += 1

        for child in self.children:
            part_string = "".join([part.name for part in child.particles()])
            child_string = "".join([child.name for child in child.children])

            if len(child.children) == 0:
                n_bonds = child.n_direct_bonds
            else:
                n_bonds = child.n_bonds

            identifier = f"{child.name}_{len(child.children)}_{child_string}_{child.n_particles}_{part_string}_{n_bonds}"

            if duplicates[identifier][1]:
                yield {
                    "level": level,
                    "parent_id": id(self),
                    "comp_id": id(child),
                    "comp": child,
                    "n_dup": duplicates[identifier][0],
                }

                for subchild in child._get_hierarchy_nodup(level + 1):
                    yield subchild
                duplicates[identifier][1] = False

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
        """Get the Compound at the top of self's hierarchy.

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
        """Return all Particles of the Compound with a specific name.

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
        """Return all Particles of the Compound with a specific element.

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
    def mass(self):
        """Return the total mass of a compound.

        If the compound contains children compouds, the total mass of all
        children compounds is returned.
        If the compound contains element information (Compound.element) then
        the mass is inferred from the elemental mass.
        If Compound.mass has been set explicitly, then it will override the
        mass inferred from Compound.element.
        If neither of a Compound's element or mass attributes have been set,
        then a mass of zero is returned.
        """
        if self._contains_only_ports():
            return self._particle_mass(self)
        else:
            particle_masses = [self._particle_mass(p) for p in self.particles()]
            if None in particle_masses:
                logger.info(
                    f"Some particle of {self} does not have mass."
                    "They will not be accounted for during this calculation."
                )
            filtered_masses = [mass for mass in particle_masses if mass is not None]
            return sum(filtered_masses) if filtered_masses else None

    @staticmethod
    def _particle_mass(particle):
        if particle._mass is not None:
            return particle._mass
        else:
            if particle.element:
                return particle.element.mass
            else:
                return None

    @mass.setter
    def mass(self, value):
        if self._contains_only_ports() is False:
            raise MBuildError(
                "Cannot set the mass of a Compound containing children compounds"
            )

        value = float(value)
        if value < 0.0:
            raise ValueError("Cannot set a mass value less than zero")
        self._mass = value

    @property
    def charge(self):
        """Return the total charge of a compound.

        If the compound contains children compouds, the total charge of all
        children compounds is returned.

        If the charge of a particle has not been explicitly set
        then the particle's charge is None, and are not used when
        calculating the total charge.
        """
        if self._contains_only_ports():
            return self._particle_charge(self)
        charges = [p._charge for p in self.particles()]
        if None in charges:
            logger.info(
                f"Some particle of {self} does not have a charge. "
                "They will not be accounted for during this calculation."
            )
        filtered_charges = [charge for charge in charges if charge is not None]
        return sum(filtered_charges) if filtered_charges else None

    @staticmethod
    def _particle_charge(particle):
        """Return charge of a Compound with no children."""
        return particle._charge

    @charge.setter
    def charge(self, value):
        if self._contains_only_ports():
            self._charge = value
        else:
            raise AttributeError(
                "charge is immutable for Compounds that are "
                "not at the bottom of the containment hierarchy."
            )

    def add(
        self,
        new_child,
        label=None,
        containment=True,
        replace=False,
        inherit_periodicity=None,
        inherit_box=False,
        check_box_size=True,
    ):
        """Add a part to the Compound.

        Note:
            This does not necessarily add the part to self.children but may
            instead be used to add a reference to the part to self.labels. See
            'containment' argument.

        Parameters
        ----------
        new_child : mb.Compound or list-like of mb.Compound
            The object(s) to be added to this Compound.
        label : str, or list-like of str, optional, default None
            A descriptive string for the part; if a list, must be the same length/shape as new_child.
        containment : bool, optional, default=True
            Add the part to self.children.
        replace : bool, optional, default=True
            Replace the label if it already exists.
        inherit_periodicity : bool, optional, default=True
            Replace the periodicity of self with the periodicity of the
            Compound being added
        inherit_box: bool, optional, default=False
            Replace the box of self with the box of the Compound being added
        check_box_size : bool, optional, default=True
            Checks and warns if compound box is smaller than its bounding box after adding new_child.
        """
        # Support batch add via lists, tuples and sets.
        # If iterable, we will first compose all the bondgraphs of individual
        # Compounds in the list for efficiency
        from mbuild.port import Port

        if isinstance(new_child, Iterable) and not isinstance(new_child, str):
            compound_list = [c for c in _flatten_list(new_child)]
            if label is not None and isinstance(label, (list, tuple)):
                label_list = [c for c in _flatten_list(label)]
                if len(label_list) != len(compound_list):
                    raise ValueError(
                        "The list-like object for label must be the same length as"
                        "the list-like object of child Compounds. "
                        f"total length of labels: {len(label_list)}, new_child: {len(new_child)}."
                    )
            temp_bond_graphs = []
            for child in compound_list:
                # create a list of bond graphs of the children to add
                if containment:
                    if child.bond_graph and not isinstance(self, Port):
                        temp_bond_graphs.append(child.bond_graph)

            # compose children bond_graphs; make sure we actually have graphs to compose
            children_bond_graph = None
            if len(temp_bond_graphs) != 0:
                children_bond_graph = nx.compose_all(temp_bond_graphs)

            if (
                temp_bond_graphs
                and not isinstance(self, Port)
                and children_bond_graph is not None
            ):
                # If anything is added at self level, it is no longer a particle
                # search for self in self.root.bond_graph and remove self
                if self.root.bond_graph.has_node(self):
                    self.root.bond_graph.remove_node(self)
                # compose the bond graph of all the children with the root
                self.root.bond_graph = nx.compose(
                    self.root.bond_graph, children_bond_graph
                )
            for i, child in enumerate(compound_list):
                child.bond_graph = None
                if label is not None:
                    self.add(
                        child,
                        label=label_list[i],
                        check_box_size=False,
                    )
                else:
                    self.add(child, check_box_size=False)

            return

        if not isinstance(new_child, Compound):
            raise ValueError(
                "Only objects that inherit from mbuild.Compound can be added "
                f"to Compounds. You tried to add '{new_child}'."
            )
        if self._mass is not None and not isinstance(new_child, Port):
            logger.info(
                f"{self} has a pre-defined mass of {self._mass}, "
                "which will be reset to zero now that it contains children "
                "compounds."
            )
            self._mass = 0

        # Create children and labels on the first add operation
        if self.children is None:
            self.children = list()
        if self.labels is None:
            self.labels = OrderedDict()

        if containment:
            if new_child.parent is not None:
                raise MBuildError(
                    f"Part {new_child} already has a parent: {new_child.parent}"
                )
            self.children.append(new_child)
            new_child.parent = self

            if new_child.bond_graph is not None and not isinstance(self, Port):
                # If anything is added at self level, it is no longer a particle
                # search for self in self.root.bond_graph and remove self
                if self.root.bond_graph.has_node(self):
                    self.root.bond_graph.remove_node(self)
                # Compose bond_graph of new child
                self.root.bond_graph = nx.compose(
                    self.root.bond_graph, new_child.bond_graph
                )

                new_child.bond_graph = None

        # Add new_part to labels. Does not currently support batch add.
        if label is None:
            label = f"{new_child.name}[$]"

        if label.endswith("[$]"):
            label = label[:-3]
            all_label = "all-" + label + "s"
            if all_label not in self.labels:
                self.labels[all_label] = []
            label_pattern = label + "[{}]"

            count = len(self.labels[all_label])
            self.labels[all_label].append(new_child)
            label = label_pattern.format(count)

        if not replace and label in self.labels:
            raise MBuildError(f'Label "{label}" already exists in {self}.')
        else:
            self.labels[label] = new_child
        new_child.referrers.add(self)

        if inherit_periodicity and isinstance(new_child, Compound):
            self.periodicity = new_child.periodicity

        # If parent has no box --> inherit child box
        # If parent has box --> keep unless inherit_box == True
        # If inherit_box == True, parent box != None, child_box == None,
        # keep parent box anyway and warn
        if self.box is None:
            if new_child.box is not None:
                self.box = new_child.box
        else:
            if inherit_box:
                if new_child.box is None:
                    logger.info(
                        "The Compound you are adding has no box but "
                        "inherit_box=True. The box of the original "
                        "Compound will remain unchanged."
                    )
                else:
                    self.box = new_child.box
            else:
                if new_child.box is not None:
                    logger.info(
                        "The Compound you are adding has a box. "
                        "The box of the parent compound will be used. Use "
                        "inherit_box = True if you wish to replace the parent "
                        "compound box with that of Compound being added."
                    )

        # Check that bounding box is within box after adding compound
        if self.box and check_box_size:
            if (
                np.array(self.box.lengths) < np.array(self.get_boundingbox().lengths)
            ).any():
                logger.warning(
                    "After adding new Compound, Compound.box.lengths < "
                    "Compound.boundingbox.lengths. There may be particles "
                    "outside of the defined simulation box"
                )

    def remove(self, objs_to_remove, reset_labels=False):
        """Remove children from the Compound cleanly.

        Parameters
        ----------
        objs_to_remove : mb.Compound or list of mb.Compound
            The Compound(s) to be removed from self
        reset_labels : bool, optional, default=False
            If True, the Compound labels will be reset
        """
        # Preprocessing and validating input type
        from mbuild.port import Port

        if not hasattr(objs_to_remove, "__iter__"):
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
        particles_to_remove = set(
            [particle for obj in objs_to_remove for particle in obj.particles()]
        )

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
                    logger.warning(f"This will remove all particles in {self}")
            return

        for particle in particles_to_remove:
            _check_if_empty(particle)

        # Remove obj from bondgraph
        for removed_part in to_remove:
            self._remove(removed_part)

        # Remove references to object
        for removed_part in to_remove:
            if removed_part.parent is not None:
                removed_part.parent.children.remove(removed_part)
            self._remove_references(removed_part)

        # Remove ghost ports
        self._prune_ghost_ports()

        # Reorder labels
        if reset_labels:
            self.reset_labels()

    def reset_labels(self):
        """Reset Compound labels so that substituents and ports are renumbered, indexed from port[0] to port[N], where N-1 is the number of ports.

        Notes
        -----
        Will renumber the labels in a given Compound. Duplicated labels are named in the format "{name}[$]", where the $ stands in for the 0-indexed
        number in the Compound hierarchy with given "name".

        i.e. self.labels.keys() = ["CH2", "CH2", "CH2"] would transform into self.labels.keys() = ["CH2[0]", "CH2[1]", "CH2[2]"]
        and
        i.e. self.labels.keys() = ["CH2[1]", "CH2[3]", "CH2[5]"] would transform into self.labels.keys() = ["CH2[0]", "CH2[1]", "CH2[2]"]

        Additonally, if it doesn't exist, duplicated labels that are numbered as above with the "[$]" will also be put into a list index.
        self.labels.keys() = ["CH2", "CH2", "CH2"] would transform into self.labels.keys() = ["CH2[0]", "CH2[1]", "CH2[2]"] as shown above, but also
        have a label of self.labels["all-CH2s"], which is a list of all CH2 children in the Compound.
        """
        new_labels = OrderedDict()
        hoisted_children = {
            key: val
            for key, val in self.labels.items()
            if (
                not isinstance(val, list)
                and val.parent is not None
                and id(self) != id(val.parent)
            )
        }
        new_labels.update(hoisted_children)
        children_list = {
            id(val): [key, val]
            for key, val in self.labels.items()
            if (not isinstance(val, list))
        }
        for child in self.children:
            label = (
                children_list[id(child)][0]
                if "[" not in children_list[id(child)][0]
                else None
            )
            if label is None:
                if "Port" in child.name:
                    label = [
                        key for key, x in self.labels.items() if id(x) == id(child)
                    ][0]
                    if "port" in label:
                        label = "port[$]"
                else:
                    label = f"{child.name}[$]"
            if label.endswith("[$]"):
                label = label[:-3]
                all_label = "all-" + label + "s"
                if all_label not in new_labels:
                    new_labels[all_label] = []
                label_pattern = label + "[{}]"

                count = len(new_labels[all_label])
                new_labels[all_label].append(child)
                label = label_pattern.format(count)
            new_labels[label] = child
        self.labels = new_labels

    def _prune_ghost_ports(self):
        """Worker for remove(). Remove all ports whose anchor has been deleted."""
        all_ports_list = list(self.all_ports())
        particles = list(self.particles())
        for port in all_ports_list:
            if port.anchor not in particles:
                self._remove(port)
                port.parent.children.remove(port)
                self._remove_references(port)

    def _remove(self, removed_part):
        """Worker for remove(). Removes bonds."""
        if self.root.bond_graph.has_node(removed_part):
            for neighbor in nx.neighbors(self.root.bond_graph.copy(), removed_part):
                self.root.remove_bond((removed_part, neighbor))
            self.root.bond_graph.remove_node(removed_part)

    def _remove_references(self, removed_part):
        """Remove labels pointing to this part and vice versa."""
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

        return [port for port in self.labels.values() if isinstance(port, Port)]

    def all_ports(self):
        """Return all Ports referenced by this Compound and its successors.

        Returns
        -------
        list of mb.Compound
            A list of all Ports referenced by this Compound and its successors
        """
        from mbuild.port import Port

        return [s for s in self.successors() if isinstance(s, Port)]

    def available_ports(self):
        """Return all unoccupied Ports referenced by this Compound.

        Returns
        -------
        list of mb.Compound
            A list of all unoccupied ports referenced by the Compound
        """
        from mbuild.port import Port

        return [p for p in self.labels.values() if isinstance(p, Port) and not p.used]

    def direct_bonds(self, graph_depth=1):
        """Return a list of particles that this particle bonds to.

        Parameters
        ----------
        graph_depth : int, default=1
            Determines how many subsequent bonded neighbors to count.
            A value of 1 returns only paricles this particle is directly bonded to.
            A value of 2 returns direct bonded neighbors,
            plus their direct bonded neighbors.

        Returns
        -------
        List of mb.Compound

        See Also
        --------
        bond_graph.edges_iter : Iterations over all edges in a BondGraph
        Compound.n_direct_bonds : Returns the number of bonds a particle contains
        """
        if graph_depth <= 0 or not isinstance(graph_depth, int):
            raise ValueError("`graph_depth` must be an integer >= 1.")
        if list(self.particles()) != [self]:
            raise MBuildError(
                "The direct_bonds method can only "
                "be used on compounds at the bottom of their hierarchy."
            )
        if not self.parent:
            return None
        # Get all nodes within n edges (graph_depth=n) using BFS
        top_ancestor = [i for i in self.ancestors()][-1]
        neighbor_dict = nx.single_source_shortest_path_length(
            top_ancestor.bond_graph, self, cutoff=graph_depth
        )
        # Exclude the source node itself
        all_neighbors = {n for n, depth in neighbor_dict.items() if depth > 0}
        return all_neighbors

    def bonds(self, return_bond_order=False):
        """Return all bonds in the Compound and sub-Compounds.

        Parameters
        ----------
        return_bond_order : bool, optional, default=False
            Return the bond order of the bond as the 3rd argument in the tuple.
            bond order is returned as a dictionary with 'bo' as the key.
            If bond order is not set, the value will be set to 'default'.

        Yields
        ------
        tuple of mb.Compound
            The next bond in the Compound

        See Also
        --------
        bond_graph.edges_iter : Iterates over all edges in a BondGraph
        Compound.n_bonds : Returns the total number of bonds in the Compound and sub-Compounds
        """
        if self.root.bond_graph:
            if self.root == self:
                return self.root.bond_graph.edges(data=return_bond_order)
            else:
                return self.root.bond_graph.subgraph(self.particles()).edges(
                    data=return_bond_order
                )
        else:
            return iter(())

    @property
    def n_direct_bonds(self):
        """Return the number of bonds a particle is directly involved in.

        This method should only be used on on compounds at the bottom
        of their hierarchy (i.e. a particle).

        Returns
        -------
        int
            The number of compounds this compound is directly bonded to.
        """
        if list(self.particles()) != [self]:
            raise MBuildError(
                "The direct_bonds method can only "
                "be used on compounds at the bottom of their hierarchy."
            )
        if self.direct_bonds(graph_depth=1):
            return len(self.direct_bonds(graph_depth=1))
        else:
            return 0

    @property
    def n_bonds(self):
        """Return the total number of bonds in the Compound.

        Returns
        -------
        int
            The number of bonds in the Compound
        """
        if list(self.particles()) == [self]:
            raise MBuildError(
                "n_bonds cannot be used on Compounds "
                "at the bottom of their hierarchy (particles). "
                "Use n_direct_bonds instead."
            )
        return sum(1 for _ in self.bonds())

    def add_bond(self, particle_pair, bond_order=None):
        """Add a bond between two Particles.

        Parameters
        ----------
        particle_pair : indexable object, length=2, dtype=mb.Compound
            The pair of Particles to add a bond between
        bond_order : float, optional, default=None
            Bond order of the bond.
            Available options include "default", "single", "double",
            "triple", "aromatic" or "unspecified"
        """
        if self.root.bond_graph is None:
            self.root.bond_graph = BondGraph()
        if bond_order is None:
            bond_order = "default"
        else:
            if not isinstance(bond_order, str) or bond_order.lower() not in [
                "default",
                "single",
                "double",
                "triple",
                "aromatic",
                "unspecified",
            ]:
                raise ValueError(
                    "Invalid bond_order given. Available bond orders are: single",
                    "double",
                    "triple",
                    "aromatic",
                    "unspecified",
                )
        self.root.bond_graph.add_edge(
            particle_pair[0], particle_pair[1], bond_order=bond_order
        )

    def generate_bonds(self, name_a, name_b, dmin, dmax):
        """Add Bonds between all pairs of types a/b within [dmin, dmax].

        Parameters
        ----------
        name_a : str
            The name of one of the Particles to be in each bond
        name_b : str
            The name of the other Particle to be in each bond
        dmin : float
            The minimum distance (in nm) between Particles for considering a bond
        dmax : float
            The maximum distance (in nm) between Particles for considering a bond
        """
        if self.box is None:
            self.box = self.get_boundingbox()
        particle_kdtree = PeriodicKDTree.from_compound(compound=self, leafsize=10)
        particle_array = np.array(list(self.particles()))
        added_bonds = list()
        for p1 in self.particles_by_name(name_a):
            nearest = self.particles_in_range(
                p1,
                dmax,
                max_particles=20,
                particle_kdtree=particle_kdtree,
                particle_array=particle_array,
            )
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

    def freud_generate_bonds(self, name_a, name_b, dmin, dmax):
        """Add Bonds between all pairs of types a/b within [dmin, dmax].

        Parameters
        ----------
        name_a : str
            The name of one of the Particles to be in each bond
        name_b : str
            The name of the other Particle to be in each bond
        dmin : float
            The minimum distance (in nm) between Particles for considering a bond
        dmax : float
            The maximum distance (in nm) between Particles for considering a bond
        """
        freud = import_("freud")
        moved_positions, freud_box = self.to_freud()
        a_indices = []
        b_indices = []
        for i, part in enumerate(self.particles()):
            if part.name == name_a:
                a_indices.append(i)
            if part.name == name_b:
                b_indices.append(i)
        # If we are looking to create bonds between the same species
        # then the indices added to a_indices and b_indices will be identical.
        # In this case we need to make sure that we don't try to bond a particle
        # to itself  (i.e., excluded_ii = True).
        # If we are looking for bonds between two different species,
        # the indices we find for a and b will be distinct, with no overlap.
        # However, the way the code is structured, we don't actually pass
        # freud the indices, but rather a list of particle positions associated with each set of indices.
        # As such, the indices that freud sees will be (0, len(a_indices)) and (0, len(b_indices)), even though
        # they represent different actually particles. Thus, to get the right behavior we
        # must not exclude particles with the same index, and thus exclude_ii = False.
        if name_a == name_b:
            exclude_ii = True
        else:
            exclude_ii = False
        aq = freud.locality.AABBQuery(freud_box, moved_positions[b_indices])

        nlist = aq.query(
            moved_positions[a_indices],
            dict(r_min=dmin, r_max=dmax, exclude_ii=exclude_ii),
        ).toNeighborList()

        part_list = [part for part in self.particles(include_ports=False)]
        for i, j in nlist[:]:
            self.add_bond((part_list[a_indices[i]], part_list[b_indices[j]]))

    def remove_bond(self, particle_pair):
        """Delete a bond between a pair of Particles.

        Parameters
        ----------
        particle_pair : indexable object, length=2, dtype=mb.Compound
            The pair of Particles to remove the bond between
        """
        from mbuild.port import Port

        if self.root.bond_graph is None or not self.root.bond_graph.has_edge(
            *particle_pair
        ):
            raise MBuildError(
                "Bond between {} and {} doesn't exist!".format(*particle_pair)
            )
            return
        self.root.bond_graph.remove_edge(*particle_pair)
        bond_vector = particle_pair[0].pos - particle_pair[1].pos
        if np.allclose(bond_vector, np.zeros(3)):
            logger.warning(
                "Particles {} and {} overlap! Ports will not be added.".format(
                    *particle_pair
                )
            )
            return
        distance = np.linalg.norm(bond_vector)
        particle_pair[0].parent.add(
            Port(
                anchor=particle_pair[0],
                orientation=-bond_vector,
                separation=distance / 2,
            ),
            "port[$]",
        )
        particle_pair[1].parent.add(
            Port(
                anchor=particle_pair[1],
                orientation=bond_vector,
                separation=distance / 2,
            ),
            "port[$]",
        )

    @property
    def pos(self):
        """Get the position of the Compound.

        If the Compound contains children, returns the center.

        The position of a Compound containing children can't be set.
        """
        if not self.children:
            return self._pos
        else:
            return self.center

    @pos.setter
    def pos(self, value):
        if not self.children:
            self._pos = value
        else:
            raise MBuildError("Can't set position of Compound with children.")

    @property
    def periodicity(self):
        """Get the periodicity of the Compound."""
        return self._periodicity

    @periodicity.setter
    def periodicity(self, periods):
        if len(list(periods)) != 3:
            raise ValueError("Periodicity must be of length 3")
        if not all([isinstance(p, bool) for p in periods]):
            raise TypeError(
                "Periodicity values must be True/False; if you are trying to "
                "set the dimensions, use Compound.box."
            )
        self._periodicity = tuple(periods)

    @property
    def box(self):
        """Get the box of the Compound.

        Ports cannot have a box.
        """
        return self._box

    @box.setter
    def box(self, box):
        if box and not isinstance(box, Box):
            raise TypeError("box must be specified as an mbuild.Box")
        if self.port_particle and box is not None:
            raise ValueError("Ports cannot have a box")
        # Make sure the box is bigger than the bounding box
        if box is not None:
            if np.asarray((box.lengths < self.get_boundingbox().lengths)).any():
                logger.warning(
                    "Compound.box.lengths < Compound.boundingbox.lengths. "
                    "There may be particles outside of the defined "
                    "simulation box."
                )
        self._box = box

    @property
    def element(self):
        """Get the element of the Compound."""
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
            arr = np.fromiter(
                itertools.chain.from_iterable(p.pos for p in self.particles()),
                dtype=float,
            )
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
                    p.pos for p in self.particles(include_ports=True)
                ),
                dtype=float,
            )
            pos = arr.reshape((-1, 3))
        return pos

    @xyz.setter
    def xyz(self, arrnx3):
        """Set the positions of the particles in the Compound, excluding Ports.

        This function does not set the position of the ports.

        Parameters
        ----------
        arrnx3 : np.ndarray, shape=(n,3), dtype=float
            The new particle positions
        """
        arrnx3 = np.array(arrnx3)
        if not self.children:
            if not arrnx3.shape[0] == 1:
                raise ValueError(
                    f"Trying to set position of {self} with more than one"
                    f"coordinate: {arrnx3}"
                )
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(self._particles(include_ports=False), arrnx3):
                atom.pos = coords

    @xyz_with_ports.setter
    def xyz_with_ports(self, arrnx3):
        """Set the positions of the particles in the Compound, including Ports.

        Parameters
        ----------
        arrnx3 : np.ndarray, shape=(n,3), dtype=float
            The new particle positions
        """
        if not self.children:
            if not arrnx3.shape[0] == 1:
                raise ValueError(
                    f"Trying to set position of {self} with more than one"
                    f"coordinate: {arrnx3}"
                )
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(self._particles(include_ports=True), arrnx3):
                atom.pos = coords

    @property
    def center(self):
        """Get the cartesian center of the Compound based on its Particles.

        Returns
        -------
        np.ndarray, shape=(3,), dtype=float
            The cartesian center of the Compound based on its Particles
        """
        if np.all(np.isfinite(self.xyz)):
            return np.mean(self.xyz, axis=0)

    @property
    def mins(self):
        """Return the mimimum x, y, z coordinate of any particle in this compound."""
        return self.xyz.min(axis=0)

    @property
    def maxs(self):
        """Return the maximum x, y, z coordinate of any particle in this compound."""
        return self.xyz.max(axis=0)

    def is_independent(self):
        """Return True if there is no bond between particles of the Compound to an external Compound."""
        if not self.parent:
            # This is the very top level, and hence have to be independent
            return True
        elif not self.root.bond_graph.edges():
            # If there is no bond in the top level, then everything is independent
            return True
        else:
            # Cover the other cases
            for particle in self.particles():
                for neigh in nx.neighbors(self.root.bond_graph, particle):
                    if neigh not in self.particles():
                        return False
            return True

    def check_for_overlap(self, excluded_bond_depth, minimum_distance=0.10):
        """Check if a compound contains overlapping particles.

        Parameters:
        -----------
        excluded_bond_depth : int, required
            The depth of bonded neighbors to exclude from overlap check.
            see Compound.direct_bonds()
        minimum_distance : float, default=0.10
            Distance (in nanometers) used as the threshold in
            determining if a pair of particles overlap.

        Notes:
        ------
        If `minimum_distance` is set larger than existing bond lengths,
        adjust the `excluded_bond_depth` parameter to excluded directly
        bonded neighbors from overlap checks.

        See Also:
        ---------
        mbuild.Compound.direct_bonds()

        Returns:
        --------
        overlapping_particles : list of tuples
            A list of particle pairs that were found within minimum_distance.
        """
        if excluded_bond_depth < 0 or not isinstance(excluded_bond_depth, int):
            raise ValueError("`excluded_bond_depth must be an integer >= 0.")

        if self.box:
            min_length = np.min(self.box.lengths)
        else:
            min_length = np.min(self.get_boundingbox().lengths)

        if minimum_distance >= min_length / 2:
            raise ValueError(
                "The minimum distance chosen is greater than or equal to "
                "half of the box length."
            )

        freud = import_("freud")
        moved_positions, freud_box = self.to_freud()
        aq = freud.locality.AABBQuery(freud_box, moved_positions)
        aq_query = aq.query(
            query_points=moved_positions,
            query_args=dict(r_min=0.0, r_max=minimum_distance, exclude_ii=True),
        )
        nlist = aq_query.toNeighborList()
        # nlist contains each pair twice, get the set
        pairs_set = set([tuple(sorted((i, j))) for i, j in nlist])
        all_particles = [p for p in self.particles()]
        overlapping_particles = []
        for i, j in pairs_set:
            # Exclude bonded neighbors that are within min distance
            if excluded_bond_depth > 0:
                i_bonds = all_particles[i].direct_bonds(graph_depth=excluded_bond_depth)
                if all_particles[j] not in i_bonds:
                    overlapping_particles.append((i, j))
            else:  # Don't exclude bonded neighbors
                overlapping_particles.append((i, j))
        return overlapping_particles

    def get_boundingbox(self, pad_box=None):
        """Compute the bounding box of the compound.

        Compute and store the rectangular bounding box of the Compound.

        Parameters
        ----------
        pad_box: Sequence, optional, default=None
            Pad all lengths or a list of lengths by a specified amount in nm.
            Acceptable values are:

                - A single float: apply this pad value to all 3 box lengths.
                - A sequence of length 1: apply this pad value to all 3 box lengths.
                - A sequence of length 3: apply these pad values to the a, b, c box lengths.


        Returns
        -------
        mb.Box
            The bounding box for this Compound.

        Notes
        -----
        Triclinic bounding boxes are supported, but only for Compounds
        that are generated from mb.Lattice's and the resulting
        mb.Lattice.populate method
        """
        # case where only 1 particle exists
        is_one_particle = False
        if self.xyz.shape[0] == 1:
            is_one_particle = True

        # are any columns all equalivalent values?
        # an example of this would be a planar molecule
        # example: all z values are 0.0
        # from: https://stackoverflow.com/a/14860884
        # steps: create mask array comparing first value in each column
        # use np.all with axis=0 to do row columnar comparision
        has_dimension = [True, True, True]
        if not is_one_particle:
            missing_dimensions = np.all(
                np.isclose(self.xyz, self.xyz[0, :], atol=1e-2),
                axis=0,
            )
            for i, truthy in enumerate(missing_dimensions):
                has_dimension[i] = not truthy

        if is_one_particle:
            v1 = np.asarray([[1.0, 0.0, 0.0]])
            v2 = np.asarray([[0.0, 1.0, 0.0]])
            v3 = np.asarray([[0.0, 0.0, 1.0]])
        else:
            v1 = np.asarray((self.maxs[0] - self.mins[0], 0.0, 0.0))
            v2 = np.asarray((0.0, self.maxs[1] - self.mins[1], 0.0))
            v3 = np.asarray((0.0, 0.0, self.maxs[2] - self.mins[2]))
        vecs = [v1, v2, v3]

        # handle any missing dimensions (planar molecules)
        for i, dim in enumerate(has_dimension):
            if not dim:
                vecs[i][i] = 0.1

        if pad_box is not None:
            if isinstance(pad_box, (int, float, str, Sequence)):
                if isinstance(pad_box, Sequence):
                    if len(pad_box) == 1:
                        padding = [float(pad_box[0])] * 3
                    elif len(pad_box) == 3:
                        padding = [float(val) for val in pad_box]
                    else:
                        raise TypeError(
                            f"Expected a Sequence of length 1 or 3 for pad_box. Provided: {len(pad_box)}"
                        )
                else:
                    pad_box = float(pad_box)
                    padding = [pad_box] * 3
            else:
                raise TypeError(
                    f"Expected a value of type: int, float, str, or Sequence, was provided: {type(pad_box)}"
                )
            for dim, val in enumerate(padding):
                vecs[dim][dim] = vecs[dim][dim] + val

        bounding_box = Box.from_vectors(vectors=np.asarray([vecs]).reshape(3, 3))
        return bounding_box

    def min_periodic_distance(self, xyz0, xyz1):
        """Vectorized distance calculation considering minimum image.

        Only implemented for orthorhombic simulation boxes.

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
        if self.box is not None:
            if np.allclose(self.box.angles, 90.0):
                d = np.where(
                    d > 0.5 * np.array(self.box.lengths),
                    np.array(self.box.lengths) - d,
                    d,
                )
            else:
                raise NotImplementedError(
                    "Periodic distance calculation is not implemented "
                    "for non-orthorhombic boxes"
                )
        else:
            """
            raise MBuildError(f'Cannot calculate minimum periodic distance. '
                              f'No Box set for {self}')
            """
            logger.warning(
                f"No Box object set for {self}, using rectangular bounding box"
            )
            self.box = self.get_boundingbox()
            if np.allclose(self.box.angles, 90.0):
                d = np.where(
                    d > 0.5 * np.array(self.box.lengths),
                    np.array(self.box.lengths) - d,
                    d,
                )
            else:
                raise NotImplementedError(
                    "Periodic distance calculation is not implemented "
                    "for non-orthorhombic boxes"
                )
        return np.sqrt((d**2).sum(axis=-1))

    def particles_in_range(
        self,
        compound,
        dmax,
        max_particles=20,
        particle_kdtree=None,
        particle_array=None,
    ):
        """Find particles within a specified range of another particle.

        Parameters
        ----------
        compound : mb.Compound
            Reference particle to find other particles in range of
        dmax : float
            Maximum distance from 'compound' to look for Particles
        max_particles : int, optional, default=20
            Maximum number of Particles to return
        particle_kdtree : mb.PeriodicKDTree, optional
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
        periodic_kdtree.PerioidicKDTree : mBuild implementation of kd-trees
        scipy.spatial.kdtree : Further details on kd-trees
        """
        if self.box is None:
            self.box = self.get_boundingbox()
        if particle_kdtree is None:
            particle_kdtree = PeriodicKDTree.from_compound(self, leafsize=10)
        _, idxs = particle_kdtree.query(
            compound.pos, k=max_particles, distance_upper_bound=dmax
        )
        idxs = idxs[idxs != self.n_particles]
        if particle_array is None:
            particle_array = np.array(list(self.particles()))
        return particle_array[idxs]

    def visualize(
        self,
        show_ports=False,
        backend="py3dmol",
        color_scheme={},
        bead_size=0.3,
    ):  # pragma: no cover
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
        bead_size : float, Optional, default=0.3
            Size of beads in visualization
        """
        viz_pkg = {
            "nglview": self._visualize_nglview,
            "py3dmol": self._visualize_py3dmol,
        }
        if run_from_ipython():
            if backend.lower() in viz_pkg:
                return viz_pkg[backend.lower()](
                    show_ports=show_ports,
                    color_scheme=color_scheme,
                    bead_size=bead_size,
                )
            else:
                raise RuntimeError(
                    f"Unsupported visualization backend ({backend}). "
                    "Currently supported backends include nglview and py3dmol"
                )

        else:
            raise RuntimeError("Visualization is only supported in Jupyter Notebooks.")

    def _visualize_py3dmol(self, show_ports=False, color_scheme={}, bead_size=0.3):
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
        bead_size : float, Optional, default=0.3
            Size of beads in visualization

        Returns
        -------
        view : py3Dmol.view
        """
        py3Dmol = import_("py3Dmol")

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
                particle.name = "UNK"
        tmp_dir = tempfile.mkdtemp()
        cloned.save(
            os.path.join(tmp_dir, "tmp.mol2"),
            include_ports=show_ports,
            overwrite=True,
        )

        view = py3Dmol.view()
        with open(os.path.join(tmp_dir, "tmp.mol2"), "r") as f:
            view.addModel(f.read(), "mol2", keepH=True)

        view.setStyle(
            {
                "stick": {"radius": bead_size * 0.6, "color": "grey"},
                "sphere": {
                    "scale": bead_size,
                    "colorscheme": modified_color_scheme,
                },
            }
        )
        view.zoomTo()

        return view

    def _visualize_nglview(self, show_ports=False, color_scheme={}, bead_size=0.3):
        """Visualize the Compound using nglview.

        Allows for visualization of a Compound within a Jupyter Notebook.

        Parameters
        ----------
        include_ports : bool, optional, default=False
            Visualize Ports in addition to Particles
        """
        nglview = import_("nglview")
        mdtraj = import_("mdtraj")  # noqa: F841
        from mdtraj.geometry.sasa import _ATOMIC_RADII

        def remove_digits(x):
            return "".join(i for i in x if not i.isdigit() or i == "_")

        for particle in self.particles():
            particle.name = remove_digits(particle.name).upper()
            if not particle.name:
                particle.name = "UNK"
        tmp_dir = tempfile.mkdtemp()
        self.save(
            os.path.join(tmp_dir, "tmp.mol2"),
            include_ports=show_ports,
            overwrite=True,
        )
        widget = nglview.show_file(os.path.join(tmp_dir, "tmp.mol2"))
        widget.clear()
        widget.add_ball_and_stick(cylinderOnly=True)
        elements = set([particle.name for particle in self.particles()])
        scale = 50.0
        for element in elements:
            try:
                widget.add_ball_and_stick(
                    f"_{element.upper()}",
                    aspect_ratio=_ATOMIC_RADII[element.title()] ** 1.5 * scale,
                )
            except KeyError:
                ids = [
                    str(i)
                    for i, particle in enumerate(self.particles())
                    if particle.name == element
                ]
                widget.add_ball_and_stick(
                    f"@{','.join(ids)}",
                    aspect_ratio=0.17**1.5 * scale,
                    color="grey",
                )
        if show_ports:
            widget.add_ball_and_stick("_VS", aspect_ratio=1.0, color="#991f00")
        overwrite_nglview_default(widget)
        return widget

    def condense(self, inplace=True):
        """Condense the hierarchical structure of the Compound to the level of molecules.

        Modify the mBuild Compound to become a Compound with 3 distinct levels in the hierarchy.
        The top level container (self), contains molecules (i.e., connected Compounds) and the
        third level represents Particles (i.e., Compounds with no children).
        If the system contains a Particle(s) without any connections to other Compounds, it will
        appear in the 2nd level (with the top level self as a parent).

        Parameter
        ---------
        inplace : bool, optional, default=True
            Option to perform the condense operation inplace or return a copy

        Return
        ------
        self : mb.Compound or None
            return a condensed Compound if inplace is False.
        """
        # temporary list of components
        comp_list = []
        comp_list_id = []
        connected_subgraph = self.root.bond_graph.connected_components()
        unbound_particles = []
        for molecule in connected_subgraph:
            if len(molecule) == 1:
                ancestors = [molecule[0]]
                unbound_particles.append(molecule[0])
            else:
                ancestors = IndexedSet(molecule[0].ancestors())
                for particle in molecule[1:]:
                    # This works because the way in which particle.ancestors is
                    # traversed, the lower level will be in the front.
                    # The intersection will be left at the end,
                    # ancestor of the first particle is used as reference.
                    # Hence, this called will return the lowest-level Compound
                    # that is a molecule
                    ancestors = ancestors.intersection(IndexedSet(particle.ancestors()))

                """Parse molecule information"""
                molecule_tag = ancestors[0]
                comp_list.append(clone(molecule_tag))
                # generate a list of particle ids within the compound
                # this will also include any particles that are part of the
                # Compound but are not bonded

                pids = [id(p) for p in molecule_tag.particles()]
                comp_list_id += pids
        # loop over particles without bonds
        # if any of the particles exist in a compound
        # don't add them, as they already exist
        for ubp in unbound_particles:
            if id(ubp) not in comp_list_id:
                comp_list.append(clone(ubp))

        if inplace:
            for child in [self.children]:
                # Need to handle the case when child is a port
                self.remove(child)
            self.add(comp_list)
        else:
            new_compound = Compound(name=self.name)
            new_compound.add(comp_list)
            return new_compound

    def flatten(self, inplace=True):
        """Flatten the hierarchical structure of the Compound.

        Modify the mBuild Compound to become a Compound where there is
        a single container (self) that contains all the particles.

        Parameter
        ---------
        inplace : bool, optional, default=True
            Option to perform the flatten operation inplace or return a copy

        Return
        ------
        self : mb.Compound or None
            return a flattened Compound if inplace is False.
        """
        ports_list = list(self.all_ports())
        children_list = list(self.children)
        particle_list = list(self.particles())
        bond_graph = self.root.bond_graph

        # Make a list of bond that involved the particles of this compound.
        # This include bonds made exist between this compound and other
        # component of the system
        new_bonds = list()
        for particle in particle_list:
            for neighbor in nx.neighbors(bond_graph, particle):
                new_bonds.append((particle, neighbor))

        # Remove all labels which refer to children in the hierarchy
        self.labels.clear()

        # Remove all the children
        if inplace:
            for child in children_list:
                # Need to handle the case when child is a port
                self.remove(child)

            # Re-add the particles and bonds
            self.add(particle_list)
            self.add(ports_list)

            for bond in new_bonds:
                self.add_bond(bond)
        else:
            comp = clone(self)
            comp.flatten(inplace=True)
            return comp
        self.reset_labels()

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
        """Adjust port locations after particles have moved.

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
        """Slightly adjust all coordinates in a Compound.

        Provides a slight adjustment to coordinates to kick them out of local
        energy minima.
        """
        xyz_init = self.xyz
        for particle in self.particles():
            particle.pos += (np.random.rand(3) - 0.5) / 100
        self._update_port_locations(xyz_init)

    def energy_minimize(
        self,
        forcefield="UFF",
        steps=1000,
        shift_com=True,
        anchor=None,
        **kwargs,
    ):
        """Perform an energy minimization on a Compound.

        Default behavior utilizes `Open Babel <http://openbabel.org/docs/dev/>`_
        to perform an energy minimization/geometry optimization on a Compound by
        applying a generic force field

        Can also utilize `OpenMM <http://openmm.org/>`_ to energy minimize after
        atomtyping a Compound using
        `Foyer <https://github.com/mosdef-hub/foyer>`_ to apply a forcefield XML
        file that contains valid SMARTS strings.

        This function is primarily intended to be used on smaller components,
        with sizes on the order of 10's to 100's of particles, as the energy
        minimization scales poorly with the number of particles.

        Parameters
        ----------
        steps : int, optional, default=1000
            The number of optimization iterations
        forcefield : str, optional, default='UFF'
            The generic force field to apply to the Compound for minimization.
            Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', 'Ghemical'.
            Please refer to the `Open Babel documentation
            <http://open-babel.readthedocs.io/en/latest/Forcefields/Overview.html>`_
            when considering your choice of force field.
            Utilizing OpenMM for energy minimization requires a forcefield
            XML file with valid SMARTS strings. Please refer to `OpenMM docs
            <http://docs.openmm.org/7.0.0/userguide/application.html#creating-force-fields>`_
            for more information.
        shift_com : bool, optional, default=True
            If True, the energy-minimized Compound is translated such that the
            center-of-mass is unchanged relative to the initial configuration.
        anchor : Compound, optional, default=None
            Translates the energy-minimized Compound such that the
            position of the anchor Compound is unchanged relative to the
            initial configuration.


        Other Parameters
        ----------------
        algorithm : str, optional, default='cg'
            The energy minimization algorithm.  Valid options are 'steep', 'cg',
            and 'md', corresponding to steepest descent, conjugate gradient, and
            equilibrium molecular dynamics respectively.
            For _energy_minimize_openbabel
        fixed_compounds : Compound, optional, default=None
            An individual Compound or list of Compounds that will have their
            position fixed during energy minimization. Note, positions are fixed
            using a restraining potential and thus may change slightly.
            Position fixing will apply to all Particles (i.e., atoms) that exist
            in the Compound and to particles in any subsequent sub-Compounds.
            By default x,y, and z position is fixed. This can be toggled by instead
            passing a list containing the Compound and an list or tuple of bool values
            corresponding to x,y and z; e.g., [Compound, (True, True, False)]
            will fix the x and y position but allow z to be free.
            For _energy_minimize_openbabel
        ignore_compounds: Compound, optional, default=None
            An individual compound or list of Compounds whose underlying particles
            will have their positions fixed and not interact with other atoms via
            the specified force field during the energy minimization process.
            Note, a restraining potential used and thus absolute position may vary
            as a result of the energy minimization process.
            Interactions of these ignored atoms can  be specified by the user,
            e.g., by explicitly setting a distance constraint.
            For _energy_minimize_openbabel
        distance_constraints: list, optional, default=None
            A list containing a pair of Compounds as a tuple or list and
            a float value specifying the target distance between the two Compounds, e.g.,:
            [(compound1, compound2), distance].
            To specify more than one constraint, pass constraints as a 2D list, e.g.,:
            [ [(compound1, compound2), distance1],  [(compound3, compound4), distance2] ].
            Note, Compounds specified here must represent individual point particles.
            For _energy_minimize_openbabel
        constraint_factor: float, optional, default=50000.0
            Harmonic springs are used to constrain distances and fix atom positions, where
            the resulting energy associated with the spring is scaled by the
            constraint_factor; the energy of this spring is considering during the minimization.
            As such, very small values of the constraint_factor may result in an energy
            minimized state that does not adequately restrain the distance/position of atoms.
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
        constraints : str, optional, default="AllBonds"
            Specify constraints on the molecule to minimize, options are:
            None, "HBonds", "AllBonds", "HAngles"
            For _energy_minimize_openmm

        References
        ----------
        If using _energy_minimize_openmm(), please cite:

        .. [Eastman2013] P. Eastman, M. S. Friedrichs, J. D. Chodera,
           R. J. Radmer, C. M. Bruns, J. P. Ku, K. A. Beauchamp, T. J. Lane,
           L.-P. Wang, D. Shukla, T. Tye, M. Houston, T. Stich, C. Klein,
           M. R. Shirts, and V. S. Pande. "OpenMM 4: A Reusable, Extensible,
           Hardware Independent Library for High Performance Molecular
           Simulation." J. Chem. Theor. Comput. 9(1): 461-469. (2013).

        If using _energy_minimize_openbabel(), please cite:

        .. [OBoyle2011] O'Boyle, N.M.; Banck, M.; James, C.A.; Morley, C.;
           Vandermeersch, T.; Hutchison, G.R. "Open Babel: An open chemical
           toolbox." (2011) J. Cheminf. 3, 33

        .. [OpenBabel] Open Babel, version X.X.X http://openbabel.org,
           (installed Month Year)

        If using the 'MMFF94' force field please also cite the following:

        .. [Halgren1996a] T.A. Halgren, "Merck molecular force field. I. Basis,
           form, scope, parameterization, and performance of MMFF94." (1996)
           J. Comput. Chem. 17, 490-519

        .. [Halgren1996b] T.A. Halgren, "Merck molecular force field. II. MMFF94
           van der Waals and electrostatic parameters for intermolecular
           interactions." (1996) J. Comput. Chem. 17, 520-552

        .. [Halgren1996c] T.A. Halgren, "Merck molecular force field. III.
           Molecular geometries and vibrational frequencies for MMFF94." (1996)
           J. Comput. Chem. 17, 553-586

        .. [Halgren1996d] T.A. Halgren and R.B. Nachbar, "Merck molecular force
           field. IV. Conformational energies and geometries for MMFF94." (1996)
           J. Comput. Chem. 17, 587-615

        .. [Halgren1996e] T.A. Halgren, "Merck molecular force field. V.
           Extension of MMFF94 using experimental data, additional computational
           data, and empirical rules." (1996) J. Comput. Chem. 17, 616-641

        If using the 'MMFF94s' force field please cite the above along with:

        .. [Halgren1999] T.A. Halgren, "MMFF VI. MMFF94s option for energy minimization
           studies." (1999) J. Comput. Chem. 20, 720-729

        If using the 'UFF' force field please cite the following:

        .. [Rappe1992] Rappe, A.K., Casewit, C.J., Colwell, K.S., Goddard, W.A.
           III, Skiff, W.M. "UFF, a full periodic table force field for
           molecular mechanics and molecular dynamics simulations." (1992)
           J. Am. Chem. Soc. 114, 10024-10039

        If using the 'GAFF' force field please cite the following:

        .. [Wang2004] Wang, J., Wolf, R.M., Caldwell, J.W., Kollman, P.A.,
           Case, D.A. "Development and testing of a general AMBER force field"
           (2004) J. Comput. Chem. 25, 1157-1174

        If using the 'Ghemical' force field please cite the following:

        .. [Hassinen2001] T. Hassinen and M. Perakyla, "New energy terms for
           reduced protein models implemented in an off-lattice force field"
           (2001) J. Comput. Chem. 22, 1229-1242

        """
        # TODO: Update mbuild tutorials to provide overview of new features
        #   Preliminary tutorials: https://github.com/chrisiacovella/mbuild_energy_minimization
        com = self.pos
        anchor_in_compound = False
        if anchor is not None:
            # check to see if the anchor exists
            # in the Compound to be energy minimized
            for succesor in self.successors():
                if id(anchor) == id(succesor):
                    anchor_in_compound = True
                    anchor_pos_old = anchor.pos

            if not anchor_in_compound:
                raise MBuildError(
                    f"Anchor: {anchor} is not part of the Compound: {self}"
                    "that you are trying to energy minimize."
                )
        self._kick()
        extension = os.path.splitext(forcefield)[-1]
        openbabel_ffs = ["MMFF94", "MMFF94s", "UFF", "GAFF", "Ghemical"]
        if forcefield in openbabel_ffs:
            self._energy_minimize_openbabel(
                forcefield=forcefield, steps=steps, **kwargs
            )
        else:
            tmp_dir = tempfile.mkdtemp()
            self.save(os.path.join(tmp_dir, "un-minimized.mol2"))

            if extension == ".xml":
                self._energy_minimize_openmm(
                    tmp_dir,
                    forcefield_files=forcefield,
                    forcefield_name=None,
                    steps=steps,
                    **kwargs,
                )
            else:
                self._energy_minimize_openmm(
                    tmp_dir,
                    forcefield_files=None,
                    forcefield_name=forcefield,
                    steps=steps,
                    **kwargs,
                )

            self.update_coordinates(os.path.join(tmp_dir, "minimized.pdb"))

        if shift_com:
            self.translate_to(com)

        if anchor_in_compound:
            anchor_pos_new = anchor.pos
            delta = anchor_pos_old - anchor_pos_new
            self.translate(delta)

    def _energy_minimize_openmm(
        self,
        tmp_dir,
        forcefield_files=None,
        forcefield_name=None,
        steps=1000,
        scale_bonds=1,
        scale_angles=1,
        scale_torsions=1,
        scale_nonbonded=1,
        constraints="AllBonds",
    ):
        """Perform energy minimization using OpenMM.

        Converts an mBuild Compound to a ParmEd Structure,
        applies a forcefield using Foyer, and creates an OpenMM System.

        Parameters
        ----------
        forcefield_files : str or list of str, optional, default=None
            Forcefield files to load
        forcefield_name : str, optional, default=None
            Apply a named forcefield to the output file using the `foyer`
            package, e.g. 'oplsaa'. `Foyer forcefields`
            <https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields>_
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
        constraints : str, optional, default="AllBonds"
            Specify constraints on the molecule to minimize, options are:
            None, "HBonds", "AllBonds", "HAngles"

        Notes
        -----
        Assumes a particular organization for the force groups
        (HarmonicBondForce, HarmonicAngleForce, RBTorsionForce, NonBondedForce)

        References
        ----------
        [Eastman2013]_
        """
        foyer = import_("foyer")

        to_parmed = self.to_parmed()
        ff = foyer.Forcefield(forcefield_files=forcefield_files, name=forcefield_name)
        to_parmed = ff.apply(to_parmed)

        import openmm.unit as u
        from openmm.app import AllBonds, HAngles, HBonds
        from openmm.app.pdbreporter import PDBReporter
        from openmm.app.simulation import Simulation
        from openmm.openmm import LangevinIntegrator

        if constraints:
            if constraints == "AllBonds":
                constraints = AllBonds
            elif constraints == "HBonds":
                constraints = HBonds
            elif constraints == "HAngles":
                constraints = HAngles
            else:
                raise ValueError(
                    f"Provided constraints value of: {constraints}.\n"
                    f'Expected "HAngles", "AllBonds" "HBonds".'
                )
            system = to_parmed.createSystem(
                constraints=constraints
            )  # Create an OpenMM System
        else:
            system = to_parmed.createSystem()  # Create an OpenMM System
        # Create a Langenvin Integrator in OpenMM
        integrator = LangevinIntegrator(
            298 * u.kelvin, 1 / u.picosecond, 0.002 * u.picoseconds
        )
        # Create Simulation object in OpenMM
        simulation = Simulation(to_parmed.topology, system, integrator)

        # Loop through forces in OpenMM System and set parameters
        for force in system.getForces():
            if type(force).__name__ == "HarmonicBondForce":
                for bond_index in range(force.getNumBonds()):
                    atom1, atom2, r0, k = force.getBondParameters(bond_index)
                    force.setBondParameters(
                        bond_index, atom1, atom2, r0, k * scale_bonds
                    )
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "HarmonicAngleForce":
                for angle_index in range(force.getNumAngles()):
                    atom1, atom2, atom3, r0, k = force.getAngleParameters(angle_index)
                    force.setAngleParameters(
                        angle_index, atom1, atom2, atom3, r0, k * scale_angles
                    )
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "RBTorsionForce":
                for torsion_index in range(force.getNumTorsions()):
                    (
                        atom1,
                        atom2,
                        atom3,
                        atom4,
                        c0,
                        c1,
                        c2,
                        c3,
                        c4,
                        c5,
                    ) = force.getTorsionParameters(torsion_index)
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
                        c5 * scale_torsions,
                    )
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "NonbondedForce":
                for nb_index in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(nb_index)
                    force.setParticleParameters(
                        nb_index, charge, sigma, epsilon * scale_nonbonded
                    )
                force.updateParametersInContext(simulation.context)

            elif type(force).__name__ == "CMMotionRemover":
                pass

            else:
                logger.warning(
                    f"OpenMM Force {type(force).__name__} is "
                    "not currently supported in _energy_minimize_openmm. "
                    "This Force will not be updated!"
                )

        simulation.context.setPositions(to_parmed.positions)
        # Run energy minimization through OpenMM
        simulation.minimizeEnergy(maxIterations=steps)
        reporter = PDBReporter(os.path.join(tmp_dir, "minimized.pdb"), 1)
        reporter.report(simulation, simulation.context.getState(getPositions=True))

    def _check_openbabel_constraints(
        self,
        particle_list,
        successors_list,
        check_if_particle=False,
    ):
        """Provide routines commonly used to check constraint inputs."""
        for part in particle_list:
            if not isinstance(part, Compound):
                raise MBuildError(f"{part} is not a Compound.")
            if id(part) != id(self) and id(part) not in successors_list:
                raise MBuildError(f"{part} is not a member of Compound {self}.")

            if check_if_particle:
                if len(part.children) != 0:
                    raise MBuildError(
                        f"{part} does not correspond to an individual particle."
                    )

    def _energy_minimize_openbabel(
        self,
        steps=1000,
        algorithm="cg",
        forcefield="UFF",
        constraint_factor=50000.0,
        distance_constraints=None,
        fixed_compounds=None,
        ignore_compounds=None,
    ):
        """Perform an energy minimization on a Compound.

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
            Valid options are 'MMFF94', 'MMFF94s', ''UFF', 'GAFF', 'Ghemical'.
            Please refer to the Open Babel documentation
            (http://open-babel.readthedocs.io/en/latest/Forcefields/Overview.html)
            when considering your choice of force field.
        fixed_compounds : Compound, optional, default=None
            An individual Compound or list of Compounds that will have their
            position fixed during energy minimization. Note, positions are fixed
            using a restraining potential and thus may change slightly.
            Position fixing will apply to all Particles (i.e., atoms) that exist
            in the Compound and to particles in any subsequent sub-Compounds.
            By default x,y, and z position is fixed. This can be toggled by instead
            passing a list containing the Compound and a list or tuple of bool values
            corresponding to x,y and z; e.g., [Compound, (True, True, False)]
            will fix the x and y position but allow z to be free.
        ignore_compounds: Compound, optional, default=None
            An individual compound or list of Compounds whose underlying particles
            will have their positions fixed and not interact with other atoms via
            the specified force field during the energy minimization process.
            Note, a restraining potential is used and thus absolute position may vary
            as a result of the energy minimization process.
            Interactions of these ignored atoms can  be specified by the user,
            e.g., by explicitly setting a distance constraint.
        distance_constraints: list, optional, default=None
            A list containing a pair of Compounds as a tuple or list and
            a float value specifying the target distance between the two Compounds, e.g.,:
            [(compound1, compound2), distance].
            To specify more than one constraint, pass constraints as a 2D list, e.g.,:
            [ [(compound1, compound2), distance1],  [(compound3, compound4), distance2] ].
            Note, Compounds specified here must represent individual point particles.
        constraint_factor: float, optional, default=50000.0
            Harmonic springs are used to constrain distances and fix atom positions, where
            the resulting energy associated with the spring is scaled by the
            constraint_factor; the energy of this spring is considering during the minimization.
            As such, very small values of the constraint_factor may result in an energy
            minimized state that does not adequately restrain the distance/position of atom(s)e.


        References
        ----------
        [OBoyle2011]_
        [OpenBabel]_

        If using the 'MMFF94' force field please also cite the following:
        [Halgren1996a]_
        [Halgren1996b]_
        [Halgren1996c]_
        [Halgren1996d]_
        [Halgren1996e]_

        If using the 'MMFF94s' force field please cite the above along with:
        [Halgren1999]_

        If using the 'UFF' force field please cite the following:
        [Rappe1992]_

        If using the 'GAFF' force field please cite the following:
        [Wang2001]_

        If using the 'Ghemical' force field please cite the following:
        [Hassinen2001]_
        """
        openbabel = import_("openbabel")
        for particle in self.particles():
            if particle.element is None:
                try:
                    particle._element = element_from_symbol(particle.name)
                except ElementError:
                    try:
                        particle._element = element_from_name(particle.name)
                    except ElementError:
                        raise MBuildError(
                            f"No element assigned to {particle}; element could not be"
                            f"inferred from particle name {particle.name}. Cannot perform"
                            "an energy minimization."
                        )
        # Create a dict containing particle id and associated index to speed up looping
        particle_idx = {
            id(particle): idx for idx, particle in enumerate(self.particles())
        }

        # A list containing all Compounds ids contained in self. Will be used to check if
        # compounds refered to in the constrains are actually in the Compound we are minimizing.
        successors_list = [id(compound) for compound in self.successors()]

        # initialize constraints
        ob_constraints = openbabel.OBFFConstraints()

        if distance_constraints is not None:
            # if a user passes single constraint as a 1-D array,
            # i.e., [(p1,p2), 2.0]  rather than [[(p1,p2), 2.0]],
            # just add it to a list so we can use the same looping code
            if len(np.array(distance_constraints, dtype=object).shape) == 1:
                distance_constraints = [distance_constraints]

            for con_temp in distance_constraints:
                p1 = con_temp[0][0]
                p2 = con_temp[0][1]

                self._check_openbabel_constraints(
                    [p1, p2], successors_list, check_if_particle=True
                )
                if id(p1) == id(p2):
                    raise MBuildError(
                        f"Cannot create a constraint between a Particle and itself: {p1} {p2} ."
                    )

                # openbabel indices start at 1
                pid_1 = particle_idx[id(p1)] + 1
                # openbabel indices start at 1
                pid_2 = particle_idx[id(p2)] + 1
                dist = (
                    con_temp[1] * 10.0
                )  # obenbabel uses angstroms, not nm, convert to angstroms

                ob_constraints.AddDistanceConstraint(pid_1, pid_2, dist)

        if fixed_compounds is not None:
            # if we are just passed a single Compound, wrap it into
            # and array so we can just use the same looping code
            if isinstance(fixed_compounds, Compound):
                fixed_compounds = [fixed_compounds]

            # if fixed_compounds is a 1-d array and it is of length 2, we need to determine whether it is
            # a list of two Compounds or if fixed_compounds[1] should correspond to the directions to constrain
            if len(np.array(fixed_compounds, dtype=object).shape) == 1:
                if len(fixed_compounds) == 2:
                    if not isinstance(fixed_compounds[1], Compound):
                        # if it is not a list of two Compounds, make a 2d array so we can use the same looping code
                        fixed_compounds = [fixed_compounds]

            for fixed_temp in fixed_compounds:
                # if an individual entry is a list, validate the input
                if isinstance(fixed_temp, list):
                    if len(fixed_temp) == 2:
                        msg1 = (
                            "Expected tuple or list of length 3 to set"
                            "which dimensions to fix motion."
                        )
                        assert isinstance(fixed_temp[1], (list, tuple)), msg1

                        msg2 = (
                            "Expected tuple or list of length 3 to set"
                            "which dimensions to fix motion, "
                            f"{len(fixed_temp[1])} found."
                        )
                        assert len(fixed_temp[1]) == 3, msg2

                        dims = [dim for dim in fixed_temp[1]]
                        msg3 = (
                            "Expected bool values for which directions are fixed."
                            f"Found instead {dims}."
                        )
                        assert all(isinstance(dim, bool) for dim in dims), msg3

                        p1 = fixed_temp[0]

                    # if fixed_compounds is defined as [[Compound],[Compound]],
                    # fixed_temp will be a list of length 1
                    elif len(fixed_temp) == 1:
                        p1 = fixed_temp[0]
                        dims = [True, True, True]

                else:
                    p1 = fixed_temp
                    dims = [True, True, True]

                all_true = all(dims)

                self._check_openbabel_constraints([p1], successors_list)

                if len(p1.children) == 0:
                    pid = particle_idx[id(p1)] + 1  # openbabel indices start at 1

                    if all_true:
                        ob_constraints.AddAtomConstraint(pid)
                    else:
                        if dims[0]:
                            ob_constraints.AddAtomXConstraint(pid)
                        if dims[1]:
                            ob_constraints.AddAtomYConstraint(pid)
                        if dims[2]:
                            ob_constraints.AddAtomZConstraint(pid)
                else:
                    for particle in p1.particles():
                        pid = (
                            particle_idx[id(particle)] + 1
                        )  # openbabel indices start at 1

                        if all_true:
                            ob_constraints.AddAtomConstraint(pid)
                        else:
                            if dims[0]:
                                ob_constraints.AddAtomXConstraint(pid)
                            if dims[1]:
                                ob_constraints.AddAtomYConstraint(pid)
                            if dims[2]:
                                ob_constraints.AddAtomZConstraint(pid)

        if ignore_compounds is not None:
            temp1 = np.array(ignore_compounds, dtype=object)
            if len(temp1.shape) == 2:
                ignore_compounds = list(temp1.reshape(-1))

            # Since the ignore_compounds can only be passed as a list
            # we can check the whole list at once before looping over it
            self._check_openbabel_constraints(ignore_compounds, successors_list)

            for ignore in ignore_compounds:
                p1 = ignore
                if len(p1.children) == 0:
                    pid = particle_idx[id(p1)] + 1  # openbabel indices start at 1
                    ob_constraints.AddIgnore(pid)

                else:
                    for particle in p1.particles():
                        pid = (
                            particle_idx[id(particle)] + 1
                        )  # openbabel indices start at 1
                        ob_constraints.AddIgnore(pid)

        mol = self.to_pybel()
        mol = mol.OBMol

        mol.PerceiveBondOrders()
        mol.SetAtomTypesPerceived()

        ff = openbabel.OBForceField.FindForceField(forcefield)
        if ff is None:
            raise MBuildError(
                f"Force field '{forcefield}' not supported for energy "
                "minimization. Valid force fields are 'MMFF94', "
                "'MMFF94s', 'UFF', 'GAFF', and 'Ghemical'."
                ""
            )
        logger.info(
            "Performing energy minimization using the Open Babel package. "
            "Please refer to the documentation to find the appropriate "
            f"citations for Open Babel and the {forcefield} force field"
        )

        if (
            distance_constraints is not None
            or fixed_compounds is not None
            or ignore_compounds is not None
        ):
            ob_constraints.SetFactor(constraint_factor)
            if ff.Setup(mol, ob_constraints) == 0:
                raise MBuildError(
                    "Could not setup forcefield for OpenBabel Optimization."
                )
        else:
            if ff.Setup(mol) == 0:
                raise MBuildError(
                    "Could not setup forcefield for OpenBabel Optimization."
                )

        if algorithm == "steep":
            ff.SteepestDescent(steps)
        elif algorithm == "md":
            ff.MolecularDynamicsTakeNSteps(steps, 300)
        elif algorithm == "cg":
            ff.ConjugateGradients(steps)
        else:
            raise MBuildError(
                "Invalid minimization algorithm. Valid options "
                "are 'steep', 'cg', and 'md'."
            )
        ff.UpdateCoordinates(mol)

        # update the coordinates in the Compound
        for i, obatom in enumerate(openbabel.OBMolAtomIter(mol)):
            x = obatom.GetX() / 10.0
            y = obatom.GetY() / 10.0
            z = obatom.GetZ() / 10.0
            self[i].pos = np.array([x, y, z])

    def save(
        self,
        filename,
        include_ports=False,
        box=None,
        overwrite=False,
        residues=None,
        **kwargs,
    ):
        """Save the Compound to a file.

        Parameters
        ----------
        filename : str
            Filesystem path in which to save the trajectory. The extension or
            prefix will be parsed and control the format. Supported extensions:
            'gsd', 'gro', 'top', 'mcf', 'pdb', 'xyz',
            'json', 'mol2', 'sdf', 'psf'. See `mbuild.conversion.save()`
            for more information about writer methods.
        include_ports : bool, optional, default=False
            Save ports contained within the compound.
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be written to the output file. If 'None', a
            bounding box is used with 0.25nm buffers at each face to avoid
            overlapping atoms.
        overwrite : bool, optional, default=False
            Overwrite if the filename already exists
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        **kwargs
            See `mbuild.conversion.save()`.
            Depending on the file extension these will be passed to either
            Parmed or GMSO backend writers
            See https://parmed.github.io/ParmEd/html/structobj/parmed.structure.
            Structure.html#parmed.structure.Structure.save and
            https://github.com/mosdef-hub/gmso/tree/main/gmso/formats

        Other Parameters
        ----------------
        ref_distance : float, optional, default=1.0
            Normalization factor used when saving to the .gsd format
            for converting distance values to reduced units.
        ref_energy : float, optional, default=1.0
            Normalization factor used when saving to the .gsd format
            for converting energy values to reduced units.
        ref_mass : float, optional, default=1.0
            Normalization factor used when saving to the .gsd format
            for converting mass values to reduced units.

        Notes
        -----
        When saving the compound as a json, only the following arguments are
        used:
        * filename
        * include_ports

        The savers used for each supported file type are:
        GMSO: .gro, .gsd, .data, .xyz, .mcf, .top
        Parmed: .mol2, .pdb, .prmtop, .cif, .crd
        PyBel: .sdf

        See Also
        --------
        mbuild.conversion.save : Main saver logic
        mbuild.formats.cassandramcf.write_mcf : Write to Cassandra MCF format
        mbuild.formats.json_formats.compound_to_json : Write to a json file
        """
        conversion.save(
            compound=self,
            filename=filename,
            include_ports=include_ports,
            box=box,
            overwrite=overwrite,
            residues=residues,
            **kwargs,
        )

    def translate(self, by):
        """Translate the Compound by a vector.

        Parameters
        ----------
        by : np.ndarray, shape=(3,), dtype=float
        """
        new_positions = _translate(self.xyz_with_ports, by)
        self.xyz_with_ports = new_positions

    def translate_to(self, pos):
        """Translate the Compound to a specific position.

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

    def spin(self, theta, around, anchor=None):
        """Rotate Compound in place around an arbitrary vector.

        Parameters
        ----------
        theta : float
            The angle by which to rotate the Compound, in radians.
        around : np.ndarray, shape=(3,), dtype=float
            The axis about which to spin the Compound.
        anchor : mb.Compound, optional, default=None (self)
            Anchor compound/particle to perform spinning.
            If the anchor is not a particle, the spin will be
            around the center of the anchor Compound.
        """
        around = np.asarray(around).reshape(3)

        if anchor:
            msg = f"{anchor} is not part of {self}."
            assert anchor in self.successors(), msg
        else:
            anchor = self
        anchor_pos = anchor.center

        self.translate(-anchor_pos)
        self.rotate(theta, around)
        self.translate(anchor_pos)

    def rotate_dihedral(self, bond, phi):
        """Rotate a dihedral about a central bond.

        Parameters
        ----------
        bond : indexable object, length=2, dtype=mb.Compound
            The pair of bonded Particles in the central bond of the dihedral
        phi : float
            The angle by which to rotate the dihedral, in radians.
        """
        nx = import_("networkx")

        # Generate a bond graph and convert to networkX
        mb_bondgraph = self.bond_graph
        G = nx.Graph(mb_bondgraph.edges())

        # Remove separate the compound in to two pieces by removing the bond
        G.remove_edge(*bond)
        assert len([i for i in nx.connected_components(G)]) == 2
        components = [G.subgraph(c).copy() for c in nx.connected_components(G)]
        component1 = components[1]  # One piece of the compound

        # Get original coordinates
        original_bond_positions = [bond[0].pos, bond[1].pos]

        # Get the vector along the bond
        bond_vec = bond[1].pos - bond[0].pos

        # Rotate the coordinates of the piece by phi about the bond vector
        xyz = np.array([p.pos for p in component1.nodes])
        transformed_xyz = _rotate(xyz, phi, bond_vec)
        for atom, coord in zip(component1.nodes, transformed_xyz):
            atom.translate_to(coord)

        # Move atoms involved in the bond to original positions
        # This is neccessary since the piece is rotated about its center
        if bond[0] in set(component1.nodes):
            trans_vec = original_bond_positions[0] - bond[0].pos
        elif bond[1] in set(component1.nodes):
            trans_vec = original_bond_positions[1] - bond[1].pos

        for atom in component1.nodes:
            atom.translate(trans_vec)

    # Interface to GMSO Topology for reading/writing mol2 files
    def from_gmso(self, topology, coords_only=False, infer_hierarchy=True):
        """Convert a GMSO Topology to mBuild Compound.

        Parameter
        ---------
        topology : gmso.Topology
            The GMSO Topology to be converted.
        compound : mb.Compound, optional, default=None
            Host mb.Compound that we are loading to.
        coords_only : bool, optional, default=False
            Set preexisting atoms in compound to coordinates given by Topology.
        infer_hierarchy : bool, optional, default=True
            If True, infer compound hierarchy from Topology residue, to be implemented.
            Pending new GMSO release.

        Returns
        -------
        compound : mb.Compound
        """
        return conversion.from_gmso(
            topology=topology,
            compound=self,
            coords_only=coords_only,
            # infer_hierarchy=infer_hierarchy,
            # TO DO: enable this with new release of GMSO
        )

    def to_gmso(self, **kwargs):
        """Create a GMSO Topology from a mBuild Compound.

        Parameters
        ----------
        compound : mb.Compound
            The mb.Compound to be converted.

        Returns
        -------
        topology : gmso.Topology
            The converted gmso Topology
        """
        return conversion.to_gmso(self, **kwargs)

    def to_hoomdsnapshot(self, **kwargs):
        """Create a HOOMD-Blue snapshot from an mBuild Compound.

        Parameters
        ----------
        compound : mb.Compound
            The mb.Compound to be converted.

        Returns
        -------
        snapshot : gsd.hoomd.Frame
            HOOMD-Blue compatible topology.
        """
        return conversion.to_hoomdsnapshot(self, **kwargs)

    def to_freud(self):
        return conversion.to_freud(self)

    # Interface to Trajectory for reading/writing .pdb and .mol2 files.
    # -----------------------------------------------------------------
    def from_trajectory(self, traj, frame=-1, coords_only=False, infer_hierarchy=True):
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

        See Also
        --------
        mbuild.conversion.from_trajectory
        """
        conversion.from_trajectory(
            traj=traj,
            compound=self,
            frame=frame,
            coords_only=coords_only,
            infer_hierarchy=True,
        )

    def to_trajectory(self, include_ports=False, chains=None, residues=None, box=None):
        """Convert to an md.Trajectory and flatten the compound.

        Parameters
        ----------
        include_ports : bool, optional, default=False
            Include all port atoms when converting to trajectory.
        chains : mb.Compound or list of mb.Compound
            Chain types to add to the topology
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be used when converting to a `Trajectory`.
            If 'None', self.box is used. If self.box is None,
            a bounding box is used with a 0.5 nm buffer in each
            dimension to avoid overlapping atoms.

        Returns
        -------
        trajectory : md.Trajectory

        See Also
        --------
        _to_topology
        """
        return conversion.to_trajectory(
            compound=self,
            include_ports=include_ports,
            chains=chains,
            residues=residues,
            box=box,
        )

    def from_parmed(self, structure, coords_only=False, infer_hierarchy=True):
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
        conversion.from_parmed(
            structure=structure,
            compound=self,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy,
        )

    def to_rdkit(self):
        """Create an RDKit RWMol from an mBuild Compound.

        Returns
        -------
        rdkit.Chem.RWmol

        Notes
        -----
        Use this method to utilzie rdkit funcitonality.
        This method only works when the mBuild compound
        contains accurate element information. As a result,
        this method is not compatible with compounds containing
        abstract particles (e.g. coarse-grained systems)

        Example
        -------
        >>> import mbuild
        >>> from rdkit.Chem import Draw

        >>> benzene = mb.load("c1cccc1", smiles=True)
        >>> benzene_rdkmol = benzene.to_rdkit()
        >>> img = Draw.MolToImage(benzene_rdkmol)

        See https://www.rdkit.org/docs/GettingStartedInPython.html

        """
        return conversion.to_rdkit(self)

    def to_parmed(
        self,
        box=None,
        title="",
        residues=None,
        include_ports=False,
        infer_residues=False,
        infer_residues_kwargs={},
    ):
        """Create a ParmEd Structure from a Compound.

        Parameters
        ----------
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be used when converting to a `Structure`.
            If 'None', self.box is used. If self.box is None,
            a bounding box is used with 0.5 nm buffer in each dimension
            to avoid overlapping atoms.
        title : str, optional, default=self.name
            Title/name of the ParmEd Structure
        residues : str of list of str, optional, default=None
            Labels of residues in the Compound. Residues are assigned by checking
            against Compound.name.
        include_ports : boolean, optional, default=False
            Include all port atoms when converting to a `Structure`.
        infer_residues : bool, optional, default=True
            Attempt to assign residues based on the number of bonds and particles in
            an object. This option is not used if `residues == None`
        infer_residues_kwargs : dict, optional, default={}
            Keyword arguments for :func:`mbuild.conversion.pull_residues`

        Returns
        -------
        parmed.structure.Structure
            ParmEd Structure object converted from self

        See Also
        --------
        mbuild.conversion.to_parmed
        parmed.structure.Structure : Details on the ParmEd Structure object
        """
        return conversion.to_parmed(
            compound=self,
            box=box,
            title=title,
            residues=residues,
            include_ports=include_ports,
            infer_residues=infer_residues,
            infer_residues_kwargs=infer_residues_kwargs,
        )

    def to_networkx(self, names_only=False):
        """Create a NetworkX graph representing the hierarchy of a Compound.

        Parameters
        ----------
        names_only : bool, optional, default=False
            Store only the names of the compounds in the graph, appended with
            their IDs, for distinction even if they have the same name. When
            set to False, the default behavior, the nodes are the compounds
            themselves.

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

    def to_pybel(
        self,
        box=None,
        title="",
        residues=None,
        include_ports=False,
        infer_residues=False,
    ):
        """Create a pybel.Molecule from a Compound.

        Parameters
        ----------
        box : mb.Box, def None
        title : str, optional, default=self.name
            Title/name of the ParmEd Structure
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        include_ports : boolean, optional, default=False
            Include all port atoms when converting to a `Structure`.
        infer_residues : bool, optional, default=False
            Attempt to assign residues based on names of children

        Returns
        -------
        pybel.Molecule

        See Also
        --------
        mbuild.conversion.to_pybel

        Notes
        -----
        Most of the mb.Compound is first converted to openbabel.OBMol
        And then pybel creates a pybel.Molecule from the OBMol
        Bond orders are assumed to be 1
        OBMol atom indexing starts at 1, with spatial dimension Angstrom
        """
        return conversion.to_pybel(
            compound=self,
            box=box,
            title=title,
            residues=residues,
            include_ports=include_ports,
        )

    def to_smiles(self, backend="pybel"):
        """Create a SMILES string from an mbuild compound.

        Parameters
        ----------
        compound : mb.Compound.
            The mbuild compound to be converted.
        backend : str, optional, default="pybel"
            Backend used to do the conversion.

        Return
        ------
        smiles_string : str
        """
        return conversion.to_smiles(self, backend)

    def from_pybel(
        self,
        pybel_mol,
        use_element=True,
        coords_only=False,
        infer_hierarchy=True,
        ignore_box_warn=False,
    ):
        """Create a Compound from a Pybel.Molecule.

        Parameters
        ----------
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

        See Also
        --------
        mbuild.conversion.from_pybel
        """
        conversion.from_pybel(
            pybel_mol=pybel_mol,
            compound=self,
            use_element=use_element,
            coords_only=coords_only,
            ignore_box_warn=ignore_box_warn,
        )

    def to_intermol(self, molecule_types=None):  # pragma: no cover
        """Create an InterMol system from a Compound.

        Parameters
        ----------
        molecule_types : list or tuple of subclasses of Compound

        Returns
        -------
        intermol_system : intermol.system.System

        See Also
        --------
        mbuild.conversion.to_intermol
        """
        return conversion.to_intermol(compound=self, molecule_types=molecule_types)

    def get_smiles(self):
        """Get SMILES string for compound.

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
        """Get item from Compound."""
        if isinstance(selection, int):
            return list(self.particles())[selection]
        if isinstance(selection, str):
            if selection not in self.labels:
                raise MBuildError(f"{self.name}['{selection}'] does not exist.")
            return self.labels.get(selection)

    def __repr__(self):
        """Compound representation."""
        descr = list("<")
        descr.append(self.name + " ")

        if self.children:
            descr.append(f"{self.n_particles} particles, ")
            descr.append(f"{self.n_bonds} bonds, ")
            if self.box is not None:
                descr.append(f"System box: {self.box}, ")
            else:
                descr.append("non-periodic, ")
        else:
            descr.append(f"pos=({np.array2string(self.pos, precision=4)}), ")
            descr.append(f"{self.n_direct_bonds} bonds, ")

        descr.append(f"id: {id(self)}>")
        return "".join(descr)

    def _clone(self, clone_of=None, root_container=None):
        """Clones compound faster than deepcopying.

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
        newone._pos = deepcopy(self._pos)
        newone.port_particle = deepcopy(self.port_particle)
        newone._box = deepcopy(self._box)
        newone._periodicity = deepcopy(self._periodicity)
        newone._charge = deepcopy(self._charge)
        newone._mass = deepcopy(self._mass)
        if hasattr(self, "index"):
            newone.index = deepcopy(self.index)

        if self.children is None:
            newone.children = None
        else:
            newone.children = list()
        # Parent should be None initially.
        newone.parent = None
        newone.labels = OrderedDict()
        newone.referrers = set()
        newone.bond_graph = None

        # Add children to clone.
        if self.children:
            for child in self.children:
                newchild = child._clone(clone_of, root_container)
                newone.children.append(newchild)
                newchild.parent = newone

        # Copy labels, except bonds with atoms outside the hierarchy.
        if self.labels:
            for label, compound in self.labels.items():
                if not isinstance(compound, list):
                    newone.labels[label] = compound._clone(clone_of, root_container)
                    compound.referrers.add(clone_of[compound])
                else:
                    # compound is a list of compounds, so we create an empty
                    # list, and add the clones of the original list elements.
                    newone.labels[label] = []
                    for subpart in compound:
                        newone.labels[label].append(
                            subpart._clone(clone_of, root_container)
                        )
                        # Referrers must have been handled already, or the will
                        # be handled

        return newone

    def _clone_bonds(self, clone_of=None):
        """Clone the bond of the source compound to clone compound."""
        newone = clone_of[self]
        newone.bond_graph = BondGraph()
        for particle in self.particles():
            newone.bond_graph.add_node(clone_of[particle])
        for c1, c2, data in self.bonds(return_bond_order=True):
            try:
                # bond order is added to the data dictionary as 'bo'
                newone.add_bond(
                    (clone_of[c1], clone_of[c2]), bond_order=data["bond_order"]
                )
            except KeyError:
                raise MBuildError(
                    "Cloning failed. Compound contains bonds to "
                    "Particles outside of its containment hierarchy."
                )


Particle = Compound


def _flatten_list(c_list):
    """Flatten a list.

    Helper function to flatten a list that may be nested, e.g. [comp1, [comp2, comp3]].
    """
    if isinstance(c_list, Iterable) and not isinstance(c_list, str):
        for c in c_list:
            if isinstance(c, Iterable) and not isinstance(c, str):
                yield from _flatten_list(c)
            else:
                yield c
