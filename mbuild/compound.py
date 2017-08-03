from __future__ import print_function, division

__all__ = ['load', 'clone', 'Compound', 'Particle']

import collections
from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
from itertools import product
import os
import sys
import tempfile
from warnings import warn

import mdtraj as md
import numpy as np
from oset import oset as OrderedSet
import parmed as pmd
from parmed.periodic_table import AtomicNum, element_by_name, Mass
import simtk.openmm.app.element as elem
from six import integer_types, string_types

from mbuild.bond_graph import BondGraph
from mbuild.box import Box
from mbuild.exceptions import MBuildError
from mbuild.formats.hoomdxml import write_hoomdxml
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.formats.gsdwriter import write_gsd
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.utils.io import run_from_ipython, import_
from mbuild.coordinate_transform import _translate, _rotate, \
    normalized_matrix, unit_vector, angle


def load(filename, relative_to_module=None, compound=None, coords_only=False,
         rigid=False, use_parmed=False, **kwargs):
    """Load a file into an mbuild compound.

    Files are read using the MDTraj package unless the `use_parmed` argument is
    specified as True. Please refer to http://mdtraj.org/1.8.0/load_functions.html
    for formats supported by MDTraj and https://parmed.github.io/ParmEd/html/
    readwrite.html for formats supported by ParmEd.

    Parameters
    ----------
    filename : str
        Name of the file from which to load atom and bond information.
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory, look for the file
        where this module is defined. This is typically used in Compound
        classes that will be instantiated from a different directory
        (such as the Compounds located in mbuild.lib).
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information into.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compoint.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    use_parmed : bool, optional, default=False
        Use readers from ParmEd instead of MDTraj.
    **kwargs : keyword arguments
        Key word arguments passed to mdTraj for loading.

    Returns
    -------
    compound : mb.Compound

    """
    # Handle mbuild *.py files containing a class that wraps a structure file
    # in its own folder. E.g., you build a system from ~/foo.py and it imports
    # from ~/bar/baz.py where baz.py loads ~/bar/baz.pdb.
    if relative_to_module:
        script_path = os.path.realpath(sys.modules[relative_to_module].__file__)
        file_dir = os.path.dirname(script_path)
        filename = os.path.join(file_dir, filename)

    if compound is None:
        compound = Compound()

    if use_parmed:
        warn("use_parmed set to True.  Bonds may be inferred from inter-particle "
             "distances and standard residue templates!")
        structure = pmd.load_file(filename, structure=True, **kwargs)
        compound.from_parmed(structure, coords_only=coords_only)
    else:
        traj = md.load(filename, **kwargs)
        compound.from_trajectory(traj, frame=-1, coords_only=coords_only)

    if rigid:
        compound.label_rigid_bodies()
    return compound


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
    boundingbox
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
                 periodicity=None, port_particle=False):
        super(Compound, self).__init__()

        if name:
            if not isinstance(name, string_types):
                raise ValueError('Compound.name should be a string. You passed '
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

        # self.add() must be called after labels and children are initialized.
        if subcompounds:
            if charge:
                raise MBuildError('Cannot set the charge of a Compound containing '
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
    def my_label(self):
       """Returns the label of the current compound, as seen by the compound's
       parent.

       The default MBuild labeling convention when building compounds
       is label = "name[{}]".format(ii) where ii is the order (zero indexed)
       which that kind of that compound/particle is added.
       """
       if self.parent:
           for lab in self.parent.labels:
               if isinstance(lab, list):
                   continue
               if self.parent[lab] is self:
                   return lab
                   break
           else:
               raise AttributeError("Developer Error")
               # revisit this error an also should I make this a property?
       else:
           warn ("Object {} is at the top of its hierarchy and thus has no label."
                 " Returning None.".format(self))
           return None


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

    def find_particle_in_path(self, within_path):
       """"
       Yields all particles that exist within the hierarchal pathway description provided
        in the parameter 'within_path'.

       A hierarchal pathway is a list or tuple containing any combination of strings,
       list/tuples, or mb.Compounds. Each element of the list/tuple either describes a series
       of subcompounds (this occurs in the instance where an inner list/tuple is passed), or describes
       one subcompound or type of subcompound (when a string is passed), or even IS a subcompound
       (in the instance where a mb.Compound object is passed). Strings correspond
       to either the names or labels of subcompounds, and list/tuples hold multiple strings that
       correspond to names and/or labels. They are used when the user wishes to describe multiple
       pathways, for example, path = [..., ["subcompound[1]", "subcompound[4]"], ...].
       This example demonstrates that the user can describe some but not all of the pathways
       that have the name "subcompound". The order of the elements in the outer list/tuple correspond
       to their position in the hierarchal pathway, where the first index is the lowest level and the
       last is the highest specified. The number of subcompounds the user can specify is unlimited so
       long as each subcompound specified lies within the hierarchy of the list/tuple element that
       follows. In the context of this function, the first index correspond to a mb.Particle.

       The idea of a pathway is similar to how one sorts through directories on a computer,
       i.e. "C:/user/username/documents" BUT since MBuild uses hierarchal pathways from lowest
       to highest, the MBuild style of writing it would be "documents/username/user/C:".
       # best hierarchal description

       EX: path =["target",
                         ["SubSubSubCompound[0]",
                          "SubSubSubCompound[3]",
                          "SubSubSubCompound[4]"],
                        "SubSubCompound",
                        "SubCompound[6]"]

       :param within_path: accepts list, tuple, or mb.particle
           If a mb.particle is provided, the function will return within_path.
           The a list/tuple is specified, each element is either a str, mb.Compound,
           list/tuple (containing any combination of strs and mb.Compounds). If
           a mb.Compound is provided the function skips to that index and ignores
           all data beyond that index. For example, given
           ["a0", "a1", <mb.Compound>, "a3"], the function would ignore the value "a3"
            and skip straight to the mb.Compound. See
           above description for more information on pathways.

       :yields
           All particles that match the hierarchal description provided.

       """
       # revisit bc of the list option. i'm worried it won't return any errors. maybe track which ones
       # it doesnt find. That would probably need to be done in find_subc_in_path
       if not isinstance(within_path, (list,tuple)):
           if not isinstance(within_path, mb.Particle):
               raise TypeError("within_path must be of type list or tuple. "
                               "User passed type: {}. within_path can also "
                               "accept a mb.Particle but in this instance "
                               "find_particle_in_path just returns within_path"
                               ".".format(type(within_path)))
           else:
               yield within_path
               return
       if not within_path:
           raise ValueError("within_path cannot be empty.")
       within_path = list(within_path)
       no_yield = True
       parti = within_path[0]
       if len(within_path) ==1:
           if isinstance(parti,mb.Particle):
               yield parti
               return
           elif isinstance(parti, str):
               for parts in self:
                   if parts.name == parti or parts.my_label == parti:
                       if no_yield:
                           no_yield = False
                       yield parts
           elif isinstance(parti, (list, tuple)):
               for parts in self:
                   if parts.name in parti or parts.my_label in parti:
                       if no_yield:
                           no_yield = False
                       yield parts
           else:
               raise TypeError("The object contained in within_path is not of "
                               "an acceptable type. Acceptable types are str, tuple,"
                               "list, mb.Particle (for only the first index), "
                               "and mb.Compound (valid for any index except the first)."
                               " User passed object of type: {}.".format(type(parti)))
       elif isinstance(parti, (list,tuple)):
           for subc in self.find_subcompounds_in_path(pathway=within_path[1:]):
               if subc:
                   for parts in subc:
                       if parts.name in parti or parts.my_label in parti:
                           if no_yield:
                               no_yield = False
                           yield parts
       else:
           for subc in self.find_subcompounds_in_path(pathway= within_path[1:]):
               if subc:
                   for parts in subc:
                       if parts.name == parti or parts.my_label == parti:
                           if no_yield:
                               no_yield = False
                           yield parts
       if no_yield:
           raise ValueError("Particle in path {} not found. Verify that "
                            "this is the correct path.".format(within_path))


    def subcompounds_by_name_or_label(self, looking_for):
       """
       Yields all the subcompounds that exhibit the specified name or label.

       Parameters:

       looking_for: accepts str
           This string will specify the name or label of the particle(s) the user wishes
           to find.

       Editors note:
       Whenever calling this function within a function make sure to add in a method to track
       if anything in looking_for was not found
       """
       if isinstance(looking_for, str):
           for parti in self.children:
               if parti.name == looking_for or parti.my_label == looking_for:
                   if parti.n_particles > 1:
                       yield parti
                   else:
                       raise ValueError("The user passed {}, the name of an atom/ particle within this "
                                       "object. Please use particles_by_name, find_particle_in_path,"
                                        "or a similar method instead.".format(parti.name))
               else:
                   if parti.n_particles > 1:
                       yield from parti.subcompounds_by_name_or_label(looking_for)
                   else:
                       yield None

       elif isinstance(looking_for, (list, tuple)) and all(looking_for):
           for l in looking_for:
               yield from self.subcompounds_by_name_or_label(looking_for=l)
       else:
           raise TypeError("looking_for must be of type str or a list/tuple of strs."
                           " User passed: {}.".format(type(looking_for)))

    def find_subcompounds_in_path(self, pathway):
       """
       yield all subcompounds that are in the specified hierarchal pathway

       :param pathway: list or tuple containing strings or mb.Compounds
           A hierarchal pathway (see below) to the desired subcompounds.

       :return: yields all particles that match the path description.
               yields None if the particle path specified doesn't exist

       A hierarchal pathway is a list or tuple containing any combination of strings,
       list/tuples, or mb.Compounds. Each element of the list/tuple either describes a series
       of subcompounds (this occurs in the instance where an inner list/tuple is passed), or describes
       one subcompound or type of subcompound (when a string is passed), or even IS a subcompound
       (in the instance where a mb.Compound object is passed). Strings correspond
       to either the names or labels of subcompounds, and list/tuples hold multiple strings that
       correspond to names and/or labels. They are used when the user wishes to describe multiple
       pathways, for example, path = [..., ["subcompound[1]", "subcompound[4]"], ...].
       This example demonstrates that the user can describe some but not all of the pathways
       that have the name "subcompound". The order of the elements in the outer list/tuple correspond
       to their position in the hierarchal pathway, where the first index is the lowest level and the
       last is the highest specified. The number of subcompounds the user can specify is unlimited so
       long as each subcompound specified lies within the hierarchy of the list/tuple element that
       follows. In the context of this function, the first index must correspond to a subcompound, not
       a mb.Particle.

       The idea of a pathway is similar to how one sorts through directories on a computer,
       i.e. "C:/user/username/documents" BUT since MBuild uses hierarchal pathways from lowest
       to highest, the MBuild style of writing it would be "documents/username/user/C:".

       EX: path =["target",
                         ["SubSubSubCompound[0]",
                          "SubSubSubCompound[3]",
                          "SubSubSubCompound[4]"],
                        "SubSubCompound",
                        "SubCompound[6]"]
       """

       if not isinstance(pathway, (list, tuple)):
           raise TypeError("Parameter pathway must be of type list or tuple. User"
                           " passed type: {}.".format(type(pathway)))
       if not pathway:
           raise ValueError("Parameter 'pathway' cannot be an empty {}.".format(type(pathway)))
       pathway = list(pathway)
       for n, ii in enumerate(pathway):
           if isinstance(ii, (str, list, tuple)):
               pass
           elif isinstance(ii, mb.Compound):
               if pathway[:n]:
                   yield from ii._which_subc(looking_for=pathway[:n])
               else:
                   yield ii
               break
           else:
               raise TypeError("pathway parameter must be either a list or tuple containing"
                               " only strings, lists, tuples or mb.Compounds. User passed {}"
                               " containing invalid type: "
                               "{} at index {}.".format(type(pathway), type(ii), n))
       else:
           yield from self._which_subc(looking_for=pathway)

    def _which_subc(self, looking_for):
       """
       refer to def find_subcompounds_in_path
       """
       shorten = len(looking_for)-1
       lf = looking_for[-1]
       if len(looking_for) > 1:
           we_ok = False
           if isinstance(lf, str):
               for subp in self.subcompounds_by_name_or_label(looking_for=lf):
                   if subp:
                       we_ok = True
                       yield from subp._which_subc(looking_for=looking_for[:shorten])
           else:
               for l in lf:
                   for subp in self.subcompounds_by_name_or_label(looking_for=l):
                       if subp:
                           we_ok = True
                           yield from subp._which_subc(looking_for=looking_for[:shorten])
           if not we_ok:
               yield None
               #raise ValueError('{} was not found within {}'.format(within, self.name))
       else:
           if isinstance(lf, str):
               yield from self.subcompounds_by_name_or_label(looking_for=lf)
           else:
               for l in lf:
                   yield from self.subcompounds_by_name_or_label(looking_for=l)

    @property
    def charge(self):
        return sum([particle._charge for particle in self.particles()])

    @charge.setter
    def charge(self, value):
        if self._contains_only_ports():
            self._charge = value
        else:
            raise AttributeError("charge is immutable for Compounds that are "
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
            raise AttributeError("rigid_id is immutable for Compounds that are "
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
            if isinstance(discrete_bodies, string_types):
                discrete_bodies = [discrete_bodies]
        if rigid_particles is not None:
            if isinstance(rigid_particles, string_types):
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
        unique_rigid_ids = sorted(set([p.rigid_id for p in self.rigid_particles()]))
        n_unique_rigid = len(unique_rigid_ids)
        if max_rigid and n_unique_rigid != max_rigid + 1:
            missing_rigid_id = (unique_rigid_ids[-1] * (unique_rigid_ids[-1] + 1))/2 - sum(unique_rigid_ids)
            for successor in self.successors():
                if successor.rigid_id is not None:
                    if successor.rigid_id > missing_rigid_id:
                        successor.rigid_id -= 1
            if self.rigid_id:
                if self.rigid_id > missing_rigid_id:
                    self.rigid_id -= 1

    def add(self, new_child, label=None, containment=True, replace=False,
            inherit_periodicity=True, reset_rigid_ids=True):
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
        reset_rigid_ids : bool, optional, default=True
            If the Compound to be added contains rigid bodies, reset the
            rigid_ids such that values remain distinct from rigid_ids
            already present in `self`. Can be set to False if attempting
            to add Compounds to an existing rigid body.

        """
        # Support batch add via lists, tuples and sets.
        if (isinstance(new_child, collections.Iterable) and
                not isinstance(new_child, string_types)):
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
            label = '{0}[$]'.format(new_child.name)

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

    def remove(self, objs_to_remove):
        """Remove children from the Compound.

        Parameters
        ----------
        objs_to_remove : mb.Compound or list of mb.Compound
            The Compound(s) to be removed from self

        """
        if not self.children:
            return

        if not hasattr(objs_to_remove, '__iter__'):
            objs_to_remove = [objs_to_remove]
        objs_to_remove = set(objs_to_remove)

        if len(objs_to_remove) == 0:
            return

        remove_from_here = objs_to_remove.intersection(self.children)
        self.children -= remove_from_here
        yet_to_remove = objs_to_remove - remove_from_here

        for removed in remove_from_here:
            for child in removed.children:
                removed.remove(child)

        for removed_part in remove_from_here:
            if removed_part.rigid_id is not None:
                for ancestor in removed_part.ancestors():
                    ancestor._check_if_contains_rigid_bodies = True
            if self.root.bond_graph and self.root.bond_graph.has_node(removed_part):
                self.root.bond_graph.remove_node(removed_part)
            self._remove_references(removed_part)

        # Remove the part recursively from sub-compounds.
        for child in self.children:
            child.remove(yet_to_remove)
            if child.contains_rigid:
                self.root._reorder_rigid_ids()

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
                return self.root.bond_graph.subgraph(self.particles()).edges_iter()
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
        particle_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
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
        if self.root.bond_graph is None or not self.root.bond_graph.has_edge(*particle_pair):
            warn("Bond between {} and {} doesn't exist!".format(*particle_pair))
            return
        self.root.bond_graph.remove_edge(*particle_pair)

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
            arr = np.fromiter(itertools.chain.from_iterable(
                particle.pos for particle in self.particles(include_ports=True)), dtype=float)
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
                raise ValueError('Trying to set position of {} with more than one'
                                 'coordinate: {}'.format(self, arrnx3))
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(self._particles(include_ports=False), arrnx3):
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
                raise ValueError('Trying to set position of {} with more than one'
                                 'coordinate: {}'.format(self, arrnx3))
            self.pos = np.squeeze(arrnx3)
        else:
            for atom, coords in zip(self._particles(include_ports=True), arrnx3):
                atom.pos = coords

    @property
    def center(self):
        """The cartesian center of the Compound based on its Particles.

        Returns
        -------
        np.ndarray, shape=(3,), dtype=float
            The cartesian center of the Compound based on its Particles

        """
        if self.xyz.any():
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

    def particles_in_range(self, compound, dmax, max_particles=20, particle_kdtree=None,
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
            particle_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        _, idxs = particle_kdtree.query(compound.pos, k=max_particles, distance_upper_bound=dmax)
        idxs = idxs[idxs != self.n_particles]
        if particle_array is None:
            particle_array = np.array(list(self.particles()))
        return particle_array[idxs]

    def visualize(self, show_ports=False):
        """Visualize the Compound using nglview.

        Allows for visualization of a Compound within a Jupyter Notebook.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Visualize Ports in addition to Particles

        """
        nglview = import_('nglview')
        if run_from_ipython():
            structure = self.to_trajectory(show_ports)
            return nglview.show_mdtraj(structure)
        else:
            raise RuntimeError('Visualization is only supported in Jupyter '
                               'Notebooks.')

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
            load(filename, compound=self, coords_only=True)
            self._update_port_locations(xyz_init)
        else:
            load(filename, compound=self, coords_only=True)

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

    def energy_minimization(self, steps=2500, algorithm='cg',
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

        for particle in self.particles():
            try:
                elem.get_by_symbol(particle.name)
            except KeyError:
                raise MBuildError("Element name {} not recognized. Cannot "
                                  "perform minimization."
                                  "".format(particle.name))

        tmp_dir = tempfile.mkdtemp()
        original = clone(self)
        self._kick()
        self.save(os.path.join(tmp_dir,'un-minimized.mol2'))
        obConversion = openbabel.OBConversion()
        obConversion.SetInAndOutFormats("mol2", "mol2")
        mol = openbabel.OBMol()

        obConversion.ReadFile(mol, os.path.join(tmp_dir, "un-minimized.mol2"))

        ff = openbabel.OBForceField.FindForceField(forcefield)
        if ff is None:
            raise MBuildError("Force field '{}' not supported for energy "
                              "minimization. Valid force fields are 'MMFF94', "
                              "'MMFF94s', 'UFF', 'GAFF', and 'Ghemical'."
                              "".format(forcefield))
        warn("Performing energy minimization using the Open Babel package. Please "
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

        obConversion.WriteFile(mol, os.path.join(tmp_dir, 'minimized.mol2'))
        self.update_coordinates(os.path.join(tmp_dir, 'minimized.mol2'))

    def save(self, filename, show_ports=False, forcefield_name=None,
             forcefield_files=None, box=None, overwrite=False, residues=None,
             references_file=None, combining_rule='lorentz', **kwargs):
        """Save the Compound to a file.

        Parameters
        ----------
        filename : str
            Filesystem path in which to save the trajectory. The extension or
            prefix will be parsed and control the format. Supported
            extensions are: 'hoomdxml', 'gsd', 'gro', 'top', 'lammps', 'lmp'
        show_ports : bool, optional, default=False
            Save ports contained within the compound.
        forcefield_file : str, optional, default=None
            Apply a forcefield to the output file using a forcefield provided
            by the `foyer` package.
        forcefield_name : str, optional, default=None
            Apply a named forcefield to the output file using the `foyer`
            package, e.g. 'oplsaa'. Forcefields listed here:
            https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields
        box : mb.Box, optional, default=self.boundingbox (with buffer)
            Box information to be written to the output file. If 'None', a
            bounding box is used with 0.25nm buffers at each face to avoid
            overlapping atoms.
        overwrite : bool, optional, default=False
            Overwrite if the filename already exists
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.
        references_file : str, optional, default=None
            Specify a filename to write references for the forcefield that is
            to be applied. References are written in BiBTeX format.
        combining_rule : str, optional, default='lorentz'
            Specify the combining rule for nonbonded interactions. Only relevant
            when the `foyer` package is used to apply a forcefield. Valid
            options are 'lorentz' and 'geometric', specifying Lorentz-Berthelot
            and geometric combining rules respectively.

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

        See Also
        --------
        formats.gsdwrite.write_gsd : Write to GSD format
        formats.hoomdxml.write_hoomdxml : Write to Hoomd XML format
        formats.lammpsdata.write_lammpsdata : Write to LAMMPS data format

        """
        extension = os.path.splitext(filename)[-1]
        if extension == '.xyz':
            traj = self.to_trajectory(show_ports=show_ports)
            traj.save(filename)
            return

        # Savers supported by mbuild.formats
        savers = {'.hoomdxml': write_hoomdxml,
                  '.gsd': write_gsd,
                  '.lammps': write_lammpsdata,
                  '.lmp': write_lammpsdata}

        try:
            saver = savers[extension]
        except KeyError:
            saver = None

        if os.path.exists(filename) and not overwrite:
            raise IOError('{0} exists; not overwriting'.format(filename))

        structure = self.to_parmed(box=box, residues=residues)
        # Apply a force field with foyer if specified
        if forcefield_name or forcefield_files:
            from foyer import Forcefield
            ff = Forcefield(forcefield_files=forcefield_files,
                            name=forcefield_name)
            structure = ff.apply(structure, references_file=references_file)
            structure.combining_rule = combining_rule

        total_charge = sum([atom.charge for atom in structure])
        if round(total_charge, 4) != 0.0:
            warn('System is not charge neutral. Total charge is {}.'
                 ''.format(total_charge))

        # Provide a warning if rigid_ids are not sequential from 0
        if self.contains_rigid:
            unique_rigid_ids = sorted(set([p.rigid_id
                                           for p in self.rigid_particles()]))
            if max(unique_rigid_ids) != len(unique_rigid_ids) - 1:
                warn("Unique rigid body IDs are not sequential starting from zero.")

        if saver:  # mBuild supported saver.
            if extension in ['.gsd', '.hoomdxml']:
                kwargs['rigid_bodies'] = [p.rigid_id for p in self.particles()]
            saver(filename=filename, structure=structure, **kwargs)
        else:  # ParmEd supported saver.
            structure.save(filename, overwrite=overwrite, **kwargs)

    def align_vectors(self, align_these, with_these, anchor_pt = None, lattice_override=False):
       """
       Given 2 sets (align_these and with_these) of 2 vectors, rotate a compound
       so that the vectors align_these point in the direction that with_these do.

       :param align_these: list-like
           The vectors to be aligned. Must represent 3D cartesian coordinates.

       :param with_these: list-like
           The vectors to serve as the end goal for align_these to be aligned with.
           Must represent 3D cartesian coordinates.

       :param anchor_pt: optional, accepts list-like, defaults to self.center
           anchor_pt is used as a way to identify a point in the compound that will remain
           unchanged after alignment. The list-like either contains 3D coordinates or the
           hierarchal pathway (see below) to a unique particle that will serve as the anchor
           point.

       :param lattice_override:
           ### revisit when i see justin

       A hierarchal pathway is a list or tuple containing any combination of strings,
       list/tuples, or mb.Compounds. Each element of the list/tuple either describes a series
       of subcompounds (this occurs in the instance where an inner list/tuple is passed), or describes
       one subcompound or type of subcompound (when a string is passed), or even IS a subcompound
       (in the instance where a mb.Compound object is passed). Strings correspond
       to either the names or labels of subcompounds, and list/tuples hold multiple strings that
       correspond to names and/or labels. They are used when the user wishes to describe multiple
       pathways, for example, path = [..., ["subcompound[1]", "subcompound[4]"], ...].
       This example demonstrates that the user can describe some but not all of the pathways
       that have the name "subcompound". The order of the elements in the outer list/tuple correspond
       to their position in the hierarchal pathway, where the first index is the lowest level and the
       last is the highest specified. The number of subcompounds the user can specify is unlimited so
       long as each subcompound specified lies within the hierarchy of the list/tuple element that
       follows. In the context of this function, the first index must correspond to a mb.Particle.

       The idea of a pathway is similar to how one sorts through directories on a computer,
       i.e. "C:/user/username/documents" BUT since MBuild uses hierarchal pathways from lowest
       to highest, the MBuild style of writing it would be "documents/username/user/C:".

       EX: path =["target",
                         ["SubSubSubCompound[0]",
                          "SubSubSubCompound[3]",
                          "SubSubSubCompound[4]"],
                        "SubSubCompound",
                        "SubCompound[6]"]

       """
       if self.made_from_lattice and not lattice_override:
           warn("This compound was made from a lattice, please use the "
                "Lattice.rotate(axis_align= True) or "
                "Lattice.rotate(miller_directions=True) methods."
                " To proceed use the Compound.align_vectors() method  with "
                "this compound, pass align_vectors's optional parameter "
                "lattice_override as True. This compound has not "
                "been changed.")
           return
       for aligner in list([align_these, with_these]):
           if not isinstance(aligner, (list,tuple)):
               if not isinstance(aligner, np.ndarry):
                   raise TypeError("Parameters align_these and with_these must be a list-like of"
                                   " list-like types.")
               else:
                   aligner = aligner.tolist
           else:
               aligner = list(aligner)
           if len(aligner) !=2:
               raise ValueError("Vector pair {} is not of length 2. Both vectors are required to "
                                "sufficienly and concisely describe a plane. If you are in"
                                " the 2D case, please pass (0,0,1) as one of your vectors.".format(aligner))
       ang_current, ang_goal = map(lambda x: angle(x[0], x[1]),
                                   [align_these, with_these])
       if not np.allclose(ang_current, ang_goal, atol= 1e-2):
           raise ValueError("The vectors specified cannot be aligned because the "
                            "angle between the vectors specified in align_these "
                            "is too different from the angle between the vectors "
                            "specified in with_these. Angles were {} and {} degrees, "
                            "respectively.".format(ang_current*180/np.pi,
                                                   ang_goal*180/np.pi))
       align_these, with_these = map(lambda x: normalized_matrix(x), [align_these, with_these])
       if not np.allclose(ang_goal,np.pi/2, atol= 5e-3):
           align_these[1], with_these[1] = map(lambda x: unit_vector(np.cross(x[0], x[1])),
                                                     [align_these, with_these])
           # this ensures that the vector pair will be orthagonal
       if anchor_pt is None:
           anchor_pt = self.center
       else:
           if isinstance(anchor_pt, np.ndarray):
               pass
           elif isinstance(anchor_pt, (tuple, list)):
               if all(isinstance(ap, (int,float)) for ap in anchor_pt):
                   anchor_pt = np.array(anchor_pt)
               else:
                   path = deepcopy(anchor_pt)
                   anchor_pt = list(self.find_particle_in_path(within_path=anchor_pt))
                   if len(anchor_pt) > 1:
                       raise MBuildError("This is not a unique anchor point. "
                                         "The hierarchal path {} is invalid."
                                         "".format(path))
                   else:
                       anchor_pt = anchor_pt[0].pos

           else:
               raise TypeError("Parameter anchor_pt must be of type list, tuple, or"
                               " np.ndarray.")

       self._align(align_these=align_these, with_these=with_these,
                   anchor_pt=anchor_pt, lattice_override=lattice_override)


    def _align(self, align_these, with_these, anchor_pt, lattice_override=False):
       """
       This alignment technique assumes that all the input methods have been checked.
       The function align_vectors() checks input and calls upon this to execute the alignment.
       See def align_vectors() for more information.
       """
       current = deepcopy(np.array(align_these))
       goal = np.array(with_these)
       self.translate(-anchor_pt)
       for ii in range(2):
           if np.allclose(current[ii], goal[ii], atol=1e-3):
               continue
           elif np.allclose(current[ii]*-1, goal[ii], atol= 1e-3):
               self. rotate(theta = np.pi, around= current[(ii+1)%2])
               current[ii]*=-1
               continue
           orthag = np.cross(current[ii], goal[ii])
           theta = abs(angle(current[ii], goal[ii]))
           current = np.array(list(_rotate(coordinates=current, around=orthag, theta=theta)))
           # look into checks here like allclose's
           current = normalized_matrix(current)
           self.rotate(theta=theta, around=orthag)
       self.translate(anchor_pt)



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
    def from_trajectory(self, traj, frame=-1, coords_only=False):
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

        """
        if coords_only:
            if traj.n_atoms != self.n_particles:
                raise ValueError('Number of atoms in {traj} does not match'
                                 ' {self}'.format(**locals()))
            atoms_particles = zip(traj.topology.atoms,
                                  self._particles(include_ports=False))
            for mdtraj_atom, particle in atoms_particles:
                particle.pos = traj.xyz[frame, mdtraj_atom.index]
            return

        atom_mapping = dict()
        for chain in traj.topology.chains:
            if traj.topology.n_chains > 1:
                chain_compound = Compound()
                self.add(chain_compound, 'chain[$]')
            else:
                chain_compound = self
            for res in chain.residues:
                for atom in res.atoms:
                    new_atom = Particle(name=str(atom.name), pos=traj.xyz[frame, atom.index])
                    chain_compound.add(new_atom, label='{0}[$]'.format(atom.name))
                    atom_mapping[atom] = new_atom

        for mdtraj_atom1, mdtraj_atom2 in traj.topology.bonds:
            atom1 = atom_mapping[mdtraj_atom1]
            atom2 = atom_mapping[mdtraj_atom2]
            self.add_bond((atom1, atom2))

        if np.any(traj.unitcell_lengths) and np.any(traj.unitcell_lengths[0]):
            self.periodicity = traj.unitcell_lengths[0]
        else:
            self.periodicity = np.array([0., 0., 0.])

    def to_trajectory(self, show_ports=False, chains=None,
                      residues=None):
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

        Returns
        -------
        trajectory : md.Trajectory

        See also
        --------
        _to_topology

        """
        atom_list = [particle for particle in self.particles(show_ports)]

        top = self._to_topology(atom_list, chains, residues)

        # Coordinates.
        xyz = np.ndarray(shape=(1, top.n_atoms, 3), dtype='float')
        for idx, atom in enumerate(atom_list):
            xyz[0, idx] = atom.pos

        # Unitcell information.
        box = self.boundingbox
        unitcell_lengths = np.empty(3)
        for dim, val in enumerate(self.periodicity):
            if val:
                unitcell_lengths[dim] = val
            else:
                unitcell_lengths[dim] = box.lengths[dim]

        return md.Trajectory(xyz, top, unitcell_lengths=unitcell_lengths,
                             unitcell_angles=np.array([90, 90, 90]))

    def _to_topology(self, atom_list, chains=None, residues=None):
        """Create a mdtraj.Topology from a Compound.

        Parameters
        ----------
        atom_list : list of mb.Compound
            Atoms to include in the topology
        chains : mb.Compound or list of mb.Compound
            Chain types to add to the topology
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.

        Returns
        -------
        top : mdtraj.Topology

        See Also
        --------
        mdtraj.Topology : Details on the mdtraj Topology object

        """
        from mdtraj.core.element import get_by_symbol
        from mdtraj.core.topology import Topology

        if isinstance(chains, string_types):
            chains = [chains]
        if isinstance(chains, (list, set)):
            chains = tuple(chains)

        if isinstance(residues, string_types):
            residues = [residues]
        if isinstance(residues, (list, set)):
            residues = tuple(residues)
        top = Topology()
        atom_mapping = {}

        default_chain = top.add_chain()
        default_residue = top.add_residue('RES', default_chain)

        compound_residue_map = dict()
        atom_residue_map = dict()
        compound_chain_map = dict()
        atom_chain_map = dict()

        for atom in atom_list:
            # Chains
            if chains:
                if atom.name in chains:
                    current_chain = top.add_chain()
                    compound_chain_map[atom] = current_chain
                else:
                    for parent in atom.ancestors():
                        if chains and parent.name in chains:
                            if parent not in compound_chain_map:
                                current_chain = top.add_chain()
                                compound_chain_map[parent] = current_chain
                                current_residue = top.add_residue('RES', current_chain)
                            break
                    else:
                        current_chain = default_chain
            else:
                current_chain = default_chain
            atom_chain_map[atom] = current_chain

            # Residues
            if residues:
                if atom.name in residues:
                    current_residue = top.add_residue(atom.name, current_chain)
                    compound_residue_map[atom] = current_residue
                else:
                    for parent in atom.ancestors():
                        if residues and parent.name in residues:
                            if parent not in compound_residue_map:
                                current_residue = top.add_residue(parent.name, current_chain)
                                compound_residue_map[parent] = current_residue
                            break
                    else:
                        current_residue = default_residue
            else:
                if chains:
                    try: # Grab the default residue from the custom chain.
                        current_residue = next(current_chain.residues)
                    except StopIteration: # Add the residue to the current chain
                        current_residue = top.add_residue('RES', current_chain)
                else: # Grab the default chain's default residue
                    current_residue = default_residue
            atom_residue_map[atom] = current_residue

            # Add the actual atoms
            try:
                elem = get_by_symbol(atom.name)
            except KeyError:
                elem = get_by_symbol("VS")
            at = top.add_atom(atom.name, elem, atom_residue_map[atom])
            at.charge = atom.charge
            atom_mapping[atom] = at

        # Remove empty default residues.
        chains_to_remove = [chain for chain in top.chains if chain.n_atoms == 0]
        residues_to_remove = [res for res in top.residues if res.n_atoms == 0]
        for chain in chains_to_remove:
            top._chains.remove(chain)
        for res in residues_to_remove:
            for chain in top.chains:
                try:
                    chain._residues.remove(res)
                except ValueError:  # Already gone.
                    pass

        for atom1, atom2 in self.bonds():
            # Ensure that both atoms are part of the compound. This becomes an
            # issue if you try to convert a sub-compound to a topology which is
            # bonded to a different subcompound.
            if all(a in atom_mapping.keys() for a in [atom1, atom2]):
                top.add_bond(atom_mapping[atom1], atom_mapping[atom2])
        return top

    def from_parmed(self, structure, coords_only=False):
        """Extract atoms and bonds from a pmd.Structure.

        Will create sub-compounds for every chain if there is more than one
        and sub-sub-compounds for every residue.

        Parameters
        ----------
        structure : pmd.Structure
            The structure to load.
        coords_only : bool
            Set preexisting atoms in compound to coordinates given by structure.

        """
        if coords_only:
            if len(structure.atoms) != self.n_particles:
                raise ValueError('Number of atoms in {structure} does not match'
                                 ' {self}'.format(**locals()))
            atoms_particles = zip(structure.atoms,
                                  self._particles(include_ports=False))
            for parmed_atom, particle in atoms_particles:
                particle.pos = np.array([parmed_atom.xx,
                                         parmed_atom.xy,
                                         parmed_atom.xz]) / 10
            return

        atom_mapping = dict()
        chain_id = None
        chains = defaultdict(list)
        for residue in structure.residues:
            chains[residue.chain].append(residue)

        for chain, residues in chains.items():
            if len(chains) > 1:
                chain_compound = Compound()
                self.add(chain_compound, chain_id)
            else:
                chain_compound = self
            for residue in residues:
                for atom in residue.atoms:
                    pos = np.array([atom.xx, atom.xy, atom.xz]) / 10
                    new_atom = Particle(name=str(atom.name), pos=pos)
                    chain_compound.add(new_atom, label='{0}[$]'.format(atom.name))
                    atom_mapping[atom] = new_atom

        for bond in structure.bonds:
            atom1 = atom_mapping[bond.atom1]
            atom2 = atom_mapping[bond.atom2]
            self.add_bond((atom1, atom2))

        if structure.box is not None:
            self.periodicity = structure.box[0:3]
        else:
            self.periodicity = np.array([0., 0., 0.])

    def to_parmed(self, box=None, title='', residues=None):
        """Create a ParmEd Structure from a Compound.

        Parameters
        ----------
        title : str, optional, default=self.name
            Title/name of the ParmEd Structure
        residues : str of list of str
            Labels of residues in the Compound. Residues are assigned by
            checking against Compound.name.

        Returns
        -------
        parmed.structure.Structure
            ParmEd Structure object converted from self

        See Also
        --------
        parmed.structure.Structure : Details on the ParmEd Structure object

        """
        structure = pmd.Structure()
        structure.title = title if title else self.name
        atom_mapping = {}  # For creating bonds below
        guessed_elements = set()

        if isinstance(residues, string_types):
            residues = [residues]
        if isinstance(residues, (list, set)):
            residues = tuple(residues)

        default_residue = pmd.Residue('RES')
        compound_residue_map = dict()
        atom_residue_map = dict()

        for atom in self.particles():
            if residues and atom.name in residues:
                current_residue = pmd.Residue(atom.name)
                atom_residue_map[atom] = current_residue
                compound_residue_map[atom] = current_residue
            elif residues:
                for parent in atom.ancestors():
                    if residues and parent.name in residues:
                        if parent not in compound_residue_map:
                            current_residue = pmd.Residue(parent.name)
                            compound_residue_map[parent] = current_residue
                        atom_residue_map[atom] = current_residue
                        break
                else:  # Did not find specified residues in ancestors.
                    current_residue = default_residue
                    atom_residue_map[atom] = current_residue
            else:
                current_residue = default_residue
                atom_residue_map[atom] = current_residue

            if current_residue not in structure.residues:
                structure.residues.append(current_residue)

            atomic_number = None
            name = ''.join(char for char in atom.name if not char.isdigit())
            try: atomic_number = AtomicNum[atom.name]
            except KeyError:
                element = element_by_name(atom.name)
                if name not in guessed_elements:
                    warn('Guessing that "{}" is element: "{}"'.format(atom, element))
                    guessed_elements.add(name)
            else:
                element = atom.name

            atomic_number = atomic_number or AtomicNum[element]
            mass = Mass[element]
            pmd_atom = pmd.Atom(atomic_number=atomic_number, name=atom.name,
                                mass=mass, charge=atom.charge)
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10  # Angstroms

            residue = atom_residue_map[atom]
            structure.add_atom(pmd_atom, resname=residue.name,
                               resnum=residue.idx)

            atom_mapping[atom] = pmd_atom

        structure.residues.claim()

        for atom1, atom2 in self.bonds():
            bond = pmd.Bond(atom_mapping[atom1], atom_mapping[atom2])
            structure.bonds.append(bond)
        # pad box with .25nm buffers
        if box is None:
            box = self.boundingbox
            box_vec_max = box.maxs.tolist()
            box_vec_min = box.mins.tolist()
            for dim, val in enumerate(self.periodicity):
                if val:
                    box_vec_max[dim] = val
                    box_vec_min[dim] = 0.0
                if not val:
                    box_vec_max[dim] += 0.25
                    box_vec_min[dim] -= 0.25
            box.mins = np.asarray(box_vec_min)
            box.maxs = np.asarray(box_vec_max)

        box_vector = np.empty(6)
        box_vector[3] = box_vector[4] = box_vector[5] = 90.0
        for dim in range(3):
            box_vector[dim] = box.lengths[dim] * 10
        structure.box = box_vector
        return structure

    def to_intermol(self, molecule_types=None):
        """Create an InterMol system from a Compound.

        Parameters
        ----------
        molecule_types : list or tuple of subclasses of Compound

        Returns
        -------
        intermol_system : intermol.system.System

        """
        from intermol.atom import Atom as InterMolAtom
        from intermol.molecule import Molecule
        from intermol.system import System
        import simtk.unit as u

        if isinstance(molecule_types, list):
            molecule_types = tuple(molecule_types)
        elif molecule_types is None:
            molecule_types = (type(self),)
        intermol_system = System()

        last_molecule_compound = None
        for atom_index, atom in enumerate(self.particles()):
            for parent in atom.ancestors():
                # Don't want inheritance via isinstance().
                if type(parent) in molecule_types:
                    # Check if we have encountered this molecule type before.
                    if parent.name not in intermol_system.molecule_types:
                        self._add_intermol_molecule_type(intermol_system, parent)
                    if parent != last_molecule_compound:
                        last_molecule_compound = parent
                        last_molecule = Molecule(name=parent.name)
                        intermol_system.add_molecule(last_molecule)
                    break
            else:
                # Should never happen if molecule_types only contains type(self)
                raise ValueError('Found an atom {} that is not part of any of '
                                 'the specified molecule types {}'.format(atom, molecule_types))

            # Add the actual intermol atoms.
            intermol_atom = InterMolAtom(atom_index + 1, name=atom.name,
                                         residue_index=1, residue_name='RES')
            intermol_atom.position = atom.pos * u.nanometers
            last_molecule.add_atom(intermol_atom)
        return intermol_system

    @staticmethod
    def _add_intermol_molecule_type(intermol_system, parent):
        """Create a molecule type for the parent and add bonds. """
        from intermol.moleculetype import MoleculeType
        from intermol.forces.bond import Bond as InterMolBond

        molecule_type = MoleculeType(name=parent.name)
        intermol_system.add_molecule_type(molecule_type)

        for index, parent_atom in enumerate(parent.particles()):
            parent_atom.index = index + 1

        for atom1, atom2 in parent.bonds():
            intermol_bond = InterMolBond(atom1.index, atom2.index)
            molecule_type.bonds.add(intermol_bond)

    def __getitem__(self, selection):
        if isinstance(selection, integer_types):
            return list(self.particles())[selection]
        if isinstance(selection, string_types):
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
        newone.periodicity = deepcopy(self.periodicity)
        newone._pos = deepcopy(self._pos)
        newone.port_particle = deepcopy(self.port_particle)
        newone._check_if_contains_rigid_bodies = deepcopy(self._check_if_contains_rigid_bodies)
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
                    newone.labels[label] = compound._clone(clone_of, root_container)
                    compound.referrers.add(clone_of[compound])
                else:
                    # compound is a list of compounds, so we create an empty
                    # list, and add the clones of the original list elements.
                    newone.labels[label] = []
                    for subpart in compound:
                        newone.labels[label].append(subpart._clone(clone_of, root_container))
                        # Referrers must have been handled already, or the will be handled

        return newone

    def _clone_bonds(self, clone_of=None):
        newone = clone_of[self]
        for c1, c2 in self.bonds():
            try:
                newone.add_bond((clone_of[c1], clone_of[c2]))
            except KeyError:
                raise MBuildError("Cloning failed. Compound contains bonds to "
                                  "Particles outside of its containment hierarchy.")


Particle = Compound
