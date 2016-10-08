from __future__ import print_function, division

__all__ = ['load', 'clone', 'Compound', 'Particle']

import collections
from collections import OrderedDict, defaultdict
from copy import deepcopy
import itertools
import os
import sys
from warnings import warn

import mdtraj as md
import nglview
import numpy as np
from oset import oset as OrderedSet
import parmed as pmd
from parmed.periodic_table import AtomicNum, element_by_name, Mass
from six import integer_types, string_types

from mbuild.bond_graph import BondGraph
from mbuild.box import Box
from mbuild.exceptions import MBuildError
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.utils.io import run_from_ipython
from mbuild.formats.hoomdxml import write_hoomdxml
from mbuild.formats.lammpsdata import write_lammpsdata


def load(filename, relative_to_module=None, compound=None, coords_only=False,
         **kwargs):
    """Load a file into an mbuild compound.

    Parameters
    ----------
    filename : str
    relative_to_module :
    compound : mb.Compound, optional
    coords_only : bool, optional
        Only load the coordinates into an existing compoint.

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

    traj = md.load(filename, **kwargs)
    compound.from_trajectory(traj, frame=-1, coords_only=coords_only)
    return compound


def clone(existing_compound, clone_of=None, root_container=None):
    """A faster alternative to deepcopying.

    Does not resolve circular dependencies. This should be safe provided
    you never try to add the top of a Compound hierarchy to a
    sub-Compound.
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
    Composite design pattern (Gamma, Erich; Richard Helm; Ralph Johnson; John
    M. Vlissides (1995). Design Patterns: Elements of Reusable Object-Oriented
    Software. Addison-Wesley. p. 395. ISBN 0-201-63361-2.), with Compound being
    the composite, and Particle playing the role of the primitive (leaf) part,
    where Particle is in fact simply an alias to the Compound class.

    Compound maintains a list of children (other Compounds contained within), and
    provides a means to tag the children with labels, so that the compounds can
    be easily looked up later. Labels may also point to objects outside the
    Compound's containment hierarchy. Compound has built-in support for copying
    and deepcopying Compound hierarchies, enumerating particles or bonds in the
    hierarchy, proximity based searches, visualization, I/O operations, and a
    number of other convenience methods.

    Parameters
    ----------
    subcompounds : Compound, optional, default=None
        One or more compounds to be added to self.
    name : str, optional, default=self.__class__.__name__
        The type of Compound.
    periodicity : np.ndarray, shape=(3,), dtype=float, optional
        The periodic lengths of the Compound in the x, y and z directions.
        Defaults to zeros which is treated as non-periodic.

    Attributes
    ----------
    name : str, optional, default=self.__class__.__name__
        The type of Compound.
    periodicity : np.ndarray, shape=(3,), dtype=float, optional
        The periodic lengths of the Compound in the x, y and z directions.
        Defaults to zeros which is treated as non-periodic.
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

        # A periodocity of zero in any direction is treated as non-periodic.
        if periodicity is None:
            self._periodicity = np.array([0.0, 0.0, 0.0])
        else:
            self._periodicity = np.asarray(periodicity)

        if pos is not None:
            self._pos = np.asarray(pos, dtype=float)
        else:
            self._pos = np.zeros(3)

        self.charge = charge

        self.parent = None
        self.children = OrderedSet()
        self.labels = OrderedDict()
        self.referrers = set()

        self.bond_graph = None
        self.port_particle = port_particle

        # self.add() must be called after labels and children are initialized.
        if subcompounds:
            self.add(subcompounds)

    def particles(self, include_ports=False):
        """ """
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
        """Yield Compounds below self in the hierarchy. """
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
        """Generate all ancestors of the Compound recursively. """
        if self.parent is not None:
            yield self.parent
            for ancestor in self.parent.ancestors():
                yield ancestor

    @property
    def root(self):
        parent = None
        for parent in self.ancestors():
            pass
        if parent is None:
            return self
        return parent

    def particles_by_name(self, name):
        for particle in self.particles():
            if particle.name == name:
                yield particle

    def add(self, new_child, label=None, containment=True, replace=False,
            inherit_periodicity=True):
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

        """
        # Support batch add via lists, tuples and sets.
        if (isinstance(new_child, collections.Iterable) and
                not isinstance(new_child, string_types)):
            for child in new_child:
                self.add(child)
            return

        if not isinstance(new_child, Compound):
            raise ValueError('Only objects that inherit from mbuild.Compound '
                             'can be added to Compounds. You tried to add '
                             '"{}".'.format(new_child))

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

    def remove(self, objs_to_remove):
        """Remove children from the Compound. """
        if not self.children:
            return

        if not hasattr(objs_to_remove, '__iter__'):
            objs_to_remove = [objs_to_remove]
        objs_to_remove = set(objs_to_remove)

        if len(objs_to_remove) == 0:
            return

        intersection = objs_to_remove.intersection(self.children)
        self.children -= intersection
        objs_to_remove -= intersection

        for removed_part in intersection:
            if self.root.bond_graph and self.root.bond_graph.has_node(removed_part):
                self.root.bond_graph.remove_node(removed_part)
            self._remove_references(removed_part)

        # Remove the part recursively from sub-compounds.
        if self.children:
            for part in self.children:
                part.remove(objs_to_remove)

    @staticmethod
    def _remove_references(removed_part):
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
            for label, part in removed_part.labels.items():
                if removed_part not in part.ancestors():
                    part.referrers.remove(removed_part)
                    labels_to_delete.append(label)
        for label in labels_to_delete:
            del removed_part.labels[label]

    def referenced_ports(self):
        """Return all Ports referenced by this Compound. """
        from mbuild.port import Port
        return [port for port in self.labels.values()
                if isinstance(port, Port)]

    def available_ports(self):
        """Return all unoccupied Ports referenced by this Compound. """
        from mbuild.port import Port
        return [port for port in self.labels.values()
                if isinstance(port, Port) and not port.used]

    def bonds(self):
        """A list of all Bonds in the Compound and sub-Compounds. """
        if self.root.bond_graph:
            if self.root == self:
                return self.root.bond_graph.edges_iter()
            else:
                return self.root.bond_graph.subgraph(self.particles()).edges_iter()
        else:
            return iter(())

    @property
    def n_bonds(self):
        """Return the number of Bonds in the Compound. """
        return sum(1 for _ in self.bonds())

    def add_bond(self, particle_pair):
        """"""
        if self.root.bond_graph is None:
            self.root.bond_graph = BondGraph()

        self.root.bond_graph.add_edge(particle_pair[0], particle_pair[1])

    def generate_bonds(self, name_a, name_b, dmin, dmax):
        """Add Bonds between all pairs of types a/b within [dmin, dmax]. """
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
                if (p2.name == name_b) and (dmin <= self.min_periodic_distance(p2.pos, p1.pos) <= dmax):
                    self.add_bond((p1, p2))
                    added_bonds.append(bond_tuple)

    def remove_bond(self, particle_pair):
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
        pos : np.ndarray, shape=(n, 3)
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
        """Return all particle coordinates in this compound including ports. """
        if not self.children:
            pos = self._pos
        else:
            arr = np.fromiter(itertools.chain.from_iterable(
                particle.pos for particle in self.particles(include_ports=True)), dtype=float)
            pos = arr.reshape((-1, 3))
        return pos

    @property
    def center(self):
        """The cartesian center of the Compound based on its Atoms. """
        if self.xyz.any():
            return np.mean(self.xyz, axis=0)

    @property
    def boundingbox(self):
        """Compute the bounding box of the compound. """
        xyz = self.xyz
        return Box(mins=xyz.min(axis=0), maxs=xyz.max(axis=0))

    def min_periodic_distance(self, xyz0, xyz1):
        """Vectorized distance calculation considering minimum image. """
        d = np.abs(xyz0 - xyz1)
        d = np.where(d > 0.5 * self.periodicity, self.periodicity - d, d)
        return np.sqrt((d ** 2).sum(axis=-1))

    def particles_in_range(self, compound, dmax, max_particles=20, particle_kdtree=None, particle_array=None):
        """Find particles within a specified range of another particle. """
        if particle_kdtree is None:
            particle_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        _, idxs = particle_kdtree.query(compound.pos, k=max_particles, distance_upper_bound=dmax)
        idxs = idxs[idxs != self.n_particles]
        if particle_array is None:
            particle_array = np.array(list(self.particles()))
        return particle_array[idxs]

    def view_hierarchy(self, show_ports=False):
        """Visualize a compound hierarchy as a tree. """
        raise NotImplementedError('Coming soon!')

    def visualize(self, show_ports=False):
        """Visualize the Compound using nglview. """
        if run_from_ipython():
            structure = self.to_trajectory(show_ports)
            return nglview.show_mdtraj(structure)
        else:
            try:
                """Visualize the Compound using imolecule. """
                import imolecule
                json_mol = self._to_json(show_ports)
                imolecule.draw(json_mol, format='json', shader='lambert',
                               drawing_type='ball and stick', camera_type='perspective',
                               element_properties=None)

            except ImportError:
                raise RuntimeError('Visualization is only supported in Jupyter '
                                                      'Notebooks.')

    def _to_json(self, show_ports=False):
        import imolecule
        atoms = list()

        for idx, particle in enumerate(self._particles(include_ports=show_ports)):
            particle.index = idx
            atoms.append({'element': particle.name,'location': list(np.asarray(particle.pos, dtype=float) * 10)})

        bonds = [{'atoms': [atom1.index, atom2.index],'order': 1} for atom1, atom2 in self.bonds()]
        output = {'name': self.name, 'atoms': atoms, 'bonds': bonds}

        # Remove the index attribute on particles.
        for idx, particle in enumerate(self.particles()):
            if not show_ports and particle.port_particle:
                continue
            del particle.index

        return imolecule.json_formatter.compress(output)

    def update_coordinates(self, filename):
        """Update the coordinates of this Compound from a file. """
        load(filename, compound=self, coords_only=True)

    def save(self, filename, show_ports=False, forcefield=None, **kwargs):
        """Save the Compound to a file.

        Parameters
        ----------
        filename : str
            Filesystem path in which to save the trajectory. The extension or
            prefix will be parsed and will control the format.
        show_ports : bool, default=False
            Save ports contained within the compound.
        forcefield : str, default=None
            Apply a forcefield to the output file using the `foyer` package.

        Other Parameters
        ----------------
        force_overwrite : bool

        """
        extension = os.path.splitext(filename)[-1]

        savers = {'.hoomdxml': self.save_hoomdxml,
                  '.gro': self.save_gromacs,
                  '.top': self.save_gromacs,
                  '.lammps': self.save_lammpsdata,
                  '.lmp': self.save_lammpsdata}

        try:
            saver = savers[extension]
        except KeyError:  # TODO: better reporting
            saver = None

        structure = self.to_parmed(show_ports, **kwargs)
        if saver:  # mBuild/InterMol supported saver.
            return saver(filename, structure, forcefield, **kwargs)
        else:  # ParmEd supported saver.
            return structure.save(filename, **kwargs)

    def save_hoomdxml(self, filename, structure, forcefield, box=None, **kwargs):
        """ """
        if forcefield:
            from foyer.forcefield import apply_forcefield
            structure = apply_forcefield(structure, forcefield=forcefield)
        if not box:
            box = self.boundingbox
            for dim, val in enumerate(self.periodicity):
                if val:
                    box.lengths[dim] = val
                    box.maxs[dim] = val
                    box.mins[dim] = 0.0
                if not val:
                    box.maxs[dim] += 0.25
                    box.mins[dim] -= 0.25
                    box.lengths[dim] += 0.5
        write_hoomdxml(structure, filename, forcefield, box, **kwargs)

    def save_gromacs(self, filename, structure, forcefield, force_overwrite=False, **kwargs):
        """ """

        # Create separate file paths for .gro and .top
        filepath, filename = os.path.split(filename)
        basename = os.path.splitext(filename)[0]
        top_filename = os.path.join(filepath, basename + '.top')
        gro_filename = os.path.join(filepath, basename + '.gro')

        if forcefield:
            from foyer.forcefield import apply_forcefield
            structure = apply_forcefield(structure, forcefield=forcefield)
        structure.save(top_filename, 'gromacs', **kwargs)
        structure.save(gro_filename, 'gro', **kwargs)

    def save_lammpsdata(self, filename, structure, forcefield, box=None, **kwargs):
        """ """
        if forcefield:
            from foyer.forcefield import apply_forcefield
            structure = apply_forcefield(structure, forcefield=forcefield)
        if not box:
            box = self.boundingbox
            for dim, val in enumerate(self.periodicity):
                if val:
                    box.lengths[dim] = val
                    box.maxs[dim] = val
                    box.mins[dim] = 0.0
                if not val:
                    box.maxs[dim] += 0.25
                    box.mins[dim] -= 0.25
                    box.lengths[dim] += 0.5

        write_lammpsdata(structure, filename, forcefield, box, **kwargs)

    # Interface to Trajectory for reading/writing .pdb and .mol2 files.
    # -----------------------------------------------------------------
    def from_trajectory(self, traj, frame=-1, coords_only=False):
        """Extract atoms and bonds from a md.Trajectory.

        Will create sub-compounds for every chain if there is more than one
        and sub-sub-compounds for every residue.

        Parameters
        ----------
        traj : md.Trajectory
            The trajectory to load.
        frame : int
            The frame to take coordinates from.
;
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

    def to_trajectory(self, show_ports=False, chain_types=None,
                      residue_types=None, **kwargs):
        """Convert to an md.Trajectory and flatten the compound.

        Parameters
        ----------
        show_ports : bool, optional, default=False
            Include all port atoms when converting to trajectory.

        Returns
        -------
        trajectory : md.Trajectory

        See also
        --------
        _to_topology

        """
        import mdtraj as md

        atom_list = [particle for particle in self.particles(show_ports)]

        top = self._to_topology(atom_list, chain_types, residue_types)

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

    def _to_topology(self, atom_list, chain_types=None, residue_types=None):
        """Create a mdtraj.Topology from a Compound.

        Parameters
        ----------
        atom_list :
        chain_types :
        residue_types :

        Returns
        -------
        top : mtraj.Topology

        """
        from mdtraj.core.element import get_by_symbol
        from mdtraj.core.topology import Topology

        if isinstance(chain_types, Compound):
            chain_types = [Compound]
        if isinstance(chain_types, (list, set)):
            chain_types = tuple(chain_types)

        if isinstance(residue_types, Compound):
            residue_types = [Compound]
        if isinstance(residue_types, (list, set)):
            residue_types = tuple(residue_types)
        top = Topology()
        atom_mapping = {}

        default_chain = top.add_chain()
        default_residue = top.add_residue('RES', default_chain)

        last_residue_compound = None
        last_chain_compound = None
        last_residue = None
        last_chain = None

        for atom in atom_list:
            # Chains
            for parent in atom.ancestors():
                if chain_types and isinstance(parent, chain_types):
                    if parent != last_chain_compound:
                        last_chain_compound = parent
                        last_chain = top.add_chain()
                        last_chain_default_residue = top.add_residue('RES', last_chain)
                        last_chain.compound = last_chain_compound
                    break
            else:
                last_chain = default_chain
                last_chain.compound = last_chain_compound

            # Residues
            for parent in atom.ancestors():
                if residue_types and isinstance(parent, residue_types):
                    if parent != last_residue_compound:
                        last_residue_compound = parent
                        last_residue = top.add_residue(parent.__class__.__name__, last_chain)
                        last_residue.compound = last_residue_compound
                    break
            else:
                if last_chain != default_chain:
                    last_residue = last_chain_default_residue
                else:
                    last_residue = default_residue
                last_residue.compound = last_residue_compound

            # Add the actual atoms
            try:
                elem = get_by_symbol(atom.name)
            except KeyError:
                elem = get_by_symbol("VS")
            at = top.add_atom(atom.name, elem, last_residue)
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
                particle.pos = structure.coordinates[parmed_atom.idx]
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
                    new_atom = Particle(name=str(atom.name), pos=structure.coordinates[atom.idx])
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

    def to_parmed(self, title='', **kwargs):
        """Create a ParmEd Structure from a Compound. """
        structure = pmd.Structure()
        structure.title = title if title else self.name
        atom_mapping = {}  # For creating bonds below
        for atom in self.particles():
            atomic_number = None
            try:
                atomic_number = AtomicNum[atom.name]
            except KeyError:
                element = element_by_name(atom.name)
                warn('Guessing that {} is element: {}'.format(atom, element))
            else:
                element = atom.name

            atomic_number = atomic_number or AtomicNum[element]
            mass = Mass[element]
            pmd_atom = pmd.Atom(atomic_number=atomic_number, name=atom.name,
                                mass=mass)
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10  # Angstroms
            structure.add_atom(pmd_atom, resname='RES', resnum=1)
            atom_mapping[atom] = pmd_atom

        for atom1, atom2 in self.bonds():
            bond = pmd.Bond(atom_mapping[atom1], atom_mapping[atom2])
            structure.bonds.append(bond)

        box = self.boundingbox
        box_vector = np.empty(6)
        box_vector[3] = box_vector[4] = box_vector[5] = 90.0
        for dim, val in enumerate(self.periodicity):
            if val:
                box_vector[dim] = val * 10
            else:
                box_vector[dim] = box.lengths[dim] * 10 + 5
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
                                 'the specified molecule types {}'.format(
                    atom, molecule_types))

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
            descr.append('pos=({: .4f},{: .4f},{: .4f}), '.format(self.pos[0], self.pos[1], self.pos[2]))

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
        newone.charge = deepcopy(self.charge)
        newone.port_particle = deepcopy(self.port_particle)
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
            newone.add_bond((clone_of[c1], clone_of[c2]))

Particle = Compound
