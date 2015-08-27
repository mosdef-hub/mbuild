from __future__ import print_function

from collections import OrderedDict, Counter
from copy import deepcopy
import itertools
import os
import sys

import imolecule
import numpy as np
import mdtraj as md
from mdtraj.core.element import get_by_symbol
from mdtraj.core.topology import Topology
from oset import oset as OrderedSet

from mbuild.atom import Atom
from mbuild.box import Box
from mbuild.bond import Bond
from mbuild.formats.mol2 import write_mol2
from mbuild.part import Part
from mbuild.periodic_kdtree import PeriodicCKDTree


__all__ = ['load', 'clone', 'Compound']


def load(filename, relative_to_module=None, frame=-1, compound=None,
         coords_only=False, **kwargs):
    """ """
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
    compound.from_trajectory(traj, frame=frame, coords_only=coords_only)
    return compound


def clone(what, clone_of=None, root_container=None):
    # we can use the good old deepcopy if we wanted to
    # return deepcopy(what)
    return what._clone(clone_of=clone_of, root_container=root_container)

class Compound(Part):
    """A building block in the mBuild hierarchy.

    Compound is the superclass of all composite building blocks in the mBuild
    hierarchy. That is, all composite building blocks must inherit from
    compound, either directly or indirectly. The design of Compound follows the
    Composite design pattern (Gamma, Erich; Richard Helm; Ralph Johnson; John
    M. Vlissides (1995). Design Patterns: Elements of Reusable Object-Oriented
    Software. Addison-Wesley. p. 395. ISBN 0-201-63361-2.), with Compound being
    the composite, and Atom playing the role of the primitive (leaf) part.

    Compound maintains a list of parts (contained Compounds, Atoms, Bonds,
    etc., that inherit from Part), and provides a means to tag the parts
    with labels, so that the parts can be easily looked up later. Labels may
    also point to objects outside the Compound's containment hierarchy.
    Compound has built-in support for copying and deepcopying Compound
    hierarchies, enumerating atoms or bonds in the hierarchy, proximity based
    searches, visualization, I/O operations, and a number of other convenience
    methods.

    Parameters
    ----------
    subcompounds : Parts, optional, default=None
        One or more parts to be added to self.
    kind : str, optional, default=self.__class__.__name__
        The type of Compound.
    periodicity : np.ndarray, shape=(3,), dtype=float, optional
        The periodic lengths of the Compound in the x, y and z directions.
        Defaults to zeros which is treated as non-periodic.

    Attributes
    ----------
    kind : str, optional, default=self.__class__.__name__
        The type of Compound.
    periodicity : np.ndarray, shape=(3,), dtype=float, optional
        The periodic lengths of the Compound in the x, y and z directions.
        Defaults to zeros which is treated as non-periodic.
    parts : OrderedSet
        Contains all child parts. Parts can be Atom, Bond or Compound - anything
        that inherits from Part.
    labels : OrderedDict
        Labels to Compound/Atom mappings. These do not necessarily need not be
        in parts.
    parent : mb.Compound
        The parent Compound that contains this part. Can be None if this
        compound is the root of the containment hierarchy.
    referrers : set
        Other compounds that reference this part with labels.

    """
    def __init__(self, subcompounds=None, kind=None, periodicity=None):
        super(Compound, self).__init__()

        if kind:
            self.kind = kind
        else:
            self.kind = self.__class__.__name__

        # A periodocity of zero in any direction is treated as non-periodic.
        if not periodicity:
            periodicity = np.array([0.0, 0.0, 0.0])
        self._periodicity = periodicity

        self.parts = OrderedSet()
        self.labels = OrderedDict()

        # self.add() must be called after labels and parts are initialized.
        if subcompounds:
            self.add(subcompounds)

    @property
    def atoms(self):
        """A list of all Atoms in the Compound and sub-Compounds.  """
        return self.atom_list_by_name(exclude_ports=True)

    def yield_atoms(self):
        """ """
        return self._yield_parts(Atom)

    @property
    def n_atoms(self):
        """Return the number of Atoms in the Compound. """
        return len(self.atoms)

    def atom_list_by_name(self, name='*', exclude_ports=True):
        """Return a list of Atoms filtered by their name.

        Parameters
        ----------
        name : str
            Return only atoms of this type. '*' indicates all.
        exclude_ports : bool
            Exclude Port particles of kind 'G' - reserved for Ports.

        Returns
        -------
        atom_list : list
            List of Atoms matching the inputs.

        """
        atom_list = []
        for atom in self.yield_atoms():
            if not (exclude_ports and atom.name == 'G'):
                if name == '*':
                    atom_list.append(atom)
                elif atom.name == name:
                    atom_list.append(atom)
        return atom_list

    @property
    def bonds(self):
        """A list of all Bonds in the Compound and sub-Compounds. """
        return self.bond_list_by_kind()

    def yield_bonds(self):
        """ """
        return self._yield_parts(Bond)

    @property
    def n_bonds(self):
        """Return the number of Bonds in the Compound. """
        return len(self.bonds)

    def bond_list_by_kind(self, kind='*'):
        """Return a list of Bonds filtered by their kind. """
        bond_list = []
        for bond in self.yield_bonds():
            if kind == '*':
                bond_list.append(bond)
            elif bond.kind == kind:
                bond_list.append(bond)
        return bond_list

    def _yield_parts(self, part_type):
        """Yield parts of a specified type in the Compound recursively. """
        for part in self.parts:
            # Parts local to the current Compound.
            if isinstance(part, part_type):
                yield part
            # Parts further down the hierarchy.
            if isinstance(part, Compound):
                for subpart in part._yield_parts(part_type):
                    yield subpart

    @property
    def periodicity(self):
        return self._periodicity

    @periodicity.setter
    def periodicity(self, periods):
        self._periodicity = np.array(periods)

    @property
    def xyz(self):
        """Return all atom coordinates in this compound. """
        arr = np.fromiter(itertools.chain.from_iterable(
            atom.pos for atom in self.yield_atoms() if atom.name != 'G'), dtype=float)
        return arr.reshape((-1, 3))

    @property
    def xyz_with_ports(self):
        """Return all atom coordinates in this compound including ports. """
        arr = np.fromiter(itertools.chain.from_iterable(
            atom.pos for atom in self.yield_atoms()), dtype=float)
        return arr.reshape((-1, 3))

    @property
    def center(self):
        """The cartesian center of the Compound based on its Atoms. """
        xyz = self.xyz
        if xyz.any():
            return np.mean(self.xyz, axis=0)
        else:  # It's a port or an empty compound.
            atoms = self.atom_list_by_name('G', exclude_ports=False)
            try:
                return sum(atom.pos for atom in atoms) / len(atoms)
            except ZeroDivisionError as err:
                print('Compound contains no atoms.')
                raise err

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

    def atoms_in_range(self, atom, dmax, max_atoms=20, atom_kdtree=None, atom_array=None):
        """"""
        if atom_kdtree is None:
            atom_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        _, idxs = atom_kdtree.query(atom.pos, k=max_atoms, distance_upper_bound=dmax)
        idxs = idxs[idxs != self.n_atoms]
        if atom_array is None:
            atom_array = np.array(self.atoms)
        return atom_array[idxs]

    def add_bonds(self, type_a, type_b, dmin, dmax, kind=None):
        """Add Bonds between all pairs of types a/b within [dmin, dmax]. """
        atom_kdtree = PeriodicCKDTree(data=self.xyz, bounds=self.periodicity)
        atom_array = np.array(self.atoms)
        added_bonds = list()
        for a1 in self.atom_list_by_name(type_a):
            nearest = self.atoms_in_range(a1, dmax, max_atoms=20,
                                          atom_kdtree=atom_kdtree,
                                          atom_array=atom_array)
            for a2 in nearest:
                bond_tuple = (a1, a2) if id(a1) < id(a2) else (a2, a1)
                if bond_tuple in added_bonds:
                    continue
                if (a2.name == type_b) and (dmin <= self.min_periodic_distance(a2.pos, a1.pos) <= dmax):
                    self.add(Bond(a1, a2, kind=kind))
                    added_bonds.append(bond_tuple)

    def add(self, new_part, label=None, containment=True, replace=False,
            inherit_periodicity=True):
        """Add a part to the Compound.

        Note:
            This does not necessarily add the part to self.parts but may
            instead be used to add a reference to the part to self.labels. See
            'containment' argument.

        Parameters
        ----------
        new_part : mb.Atom, mb.Bond or mb.Compound
            The object to be added to this Compound.
        label : str, optional
            A descriptive string for the part.
        containment : bool, optional, default=True
            Add the part to self.parts.
        replace : bool, optional, default=True
            Replace the label if it already exists.

        """
        # Support batch add via lists, tuples and sets.
        if isinstance(new_part, (list, tuple, set)):
            for part in new_part:
                self.add(part)
            return

        if not isinstance(new_part, Part):
            raise ValueError('Only objects that inherit from mbuild.Part '
                             'can be added to Compounds. You tried to add '
                             '{}'.format(new_part))

        if containment:
            if new_part.parent is not None:
                raise ValueError('Part {} already has a parent: {}'.format(
                    new_part, new_part.parent))
            self.parts.add(new_part)
            new_part.parent = self

        # Add new_part to labels. Does not currently support batch add.
        if label is None:
            label = '_{0}[$]'.format(new_part.__class__.__name__)

        if label is not None:
            if label.endswith('[$]'):
                label = label[:-3]
                if label not in self.labels:
                    self.labels[label] = []
                label_pattern = label + '[{}]'

                count = len(self.labels[label])
                self.labels[label].append(new_part)
                label = label_pattern.format(count)

            if not replace and label in self.labels:
                raise Exception('Label {0} already exists in {1}'.format(label, self))
            else:
                self.labels[label] = new_part
        new_part.referrers.add(self)

        if (inherit_periodicity and isinstance(new_part, Compound) and
                new_part.periodicity.any()):
            self.periodicity = new_part.periodicity

    def remove(self, objs_to_remove):
        """Remove parts (Atom, Bond or Compound) from the Compound. """
        if not isinstance(objs_to_remove, (list, tuple, set)):
            objs_to_remove = [objs_to_remove]
        objs_to_remove = set(objs_to_remove)

        if len(objs_to_remove) == 0:
            return

        intersection = objs_to_remove.intersection(self.parts)
        self.parts -= intersection
        objs_to_remove -= intersection

        for removed_part in intersection:
            self._remove_bonds(removed_part)
            self._remove_references(removed_part)

        # Remove the part recursively from sub-compounds.
        for part in self.parts:
            if isinstance(part, Compound) and len(objs_to_remove) > 0:
                part.remove(objs_to_remove)

    @staticmethod
    def _remove_bonds(removed_part):
        """If removing an atom, make sure to remove the bonds it's part of. """
        if isinstance(removed_part, Atom):
            for bond in removed_part.bonds:
                bond.other_atom(removed_part).bonds.remove(bond)
                if bond.parent is not None:
                    bond.parent.remove(bond)

    @staticmethod
    def _remove_references(removed_part):
        """Remove labels pointing to this part and vice versa. """
        removed_part.parent = None

        # Remove labels in the hierarchy pointing to this part.
        referrers_to_remove = set()
        for referrer in removed_part.referrers:
            if removed_part not in referrer.ancestors():
                for label, referred_part in referrer.labels.items():
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
        return [port for port in self.labels.values() if isinstance(port, Port)]

    def update_coordinates(self, filename):
        """ """
        load(filename, compound=self, coords_only=True)

    def save(self, filename, show_ports=False, forcefield=None, **kwargs):
        """Save the Compound to a file.

        Parameters
        ----------
        filename : str
            Filesystem path in which to save the trajectory. The extension or
            prefix will be parsed and will control the format.

        Other Parameters
        ----------------
        force_overwrite : bool

        """
        extension = os.path.splitext(filename)[-1]

        savers = {'.hoomdxml': self.save_hoomdxml,
                  '.gro': self.save_gromacs,
                  '.top': self.save_gromacs,
                  '.mol2': self.save_mol2,
                  '.lammps': self.save_lammpsdata,
                  '.lmp': self.save_lammpsdata,
                  }

        try:
            saver = savers[extension]
        except KeyError:  # TODO: better reporting
            saver = None

        if (not saver or extension == '.mol2') and forcefield:
            ff_formats = ', '.join(set(savers.keys()) - set(['.mol2']))
            raise ValueError('The only supported formats with forcefield'
                             'information are: {0}'.format(ff_formats))

        if saver:  # mBuild/InterMol supported saver.
            traj = self.to_trajectory(show_ports=show_ports, **kwargs)
            return saver(filename, traj)
        else:  # MDTraj supported saver.
            traj = self.to_trajectory(show_ports=show_ports, **kwargs)
            return traj.save(filename, **kwargs)

    def save_mol2(self, filename, traj, **kwargs):
        """ """
        write_mol2(filename, traj)

    def save_hoomdxml(self, filename, traj, force_overwrite=True, **kwargs):
        """ """
        raise NotImplementedError('Interface to InterMol missing')

    def save_gromacs(self, filename, traj, force_overwrite=True, **kwargs):
        """ """
        raise NotImplementedError('Interface to InterMol missing')

    def save_lammpsdata(self, filename, traj, force_overwrite=True, **kwargs):
        """ """
        raise NotImplementedError('Interface to InterMol missing')

    # Visualization
    # -------------
    def view_hierarchy(self, show_ports=False):
        """Visualize a compound hierarchy as a tree.

        A tree is constructed from the compound hierarchy with self as the root.
        The tree is then rendered in a web browser window using D3.js.

        Note
        ------
        Portions of this code are adapted from https://gist.github.com/mbostock/4339083.
        """
        raise NotImplementedError('To be replaced with igraph')

        # try:
        #     import urlparse
        # except ImportError:
        #     import urllib.parse as urlparse
        # try:
        #     from urllib import pathname2url
        # except ImportError:
        #     from urllib.request import pathname2url
        #
        # try:
        #     import networkx as nx
        # except ImportError:
        #     raise ImportError('Networkx is required to visualize the compound hierarchy.')
        # from networkx.readwrite import json_graph
        # from mbuild.utils.visualization import d3_tree_template
        #
        # tempdir = tempfile.mkdtemp(prefix='mbuild_view_hierarchy_')
        #
        # compound_tree = nx.DiGraph()
        # compound_tree.add_node(self.kind)
        # compound_frequency = Counter([self.kind])
        # for sub_compound in self._yield_parts(Compound):
        #     if not show_ports and sub_compound.kind in ["Port", "subport"]:
        #         continue
        #     compound_frequency[sub_compound.kind] += 1
        #     compound_tree.add_node(sub_compound.kind)
        #     if sub_compound.parent:
        #         compound_tree.add_edge(sub_compound.parent.kind, sub_compound.kind)
        #
        # labels = {"'children'": '"children"'}
        # for compound in compound_tree:
        #     node_key = "'{}'".format(compound)
        #     labels[node_key] = '"{} {:d}"'.format(compound, compound_frequency[compound])
        #
        # json_template = json_graph.tree_data(compound_tree, self.kind,
        #                                      dict(id="name", children="children"))
        # json_template = str(json_template)
        # for label in labels:
        #     json_template = json_template.replace(label, labels[label])
        #
        # # generate png image from the compound
        # self.save_png(os.path.join(tempdir, 'visualize_{}.png'.format(self.kind)), show_ports=show_ports)
        # sub_compounds_dict = {labels["'{}'".format(self.kind)]:"visualize_{}".format(self.kind)}
        #
        # # generate png image from all subcompounds
        # for sub_compound in self._yield_parts(Compound):
        #     if not show_ports and sub_compound.kind in ["Port", "subport"]:
        #         continue
        #     if labels["'{}'".format(sub_compound.kind)] not in sub_compounds_dict.keys():
        #         if not show_ports and sub_compound.kind in ["Port", "subport"]:
        #             continue
        #         sub_compound.save_png(os.path.join(tempdir, 'visualize_{}.png'.format(sub_compound.kind)), show_ports=show_ports)
        #
        #         label = labels["'{}'".format(sub_compound.kind)]
        #         sub_compounds_dict[label] = "visualize_{}".format(sub_compound.kind)
        #
        # for key in sub_compounds_dict:
        #     filename = sub_compounds_dict[key]
        #     json_template = json_template.replace("'name': {}".format(key),
        #                                           '"name": {0}, "icon": "{1}"'
        #                                           .format(key, filename+'.png'))
        #
        # html = d3_tree_template % str(json_template)
        # html_file = os.path.join(tempdir, 'view_hierarchy.html')
        # with open(html_file, 'w') as the_file:
        #     the_file.write(html)
        #
        # webbrowser.open(urlparse.urljoin('file:', pathname2url(os.path.join(tempdir, html_file))))
        # # and leave all temp files in place...

    def visualize(self, show_ports=False, shader='lambert',
                  drawing_type='ball and stick', camera_type='perspective'):
        """Visualize the Compound using imolecule. """
        # try:
        #     __IPYTHON__
        # except NameError:
        #     from imolecule import viewer
        #     viewer.visualize(self._to_json(show_ports=show_ports), title=self.kind)
        # else:
        json_mol = self._to_json(show_ports)
        imolecule.draw(json_mol, format='json', shader=shader,
                       drawing_type=drawing_type, camera_type=camera_type)

    def save_png(self, image_filename, show_ports=False):
        raise NotImplementedError('To be replaced by imolecule')

    def _to_json(self, show_ports=False):
        atoms = list()
        for idx, atom in enumerate(self.atom_list_by_name(exclude_ports=not show_ports)):
            atom.index = idx
            atoms.append({'element': atom.name,
                          'location': list(atom.pos * 10)})

        bonds = [{'atoms': [bond.atom1.index, bond.atom2.index],
                  'order': 1}
                 for bond in self.yield_bonds()]
        output = {'name': self.kind, 'atoms': atoms, 'bonds': bonds}
        return imolecule.json_formatter.compress(output)

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

        """
        if coords_only:
            if traj.n_atoms != self.n_atoms:
                raise ValueError('Number of atoms in {traj} does not match {self}'.format(**locals()))
            for mdtraj_atom, mbuild_atom in zip(traj.topology.atoms, self.atoms):
                mbuild_atom.pos = traj.xyz[frame, mdtraj_atom.index]
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
                    new_atom = Atom(str(atom.name), traj.xyz[frame, atom.index])
                    chain_compound.add(new_atom, label='{0}[$]'.format(atom.name))
                    atom_mapping[atom] = new_atom

        for a1, a2 in traj.topology.bonds:
            atom1 = atom_mapping[a1]
            atom2 = atom_mapping[a2]
            self.add(Bond(atom1, atom2))

        if np.any(traj.unitcell_lengths) and np.any(traj.unitcell_lengths[0]):
            self.periodicity = traj.unitcell_lengths[0]
        else:
            self.periodicity = np.array([0., 0., 0.])

    def to_trajectory(self, show_ports=False, chain_types=None, residue_types=None):
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
        mbuild.topology

        """
        exclude = not show_ports
        atom_list = self.atom_list_by_name('*', exclude_ports=exclude)
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

        for bond in self.bonds:
            a1 = bond.atom1
            a2 = bond.atom2
            # Ensure that both atoms are part of the compound. This becomes an
            # issue if you try to convert a sub-compound to a topology which is
            # bonded to a different subcompound.
            if all(a in atom_mapping.keys() for a in [a1, a2]):
                top.add_bond(atom_mapping[a1], atom_mapping[a2])
        return top

    # Interface to InterMol for writing fully parameterized systems.
    # --------------------------------------------------------------
    def _to_intermol(self, molecule_types=None):
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

        if isinstance(molecule_types, list):
            molecule_types = tuple(molecule_types)
        elif molecule_types is None:
            molecule_types = (type(self),)
        intermol_system = System()

        last_molecule_compound = None
        for atom_index, atom in enumerate(self.atoms):
            for parent in atom.ancestors():
                # Don't want inheritance via isinstance().
                if type(parent) in molecule_types:
                    # Check if we have encountered this molecule type before.
                    if parent.kind not in intermol_system.molecule_types:
                        self._add_intermol_molecule_type(intermol_system, parent)
                    if parent != last_molecule_compound:
                        last_molecule_compound = parent
                        last_molecule = Molecule(name=parent.kind)
                        intermol_system.add_molecule(last_molecule)
                    break
            else:
                # Should never happen if molecule_types only contains type(self)
                raise ValueError('Found an atom {} that is not part of any of '
                                 'the specified molecule types {}'.format(atom, molecule_types))

            # Add the actual intermol atoms.
            intermol_atom = InterMolAtom(atom_index + 1, name=atom.name,
                                         residue_index=1, residue_name='RES')
            intermol_atom.position = atom.pos
            last_molecule.add_atom(intermol_atom)
        return intermol_system

    @staticmethod
    def _add_intermol_molecule_type(intermol_system, parent):
        """Create a molecule type for the parent and add bonds. """
        from intermol.moleculetype import MoleculeType
        from intermol.forces.bond import Bond as InterMolBond

        molecule_type = MoleculeType(name=parent.kind)
        intermol_system.add_molecule_type(molecule_type)

        for index, parent_atom in enumerate(parent.atoms):
            parent_atom.index = index + 1

        for bond in parent.bonds:
            intermol_bond = InterMolBond(bond.atom1.index, bond.atom2.index)
            molecule_type.bonds.add(intermol_bond)

    # Magic
    # -----
    def __getattr__(self, attr):
        assert 'labels' != attr, ('Compound __init__ never called. Make '
                                  'sure to call super().__init__() in the '
                                  '__init__ method of your class.')
        if attr in self.labels:
            return self.labels[attr]
        else:
            raise AttributeError("'{}' object has no attribute '{}'".format(
                self, attr))

    def __repr__(self):
        descr = ['<{:s}, {:d} atoms, {:d} bonds, '.format(
            self.kind, self.n_atoms, self.n_bonds
        )]
        if any(self.periodicity):
            descr.append('periodicity: {}'.format(self.periodicity))
        else:
            descr.append('non-periodic')
        descr.append('; ID: {}>'.format(id(self)))
        return ''.join(descr)


    def _clone(self, clone_of=None, root_container=None):
        from copy import deepcopy
        # make the current container the root container if it's not yet specified
        if not root_container:
            root_container=self

        # create the clone_of dict if it's None
        if not clone_of:
            clone_of=dict()

        # if this compound has been cloned, return it
        if self in clone_of:
            return clone_of[self]

        # else we make a new clone

        cls = self.__class__
        newone = cls.__new__(cls)

        # remember that we're cloning the new one of of self
        clone_of[self] = newone

        # First copy those attributes that don't need deepcopying.

        newone.kind = self.kind

        newone.periodicity = deepcopy(self.periodicity)

        # Create empty containers.
        newone.parts = OrderedSet()
        newone.labels = OrderedDict()
        newone.referrers = set()

        # parent should be None initially
        newone.parent = None

        # Add parts to clone
        for part in self.parts:
            if isinstance(part, Bond) and part.has_atoms_outside_of(root_container):
                # ignore bonds with atoms outside the hierarchy.
                continue
            else:
                newpart = clone(part, clone_of, root_container)
                newone.parts.add(newpart)
                newpart.parent = newone

        # Copy labels, except bonds with atoms outside the hierarchy
        for label, part in self.labels.items():
            if isinstance(part, Bond) and part.has_atoms_outside_of(root_container):
                # it's a bond that has atoms outside the current containment hierarchy, so we skip it
                continue
            else:
                if not isinstance(part, list):
                    newone.labels[label] = clone(part, clone_of, root_container)
                    part.referrers.add(clone_of[part])
                else:
                    # part is a list of parts, so we create an empty list, and add the clones of the original list elements
                    newone.labels[label] = []
                    for p in part:
                        newone.labels[label].append(clone(p, clone_of, root_container))
                        # referrers must have been handled already, or the will be handled
        return newone


    def __deepcopy__(self, memo):

        cls = self.__class__
        newone = cls.__new__(cls)
        if len(memo) == 0:
            memo[0] = self
        memo[id(self)] = newone

        # First copy those attributes that don't need deepcopying.
        newone.kind = deepcopy(self.kind, memo)
        newone.periodicity = deepcopy(self.periodicity, memo)

        # Create empty containers.
        newone.parts = OrderedSet()
        newone.labels = OrderedDict()
        newone.referrers = set()

        # Copy the parent of everyone, except topmost Compound being deepcopied.
        if memo[0] == self:
            newone.parent = None
        else:
            newone.parent = deepcopy(self.parent, memo)

        # Copy parts, except bonds with atoms outside the hierarchy.
        for part in self.parts:
            if isinstance(part, Bond):
                if memo[0] in part.atom1.ancestors() and memo[0] in part.atom2.ancestors():
                    newone.parts.add(deepcopy(part, memo))
            else:
                newone.parts.add(deepcopy(part, memo))

        # Copy labels, except bonds with atoms outside the hierarchy
        for k, v in self.labels.items():
            if isinstance(v, Bond):
                if memo[0] in v.atom1.ancestors() and memo[0] in v.atom2.ancestors():
                    newone.labels[k] = deepcopy(v, memo)
                    newone.labels[k].referrers.add(newone)
            else:
                newone.labels[k] = deepcopy(v, memo)
                if not isinstance(newone.labels[k], list):
                    newone.labels[k].referrers.add(newone)

        # Copy referrers that do not point out of the hierarchy.
        for r in self.referrers:
            if memo[0] in r.ancestors():
                newone.referrers.add(deepcopy(r, memo))

        return newone

