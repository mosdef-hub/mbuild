from mdtraj.core.topology import Topology as MDTrajTopology
from mdtraj.core.topology import Atom as MDTrajAtom


__all__ = ['Topology']


class Topology(MDTrajTopology):
    """Derivative of MDTraj's Topology class with additional functionality. """

    def __init__(self):
        super(Topology, self).__init__()

        self.forcefield = None
        self._angles = []
        self._dihedrals = []
        self._ff_bonds = set()
        self._ff_angles = set()
        self._ff_dihedrals = set()
        self._ff_impropers = set()
        self._ff_pairs = set()

    @property
    def angles(self):
        """Iterator over all angles (tuple of three atoms). """
        return iter(self._angles)

    @property
    def n_angles(self):
        """Get the number of angles in the Topology"""
        return len(self._angles)

    @property
    def dihedrals(self):
        """Iterator over all dihedrals (tuple of four atoms). """
        return iter(self._dihedrals)

    @property
    def n_dihedrals(self):
        """Get the number of dihedrals in the Topology. """
        return len(self._dihedrals)

    @property
    def ff_bonds(self):
        return iter(self._ff_bonds)

    @property
    def ff_angles(self):
        return iter(self._ff_angles)

    @property
    def ff_dihedrals(self):
        return iter(self._ff_dihedrals)

    @property
    def ff_impropers(self):
        return iter(self._ff_impropers)

    @property
    def ff_pairs(self):
        return iter(self._ff_pairs)

    @property
    def n_ff_bonds(self):
        return len(self._ff_bonds)

    @property
    def n_ff_angles(self):
        return len(self._ff_angles)

    @property
    def n_ff_dihedrals(self):
        return len(self._ff_dihedrals)

    @property
    def n_ff_impropers(self):
        return len(self._ff_impropers)

    def add_ff_bond(self, bond):
        """ """
        self._ff_bonds.add(bond)

    def add_ff_angle(self, angle):
        """ """
        self._ff_angles.add(angle)

    def add_ff_dihedral(self, dihedral):
        """ """
        self._ff_dihedrals.add(dihedral)

    def add_ff_improper(self, improper):
        """ """
        self._ff_impropers.add(improper)

    def gen_pairs(self, n_excl=4):
        if n_excl == 4:
            for dihedral in self.ff_dihedrals:
                self._ff_pairs.add((dihedral.atom1, dihedral.atom4))
        else:
            raise ValueError('Unsupported number of pair exclusions.')


class Atom(MDTrajAtom):
    """Derivative of MDTraj's Atom class with additional functionality.

    Attributes
    ----------
    name : str
        The name of the Atom
    element : mdtraj.element.Element
        The element of the Atoms
    index : int
        The index of the Atom within its Topology
    residue : mdtraj.topology.Residue
        The Residue this Atom belongs to

    """
    def __init__(self, name, element, index, residue, serial=None):
        super(Atom, self).__init__(name, element, index, residue, serial=serial)
        self._neighbors = []

    @property
    def neighbors(self):
        """Return a list of all neighboring Atoms. """
        return self._neighbors

    @neighbors.setter
    def neighbors(self, bonded_atoms):
        self._neighbors = bonded_atoms
