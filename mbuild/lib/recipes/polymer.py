import itertools as it

from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.coordinate_transform import force_overlap
from mbuild.utils.validation import assert_port_exists
from mbuild.lib.atoms import H
from mbuild import clone


__all__ = ['Polymer']


class Polymer(Compound):
    """Connect one or more components in a specified sequence.

    Parameters
    ----------
    monomers : mb.Compound or list of mb.Compound
        The compound(s) to replicate.
    n : int
        The number of times to replicate the sequence.
    sequence : str, optional, default='A'
        A string of characters where each unique character represents one
        repetition of a monomer. Characters in `sequence` are assigned to
        monomers in the order assigned by the built-in `sorted()`.
    port_labels : 2-tuple of strs, optional, default=('up', 'down')
        The names of the two ports to use to connect copies of proto.

    """
    def __init__(self):
        super(Polymer, self).__init__()
        self.monomers = []
        self.port_labels = []
        self.end_groups = []

    def build(self, n, sequence='A'):
        if n < 1:
            raise ValueError('n must be 1 or more')
        if isinstance(self.monomers, Compound):
            self.monomers = (self.monomers,)
        for monomer in self.monomers:
            for label in self.port_labels:
                assert_port_exists(label, monomer)

        unique_seq_ids = sorted(set(sequence))
        if len(self.monomers) != len(unique_seq_ids):
            raise ValueError('Number of monomers passed to `Polymer` class must'
                             ' match number of unique entries in the specified'
                             ' sequence.')

        # 'A': monomer_1, 'B': monomer_2....
        seq_map = dict(zip(unique_seq_ids, self.monomers))

        last_part = None
        for n_added, seq_item in enumerate(it.cycle(sequence)):
            this_part = clone(seq_map[seq_item])
            self.add(this_part, 'monomer[$]')
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that it's bottom port is rotated
                # and translated to the last part's top port.
                force_overlap(this_part,
                              this_part.labels[self.port_labels[1]],
                              last_part.labels[self.port_labels[0]])
            last_part = this_part
            if n_added == n * len(sequence) - 1:
                break

        # Hoist the last part's top port to be the top port of the polymer.
        self.add(last_part.labels[self.port_labels[0]], self.port_labels[0], containment=False)

        # Hoist the first part's bottom port to be the bottom port of the polymer.
        self.add(first_part.labels[self.port_labels[1]], self.port_labels[1], containment=False)

        # Add the end groups
        head = self['monomer[0]'] # First monomer group
        tail = self['monomer[{}]'.format(n - 1)] # Last monomer group
        for label in self.port_labels:
            if not head[label].used:
                head_port = head[label]
            if not tail[label].used:
                tail_port = tail[label]

        if not self.end_groups:
            self.end_groups.extend([H(), H()])
        for compound in self.end_groups:
            self.add(compound)
        
        # Update head_port and tail_port orientation and separation before using force_overlap
        # Goal is for the orientation and separation in head and tail ports to match
        # what was given in the ports created in the self.end_groups compounds

        head_port.orientation = -self.end_groups[0].labels['up'].orientation
        head_port.separation = self.end_groups[0].labels['up'].separation

        tail_port.orientation = -self.end_groups[1].labels['up'].orientation
        tail_port.separation = self.end_groups[1].labels['up'].separation

        force_overlap(self.end_groups[0],
                     self.end_groups[0].labels['up'],
                     head_port
                     )
        force_overlap(self.end_groups[1],
                     self.end_groups[1].labels['up'],
                     tail_port
                     )

        
    def add_monomer(self, monomer, bonding_indices, separation,
                    port_labels=['A', 'B'], orientation=None,
                    replace=True):
        """
        Add an mBuild compound to self.monomers which will be used to build the polymer.
        Call this function for each unique monomer to be used in the polymer.

        
        Parameters
        ----------
        monomer : mb.Compound
            A compound of an individual monomer
        bonding_indices : list of int of length 2
            The particle indicies of monomer that represent the polymer
            bonding sites. You can specify the indices of particles that will
            be replaced by the polymer bond, or indices of particles that act
            as the bonding sites. See the 'replace' parameter notes.
        separation : float, units nm
            The bond length desired at the monomer-monomer bonding site.
            (separation / 2) is used to set the length of each port
        port_labels : list of str of length 2, default=['A', 'B']
            Labels given to the two ports added to monomer.
            Ex.) ['head', 'tail'] or ['A', 'B']
            The same port labels must be used for any subsequent
            monomer created using add_monomer()
        orientation : array-like, shape=(3,), default=None
            Vector along which to orient the port
            If replace = True, then the orientation of the bond
            between the particle being removed and the anchor particle
            is used.
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
        if self.port_labels:
            if not sorted(set(port_labels)) == sorted(set(self.port_labels)):
                raise ValueError("The port labels given for each" +
                                "monomer must match. The previous" +
                                "port labels used were {}".format(self.port_labels)
                                )
        else:
            self.port_labels.extend(port_labels)

        for idx, label in zip(bonding_indices, port_labels):
            _add_port(monomer, label, idx, separation, orientation, replace)

        self.monomers.append(monomer)

    def add_end_groups(self, compound, bond_index, separation, orientation=None, replace=True):
        """
        """
        compound_2 = clone(compound)
        _add_port(compound, 'up', bond_index, separation, orientation, replace)
        _add_port(compound_2, 'up', bond_index, separation, orientation, replace)
        self.end_groups.extend([compound, compound_2])


def _add_port(compound, label, atom_idx, separation, orientation=None, replace=True):
    """
    """
    if replace:
        atom_bonds = [bond for bond in compound.bonds() if compound[atom_idx] in bond][0]
        anchor_particle = [p for p in atom_bonds if p != compound[atom_idx]][0]
        orientation = compound[atom_idx].pos - anchor_particle.pos
        compound.remove(compound[atom_idx])
    else:
        anchor_particle = compound[atom_idx]
        
    port = Port(anchor = anchor_particle,
                orientation=orientation,
                separation=separation/2
                )
    compound.add(port, label=label)

