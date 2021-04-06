import itertools as it

import numpy as np

from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.coordinate_transform import force_overlap
from mbuild.lib.atoms import H
from mbuild.utils.validation import assert_port_exists
from mbuild import clone


__all__ = ["Polymer"]


class Polymer(Compound):
    """Connect one or more components in a specified sequence.

    Attributes
    ----------
    monomers : List of mbuild.Compounds
        The compound(s) to replicate. Add to this list using the add_monomers
        method.
    end_groups : List of mbuild.Compounds
        The compound to cap the end of the polymer. Add to this list using the
        add_end_groups method.

    Methods
    -------
    add_monomer(monomer, indices, separation, port_labels, orientation, replace)
        Use to add a monomer compound to Polymer.monomers

    add_end_groups(compound, index, separation, orientation, replace)
        Use to add an end group compound to Polymer.end_groups

    build(n, sequence)
        Use to create a single polymer compound. This method uses the compounds
        created by calling the add_monomer and add_end_group methods.
    """

    def __init__(self, monomers=None, end_groups=None):
        super(Polymer, self).__init__()
        self._monomers = monomers or []
        self._end_groups = end_groups or [None, None]
        if len(self._end_groups) != 2:
            raise ValueError(
                "Please provide two end groups;"
                f" you provided {len(self._end_groups)}"
            )
        self._port_labels = ["up", "down"]
        self._headtail = [None, None]

    @property
    def monomers(self):
        return self._monomers

    @property
    def end_groups(self):
        return self._end_groups

    def build(self, n, sequence="A", add_hydrogens=True):
        """Connect one or more components in a specified sequence.

        Parameters
        ----------
        monomers : mbuild.Compound or list of mbuild.Compound
            The compound(s) to replicate.
        n : int
            The number of times to replicate the sequence.
        sequence : str, optional, default 'A'
            A string of characters where each unique character represents one
            repetition of a monomer. Characters in `sequence` are assigned to
            monomers in the order they appear in `Polymer.monomers`.
        add_hydrogens : bool, default True
            If True and an end group compound is None, then the head or tail
            of the polymer will be capped off with hydrogen atoms. If end group
            compounds exist, then they will be used.
            If False and an end group compound is None, then the head or tail
            port will be exposed in the polymer.
        """
        if n < 1:
            raise ValueError("n must be 1 or more")
        n_monomers = n * len(sequence)

        for monomer in self._monomers:
            for label in self._port_labels:
                assert_port_exists(label, monomer)

        unique_seq_ids = sorted(set(sequence))
        if len(self._monomers) != len(unique_seq_ids):
            raise ValueError(
                "Number of monomers passed to `Polymer` class must"
                " match number of unique entries in the specified"
                " sequence."
            )

        # 'A': monomer_1, 'B': monomer_2....
        seq_map = dict(zip(unique_seq_ids, self._monomers))
        last_part = None
        for n_added, seq_item in enumerate(it.cycle(sequence)):
            this_part = clone(seq_map[seq_item])
            self.add(this_part, "monomer[$]")
            if last_part is None:
                first_part = this_part
            else:
                # Transform this part, such that its bottom port is rotated
                # and translated to the last parts top port.
                force_overlap(
                    this_part,
                    this_part.labels[self._port_labels[0]],
                    last_part.labels[self._port_labels[1]],
                )
            last_part = this_part
            if n_added == n * len(sequence) - 1:
                break

        # Add the end groups
        head = self["monomer[0]"]  # First monomer
        tail = self["monomer[{}]".format(n_monomers - 1)]  # Last monomer
        if not head["up"].used:
            head_port = head["up"]
        else:
            head_port = None
        if not tail["down"].used:
            tail_port = tail["down"]
        else:
            tail_port = None

        head_tail = [head_port, tail_port]

        for i, compound in enumerate(self._end_groups):
            if compound is not None:
                if self._headtail[i] is not None:
                    head_tail[i].update_separation(self._headtail[i])
                self.add(compound)
                force_overlap(compound, compound.labels["up"], head_tail[i])
            else:
                if add_hydrogens:
                    hydrogen = H()
                    # Defaut to 1/2 H-C bond len
                    head_tail[i].update_separation(0.0547)
                    hydrogen['up'].update_separation(0.0547)
                    self.add(hydrogen)
                    force_overlap(hydrogen, hydrogen["up"], head_tail[i])
                else:
                    # if None, hoist port to polymer level
                    self.add(
                            head_tail[i],
                            self._port_labels[i],
                            containment=False
                            )

        for port in self.all_ports():
            if port not in self.available_ports():
                self.remove(port)

    def add_monomer(
        self,
        compound,
        indices,
        separation=None,
        orientation=[None, None],
        replace=True,
    ):
        """Add an mBuild compound to self.monomers which will be used to build
        the polymer. Call this function for each unique monomer to be used in
        the polymer.

        Parameters
        ----------
        compound : mbuild.Compound
            A compound of the individual monomer
        indices : list of int of length 2
            The particle indicies of compound that represent the polymer
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
        port_labels = ["up", "down"]
        comp = clone(compound)

        for idx, label, orientation in zip(indices, port_labels, orientation):
            _add_port(comp, label, idx, separation, orientation, replace)
        if replace:
            remove_atom1 = comp[indices[0]]
            remove_atom2 = comp[indices[1]]
            comp.remove(remove_atom1)
            comp.remove(remove_atom2)
        self._monomers.append(comp)

    def add_end_groups(
        self,
        compound,
        index,
        separation=None,
        orientation=None,
        replace=True,
        label="head",
        duplicate=True,
    ):
        """
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
        replace : Bool, default True
            If True, then the particle identified by `index` will be removed
            and ports are added to the particle it was initially bonded to.
            Only use `replace=True` in the case that index points to a hydrogen
            atom bonded to the desired bonding site particle.
            If False, then the particle identified by `index` will have a port
            added and no particle is removed from the end group compound.
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
        separation = _add_port(
            comp, "up", index, separation, orientation, replace
        )
        if replace:
            comp.remove(comp[index])
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


def _add_port(compound, label, idx, separation, orientation=None, replace=True):
    """Add the ports to the compound at the specified particle index, either
    using that particle as the anchor or replacing it entirely.
    """
    if replace:
        atom_bonds = [b for b in compound.bonds() if compound[idx] in b][0]
        anchor = [p for p in atom_bonds if p != compound[idx]][0]
        if orientation is None:
            orientation = compound[idx].pos - anchor.pos
        if separation is None:
            separation = np.linalg.norm(compound[idx].pos - anchor.pos)
    else:
        anchor = compound[idx]

    port = Port(
        anchor=anchor, orientation=orientation, separation=separation / 2,
    )
    compound.add(port, label=label)
    return separation
