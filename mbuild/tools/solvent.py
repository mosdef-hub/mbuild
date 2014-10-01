from __future__ import division

from copy import deepcopy

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.periodic_kdtree import PeriodicCKDTree
from mbuild.coordinate_transform import translate


def vdw_radius(atomic_number):
    return 0.2


def solvent_box(solvent, box):
    """Fill a box with solvent.

    Args:
        solvent (Compound):
        box (Box):

    """
    if np.all(solvent.periodicity == 0):
        solvent.periodicity = solvent.boundingbox().lengths
    solvent_box = Box(lengths=solvent.periodicity)
    num_replicas = np.ceil(box.lengths / solvent_box.lengths)
    num_replicas = num_replicas.astype('int')
    print num_replicas
    compound = Compound()
    for xi in range(num_replicas[0]):
        for yi in range(num_replicas[1]):
            for zi in range(num_replicas[2]):
                temp_solvent = deepcopy(solvent)
                translate(temp_solvent, -solvent_box.mins + box.mins + np.array([xi, yi, zi])*solvent_box.lengths)

                # Remove atoms outside the host's box and anything bonded to them.
                guest_atoms = list()
                guest_atom_pos_list = list()
                for atom in temp_solvent.atoms():
                    guest_atoms.append(atom)
                    guest_atom_pos_list.append(atom.pos)
                guest_atom_pos_list = np.array(guest_atom_pos_list)

                atoms_to_remove = set()
                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < box.mins, axis=1),
                        np.any(guest_atom_pos_list > box.maxs, axis=1)))[0]
                for ai in atom_indicies:
                    atoms_to_remove.add(guest_atoms[ai])
                    atoms_to_remove.update(guest_atoms[ai].bonded_atoms())
                temp_solvent.remove(atoms_to_remove)
                compound.add(temp_solvent, "solvent_{}_{}_{}".format(xi, yi, zi))
    return compound


def solvate(host_compound, guest_compound, host_box, guest_box, overlap=vdw_radius):
    """Solvate a Compound in a Box.

    Typical usage of this function would be to use a pre-equilibrated solvent
    system as the guest which will then be replicated to fill the host's Box.
    All Compounds that have overlapping Atoms with the host are removed.

    Args:
        host_compound (Compound):
        guest_compound (Compound):
        host_box (Box):
        guest_box (Box):
        overlap (function): A function which translate atoms to overlap radii.

    """
    assert isinstance(host_box, Box)
    assert isinstance(guest_box, Box)

    # TODO: we may want to make sure that the axes of the two boxes line up

    # Replicate the guest so that it fills or overfills the host box.
    host_atom_list = [atom for atom in host_compound.atoms() if atom.kind != 'G']
    host_atom_pos_list = [atom.pos for atom in host_atom_list]
    kdtree = PeriodicCKDTree(host_atom_pos_list)

    num_replicas = np.ceil(host_box.lengths / guest_box.lengths)
    num_replicas = num_replicas.astype('int')
    print num_replicas
    for xi in range(num_replicas[0]):
        for yi in range(num_replicas[1]):
            for zi in range(num_replicas[2]):
                guest = deepcopy(guest_compound)
                translate(guest, -guest_box.mins + host_box.mins + np.array([xi, yi, zi])*guest_box.lengths)

                # Remove atoms outside the host's box and anything bonded to them.
                guest_atoms = list()
                guest_atom_pos_list = list()
                for atom in guest.atoms():
                    guest_atoms.append(atom)
                    guest_atom_pos_list.append(atom.pos)
                guest_atom_pos_list = np.array(guest_atom_pos_list)

                atoms_to_remove = set()
                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < host_box.mins, axis=1),
                        np.any(guest_atom_pos_list > host_box.maxs, axis=1)))[0]
                for ai in atom_indicies:
                    atoms_to_remove.add(guest_atoms[ai])
                    atoms_to_remove.update(guest_atoms[ai].bonded_atoms())
                guest.remove(atoms_to_remove)

                # Remove overlapping atoms and anything bonded to them.
                atoms_to_remove = set()
                for guest_atom in guest.atoms():
                    _, neighbors = kdtree.query(guest_atom.pos, k=10)
                    for host_atom_idx in neighbors:
                        if host_atom_idx < len(host_atom_list):
                            host_atom = host_atom_list[host_atom_idx]
                            if host_compound.min_periodic_distance(host_atom, guest_atom) < (overlap(host_atom) + overlap(guest_atom)):
                                atoms_to_remove.add(guest_atom)
                                atoms_to_remove.update(guest_atom.bonded_atoms())
                guest.remove(atoms_to_remove)
                host_compound.add(guest, "guest_{}_{}_{}".format(xi, yi, zi))
