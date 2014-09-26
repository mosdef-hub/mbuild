from __future__ import division

from copy import deepcopy

import numpy as np

from box import Box
from bond import Bond
from periodic_kdtree import PeriodicCKDTree
from coordinate_transform import translate


def vdw_radius(atomic_number):
    return 0.2


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
    assert np.all(guest_compound.periodicity == 0)

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
                #guest_atoms = guest.atom_list_by_kind('*')
                #guest_atoms = [atom for atom in guest.atoms()]
                #guest_atom_pos_list = np.array([atom.pos for atom in guest_atoms])
                guest_atoms = list()
                guest_atom_pos_list = list()
                for atom in guest.atoms():
                    guest_atoms.append(atom)
                    guest_atom_pos_list.append(atom.pos)
                guest_atom_pos_list = np.array(guest_atom_pos_list)

                atoms_to_remove = set()
                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < host_box.mins, axis=1),
                        np.any(guest_atom_pos_list > host_box.maxes, axis=1)))[0]
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


# TODO: Move and incorporate into forcefield filtering functions.
def add_bond(compound, type_A, type_B, dmin, dmax, kind=None):
    """Ai-Bj distance is in [dmin, dmax] => add bond A1xB(Ai,Bj) (symmetric)."""

    for a1 in compound.atom_list_by_kind(type_A):
        nearest = compound.atoms_in_range(a1.pos, dmax)
        for a2 in nearest:
            if (a2.kind==type_B) and (dmin <= compound.min_periodic_distance(a2.pos, a1.pos) <= dmax):
                compound.add(Bond(a1, a2, kind=kind))


def add_angle(traj, type_A, type_B, type_C):
    """ """
    for ab in traj.bonds_by_atom_type(type_A, type_B):
        nearest = traj.neighbor_bonds(ab)
        for bc in nearest:
            bc = traj.order_bond(bc, type_B, type_C)
            if bc is None:
                continue

            temp_ang = None
            if ab[1] == bc[0]:
                temp_ang = (ab[0], ab[1], bc[1])
            if ab[0] == bc[1]:
                temp_ang = (ab[0], ab[1], bc[1])
            if temp_ang is None:
                continue

            # to be improved...
            # handle an edge case where we're adding the same angle multiple times
            if type_A == type_B == type_C:
                if ab[0].index < bc[1].index:
                    continue
            traj.top.add_angle(*temp_ang)


def compute_center_of_mass(traj, masses=None):
    """Compute the center of mass for each frame.

    Note: This function alters the equivalent mdtraj by adding the option to
        provide an array of masses.

    Args:
        traj (Trajectory): Trajectory to compute center of mass for.
        masses (np.ndarray, optional):

    Returns:
        com (np.ndarray, shape=(n_frames, 3)): Coordinates of the center of mass
            for each frame.

    """

    com = np.zeros((traj.n_frames, 3))
    if masses is None:
        try:
            masses = np.array([a.element.mass for a in traj.top.atoms])
        except:
            masses = np.ones(traj.n_atoms)
    else:
        assert masses.shape[0] == traj.xyz[0].shape[0]

    masses /= masses.sum()

    for i, x in enumerate(traj.xyz):
        com[i, :] = x.astype('float64').T.dot(masses)
    return com


def compute_inertia_tensor(traj, masses=None):
    """Calculates the moment of inertia tensor for each frame.

    Function adapted from MDAnalysis (https://code.google.com/p/mdanalysis/)

    Args:
        traj (Trajectory): Trajectory to compute center of mass for.
        masses (np.ndarray, optional):

    Returns:
        I (np.ndarray, shape=(n_frames, 3, 3): Moment of inertia tensor for each
            frame.
    """

    if not masses:
        try:
            masses = np.array([a.element.mass for a in traj.top.atoms])
        except:
            masses = np.ones(traj.n_atoms)
    else:
        assert masses.shape[0] == traj.xyz[0].shape[0]

    com = compute_center_of_mass(traj, masses)
    I = np.zeros(shape=(traj.n_frames, 3, 3))
    for n, frame in enumerate(traj.xyz):
        recenteredpos = frame - com[n]
        values = zip(masses, recenteredpos)
        Ixx = reduce(lambda t,a: t+a[0]*(a[1][1]*a[1][1]+a[1][2]*a[1][2]), values, 0.)
        Iyy = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][2]*a[1][2]), values, 0.)
        Izz = reduce(lambda t,a: t+a[0]*(a[1][0]*a[1][0]+a[1][1]*a[1][1]), values, 0.)
        Ixy = Iyx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][1], values, 0.)
        Ixz = Izx = -1*reduce(lambda t,a: t+a[0]*a[1][0]*a[1][2], values, 0.)
        Iyz = Izy = -1*reduce(lambda t,a: t+a[0]*a[1][1]*a[1][2], values, 0.)
        I[n] = np.array([[Ixx, Ixy, Ixz], [Iyx, Iyy, Iyz], [Izx, Izy, Izz]])
    return I

if __name__ == "__main__":
    print "hello"



