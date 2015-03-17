from __future__ import division

from copy import deepcopy
from distutils.spawn import find_executable
from subprocess import Popen, PIPE
import tempfile

import numpy as np

from mbuild.box import Box
from mbuild.components.small_groups.h2o import H2O
from mbuild.coordinate_transform import translate
from mbuild.periodic_kdtree import PeriodicCKDTree


__all__ = ['solvent_box', 'solvate']


def vdw_radius(atomic_number):
    return 0.2


def solvent_box(solvent, box):
    """Fill a box with solvent.

    Parameters
    ----------
    solvent : Compound
    box : Box

    """

    if np.all(solvent.periodicity == 0):
        solvent.periodicity = solvent.boundingbox.lengths
    num_replicas = np.ceil(box.lengths / solvent.periodicity)   # is this right?
    num_replicas = num_replicas.astype('int')
    print(num_replicas)

    from mbuild.compound import Compound
    compound = Compound()
    translate(solvent, -solvent.boundingbox.mins + box.mins)
    for xi in range(num_replicas[0]):   # x replicas
        for yi in range(num_replicas[1]):   # y replicas
            for zi in range(num_replicas[2]):   # z replicas
                temp_solvent = deepcopy(solvent)   # copy of solvent with box
                # translate solvent box
                translate(temp_solvent, 
                        np.array([xi, yi, zi])*solvent.periodicity)

                # Remove atoms outside the host's box and anything bonded to them.
                guest_atoms = list()   # atoms in box
                guest_atom_pos_list = list()   # positions of atoms in box
                for atom in temp_solvent.yield_atoms():
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


def solvate(solute, solvent, box, overlap=0.1):
    """Solvate a compound in a box of solvent.

    This function will try to use the solvate tool from gromacs 5.0+. If the
    'gmx' executable is not on your path, it will default to a (very) slow
    python implementation. These two approaches will NOT produce identical
    results - use at your own risk.

    Parameters
    ----------
    solute : mb.Compound
    solvent : mb.Compound
    box : mb.Box
    overlap : float

    """

    GMX = find_executable('gmx')
    if not GMX:  # Use the slow python version.
        return _solvate(solute, solvent, box, overlap=overlap)

    # Temporary files used by gmx solvate.
    solute_file = tempfile.mkstemp(suffix='.pdb')[1]
    translate(solute, -solute.boundingbox.mins)
    solute.save(solute_file)
    if solvent == 'water':
        solvent = H2O()
        solvent_file = 'spc216.gro'
    else:
        #solvent_file = tempfile.mkstemp(suffix='.pdb')[1]
        solvent_file = 'solvent.pdb'
        solvent.save(solvent_file)
    solvated_file = tempfile.mkstemp(suffix='.pdb')[1]
    box_lengths = ' '.join("{0:.3f}".format(f) for f in box.lengths)

    # Call gmx solvate.
    gmx_solvate = 'gmx solvate -cp {0} -cs {1} -box {2} -o {3}'.format(
        solute_file, solvent_file, box_lengths, solvated_file)
    proc = Popen(gmx_solvate, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    out, err = proc.communicate()

    # Figure out how many solvent molecules were added...
    n_solvent_line = [line for line in err.splitlines() if 'Number of SOL molecules' in line][0]
    n_solvent = int(n_solvent_line.split()[-1])

    # ...add the appropriate topology information...
    for _ in range(n_solvent):
        solute.add(deepcopy(solvent), 'SOL[$]')

    # ...and the coordinates.
    solute.update_coordinates(solvated_file)
    return solute


def _solvate(solute, solvent, box, overlap=0.1):
    """Slow, fairly naive python implementation. """
    # Replicate the guest so that it fills or overfills the host box.
    host_atom_list = [atom for atom in solute.yield_atoms() if atom.name != 'G']
    host_atom_pos_list = [atom.pos for atom in host_atom_list]
    kdtree = PeriodicCKDTree(host_atom_pos_list)

    if solvent.periodicity.all():
        solvent_box = Box(lengths=solvent.periodicity)
    else:
        solvent_box = solvent.boundingbox
    num_replicas = np.ceil(box.lengths / solvent_box.lengths)
    num_replicas = num_replicas.astype('int')
    for xi in range(num_replicas[0]):
        for yi in range(num_replicas[1]):
            for zi in range(num_replicas[2]):
                guest = deepcopy(solvent)
                translate(guest, -solvent_box.mins + box.mins + np.array([xi, yi, zi]) * solvent_box.lengths)

                # Remove atoms outside the host's box and anything bonded to them.
                guest_atoms = list()
                guest_atom_pos_list = list()
                for atom in guest.yield_atoms():
                    guest_atoms.append(atom)
                    guest_atom_pos_list.append(atom.pos)
                guest_atom_pos_list = np.array(guest_atom_pos_list)

                atoms_to_remove = set()
                atom_indicies = np.where(np.logical_or(np.any(guest_atom_pos_list < box.mins, axis=1),
                        np.any(guest_atom_pos_list > box.maxs, axis=1)))[0]
                for ai in atom_indicies:
                    atoms_to_remove.add(guest_atoms[ai])
                    atoms_to_remove.update(guest_atoms[ai].bonded_atoms())
                guest.remove(atoms_to_remove)

                # Remove overlapping atoms and anything bonded to them.
                atoms_to_remove = set()
                for guest_atom in guest.yield_atoms():
                    _, neighbors = kdtree.query(guest_atom.pos, k=10)
                    for host_atom_idx in neighbors:
                        if host_atom_idx < len(host_atom_list):
                            host_atom = host_atom_list[host_atom_idx]
                            if solute.min_periodic_distance(host_atom, guest_atom) < (2 * overlap):
                                atoms_to_remove.add(guest_atom)
                                atoms_to_remove.update(guest_atom.bonded_atoms())
                guest.remove(atoms_to_remove)
                solute.add(guest, "guest_{}_{}_{}".format(xi, yi, zi))
    return solute
