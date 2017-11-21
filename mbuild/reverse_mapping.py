import pdb
import numpy as np
import mbuild as mb
import pdb

from parmed.periodic_table import Mass
from collections import OrderedDict

__all__ = ['reverse_map']

def reverse_map(coarse_grained, mapping_moieties, minimize_energy=True,
        forcefield='UFF', steps=1000):
    """ Reverse map an mb.Compound

    Parameters
    ---------
    coarse_grained : mb.Compound
        original structure
    mapping_moieties : dictionary
        Relate CG bead names to finer-detailed mbuild Compound
    minimize_energy : boolean, optional, default=True
        Perform energy minimization on reverse-mapped compound

    **kwargs : keyword arguments
        Key word arguments passed to energy_minimization

    """
    # Get molecular information through bonding 
    molecules = coarse_grained.bond_graph.connected_components()
    
    aa_system = mb.Compound()
    # CG to AA relates the CG bead to its AA representation
    cg_to_aa = OrderedDict()

    # For each bead, replace it with the appropriate mb compound
    for molecule in molecules:
        new_molecule =  mb.Compound()
        for bead in molecule:
            new_atom = mapping_moieties[bead.name]()
            cg_to_aa[bead] = new_atom
            new_atom.translate(bead.pos)
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)
    # Go back and include bonds
    for p_i, p_j in coarse_grained.bonds():
        mb.force_overlap(cg_to_aa[p_i],
                from_positions=cg_to_aa[p_i].available_ports()[0],
                to_positions=cg_to_aa[p_j].available_ports()[0])
    # Iterative energy minimization
    # Energy minimize each molecule separately,
    # compute RMSD, check tolerance, translate back to CG representation
    if minimize_energy:
        for molecule in aa_system.children:
            molecule.energy_minimization(forcefield=forcefield, steps=steps)
        #aa_system = _energy_minimize_loop(aa_system, cg_to_aa, n_iter=10,)
    return aa_system

def _energy_minimize_loop(aa_system, cg_to_aa, n_iter=10, rel_tol=1e-4,
        forcefield='UFF'):
    """ Minimize reverse-mapped structure according to rmsd

    aa_system : mb.Compound()
    cg_to_aa : OrderedDict()
        Relates CG bead (mb.Compound) to its AA representation (mb.Compound)
    n_iter : int
        Number of iterations for EM loop
    rel_tol : float
        Relative tolerance for RMSD comparisons 
    """
    
    # Loop molecule by molecule for each EM 
    for molecule in aa_system.children:
        loop_counter = 0
        rel_err_i = 100
        old_rmsd = _compute_rmsd(cg_to_aa)
        print("Initial RMSD: {old_rmsd}".format(**locals()))


        # n_iterations or rmsd tolerance
        while loop_counter < n_iter and rel_err_i > rel_tol:
            # Translate AA particles back to CG position
            for cg_particle, aa_particles in cg_to_aa.items():
                aa_particles.translate_to(cg_particle.pos)

            # Minimize energy
            # Reduce iterations here?
            #aa_system.energy_minimization()
            molecule.energy_minimization(forcefield=forcefield)

            # Measure RMSD
            new_rmsd = _compute_rmsd(cg_to_aa)
            print("RMSD ({loop_counter}): {new_rmsd}".format(**locals()))

            # Compute relative error
            rel_err_i = abs(new_rmsd - old_rmsd)/old_rmsd


            # While loop things
            old_rmsd = new_rmsd
            loop_counter+=1

    return aa_system

def _compute_rmsd(cg_to_aa):
    """ Compute RMSD via centers of masses of beads"""
    rmsd = 0
    for cg_particle, aa_particles in cg_to_aa.items():
        cg_com = cg_particle.pos
        aa_com = _compute_center_of_mass(aa_particles.children)
        rmsd += sum([(cg_com[i] - aa_com[i])**2 for i in range(3)])
    return rmsd

def _compute_center_of_mass(particles):
    """ Compute center of mass"""
    masses = [Mass[particle.name[1:2]] for particle in particles]
    total_mass = sum(masses)
    com = np.ndarray(3)
    for i in range(3):
        com[i] = sum([particle.pos[i]*mass/total_mass 
            for particle, mass in zip(particles, masses)])
    return com 


