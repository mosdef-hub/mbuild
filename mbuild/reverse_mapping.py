from collections import OrderedDict

from warnings import warn

from mbuild.compound import Compound
from mbuild.compound import clone
from mbuild.coordinate_transform import force_overlap

import pdb

__all__ = ['reverse_map']

def reverse_map(coarse_grained, mapping_moieties, target_structure=None):
    """ Reverse map an mb.Compound

    Parameters
    ---------
    coarse_grained : mb.Compound
        original structure
    mapping_moieties : dictionary
        Relate CG bead names to finer-detailed mbuild Compound
    target_structure : mb.Compound, optional, default=False
        A target atomistic structure which can be used to reconstruct 
        bonding.
        Bond network in the reverse-mapped structure will be completely
        overridden by the bond network from the target atomistic structure
        Care must be taken that atom indices match perfectly
    
    """
    # Get molecular information through bonding 
    # molecules is a list where each item is a list of bonded components
    molecules = coarse_grained.bond_graph.connected_components()
    
    aa_system = Compound()
    # CG to AA relates the CG bead to its AA representation
    cg_to_aa = OrderedDict()

    # For each bead, replace it with the appropriate mb compound
    # Iterate through each molecule (set of particles that are bonded together)
    for molecule in molecules:
        new_molecule = Compound()
        # Rather than sort through the molecule, which may be unsorted
        # Look at the parent's particles, which will be sorted
        for bead in molecule[0].parent.particles():
            new_atom = clone(mapping_moieties[bead.name])
            cg_to_aa[bead] = new_atom
            new_atom.translate(bead.pos)
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)

    # Go back and include bonds
    if target_structure:
        # If a target atomistic structure is provided, just its bond graph
        # to the reverse-mapped structure
        aa_system.root.bond_graph = None
        target_traj = target_structure.to_trajectory()

        for (i,j) in target_traj.topology.bonds:
            aa_system.add_bond([aa_system[i.index], aa_system[j.index]])
        
    else:
        # If no target atomistic structure is provided, look at each molecule, 
        # working inwards from the ends of the molecule
        
        cg_bonds = list(coarse_grained.bonds())
        # Repeatedly iterate through the coarse grained bonds, but only bond
        # particles that have a certain number of available ports
        while len(cg_bonds) > 0:
            for p_i, p_j in cg_bonds:
                if 0 < len(cg_to_aa[p_i].available_ports()) <= 1 or \
                    0 < len(cg_to_aa[p_j].available_ports()) <= 1:
                            p_i_port, p_j_port = _find_matching_ports(cg_to_aa[p_i], 
                                cg_to_aa[p_j])
                            force_overlap(cg_to_aa[p_i], from_positions=p_i_port, 
                                to_positions=p_j_port)
                            cg_bonds.remove((p_i, p_j))


    # Put molecules back after energy minimization
    for cg_particle, aa_particles in cg_to_aa.items():
        aa_particles.translate_to(cg_particle.pos)

    return aa_system



def _find_matching_ports(i, j):
    """ Find corresponding ports on two mBuild compounds"""

    i_ports = i.available_ports()
    j_ports = j.available_ports()
    i_port_names = [p.name for p in i.available_ports()]
    j_port_names = [p.name for p in j.available_ports()]
    common_name = list(set(i_port_names).intersection(j_port_names))
    if len(common_name) != 1:
        warn("{} ports were found with corresponding names for"
                " particles {} and {}".format(len(common_name), i,j))
    i_port = [p for p in i.available_ports() if p.name == common_name[0]]
    j_port = [p for p in j.available_ports() if p.name == common_name[0]]
    #for j_port in j_ports:
        #jif j_port.name == i_port.name:
            #return i_port, j_port
    return i_port[0], j_port[0]


