from collections import OrderedDict

from mbuild.compound import Compound
from mbuild.compound import clone
from mbuild.coordinate_transform import force_overlap
from mbuild.exceptions import MBuildError
from warnings import warn

import pdb

__all__ = ['coarse_grain', 'reverse_map']


def coarse_grain(real_thing, memo=None, particle_classes=None):
    if memo is None:
        memo = OrderedDict()

    if particle_classes is None:
        particle_classes = []

    proxy = _create_proxy_compounds(real_thing, memo, particle_classes)
    _create_proxy_bonds(real_thing, memo, particle_classes)
    _create_proxy_labels(real_thing, memo)

    return proxy


def reverse_map(coarse_grained, mapping_moieties, target_structure=None,
        energy_minimize=True, n_loops = 10, **kwargs):
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
    minimize_energy : boolean, optional, default=True
        Perform energy minimization on reverse-mapped compound
    n_loops : int, optional, default=True
        Number of energy minimization loops to perform

    **kwargs : keyword arguments
        Key word arguments passed to energy_minimization

    """
    # Get molecular information through bonding 
    molecules = coarse_grained.bond_graph.connected_components()
    
    aa_system = Compound()
    # CG to AA relates the CG bead to its AA representation
    cg_to_aa = OrderedDict()

    # For each bead, replace it with the appropriate mb compound
    for molecule in molecules:
        new_molecule =  Compound()
        for bead in molecule:
            new_atom = clone(mapping_moieties[bead.name])
            cg_to_aa[bead] = new_atom
            new_atom.translate(bead.pos)
            new_molecule.add(new_atom)
        aa_system.add(new_molecule)

    # Go back and include bonds
    if target_structure:
        aa_system.root.bond_graph = None
        target_traj = target_structure.to_trajectory()

        for (i,j) in target_traj.topology.bonds:
            aa_system.add_bond([aa_system[i.index], aa_system[j.index]])
        
    else:
        # What if i look at bonding particle-by-particle?
        for p_i, p_j in coarse_grained.bonds():
            # So p_i is supposed to bond to p_i_partners
            p_i_bonds = [bond for bond in coarse_grained.bonds() if 
                    p_i in bond]
            # Sort through each of the bonding partners to find the correct matches
            bonding_candidates = _find_all_matching_ports(cg_to_aa,
                    p_i, p_i_bonds)
            # Now that we know all the bonding port candidates, 
            # Find the ports with the fewest common names
            re_arranged = sorted(bonding_candidates, key=lambda x: len(x[2]))

            # Take the first bonding pair (bonding partners with the fewest
            # common names, and bond them
            i_port = [port for port in cg_to_aa[bonding_candidates[0][0]].available_ports() if port.name==bonding_candidates[0][2]]
            j_port = [port for port in cg_to_aa[bonding_candidates[0][1]].available_ports() if port.name==bonding_candidates[0][2]]
            if len(i_port) ==0 or len(j_port ==0):
                print("Bonding issue with {}".format(p_i))
            force_overlap(cg_to_aa[p_i], from_positions=i_port[0], to_positions=j_port[0])





        # This is the stuff that was working
        #for p_i, p_j in coarse_grained.bonds():
        #    p_i_port, p_j_port = _find_matching_ports(cg_to_aa[p_i], 
        #            cg_to_aa[p_j])
        #    print("{} and {}: {} and {}".format(p_i.name, p_j.name,
        #        p_i_port.name, p_j_port.name))
        #    force_overlap(cg_to_aa[p_i], from_positions=p_i_port, 
        #            to_positions=p_j_port)

    # Put molecules back after energy minimization
    for cg_particle, aa_particles in cg_to_aa.items():
        aa_particles.translate_to(cg_particle.pos)

    # Iterative energy minimization
    # Energy minimize each molecule separately,
    if energy_minimize:
        for i in range(n_loops):
           # Put molecules back after energy minimization
           for cg_particle, aa_particles in cg_to_aa.items():
                aa_particles.translate_to(cg_particle.pos) 
           for molecule in aa_system.children:
                molecule.energy_minimize(**kwargs)
            
    return aa_system

def _find_all_matching_ports(cg_to_aa, p_i, p_i_bonds):
    """ Determine the ports that should be bonded to each other"""

    bonding_candidates =[]
    for bond in p_i_bonds:
        if p_i != bond[0]:
            p_j = bond[0]
        else:
            p_j = bond[1]
        #i_ports = p_i.available_ports()
        #j_ports = p_j.available_ports()
        i_port_names = [p.name for p in cg_to_aa[p_i].available_ports()]
        j_port_names = [p.name for p in cg_to_aa[p_j].available_ports()]
        common_names = list(set(i_port_names).intersection(j_port_names))
        bonding_candidates.append([p_i, p_j, common_names])
    return bonding_candidates
        



def _find_matching_ports(i, j):
    """ Find corresponding ports on two mBuild compounds"""
    def _sort_by_name(port):
        return port.name

    #i_port = i_ports[0]
    #j_ports = sorted(j.available_ports(), key=_sort_by_name)

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




class Proxy(Compound):

    def __init__(self, compound):
        if compound.name == 'G':
            name = 'G'
        else:
            name = compound.name + ' (proxy) '
        super(Proxy, self).__init__(name=name)

        self.wrapped = compound
        self.children = None
        self.labels = None
        self.parent = None
        self.referrers = set()
        self.index = None
        self.bond_graph = None

    def proxy_for(self):
        if hasattr(self.wrapped, 'wrapped'):
            return self.wrapped.proxy_for()
        else:
            return self.wrapped.__class__

    @property
    def pos(self):
        return self.wrapped.pos

    @pos.setter
    def pos(self, value):
        self.wrapped.pos = value

    def __getattr__(self, attr):
        return getattr(self.wrapped, attr)


def is_leaf(what):
    return hasattr(what, 'parts') and not what.children


def _create_proxy_compounds(real_thing, memo, particle_classes):
    proxy = Proxy(real_thing)
    memo[real_thing] = proxy

    if not type(real_thing) in particle_classes:
        if not is_leaf(real_thing):  # Recurse only if it has parts.
            # Recursively create proxies for parts.
            for part in real_thing.children:
                part_proxy = _create_proxy_compounds(part, memo,
                                                     particle_classes)
                proxy.add(part_proxy)

    return proxy


def _proxy_of(real_thing, memo):
    if real_thing in memo:
        return memo[real_thing]
    else:
        return _proxy_of(real_thing.parent, memo)


def _create_proxy_bonds(real_thing, memo, leaf_classes):
    proxy = memo[real_thing]

    if type(real_thing) in leaf_classes or is_leaf(real_thing):
        # It is a leaf of the proxy, so we don't recurse.
        pass
    else:
        for part in real_thing.children:
            _create_proxy_bonds(part, memo, leaf_classes)

    # Check if there's a contained bond that needs to be added to the proxy.
    if hasattr(real_thing, 'bonds'):
        for a1, a2 in real_thing.bonds():
            pa1 = _proxy_of(a1, memo)
            pa2 = _proxy_of(a2, memo)
            if pa1 != pa2:  # Do not add internal bonds.
                proxy.add_bond((pa1, pa2))


def _create_proxy_labels(real_thing, memo):
    if not is_leaf(real_thing):
        for label, part in real_thing.labels.items():
            if isinstance(part, list):
                # TODO support lists with labels
                continue
            if part in memo:
                memo[real_thing].labels[label] = memo[part]

        for part in real_thing.children:
            _create_proxy_labels(part, memo)


