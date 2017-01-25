from __future__ import division

__all__ = ['write_forcefield']


import re
import json
import parmed as pmd
import numpy as np
from parmed.periodic_table import AtomicNum
from .hoomdxml import RB_to_OPLS
from collections import OrderedDict

def write_forcefield(structure, filename, ref_distance=1.0, ref_energy=1.0):
    """Output force field information in JSON format.

    Parameters
    ----------
    structure : parmed.GromacsTopologyFile
        Parmed structure object
    filename : str
        Path of the output file.
    ref_distance : float, default=1.0
        Reference distance for conversion to reduced units
    ref_energy : float, default=1.0
        Reference energy for conversion to reduced units
    """

    params = pmd.ParameterSet.from_structure(structure)
    ff_data = OrderedDict()

    styles = OrderedDict()

    with open(filename, 'w') as f:

        # Pair
        charges = np.any([atom.charge for atom in structure.atoms])
        if charges:
            styles['pair'] = 'lj/coul'
        else:
            styles['pair'] = 'lj'

        pair_data = dict()
        for key in params.atom_types.items():
            temp_dict = OrderedDict()
            temp_dict['element'] = {enum:ename for ename,enum in AtomicNum.items()}[key[1].atomic_number]
            temp_dict['epsilon'] = round(key[1].epsilon,3) / ref_energy
            temp_dict['sigma'] = round(key[1].sigma,3) / ref_distance
            pair_data[key[0]] = temp_dict
        pair_data = OrderedDict(sorted(pair_data.items(), key=lambda(key,value):int(key.split('_')[1])))

        # Bonds
        bonds = [bond for bond in structure.bonds]
        if bonds:
            styles['bond'] = 'harmonic'
            bond_data = OrderedDict()
            for key in structure.bond_types:
                temp_dict = OrderedDict()
                temp_dict['k'] = round(key.k*2,3) * ((ref_distance**2)/ref_energy)
                temp_dict['r0'] = round(key.req,3) / ref_distance
                bond_data[str(key.idx)] = temp_dict

        # Angles
        angles = [angle for angle in structure.angles]
        if angles:
            styles['angle'] = 'harmonic'
            angle_data = OrderedDict()
            for key in structure.angle_types:
                temp_dict = OrderedDict()
                temp_dict['k'] = round(key.k*2,3) / ref_energy
                temp_dict['t0'] = round(key.theteq,3) * (np.pi/180)
                angle_data[str(key.idx)] = temp_dict

        # Dihedrals
        dihedrals = [dihedral for dihedral in structure.rb_torsions]
        if dihedrals:
            styles['dihedral'] = 'opls'
            dihedrals = structure.rb_torsion_types
            for i,dihedral in enumerate(dihedrals):
                for j in range(i+1, len(dihedrals)):
                    dihedral2 = dihedrals[j]
                    if dihedral == dihedral2:
                        dihedrals[j] = dihedral
            dihedrals = list(set(dihedrals))
            dihedrals_opls = [RB_to_OPLS(dihedral.c0,
                                         dihedral.c1,
                                         dihedral.c2,
                                         dihedral.c3,
                                         dihedral.c4,
                                         dihedral.c5) for dihedral in dihedrals]
            dihedral_data = OrderedDict()
            for idx,key in enumerate(dihedrals_opls):
                temp_dict = OrderedDict()
                temp_dict['k1'] = round(key[0],3) / ref_energy
                temp_dict['k2'] = round(key[1],3) / ref_energy
                temp_dict['k3'] = round(key[2],3) / ref_energy
                temp_dict['k4'] = round(key[3],3) / ref_energy
                dihedral_data[str(idx)] = temp_dict
            
        ff_data['styles'] = styles
        ff_data['pair_coeffs'] = pair_data
        ff_data['bond_coeffs'] = bond_data
        ff_data['angle_coeffs'] = angle_data
        ff_data['dihedral_coeffs'] = dihedral_data

        json.dump(ff_data,f,indent=4)
