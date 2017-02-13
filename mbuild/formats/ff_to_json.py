from __future__ import division

__all__ = ['write_forcefield']


import re
import json
import parmed as pmd
import numpy as np
from parmed.periodic_table import AtomicNum
from .hoomdxml import RB_to_OPLS
from collections import OrderedDict
from oset import oset as OrderedSet

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
        pair_data = OrderedDict(sorted(pair_data.items(),
                                       key=lambda p: int(p[0].split('_')[1])))

        # Bonds
        bonds = [bond for bond in structure.bonds]
        if bonds:
            styles['bond'] = 'harmonic'
            unique_bond_types = dict(enumerate(OrderedSet([(round(bond.type.k,3),
                                                            round(bond.type.req,3)) for bond in structure.bonds])))
            bond_data = OrderedDict()
            for idx,key in unique_bond_types.items():
                temp_dict = OrderedDict()
                temp_dict['k'] = round(key[0]*2,3) * ((ref_distance**2)/ref_energy)
                temp_dict['r0'] = round(key[1],3) / ref_distance
                bond_data[str(idx)] = temp_dict

        # Angles
        angles = [angle for angle in structure.angles]
        if angles:
            styles['angle'] = 'harmonic'
            unique_angle_types = dict(enumerate(OrderedSet([(round(angle.type.k,3),
                                                    round(angle.type.theteq,3)) for angle in structure.angles])))
            angle_data = OrderedDict()
            for idx,key in unique_angle_types.items():
                temp_dict = OrderedDict()
                temp_dict['k'] = round(key[0]*2,3) / ref_energy
                temp_dict['t0'] = round(key[1],3) * (np.pi/180)
                angle_data[str(idx)] = temp_dict

        # Dihedrals
        dihedrals = [dihedral for dihedral in structure.rb_torsions]
        if dihedrals:
            styles['dihedral'] = 'opls'

            unique_dihedral_types = dict(enumerate(OrderedSet([(round(dihedral.type.c0,3),
                                                    round(dihedral.type.c1,3),
                                                    round(dihedral.type.c2,3),
                                                    round(dihedral.type.c3,3),
                                                    round(dihedral.type.c4,3),
                                                    round(dihedral.type.c5,3),
                                                    round(dihedral.type.scee,1),
                                                    round(dihedral.type.scnb,1)) for dihedral in structure.rb_torsions])))
            dihedrals_opls = [RB_to_OPLS(y[0],y[1],y[2],y[3],y[4],y[5]) for x,y in unique_dihedral_types.items()]
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
        if bonds:
            ff_data['bond_coeffs'] = bond_data
        if angles:
            ff_data['angle_coeffs'] = angle_data
        if dihedrals:
            ff_data['dihedral_coeffs'] = dihedral_data

        json.dump(ff_data,f,indent=4)
