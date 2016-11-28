from __future__ import division

__all__ = ['write_forcefield']


import json
import parmed as pmd
from parmed.periodic_table import AtomicNum
from .hoomdxml import RB_to_OPLS
from collections import OrderedDict

def write_forcefield(structure, filename):

    params = pmd.ParameterSet.from_structure(structure)

    ff_data = OrderedDict()
    ff_data['styles'] = OrderedDict([('pair','lj'),
                                     ('bond','harmonic'),
                                     ('angle','harmonic'),
                                     ('dihedral','opls')])

    with open(filename, 'w') as f:

        # Pair
        pair_data = OrderedDict()
        for key in params.atom_types.items():
            temp_dict = OrderedDict()
            temp_dict['element'] = {enum:ename for ename,enum in AtomicNum.items()}[key[1].atomic_number]
            temp_dict['epsilon'] = round(key[1].epsilon,3)
            temp_dict['sigma'] = round(key[1].sigma,3)
            pair_data[key[0]] = temp_dict

        ff_data['pair_coeffs'] = pair_data

        # Bonds
        bond_data = OrderedDict()
        for key in structure.bond_types:
            temp_dict = OrderedDict()
            temp_dict['k'] = round(key.k,3)
            temp_dict['req'] = round(key.req,3)
            bond_data[str(key.idx)] = temp_dict

        ff_data['bond_coeffs'] = bond_data

        # Angles
        angle_data = OrderedDict()
        for key in structure.angle_types:
            temp_dict = OrderedDict()
            temp_dict['k'] = round(key.k,3)
            temp_dict['theteq'] = round(key.theteq,3)
            angle_data[str(key.idx)] = temp_dict

        ff_data['angle_coeffs'] = angle_data

        # Dihedrals
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
            temp_dict['k0'] = round(key[0],3)
            temp_dict['k1'] = round(key[1],3)
            temp_dict['k2'] = round(key[2],3)
            temp_dict['k3'] = round(key[3],3)
            dihedral_data[str(idx)] = temp_dict
            
        ff_data['dihedral_coeffs'] = dihedral_data

        json.dump(ff_data,f,indent=4)
