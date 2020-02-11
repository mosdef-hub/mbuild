from collections import OrderedDict
from warnings import warn
import itertools as it

import numpy as np
from parmed.parameters import ParameterSet

from mbuild import Box
from mbuild.utils.conversion import RB_to_OPLS
from mbuild.utils.sorting import natural_sort
from scipy.constants import epsilon_0

__all__ = ['write_lammpsdata']

def write_lammpsdata(structure, filename, atom_style='full', 
                    unit_style='real',
                    detect_forcefield_style=True, nbfix_in_data_file=True,
                    use_urey_bradleys=False,
                    use_rb_torsions=True, use_dihedrals=False):
    """Output a LAMMPS data file.
    
    Outputs a LAMMPS data file in the 'full' atom style format. Default units are 
    'real' units. See http://lammps.sandia.gov/doc/atom_style.html for
    more information on atom styles.

    Parameters
    ----------
    structure : parmed.Structure
        ParmEd structure object
    filename : str
        Path of the output file
    atom_style: str
        Defines the style of atoms to be saved in a LAMMPS data file. The following atom
        styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more
        information on atom styles.
    unit_style: str
        Defines to unit style to be save in a LAMMPS data file.  Defaults to 'real' units.
        Current styles are supported: 'real', 'lj'
        see https://lammps.sandia.gov/doc/99/units.html for more information
        on unit styles
    detect_forcefield_style: boolean
        If True, format lammpsdata parameters based on the contents of 
        the parmed Structure
    use_urey_bradleys: boolean
        If True, will treat angles as CHARMM-style angles with urey bradley terms
        while looking for `structure.urey_bradleys`
    use_rb_torsions:
        If True, will treat dihedrals OPLS-style torsions while looking for
        `structure.rb_torsions`
    use_dihedrals:
        If True, will treat dihedrals as CHARMM-style dihedrals while looking for 
        `structure.dihedrals`

    Notes
    -----
    See http://lammps.sandia.gov/doc/2001/data_format.html for a full description
    of the LAMMPS data format. Currently the following sections are supported (in
    addition to the header): *Masses*, *Nonbond Coeffs*, *Bond Coeffs*, *Angle
    Coeffs*, *Dihedral Coeffs*, *Atoms*, *Bonds*, *Angles*, *Dihedrals*, *Impropers*
    OPLS and CHARMM forcefield styles are supported, AMBER forcefield styles are NOT

    Some of this function has beed adopted from `mdtraj`'s support of the LAMMPSTRJ
    trajectory format. See https://github.com/mdtraj/mdtraj/blob/master/mdtraj/formats/lammpstrj.py for details.

    """

    if atom_style not in ['atomic', 'charge', 'molecular', 'full']:
        raise ValueError('Atom style "{}" is invalid or is not currently supported'.format(atom_style))

    # Check if structure is paramterized
    if unit_style == 'lj':
        if any([atom.sigma for atom in structure.atoms]) is None:
           raise ValueError('LJ units specified but one or more atoms has undefined LJ parameters.') 

    xyz = np.array([[atom.xx,atom.xy,atom.xz] for atom in structure.atoms])

    forcefield = True
    if structure[0].type == '':
        forcefield = False

    """
    Note:
    -----
    unique_types : a sorted list of unique atomtypes for all atoms in the structure.
        Defined by:
            atomtype : atom.type
    unique_bond_types: an enumarated OrderedDict of unique bond types for all bonds in the structure.
        Defined by bond parameters and component atomtypes, in order:
            k : bond.type.k
            req : bond.type.req
            atomtypes : sorted((bond.atom1.type, bond.atom2.type))
    unique_angle_types: an enumerated OrderedDict of unique angle types for all angles in the structure.
        Defined by angle parameters and component atomtypes, in order:
            k : angle.type.k
            theteq : angle.type.theteq
            vertex atomtype: angle.atom2.type
            atomtypes: sorted((bond.atom1.type, bond.atom3.type))
    unique_dihedral_types: an enumerated OrderedDict of unique dihedrals type for all dihedrals in the structure.
        Defined by dihedral parameters and component atomtypes, in order:
            c0 : dihedral.type.c0
            c1 : dihedral.type.c1
            c2 : dihedral.type.c2
            c3 : dihedral.type.c3
            c4 : dihedral.type.c4
            c5 : dihedral.type.c5
            scee : dihedral.type.scee
            scnb : dihedral.type.scnb
            atomtype 1 : dihedral.atom1.type
            atomtype 2 : dihedral.atom2.type
            atomtype 3 : dihedral.atom3.type
            atomtype 4 : dihedral.atom4.type
    """
    if forcefield:
        types = [atom.type for atom in structure.atoms]
    else:
        types = [atom.name for atom in structure.atoms]

    unique_types = list(set(types))
    unique_types.sort(key=natural_sort)

    charges = np.array([atom.charge for atom in structure.atoms])

    # Convert coordinates to LJ units
    if unit_style == 'lj':
        # Get sigma, mass, and epsilon conversions by finding maximum of each
        sigma_conversion_factor = np.max([atom.sigma for atom in structure.atoms])
        epsilon_conversion_factor = np.max([atom.epsilon for atom in structure.atoms])
        mass_conversion_factor = np.max([atom.mass for atom in structure.atoms])

        xyz = xyz / sigma_conversion_factor
        charges = (charges*1.6021e-19) / np.sqrt(4*np.pi*(sigma_conversion_factor*1e-10)*
          (epsilon_conversion_factor*4184)*epsilon_0)
        charges[np.isinf(charges)] = 0 
        # TODO: FIX CHARGE UNIT CONVERSION
    else:
        sigma_conversion_factor = 1
        epsilon_conversion_factor = 1
        mass_conversion_factor = 1

    # Internally use nm
    box = Box(lengths=np.array([0.1 * val for val in structure.box[0:3]]),
              angles=structure.box[3:6])
    # Divide by conversion factor
    box.maxs /= sigma_conversion_factor
    
    # Lammps syntax depends on the functional form
    # Infer functional form based on the properties of the structure
    if detect_forcefield_style:
        # Check angles
        if len(structure.urey_bradleys) > 0 :
            print("Urey bradley terms detected, will use angle_style charmm")
            use_urey_bradleys = True
        else:
            print("No urey bradley terms detected, will use angle_style harmonic")
            use_urey_bradleys = False

        # Check dihedrals
        if len(structure.rb_torsions) > 0:
            print("RB Torsions detected, will use dihedral_style opls")
            use_rb_torsions = True
        else:
            use_rb_torsions = False
        if len(structure.dihedrals) > 0:
            print("Charmm dihedrals detected, will use dihedral_style charmm")
            use_dihedrals = True
        else:
            use_dihedrals = False
    if use_rb_torsions and use_dihedrals:
        raise ValueError("Multiple dihedral styles detected, check your "
                         "Forcefield XML and structure")

    # Check impropers
    for dihedral in structure.dihedrals:
        if dihedral.improper:
            raise ValueError("Amber-style impropers are currently not supported")

    bonds = [[bond.atom1.idx+1, bond.atom2.idx+1] for bond in structure.bonds]
    angles = [[angle.atom1.idx+1,
               angle.atom2.idx+1,
               angle.atom3.idx+1] for angle in structure.angles]
    if use_rb_torsions:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure.rb_torsions]
    elif use_dihedrals:
        dihedrals = [[dihedral.atom1.idx+1,
                      dihedral.atom2.idx+1,
                      dihedral.atom3.idx+1,
                      dihedral.atom4.idx+1] for dihedral in structure.dihedrals]
    else:
        dihedrals = []
    impropers = [[improper.atom1.idx+1,
                  improper.atom2.idx+1,
                  improper.atom3.idx+1,
                  improper.atom4.idx+1] for improper in structure.impropers]


    if bonds :
        if len(structure.bond_types) == 0:
            bond_types = np.ones(len(bonds),dtype=int)
        else:
            bond_types, unique_bond_types = _get_bond_types(structure,
                    bonds, sigma_conversion_factor,
                    epsilon_conversion_factor)

    if angles:
        angle_types, unique_angle_types = _get_angle_types(structure,
                use_urey_bradleys, sigma_conversion_factor,
                epsilon_conversion_factor)

    if dihedrals:
        dihedral_types, unique_dihedral_types = _get_dihedral_types(
                structure, use_rb_torsions, use_dihedrals,
                epsilon_conversion_factor)
            
    if impropers:
        improper_types, unique_improper_types = _get_impropers(structure,
                epsilon_conversion_factor)
    

    with open(filename, 'w') as data:
        data.write(filename+' - created by mBuild; units = {}\n\n'.format(
            unit_style))
        data.write('{:d} atoms\n'.format(len(structure.atoms)))
        if atom_style in ['full', 'molecular']:
            data.write('{:d} bonds\n'.format(len(bonds)))
            data.write('{:d} angles\n'.format(len(angles)))
            data.write('{:d} dihedrals\n'.format(len(dihedrals)))
            data.write('{:d} impropers\n\n'.format(len(impropers)))

        data.write('{:d} atom types\n'.format(len(set(types))))
        if atom_style in ['full', 'molecular']:
            if bonds:
                data.write('{:d} bond types\n'.format(len(set(bond_types))))
            if angles:
                data.write('{:d} angle types\n'.format(len(set(angle_types))))
            if dihedrals:
                data.write('{:d} dihedral types\n'.format(len(set(dihedral_types))))
            if impropers:
                data.write('{:d} improper types\n'.format(len(set(improper_types))))


        data.write('\n')
        # Box data
        if np.allclose(box.angles, np.array([90, 90, 90])):
            for i,dim in enumerate(['x','y','z']):
                data.write('{0:.6f} {1:.6f} {2}lo {2}hi\n'.format(
                    10.0 * box.mins[i],
                    10.0 * box.maxs[i],
                    dim))
        else:
            a, b, c = 10.0 * box.lengths
            alpha, beta, gamma = np.radians(box.angles)

            lx = a
            xy = b * np.cos(gamma)
            xz = c * np.cos(beta)
            ly = np.sqrt(b**2 - xy**2)
            yz = (b*c*np.cos(alpha) - xy*xz) / ly
            lz = np.sqrt(c**2 - xz**2 - yz**2)

            xlo, ylo, zlo = 10.0 * box.mins
            xhi = xlo + lx
            yhi = ylo + ly
            zhi = zlo + lz

            xlo_bound = xlo + np.min([0.0, xy, xz, xy+xz])
            xhi_bound = xhi + np.max([0.0, xy, xz, xy+xz])
            ylo_bound = ylo + np.min([0.0, yz])
            yhi_bound = yhi + np.max([0.0, yz])
            zlo_bound = zlo
            zhi_bound = zhi

            data.write('{0:.6f} {1:.6f} xlo xhi\n'.format(
                xlo_bound, xhi_bound))
            data.write('{0:.6f} {1:.6f} ylo yhi\n'.format(
                ylo_bound, yhi_bound))
            data.write('{0:.6f} {1:.6f} zlo zhi\n'.format(
                zlo_bound, zhi_bound))
            data.write('{0:.6f} {1:.6f} {2:6f} xy xz yz\n'.format(
                xy, xz, yz))

        # Mass data
        masses = np.array([atom.mass for atom in structure.atoms]) / mass_conversion_factor
        mass_dict = dict([(unique_types.index(atom_type)+1,mass) for atom_type,mass in zip(types,masses)])

        data.write('\nMasses\n\n')
        for atom_type,mass in mass_dict.items():
            data.write('{:d}\t{:.6f}\t# {}\n'.format(atom_type,mass,unique_types[atom_type-1]))

        if forcefield:
            epsilons = np.array([atom.epsilon for atom in structure.atoms]) / epsilon_conversion_factor
            sigmas = np.array([atom.sigma for atom in structure.atoms]) / sigma_conversion_factor
            forcefields = [atom.type for atom in structure.atoms]
            epsilon_dict = dict([(unique_types.index(atom_type)+1,epsilon) for atom_type,epsilon in zip(types,epsilons)])
            sigma_dict = dict([(unique_types.index(atom_type)+1,sigma) for atom_type,sigma in zip(types,sigmas)])
            forcefield_dict = dict([(unique_types.index(atom_type)+1,forcefield) for atom_type,forcefield in zip(types,forcefields)])
            


            # Modified cross-interactions
            if structure.has_NBFIX():
                params = ParameterSet.from_structure(structure)
                # Sort keys (maybe they should be sorted in ParmEd)
                new_nbfix_types = OrderedDict()
                for key, val in params.nbfix_types.items():
                    sorted_key = tuple(sorted(key))
                    if sorted_key in new_nbfix_types:
                        warn('Sorted key matches an existing key')
                        if new_nbfix_types[sorted_key]:
                            warn('nbfixes are not symmetric, overwriting old nbfix')
                    new_nbfix_types[sorted_key] = params.nbfix_types[key]
                params.nbfix_types = new_nbfix_types
                warn('Explicitly writing cross interactions using mixing rule: {}'.format(
                    structure.combining_rule))
                coeffs = OrderedDict()
                for combo in it.combinations_with_replacement(unique_types, 2):
                    # Attempt to find pair coeffis in nbfixes
                    if combo in params.nbfix_types:
                        type1 = unique_types.index(combo[0])+1
                        type2 = unique_types.index(combo[1])+1
                        rmin = params.nbfix_types[combo][0] # Angstrom OR lj units
                        epsilon = params.nbfix_types[combo][1] # kcal OR lj units
                        sigma = rmin/2**(1/6)
                        coeffs[(type1, type2)] = (round(sigma, 8), round(epsilon, 8))
                    else:
                        type1 = unique_types.index(combo[0]) + 1
                        type2 = unique_types.index(combo[1]) + 1
                        # Might not be necessary to be this explicit
                        if type1 == type2:
                            sigma = sigma_dict[type1]
                            epsilon = epsilon_dict[type1]
                        else:
                            if structure.combining_rule == 'lorentz':
                                sigma = (sigma_dict[type1]+sigma_dict[type2])*0.5
                            elif structure.combining_rule == 'geometric':
                                sigma = (sigma_dict[type1]*sigma_dict[type2])**0.5
                            else:
                                raise ValueError('Only lorentz and geometric combining rules are supported')
                            epsilon = (epsilon_dict[type1]*epsilon_dict[type2])**0.5
                        coeffs[(type1, type2)] = (round(sigma, 8), round(epsilon, 8))
                if nbfix_in_data_file:
                    data.write('\nPairIJ Coeffs # modified lj\n')
                    data.write('# type1 type2 \tepsilon (kcal/mol) \tsigma (Angstrom)\n')
                    for (type1, type2), (sigma, epsilon) in coeffs.items():
                        data.write('{0} \t{1} \t{2} \t\t{3}\t\t# {4}\t{5}\n'.format(
                            type1, type2, epsilon, sigma, forcefield_dict[type1], forcefield_dict[type2]))
                else:
                    data.write('\nPair Coeffs # lj\n\n')
                    for idx,epsilon in epsilon_dict.items():
                        data.write('{}\t{:.5f}\t{:.5f}\n'.format(idx,epsilon,sigma_dict[idx]))
                    print('Copy these commands into your input script:\n')
                    print('# type1 type2 \tepsilon (kcal/mol) \tsigma (Angstrom)\n')
                    for (type1, type2), (sigma, epsilon) in coeffs.items():
                        print('pair_coeff\t{0} \t{1} \t{2} \t\t{3} \t\t# {4} \t{5}'.format(
                            type1, type2, epsilon, sigma,forcefield_dict[type1],forcefield_dict[type2]))

            # Pair coefficients
            else:
                data.write('\nPair Coeffs # lj \n')
                if unit_style == 'real':
                    data.write('#\tepsilon (kcal/mol)\t\tsigma (Angstrom)\n')
                elif unit_style == 'lj':
                    data.write('#\treduced_epsilon \t\treduced_sigma \n')
                for idx,epsilon in epsilon_dict.items():
                    data.write('{}\t{:.5f}\t\t{:.5f}\t\t# {}\n'.format(idx,epsilon,sigma_dict[idx],forcefield_dict[idx]))

            # Bond coefficients
            if bonds:
                data.write('\nBond Coeffs # harmonic\n')
                if unit_style == 'real':
                    data.write('#\tk(kcal/mol/angstrom^2)\t\treq(angstrom)\n')
                elif unit_style == 'lj':
                    data.write('#\treduced_k\t\treduced_req\n')
                for params,idx in unique_bond_types.items():
                    data.write('{}\t{}\t\t{}\t\t# {}\t{}\n'.format(idx,params[0],params[1],params[2][0],params[2][1]))

            # Angle coefficients
            if angles:
                if use_urey_bradleys:
                    data.write('\nAngle Coeffs # charmm\n')
                    data.write('#\tk(kcal/mol/rad^2)\t\ttheteq(deg)\tk(kcal/mol/angstrom^2)\treq(angstrom)\n')
                    for params,idx in unique_angle_types.items():
                        data.write('{}\t{}\t{:.5f}\t{:.5f}\t{:.5f}\n'.format(idx,*params))

                else:
                    data.write('\nAngle Coeffs # harmonic\n')
                    data.write('#\treduced_k\t\ttheteq(deg)\n')
                    for params,idx in unique_angle_types.items():
                        data.write('{}\t{}\t\t{:.5f}\t# {}\t{}\t{}\n'.format(idx,params[0],params[1],
                                                                             params[3][0],params[2],params[3][1]))

            # Dihedral coefficients
            if dihedrals:
                if use_rb_torsions:
                    data.write('\nDihedral Coeffs # opls\n')
                    if unit_style == 'real':
                        data.write('#\tf1(kcal/mol)\tf2(kcal/mol)\tf3(kcal/mol)\tf4(kcal/mol)\n')
                    elif unit_style == 'lj':
                        data.write('#\tf1\tf2\tf3\tf4 (all lj reduced units)\n')
                    for params,idx in unique_dihedral_types.items():
                        opls_coeffs = RB_to_OPLS(params[0],
                                                 params[1],
                                                 params[2],
                                                 params[3],
                                                 params[4],
                                                 params[5])
                        data.write('{}\t{:.5f}\t{:.5f}\t\t{:.5f}\t\t{:.5f}\t# {}\t{}\t{}\t{}\n'.format(idx,opls_coeffs[0],
                                                                                                       opls_coeffs[1],
                                                                                                       opls_coeffs[2],
                                                                                                       opls_coeffs[3],
                                                                                                       params[8],params[9],
                                                                                                       params[10],params[11]))
                elif use_dihedrals:
                    data.write('\nDihedral Coeffs # charmm\n')
                    data.write('#k, n, phi, weight\n')
                    for params, idx in unique_dihedral_types.items():
                        data.write('{}\t{:.5f}\t{:d}\t{:d}\t{:.5f}\t# {}\t{}\t{}\t{}\n'.format(idx, params[0],
                                                                                                params[1], params[2],
                                                                                                params[3], params[6],
                                                                                                params[7], params[8], params[9]))

            # Improper coefficients
            if impropers:
                data.write('\nImproper Coeffs # harmonic\n')
                data.write('#k, phi\n')
                for params,idx in unique_improper_types.items():
                    data.write('{}\t{:.5f}\t{:.5f}\t# {}\t{}\t{}\t{}\n'.format(idx, params[0],
                                                                                params[1], params[2],
                                                                                params[3], params[4],
                                                                                params[5]))

        # Atom data
        data.write('\nAtoms\n\n')
        if atom_style == 'atomic':
            atom_line = '{index:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'charge':
            if unit_style == 'real':
                atom_line = '{index:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
            elif unit_style == 'lj':
                atom_line = '{index:d}\t{type_index:d}\t{charge:.4ef}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'molecular':
            atom_line = '{index:d}\t{zero:d}\t{type_index:d}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
        elif atom_style == 'full':
            if unit_style == 'real':
                atom_line ='{index:d}\t{zero:d}\t{type_index:d}\t{charge:.6f}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'
            elif unit_style == 'lj':
                atom_line ='{index:d}\t{zero:d}\t{type_index:d}\t{charge:.4e}\t{x:.6f}\t{y:.6f}\t{z:.6f}\n'

        for i,coords in enumerate(xyz):
            data.write(atom_line.format(
                index=i+1,type_index=unique_types.index(types[i])+1,
                zero=structure.atoms[i].residue.idx,charge=charges[i],
                x=coords[0],y=coords[1],z=coords[2]))

        if atom_style in ['full', 'molecular']:
            # Bond data
            if bonds:
                data.write('\nBonds\n\n')
                for i,bond in enumerate(bonds):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,bond_types[i],bond[0],bond[1]))

            # Angle data
            if angles:
                data.write('\nAngles\n\n')
                for i,angle in enumerate(angles):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,angle_types[i],angle[0],angle[1],angle[2]))

            # Dihedral data
            if dihedrals:
                data.write('\nDihedrals\n\n')
                for i,dihedral in enumerate(dihedrals):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,dihedral_types[i],dihedral[0],
                        dihedral[1],dihedral[2],dihedral[3]))
            # Dihedral data
            if impropers:
                data.write('\nImpropers\n\n')
                for i,improper in enumerate(impropers):
                    data.write('{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n'.format(
                        i+1,improper_types[i],improper[2],
                        improper[1],improper[0],improper[3]))

def _get_bond_types(structure, bonds, sigma_conversion_factor, 
        epsilon_conversion_factor):
    unique_bond_types = dict(enumerate(set([(round(bond.type.k*(
        sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                             round(bond.type.req/sigma_conversion_factor,3),
                                             tuple(sorted((bond.atom1.type,bond.atom2.type)))
                                             ) for bond in structure.bonds])))
    unique_bond_types = OrderedDict([(y,x+1) for x,y in unique_bond_types.items()])
    bond_types = [unique_bond_types[(round(bond.type.k*
        (sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                     round(bond.type.req/sigma_conversion_factor,3),
                                     tuple(sorted((bond.atom1.type,bond.atom2.type)))
                                     )] for bond in structure.bonds]
    return bond_types, unique_bond_types

def _get_angle_types(structure, use_urey_bradleys,
        sigma_conversion_factor, epsilon_conversion_factor):
    if use_urey_bradleys:
        charmm_angle_types = []
        for angle in structure.angles:
            ub_k = 0
            ub_req = 0
            for ub in structure.urey_bradleys:
                if (angle.atom1, angle.atom3) == (ub.atom1, ub.atom2):
                    ub_k = ub.type.k
                    ub_req = ub.type.req
            charmm_angle_types.append((round(angle.type.k*(
                sigma_conversion_factor**2/epsilon_conversion_factor),3), 
                                       round(angle.type.theteq,3),
                                       round(ub_k/epsilon_conversion_factor, 3),
                                       round(ub_req, 3),
                                       tuple(sorted((angle.atom1.type,angle.atom3.type)))))

        unique_angle_types = dict(enumerate(set(charmm_angle_types)))
        unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[ub_info] for ub_info in charmm_angle_types]

    else:
        unique_angle_types = dict(enumerate(set([(round(angle.type.k*(
            sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                                  round(angle.type.theteq,3),
                                                  angle.atom2.type,
                                                  tuple(sorted((angle.atom1.type,angle.atom3.type)))
                                                  ) for angle in structure.angles])))
        unique_angle_types = OrderedDict([(y,x+1) for x,y in unique_angle_types.items()])
        angle_types = [unique_angle_types[(round(angle.type.k*(
            sigma_conversion_factor**2/epsilon_conversion_factor),3),
                                           round(angle.type.theteq,3),
                                           angle.atom2.type,
                                           tuple(sorted((angle.atom1.type,angle.atom3.type)))
                                           )] for angle in structure.angles]

    return angle_types, unique_angle_types

def _get_dihedral_types(structure, use_rb_torsions, use_dihedrals,
         epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    if use_rb_torsions:
        unique_dihedral_types = dict(enumerate(set([(round(dihedral.type.c0*lj_unit,3),
                                                     round(dihedral.type.c1*lj_unit,3),
                                                     round(dihedral.type.c2*lj_unit,3),
                                                     round(dihedral.type.c3*lj_unit,3),
                                                     round(dihedral.type.c4*lj_unit,3),
                                                     round(dihedral.type.c5*lj_unit,3),
                                                     round(dihedral.type.scee,1),
                                                     round(dihedral.type.scnb,1),
                                                     dihedral.atom1.type, dihedral.atom2.type,
                                                     dihedral.atom3.type, dihedral.atom4.type
                                                     ) for dihedral in structure.rb_torsions])))
        unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])
        dihedral_types = [unique_dihedral_types[(round(dihedral.type.c0*lj_unit,3),
                                                 round(dihedral.type.c1*lj_unit,3),
                                                 round(dihedral.type.c2*lj_unit,3),
                                                 round(dihedral.type.c3*lj_unit,3),
                                                 round(dihedral.type.c4*lj_unit,3),
                                                 round(dihedral.type.c5*lj_unit,3),
                                                 round(dihedral.type.scee,1),
                                                 round(dihedral.type.scnb,1),
                                                 dihedral.atom1.type, dihedral.atom2.type,
                                                 dihedral.atom3.type, dihedral.atom4.type
                                                 )] for dihedral in structure.rb_torsions]
    elif use_dihedrals:
        charmm_dihedrals = []
        structure.join_dihedrals()
        for dihedral in structure.dihedrals:
            if not dihedral.improper:
                weight = 1 / len(dihedral.type)
                for dih_type in dihedral.type:
                    charmm_dihedrals.append((round(dih_type.phi_k*lj_unit,3),
                                             int(round(dih_type.per,0)),
                                             int(round(dih_type.phase,0)),
                                             round(weight, 4),
                                             round(dih_type.scee,1),
                                             round(dih_type.scnb,1),
                                             dihedral.atom1.type, dihedral.atom2.type,
                                             dihedral.atom3.type, dihedral.atom4.type))

        unique_dihedral_types = dict(enumerate(set(charmm_dihedrals)))
        unique_dihedral_types = OrderedDict([(y,x+1) for x,y in unique_dihedral_types.items()])
        dihedral_types = [unique_dihedral_types[dihedral_info] for dihedral_info in charmm_dihedrals]

    return dihedral_types, unique_dihedral_types

def _get_impropers(structure, epsilon_conversion_factor):
    lj_unit = 1 / epsilon_conversion_factor
    unique_improper_types = dict(enumerate(set([(round(improper.type.psi_k*lj_unit,3),
                                                 round(improper.type.psi_eq,3),
                                                 improper.atom1.type, improper.atom2.type,
                                                 improper.atom3.type, improper.atom4.type) for improper in structure.impropers])))
    unique_improper_types = OrderedDict([(y,x+1) for x,y in unique_improper_types.items()])
    improper_types = [unique_improper_types[(round(improper.type.psi_k*lj_unit,3),
                                             round(improper.type.psi_eq,3),
                                             improper.atom1.type, improper.atom2.type,
                                             improper.atom3.type, improper.atom4.type)] for improper in structure.impropers]

    return improper_types, unique_improper_types
