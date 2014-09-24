from future.builtins import range

import operator

from mdtraj.utils import in_units_of


def save_lammps(traj, step=-1, filename='data.mbuild', unit_set='real'):
    """Output a Trajectory as a LAMMPS data file.

    Args:
        traj (Trajectory): The Trajectory to be output.
        step (int, optional): The frame in traj to save. Default is last frame.
        filename (str, optional): Path of the output file.
        unit_set (str, optional): The LAMMPS unit set to write the data file in.

    """
    _radians_unit = 'radians'
    _degrees_unit = 'degrees'
    if unit_set == 'real':
        _distance_unit = 'angstroms'
        _velocity_unit = 'angstroms/femtosecond'
        _energy_unit = 'kilocalories/mole'
        _mass_unit = 'grams/mole'
        _charge_unit = 'elementary_charge'
        _mole_unit = 'mole'
    else:
        raise Exception("Unsupported unit set specified: {0}".format(unit_set))

    directives_to_write = list()
    mass_list = list()
    mass_list.append('\n')
    mass_list.append('Masses\n')
    mass_list.append('\n')

    atom_list = list()
    atom_list.append('\n')
    atom_list.append('Atoms\n')
    atom_list.append('\n')

    numeric_types = dict()
    atom_type_n = 1
    for chain in traj.top.chains:
        for atom in chain.atoms:
            if atom.name not in numeric_types:
                numeric_types[atom.name] = atom_type_n
                mass = in_units_of(atom.element.mass, 'grams/moles', _mass_unit)
                mass_list.append('{0:d} {1:8.4f}\n'.format(atom_type_n, mass))
                atom_type_n += 1
            x, y, z = in_units_of(traj.xyz[step][atom.index], 'nanometers',
                                  _distance_unit)
            entry = '{0:-d} {1:-d} {2:-d} {3:5.8f} {4:8.5f} {5:8.5f} {6:8.5f}  # {7}\n'.format(
                atom.index + 1, chain.index + 1, numeric_types[atom.name],
                0.0, x, y, z, atom.name)
            atom_list.append(entry)

    directives_to_write.append(mass_list)
    directives_to_write.append(atom_list)
    print '(string: numeric) types for atoms'
    print sorted(numeric_types.iteritems(), key=operator.itemgetter(1))

    optional_directives = [('bonds', 2), ('angles', 3),
                           ('dihedrals', 4), ('impropers', 4)]
    number_of_terms = dict()
    for directive, n_terms in optional_directives:
        number_of_terms[directive] = 0
        numeric_types = dict()
        if getattr(traj.top, '_ff_{0}'.format(directive)):
            list_of_terms = list()
            list_of_terms.append('\n')
            list_of_terms.append('{0}\n'.format(directive.title()))
            list_of_terms.append('\n')

            term_type_n = 1
            for term_n, term in enumerate(getattr(traj.top, 'ff_{0}'.format(directive))):
                if term.kind not in numeric_types:
                    numeric_types[term.kind] = term_type_n
                    term_type_n += 1
                entry = '{0:-d} {1:d} '.format(term_n + 1, numeric_types[term.kind])
                for n in range(n_terms):
                    entry += '{0:d} '.format(getattr(term, 'atom{0}'.format(n + 1)).index + 1)
                entry += ' # {0}\n'.format(term.kind)
                list_of_terms.append(entry)
            number_of_terms[directive] = term_n + 1
            directives_to_write.append(list_of_terms)
        print '(string: numeric types) for {0}'.format(directive)
        print sorted(numeric_types.iteritems(), key=operator.itemgetter(1))

    with open(filename, 'w') as f:
        n_bonds = number_of_terms['bonds']
        n_angles = number_of_terms['angles']
        n_dihedrals = number_of_terms['dihedrals']
        n_impropers = number_of_terms['impropers']

        f.write('{0} atoms\n'.format(traj.n_atoms))
        f.write('{0} bonds\n'.format(n_bonds))
        f.write('{0} angles\n'.format(n_angles))
        f.write('{0} dihedrals\n'.format(n_dihedrals))
        f.write('{0} impropers\n'.format(n_impropers))
        f.write('\n')

        box = traj.boundingbox(step)
        box.mins = in_units_of(box.mins, 'nanometers', _distance_unit)
        box.maxes = in_units_of(box.mins, 'nanometers', _distance_unit)
        f.write(
            '{0:10.6f} {1:10.6f} xlo xhi\n'.format(box.mins[0], box.maxes[0]))
        f.write(
            '{0:10.6f} {1:10.6f} ylo yhi\n'.format(box.mins[1], box.maxes[1]))
        f.write(
            '{0:10.6f} {1:10.6f} zlo zhi\n'.format(box.mins[2], box.maxes[2]))

        for directive in directives_to_write:
            for entry in directive:
                f.write(entry)


if __name__ == "__main__":
    from mbuild.examples.ethane.ethane import Ethane

    ethane = Ethane()
    ethane = ethane.to_trajectory()
    ethane.top.find_forcefield_terms()
    save_lammps(ethane, filename='data.ethane')
