__author__ = 'sallai'

import os.path

from compound import *
import unit as units

class Xyz(Compound):

    def __init__(self, path, ctx={}, labels=False, cwd="", distance_unit=units.angstroms):
        super(Xyz, self).__init__(ctx=ctx)
        self.DIST = distance_unit

        filename = os.path.join(cwd, path)
        with open(filename, 'r') as f:
            num_atoms = int(f.readline())
            comment = f.readline()
            for line in f:
                fields = line.split()
                atom_type = fields[0]
                #x = float(fields[1]) * self.DIST
                #y = float(fields[2]) * self.DIST
                #z = float(fields[3]) * self.DIST
                #atom_pos = (x, y, z)
                atom_pos = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if labels:
                    self.add(Atom(kind=atom_type, pos=atom_pos), atom_type + '_#')
                else:
                    self.add(Atom(kind=atom_type, pos=atom_pos))

    @staticmethod
    def save(compound, path, cwd="", print_ports=False):
        """
        Save into an xyz file
        :param fn: file name
        :param print_ports: if False, ghost points are not written
        """

        fn = os.path.join(cwd, path)

        with open(fn, 'w') as f:
            if print_ports:
                f.write(str(compound.atoms().__len__()) + '\n\n')
            else:
                i = 0
                for value in compound.atoms():
                    if value.kind != 'G':
                        i += 1
                f.write(str(i) + '\n\n')
            for value in compound.atoms():
                if print_ports:
                    f.write(value.kind + '\t' +
                            str(value.pos[0]) + '\t' +
                            str(value.pos[1]) + '\t' +
                            str(value.pos[2]) + '\n')
                else:
                    if value.kind != 'G':
                        f.write(value.kind + '\t' +
                                str(value.pos[0]) + '\t' +
                                str(value.pos[1]) + '\t' +
                                str(value.pos[2]) + '\n')


if __name__ == "__main__":
    m = Xyz(path="wip/xml/c60.xyz")
    print [a for a in m.atoms()]

    print m.boundingbox()

    from treeview import TreeView
    TreeView(m).show()
