__author__ = 'sallai'

from .treeview import TreeView
from .compound import *
import os.path

class Xyz(Compound):

    def __init__(self, path, ctx={}, labels=False, cwd=""):
        super(Xyz, self).__init__(ctx=ctx)

        fn = os.path.join(cwd, path)

        f = open(fn)
        num_atoms = int(f.readline())
        comment = f.readline()
        for line in f:
            split_line = line.split()
            atom_type = split_line[0]
            atom_pos = (float(split_line[1]),float(split_line[2]),float(split_line[3]))
            if labels:
                self.add(Atom(kind=atom_type, pos=atom_pos), atom_type + '_#')
            else:
                self.add(Atom(kind=atom_type, pos=atom_pos))
        f.close()

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
    m = Xyz(path="c60.xyz")
    print [a for a in m.atoms()]

    print m.boundingbox()

    TreeView(m).show()
