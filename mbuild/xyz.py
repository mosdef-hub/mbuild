__author__ = 'sallai'

import time
import os.path
import sys
import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound

def load_xyz(filename, component=None, labels=False):

        if component is None:
            component = Compound()

        with open(filename, 'r') as f:
            num_atoms = int(f.readline())
            comment = f.readline()
            for line in f:
                fields = line.split()
                atom_type = fields[0]
                atom_pos = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                if labels:
                    component.add(Atom(kind=atom_type, pos=atom_pos), atom_type + '_#')
                else:
                    component.add(Atom(kind=atom_type, pos=atom_pos))

        return component


def write_xyz(compound, fn, print_ports=False):

    print "    Writing to '{0}'...".format(fn)
    start = time.time()
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
                    atom_name = value.kind
                    f.write(atom_name + '\t' +
                            str(value.pos[0]) + '\t' +
                            str(value.pos[1]) + '\t' +
                            str(value.pos[2]) + '\n')
    print "    Done. ({0:.2f} s)".format(time.time() - start)

if __name__ == "__main__":
    # Look for data file in same directory as this python module.
    current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
    new_path = os.path.join(current_dir, "../wip/resources/c60.xyz")


    m = load_xyz(new_path)
    print [a for a in m.atoms()]

    print m.boundingbox()

    # from treeview import TreeView
    # TreeView(m).show()
    from mbuild.plot import Plot
    Plot(m).show()
