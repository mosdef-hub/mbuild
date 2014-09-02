from mbuild.flat_compound import FlatCompound
from mbuild.plugins.system import FlatCompound

__author__ = 'sallai'

import time
import os.path
import sys

import numpy as np

from mbuild.atom import Atom
from mbuild.compound import Compound

def load_xyz(filename):
    """ """

    coords = None
    types = None

    with open(filename, 'r') as f:
        n_atoms = int(f.readline())
        comment = f.readline()
        coords = np.ndarray(shape=(n_atoms,3), dtype='float')
        types = np.empty(n_atoms, dtype='string')
        idx = 0
        for line in f:
            fields = line.split()
            types[idx] = fields[0]
            # coords[idx,0] = float(fields[1])
            # coords[idx,1] = float(fields[2])
            # coords[idx,2] = float(fields[3])
            coords[idx] = np.array([float(fields[1]),float(fields[2]),float(fields[3])])
            idx += 1

    return FlatCompound(coords=coords, types=types)


def write_xyz(system, fn):
    """ """
    print "    Writing to '{0}'...".format(fn)
    start = time.time()
    with open(fn, 'w') as f:
        f.write(str(system.n_atoms) + '\n\n')
        for idx, kind in enumerate(system.types):
                f.write(kind + '\t' +
                        str(system.coords[idx, 0]) + '\t' +
                        str(system.coords[idx, 1]) + '\t' +
                        str(system.coords[idx, 2]) + '\n')
    print "    Done. ({0:.2f} s)".format(time.time() - start)

if __name__ == "__main__":
    # Look for data file in same directory as this python module.
    current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
    new_path = os.path.join(current_dir, "../../wip/resources/c60.xyz")

    sys = load_xyz(new_path)
    print sys.coords

    m = sys.to_compound()

    # # from treeview import TreeView
    # # TreeView(m).show()
    from mbuild.plot import Plot
    Plot(m).show()
