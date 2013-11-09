__author__ = 'sallai'
from mbuild.compound import *


class Xyz(Compound):

    @classmethod
    def create(cls, fn, label=None):
        m = super(Xyz, cls).create(label)
        f = open(fn)
        num_atoms = int(f.readline())
        comment = f.readline()
        atom_cnt = {}
        for line in f:
            split_line = line.split()
            atom_type = split_line[0]
            atom_pos = (float(split_line[1]),float(split_line[2]),float(split_line[3]))
            if atom_type not in atom_cnt.keys():
                atom_label = atom_type + '0'
                atom_cnt[atom_type] = 1
            else:
                atom_label = atom_type + str(atom_cnt[atom_type])
                atom_cnt[atom_type] += 1

            m.add(Atom.create(atom_type, atom_pos), atom_label)
        return m

if __name__ == "__main__":
    m = Xyz.create("c60.xyz")
    # # print ethane
    print [label for label,a in m.atoms()]

    print m.boundingbox()

    m.plot(labels=False)
