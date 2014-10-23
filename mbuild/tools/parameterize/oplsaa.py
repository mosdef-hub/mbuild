from collections import defaultdict
import mdtraj
from networkx.algorithms.isomorphism.vf2userfunc import GraphMatcher
from mbuild.examples.alkane.alkane import Alkane
from mbuild.examples.ethane.methyl import Methyl



def opls_atomtypes(compound):

    for atom in compound.yield_atoms():
        if atom.kind == 'G':
            continue
        elif atom.kind == 'C':
            carbon(atom)
        elif atom.kind == 'H':
            hydrogen(atom)
        else:
            print ('atom kind {} not supported'.format(atom.kind))




def carbon(atom):
    valency = len(atom.bonds)

    assert valency<5, 'carbon with valency {}'.format(valency)

    if valency == 4:
        # alkanes

        # CH3
        if neighbor_types(atom)['H'] == 3 and neighbor_types(atom)['C'] == 1:
            atom.extras['opls_type'] = '135'

        # CH2
        elif neighbor_types(atom)['H'] == 2 and neighbor_types(atom)['C'] == 2:
            atom.extras['opls_type'] = '136'

        # CH
        elif neighbor_types(atom)['H'] == 1 and neighbor_types(atom)['C'] == 3:
            atom.extras['opls_type'] = '137'

        # CH4
        elif neighbor_types(atom)['H'] == 4:
            atom.extras['opls_type'] = '138'

        # C
        elif neighbor_types(atom)['C'] == 4:
            atom.extras['opls_type'] = '139'

        else:
            print "unidentified {}-valent carbon with bonds: {}".format(valency, atom.bonds)

    else:
        print "unidentified {}-valent carbon with bonds: {}".format(valency, atom.bonds)


def hydrogen(atom):
    valency = len(atom.bonds)

    assert valency<2, 'hydrogen with valency {}'.format(valency)

    if neighbor_types(atom)['C'] == 1:
        atom.extras['opls_type'] = '140'

    else:
        print "unidentified {}-valent hydrogen with bonds: {}".format(valency, atom.bonds)


neighbor_types_map = {}

def neighbor_types(atom):

    if atom in neighbor_types_map:
        return neighbor_types_map[atom]
    else:

        rval = defaultdict(int)
        for b in atom.bonds:
            kind = b.other_atom(atom).kind
            rval[kind] += 1

        neighbor_types_map[atom] = rval

    return rval




if __name__ == "__main__":
    alkane = Alkane(n=3)
    opls_atomtypes(alkane)

    for atom in alkane.yield_atoms():
        print "atom kind={} opls_type={}".format(atom.kind, atom.extras['opls_type'])





