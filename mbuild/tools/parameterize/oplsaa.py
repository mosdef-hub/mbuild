from collections import defaultdict
from copy import deepcopy, copy
import os
import mdtraj
from networkx.algorithms.isomorphism.vf2userfunc import GraphMatcher
from mbuild.compound import Compound
from mbuild.examples.alkane.alkane import Alkane
from mbuild.examples.ethane.methyl import Methyl
from pkg_resources import resource_filename
from mbuild.orderedset import OrderedSet
import sys


# map opls ids to the functions that check for them
rules = dict()

def opls_atomtypes(compound):

    for fn, fcn in sys.modules[__name__].__dict__.items():
        if fn.startswith('opls_'):
            rules[fn.split("_")[1]] = fcn

    for atom in compound.yield_atoms():
        prepare(atom)

    max_iter = 10
    iter_cnt = 0
    while(True):

        print ("Iteration {}".format(iter_cnt))
        old_len = 0
        new_len = 0
        for atom in compound.yield_atoms():

            old_len += len(atom.opls_whitelist)
            old_len += len(atom.opls_blacklist)

            if atom.kind == 'G':
                continue
            elif atom.kind == 'C':
                carbon(atom)
            elif atom.kind == 'H':
                hydrogen(atom)
            else:
                print ('atom kind {} not supported'.format(atom.kind))

            new_len += len(atom.opls_whitelist)
            new_len += len(atom.opls_blacklist)

        if old_len == new_len:
            break

        iter_cnt += 1
        if max_iter == iter_cnt:
            print("reached maximum iteration count")
            break





def prepare(atom):
    atom.extras['opls_whitelist'] = OrderedSet()
    atom.extras['opls_blacklist'] = OrderedSet()


def whitelist(atom, rule):
    if isinstance(rule, (list, tuple, set)):
        for r in rule:
            atom.opls_whitelist.add(str(r))
    else:
        atom.opls_whitelist.add(str(rule))

def blacklist(atom, rule):
    if isinstance(rule, (list, tuple, set)):
        for r in rule:
            atom.opls_blacklist.add(str(r))
    else:
        atom.opls_blacklist.add(str(rule))


# alkanes
def opls_135(atom):
    # CH3
    if neighbor_types(atom)['H'] == 3 and neighbor_types(atom)['C'] == 1:
        whitelist(atom, 135)
        blacklist(atom, [136, 137, 138, 139])


def opls_136(atom):
    # CH2
    if neighbor_types(atom)['H'] == 2 and neighbor_types(atom)['C'] == 2:
        whitelist(atom, 136)
        blacklist(atom, [135, 137, 138, 139])

def opls_137(atom):
    # CH
    if neighbor_types(atom)['H'] == 1 and neighbor_types(atom)['C'] == 3:
        whitelist(atom, 137)
        blacklist(atom, [135, 136, 138, 139])

def opls_138(atom):
    # CH4
    if neighbor_types(atom)['H'] == 4:
        whitelist(atom, 138)
        blacklist(atom, [135, 136, 137, 139])

def opls_139(atom):
    # C
    if neighbor_types(atom)['C'] == 4:
        whitelist(atom, 139)
        blacklist(atom, [135, 136, 137, 138])

def opls_140(atom):
    if neighbor_types(atom)['C'] == 1:
        whitelist(atom, 140)


# alkenes
def opls_141(atom):
    # opls_141   12.01100  ; alkene C (R2-C=)
    # TODO: check notation
    if neighbor_types(atom)['C'] == 3:
        whitelist(atom, 141)
        blacklist(atom, [142, 143])

def opls_142(atom):
    # opls_142   12.01100  ; alkene C (RH-C=)
    if neighbor_types(atom)['C'] == 2 and neighbor_types(atom)['H'] == 1:
        whitelist(atom, 142)
        blacklist(atom, [141, 143])

def opls_143(atom):
    # opls_143   12.01100  ; alkene C (H2-C=)
    if neighbor_types(atom)['C'] == 1 and neighbor_types(atom)['H'] == 2:
        whitelist(atom, 143)
        blacklist(atom, [141, 142])

def opls_144(atom):
    # opls_144    1.00800  ; alkene H (H-C=)

    if neighbor_types(atom)['C'] == 1:
        c = atom.neighbors[0]
        rule_ids = set(['141', '142', '143'])
        rule_ids.intersection_update(c.opls_whitelist)
        rule_ids.difference_update(c.opls_blacklist)

        if rule_ids:
            whitelist(atom, 144)
            blacklist(atom, 140)

# benzene
def opls_145(atom):
    # opls_145   12.01100  ; Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl

    if neighbor_types(atom)['C'] == 2 and neighbor_types(atom)['H'] == 1:
        r = Rings(atom, 6)

        # 2 rings, because we count the traversal in both directions
        if len(r.rings) == 2:
            for c in r.rings[0]:
                if not (c.kind == 'C' and len(c.neighbors) == 3):
                    break
            else:
                whitelist(atom, 145)
                # blacklist alkene carbons (with valency 3)
                blacklist(atom, [141, 142, 143])


def opls_145B(atom):
    #opls_145B  12.01100  ; Biphenyl C1

    if neighbor_types(atom)['C'] == 3:
        r = Rings(atom, 6)

        # 2 rings, because we count the traversal in both directions
        if len(r.rings) == 2:
            for c in r.rings[0]:
                if not (c.kind == 'C' and len(c.neighbors) == 3):
                    break
            else:
                for neighbor in atom.neighbors:
                    if neighbor not in r.rings[0]:
                        r = Rings(neighbor, 6)

                        # 2 rings, because we count the traversal in both directions
                        if len(r.rings) == 2:
                            for c in r.rings[0]:
                                if not (c.kind == 'C' and len(c.neighbors) == 3):
                                    break
                            else:
                                whitelist(atom, '145B')
                                # blacklist alkene carbons (with valency 3)
                                blacklist(atom, [141, 142, 143])
                                # blacklist benzene carbon
                                blacklist(atom, 145)

def opls_146(atom):
    #  opls_146    1.00800  ; Benzene H - 12 site.

    if neighbor_types(atom)['C'] == 1:
        c = atom.neighbors[0]
        rule_ids = set(['145'])
        rule_ids.intersection_update(c.opls_whitelist)
        rule_ids.difference_update(c.opls_blacklist)

        if rule_ids:
            whitelist(atom, 146)
            blacklist(atom, 144)
            blacklist(atom, 140)


def run_rule(atom, rule_id):
    if not rule_id in atom.opls_blacklist and not rule_id in atom.opls_whitelist:
        try:
            rule_fn = rules[rule_id]
        except KeyError:
            raise KeyError('rule {} not implemented'.format(rule_id))

        rule_fn(atom)


def carbon(atom):
    valency = len(atom.bonds)

    assert valency<5, 'carbon with valency {}'.format(valency)

    if valency == 4:
        for rule_id in ['135','136','137','138','139']:
            run_rule(atom, rule_id)
    elif valency == 3:
        for rule_id in ['141','142','143', '145', '145B']:
            run_rule(atom, rule_id)
    else:
        print "found no rules for {}-valent carbon".format(valency)




def hydrogen(atom):
    valency = len(atom.bonds)

    assert valency<2, 'hydrogen with valency {}'.format(valency)

    if valency == 1:
        for rule_id in ['140', '144', '146']:
            run_rule(atom, rule_id)

    else:
        print "found no rules for {}-valent hydrogen".format(valency)


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



def get_fn(name):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``mbuild/testing/reference``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Args:
        name (str): Name of the file to load (with respect to the reference/ folder).

    """
    fn = resource_filename('mbuild', os.path.join('..', 'opls_validation', name))

    if not os.path.exists(fn):
        raise ValueError('Sorry! %s does not exists. If you just '
                         'added it, you\'ll have to re install' % fn)

    return fn


class Rings(object):
    """ """
    def __init__(self, atom, ring_length):
        """ """
        self.rings = list()
        self.current_path = list()
        self.branch_points = OrderedSet()
        self.ring_length = ring_length
        self.find_rings(atom)

    def find_rings(self, atom):
        self.current_path = list()
        self.current_path.append(atom)
        self.step(atom)

    def step(self, atom):
        neighbors = atom.neighbors
        if len(neighbors) > 1:
            if len(neighbors) > 2:
                self.branch_points.add(atom)
            for n in neighbors:
                # Check to see if we found a ring.
                if len(self.current_path) > 2 and n == self.current_path[0]:
                    self.rings.append(copy(self.current_path))
                # Prevent stepping backwards.
                elif n in self.current_path:
                    continue
                else:
                    if len(self.current_path) < self.ring_length:
                        # Take another step.
                        self.current_path.append(n)
                        self.step(n)
                    else:
                        # Reached max length.
                        continue
            else:
                # Finished looping over all neighbors.
                del self.current_path[-1]
                if atom in self.branch_points:
                    self.branch_points.discard(atom)
        else:
            # Found a dead end.
            del self.current_path[-1]


if __name__ == "__main__":


    # m = Alkane(n=3)
    # opls_atomtypes(alkane)
    #
    # for atom in alkane.yield_atoms():
    #     print "atom kind={} opls_type={}".format(atom.kind, atom.extras['opls_type'])



    m = Compound()
    # m.append_from_file(get_fn('isopropane.pdb'))
    # m.append_from_file(get_fn('cyclohexane.pdb'))
    # m.append_from_file(get_fn('neopentane.pdb'))
    # m.append_from_file(get_fn('benzene.pdb'))
    # m.append_from_file(get_fn('1-propene.pdb'))
    m.append_from_file(get_fn('biphenyl.pdb'))
    opls_atomtypes(m)

    for atom in m.atoms:
        print "atom kind={} opls_whitelist={}  opls_black={}".format(atom.kind, atom.opls_whitelist, atom.opls_blacklist)




