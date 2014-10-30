from collections import defaultdict
from copy import copy
import os
from warnings import warn
from pkg_resources import resource_filename
import sys

from mbuild.compound import Compound
from mbuild.orderedset import OrderedSet

# Map opls ids to the functions that check for them.
rules = dict()

# Globally maintained neighbor information (see `neighbor_types()`).
neighbor_types_map = {}


def opls_atomtypes(compound):
    """Determine OPLS-aa atomtypes for all atoms in `compound`.

    This is where everything is orchestrated and the outer iteration happens.


    TODO: look into factoring out functions for different rings (see 145)
    """
    # Build a map to all of the supported opls_* functions.
    for fn, fcn in sys.modules[__name__].__dict__.items():
        if fn.startswith('opls_'):
            rules[fn.split("_")[1]] = fcn

    # Add white- and blacklists to all atoms.
    for atom in compound.yield_atoms():
        prepare(atom)

    max_iter = 10
    for iter_cnt in range(max_iter):
        #print ("Iteration {}".format(iter_cnt))

        # For comparing the lengths of the white- and blacklists.
        old_len = 0
        new_len = 0
        for atom in compound.yield_atoms():
            old_len += len(atom.opls_whitelist)
            old_len += len(atom.opls_blacklist)

            if atom.kind == 'G':  # Ignore Ports.
                continue
            elif atom.kind == 'C':
                carbon(atom)
            elif atom.kind == 'H':
                hydrogen(atom)
            elif atom.kind == 'O':
                oxygen(atom)
            else:
                warn("Atom kind '{}' not supported".format(atom.kind))

            new_len += len(atom.opls_whitelist)
            new_len += len(atom.opls_blacklist)

        # Nothing changed, we're done!
        if old_len == new_len:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")


def prepare(atom):
    """Add white- and blacklists to atom. """
    atom.extras['opls_whitelist'] = OrderedSet()
    atom.extras['opls_blacklist'] = OrderedSet()


def run_rule(atom, rule_id):
    """Execute the rule function for a specified OPLS-aa atomtype. """
    if not rule_id in atom.opls_blacklist and not rule_id in atom.opls_whitelist:
        try:
            rule_fn = rules[str(rule_id)]
        except KeyError:
            raise KeyError('Rule for {} not implemented'.format(rule_id))
        rule_fn(atom)


def neighbor_types(atom):
    """Maintains number of neighbors of each element type for all atoms.

    The dict maintained is `neighbor_types_map` and is organized as follows:
        atom: defaultdict{element: number of neighbors of that element type}

    E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
        Atom: {'C': 3, 'H': 1}
    """
    if atom in neighbor_types_map:
        return neighbor_types_map[atom]
    else:
        rval = defaultdict(int)
        for b in atom.bonds:
            kind = b.other_atom(atom).kind
            rval[kind] += 1
        neighbor_types_map[atom] = rval
    return rval


def check_atom(neighbor, input_rule_ids):
    """Ensure that atom is valid candidate.

    Checks that every rule in `rule_ids` is in the white- and not the blacklist.
    """
    rule_ids = set()
    if isinstance(input_rule_ids, (list, tuple, set)):
        for r in input_rule_ids:
            rule_ids.add(str(r))
    else:
        rule_ids.add(str(input_rule_ids))
    rule_ids.intersection_update(neighbor.opls_whitelist)
    rule_ids.difference_update(neighbor.opls_blacklist)
    return rule_ids


def whitelist(atom, rule):
    """Whitelist an OPLS-aa atomtype for an atom. """
    if isinstance(rule, (list, tuple, set)):
        for r in rule:
            atom.opls_whitelist.add(str(r))
    else:
        atom.opls_whitelist.add(str(rule))


def blacklist(atom, rule):
    """Blacklist an OPLS-aa atomtype for an atom. """
    if isinstance(rule, (list, tuple, set)):
        for r in rule:
            atom.opls_blacklist.add(str(r))
    else:
        atom.opls_blacklist.add(str(rule))


def get_opls_fn(name):
    """Get the full path to a file used to validate the OPLS-aa atomtyper.

    In the source distribution, these files are in ``opls_validation``.

    Args:
        name (str): Name of the file to load.
    """
    fn = resource_filename('mbuild',
                           os.path.join('..', 'opls_validation', name))
    if not os.path.exists(fn):
        raise ValueError('Sorry! {} does not exists. If you just '
                         'added it, you\'ll have to re install'.format(fn))
    return fn


def no_pattern(atom, valency):
    warn("No connectivity patterns matched {} with valency {} and neighbors "
         "{}.".format(atom, valency, neighbor_types(atom).items()))


def no_rule(atom, valency):
    warn("Found no rules for {}-valent {}.".format(valency, atom))


class Rings(object):
    """Find all rings of a specified length that the atom is a part of.

    Note: Finds each ring twice because the graph is traversed in both directions.
    """
    def __init__(self, atom, ring_length):
        """Initialize a ring bearer. """
        self.rings = list()
        self.current_path = list()
        self.branch_points = OrderedSet()
        self.ring_length = ring_length
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

#-----------------------------------------------#
# Filters for each element to break up the code #
#-----------------------------------------------#


def carbon(atom):
    valency = len(atom.bonds)
    assert valency < 5, 'Found carbon with valency {}.'.format(valency)

    if valency == 4:
        if neighbor_types(atom)['H'] == 4:
            for rule_id in [138]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['H'] == 3 and neighbor_types(atom)['C'] == 1:
            for rule_id in [135, 148]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['H'] == 2 and neighbor_types(atom)['C'] == 2:
            for rule_id in [136, 149, 218]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['H'] == 1 and neighbor_types(atom)['C'] == 3:
            for rule_id in [137]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['C'] == 4:
            for rule_id in [139]:
                run_rule(atom, rule_id)
        elif (neighbor_types(atom)['C'] == 1 and
              neighbor_types(atom)['O'] == 1 and
              neighbor_types(atom)['H'] == 2):
            for rule_id in [218]:
                run_rule(atom, rule_id)
        else:
            no_pattern(atom, valency)
    elif valency == 3:
        if neighbor_types(atom)['H'] == 2 and neighbor_types(atom)['C'] == 1:
            for rule_id in [143]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['H'] == 1 and neighbor_types(atom)['C'] == 2:
            for rule_id in [142, 145]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['C'] == 3:
            for rule_id in [141, 145, '145B', 221]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['C'] >= 2:
            for rule_id in [145]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['C'] >= 1 and neighbor_types(atom)['O'] == 1:
            for rule_id in [232]:
                run_rule(atom, rule_id)
        else:
            no_pattern(atom, valency)
    else:
        no_rule(atom, valency)


def hydrogen(atom):
    valency = len(atom.bonds)
    assert valency < 2, 'Found hydrogen with valency {}'.format(valency)

    if valency == 1:
        if neighbor_types(atom)['C'] == 1:
            for rule_id in [140, 144, 146, 279]:
                run_rule(atom, rule_id)
        elif neighbor_types(atom)['O'] == 1:
            for rule_id in [155]:
                run_rule(atom, rule_id)
        else:
            no_pattern(atom, valency)
    else:
        no_rule(atom, valency)


def oxygen(atom):
    valency = len(atom.bonds)
    # TODO: check if OPLS has parameters for things like hydronium
    assert valency < 3, 'Found oxygen with valency {}.'.format(valency)

    if valency == 2:
        if neighbor_types(atom)['H'] == 1 and neighbor_types(atom)['C'] == 1:
            for rule_id in [154]:
                run_rule(atom, rule_id)
        else:
            no_pattern(atom, valency)
    elif valency == 1:
        if neighbor_types(atom)['C'] == 1:
            for rule_id in [278]:
                run_rule(atom, rule_id)
        else:
            no_pattern(atom, valency)
    else:
        no_rule(atom, valency)

#---------------------------------------------------------#
# Filters for some specific patterns to break up the code #
#---------------------------------------------------------#


def benzene(atom):
    """Check if atom is part of a single benzene ring. """
    benzene = Rings(atom, 6).rings
    # 2 rings, because we count the traversal in both directions.
    if len(benzene) == 2:
        for c in benzene[0]:
            if not (c.kind == 'C' and len(c.neighbors) == 3):
                break
        else:
            return benzene[0]  # Only return one direction of the ring.
    return False

#----------------#
# House of rules #
#----------------#


# Alkanes
def opls_135(atom):
    # alkane CH3
    whitelist(atom, 135)


def opls_136(atom):
    # alkane CH2
    whitelist(atom, 136)


def opls_137(atom):
    # alkane CH
    whitelist(atom, 137)


def opls_138(atom):
    # alkane CH4
    whitelist(atom, 138)


def opls_139(atom):
    # alkane C
    whitelist(atom, 139)


def opls_140(atom):
    # alkane H
    whitelist(atom, 140)


# Alkenes
def opls_141(atom):
    # alkene C (R2-C=)
    whitelist(atom, 141)


def opls_142(atom):
    # alkene C (RH-C=)
    whitelist(atom, 142)


def opls_143(atom):
    # alkene C (H2-C=)
    whitelist(atom, 143)


def opls_144(atom):
    # alkene H (H-C=)
    # Make sure that the carbon is an alkene carbon.
    rule_ids = [141, 142, 143]
    if check_atom(atom.neighbors[0], rule_ids):
        whitelist(atom, 144)
        blacklist(atom, 140)


def opls_145(atom):
    # Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl
    if benzene(atom):
        whitelist(atom, 145)
        blacklist(atom, [141, 142])


def opls_145B(atom):
    # Biphenyl C1
    # Store for checking for neighbors outside the first ring.
    ring_one = benzene(atom)
    if ring_one:
        for neighbor in atom.neighbors:
            if neighbor not in ring_one:
                if benzene(neighbor):
                    whitelist(atom, '145B')
                    blacklist(atom, [141, 145])


def opls_146(atom):
    # Benzene H - 12 site.
    if check_atom(atom.neighbors[0], 145):
        whitelist(atom, 146)
        blacklist(atom, [140, 144])


#def opls_147
    # Napthalene fusion C (C9)


def opls_148(atom):
    # C: CH3, toluene
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            if check_atom(neighbor, 145):
                whitelist(atom, 148)
                blacklist(atom, 135)


def opls_149(atom):
    # C: CH2, ethyl benzene
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            if check_atom(neighbor, 145):
                whitelist(atom, 149)
                blacklist(atom, 136)


def opls_154(atom):
    # all-atom O: mono alcohols
    whitelist(atom, 154)


def opls_155(atom):
    # all-atom H(O): mono alcohols, OP(=O)2
    if check_atom(atom.neighbors[0], 154):
        whitelist(atom, 155)


def opls_218(atom):
    # C in CH2OH - benzyl alcohols
    benzene_carbon = False
    alcohol_oxygen = False
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            benzene_carbon = check_atom(neighbor, 145)
        if neighbor.kind == 'O':
            alcohol_oxygen = check_atom(neighbor, 154)
    if benzene_carbon and alcohol_oxygen:
        whitelist(atom, 218)


def opls_221(atom):
    # C(CH2OH)   - benzyl alcohols
    if check_atom(atom, 145):  # Already identified as part of benzene.
        for neighbor in atom.neighbors:
            if check_atom(neighbor, 218):
                whitelist(atom, 221)
                # Blacklist 3-valent carbons w/ 3 carbon neighbors.
                blacklist(atom, [141, 145, '145B'])


def opls_232(atom):
    # C: C=0 in benzaldehyde, acetophenone (CH)
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            benzene_carbon = check_atom(neighbor, 145)
        if neighbor.kind == 'O':
            aldehyde_oxygen = check_atom(neighbor, 278)
    if benzene_carbon and aldehyde_oxygen:
        whitelist(atom, 232)


def opls_278(atom):
    # AA O: aldehyde
    whitelist(atom, 278)


def opls_279(atom):
    # AA H-alpha in aldehyde & formamide
    if check_atom(atom.neighbors[0], [232, 277]):
        whitelist(atom, 279)
        blacklist(atom, [140, 144, 146])


if __name__ == "__main__":
    import pdb

    # m = Alkane(n=3)
    # m = Compound.load(get_opls_fn('isopropane.pdb'))
    # m = Compound.load(get_opls_fn('cyclohexane.pdb'))
    # m = Compound.load(get_opls_fn('neopentane.pdb'))
    # m = Compound.load(get_opls_fn('benzene.pdb'))
    m = Compound.load(get_opls_fn('1-propene.pdb'))
    # m = Compound.load(get_opls_fn('biphenyl.pdb'))

    opls_atomtypes(m)

    for atom in m.atoms:
        print "Atom kind={}, opls_whitelist={},  opls_blacklist={}".format(
            atom.kind, atom.opls_whitelist, atom.opls_blacklist)




