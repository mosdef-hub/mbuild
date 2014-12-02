from collections import defaultdict
from itertools import combinations_with_replacement
import sys
from warnings import warn

import matplotlib.pyplot as plt
import networkx as nx

from mbuild.orderedset import OrderedSet

# Map ids to the functions that check for them.
rule_number_to_rule = dict()
rule_map = dict()

# Globally maintained neighbor information (see `neighbor_types()`).
neighbor_types_map = {}

#
neighbor_whitelist_map = {}


def find_atomtypes(compound, forcefield='OPLS-AA', debug=True):
    """Determine atomtypes for all atoms in `compound`. """
    if forcefield == 'OPLS-AA':
        import oplsaa
        # Build a map to all of the supported opls_* functions.
        for fn, fcn in sys.modules[oplsaa.__name__].__dict__.items():
            if fn.startswith('opls_'):
                rule_number_to_rule[fn.split("_")[1]] = fcn

    elif forcefield == 'UFF':
        import uff
        # Build a map to all of the supported opls_* functions.
        for fn, fcn in sys.modules[uff.__name__].__dict__.items():
            if fn.startswith('uff_'):
                rule_number_to_rule[fn.split("_")[1]] = fcn

    build_rule_map()
    if debug:
        sanitize()
    prepare_compound(compound)
    iterate_rules(compound, max_iter=10)
    resolve_atomtypes(compound)


def build_rule_map():
    """Build up a tree of element types-->neighbor counts-->rules. """
    for rule_number, rule in rule_number_to_rule.items():
        decorators = get_decorator_objects_by_type(rule, RuleDecorator)
        element_type = None
        neighbor_count = None
        for dec in decorators:
            if isinstance(dec, Element):
                element_type = dec.element_type
            if isinstance(dec, NeighborCount):
                neighbor_count = dec.count

        if not element_type:
            warn('Rule {} has no element type'.format(rule_number))
        if not neighbor_count:
            warn('Rule {} has no neighbor count'.format(rule_number))

        if element_type not in rule_map:
            rule_map[element_type] = dict()
        if neighbor_count not in rule_map[element_type]:
            rule_map[element_type][neighbor_count] = []
        rule_map[element_type][neighbor_count].append(rule_number)


def prepare_compound(compound):
    """Add white- and blacklists to each `atom` in `compound`. """
    for atom in compound.yield_atoms():
        atom.extras['whitelist'] = OrderedSet()
        atom.extras['blacklist'] = OrderedSet()


def iterate_rules(compound, max_iter=10):
    """Iteratively run all the rules until the white- and backlists converge.

    NOTE: We only add and never remove from the white- and blacklists so we
          don't need to check the lengths for each individual atom.

    Args:
        compound:
        max_iter:
    """
    for iter_cnt in range(max_iter):
        # For comparing the lengths of the white- and blacklists.
        old_len = 0
        new_len = 0
        for atom in compound.yield_atoms():
            old_len += len(atom.whitelist)
            old_len += len(atom.blacklist)

            if atom.kind == 'G':  # Ignore Ports.
                continue

            # Only run rules with matching element and neighbor counts.
            if atom.kind in rule_map:
                if len(atom.neighbors) in rule_map[atom.kind]:
                    for rule in rule_map[atom.kind][len(atom.neighbors)]:
                        run_rule(atom, rule)
                else:
                    warn("No rule for {}-neighbor '{}' atom".format(
                        len(atom.neighbors), atom.kind))
            else:
                warn("No rule for atom kind '{}'".format(atom.kind))

            new_len += len(atom.whitelist)
            new_len += len(atom.blacklist)

        # Nothing changed, we're done!
        if old_len == new_len:
            break
    else:
        warn("Reached maximum iterations. Something probably went wrong.")


def run_rule(atom, rule_id):
    """Execute the rule function for a specified atomtype. """
    if rule_id not in atom.whitelist:
        try:
            rule_fn = rule_number_to_rule[str(rule_id)]
        except KeyError:
            raise KeyError('Rule for {} not implemented'.format(rule_id))
        rule_fn(atom)


def resolve_atomtypes(compound):
    """Determine the final atomtype from the white- and blacklists."""
    for i, atom in enumerate(compound.atoms):
        atomtype = atom.whitelist - atom.blacklist
        atomtype = [a for a in atomtype]

        if len(atomtype) == 1:
            atom.extras['atomtype'] = [atomtype[0]]
        else:
            warn("CHECK YOUR TOPOLOGY. Found multiple or no types for atom {0} ({1}): {2}.".format(
                    i, atom.kind, atomtype))
            atom.extras['atomtype'] = ', '.join(atomtype)


def neighbor_types(atom):
    """Returns the number of neighbors of each element type for an `atom`.

    The dict maintained is `neighbor_types_map` and is organized as follows:
        atom: defaultdict{element: number of neighbors of that element type}
    E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
        Atom: {'C': 3, 'H': 1}

    If the queried `atom` is not already in `neighbor_types_map`, it entry will
    be added.
    """
    if atom not in neighbor_types_map:
        neighbors = defaultdict(int)
        for neighbor in atom.neighbors:
            kind = neighbor.kind
            neighbors[kind] += 1
        neighbor_types_map[atom] = neighbors
    return neighbor_types_map[atom]


def neighbor_whitelist_types(atom):
    """Returns the number of neighbors of each element type for an `atom`.

    The dict maintained is `neighbor_types_map` and is organized as follows:
        atom: defaultdict{element: number of neighbors of that element type}
    E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
        Atom: {'C': 3, 'H': 1}

    If the queried `atom` is not already in `neighbor_types_map`, it entry will
    be added.
    """
    if atom in neighbor_whitelist_map:
        return neighbor_whitelist_map[atom]
    else:
        return dict()


def append_neighbor_whitelist(atom, whitelist_type):
    """Returns the number of neighbors of each element type for an `atom`.

    The dict maintained is `neighbor_types_map` and is organized as follows:
        atom: defaultdict{element: number of neighbors of that element type}
    E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
        Atom: {'C': 3, 'H': 1}

    If the queried `atom` is not already in `neighbor_types_map`, it entry will
    be added.
    """
    for neighbor in atom.neighbors:
        if neighbor not in neighbor_whitelist_map:
            neighbor_whitelist_map[neighbor] = defaultdict(int)
        neighbor_whitelist_map[neighbor][whitelist_type] += 1


def check_atom(atom, input_rule_ids):
    """Check if any of the rules in `input_rule_ids` are in the whitelist.

    This means that the atom was once identified as being elligible for at least
    one of these rules. This can be useful for checking, e.g. if a carbon was
    ever identified as being part of a benzene ring.
    """
    rule_ids = set()
    if isinstance(input_rule_ids, (list, tuple, set)):
        for rule in input_rule_ids:
            rule_ids.add(str(rule))
    else:
        rule_ids.add(str(input_rule_ids))

    for rule in rule_ids:
        if rule in atom.whitelist:
            return True

# -------------------- #
# Decorators for rules #
# -------------------- #


class RuleDecorator(object):
    pass


class Element(RuleDecorator):
    def __init__(self, element_type):
        self.element_type = element_type

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if atom.kind == self.element_type:
                return f(atom)

        return wrapped


class InWhitelist(RuleDecorator):
    def __init__(self, element_type):
        self.element_type = element_type

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if self.element_type in atom.whitelist:
                return f(atom)

        return wrapped


class NeighborCount(RuleDecorator):
    def __init__(self, count):
        self.count = count

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if len(atom.neighbors) == self.count:
                return f(atom)
        return wrapped


class NeighborsBase(RuleDecorator):
    def __init__(self, neighbor_type, count):
        neighbor_type = str(neighbor_type)
        self.neighbor_type = neighbor_type
        self.count = count

    def match_count(self, atom):
        if self.neighbor_type in neighbor_types(atom):
            match_count = neighbor_types(atom)[self.neighbor_type]
        elif self.neighbor_type in neighbor_whitelist_types(atom):
            match_count = neighbor_whitelist_types(atom)[self.neighbor_type]
        else:
            match_count = 0
        return match_count


class NeighborsExactly(NeighborsBase):
    def __init__(self, neighbor_type, count):
        super(NeighborsExactly, self).__init__(neighbor_type, count)

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if self.match_count(atom) == self.count:
                return f(atom)
        return wrapped


class NeighborsAtLeast(NeighborsBase):
    def __init__(self, neighbor_type, count):
        super(NeighborsAtLeast, self).__init__(neighbor_type, count)

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if self.match_count(atom) >= self.count:
                return f(atom)
        return wrapped


class NeighborsAtMost(NeighborsBase):
    def __init__(self, neighbor_type, count):
        super(NeighborsAtMost, self).__init__(neighbor_type, count)

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if self.match_count(atom) <= self.count:
                return f(atom)
        return wrapped


class Whitelist(RuleDecorator):
    def __init__(self, rule_number):
        if isinstance(rule_number, (list, tuple, set)):
            warn("Rules should only whitelist themselves.")
        else:
            self.rule_number = str(rule_number)

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if f(atom):
                self.whitelist(atom)
                return True
        return wrapped

    def whitelist(self, atom):
        """Whitelist an OPLS-aa atomtype for an atom. """
        atom.whitelist.add(str(self.rule_number))
        append_neighbor_whitelist(atom, self.rule_number)


class Blacklist(RuleDecorator):
    def __init__(self, rule_numbers):
        if isinstance(rule_numbers, (list, tuple, set)):
            self.rule_numbers = list(map(str, rule_numbers))
            self.rule_numbers.sort()
        else:
            self.rule_numbers = [str(rule_numbers)]

    def __call__(self, f):
        # this must be called 'wrapped'
        def wrapped(atom):
            if f(atom):
                self.blacklist(atom)
                return True
        return wrapped

    def blacklist(self, atom):
        """Blacklist an OPLS-aa atomtype for an atom. """
        if isinstance(self.rule_numbers, (list, tuple, set)):
            for rule in self.rule_numbers:
                atom.blacklist.add(str(rule))
        else:
            atom.blacklist.add(str(self.rule_numbers))

# ------------------------------------- #
# Sanitization and associated functions #
# ------------------------------------- #


def sanitize():
    """Analyze all rules for possible inconsistencies.

    This function serves primarily as a tool for developers who intend to add
    new rules or modify existing ones. Ideally, it will help you identify and
    correct logical inconsistencies as early as possible. Additionally, it
    suggests other rules that you may want to consider blacklisting.
    """

    # Find all elements currently supported by rules.
    supported_elements = set()
    for rule_number, rule in rule_number_to_rule.items():
        decorators = get_decorator_objects_by_type(rule, RuleDecorator)

        element_type = None
        for dec in decorators:
            if isinstance(dec, Element):
                check_duplicate_element(element_type, rule_number)
                element_type = dec.element_type
                supported_elements.add(element_type)
    supported_elements = list(supported_elements)
    supported_elements.sort()

    # Find all elements and combinations of neighbor types that have a rule.
    # Rule matches is structured as follows:
    #   key: (element, (neighbor element 1, neighbor element 2, etc..))
    #   value: set(rule numbers)
    # Example entry (from time of writing this comment):
    #   ('C', ('C', 'C', 'H')): set(['145', '142'])
    rule_matches = dict()
    for rule_number, rule in rule_number_to_rule.items():
        decorators = get_decorator_objects_by_type(rule, RuleDecorator)

        element_type = None
        neighbor_count = None
        for dec in decorators:
            if isinstance(dec, Element):
                element_type = dec.element_type
            if isinstance(dec, NeighborCount):
                check_duplicate_neighbor_count(neighbor_count, rule_number)
                neighbor_count = dec.count
        # All POSSIBLE combinations of elements and neighbors.
        all_patterns = set(combinations_with_replacement(supported_elements, neighbor_count))

        # Remove the ones that don't actually have a rule.
        removed_patterns = set()
        for dec in decorators:
            if isinstance(dec, NeighborsExactly):
                for pattern in all_patterns:
                    if not pattern.count(dec.neighbor_type) == dec.count:
                        removed_patterns.add(pattern)
            elif isinstance(dec, NeighborsAtLeast):
                for pattern in all_patterns:
                    if not pattern.count(dec.neighbor_type) >= dec.count:
                        removed_patterns.add(pattern)
            elif isinstance(dec, NeighborsAtMost):
                for pattern in all_patterns:
                    if not pattern.count(dec.neighbor_type) <= dec.count:
                        removed_patterns.add(pattern)
        all_patterns.difference_update(removed_patterns)

        for pattern in all_patterns:
            if (element_type, pattern) not in rule_matches:
                rule_matches[(element_type, pattern)] = set([rule_number])
            else:
                rule_matches[(element_type, pattern)].add(rule_number)

    # Build directed graphs showing which rules blacklist each other.
    for key, rules in rule_matches.items():
        # Only consider patterns matched by multiple rules.
        if len(rules) < 2:
            continue

        element_type, pattern = key
        G = nx.DiGraph()
        for rule_number in rules:
            G.add_node(rule_number)
            blacklisted_rules = set()
            decorators = get_decorator_objects_by_type(rule_number_to_rule[rule_number], RuleDecorator)

            for dec in decorators:
                if isinstance(dec, Blacklist):
                    blacklisted_rules.update(dec.rule_numbers)
            for blacklisted_rule in blacklisted_rules:
                G.add_edge(rule_number, blacklisted_rule)

        # Check if graph is connected.
        if not nx.is_connected(G.to_undirected()):
            draw_rule_graph('unconnected', G, element_type, pattern)

        # Check if DAG.
        if not nx.is_directed_acyclic_graph(G):
            draw_rule_graph('not_DAG', G, element_type, pattern)

        # Check if there are multiple sinks. This is not necessarily incorrect.
        sinks = []
        for node in G.nodes():
            if len(nx.descendants(G, node)) == 0:
                sinks.append(node)
        if len(sinks) > 1:
            # Cases with multiple sinks that we can safely ignore.
            # TODO: Provide more robust way to add approved exceptions.
            if '141' in sinks and '142' in sinks and '145' in G.nodes():
                # This case occurs because 145 requires AT LEAST 2 carbon
                # neighbors (3 total) and thus can match multiple patterns.
                continue
            draw_rule_graph('multiple_sinks', G, element_type, pattern,
                            sinks=sinks)

        # Check if there are multiple sources. This is not necessarily incorrect.
        sources = []
        for node in G.nodes():
            if len(nx.ancestors(G, node)) == 0:
                sources.append(node)
        if len(sources) > 1:
            draw_rule_graph('multiple_sources', G, element_type, pattern,
                            sources=sources)


def get_decorator_objects_by_type(decorated_function, decorator_type):
    """Find all decorators of a particular type on a function.

    Args:
        decorated_function:
        decorator_type:

    Returns:
        rval:
    """
    decorators = []
    # Find an object of decorator_type in the function's closure. There should
    # be only one.
    for cell in decorated_function.func_closure:
        closure_entry = cell.cell_contents
        if isinstance(closure_entry, decorator_type):
            decorators.append(closure_entry)
            break
    # Find a function called `wrapper` in the function's closure, and recurse on that.
    for cell in decorated_function.func_closure:
        closure_entry = cell.cell_contents
        if hasattr(closure_entry, '__name__') and closure_entry.__name__ is "wrapped":
            wrapped_decorator_objects = get_decorator_objects_by_type(closure_entry, decorator_type)
            decorators += wrapped_decorator_objects
            break
    return decorators


def check_duplicate_element(element_type, rule_number):
    assert element_type is None, ("Duplicate element type decorators on rule "
                                  "{}".format(rule_number))


def check_duplicate_neighbor_count(neighbor_count, rule_number):
    assert neighbor_count is None, ("Duplicate neighbor count decorators on "
                                    "rule {}".format(rule_number))


def draw_rule_graph(issue, G, element, pattern, sinks=None, sources=None):
    """
    Args:
        issue:
        G:
        element:
        pattern:
        sinks:
    """

    # TODO: color unconnected nodes
    colors = []
    for node in G.nodes_iter():
        if sinks and node in sinks:
            colors.append('blue')
        elif sources and node in sources:
            colors.append('blue')
        else:
            colors.append('red')

    nx.draw_circular(G, node_size=1000, node_color=colors)
    fig_name = '{}-element_{}-pattern_{}.png'.format(
        issue, element, ''.join(pattern))
    plt.savefig(fig_name)
    plt.clf()

    if issue == 'unconnected':
        phrase = 'is not connected'
    elif issue == 'not_DAG':
        phrase = 'is not a DAG'
    elif issue == 'multiple_sinks':
        assert sinks is not None
        phrase = 'has multiple sinks: {}'.format(sinks)
    elif issue == 'multiple_sources':
        assert sources is not None
        phrase = 'has multiple sources: {}'.format(sources)

    warn("{} connected to {} {}. See '{}'".format(
        element, pattern, phrase, fig_name))
