from collections import defaultdict
from itertools import combinations_with_replacement
import sys
from warnings import warn

import matplotlib.pyplot as plt
import networkx as nx

from mbuild.orderedset import OrderedSet

# Map rule ids to the functions that check for them (see `find_atomtypes()`).
rule_number_to_rule = dict()

# Organizes the rules (see `build_rule_map()`).
rule_map = dict()

# Global neighbor information (see `neighbor_types()`).
neighbor_types_map = dict()

# Global neighbor whitelist information (see `neighbor_whitelist_types()`).
neighbor_whitelist_map = dict()

# Used for more descriptive output when sanitizing rules.
rule_number_to_doc_string = dict()


def find_atomtypes(compound, forcefield='OPLS-AA', debug=True):
    """Determine atomtypes for all atoms in `compound`. """
    if forcefield.lower() == 'opls-aa':
        import mbuild.tools.parameterize.oplsaa as oplsaa
        # Build a map to all of the supported opls_* functions.
        for func_name, func in sys.modules[oplsaa.__name__].__dict__.items():
            if func_name.startswith('opls_'):
                rule_number_to_rule[func_name.split("_")[1]] = func

    elif forcefield == 'uff':
        import mbuild.tools.parameterize.uff as uff
        # Build a map to all of the supported uff_* functions.
        for func_name, func in sys.modules[uff.__name__].__dict__.items():
            if func_name.startswith('uff_'):
                rule_number_to_rule[func_name.split("_")[1]] = func

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
    """Add white- and blacklists to each `atom` in `compound`.

    The use of ordered sets is not strictly necessary but it helps when
    debugging because it shows the order in which rules are added.
    """
    for atom in compound.yield_atoms():
        atom.extras['whitelist'] = OrderedSet()
        atom.extras['blacklist'] = OrderedSet()


def iterate_rules(compound, max_iter=10):
    """Iteratively run all the rules until the white- and backlists converge.

    Args:
        compound (Compound): The Compound being parameterized.
        max_iter (int, optional): The maximum number of iterations.
    """
    for _ in range(max_iter):
        # For comparing the lengths of the white- and blacklists.
        old_len = 0
        new_len = 0
        for atom in compound.yield_atoms():
            # We only add and never remove from the white- and blacklists so we
            # do not need to check the lengths for each individual atom.
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
    """Determine the final atomtypes from the white- and blacklists."""
    for i, atom in enumerate(compound.atoms):
        atomtype = atom.whitelist - atom.blacklist
        atomtype = [a for a in atomtype]

        if len(atomtype) == 1:
            atom.extras['atomtype'] = [atomtype[0]]
        else:
            warn("CHECK YOUR TOPOLOGY. Found multiple or no types for atom "
                 "{0} ({1}): {2}.".format(i, atom.kind, atomtype))
            atom.extras['atomtype'] = ', '.join(atomtype)


def neighbor_element_types(atom):
    """Returns the number of neighbors of each element type for an `atom`.

    The dict maintained is `neighbor_types_map` and is organized as follows:
        atom: defaultdict{element: number of neighbors of that element type}
    E.g. for an atom with 3 carbon and 1 hydrogen neighbors:
        Atom: {'C': 3, 'H': 1}

    If the queried `atom` is not already in `neighbor_types_map`, its entry will
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
    """Returns the number of neighbors of each type of whitelisted rule.

    The dict maintained is `neighbor_whitelist_map` and is organized as follows:
        atom: defaultdict{rule id: number of neighbors with rule id whitelisted}
    E.g. for a carbon atom in the middle of an alkane chain:
        Atom: {'136': 2, '140': 2}

    If the queried `atom` is not already in `neighbor_whitelist_map`, an empty
    dict is returned. This behavior differs from `neighbor_element_types()`
    because whitelisted rules are not known a priori. New entries are only
    added when the @Whitelist decorator is used.

    See Also `increment_neighbor_whitelist()`.
    """
    if atom in neighbor_whitelist_map:
        return neighbor_whitelist_map[atom]
    else:
        return dict()


def increment_neighbor_whitelist(atom, whitelist_type):
    """Increment the counts in an atom's `neighbor_whitelist_map` entry. """
    for neighbor in atom.neighbors:
        if neighbor not in neighbor_whitelist_map:
            neighbor_whitelist_map[neighbor] = defaultdict(int)
        neighbor_whitelist_map[neighbor][whitelist_type] += 1


def check_atom(atom, input_rule_ids):
    """Check if any of the rules in `input_rule_ids` are in the whitelist.

    This means that the atom was once identified as being elligible for at
    least one of these rules. This can be useful for checking, e.g., if a carbon
    was ever identified as being part of a benzene ring.
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
    """Base class for rule decorators. """
    @staticmethod
    def extract_doc_string(func):
        if func.func_doc:
            rule_number = func.func_name.split("_")[1]
            rule_number_to_doc_string[rule_number] = func.func_doc


class Element(RuleDecorator):
    """Designate the element that this rule applies to. """
    def __init__(self, element_type):
        self.element_type = element_type

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if atom.kind == self.element_type:
                return func(atom)
        return wrapped


class InWhitelist(RuleDecorator):
    """Checks if this atom already has different rule whitelisted. """
    def __init__(self, element_type):
        self.element_type = element_type

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if self.element_type in atom.whitelist:
                return func(atom)
        return wrapped


class NeighborCount(RuleDecorator):
    """Designate the number of neighbors an atom must have for this rule. """
    def __init__(self, count):
        self.count = count

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if len(atom.neighbors) == self.count:
                return func(atom)
        return wrapped


class NeighborsBase(RuleDecorator):
    """Super class for rule decorators that count types of neighbors. """
    def __init__(self, neighbor_type, count):
        neighbor_type = str(neighbor_type)
        self.neighbor_type = neighbor_type
        self.count = count

    def match_count(self, atom):
        if self.neighbor_type in neighbor_element_types(atom):
            match_count = neighbor_element_types(atom)[self.neighbor_type]
        elif self.neighbor_type in neighbor_whitelist_types(atom):
            match_count = neighbor_whitelist_types(atom)[self.neighbor_type]
        else:
            match_count = 0
        return match_count


class NeighborsExactly(NeighborsBase):
    """Designate the rule's exact number of neighbors of a specific type.

    The "specific type" can either be an element or a whitelisted rule.
    """
    def __init__(self, neighbor_type, count):
        super(NeighborsExactly, self).__init__(neighbor_type, count)

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if self.match_count(atom) == self.count:
                return func(atom)
        return wrapped


class NeighborsAtLeast(NeighborsBase):
    """Designate the rule's minimum number of neighbors of a specific type.

    The "specific type" can either be an element or a whitelisted rule.
    """
    def __init__(self, neighbor_type, count):
        super(NeighborsAtLeast, self).__init__(neighbor_type, count)

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if self.match_count(atom) >= self.count:
                return func(atom)
        return wrapped


class NeighborsAtMost(NeighborsBase):
    """Designate the rule's maximum number of neighbors of a specific type.

    The "specific type" can either be an element or a whitelisted rule.
    """
    def __init__(self, neighbor_type, count):
        super(NeighborsAtMost, self).__init__(neighbor_type, count)

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if self.match_count(atom) <= self.count:
                return func(atom)
        return wrapped


class Whitelist(RuleDecorator):
    """Whitelist an atomtype for an atom. """
    def __init__(self, rule_number):
        if isinstance(rule_number, (list, tuple, set)):
            warn("Rules should only whitelist themselves.")
        else:
            self.rule_number = str(rule_number)

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if func(atom):
                self.whitelist(atom)
                return True
        return wrapped

    def whitelist(self, atom):
        """Add the rule to the atom's whitelist. """
        atom.whitelist.add(str(self.rule_number))
        increment_neighbor_whitelist(atom, self.rule_number)


class Blacklist(RuleDecorator):
    """Blacklist an atomtype for an atom. """
    def __init__(self, rule_numbers):
        if isinstance(rule_numbers, (list, tuple, set)):
            self.rule_numbers = list(map(str, rule_numbers))
            self.rule_numbers.sort()
        else:
            self.rule_numbers = [str(rule_numbers)]

    def __call__(self, func):
        self.extract_doc_string(func)

        def wrapped(atom):  # this must be called 'wrapped'
            if func(atom):
                self.blacklist(atom)
                return True
        return wrapped

    def blacklist(self, atom):
        """Add the rule to the atom's blacklist. """
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
    supported_elements = find_all_supported_elements()
    rule_matches = find_all_rule_matches(supported_elements)

    # Build directed graphs showing which rules blacklist each other.
    for key, rules in rule_matches.items():
        # Only consider patterns matched by multiple rules.
        if len(rules) < 2:
            continue

        element_type, pattern = key
        graph = nx.DiGraph()
        for rule_number in rules:
            graph.add_node(rule_number)
            blacklisted_rules = set()
            decorators = get_decorator_objects_by_type(
                rule_number_to_rule[rule_number], RuleDecorator)

            for dec in decorators:
                if isinstance(dec, Blacklist):
                    blacklisted_rules.update(dec.rule_numbers)
            for blacklisted_rule in blacklisted_rules:
                graph.add_edge(rule_number, blacklisted_rule)

        if not nx.is_connected(graph.to_undirected()):
            draw_rule_graph('unconnected', graph, element_type, pattern)

        if not nx.is_directed_acyclic_graph(graph):
            draw_rule_graph('not_DAG', graph, element_type, pattern)

        # Check for multiple sinks. This is not necessarily incorrect.
        sinks = []
        for node in graph.nodes():
            if len(nx.descendants(graph, node)) == 0:
                sinks.append(node)
        if len(sinks) > 1:
            # Cases with multiple sinks that we can safely ignore.
            # TODO: Provide more robust way to add approved exceptions.
            if '141' in sinks and '142' in sinks and '145' in graph.nodes():
                # This case occurs because 145 requires AT LEAST 2 carbon
                # neighbors (3 total) and thus can match multiple patterns.
                continue
            draw_rule_graph('multiple_sinks', graph, element_type, pattern,
                            sinks=sinks)

        # Check for multiple sources. This is not necessarily incorrect.
        sources = []
        for node in graph.nodes():
            if len(nx.ancestors(graph, node)) == 0:
                sources.append(node)
        if len(sources) > 1:
            draw_rule_graph('multiple_sources', graph, element_type, pattern,
                            sources=sources)


def find_all_supported_elements():
    """Find all elements currently supported by rules.

    Returns:
        supported_elements (list): A sorted list of all the supported elements.
    """
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
    return supported_elements


def find_all_rule_matches(supported_elements):
    """Find all elements and combinations of neighbor types that have a rule.

    Returns:
        rule_matches (dict): All patterns and rules that apply. The dictionary
            is structured as follows:
                key: (element, (neighbor element 1, neighbor element 2, etc..))
                value: set(rule numbers)
            Example entry (from time of writing this comment):
                ('C', ('C', 'C', 'H')): set(['145', '142'])
    """
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
        all_patterns = set(combinations_with_replacement(supported_elements,
                                                         neighbor_count))

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
                rule_matches[(element_type, pattern)] = {rule_number}
            else:
                rule_matches[(element_type, pattern)].add(rule_number)
    return rule_matches


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
    # Find a function called `wrapper` in the function's closure, and recurse.
    for cell in decorated_function.func_closure:
        closure_entry = cell.cell_contents
        if hasattr(closure_entry, '__name__') and closure_entry.__name__ is "wrapped":
            wrapped_decorator_objects = get_decorator_objects_by_type(
                closure_entry, decorator_type)
            decorators += wrapped_decorator_objects
            break
    return decorators


def check_duplicate_element(element_type, rule_number):
    """Used to catch rules with multiple @Element decorators.

    Args:
        element_type (str): The element type in the decorator. This is only
            used to check if it has been assigned previously.
        rule_number (int): The decorated rule in question.
    """
    assert element_type is None, ("Duplicate element type decorators on rule "
                                  "{}".format(rule_number))


def check_duplicate_neighbor_count(neighbor_count, rule_number):
    """Used to catch rules with multiple @NeighborCount decorators.

    Args:
        neighbor_count (int): The neighbor count in the decorator. This is only
            used to check if it has been assigned previously.
        rule_number (int): The decorated rule in question.
    """
    assert neighbor_count is None, ("Duplicate neighbor count decorators on "
                                    "rule {}".format(rule_number))


def draw_rule_graph(issue, graph, element, pattern, sinks=None, sources=None):
    """Visualize the blacklist graph for a pattern with a logical inconsitency.

    Args:
        issue (str): Descriptive string of the logical inconsistency. Used for
            labeling purposes.
        graph (nx.DiGraph): The blacklist graph.
        element (str): The type of element that the conflicting rules describe.
        pattern (tuple of str): The immediate neighbors of the element in
            question.
        sinks (list of nodes, optional): Nodes that are sinks in the graph.
        sources (list of nodes, optional): Nodes that are sources in the graph.
    """
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
    else:
        warn("Can't draw rule graph. Unrecognized issue: {}".format(issue))
        return

    # TODO: color unconnected nodes
    colors = []
    for node in graph.nodes_iter():
        if sinks and node in sinks:
            colors.append('blue')
        elif sources and node in sources:
            colors.append('blue')
        else:
            colors.append('red')

    labels = {rule: '{}\n{}'.format(rule, doc) for rule, doc in
              rule_number_to_doc_string.iteritems() if rule in graph.nodes()}

    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=1000, node_color=colors)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos, labels=labels)

    # TODO: Prettier way to show the doc strings on the blacklist graph.
    plt.axis('off')
    xmin, xmax = plt.xlim()
    plt.xlim(xmin*1.5, xmax*1.5)

    fig_name = '{}-element_{}-pattern_{}.png'.format(
        issue, element, ''.join(pattern))
    plt.tight_layout()
    plt.savefig(fig_name)
    plt.clf()

    warn("{} connected to {} {}. See '{}'".format(
        element, pattern, phrase, fig_name))
