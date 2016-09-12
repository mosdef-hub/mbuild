from collections import OrderedDict, defaultdict
from oset import oset

class OrderedDefaultDict(OrderedDict, defaultdict):
    def __init__(self, default_factory=None, *args, **kwargs):
        #in python3 you can omit the args to super
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)
        self.default_factory = default_factory

class BondGraph(object):
    """A graph-like object used to store and manipulate bonding information.

    `BondGraph` is designed to mimic the API and partial functionality of
     NetworkX's `Graph` data structure.

    """
    def __init__(self):
        self._data = OrderedDefaultDict(oset)

    def add_node(self, node):
        if not self.has_node(node):
            self._data[node] = oset()

    def remove_node(self, node):
        for other_node in self.nodes():
            if node in self._data[other_node]:
                self.remove_edge(node, other_node)

    def has_node(self, node):
        return node in self._data

    def nodes(self):
        return [node for node in self._data]

    def nodes_iter(self):
        for node in self._data:
            yield node

    def number_of_nodes(self):
        return sum(1 for _ in self.nodes_iter())

    def add_edge(self, node1, node2):
        self._data[node1].add(node2)
        self._data[node2].add(node1)

    def remove_edge(self, node1, node2):
        if self.has_node(node1) and self.has_node(node2):
            self._data[node1].remove(node2)
            self._data[node2].remove(node1)
            if not self._data[node1]:
                del self._data[node1]
            if not self._data[node2]:
                del self._data[node2]
        else:
            raise ValueError('There is no edge between {} and {}'.format(
                node1, node2))

    def has_edge(self, node1, node2):
        if self.has_node(node1):
            return node2 in self._data[node1]

    def edges(self):
        edges = oset()
        for node, neighbors in self._data.items():
            for neighbor in neighbors:
                bond = (node, neighbor) if id(node) < id(neighbor) else (neighbor, node)
                edges.add(bond)
        return list(edges)

    def edges_iter(self):
        for edge in self.edges():
            yield edge

    def number_of_edges(self):
        return sum(1 for _ in self.edges())

    def neighbors(self, node):
        if self.has_node(node):
            return [neighbor for neighbor in self._data[node]]
        else:
            return []

    def neighbors_iter(self, node):
        if self.has_node(node):
            return (neighbor for neighbor in self._data[node])
        else:
            return iter(())

    def compose(self, graph):
        for node, neighbors in graph._data.items():
            if self.has_node(node):
                self._data[node].union(neighbors)
            elif neighbors:
                self._data[node] = neighbors

    def subgraph(self, nodes):
        new_graph = BondGraph()
        nodes = list(nodes)
        for node in nodes:
            if node not in self._data:
                continue
            for neighbor in self._data[node]:
                if neighbor in nodes:
                    new_graph.add_edge(node, neighbor)
        return new_graph

    def connected_components(self):
        """ """
        def go_deeper(current_component, n):
            """ """
            current_component.add(n)
            neighbors = self._data[n]
            for neighbor in neighbors:
                if neighbor not in current_component:
                    go_deeper(current_component, neighbor)

        components = []
        for node in self._data:
            # Is the node in another component?
            for component in components:
                if node in component:
                    current_component = component
                    break
            else:  # We're in a new component.
                current_component = oset()
                components.append(current_component)

            go_deeper(current_component, node)

        return [list(component) for component in components]


