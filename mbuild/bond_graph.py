from collections import defaultdict


class BondGraph(object):
    def __init__(self):
        self._data = defaultdict(set)

    def add_node(self, node):
        if not self.has_node(node):
            self._data[node] = set()

    def remove_node(self, node):
        del self._data[node]

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
        edges = set()
        for node, neighbors in self._data.items():
            for neighbor in neighbors:
                if id(node) < id(neighbor):
                    edges.add((node, neighbor))
                else:
                    edges.add((neighbor, node))
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
        components = []
        for node, neighbors in self._data.items():
            for component in components:
                if node in component:  # Is the node already present?
                    current_component = component
                    break
            else:
                current_component = set()
                components.append(current_component)

            current_component.add(node)
            current_component.union(neighbors)

        return [list(component) for component in components]


