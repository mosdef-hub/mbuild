"""Bond graph for mBuild Compounds.

NOTE: The functions 'connected_components' and '_bfs' have been obtained (with
modifications) from the Networkx Python package, which is distributed under the
following BSD license:

Copyright (C) 2004-2016, NetworkX Developers
Aric Hagberg <hagberg@lanl.gov>
Dan Schult <dschult@colgate.edu>
Pieter Swart <swart@lanl.gov>
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

 * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

 * Neither the name of the NetworkX Developers nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from collections import defaultdict

from mbuild.utils.orderedset import OrderedSet


class BondGraph(object):
    """A graph-like object used to store and manipulate bonding information.

    `BondGraph` is designed to mimic the API and partial functionality of
     NetworkX's `Graph` data structure.

    """

    def __init__(self):
        self._adj = defaultdict(OrderedSet)

    def add_node(self, node):
        """Add a node to the bond graph."""
        if not self.has_node(node):
            self._adj[node] = OrderedSet()

    def remove_node(self, node):
        """Remove a node from the bond graph."""
        adj = self._adj
        for other_node in self.nodes():
            if node in adj[other_node]:
                self.remove_edge(node, other_node)

    def has_node(self, node):
        """Determine whether the graph contains a node."""
        return node in self._adj

    def nodes(self):
        """Return all nodes of the bond graph."""
        return [node for node in self._adj]

    def nodes_iter(self):
        """Iterate through the nodes."""
        for node in self._adj:
            yield node

    def number_of_nodes(self):
        """Get the number of nodes in the graph."""
        return sum(1 for _ in self.nodes_iter())

    def add_edge(self, node1, node2):
        """Add an edge to the bond graph."""
        self._adj[node1].add(node2)
        self._adj[node2].add(node1)

    def remove_edge(self, node1, node2):
        """Remove an edge from the bond graph."""
        adj = self._adj
        if self.has_node(node1) and self.has_node(node2):
            adj[node1].remove(node2)
            adj[node2].remove(node1)
            if not adj[node1]:
                del adj[node1]
            if not adj[node2]:
                del adj[node2]
        else:
            raise ValueError(
                "There is no edge between {} and {}".format(node1, node2)
            )

    def has_edge(self, node1, node2):
        """Determine whether the graph contains an edge."""
        if self.has_node(node1):
            return node2 in self._adj[node1]

    def edges(self):
        """Return all edges in the bond graph."""
        edges = OrderedSet()
        for node, neighbors in self._adj.items():
            for neighbor in neighbors:
                bond = (
                    (node, neighbor)
                    if self.nodes().index(node) > self.nodes().index(neighbor)
                    else (neighbor, node)
                )
                edges.add(bond)
        return list(edges)

    def edges_iter(self):
        """Iterate through the edges in the bond graph."""
        for edge in self.edges():
            yield edge

    def number_of_edges(self):
        """Get the number of edges in the graph."""
        return sum(1 for _ in self.edges())

    def neighbors(self, node):
        """Get all neighbors of the given node."""
        if self.has_node(node):
            return [neighbor for neighbor in self._adj[node]]
        else:
            return []

    def neighbors_iter(self, node):
        """Iterate through the neighbors of the given node."""
        if self.has_node(node):
            return (neighbor for neighbor in self._adj[node])
        else:
            return iter(())

    def compose(self, graph):
        """Compose this graph with the given graph."""
        adj = self._adj
        for node, neighbors in graph._adj.items():
            if self.has_node(node):
                (adj[node].add(neighbor) for neighbor in neighbors)
            else:
                # Add new node even if it has no bond/neighbor
                adj[node] = neighbors

    def subgraph(self, nodes):
        """Return a subgraph view of the subgraph induced on given nodes."""
        new_graph = BondGraph()
        nodes = list(nodes)
        adj = self._adj
        for node in nodes:
            if node not in adj:
                continue
            for neighbor in adj[node]:
                if neighbor in nodes:
                    new_graph.add_edge(node, neighbor)
        return new_graph

    def connected_components(self):
        """Generate connected components."""
        seen = set()
        components = []
        for v in self.nodes():
            if v not in seen:
                c = set(self._bfs(v))
                components.append(list(c))
                seen.update(c)

        return components

    def _bfs(self, source):
        seen = set()
        nextlevel = {source}
        while nextlevel:
            thislevel = nextlevel
            nextlevel = set()
            for v in thislevel:
                if v not in seen:
                    yield v
                    seen.add(v)
                    nextlevel.update(self.neighbors(v))
