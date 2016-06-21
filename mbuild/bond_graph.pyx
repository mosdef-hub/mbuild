import numpy as np
cimport numpy as np

very_long = np.int64 # runtime type
ctypedef np.int64_t very_long_ct # compile-time type
ctypedef Py_ssize_t index_size

cdef class BondGraph:
    cdef index_size number_of_nodes
    cdef index_size max_number_of_nodes
    cdef index_size max_adjacency_list_length

    def __init__(self):
        """An alternative to networkx graphs.
        """
        self.number_of_nodes = 0
        self.max_number_of_nodes = 20
        self.max_adjacency_list_length = 10
        
        nodes = np.zeros(
                shape=(self.max_number_of_nodes, self.max_adjacency_list_length),
                dtype=very_long)
        node_id_list = np.zeros(shape=(self.max_number_of_nodes), dtype=very_long)
        node_list = np.zeros(shape=(self.max_number_of_nodes), dtype=object)

        cdef very_long_ct [:,:] nodes_view = self.nodes
        cdef very_long_ct [:] node_id_list_view = self.node_id_list
        cdef very_long_ct [:] node_list_view = self.node_list

    cdef _grow_nodes(self, index_size new_size):
        """Grow the nodes array and the maximum number of nodes to new_size, or
        by a factor of two, if new_size is 0. If new_size is less than
        max_number_of_nodes, the nodes array will be shrunk rather than grown.
        Note: new_size must be greater than number_of_nodes, otherwise data will
        be lost.
        """
        cdef index_size old_max_size = self.max_number_of_nodes
        cdef np.ndarray[very_long_ct, ndim=2] empty_id_array
        cdef np.ndarray[very_long_ct, ndim=1] empty_id_list
        cdef np.ndarray empty_node_list

        if new_size != 0:
            self.max_number_of_nodes = new_size
        else:
            self.max_number_of_nodes *= 2
            self.max_number_of_nodes += 1

        # growing
        if old_max_size < self.max_number_of_nodes:
            # growing nodes
            empty_id_array = np.zeros(
                    shape=(self.max_number_of_nodes-old_max_size, self.max_adjacency_list_length),
                    dtype=very_long)
            self.nodes = np.concatenate((self.nodes, empty_id_array), axis=0)
            # growing node_id_list
            empty_id_list = np.zeros(
                    shape=(self.max_number_of_nodes-old_max_size),
                    dtype=very_long)
            self.node_id_list = np.concatenate(
                    (self.node_id_list, empty_id_list),
                    axis=0)
            # growing node_list
            empty_node_list = np.zeros(
                    shape=(self.max_number_of_nodes-old_max_size),
                    dtype=very_long)
            self.node_list = np.concatenate(
                    (self.node_list, empty_node_list),
                    axis=0)
        # shrinking array
        elif old_max_size > self.max_number_of_nodes:
            self.nodes = self.nodes[0:self.max_number_of_nodes, :]
            self.node_id_list = self.node_id_list[0:self.max_number_of_nodes]
            self.node_list = self.node_list[0:self.max_number_of_nodes]

    cdef _grow_adjacency_list(self, index_size new_size):
        """Grow the nodes array and the maximum adjacency list length to
        new_size, or by a factor of two, if new_size is 0. new_size must be
        greater than max_adjacency_list_length - attempts to shrink the array
        will do nothing.
        """
        cdef index_size old_max_size = self.max_adjacency_list_length
        cdef np.ndarray[very_long_ct, ndim=2] empty_array

        if new_size != 0:
            self.max_adjacency_list_length = new_size
        else:
            self.max_adjacency_list_length *= 2
            self.max_adjacency_list_length += 1

        # growing array
        if old_max_size < self.max_adjacency_list_length:
            empty_array = np.zeros(
                    shape=(self.max_number_of_nodes, self.max_adjacency_list_length-old_max_size),
                    dtype=very_long)
            self.nodes = np.concatenate((self.nodes, empty_array), axis=1)
        else:
            self.max_adjacency_list_length = old_max_size

    cdef _get_node(self, very_long_ct node_id):
        """Given a node's id (as is stored in self.nodes), finds and returns the
        object that id points to, the object BondGraph was asked to store.
        """
        cdef index_size index = self.node_id_list.index(node_id)
        return self.node_list[index]

    def has_node(self, check_node):
        """Check if check_node is in BondGraph.
        Return true if check_node is in BondGraph; false otherwise.
        """
        if id(check_node) in self.nodes[0:self.number_of_nodes,0]:
            return True
        return False

        # also can check if check_node is in self.node_list TODO: check which is faster

    def _find_node(self, check_node):
        """Finds a node in BondGraph, and returns its row's index.
        If node doesn't exist, returns -1.
        """
        if not (isinstance(check_node, np.int64) or isinstance(check_node, int)):
            check_node = id(check_node)

        index, = np.where(self.nodes[0:self.number_of_nodes,0]==check_node)
        if len(index) == 0:
            return -1
        return index[0]

    def add_node(self, new_node):
        """Add new_node to BondGraph.
        If new_node is already in BondGraph, does nothing.
        """
        if self.has_node(new_node):
            return
        self._add_node(new_node)

    def _add_node(self, new_node):
        """Add new_node to BondGraph.
        Assumes new_node is not already in BondGraph.
        """
        if self.number_of_nodes == self.max_number_of_nodes:
            self._grow_nodes(0)

        self.nodes[self.number_of_nodes,0] = id(new_node)
        self.number_of_nodes += 1

        self.node_list[self.number_of_nodes] = new_node
        self.node_id_list[self.number_of_nodes] = id(new_node)

    def has_edge(self, node1, node2):
        """Check if edge between node1 and node2 is in BondGraph.
        Return true if edge exists; false otherwise.
        """
        index = self._find_node(node1)
        if index == -1:
            return False
        if id(node2) in self.nodes[index,:]:
            return True
        return False

    def add_edge(self, node1, node2):
        """Adds edge between node1 and node2
        If node1 or node2 aren't in BondGraph, they are created.
        If an edge between node1 and node2 already exists, does nothing.
        """
        if self.has_edge(node1, node2):
            return

        index1 = self._find_node(node1)
        if index1 == -1:
            self._add_node(node1)
            index1 = self.number_of_nodes-1
        index2 = self._find_node(node2)
        if index2 == -1:
            self._add_node(node2)
            index2 = self.number_of_nodes-1

        for i in range(1, self.max_adjacency_list_length):
            if self.nodes[index1,i] == 0:
                self.nodes[index1,i] = id(node2)
                break
            elif i == self.max_adjacency_list_length-1:
                self._grow_adjacency_list(0)
                self.nodes[index1,i+1] = id(node2)

        for i in range(1, self.max_adjacency_list_length):
            if self.nodes[index2,i] == 0:
                self.nodes[index2,i] = id(node1)
                break
            elif i == self.max_adjacency_list_length-1:
                self._grow_adjacency_list(0)
                self.nodes[index2,i+1] = id(node1)

    def remove_edge(self, node1, node2):
        """Remove edge between node1 and node2 from BondGraph.
        If edge does not exist, does nothing.
        """
        if not self.has_edge(node1, node2):
            return
        index1 = self._find_node(node1)
        index2 = self._find_node(node2)

        replace = False
        for i in range(1, self.max_adjacency_list_length):
            if replace:
                self.nodes[index1,i-1] = self.nodes[index1,i]
            if self.nodes[index1,i] == 0:
                break
            if not replace and id(node2) == self.nodes[index1,i]:
                replace = True
                self.nodes[index1,i] = 0

        replace = False
        for i in range(1, self.max_adjacency_list_length):
            if replace:
                self.nodes[index2,i-1] = self.nodes[index2,i]
            if self.nodes[index2,i] == 0:
                break
            if not replace and id(node1) == self.nodes[index2,i]:
                replace = True
                self.nodes[index2,i] = 0

    def remove_node(self, del_node):
        """Remove del_node from BondGraph.
        If del_node is not in BondGraph, does nothing.
        If del_node is a part of any edges, those edges will be removed.
        """
        index = self._find_node(del_node)
        if index == -1:
            return

        for i in range(self.max_adjacency_list_length-1, 0, -1): # Iterating backwards to avoid cost of back-shifting every time
            if self.nodes[index,i] != 0:
                self.remove_edge(del_node, self._get_node(self.nodes[index,i])) #Make a _remove_edge() method akin to _add_node()

        for i in range(index+1, self.number_of_nodes):
            self.nodes[i-1] = self.nodes[i]
        self.number_of_nodes -= 1

        index = self.node_id_list.index(id(del_node))
        self.node_id_list = np.delete(self.node_id_list, index)
        self.node_list = np.delete(self.node_list, index)

    def compose(self, graph2):
        """Append a bond graph to the end of this one.
        """
        if self.max_number_of_nodes != self.number_of_nodes:
            self._grow_nodes(self.number_of_nodes)

        if self.max_adjacency_list_length < graph2.max_adjacency_list_length:
            self._grow_adjacency_list(graph2.max_adjacency_list_length)
        elif self.max_adjacency_list_length > graph2.max_adjacency_list_length:
            graph2._grow_adjacency_list(self.max_adjacency_list_length)

        # If there's an intersection between the graphs
        intersect = np.in1d(graph2.nodes[0:graph2.number_of_nodes,0], self.nodes[:,0])
        if intersect.any():
            index_list, = np.where(intersect)
            duplicates = len(index_list)

            self.nodes = np.concatenate((self.nodes, graph2.nodes[0:index_list[0],:]))
            self.node_list = np.concatenate((self.node_list, graph2.node_list[0:index_list[0]]))
            self.node_id_list = np.concatenate((self.node_id_list, graph2.node_id_list[0:index_list[0]]))
            for i in range(duplicates-1):
                self.nodes = np.concatenate((self.nodes, graph2.nodes[index_list[i]+1:index_list[i+1],:]))
                self.node_list = np.concatenate((self.node_list, graph2.node_list[index_list[i]+1:index_list[i+1]]))
                self.node_id_list = np.concatenate((self.node_id_list, graph2.node_id_list[index_list[i]+1:index_list[i+1]]))
            self.nodes = np.concatenate((self.nodes, graph2.nodes[index_list[-1]+1:,:]))
            self.node_list = np.concatenate((self.node_list, graph2.node_list[index_list[-1]+1:]))
            self.node_id_list = np.concatenate((self.node_id_list, graph2.node_id_list[index_list[-1]+1:]))

            self.number_of_nodes += graph2.number_of_nodes - duplicates
            self.max_number_of_nodes += graph2.max_number_of_nodes - duplicates

            # Add edges from duplicate nodes
            for i in index_list:
                for j in range(1, graph2.max_adjacency_list_length):
                    if graph2.nodes[i,j] == 0:
                        break
                    node1 = graph2.nodes[i,0]
                    node2 = graph2.nodes[i,j]
                    index1 = self._find_node(node1)
                    index2 = self._find_node(node2)

                    for k in range(1, self.max_adjacency_list_length):
                        if self.nodes[index1,k] == node2:
                            break
                        elif self.nodes[index1,k] == 0:
                            self.nodes[index1,k] = node2
                            break
                        elif i == self.max_adjacency_list_length-1:
                            self._grow_adjacency_list(0)
                            self.nodes[index1,k+1] = node2

                    for k in range(1, self.max_adjacency_list_length):
                        if self.nodes[index2,k] == node1:
                            break
                        elif self.nodes[index2,k] == 0:
                            self.nodes[index2,k] = node1
                            break
                        elif i == self.max_adjacency_list_length-1:
                            self._grow_adjacency_list(0)
                            self.nodes[index2,k+1] = node1


        else: # no intersection
            self.nodes = np.concatenate((self.nodes, graph2.nodes))
            self.node_list = np.concatenate((self.node_list, graph2.node_list))
            self.node_id_list = np.concatenate((self.node_id_list, graph2.node_id_list))

            self.number_of_nodes += graph2.number_of_nodes
            self.max_number_of_nodes += graph2.max_number_of_nodes

    def edges_iter(self, yield_nodes):
        """Return an iterator over the edges between nodes in yield_nodes, a generator.
        Edges are returned as tuples.
        """
        node_list = []
        edge_list = []
        for node in yield_nodes:
            index = self._find_node(node)
            if index == -1:
                continue
            node_list.append(id(node))

            for i in range(1, self.max_adjacency_list_length):
                if self.nodes[index,i] == 0:
                    break
                if self.nodes[index,i] in node_list:
                    node1 = self._get_node(self.nodes[index,i])
                    node2 = self._get_node(self.nodes[index,0])
                    edge_list.append((node1, node2))

        return iter(edge_list)

    def neighbors(self, check_node):
        """Return a list of neighbors connected to check_node.
        Returns an empty list if check_node is not in BondGraph
        """
        neighbor_list = []

        index = self._find_node(check_node)
        if index == -1:
            return neighbor_list

        for i in range(1, self.max_adjacency_list_length):
            if self.nodes[index,i] == 0:
                break
            neighbor_list.append(self._get_node(self.nodes[index,i]))
        return neighbor_list

#
# The remaining methods serve mainly for testing, debuging, and other diagnostic
# purposes; and at the time of this writing are not used in any mbuild.Compound
# functionality. This author recommends against accessing the bond graph except
# via the compound class; and for that purpose the above methods are currently
# sufficient

    def number_of_edges(self):
        """Returns the number of edges in BondGraph. Does not count edges twice.
        """
        edges = 0
        for i in range(self.number_of_nodes):
            for j in range(1, self.max_adjacency_list_length):
                if self.nodes[i,j] == 0:
                    break
                edges += 1

        return edges/2 # Since this is not a digraph, this will give the real number of edges

    def __len__(self):
        """Returns number of nodes. Used by python's built-in function len()
        """
        return self.number_of_nodes

    def is_directed(self):
        """Returns True if BondGraph is directed.
        BondGraph is always undirected, so this will always return False.
        """
        return False

    def is_multigraph(self):
        """Returns True if BondGraph is a multigraph.
        BondGraph is never a multigraph, so this will always return False.
        """
        return False

    def connected_components(self):
        """Returns a generator of lists of nodes for each component of BondGraph
        """
        already_used = []
        for i in range(self.number_of_nodes):
            if self._get_node(self.nodes[i,0]) not in already_used:
                component_list = []
                component_list.append(self._get_node(self.nodes[i,0]))
                self._component(i, component_list)
                already_used.extend(component_list)
                yield component_list
                if len(already_used) == self.number_of_nodes:
                    break

    def _component(self, index, component_list=[]):
        """A helper function for connected_components()
        Returns a list of nodes an index is connected to, including itself, with
        an option to build onto an existing list component_list
        """
        for i in range(1, self.max_adjacency_list_length):
            if self.nodes[index,i] == 0:
                break
            current_node = self._get_node(self.nodes[index,i])
            if current_node not in component_list:
                component_list.append(current_node)
                new_index = self._find_node(current_node)
                self._component(new_index, component_list)
