import numpy as np
cimport numpy as np
import ctypes

very_long = np.int64 # runtime type
ctypedef np.int64_t very_long_ct # compile-time type
ctypedef Py_ssize_t index_size

cdef class BondGraph:
    cdef public index_size number_of_nodes
    cdef public index_size max_number_of_nodes
    cdef public index_size max_adjacency_list_length
    cdef public object nodes

    def __init__(self):
        """An efficiency-oriented alternative to networkx graphs for keeping
        track of bonds between atoms.
        """
        self.number_of_nodes = 0
        self.max_number_of_nodes = 10
        self.max_adjacency_list_length = 10

        # self.nodes is the primary container of this class. It is an array of
        # very long integers, (primarily) holding object IDs, keeping track of
        # which nodes are connected.
        # Every row is for a different node. The first element of a row is the
        # node ID stored in that row. The second element of a row is the number
        # of non-garbage elements in that row; that is, the number of edges that
        # row's node is connected to, plus two. Every subsequent element is the
        # ID of every node that node is connected to. Any elements beyond the
        # number of edges stored are garbage and should not be looked up.
        self.nodes = np.zeros(
                shape=(self.max_number_of_nodes, self.max_adjacency_list_length),
                dtype=very_long)

        # cdef very_long_ct [:,:] nodes_view = self.nodes # I don't know how to use this

    cdef void _grow_nodes(self, index_size new_size):
        """Grow the nodes array and the maximum number of nodes to new_size, or
        by a factor of two, if new_size is 0. If new_size is less than
        max_number_of_nodes, the nodes array will be shrunk rather than grown.
        Note: new_size must be greater than number_of_nodes, otherwise data will
        be lost.
        """
        cdef index_size old_max_size = self.max_number_of_nodes
        cdef np.ndarray[very_long_ct, ndim=2] empty_array

        if new_size != 0:
            self.max_number_of_nodes = new_size
        else:
            self.max_number_of_nodes *= 2
            self.max_number_of_nodes += 1

        # growing array
        if old_max_size < self.max_number_of_nodes:
            empty_array = np.zeros(
                    shape=(self.max_number_of_nodes-old_max_size, self.max_adjacency_list_length),
                    dtype=very_long)
            self.nodes = np.concatenate((self.nodes, empty_array), axis=0)

        # shrinking array
        elif old_max_size > self.max_number_of_nodes:
            self.nodes = self.nodes[0:self.max_number_of_nodes, :]

    cpdef void _grow_adjacency_list(self, index_size new_size):
        """Grow the nodes array and the maximum adjacency list length to
        new_size, or by a factor of two, if new_size is 0. new_size must be
        greater than max_adjacency_list_length - attempts to shrink the array
        will do nothing.
        This method can be called from cython or python, so when a BondGraph
        object is passed to a BondGraph method, as with compose(), this
        method can be called on that object. This comes with the price of slight
        overhead.
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

    cdef object _get_node(self, very_long_ct node_id):
        """Given a node's id (as is stored in self.nodes), finds and returns the
        object that id points to, the object BondGraph was asked to store.
        """
        return ctypes.cast(node_id, ctypes.py_object).value

    cdef bint _has_node(self, very_long_ct check_node):
        """Check if object ID check_node is in BondGraph.
        Return true if check_node is in BondGraph; false otherwise.
        """
        if check_node in self.nodes[0:self.number_of_nodes,0]:
            return True
        return False

    def has_node(self, check_node):
        """Check if check_node is in BondGraph.
        Return true if check_node is in BondGraph; false otherwise.
        """
        return self._has_node(id(check_node))

    cdef index_size _find_node(self, very_long_ct node_id):
        """Finds a node ID in BondGraph, and returns its row's index.
        If node doesn't exist, returns -1.
        """
        index, = np.where(self.nodes[0:self.number_of_nodes,0]==node_id)
        if len(index) == 0:
            return -1
        return index[0]

    def add_node(self, new_node):
        """Add new_node to BondGraph.
        If new_node is already in BondGraph, does nothing.
        """
        cdef very_long_ct new_node_id = id(new_node)
        if self._has_node(new_node_id):
            return
        self._add_node(new_node_id)

    cdef void _add_node(self, very_long_ct new_node):
        """Add object ID new_node to BondGraph.
        Assumes new_node is not already in BondGraph.
        """
        if self.number_of_nodes == self.max_number_of_nodes:
            self._grow_nodes(0)

        self.nodes[self.number_of_nodes,0] = new_node
        self.nodes[self.number_of_nodes,1] = 2
        self.number_of_nodes += 1

    def has_edge(self, node1, node2):
        """Check if edge between node1 and node2 is in BondGraph.
        Return true if edge exists; false otherwise.
        """
        return self._has_edge(id(node1), id(node2))

    cdef bint _has_edge(self, very_long_ct node1, very_long_ct node2):
        """Check if edge between object IDs node1 and node2 is in BondGraph.
        Return true if edge exists; false otherwise.
        """
        cdef index_size index = self._find_node(node1)
        if index == -1:
            return False

        cdef index_size i
        for i in range(2, self.nodes[index,1]):
            if self.nodes[index,i]==node2:
                return True
        return False

    def add_edge(self, node1, node2):
        """Adds edge between node1 and node2
        If node1 or node2 aren't in BondGraph, they are created.
        If an edge between node1 and node2 already exists, does nothing.
        """
        cdef very_long_ct node1_id = id(node1)
        cdef very_long_ct node2_id = id(node2)
        cdef index_size index1, index2
        cdef very_long_ct node1_num_edges, node2_num_edges

        if self._has_edge(node1_id, node2_id):
            return

        index1 = self._find_node(node1_id)
        if index1 == -1:
            self._add_node(node1_id)
            index1 = self.number_of_nodes-1
        index2 = self._find_node(node2_id)
        if index2 == -1:
            self._add_node(node2_id)
            index2 = self.number_of_nodes-1

        node1_num_edges = self.nodes[index1,1]
        node2_num_edges = self.nodes[index2,1]

        if(node1_num_edges == self.max_adjacency_list_length or
           node2_num_edges == self.max_adjacency_list_length):
            self._grow_adjacency_list(0)

        self.nodes[index1, node1_num_edges] = self.nodes[index2,0]
        self.nodes[index2, node2_num_edges] = self.nodes[index1,0]
        self.nodes[index1,1] += 1
        self.nodes[index2,1] += 1

    def remove_edge(self, node1, node2):
        """Remove edge between node1 and node2 from BondGraph.
        If edge does not exist, does nothing.
        """
        self._remove_edge(id(node1), id(node2))

    cdef _remove_edge(self, very_long_ct node1, very_long_ct node2):
        """Remove edge between object IDs node1 and node2 from BondGraph.
        If edge does not exist, does nothing.
        """
        cdef bint replace
        cdef index_size index1, index2
        cdef index_size i

        if not self._has_edge(node1, node2):
            return
        index1 = self._find_node(node1)
        index2 = self._find_node(node2)

        replace = False
        self.nodes[index1,1] -= 1
        for i in range(2, self.nodes[index1,1]):
            if not replace and node2 == self.nodes[index1,i]:
                replace = True
            if replace:
                self.nodes[index1,i] = self.nodes[index1,i+1]

        replace = False
        self.nodes[index2,1] -= 1
        for i in range(2, self.nodes[index2,1]):
            if not replace and node1 == self.nodes[index2,i]:
                replace = True
            if replace:
                self.nodes[index2,i] = self.nodes[index2,i+1]

    def remove_node(self, del_node):
        """Remove del_node from BondGraph.
        If del_node is not in BondGraph, does nothing.
        If del_node is a part of any edges, those edges will be removed.
        """
        cdef very_long_ct del_node_id = id(del_node)
        cdef index_size index = self._find_node(del_node_id)
        cdef index_size i

        if index == -1:
            return

        for i in range(self.nodes[index,1]-1, 1, -1): # Iterating backwards to avoid cost of back-shifting every time
            self._remove_edge(del_node_id, self.nodes[index,i])

        for i in range(index+1, self.number_of_nodes):
            self.nodes[i-1] = self.nodes[i]
        self.number_of_nodes -= 1

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
        #cdef np.ndarray[np.uint8_t, ndim=1] intersect # I can't seem to get this to work
        cdef np.ndarray[index_size, ndim=1] index_list
        cdef index_size duplicates
        cdef very_long_ct node1, node2
        cdef index_size index1, index2
        cdef very_long_ct node1_num_edges, node2_num_edges
        cdef index_size i, j, k
        cdef bint already_has_edge
        intersect = np.in1d(graph2.nodes[0:graph2.number_of_nodes,0], self.nodes[:,0])
        if intersect.any():
            index_list, = np.where(intersect)
            duplicates = len(index_list)

            self.nodes = np.concatenate((self.nodes, graph2.nodes[0:index_list[0],:]))
            for i in range(duplicates-1):
                self.nodes = np.concatenate((self.nodes, graph2.nodes[index_list[i]+1:index_list[i+1],:]))
            self.nodes = np.concatenate((self.nodes, graph2.nodes[index_list[-1]+1:,:]))

            self.number_of_nodes += graph2.number_of_nodes - duplicates
            self.max_number_of_nodes += graph2.max_number_of_nodes - duplicates

            # Add edges from duplicate nodes
            for i in index_list:
                for j in range(2, graph2.nodes[i,1]):
                    node1 = graph2.nodes[i,0]
                    node2 = graph2.nodes[i,j]
                    index1 = self._find_node(node1)
                    index2 = self._find_node(node2)
                    node1_num_edges = self.nodes[index1,1]
                    node2_num_edges = self.nodes[index2,1]

                    already_has_edge = False
                    for k in range(2, node1_num_edges):
                        if self.nodes[index1,k]==node2:
                            already_has_edge = True
                            break
                    if not already_has_edge:
                        if node1_num_edges == self.max_adjacency_list_length:
                            self._grow_adjacency_list(0)

                        self.nodes[index1, node1_num_edges] = self.nodes[index2,0]
                        self.nodes[index1,1] += 1

                    already_has_edge = False
                    for k in range(2, node2_num_edges):
                        if self.nodes[index2,k]==node1:
                            already_has_edge = True
                            break
                    if not already_has_edge:
                        if node2_num_edges == self.max_adjacency_list_length:
                            self._grow_adjacency_list(0)

                        self.nodes[index2, node2_num_edges] = self.nodes[index1,0]
                        self.nodes[index2,1] += 1

        else: # no intersection
            self.nodes = np.concatenate((self.nodes, graph2.nodes))

            self.number_of_nodes += graph2.number_of_nodes
            self.max_number_of_nodes += graph2.max_number_of_nodes

    def edges_iter(self, yield_nodes):
        """Return an iterator over the edges between nodes in yield_nodes, a generator.
        Edges are returned as tuples.
        """
        node_set = set()
        edge_list = []
        cdef object node
        cdef very_long_ct node_id
        cdef index_size index

        for node in yield_nodes:
            node_id = id(node)
            index = self._find_node(node_id)
            if index == -1:
                continue
            node_set.add(node_id)

            for i in range(2, self.nodes[index,1]):
                if self.nodes[index,i] in node_set:
                    node1 = self._get_node(self.nodes[index,i])
                    node2 = self._get_node(self.nodes[index,0])
                    edge_list.append((node1, node2))

        return iter(edge_list)

    def neighbors(self, check_node):
        """Return a list of neighbors connected to check_node.
        Returns an empty list if check_node is not in BondGraph
        """
        cdef index_size index = self._find_node(id(check_node))
        neighbor_list = []

        if index == -1:
            return neighbor_list

        cdef index_size i
        for i in range(2, self.nodes[index,1]):
            neighbor_list.append(self._get_node(self.nodes[index,i]))
        return neighbor_list

# The remaining methods serve mainly for testing, debuging, and other diagnostic
# purposes; and at the time of this writing are not used in any mbuild.Compound
# functionality. This author recommends against accessing the bond graph except
# via the compound class; and for that purpose the above methods are currently
# sufficient

    def number_of_edges(self):
        """Returns the number of edges in BondGraph. Does not count edges twice.
        """
        cdef index_size edges = 0
        cdef index_size i
        for i in range(self.number_of_nodes):
            edges += self.nodes[i,1]-2

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
        cdef index_size i

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
        cdef object current_node
        cdef index_size new_index
        cdef index_size i

        for i in range(2, self.nodes[index,1]):
            current_node = self._get_node(self.nodes[index,i])
            if current_node not in component_list:
                component_list.append(current_node)
                new_index = self._find_node(id(current_node))
                self._component(new_index, component_list)
