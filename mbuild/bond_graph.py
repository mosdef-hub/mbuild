#TODO: Better looping; cythonizing

import numpy as np

class BondGraph():
    def __init__(self):
        """An alternative to networkx graphs.
        """
        self.numberOfNodes = 0
        self.maxNumberOfNodes = 20
        self.maxAdjListLen = 20 # Maximum length of an adjacency list, including initial node
        arrayShape = (self.maxNumberOfNodes, self.maxAdjListLen)
        self.nodes = np.ndarray(shape=arrayShape, dtype=object)

    def _growNodes(self, newSize=None):
        """Grow the nodes array and the maximum number of nodes, by a factor of
        two, or to newSize. If newSize is less than maxNumberOfNodes, the nodes
        array will be shrunk rather than grown. Note: newSize must be greater
        than numberOfNodes, otherwise data will be lost.
        """
        oldMaxSize = self.maxNumberOfNodes

        if newSize is not None:
            self.maxNumberOfNodes = newSize
        else:
            self.maxNumberOfNodes *= 2
            self.maxNumberOfNodes += 1

        if oldMaxSize < self.maxNumberOfNodes: # growing array
            arrayShape = (self.maxNumberOfNodes-oldMaxSize, self.maxAdjListLen)
            emptyArray = np.ndarray(shape=arrayShape, dtype=object)
            self.nodes = np.concatenate((self.nodes, emptyArray), axis=0)
        elif oldMaxSize > self.maxNumberOfNodes: # shrinking array
            self.nodes = self.nodes[0:self.maxNumberOfNodes, :]

    def _growAdjList(self, newSize=None):
        """Grow the nodes array and the maximum adjacency list length, by a
        factor of two, or to newSize. Note: newSize must be greater than
        maxAdjListLen - attempts to shrink the array will do nothing.
        """
        oldMaxSize = self.maxAdjListLen

        if newSize is not None:
            self.maxAdjListLen = newSize
        else:
            self.maxAdjListLen *= 2
            self.maxAdjListLen += 1

        if oldMaxSize < self.maxAdjListLen: # growing array
            arrayShape = (self.maxNumberOfNodes, self.maxAdjListLen-oldMaxSize)
            emptyArray = np.ndarray(shape=arrayShape, dtype=object)
            self.nodes = np.concatenate((self.nodes, emptyArray), axis=1)
        else:
            self.maxAdjListLen = oldMaxSize

    def has_node(self, checkNode):
        """Check if checkNode is in BondGraph.
        Return true if checkNode is in BondGraph; false otherwise.
        """
        for i in range(self.numberOfNodes):
            if checkNode == self.nodes[i][0]:
                return True
        return False

    def _find_node(self, checkNode):
        """Finds a node in BondGraph, and returns its row's index.
        If node doesn't exist, returns -1.
        """
        for i in range(self.numberOfNodes):
            if checkNode == self.nodes[i][0]:
                return i
        return -1

    def add_node(self, newNode):
        """Add newNode to BondGraph.
        If newNode is already in BondGraph, does nothing.
        """
        if has_node(newNode):
            return
        self._add_node(newNode)

    def _add_node(self, newNode):
        """Add newNode to BondGraph.
        Assumes newNode is not already in BondGraph.
        """
        if self.numberOfNodes == self.maxNumberOfNodes:
            self._growNodes()
        try:
            self.nodes[self.numberOfNodes][0] = newNode
        except IndexError:
            print("size: ", self.nodes.shape)
        self.numberOfNodes += 1

    def has_edge(self, node1, node2):
        """Check if edge between node1 and node2 is in BondGraph.
        Return true if edge exists; false otherwise.
        """
        for i in range(0, self.numberOfNodes):
            if node1 == self.nodes[i][0]:
                for j in range(1,self.maxAdjListLen):
                    if node2 == self.nodes[i][j]:
                        return True
                return False
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
            index1 = self.numberOfNodes-1
        index2 = self._find_node(node2)
        if index2 == -1:
            self._add_node(node2)
            index2 = self.numberOfNodes-1

        for i in range(1, self.maxAdjListLen):
            if self.nodes[index1][i] is None or self.nodes[index1][i] is 0:
                self.nodes[index1][i] = node2
                break
            elif i == self.maxAdjListLen-1:
                self._growAdjList()
                self.nodes[index1][i+1] = node2

        for i in range(1, self.maxAdjListLen):
            if self.nodes[index2][i] is None or self.nodes[index2][i] is 0:
                self.nodes[index2][i] = node1
                break
            elif i == self.maxAdjListLen-1:
                self._growAdjList()
                self.nodes[index2][i+1] = node1

    def remove_edge(self, node1, node2):
        """Remove edge between node1 and node2 from BondGraph.
        If edge does not exist, does nothing.
        """
        if not self.has_edge(node1, node2):
            return
        index1 = self._find_node(node1)
        index2 = self._find_node(node2)

        replace = False
        for i in range(1, self.maxAdjListLen):
            if replace:
                self.nodes[index1][i-1] = self.nodes[index1][i]
            if self.nodes[index1][i] == None or self.nodes[index1][i] == 0:
                break
            if node2 == self.nodes[index1][i]:
                replace = True
                self.nodes[index1][i] = None

        replace = False
        for i in range(1, self.maxAdjListLen):
            if replace:
                self.nodes[index2][i-1] = self.nodes[index2][i]
            if self.nodes[index2][i] == None or self.nodes[index2][i] == 0:
                break
            if node1 == self.nodes[index2][i]:
                replace = True
                self.nodes[index2][i] = None

    def remove_node(self, delNode):
        """Remove delNode from BondGraph.
        If delNode is not in BondGraph, does nothing.
        If delNode is a part of any edges, those edges will be removed.
        """
        index = self._find_node(delNode)
        if index == -1:
            return

        for i in range(self.maxAdjListLen-1, 0, -1): # Iterating backwards to avoid cost of back-shifting every time
            if self.nodes[index][i] is not None and self.nodes[index][i] is not 0:
                self.remove_edge(delNode, self.nodes[index][i]) #Make a _remove_edge() method akin to _add_node()

        for i in range(index+1, self.numberOfNodes):
            self.nodes[i-1] = self.nodes[i]
        self.numberOfNodes -= 1
        #self._growNodes(self.maxAdjListLen-1) #I'm not sure what this is for

    def compose(self, graph2):
        """Append a bond graph to the end of this one.
        Note: the two bond graphs are expected to have no nodes in common.
        Note: graph2 is modified in this method, so don't pass a graph as an
        argument that will be used later - this is for bringing two graphs
        together permanently.
        """
        if self.maxNumberOfNodes != self.numberOfNodes:
            self._growNodes(self.numberOfNodes)

        if self.maxAdjListLen < graph2.maxAdjListLen:
            self._growAdjList(graph2.maxAdjListLen)
        elif self.maxAdjListLen > graph2.maxAdjListLen:
            graph2._growAdjList(self.maxAdjListLen)

        self.nodes = np.concatenate((self.nodes, graph2.nodes))
        self.numberOfNodes += graph2.numberOfNodes
        self.maxNumberOfNodes += graph2.maxNumberOfNodes

    def edges_iter(self, yieldNodes):
        """Return an iterator over the edges between nodes in yieldNodes, a generator.
        Edges are returned as tuples.
        """
        nodeList = []
        edgeList = []
        for node in yieldNodes:
            index = self._find_node(node)
            if index == -1:
                continue
            nodeList.append(node)

            for i in range(1, self.maxAdjListLen):
                if self.nodes[index][i] == None or self.nodes[index][i] == 0:
                    break
                if self.nodes[index][i] in nodeList:
                    edgeList.append((self.nodes[index][0], self.nodes[index][i]))

        return iter(edgeList)

    def neighbors(self, checkNode):
        """Return a list of neighbors connected to checkNode.
        Returns an empty list if checkNode is not in BondGraph
        """
        neighborList = []

        index = self._find_node(checkNode)
        if index == -1:
            return neighborList

        for i in range(1, self.maxAdjListLen):
            if self.nodes[index][i] == None or self.nodes[index][i] == 0:
                break
            neighborList.append(self.nodes[index][i])
        return neighborList

# The remaining methods serve mainly for testing, debugging, and other diagnostic
# purposes; and at the time of this writing are not used in any mbuild.Compound
# functionality. This author recommends against accessing the bond graph except
# via the compound class; and for that purpose the above methods are sufficient

    def number_of_edges(self):
        """Returns the number of edges in BondGraph. Does not count edges twice.
        """
        edges = 0
        for i in range(self.numberOfNodes):
            for j in range(1, self.maxAdjListLen):
                if self.nodes[i][j] is None or self.nodes[i][j] is 0:
                    break
                edges += 1

        return edges/2 # Since this is not a digraph, this will give the real number of edges

    def __len__(self):
        """Returns number of nodes. Used by python's built-in function len()
        """
        return self.numberOfNodes

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
        alreadyUsed = []
        for i in range(self.numberOfNodes):
            if self.nodes[i][0] not in alreadyUsed:
                componentList = []
                self._component(i, componentList)
                alreadyUsed.extend(componentList)
                yield componentList

    def _component(self, index, componentList=[]):
        """A helper function for connected_components()
        Returns a list of nodes index is connected to, with an option to build
        on an existing list componentList
        """
        for i in range(1, self.maxAdjListLen):
            if self.nodes[index][i] is None or self.nodes[index][i] is 0:
                break
            if self.nodes[index][i] not in componentList:
                componentList.append(self.nodes[index][i])
                newIndex = self._find_node(self.nodes[index][i])
                self._component(newIndex, componentList)
