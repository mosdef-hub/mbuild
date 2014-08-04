from __future__ import division
from itertools import ifilter
from copy import deepcopy
from numpy import pi
import sys
import numpy as np
from examples.pmpc.surface import Surface
from examples.pmpc.brush import Brush
from mbuild.bond import Bond
from mbuild.coordinate_transform import equivalence_transform
from mbuild.plot import Plot
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tiled_compound import TiledCompound

class BrushLayer(Compound):
    """
    """

    def __init__(self, tile_x=1, tile_y=1, chain_length=4,
            alpha=pi/4, mask=None):
        """
        """
        super(BrushLayer, self).__init__()

        surface = Surface()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')

        brush_proto = Brush(chain_length=chain_length, alpha=alpha)
        bbmin, bbmax, bbsize = self.boundingbox(excludeG=False)
        mask = mask * bbsize + bbmin

        n_ports = sum(isinstance(part, Port) for part in self.tiled_surface.labels.values())
        # print n_ports

        port_pos = np.empty((n_ports,3))
        port_list = []
        for pidx, port in enumerate(ifilter(lambda x: isinstance(x, Port), self.tiled_surface.labels.values())):
            port_pos[pidx, :] = port.middle.pos
            port_list.append(port)


        for mp in mask:
            closest_point_idx = np.argmin(tc.min_periodic_distance(mp, port_pos))
            closest_port = port_list[closest_point_idx]
            brush = deepcopy(brush_proto)
            equivalence_transform(brush, brush.port, closest_port)
            self.add(brush)
            self.add(Bond(brush.port.SI_1, closest_port.O))
            port_pos[closest_point_idx,:] = np.array([np.inf, np.inf, np.inf])

if __name__ == "__main__":
    import time
    profile = False
    if len(sys.argv) > 1 and sys.argv[1] == 'profile':
        profile = True

    print "Generating model..."
    start = time.time()

    # # random mask
    n_chains = 1
    mask = np.random.random((n_chains, 3))
    mask[:, 2] = 0

    # mask = np.array([[.3, .3, 0], [.5, .7, 0], [.7, .3, 0]])
    # mask = np.array([[.3, .3, 0], [.7, .3, 0]])
    # mask = np.array([[.3, .3, 0]])

    """
    # grid mask
    n = 10
    mask = np.zeros(shape=(n*n, 3), dtype=float)
    for i in range(n):
        for j in range(n):
            mask[i*n + j, 0] = i
            mask[i*n + j, 1] = j
    mask[:,0] = mask[:,0] / np.max(mask[:,0])
    mask[:,1] = mask[:,1] / np.max(mask[:,1])
    """

    m = BrushLayer(chain_length=5, alpha=pi/4, mask=mask, tile_x=20, tile_y=20)
    print "Done. ({0:.2f} s)".format(time.time() - start)


    #print "Visualizing..."
    #from mbuild.plot import Plot
    #from mayavi import mlab

    print len([a for a in m.atoms()])

    # Plot(m, bonds=True, angles=False, dihedrals=False, periodic_bonds=True).show()

    # from mbuild.treeview import TreeView
    # tv = TreeView(m)
    # tv.show()

