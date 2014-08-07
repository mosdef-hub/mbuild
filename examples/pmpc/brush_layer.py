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
from mbuild.tools import *

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

        apply_mask(self.tiled_surface, brush_proto, mask)




if __name__ == "__main__":
    import time
    profile = False
    if len(sys.argv) > 1 and sys.argv[1] == 'profile':
        profile = True

    print "Generating model..."
    start = time.time()

    # mask = random_mask_2d(4)
    mask = grid_mask_2d(3,3)

    print mask

    m = BrushLayer(chain_length=5, alpha=pi/4, mask=mask, tile_x=2, tile_y=2)
    print "Done. ({0:.2f} s)".format(time.time() - start)


    #print "Visualizing..."
    #from mbuild.plot import Plot
    #from mayavi import mlab

    print len([a for a in m.atoms()])

    Plot(m, bonds=True, angles=False, dihedrals=False, periodic_bonds=True).show()

    # from mbuild.treeview import TreeView
    # tv = TreeView(m)
    # tv.show()

