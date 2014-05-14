from __future__ import division
from itertools import ifilter
from copy import deepcopy
from numpy import pi
import sys

import numpy as np

from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tiled_compound import TiledCompound

from mbuild.components.brush import Brush
from mbuild.components.surface import Surface

import mbuild.unit as units
from mbuild.ff.opls_rules import OplsRules
from mbuild.ff.opls_forcefield import OplsForceField

class BrushLayer(Compound):
    """
    """

    def __init__(self, ctx=None, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4, mask=None):
        """
        """
        if not ctx:
            ctx = {}
        super(BrushLayer, self).__init__(ctx=ctx)

        surface = Surface(ctx=ctx)
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')

        brush_proto = Brush(chain_length=chain_length, alpha=alpha)
        bbmin, bbmax, bbsize = self.boundingbox(excludeG=False)
        mask = mask * bbsize + bbmin

        n_ports = sum(isinstance(part, Port) for part in self.tiled_surface.references.values())
        port_pos = np.empty((n_ports,3))
        port_list = []
        for pidx, port in enumerate(ifilter(lambda x: isinstance(x, Port), self.tiled_surface.references.values())):
            port_pos[pidx, :] = port.middle.pos
            port_list.append(port)

        for mp in mask:
            closest_point_idx = np.argmin(tc.min_periodic_distance(mp, port_pos))
            closest_port = port_list[closest_point_idx]
            brush = deepcopy(brush_proto)
            brush.transform([(brush.port, closest_port)])
            self.add(brush)
            port_pos[closest_point_idx,:] = np.array([np.inf, np.inf, np.inf])

if __name__ == "__main__":
    import time
    profile = False
    if len(sys.argv) > 1 and sys.argv[1] == 'profile':
        profile = True

    print "Generating model..."
    start = time.time()

    # # random mask
    mask = np.random.random((5, 3))
    mask[:, 2] = 0
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

    m = BrushLayer(chain_length=1, alpha=pi/4, mask=mask, tile_x=1, tile_y=1)
    print "Done. ({0:.2f} s)".format(time.time() - start)

    print "Loading and pruning forcefield..."
    start = time.time()
    ff = OplsForceField()

    # TODO: add real parameters
    ff.add_atom_type(
            opls_type     = 'Si',
            bond_type     = 'SI',
            atomic_number = 14,
            mass          = 28.085 * units.amu,
            charge        = 0.84 * units.elementary_charge,
            sigma         = 3.5 * units.angstroms,
            epsilon       = 4.0 * units.kilojoules_per_mole)

    ff.add_atom_type(
            opls_type     = 'O',
            bond_type     = 'OS',
            atomic_number = 8,
            mass          = 16.0 * units.amu,
            charge        = -0.42 * units.elementary_charge,
            sigma         = 3.5 * units.angstroms,
            epsilon       = 4.0 * units.kilojoules_per_mole)
    """
    ff.add_atom_type(
            opls_type     = 'H',
            bond_type     = 'HO',
            atomic_number = 1,
            mass          = 1 * units.amu,
            charge        = 0.2 * units.elementary_charge,
            sigma         = 0.0 * units.angstroms,
            epsilon       = 0.0 * units.kilojoules_per_mole)
    """
    ff = ff.prune(m)
    print "Done. ({0:.2f} s)".format(time.time() - start)

    if profile:
        import cProfile, pstats, StringIO
        pr = cProfile.Profile()
        pr.enable()

    print "Generating topology..."
    start = time.time()
    rules = OplsRules(m, ff)
    rules.execute(verbose=False)
    print "Done. ({0:.2f} s)".format(time.time() - start)

    if profile:
        pr.disable()
        s = StringIO.StringIO()
        sortby = 'cumulative'
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print s.getvalue()

    print "\nNumber of atoms: {0}".format(len(m.getAtomListByKind('*')))
    print "Number of bonds: {0}".format(len(m.bonds))
    print "Number of angles: {0}".format(len(m.angles))
    print "Number of dihedrals: {0}".format(len(m.dihedrals))

    print "Saving..."
    from mbuild.lammps import Lammps
    from mbuild.xyz import Xyz
    start = time.time()
    Xyz.save(m, 'brush_layer.xyz', ff=ff)
    Lammps.save(m, ff, 'brush_layer.lmp')
    print "Done. ({0:.2f} s)".format(time.time() - start)

    #print "Visualizing..."
    #from mbuild.plot import Plot
    #Plot(m, bonds=True, angles=False, dihedrals=True).show()

    #from mbuild.treeview import TreeView
    #tv = TreeView(m)
    #tv.show()
