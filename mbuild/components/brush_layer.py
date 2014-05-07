from __future__ import division
import time
from copy import deepcopy
from numpy import pi
import sys

from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.xyz import Xyz
from mbuild.lammps import Lammps

from mbuild.ff.opls_rules import OplsRules
from mbuild.ff.opls_forcefield import OplsForceField
from mbuild.tiled_compound import TiledCompound
from mbuild.treeview import TreeView
import mbuild.unit as units

from mbuild.components.brush import Brush
from mbuild.components.surface import Surface


class BrushLayer(Compound):
    """
    """

    def __init__(self, ctx=None, tile_x=1, tile_y=1, chain_length=4, alpha=pi / 4, coverage=1):
        """
        """

        if not ctx: ctx = {}

        super(BrushLayer, self).__init__(ctx=ctx)

        surface = Surface(ctx=ctx)

        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")

        self.add(tc, 'tiled_surface')

        chains_on_surface = 0
        brush_proto = Brush(chain_length=chain_length, alpha=alpha)

        n_ports = sum(isinstance(part, Port) for part in self.tiled_surface.references.values())
        for port in self.tiled_surface.references.values():
            current_coverage = (chains_on_surface / n_ports ) * 100
            # Build a pMPC brush.
            if isinstance(port, Port) and current_coverage <  coverage:
                brush = deepcopy(brush_proto)
                brush.transform([(brush.port, port)])
                self.add(brush)
                chains_on_surface += 1
            elif current_coverage >= coverage:
                break

if __name__ == "__main__":
    profile = False
    if len(sys.argv) > 1 and sys.argv[1] == 'profile':
        profile = True

    print "Generating model..."
    start = time.time()
    m = BrushLayer(chain_length=1, alpha=pi/4, coverage=3, tile_x=1, tile_y=1)
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


    print "Saving..."
    start = time.time()
    Xyz.save(m, 'brush_layer.xyz')
    Lammps.save(m, ff, 'brush_layer.lmp')
    print "Done. ({0:.2f} s)".format(time.time() - start)

    print "Visualizing..."
    from mbuild.plot import Plot
    Plot(m, bonds=True, angles=True).show()
    #
    # tv = TreeView(m)
    # tv.show()
