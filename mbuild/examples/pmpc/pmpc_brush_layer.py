from __future__ import division

from numpy import pi

from mbuild.compound import Compound
from mbuild.tools.mask import apply_mask, random_mask_2d
from mbuild.tools.parameterize.atomtyper import find_atomtypes
from mbuild.tools.tiled_compound import TiledCompound

from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.components.atoms.H import H
from pmpc import PMPC


class PMPCLayer(Compound):
    """Create a layer of grafted pMPC brushes on a beta-cristobalite surface."""
    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=4, alpha=pi/4):
        super(PMPCLayer, self).__init__()

        surface = Betacristobalite()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')

        pmpc = PMPC(chain_length=chain_length, alpha=alpha)
        hydrogen = H()

        apply_mask(self.tiled_surface, pmpc, mask, backfill=hydrogen)


if __name__ == "__main__":
    mask = random_mask_2d(1)
    pmpc_layer = PMPCLayer(mask=mask, chain_length=3, alpha=pi/4, tile_x=1, tile_y=1)
    pmpc_layer.visualize(show_ports=True)
    # find_atomtypes(pmpc_layer)
    #
    # pmpc_layer = pmpc_layer.to_trajectory(chain_types=[Betacristobalite, PMPC],
    #                                       residue_types=[MPC])
    # pmpc_layer.topology.find_forcefield_terms()
    # pmpc_layer.save(filename='brush_layer.top')
