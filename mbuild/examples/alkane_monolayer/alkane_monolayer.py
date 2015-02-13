from mbuild.compound import Compound
from mbuild.tools.tiled_compound import TiledCompound
from mbuild.tools.mask import apply_mask

from alkylsilane import AlkylSilane
from mbuild.components.surfaces.betacristobalite import Betacristobalite
from mbuild.components.atoms.H import H


class AlkaneMonolayer(Compound):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=10):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Args:
            tile_x (int): Number of times to replicate substrate in x-direction.
            tile_y (int): Number of times to replicate substrate in y-direction.
            chain_length (int): Number of carbon atoms per chain.
            mask (np.ndarray): A 2D array of binding locations.
        """
        super(AlkaneMonolayer, self).__init__()

        surface = Betacristobalite()
        # Replicate the surface.
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")
        self.add(tc, 'tiled_surface')

        alkylsilane = AlkylSilane(chain_length)
        hydrogen = H()

        # Attach chains to specified binding sites. Other sites get a hydrogen.
        apply_mask(self.tiled_surface, alkylsilane, mask, backfill=hydrogen)


def main():
    from mbuild.tools.mask import grid_mask_2d
    mask = grid_mask_2d(8, 8)  # Evenly spaced, 2D grid of points.
    monolayer = AlkaneMonolayer(chain_length=10, mask=mask)

    monolayer = monolayer.to_trajectory()  # Convert from mbuild to mdtraj
    monolayer.topology.find_forcefield_terms()  # Find angles/dihedrals from bonds.

    monolayer.save(filename='data.c10-n64')  # Print a LAMMPS data file.
    monolayer.save(filename='c10-n64.pdb')  # Print a .pdb file.

if __name__ == "__main__":
    main()

