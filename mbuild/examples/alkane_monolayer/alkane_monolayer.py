from mbuild.atom import Atom
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tools.tiled_compound import TiledCompound
from mbuild.tools.mask import apply_mask

from alkylsilane import AlkylSilane
from surface import Surface


class AlkaneMonolayer(Compound):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, tile_x=1, tile_y=1, chain_length=10, mask=None):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Args:
            tile_x (int): Number of times to replicate substrate in x-direction.
            tile_y (int): Number of times to replicate substrate in y-direction.
            chain_length (int): Number of carbon atoms per chain.
            mask (np.ndarray): A 2D array of binding locations.
        """
        super(AlkaneMonolayer, self).__init__()

        surface = Surface()
        # Replicate the surface
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")

        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        alkylsilane = AlkylSilane(chain_length)
        # Create a hydrogen Atom with a Port. See issue #21 on github.
        hydrogen = Compound()
        hydrogen.add(Atom('H'), 'H')
        hydrogen.add(Port(anchor=hydrogen.H), 'down')

        # If no binding sites are provided, attach chains to all sites.
        if mask is None:
            from mbuild.tools.mask import random_mask_2d
            mask = random_mask_2d(tile_x * 10 + tile_y * 10)

        # Attach chains to specified binding sites. Other sites get a hydrogen.
        apply_mask(self.tiled_surface, alkylsilane, mask,
                   guest_port_name='down', backfill=hydrogen)


def main():
    from mbuild.tools.mask import grid_mask_2d
    mask = grid_mask_2d(8, 8)  # Evenly spaced, 2D grid of points.
    monolayer = AlkaneMonolayer(chain_length=10, mask=mask)
    monolayer.visualize()  # Take a peak (if you have VMD installed).

    monolayer = monolayer.to_trajectory()  # Convert from mbuild to mdtraj
    monolayer.topology.find_forcefield_terms()  # Find angles/dihedrals from bonds.

    monolayer.save(filename='data.c10-n64')  # Print a LAMMPS data file.
    monolayer.save(filename='c10-n64.pdb')  # Print a .pdb file.

if __name__ == "__main__":
    main()

