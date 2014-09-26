from __future__ import division

from mbuild.atom import Atom
from mbuild.formats.lammps import save_lammps
from mbuild.port import Port
from mbuild.compound import Compound
from mbuild.tiled_compound import TiledCompound
from mbuild.plugins.mask import apply_mask

from mbuild.examples.alkane_monolayers.alkylsilane import AlkylSilane
from mbuild.examples.alkane_monolayers.surface import Surface


class AlkaneMonolayer(Compound):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, tile_x=1, tile_y=1, chain_length=10, mask=None):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Args:
            tile_x (int):
            tile_y (int):
            chain_length (int):
            mask (np.ndarray):
        """
        super(AlkaneMonolayer, self).__init__()

        surface = Surface()
        tc = TiledCompound(surface, tile_x, tile_y, 1, kind="tiled_surface")

        self.add(tc, 'tiled_surface')
        self.periodicity = tc.periodicity

        alkylsilane = AlkylSilane(chain_length)
        # TODO: Library of all elements with ports attached.
        hydrogen = Compound('H')
        hydrogen.add(Atom('H'), 'H')
        hydrogen.add(Port(anchor=hydrogen.H), 'down')
        apply_mask(self.tiled_surface, alkylsilane, mask,
                   guest_port_name='down', backfill=hydrogen)


if __name__ == "__main__":
    from mbuild.plugins.mask import grid_mask_2d
    mask = grid_mask_2d(10, 10)  # Evenly spaced, 2D grid of points.
    monolayer = AlkaneMonolayer(chain_length=6, mask=mask)
    monolayer = monolayer.to_trajectory()  # Convert from mBuild to mdtraj

    monolayer.top.find_forcefield_terms()  # Find angles/dihedrals from bonds.

    save_lammps(monolayer, filename='data.c6-100')  # Print a LAMMPS data file.

    print monolayer.n_atoms
    print monolayer.top.n_ff_bonds
    print monolayer.top.n_ff_angles
    print monolayer.top.n_ff_dihedrals



