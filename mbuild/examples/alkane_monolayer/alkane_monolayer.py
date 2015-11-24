import mbuild as mb

from mbuild.lib.surfaces import Betacristobalite
from mbuild.lib.atoms import H
from mbuild.examples.alkane_monolayer.alkylsilane import AlkylSilane


class AlkaneMonolayer(mb.Monolayer):
    """An akylsilane monolayer on beta-cristobalite. """

    def __init__(self, mask, tile_x=1, tile_y=1, chain_length=10):
        """Create an alkylsilane monolayer on beta-cristobalite.

        Parameters
        ----------
        mask : np.ndarray, shape=(n, 3), optional, default=None
            An array of planar binding locations. If not provided, the entire
            surface will be filled with `chain`.
        tile_x : int, optional, default=1
            Number of times to replicate substrate in x-direction.
        tile_y : int, optional, default=1
            Number of times to replicate substrate in y-direction.
        chain_length : int, optional, default=10
            Number of carbon atoms per chain.
        """
        surface = Betacristobalite()
        alkylsilane = AlkylSilane(chain_length)
        hydrogen = H()
        super(AlkaneMonolayer, self).__init__(surface, alkylsilane, hydrogen,
                                              mask, tile_x, tile_y)


def main():
    mask = mb.grid_mask_2d(8, 8)  # Evenly spaced, 2D grid of points.
    monolayer = AlkaneMonolayer(chain_length=10, mask=mask)
    monolayer.save(filename='c10-n64.pdb', show_ports=True)
    monolayer.visualize()

if __name__ == "__main__":
    main()

