from copy import deepcopy

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x
from mbuild.plugins.mask import grid_mask_2d
from mbuild.tools import solvate


class Bilayer(Compound):
    """ """
    def __init__(self, lipid, n_lipids_x=10, n_lipids_y=10, apl=1.0,
                 solvent=None, host_box=None, guest_box=None):
        """

        Args:
            lipid (Compound):
            n_lipids_x (int): Number of lipids in the x-direction per layer.
            n_lipids_y (int): Number of lipids in the x-direction per layer.
            apl (float): Area per lipid.
            solvent (Compound): Compound to solvate the bilayer with. Typically,
                a pre-equilibrated box of solvent.
            host_box (Box, optional):
            guest_box (Box, optional):
        """
        super(Bilayer, self).__init__()

        mask = grid_mask_2d(n_lipids_x, n_lipids_y)
        mask *= np.sqrt(apl * n_lipids_x * n_lipids_y)

        #rotate_around_x(lipid, np.pi)
        for point in mask:
            top_lipid = deepcopy(lipid)
            # TODO: figure out labeling
            translate(top_lipid, -top_lipid.C[32] + np.asarray([0, 0, 0.1]))
            translate(top_lipid, point)
            self.add(top_lipid)

        self.save('bilayer.pdb')

        for point in mask:
            bot_lipid = deepcopy(lipid)
            translate(bot_lipid, -bot_lipid.C[32] + np.asarray([0, 0, 0.1]))
            rotate_around_x(bot_lipid, np.pi)
            translate(bot_lipid, point)
            self.add(bot_lipid)
        self.save('bilayer2.pdb')

        if not host_box:
            host_box = self.boundingbox()
            host_box.lengths[0] += 0.5 * np.sqrt(apl)
            host_box.lengths[1] += 0.5 * np.sqrt(apl)
            host_box.lengths[2] *= 1.5
            host_box = Box(host_box.lengths)
        print("Bilayer: {} Box: {}{}".format(self.n_atoms(), host_box, host_box.lengths))
        print("Water: {} Box: {}{}".format(solvent.n_atoms(), guest_box, guest_box.lengths))

        solvate(ecerns, solvent, host_box, guest_box)

if __name__ == "__main__":
    from mbuild.examples.bilayer.eceramidens import ECeramideNS
    from mbuild.trajectory import Trajectory
    from mbuild.formats.hoomdxml import save_hoomdxml
    from mbuild.testing.tools import get_fn

    ecerns = ECeramideNS()

    #water = Trajectory.load(get_fn('single_water.hoomdxml'))
    water = Trajectory.load(get_fn('spc216.pdb'))
    #water_box = Box(water.unitcell_lengths[0])
    water = water.to_compound()
    water_box = water.boundingbox()

    bilayer = Bilayer(ecerns, n_lipids_x=5, n_lipids_y=5, apl=0.5,
                      solvent=water, guest_box=water_box)

    bilayer = bilayer.to_trajectory()
    bilayer.topology.load_ff_bonds()
    save_hoomdxml(bilayer, filename='bilayer.xml')
    #import os
    #os.system('vmd -hoomd bilayer.xml')

    #from mbuild.plot import Plot
    #Plot(bilayer, verbose=True, atoms=True, bonds=True).show()

    #from mbuild.plot import Plot
    #Plot(water, verbose=True, atoms=True, bonds=True).show()
    #water.append_from_file(get_fn('spc216.pdb'))
    #water.periodicity = water.boundingbox().lengths
    #water = TiledCompound(water, 6, 6, 9)

