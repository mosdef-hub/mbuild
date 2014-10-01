from copy import deepcopy

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x
from mbuild.tools.mask import grid_mask_2d
from mbuild.tools.solvent import solvent_box


class Bilayer(Compound):
    """ """
    def __init__(self, lipid, ref_atom=0, n_lipids_x=10, n_lipids_y=10, apl=1.0,
                 solvent=None, host_box=None, spacing_z=0.5):
        """

        Args:
            lipid (Compound):
            ref_atom (int): Index of the atom in lipid to form the interface
                between the bilayers at. Typically an atom at bottom of lipid.
            n_lipids_x (int): Number of lipids in the x-direction per layer.
            n_lipids_y (int): Number of lipids in the x-direction per layer.
            apl (float): Area per lipid.
            solvent (Compound): Compound to solvate the bilayer with. Typically,
                a pre-equilibrated box of solvent.
            host_box (Box, optional):
            spacing_z (float, optional):
        """
        super(Bilayer, self).__init__()

        mask = grid_mask_2d(n_lipids_x, n_lipids_y)
        mask *= np.sqrt(apl * n_lipids_x * n_lipids_y)

        spacing = np.array([0, 0, spacing_z])
        for point in mask:
            top_lipid = deepcopy(lipid)
            translate(top_lipid, -top_lipid.atom[ref_atom] + spacing)
            translate(top_lipid, point)
            self.add(top_lipid)

            bot_lipid = deepcopy(lipid)
            translate(bot_lipid, -bot_lipid.atom[ref_atom] + spacing)
            rotate_around_x(bot_lipid, np.pi)
            translate(bot_lipid, point)
            self.add(bot_lipid)

        if host_box is None:
            host_box = self.boundingbox()
            host_box.mins -= np.array([0.5*np.sqrt(apl), 0.5*np.sqrt(apl), 0])
            host_box.maxs += np.array([0.5*np.sqrt(apl), 0.5*np.sqrt(apl), 0])
            host_box.lengths = host_box.lengths * np.array([1.0, 1.0, 3.0])

        lipid_box = self.boundingbox()
        top_box = Box(mins=[host_box.mins[0], host_box.mins[1], lipid_box.maxs[2]],
                      maxs=[host_box.maxs[0], host_box.maxs[1], host_box.maxs[2]])
        bot_box = Box(mins=[host_box.mins[0], host_box.mins[1], host_box.mins[2]],
                      maxs=[host_box.maxs[0], host_box.maxs[1], lipid_box.mins[2]])

        top_solvent = solvent_box(solvent, top_box)
        self.add(top_solvent)
        bottom_solvent = solvent_box(solvent, bot_box)
        self.add(bottom_solvent)


def main():
    from mbuild.trajectory import Trajectory
    from mbuild.testing.tools import get_fn


    water = Trajectory.load(get_fn('water.hoomdxml'))
    water = water.to_compound()
    ecerns = Trajectory.load(get_fn('ecer2.hoomdxml'))
    ecerns = ecerns.to_compound()

    bilayer = Bilayer(ecerns, n_lipids_x=5, n_lipids_y=5, apl=1.4,
                      solvent=water, ref_atom=0, spacing_z=1.0)

    bilayer = bilayer.to_trajectory()
    bilayer.topology.load_ff_bonds()
    bilayer.save(filename='bilayer.hoomdxml')

if __name__ == "__main__":
    main()
