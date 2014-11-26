from copy import deepcopy

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x
from mbuild.tools.mask import grid_mask_2d
from mbuild.tools.solvent import solvent_box
import pdb


class Bilayer(Compound):
    """Create a lipid bilayer and add solvent above and below. """
    def __init__(self, lipid, ref_atom=0, n_lipids_x=10, n_lipids_y=10, 
                 area_per_lipid=1.0, solvent=None, lipid_box=None, 
                 spacing_z=0.5, solvent_per_lipid=None, n_solvent=None):
        """
        Note: This example is still pretty immature as it only represents some
        brief scratch work. Feel free to flesh it out!

        Args:
            lipid (Compound):
            ref_atom (int): Index of the atom in lipid to form the interface
            between the bilayers at. Typically an atom at bottom of lipid.
            n_lipids_x (int): Number of lipids in the x-direction per layer.
            n_lipids_y (int): Number of lipids in the x-direction per layer.
            area_per_lipid (float): Area per lipid.
            solvent (Compound): Compound to solvate the bilayer with. Typically,
            a pre-equilibrated box of solvent.
            lipid_box (Box, optional): A Box containing the lipids where no
            solvent will be added.
            spacing_z (float, optional): Amount of space to add between opposing
            monolayers.
            solvent_per_lipid (int, optional): Number of solvent molecules per lipid
            n_solvent(int, optional): *Total* number of solvent molecules
        """
        super(Bilayer, self).__init__()

        # specify exactly one of either solvent_per_lipid or n_solvent
        # and use to get number of lipids per layer
        solvent_per_layer = 0
        assert not (solvent_per_lipid is None and n_solvent is None)
        if solvent_per_lipid is not None:
            assert n_solvent is None
            solvent_per_layer = n_lipids_x * n_lipids_y * solvent_per_lipid
        elif n_solvent is not None:
            assert solvent_per_lipid is None
            solvent_per_layer = n_solvent / 2

        mask = grid_mask_2d(n_lipids_x, n_lipids_y)
        mask *= np.sqrt(area_per_lipid * n_lipids_x * n_lipids_y)

        spacing = np.array([0, 0, spacing_z])
        for point in mask:
            top_lipid = deepcopy(lipid)
            translate(top_lipid, -top_lipid.atoms[ref_atom] + spacing)
            translate(top_lipid, point)
            self.add(top_lipid)

            bot_lipid = deepcopy(lipid)
            translate(bot_lipid, -bot_lipid.atoms[ref_atom] + spacing)
            rotate_around_x(bot_lipid, np.pi)
            translate(bot_lipid, point)
            self.add(bot_lipid)

        if lipid_box is None:
            lipid_box = self.boundingbox()
            # add buffer around lipid box
            lipid_box.mins -= np.array([0.5*np.sqrt(area_per_lipid),
                0.5*np.sqrt(area_per_lipid), 0.5*np.sqrt(area_per_lipid)])
            lipid_box.maxs += np.array([0.5*np.sqrt(area_per_lipid),
                0.5*np.sqrt(area_per_lipid), 0.5*np.sqrt(area_per_lipid)])
            #lipid_box.lengths = lipid_box.lengths * np.array([1.0, 1.0, 3.0])

        #self.save('lipids.hoomdxml')
        
        # add buffer space between lipid and solvent boxes
        # figure out size of solvent boxes based on number of lipids
        solvent_number_density = solvent.n_atoms / np.prod(solvent.periodicity)
        water_box_z = solvent_per_layer / (lipid_box.lengths[0]
                * lipid_box.lengths[1] * solvent_number_density)


        top_box = Box(mins=[lipid_box.mins[0], lipid_box.mins[1], lipid_box.maxs[2]],
                      maxs=[lipid_box.maxs[0], lipid_box.maxs[1], 
                          lipid_box.maxs[2]+water_box_z])
        bot_box = Box(mins=[lipid_box.mins[0], lipid_box.mins[1], 
                        lipid_box.mins[2]-water_box_z],
                      maxs=[lipid_box.maxs[0], lipid_box.maxs[1], lipid_box.mins[2]])

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

    bilayer = Bilayer(ecerns, n_lipids_x=5, n_lipids_y=5,
            area_per_lipid=1.4, solvent=water, ref_atom=0, 
            spacing_z=0.5, n_solvent=3000)

    #test_water_box = Box(mins=[-5, 6, -13], maxs=[15, 9, -8])
    #test_water = solvent_box(water, test_water_box)
    #test_water.save(filename='water_test.hoomdxml')

    bilayer = bilayer.to_trajectory()
    bilayer.topology.load_ff_bonds()
    bilayer.save(filename='bilayer.hoomdxml')

if __name__ == "__main__":
    main()
