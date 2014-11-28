import warnings
from copy import deepcopy
from random import seed, randint

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x, rotate_around_y
from mbuild.tools.mask import grid_mask_2d
from mbuild.tools.solvent import solvent_box


class Bilayer(Compound):
    """Create a lipid bilayer and add solvent above and below. """
    def __init__(self, lipids, ref_atoms, n_lipids_x=10, n_lipids_y=10, 
                 area_per_lipid=1.0, solvent=None, lipid_box=None, 
                 spacing_z=0.5, solvent_per_lipid=None, n_solvent=None,
                 random_seed=12345):
        """
        Note: This example is still pretty immature as it only represents some
        brief scratch work. Feel free to flesh it out!

        Args:
            lipids (list): list of tuples in format (lipid, frac) where frac is the
            fraction of that lipid in the bilayer (lipid is a Compound)
            ref_atoms (int): Indices of the atom in lipids to form the interface (one for
            each lipid in lipids) (i.e., this atom is shifted to the 'interface' level)
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
            random_seed: seed for random number generator for filling in lipids
        """
        super(Bilayer, self).__init__()

        #TODO: right now this doesn't check to make sure the fractions of each lipid
        # adds up to 1.0. The last value is just what's left after all of the other ones
        # have been added. At the very least, this should check to make sure the given
        # fractions don't add up to > 1.
        frac_sum = 0
        for lipid in lipids:
            frac_sum += lipid[1]
        assert frac_sum <= 1.0
        if frac_sum != 1.0:
            new = 1.0 - frac_sum + lipids[-1][1]
            warn_message = 'Warning: lipid fractions do not add up to 1.0. '
            warn_message += 'Adjusting fraction of lipids[-1] to %f.' % new
            warnings.warn(warn_message)
        # still, if frac_sum < 1.0, then lipids[-1] is not going to be in the same
        # abundance specified by the input
        

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

        # parse lipids to figure out how many of each lipid to add
        # TODO: make simpler for single-component bilayers
        assert len(ref_atoms) == len(lipids)
        n_lipids_per_layer = n_lipids_x * n_lipids_y
        n_each_lipid_per_layer = []
        for lipid in lipids[:-1]:
            n_each_lipid_per_layer.append(int(round(lipid[1] * n_lipids_per_layer)))
        n_each_lipid_per_layer.append(
                n_lipids_per_layer - sum(n_each_lipid_per_layer))
        assert len(n_each_lipid_per_layer) == len(lipids)

        # create points and stretch to fit in box
        mask = grid_mask_2d(n_lipids_x, n_lipids_y)
        mask *= np.sqrt(area_per_lipid * n_lipids_x * n_lipids_y)

        # now give a lipid label to each point
        # TODO: this seems unnecessarily ugly, make more elegant
        labels = np.zeros([mask.shape[0]], dtype=np.uint8)
        seed(random_seed)
        while not labels.all():
            which_point = randint(0, labels.shape[0]-1)   # random index of labels
            if not labels[which_point]:   # if it hasn't been specified
                which_lipid = randint(0, len(n_each_lipid_per_layer)-1)
                if n_each_lipid_per_layer[which_lipid] > 0:
                    labels[which_point] = which_lipid + 1   # cannot be 0 for np.all
                    n_each_lipid_per_layer[which_lipid] -= 1
                else:
                    labels[which_point] = 0

        spacing = np.array([0, 0, spacing_z])
        for i, point in enumerate(mask):
            top_lipid = deepcopy(lipids[labels[i]-1][0])
            translate(top_lipid, -top_lipid.atoms[ref_atoms[labels[i]-1]] + spacing)
            translate(top_lipid, point)
            self.add(top_lipid)

            bot_lipid = deepcopy(lipids[labels[i]-1][0])
            translate(bot_lipid, -bot_lipid.atoms[ref_atoms[labels[i]-1]] + spacing)
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
    chol = Trajectory.load(get_fn('cg-chol.hoomdxml'))
    chol = chol.to_compound()
    rotate_around_x(chol, -135.0*np.pi/180)
    rotate_around_y(chol, -45.0*np.pi/180)
    lipids = [(ecerns, 0.5), (chol, 0.5)] 

    bilayer = Bilayer(lipids, n_lipids_x=10, n_lipids_y=10,
            area_per_lipid=1.4, solvent=water, ref_atoms=[0, 6], 
            spacing_z=0.5, solvent_per_lipid=10)

    #test_water_box = Box(mins=[-5, 6, -13], maxs=[15, 9, -8])
    #test_water = solvent_box(water, test_water_box)
    #test_water.save(filename='water_test.hoomdxml')

    bilayer = bilayer.to_trajectory()
    bilayer.topology.load_ff_bonds()
    bilayer.save(filename='bilayer.hoomdxml')

if __name__ == "__main__":
    main()
