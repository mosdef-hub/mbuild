import warnings
from copy import deepcopy
from random import seed, randint, shuffle

import numpy as np

from mbuild.box import Box
from mbuild.compound import Compound
from mbuild.coordinate_transform import translate, rotate_around_x, rotate_around_y
from mbuild.coordinate_transform import rotate_around_z
from mbuild.tools.mask import grid_mask_2d
from mbuild.tools.solvent import solvent_box


class Bilayer(Compound):
    """Create a lipid bilayer and add solvent above and below. """
    def __init__(self, lipids, ref_atoms, n_lipids_x=10, n_lipids_y=10, 
                 area_per_lipid=1.0, solvent=None, lipid_box=None, 
                 spacing_z=0.5, solvent_per_lipid=None, n_solvent=None,
                 random_seed=12345, mirror=True):
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
            mirror (bool): make top and bottom layers mirrors of each other if True
        """
        super(Bilayer, self).__init__()

        # set constants from constructor call
        self.lipids = lipids
        self.solvent = solvent
        self.apl = area_per_lipid
        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y
        self.ref_atoms = ref_atoms
        self.spacing_z = spacing_z
        self.mirror = mirror
        self.random_seed = random_seed
        self.n_solvent = n_solvent
        self.solvent_per_lipid = solvent_per_lipid
        self._lipid_box = lipid_box

        # why do i have to initialize these 2?
        self._n_each_lipid_per_layer = []
        self._solvent_per_layer = None

        # a few calculations to figure things out
        self.n_lipids_per_layer = self.n_lipids_x * self.n_lipids_y
        self.spacing = np.array([0, 0, spacing_z])
        self.mask = grid_mask_2d(n_lipids_x, n_lipids_y)   # lipid locations
        self.mask *= np.sqrt(self.apl * self.n_lipids_per_layer)

        # safety checks
        self.check_fractions()
        self.check_ref_atoms()

        # containers for lipids and solvent
        self.lipid_components = Compound()
        self.solvent_components = Compound()

        # assemble the lipid layers
        # TODO(tim): random number seed here?
        seed(self.random_seed)   
        top_layer, top_lipid_labels = self.create_layer()
        self.lipid_components.add(top_layer)
        if self.mirror == True:
            bottom_layer, bottom_lipid_labels = self.create_layer(
                    lipid_labels=top_lipid_labels,
                    flip_orientation=True)
        else:
            bottom_layer, bottom_lipid_labels = self.create_layer(
                    flip_orientation=True)
        self.lipid_components.add(bottom_layer)

        # solvate the lipids
        self.solvate_bilayer()

        # add everything to the big list
        self.add(self.lipid_components)
        self.add(self.solvent_components)
        # TODO(tim): shift everything so that the lipids are centered in the box?

    @property 
    def solvent_per_layer(self):
        """
        Figure out the number of solvent molecules per single layer.

        """
        if self._solvent_per_layer:
            return self._solvent_per_layer

        assert not (self.solvent_per_lipid is None and self.n_solvent is None)
        if self.solvent_per_lipid is not None:
            assert self.n_solvent is None
            self._solvent_per_layer = self.n_lipids_per_layer * self.solvent_per_lipid
        elif self.n_solvent is not None:
            assert self.solvent_per_lipid is None
            self._solvent_per_layer = self.n_solvent / 2
        return self._solvent_per_layer

    def check_fractions(self):
        frac_sum = 0
        for lipid in self.lipids:
            frac_sum += lipid[1]
        assert frac_sum == 1.0, 'Bilayer builder error: Lipid fractions do not add up to 1.'

    def check_ref_atoms(self):
        assert len(self.ref_atoms) == len(self.lipids)

    @property
    def n_each_lipid_per_layer(self):
        import pdb
        if self._n_each_lipid_per_layer:
            return self._n_each_lipid_per_layer

        self._n_each_lipid_per_layer = []
        for lipid in self.lipids[:-1]:
            self._n_each_lipid_per_layer.append(
                    int(round(lipid[1] * self.n_lipids_per_layer)))
        # TODO: give warning if frac*n different than actual
        # rounding errors may make this off by 1, so just do total - whats_been_added
        self._n_each_lipid_per_layer.append(
                self.n_lipids_per_layer - sum(self._n_each_lipid_per_layer))
        assert len(self._n_each_lipid_per_layer) == len(self.lipids)
        return self._n_each_lipid_per_layer

    def create_layer(self, lipid_labels=None, flip_orientation=False):
        """
        Args:
            top (bool): Top (no rotation) or bottom (rotate about x) layer
        """
        layer = Compound()
        if not lipid_labels:
            lipid_labels = range(self.n_lipids_per_layer)
            shuffle(lipid_labels)
        lipids_placed = 0
        for i, n_of_lipid_type in enumerate(self.n_each_lipid_per_layer):
            current_type = self.lipids[i][0]
            for n_this_lipid_type in range(n_of_lipid_type):
                new_lipid = deepcopy(current_type)
                random_index = lipid_labels[lipids_placed]
                position = self.mask[random_index]

                # Zero and space in z-direction
                translate(
                        new_lipid,
                        -new_lipid.atoms[self.ref_atoms[i]] + self.spacing)
                # Move to point on mask
                if flip_orientation == True:
                    t_com = new_lipid.center_of_mass()
                    t_com[2] = 0.0
                    # TODO(tim): function for this? 
                    # e.g., rotate_around_x_keep_com(compound, bool(3))
                    translate(new_lipid, -t_com)
                    rotate_around_x(new_lipid, np.pi)
                    translate(new_lipid, t_com)
                translate(new_lipid, position)
                layer.add(new_lipid)
                lipids_placed += 1
        return layer, lipid_labels

    @property
    def lipid_box(self):
        if self._lipid_box:
            return self._lipid_box
        else:
            self._lipid_box = self.lipid_components.boundingbox()
            # add buffer around lipid box
            self._lipid_box.mins -= np.array(
                    [0.5*np.sqrt(self.apl),
                     0.5*np.sqrt(self.apl),
                     0.5*np.sqrt(self.apl)])
            self._lipid_box.maxs += np.array(
                    [0.5*np.sqrt(self.apl),
                     0.5*np.sqrt(self.apl),
                     0.5*np.sqrt(self.apl)])
            return self._lipid_box

    def solvate_bilayer(self):
        solvent_number_density = self.solvent.n_atoms / np.prod(self.solvent.periodicity)
        import pdb
        water_box_z = self.solvent_per_layer / (self.lipid_box.lengths[0]
                * self.lipid_box.lengths[1] * solvent_number_density)

        bilayer_solvent_box = Box(mins=[self.lipid_box.mins[0],
                                self.lipid_box.mins[1], 
                                self.lipid_box.maxs[2]],
                              maxs=[self.lipid_box.maxs[0],
                                    self.lipid_box.maxs[1], 
                                    self.lipid_box.maxs[2] + 2 * water_box_z])
        self.solvent_components.add(
                solvent_box(self.solvent, bilayer_solvent_box))

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

    import cProfile, pstats, StringIO
    #pr=cProfile.Profile()
    #pr.enable()
    bilayer = Bilayer(lipids, n_lipids_x=15, n_lipids_y=15,
            area_per_lipid=1.4, solvent=water, ref_atoms=[1, 6], 
            spacing_z=0.7, solvent_per_lipid=20, mirror=False)
    """
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'tottime'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()
    """

    #test_water_box = Box(mins=[-5, 6, -13], maxs=[15, 9, -8])
    #test_water = solvent_box(water, test_water_box)
    #test_water.save(filename='water_test.hoomdxml')

    bilayer = bilayer.to_trajectory()
    bilayer.topology.load_ff_bonds()
    bilayer.save(filename='bilayer.hoomdxml')
    import os
    os.system('vmd_MACOSXX86 -e vis.vmd')

if __name__ == "__main__":
    main()
