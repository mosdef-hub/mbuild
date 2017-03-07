# -*- coding: utf-8 -*-


# -- ==bilayer== --
from random import seed, shuffle

import numpy as np

import mbuild as mb
from mbuild import clone

class Bilayer(mb.Compound):
    """Create a lipid bilayer and add solvent above and below.

     This example is still pretty immature as it only represents some brief
     scratch work. Feel free to flesh it out!

    Attributes
    ----------
    lipids : list
        List of tuples in format (lipid, frac) where frac is the fraction of
        that lipid in the bilayer (lipid is a Compound)
    ref_atoms : int
        Indices of the atom in lipids to form the interface, one for each lipid
        in lipids (i.e., this atom is shifted to the 'interface' level)
    n_lipids_x : int
        Number of lipids in the x-direction per layer.
    n_lipids_y : int
        Number of lipids in the x-direction per layer.
    area_per_lipid : float
        Area per lipid.
    solvent : Compound
        Compound to solvate the bilayer with. Typically, a pre-equilibrated box
        of solvent.
    lipid_box : Box, optional
        A Box containing the lipids where no solvent will be added.
    spacing_z : float, optional
        Amount of space to add between opposing monolayers.
    solvent_per_lipid : int, optional, default=
        Number of solvent molecules per lipid
    n_solvent : int, optional, default=None
        *Total* number of solvent molecules
    random_seed : int, optional, default=12345
        Seed for random number generator for filling in lipids.
    mirror : bool, optional, default=True
        Make top and bottom layers mirrors of each other.

    """
    def __init__(self, lipids, ref_atoms, n_lipids_x=10, n_lipids_y=10, 
                 area_per_lipid=1.0, solvent=None, lipid_box=None, 
                 spacing_z=0.5, solvent_per_lipid=None, n_solvent=None,
                 random_seed=12345, mirror=True):
        super(Bilayer, self).__init__()

        # Santitize inputs.
        if sum([lipid[1] for lipid in lipids]) != 1.0:
            raise ValueError('Lipid fractions do not add up to 1.')
        assert len(ref_atoms) == len(lipids)

        self.lipids = lipids
        self.ref_atoms = ref_atoms
        self._lipid_box = lipid_box

        # 2D Lipid locations.
        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y
        self.apl = area_per_lipid
        self.n_lipids_per_layer = self.n_lipids_x * self.n_lipids_y
        self.pattern = mb.Grid2DPattern(n_lipids_x, n_lipids_y)
        self.pattern.scale(np.sqrt(self.apl * self.n_lipids_per_layer))

        # Solvent info.
        self.solvent = solvent
        self.n_solvent = n_solvent
        self.solvent_per_lipid = solvent_per_lipid

        # Other inputs.
        self.spacing = np.array([0, 0, spacing_z])
        self.random_seed = random_seed
        self.mirror = mirror

        self._number_of_each_lipid_per_layer = []
        self._solvent_per_layer = None

        # Containers for lipids and solvent.
        self.lipid_components = mb.Compound()
        self.solvent_components = mb.Compound()

        # Assemble the lipid layers
        seed(self.random_seed)
        top_layer, top_lipid_labels = self.create_layer()
        self.lipid_components.add(top_layer)
        if self.mirror == True:
            bottom_layer, bottom_lipid_labels = self.create_layer(lipid_indices=top_lipid_labels, flip_orientation=True)
        else:
            bottom_layer, bottom_lipid_labels = self.create_layer(flip_orientation=True)
        self.lipid_components.add(bottom_layer)

        # solvate the lipids
        #self.solvate_bilayer()  # TODO: needs fixing

        # add everything to the big list
        self.add(self.lipid_components)
        self.add(self.solvent_components)
        print(self.number_of_each_lipid_per_layer)
        # TODO(tim): shift everything so that the lipids are centered in the box?

    def create_layer(self, lipid_indices=None, flip_orientation=False):
        """Create a monolayer of lipids.

        Parameters
        ----------
        lipid_indices : list, optional, default=None
            A list of indices associated with each lipid in the layer.
        flip_orientation : bool, optional, default=False
            Flip the orientation of the layer with respect to the z-dimension.

        """
        layer = mb.Compound()
        if not lipid_indices:
            lipid_indices = list(range(self.n_lipids_per_layer))
            shuffle(lipid_indices)

        for n_type, n_of_lipid_type in enumerate(self.number_of_each_lipid_per_layer):
            current_type = self.lipids[n_type][0]
            for n_this_type, n_this_lipid_type in enumerate(range(n_of_lipid_type)):
                lipids_placed = n_type + n_this_type
                new_lipid = clone(current_type)
                random_index = lipid_indices[lipids_placed]
                position = self.pattern[random_index]

                # Zero and space in z-direction
                particles = list(new_lipid.particles())
                ref_atom = self.ref_atoms[n_type]
                new_lipid.translate(-particles[ref_atom].pos + self.spacing)

                # Move to point on pattern
                if flip_orientation == True:
                    center = new_lipid.center
                    center[2] = 0.0
                    new_lipid.translate(-center)
                    new_lipid.rotate(np.pi, [1, 0, 0])
                    new_lipid.translate(center)
                new_lipid.translate(position)
                layer.add(new_lipid)
        return layer, lipid_indices

    def solvate_bilayer(self):
        """Solvate the constructed bilayer. """
        solvent_number_density = self.solvent.n_particles / np.prod(self.solvent.periodicity)

        lengths = self.lipid_box.lengths
        water_box_z = self.solvent_per_layer / (lengths[0] * lengths[1] * solvent_number_density)

        mins = self.lipid_box.mins
        maxs = self.lipid_box.maxs
        bilayer_solvent_box = mb.Box(mins=[mins[0], mins[1], maxs[2]],
                                     maxs=[maxs[0], maxs[1], maxs[2] + 2 * water_box_z])

        self.solvent_components.add(mb.fill_box(self.solvent, bilayer_solvent_box))

    @property
    def solvent_per_layer(self):
        """Determine the number of solvent molecules per single layer.  """
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

    @property
    def number_of_each_lipid_per_layer(self):
        """The number of each lipid per layer. """
        if self._number_of_each_lipid_per_layer:
            return self._number_of_each_lipid_per_layer

        for lipid in self.lipids[:-1]:
            self._number_of_each_lipid_per_layer.append(int(round(lipid[1] * self.n_lipids_per_layer)))

        # TODO: give warning if frac * n different than actual
        # Rounding errors may make this off by 1, so just do total - whats_been_added.
        self._number_of_each_lipid_per_layer.append(self.n_lipids_per_layer - sum(self._number_of_each_lipid_per_layer))
        assert len(self._number_of_each_lipid_per_layer) == len(self.lipids)
        return self._number_of_each_lipid_per_layer

    @property
    def lipid_box(self):
        """The box containing all of the lipids. """
        if self._lipid_box:
            return self._lipid_box
        else:
            self._lipid_box = self.lipid_components.boundingbox
            # Add buffer around lipid box.
            self._lipid_box.mins -= np.array([0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl)])
            self._lipid_box.maxs += np.array([0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl)])
            return self._lipid_box

# -- ==bilayer== --
