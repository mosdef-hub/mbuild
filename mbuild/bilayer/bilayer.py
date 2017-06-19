from random import shuffle, random, uniform

import numpy as np

import mbuild as mb
from mbuild import clone
from mbuild.lib.moieties import H2O
from mbuild.prototypes import DSPC, ALC, FFA, ISIS, HDHD
from mbuild.UA_molecules import DSPCUA, DMPCUA, DPPCUA, FFAUA, ISISUA, HDHDUA, ALCUA
from mbuild.lib.cg_molecules import ECer2, UCer2, Chol, FFAC16, FFAC20, FFAC24, Water

class Bilayer(mb.Compound):
    """The Bilayer Builder creates a lipid bilayer, solvates it, and stores it as an mBuild Compound. 
    Because the Bilayer Builder does not check for the orientation of the mBuild molecules the user 
    provides from a previously created class or file, the user must ensure the following things are true 
    about the lipids they input to the builder:
        - The lipid's longitudinal axis lies on the z-axis of the Cartesian coordinate system
        - The lipid's reference atom is located at the origin (see ref_atoms below)
        - The rest of the lipid is pointing in the negative z direction
        
    The user may input the fraction of each lipid in the bilayer, as well as a number of important bilayer 
    properties such as area per lipid and tilt angle.

    Parameters
    ----------
    lipids : list
        List of tuples in format (lipid, frac, z-offset) where frac is the fraction of
        that lipid in the bilayer (lipid is a Compound), amd z-offset is the distance in nanometers
        from the headgroup to the lipid-water interface
    ref_atoms : int
        Indices of the atom in lipids to form the solvent interface. By definition, this atom denotes
        the headgroup of the lipid.
    n_lipids_x : int
        Number of lipids in the x-direction per layer.
    n_lipids_y : int
        Number of lipids in the y-direction per layer.
    area_per_lipid : float, optional, default=0.3
        Area per lipid in nanometers squared.
    tilt_angle : float, optional, default=0.0
        Tilt angle of each lipid with respect to the y-axis, in degrees.
    max_tail_randomization : float, optional, default=0.0
        For a set tilt angle, the maximum amount each lipid may be randomly spun around the z-axis.
    solvent : Compound, optional, default=H2O()
        Compound to solvate the bilayer with. Typically, a pre-equilibrated box
        of solvent.
    lipid_box : Box, optional
        A Box containing the lipids where no solvent will be added.
    spacing_z : float, optional, default=0.2
        Amount of space to add between opposing monolayers, in nanometers.
    solvent_per_lipid : int, optional, default=
        Number of solvent molecules per lipid
    n_solvent : int, optional, default=None
        *Total* number of solvent molecules
    mirror : bool, optional, default=True
        Make top and bottom layers mirrors of each other. 
    cross_tilt : bool, optional, default=False
        Gives the bilayer a cross-tilting property
    unit_conversion : int, optional, default=1
        Adjust for unit conversions if necessary
        
    Attributes
    ----------
    n_lipids_x : int
        Number of lipids in the x-direction per layer.
    n_lipids_y : int
        Number of lipids in the y-direction per layer.
    lipids : list
        List of tuples in format (lipid, frac, z-offset) where frac is the fraction of
        that lipid in the bilayer (lipid is a Compound), amd z-offset is the distance in nanometers
        from the headgroup to the lipid-water interface
    ref_atoms : int
        Indices of the atom in lipids to form the solvent interface
    lipid_box : Box, optional
        A Box containing the lipids where no solvent will be added.
    solvent : Compound
        Compound to solvate the bilayer with. Typically, a pre-equilibrated box
        of solvent.
    solvent_per_lipid : int, optional, default=
        Number of solvent molecules per lipid
    unit_conversion : int, optional, default=1
        Adjust for unit conversions if necessary
    """
    
    def __init__(self, lipids, ref_atoms, n_lipids_x=5, n_lipids_y=5, 
                 area_per_lipid=0.3, tilt_angle = 0.0, 
                 max_tail_randomization = 0, solvent=H2O(), lipid_box=None,
                 spacing_z=0.4, solvent_per_lipid=None, mirror=True, 
                 cross_tilt=False, unit_conversion=1):
        super(Bilayer, self).__init__()

        self._sanitize_inputs(lipids, ref_atoms)

        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y
        self.lipids = lipids
        self.ref_atoms = ref_atoms
        self._lipid_box = lipid_box
        self.unit_conversion = unit_conversion

        # Set the 2D lipid locations on a grid pattern
        self._set_grid_pattern(area_per_lipid)
        
        #Set other important geometric attributes of the lipids
        self.tilt = tilt_angle * np.pi / 180
        self.z_orientation = max_tail_randomization

        # Solvent attributes.
        self.solvent = solvent
        self.solvent_per_lipid = solvent_per_lipid
        
        # Private attributes for getter methods 
        self._number_of_each_lipid_per_layer = []
        self._solvent_per_layer = None
        
        # Create lipid leaflets and add them to the bilayer Compound structure.
        lipid_bilayer = mb.Compound()
        solvent_components = mb.Compound()
        top_leaflet, top_lipid_labels = self.create_layer()
        lipid_bilayer.add(top_leaflet, label = 'top_leaflet')
        if mirror:
            bottom_leaflet, bottom_lipid_labels = self.create_layer(lipid_indices=top_lipid_labels)
        else:
            bottom_leaflet, bottom_lipid_labels = self.create_layer()
        lipid_bilayer.add(bottom_leaflet, label = 'bottom_leaflet')
        
        # Translate and rotate the bottom leaflet
        self._translate_bottom_leaflet(lipid_bilayer, spacing_z, cross_tilt)
        
        # solvate the lipids
        solvent_components = self.solvate_bilayer(lipid_bilayer, 
                solvent_components)

        # add everything to the overall Compound structure
        self.add(lipid_bilayer, label = 'lipid_bilayer')
        self.add(solvent_components, label='solvent')
        
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
            lipid_indices = list(range(self.lipids_per_leaflet))
            shuffle(lipid_indices)
            
        n_previous_types = 0
        for n_type, n_of_lipid_type in enumerate(self.number_of_each_lipid_per_leaflet):
            current_type = self.lipids[n_type][0]
            for n_this_type in range(n_of_lipid_type):
                lipids_placed = n_previous_types + n_this_type
                new_lipid = clone(current_type)
                random_index = lipid_indices[lipids_placed]
                position = self.pattern[random_index]

                # Zero and space in z-direction, apply necessary geometric transformations
                particles = list(new_lipid.particles())
                ref_atom = self.ref_atoms[n_type]
                new_lipid.spin(-self.tilt, [0, 1, 0])
                new_lipid.spin(theta = uniform(-self.z_orientation, 
                    self.z_orientation) * np.pi/180, around = [0, 0, 1])
                new_lipid.translate(-particles[ref_atom].pos  + 
                        [0, 0, self.lipids[n_type][2]])
                
                # Move to point on pattern and add the lipid to the leaflet Compound
                if flip_orientation == True:
                    new_lipid.spin(np.pi, [1, 0, 0]) 
                new_lipid.translate(position)
                layer.add(new_lipid)
                
            n_previous_types += n_of_lipid_type
        return layer, lipid_indices
    
    def _translate_bottom_leaflet(self, lipid_bilayer, spacing_z, 
            cross_tilt):
        mins = np.amin(lipid_bilayer['top_leaflet'].xyz, axis=0)
        z_min = mins[2]
        spacing_z *= self.unit_conversion
        lipid_bilayer['bottom_leaflet'].translate([0 , 0, 
            (z_min - spacing_z)])
        lipid_bilayer['bottom_leaflet'].spin(np.pi, [1, 0, 0])
        if not cross_tilt:
            lipid_bilayer['bottom_leaflet'].spin(np.pi, [0, 0, 1])

    def solvate_bilayer(self, lipid_bilayer, solvent_components):
        """Solvate the constructed bilayer. 
        
        Parameters
        ----------
        lipid_bilayer : mb.Compound
            An mBuild compound containing both leaflets of the bilayer
        solvent_components : mb.Compound
            A (temporarily) empty Compound that will serve as the container
            for the solvated boxes
        
        Returns
        -------
        solvent_components : mb.Compound
            The container for the solvated boxes created in this function
        """
        water_volume = .0299 * np.power(self.unit_conversion, 3) 
        #TODO: Support other solvents' volumes in a user-friendly way
        
        #Calculate box dimension
        box_volume = water_volume * self.solvent_per_layer
        x_min_box = min(lipid_bilayer.xyz[: , 0])
        x_max_box = max(lipid_bilayer.xyz[: , 0])
        y_min_box = min(lipid_bilayer.xyz[: , 1])
        y_max_box = max(lipid_bilayer.xyz[: , 1])
        
        box_area = (x_max_box - x_min_box) * (y_max_box - y_min_box)
        box_height = box_volume / box_area
        
        z_min_top = max(lipid_bilayer.xyz[: , 2])
        z_max_top = z_min_top + box_height
        
        z_max_bottom = min(lipid_bilayer.xyz[: , 2])
        z_min_bottom = z_max_bottom - box_height

        
        top_solvent_box = mb.Box(mins=[x_min_box, y_min_box, z_min_top],
                                     maxs=[x_max_box, y_max_box, z_max_top])
        
        bottom_solvent_box = mb.Box(mins=[x_min_box, y_min_box, 
            z_min_bottom], 
            maxs=[x_max_box, y_max_box, z_max_bottom])

        solvent_components.add(mb.fill_region(compound = [self.solvent, 
            self.solvent],
            n_compounds = [self.solvent_per_layer, self.solvent_per_layer],
            region=[top_solvent_box, bottom_solvent_box]))
        
        return solvent_components
            
    @property
    def solvent_per_layer(self):
        """Determine the number of solvent molecules per single layer.  """
        return self.n_lipids_x * self.n_lipids_y * self.solvent_per_lipid

    @property
    def number_of_each_lipid_per_leaflet(self):
        """The number of each lipid per leaflet. """
        if self._number_of_each_lipid_per_layer:
            return self._number_of_each_lipid_per_layer

        for lipid in self.lipids[:-1]:
            self._number_of_each_lipid_per_layer.append(int(round(lipid[1] * self.lipids_per_leaflet)))

        # TODO: give warning if frac * n different than actual
        # Rounding errors may make this off by 1, so just do total - whats_been_added.
        self._number_of_each_lipid_per_layer.append(self.lipids_per_leaflet 
                - sum(self._number_of_each_lipid_per_layer))
        assert len(self._number_of_each_lipid_per_layer) == len(self.lipids)
        return self._number_of_each_lipid_per_layer

    @property
    def lipid_box(self):
        """The box containing all of the lipids. """
        if self._lipid_box:
            return self._lipid_box
        else:
            self._lipid_box = self.lipid_bilayer.boundingbox
            # Add buffer around lipid box.
            self._lipid_box.mins -= np.array([0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl)])
            self._lipid_box.maxs += np.array([0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl),
                                              0.5*np.sqrt(self.apl)])
            return self._lipid_box
        
    def _sanitize_inputs(self, lipids, ref_atoms):
        """Check for proper inputs
    
        Ensure that the user's lipid fractions add up to 1, 
        or raise a ValueError.
    
        Ensure that the user has input the same number of reference
        atoms and lipid components, or raise an AssertionError."""
    
        if sum([lipid[1] for lipid in lipids]) != 1.0:
            raise ValueError('Lipid fractions do not add up to 1.')
        assert len(ref_atoms) == len(lipids), "Please provide one reference atom for each lipid"
        
    def _set_grid_pattern(self, area_per_lipid):
        """Utilize an mBuild 2DGridPattern to create the scaffold of points that the lipids will be laid onto"""
        
        self.lipids_per_leaflet = self.n_lipids_x * self.n_lipids_y
        pattern = mb.Grid2DPattern(self.n_lipids_x, self.n_lipids_y)
        pattern.scale(np.sqrt(area_per_lipid * self.lipids_per_leaflet) *
                self.unit_conversion)
        
        self.pattern = pattern
    
    
if __name__ == '__main__':
    lipids = [(DMPCUA(), 0.5, 0), (FFAUA(16, hcap = True), 0.25, -.25), 
            (ISISUA(), 0.25, -.2)]
    bilayer = Bilayer(lipids, n_lipids_x = 8, n_lipids_y = 8, 
            ref_atoms=[0,2,1], tilt_angle = 20, max_tail_randomization =15,
            solvent_per_lipid=20) 
    bilayer.save('bilayer_test.mol2', overwrite=True)
