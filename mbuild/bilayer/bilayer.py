from random import shuffle, uniform
import numpy as np
import math
import argparse

import mbuild as mb
from mbuild import clone
from mbuild.prototypes import DSPC, ALC, FFA, ISIS, HDHD, H2O
from mbuild.UA_molecules import DSPCUA, DMPCUA, DPPCUA, FFAUA, ISISUA, ALCUA
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
    properties such as area per lipid, bilayer size, and tilt angle.

    Parameters
    ----------
    lipids : list of tuple
        List of tuples in format (lipid, frac, z-offset) where frac is
        the fraction of that lipid in the bilayer (lipid is a
        Compound), and z-offset is the distance in nanometers from the
        headgroup to the lipid-water interface
    args : NameSpace
        The parser passed by the argparse module; allows for command-line integration with
        the Bilayer builder.
    max_tail_randomization : float, optional, default=0.0
        For a set tilt angle, the maximum amount each lipid may be
        randomly spun around the z-axis.
    solvent : Compound, optional, default=H2O()
        Compound to solvate the bilayer with. Typically, a
        pre-equilibrated box of solvent.
    lipid_box : Box, optional
        A Box containing the lipids where no solvent will be added.
    unit_conversion : int, optional, default=1
        Adjust for unit conversions if necessary
        
    Attributes
    ----------
    n_lipids_x : int
        Number of lipids in the x-direction per layer.
    n_lipids_y : int
        Number of lipids in the y-direction per layer.
    lipids : list
        List of tuples in format (lipid, frac, z-offset, ref_atom) where frac is
        the fraction of that lipid in the bilayer (lipid is a
        Compound), and z-offset is the distance in nanometers from the
        headgroup to the lipid-water interface
    ref_atoms : int
        Indices of the atom in lipids to form the solvent interface
    lipid_box : Box, optional
        A Box containing the lipids where no solvent will be added.
    solvent : Compound
        Compound to solvate the bilayer with. Typically, a
        pre-equilibrated box
        of solvent.
    solvent_per_lipid : int, optional, default=0
        Number of solvent molecules per lipid
    unit_conversion : int, optional, default=1
        Adjust for unit conversions if necessary
    """
    
    def __init__(self, lipids, args, itp_path, max_tail_randomization=0, solvent=H2O(),
                 lipid_box=None, unit_conversion=1, filename=None):
        super(Bilayer, self).__init__()

        self.n_lipids_x = args.n_lipids_x
        self.n_lipids_y = args.n_lipids_y
        self.lipids = lipids
        if args.apl is None:
            area_per_lipid = uniform(0.25 + (0.3 * (args.DSPC + args.DPPC + args.DMPC)),
                                     0.35 + (0.3 * (args.DSPC + args.DPPC + args.DMPC)))
        else:
            area_per_lipid = args.apl

        self._sanitize_inputs()

        self.ref_atoms = []
        for lipid in lipids:
            self.ref_atoms.append(lipid[3])
        assert len(self.ref_atoms) == len(lipids)
        self._lipid_box = lipid_box
        self.unit_conversion = unit_conversion
        if filename:
            self.filename = filename
        else:
            self.filename = ''
            components = [str(lipid[0].name) + '_' + str(lipid[1]) + '_' for lipid in lipids]
            for component in components:
                self.filename += component
            self.filename += 'UA'
            # self.filename += 'AA' TODO: Make a command line flag
            # self.filename += 'CG' TODO: Same as above

        # Set the 2D lipid locations on a grid pattern
        self._set_grid_pattern(area_per_lipid)
        z_spacing = args.spacing

        # Set other important geometric attributes of the lipids
        self.tilt = args.tilt_angle * np.pi / 180

        self.z_orientation = max_tail_randomization

        # Solvent attributes.
        self.solvent = solvent
        self.solvent_per_lipid = args.solvate
        
        # Private attributes for getter methods 
        self._number_of_each_lipid_per_layer = []
        self._solvent_per_layer = None
        
        # Path to .itp files
        self.itp_path = itp_path

        # Write topology file header
        print('Writing <{0}> ...'.format(self.filename + '.top'))
        top_file = self.write_top_header(self.filename)

        # Create lipid leaflets and add them to the bilayer compound 
        # structure.
        lipid_bilayer = mb.Compound()
        solvent_components = mb.Compound()
        top_file, top_leaflet, top_lipid_labels = self.create_layer(top_file)
        lipid_bilayer.add(top_leaflet, label='top_leaflet')
        if args.mirror:
            top_file, bottom_leaflet, bottom_lipid_labels = self.create_layer(top_file, lipid_indices=top_lipid_labels)
        else:
            top_file, bottom_leaflet, bottom_lipid_labels = self.create_layer(top_file)
        lipid_bilayer.add(bottom_leaflet, label='bottom_leaflet')
        
        # Translate and rotate the bottom leaflet
        self._translate_bottom_leaflet(lipid_bilayer, z_spacing, args)

        # Check to make sure leaflets have been spaced correctly
        top_z_min = np.amin(lipid_bilayer['top_leaflet'].xyz, axis=0)[2]
        bottom_z_max = np.amax(lipid_bilayer['bottom_leaflet'].xyz, axis=0)[2]
        space = top_z_min - bottom_z_max
        assert math.isclose(space, args.spacing), \
            "Desired spacing: {}, Actual spacing: {}".format(args.spacing, space)
        
        # solvate the lipids
        if self.solvent_per_lipid > 0:
            top_file, solvent_components = self.solvate_bilayer(top_file,
                                                                lipid_bilayer, solvent_components)

        # recenter the bottom leaflet over the bottom water box
        self._post_translate(solvent_components, lipid_bilayer)

        # close the completed topology file
        top_file.close()

        # add everything to the overall Compound structure
        self.add(lipid_bilayer, label='lipid_bilayer')
        if self.solvent_per_lipid > 0:
            self.add(solvent_components, label='solvent')

        # adjust the periodicity of the system, center the system in the new box
        water_box = (np.amax(self['solvent'].xyz, axis=0) - np.amin(self['solvent'].xyz, axis=0)) + 0.5
        self.periodicity = water_box
        distance = -self['solvent'].xyz.min(axis=0) + 0.25
        self.translate(distance)  # Bring the bilayer inside the bounding box for the .gro file

    def create_layer(self, top_file, lipid_indices=None):
        """Create a monolayer of lipids.

        Parameters
        ----------
        top_file : file
            A topology file into which molecule information will be
            written for simulation
        lipid_indices : list, optional, default=None
            A list of indices associated with each lipid in the layer.

        Returns
        -------
        layer : mb.Compound
            An mBuild compound containing the entire leaflet of
            molecules
        lipid indices : list
            The indices for each lipid corresponding to its location
            on the grid; necessary if user desires a mirror effect between
            the top and bottom leaflets
        
        """
        layer = mb.Compound()
        if not lipid_indices:
            lipid_indices = list(range(self.lipids_per_leaflet))
            shuffle(lipid_indices)
            
        n_previous_types = 0
        for n_type, n_of_lipid_type in enumerate(
                self.number_of_each_lipid_per_leaflet):
            current_type = self.lipids[n_type][0]
            
            # Write to topology file
            if self.number_of_each_lipid_per_leaflet[n_type] != 0:
                top_file.write("{:<10s}{:<10d}\n".format(
                    self.lipids[n_type][0].name,
                    self.number_of_each_lipid_per_leaflet[n_type]))

            # Loop through the leaflet's entire quantity of one lipid
            # type
            for n_this_type in range(n_of_lipid_type):
                lipids_placed = n_previous_types + n_this_type
                new_lipid = clone(current_type)
                random_index = lipid_indices[lipids_placed]
                position = self.pattern[random_index]

                # Zero and space in z-direction, apply necessary geometric
                # transformations
                particles = list(new_lipid.particles())
                ref_atom = self.ref_atoms[n_type]
                new_lipid.spin(-self.tilt, [0, 1, 0])
                new_lipid.spin(theta=uniform(-self.z_orientation,
                                             self.z_orientation) * np.pi/180, around=[0, 0, 1])
                new_lipid.translate(-particles[ref_atom].pos + [0, 0, self.lipids[n_type][2]])
                
                # Move to point on pattern and add the lipid to the
                # leaflet compound
                new_lipid.translate(position)
                layer.add(new_lipid)
                
            n_previous_types += n_of_lipid_type
        return top_file, layer, lipid_indices
    
    def _translate_bottom_leaflet(self, lipid_bilayer, spacing_z, args):

        spacing_z *= self.unit_conversion
        lipid_bilayer['bottom_leaflet'].spin(np.pi, [1, 0, 0])
        top_z_min = np.amin(lipid_bilayer['top_leaflet'].xyz, axis=0)[2]
        bottom_z_max = np.amax(lipid_bilayer['bottom_leaflet'].xyz, axis=0)[2]
        space = top_z_min - bottom_z_max
        lipid_bilayer['bottom_leaflet'].translate([0, 0, (space - spacing_z)])
        if not args.cross_tilt:
            lipid_bilayer['bottom_leaflet'].spin(np.pi, [0, 0, 1])

    def solvate_bilayer(self, top_file, lipid_bilayer, solvent_components):
        """Solvate the constructed bilayer. 
        
        Parameters
        ----------
        top_file : file
            Topology file into which solvent information is written
        lipid_bilayer : mb.Compound
            An mBuild compound containing both leaflets of the bilayer
        solvent_components : mb.Compound
            A (temporarily) empty Compound that will serve as the container
            for the solvated boxes
        
        Returns
        -------
        top_file : file
            Updated topology file
        solvent_components : mb.Compound
            The container for the solvated boxes created in this function
        """
        water_volume = .04 * np.power(self.unit_conversion, 3)
        # TODO: Support other solvents' volumes in a user-friendly way
        
        # Calculate box dimension
        box_volume = water_volume * self.solvent_per_layer
        x_min_box = min(self.pattern[:, 0]) # lipid_bilayer.xyz[:, 0])
        x_max_box = max(self.pattern[:, 0]) # .xyz[:, 0])
        y_min_box = min(self.pattern[:, 1]) # lipid_bilayer.xyz[:, 1])
        y_max_box = max(self.pattern[:, 1]) # lipid_bilayer.xyz[:, 1])
        
        box_area = (x_max_box - x_min_box) * (y_max_box - y_min_box)
        box_height = box_volume / box_area
        
        z_min_top = max(lipid_bilayer.xyz[:, 2]) + 0.3
        z_max_top = z_min_top + box_height
        
        z_max_bottom = min(lipid_bilayer.xyz[:, 2]) - 0.3
        z_min_bottom = z_max_bottom - box_height

        top_solvent_box = mb.Box(mins=[x_min_box, y_min_box, z_min_top],
                                 maxs=[x_max_box, y_max_box, z_max_top])
        
        bottom_solvent_box = mb.Box(mins=[x_min_box, y_min_box, z_min_bottom],
                                    maxs=[x_max_box, y_max_box, z_max_bottom])

        solvent_components.add(mb.fill_region(compound=[self.solvent, self.solvent],
                                              n_compounds=[self.solvent_per_layer, self.solvent_per_layer],
                                              region=[top_solvent_box, bottom_solvent_box]))

        # Write solvent components to topology file
        top_file.write("{:<10s}{:<10d}\n".format('SOL',
                                                 self.solvent_per_layer * 2))

        return top_file, solvent_components

    def _post_translate(self, solvent_components, lipid_bilayer):
        distance = np.amax(solvent_components.xyz, axis=0)[0] - \
                   np.amax(lipid_bilayer['bottom_leaflet'].xyz, axis=0)[0]
        lipid_bilayer['bottom_leaflet'].translate([distance, 0, 0])
            
    def write_top_header(self, filename='bilayer_test'):
        """ Generate topology file header based on existing .itp files

        Parameters
        ----------
        filename : str, optional, default = default
            Filename for the topology file

        Returns
        -------
        top_file
            The open topology file to be written into by the rest of
            the Bilayer Builder

        NOTE: Function assumes a particular path for itp files; use
        constant ITP_PATH to customize

        """
        # TODO: create flag for user to input a custom .itp path

        top_file = open(filename + '.top', 'w')

        # Include statements to locate .itp files for .top creation
        top_file.write(";#include \"{}ff.itp\"\n".format(self.itp_path))
        top_file.write("#include \"{}ff_b.itp\"\n".format(self.itp_path))
        top_file.write("; SPC water topology \n")
        top_file.write(";#include \"{}spc.itp\"\n".format(self.itp_path))
        top_file.write("#include \"{}spc_b.itp\"\n".format(self.itp_path))
        top_file.write("\n[ system]\n")
        top_file.write("United-atom bilayer system\n")
        # top_file.write("All-atom bilayer system\n")
        # top_file.write("Coarse-grained bilayer system\n")
        # TODO: Create flags for above system options
        top_file.write("\n[ molecules ]\n")

        return top_file

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

    def _sanitize_inputs(self):
        """Check for proper inputs
    
        Ensure that the user's lipid fractions add up to 1, 
        or raise a ValueError.
        """
    
        if sum([lipid[1] for lipid in self.lipids]) != 1.0:
            raise ValueError('Lipid fractions do not add up to 1.')
        if not (self.n_lipids_x > 0 and self.n_lipids_y > 0):
            raise ValueError('Bilayer dimensions must be greater than 0')
        if len(self.lipids) > (self.n_lipids_x * self.n_lipids_y):
            raise ValueError('Number of lipids provided exceeds number of available grid spaces')

    def _set_grid_pattern(self, area_per_lipid):
        """Utilize an mBuild 2DGridPattern to create the scaffold of points that the lipids will be laid onto"""
        
        self.lipids_per_leaflet = self.n_lipids_x * self.n_lipids_y
        pattern = mb.Grid2DPattern(self.n_lipids_x, self.n_lipids_y)
        pattern.scale(np.sqrt(area_per_lipid * self.lipids_per_leaflet) * self.unit_conversion)
        
        self.pattern = pattern
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Welcome to the Bilayer Builder!',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog='For more robust descriptions of variables and '
                                            'functionality of the Bilayer Builder, see the '
                                            'README file.')
    
    parser.add_argument('n_lipids_x', type=int, help='The number of lipids in the x direction')
    parser.add_argument('n_lipids_y', type=int, help='The number of lipids in the y direction')
    frac = parser.add_argument_group('Component Fractions', 'The fraction of the bilayer containing each component')
    frac.add_argument('--DSPC', type=float, default=1.0)
    frac.add_argument('--DPPC', type=float, default=0.0)
    frac.add_argument('--DMPC', type=float, default=0.0)
    frac.add_argument('--acd12', type=float, default=0.0)
    frac.add_argument('--acd14', type=float, default=0.0)
    frac.add_argument('--acd16', type=float, default=0.0)
    frac.add_argument('--acd18', type=float, default=0.0)
    frac.add_argument('--acd20', type=float, default=0.0)
    frac.add_argument('--acd22', type=float, default=0.0)
    frac.add_argument('--acd24', type=float, default=0.0)
    frac.add_argument('--alc12', type=float, default=0.0)
    frac.add_argument('--alc14', type=float, default=0.0)
    frac.add_argument('--alc16', type=float, default=0.0)
    frac.add_argument('--alc18', type=float, default=0.0)
    frac.add_argument('--alc20', type=float, default=0.0)
    frac.add_argument('--alc22', type=float, default=0.0)
    frac.add_argument('--alc24', type=float, default=0.0)
    frac.add_argument('--ISIS', type=float, default=0.0)
    geometry = parser.add_argument_group('Bilayer Geometry', 'Important System-Level Geometric Specifications')
    geometry.add_argument('-a', '--apl', type=float, help='The area per lipid of the bilayer')
    geometry.add_argument('-t', '--tilt_angle', type=float, default=uniform(20, 30),
                          help='The tilt angle of the lipids')
    geometry.add_argument('-z', '--spacing', type=float, default=0.4,
                          help='The spacing in between the two bilayer leaflets')
    geometry.add_argument('-m', '--mirror', action='store_true', default=False,
                          help='Makes the top and bottom leaflets mirror images of one another')
    geometry.add_argument('-c', '--cross_tilt', action='store_true', default=False,
                          help='Induces a cross-tilt in the bilayer')
    parser.add_argument('-s', '--solvate', type=int, default=0, help='Solvates the bilayer with a specified number '
                                                                     ' of solvent molecules per lipid')
    parser.add_argument('-u', '--units', type=float, default=1.0, help='Input any unit conversion that is necessary '
                                                                       'here (mBuild default units are in nanometers)')
    args = parser.parse_args()

    lipids = [(DSPCUA(), args.DSPC, 0, 0),
              (DPPCUA(), args.DPPC, -0.3, 0),
              (DMPCUA(), args.DMPC, -0.5, 0),
              (ALCUA(12), args.alc12, -0.2, 13),
              (FFAUA(12, ester=False), args.acd12, -0.2, 13),
              (ALCUA(14), args.alc14, -0.2, 15),
              (FFAUA(14, ester=False), args.acd14, -0.2, 15),
              (ALCUA(16), args.alc16, -0.4, 17),
              (FFAUA(16, ester=False), args.acd16, -0.4, 17),
              (ALCUA(18), args.alc18, -0.4, 19),
              (FFAUA(18, ester=False), args.acd18, -0.4, 19),
              (ALCUA(20), args.alc20, -0.5, 21),
              (FFAUA(20, ester=False), args.acd20, -0.5, 21),
              (ALCUA(22), args.alc22, -0.5, 23),
              (FFAUA(22, ester=False), args.acd22, -0.5, 23),
              (ALCUA(24), args.alc24, -0.4, 25),
              (FFAUA(24, ester=False), args.acd24, -0.4, 25),
              (ISISUA(), args.ISIS, -0.4, 20)]

    # Remove all lipid tuples whose fractions are 0
    lipids_to_pop = []
    for index, lipid in enumerate(lipids):
        if lipid[1] == 0.0:
            lipids_to_pop.append(index)
    for index in sorted(lipids_to_pop, reverse=True):
        lipids.pop(index)

    bilayer = Bilayer(lipids, args, itp_path="/home/loganguy/builds/setup/FF/gromos53a6/",
                      max_tail_randomization=30)

    print('Writing <{0}.gro> ...'.format(bilayer.filename))
    bilayer.save(bilayer.filename + '.gro',
                 residues=['DSPC', 'FFA12', 'ALC12', 'FFA14', 'ALC14', 'FFA16', 'ALC16',
                           'FFA18', 'ALC18', 'FFA20', 'ALC20', 'FFA22', 'ALC22', 'FFA24',
                           'ALC24', 'ISIS', 'HOH'], overwrite=True)
    print('Creating <{0}.mol2> ...'.format(bilayer.filename))
    bilayer.save(bilayer.filename + '.mol2', overwrite=True)
