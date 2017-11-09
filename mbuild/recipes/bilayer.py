from random import shuffle, uniform
import numpy as np
import argparse
import os

import mbuild as mb
from mbuild.compound import Compound
from mbuild.box import Box
from mbuild.packing import *
from mbuild.pattern import Grid2DPattern
from mbuild.exceptions import MBuildError
from mbuild import clone

# for coarse-grained bilayers
# from mbuild.lib.cg_molecules import ECer2, UCer2, Chol, FFAC16, FFAC20, FFAC24, Water

__all__ = ['Bilayer']


class Bilayer(Compound):
    """
    The Bilayer Builder creates a lipid bilayer, solvates it, and stores it as an mBuild Compound.
    Because the Bilayer Builder does not check for the orientation of the mBuild molecules the user 
    provides from a previously created class or file, the user must ensure the following things are true 
    about the lipids they input to the builder:
        - The lipid's longitudinal axis lies on the z-axis of the Cartesian coordinate system
        - The lipid's reference atom is located at the origin (see ref_atoms below)
        - The rest of the lipid is pointing in the negative z direction
        
    The user may input the fraction of each lipid in the bilayer, as well as a number of important bilayer 
    properties such as area per lipid, bilayer size, and tilt angle. For a more robust description of the
    various features of the Bilayer Builder, see the bilayer_tutorial.ipynb file, located in
    mbuild/recipes/bilayer/bilayer_tutorial.ipynb.

    Parameters
    ----------
    lipids : tuple or list of tuples
        List of tuples in format (lipid, frac, z-offset, ref_atom) where frac is
        the fraction of that lipid in the bilayer (lipid is a
        Compound), and z-offset is the distance in nanometers from the
        headgroup to the lipid-water interface
    n_lipids_x : int, optional
        Number of lipids in the x-direction
    n_lipids_y : int, optional
        Number of lipids in the y-direction
    itp_path : string, optional
        Absolute directory path to the itp files that will be used to create the topology file
    area_per_lipid : float, optional
        The area per lipid of the bilayer in nanometers squared
    tilt_angle : float, optional
        Random tilt angle applied equally to all lipids (tilted with respect to the y-axis)
    z_spacing : float, optional
        The distance between the two bilayer leaflets. A negative distance leads to interdigitation
    max_tail_randomization : float, optional, default=25.0
        For a set tilt angle, the maximum amount each lipid may be randomly spun around the z-axis
    mirror : boolean, optional, default=False
        Indicates whether or not the bottom leaflet will be a mirror image of the top leaflet
    cross_tilt : boolean, optional, default=False
        Indicates whether or not the bilayer will have a cross-tilt
    solvent : Compound, optional, default=H2O()
        Compound to solvate the bilayer with
    solvent_per_lipid : int, optional, default=20
        The ratio of solvent particles to lipid particles.
        A value of zero indicates that the bilayer will not be solvated.
    unit_conversion : int, optional, default=1
        Adjust for unit conversions if necessary
    filename : str, optional, default set in the constructor to a combination of lipid components
        and their compositions separated by a _ character
    make_files : bool, optional, default=False
        Indicates whether output files for simulation engines will be created
        
    Attributes
    ----------
    lipids : list of tuples
        List of tuples in format (lipid, frac, z-offset, ref_atom) where frac is
        the fraction of that lipid in the bilayer (lipid is a
        Compound), and z-offset is the distance in nanometers from the
        headgroup to the lipid-water interface
    ref_atoms : int
        Indices of the atom in lipids to form the solvent interface
    n_lipids_x : int
        Number of lipids in the x-direction per layer.
    n_lipids_y : int
        Number of lipids in the y-direction per layer.
    area_per_lipid : float
        The area per lipid of the bilayer in nanometers squared
    tilt : float, optional
        Random tilt angle applied equally to all lipids (tilted with respect to the y-axis)
    z_orientation : float
        For a set tilt angle, the maximum amount each lipid may be randomly spun around the z-axis
    cross_tilt : boolean
        Indicates whether or not the bilayer will have a cross-tilt
    solvent : Compound
        Compound to solvate the bilayer with. Typically, a
        pre-equilibrated box
        of solvent.
    solvent_per_lipid : int
        Number of solvent molecules per lipid
    itp_path : string
        Absolute directory path to the itp files that will be used to create the topology file
    filename : str, optional, default set in the constructor to a combination of lipid components
        and their compositions separated by a _ character
    make_files : bool, optional, default=False
        Indicates whether output files for simulation engines will be created
    unit_conversion : int
        Adjust for unit conversions if necessary
    """
    
    def __init__(self, lipids, n_lipids_x=8, n_lipids_y=8, itp_path="/home/loganguy/builds/setup/FF/gromos53a6/",
                 area_per_lipid=None, tilt_angle=uniform(20, 30), z_spacing=-0.2, max_tail_randomization=25,
                 mirror=False, cross_tilt=False, solvent=None, solvent_per_lipid=20, unit_conversion=1,
                 filename='', make_files=False):
        super(Bilayer, self).__init__()
        print('Setting bilayer attributes and system information...')

        # Ensure the user has entered plausible inputs
        self._sanitize_inputs(lipids,  n_lipids_x, n_lipids_y, tilt_angle, max_tail_randomization,
                              z_spacing, cross_tilt, solvent, solvent_per_lipid, itp_path, make_files,
                              filename, unit_conversion)

        self.ref_atoms = []
        print('Calculating proper area per lipid...')
        area_per_lipid = self.calculate_apl(area_per_lipid)

        # Private attributes for getter methods
        self._number_of_each_lipid_per_layer = []
        self._solvent_per_layer = None

        # Set the 2D lipid locations on a grid pattern
        self._set_grid_pattern(area_per_lipid)

        # Write topology file header
        top_file = None
        if self.make_files is True:
            print('Writing <{0}> ...'.format(self.filename + '.top'))
            top_file = self.write_top_header(self.filename)

        # Create lipid leaflets and add them to the bilayer compound structure
        lipid_bilayer = Compound()
        solvent_components = Compound()
        print('Creating top leaflet...')
        top_file, top_leaflet, top_lipid_labels = self.create_layer(top_file=top_file)
        lipid_bilayer.add(top_leaflet, label='top_leaflet')
        print('Creating bottom leaflet...')
        if mirror:
            top_file, bottom_leaflet, bottom_lipid_labels = self.create_layer(top_file, lipid_indices=top_lipid_labels)
        else:
            top_file, bottom_leaflet, bottom_lipid_labels = self.create_layer(top_file)
        lipid_bilayer.add(bottom_leaflet, label='bottom_leaflet')
        
        # Translate and rotate the bottom leaflet
        self._translate_bottom_leaflet(lipid_bilayer, z_spacing)

        # solvate the lipids
        if self.solvent_per_lipid > 0:
            print('Solvating the bilayer...')
            top_file, solvent_components = self.solvate_bilayer(top_file, lipid_bilayer, solvent_components)

        # recenter the bottom leaflet over the bottom water box
        if self.solvent_per_lipid > 0 and not self.cross_tilt:
            self._post_translate(solvent_components, lipid_bilayer)

        # close the completed topology file
        if self.make_files is True:
            top_file.close()

        # add everything to the overall Compound structure
        print('Adding the system components to mb.Compound...')
        self.add(lipid_bilayer, label='lipid_bilayer')
        if self.solvent_per_lipid > 0:
            self.add(solvent_components, label='solvent')

        # adjust the periodicity of the system, center the system in the new box
        if self.solvent_per_lipid > 0:
            water_box = (np.amax(self['solvent'].xyz, axis=0) - np.amin(self['solvent'].xyz, axis=0)) + 0.5
            self.periodicity = water_box
            distance = -self['solvent'].xyz.min(axis=0) + 0.25
            self.translate(distance)  # Bring the bilayer inside the bounding box for the .gro file

        # Create .gro file and .mol2 file
        if self.make_files is True:
            print('Writing <{0}.gro> ...'.format(self.filename))
            self.save(self.filename + '.gro',
                      residues=['DSPC', 'FFA12', 'ALC12', 'FFA14', 'ALC14', 'FFA16', 'ALC16',
                                'FFA18', 'ALC18', 'FFA20', 'ALC20', 'FFA22', 'ALC22', 'FFA24',
                                'ALC24', 'ISIS', 'HOH'], overwrite=True)
            print('Creating <{0}.mol2> ...'.format(self.filename))
            self.save(self.filename + '.mol2', overwrite=True)
        print('Bilayer construction complete.')

    def create_layer(self, top_file=None, lipid_indices=None):
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
        layer = Compound()
        if not lipid_indices:
            lipid_indices = list(range(self.lipids_per_leaflet))
            shuffle(lipid_indices)
            
        n_previous_types = 0
        for n_type, n_of_lipid_type in enumerate(
                self.number_of_each_lipid_per_leaflet):
            current_type = self.lipids[n_type][0]
            
            # Write to topology file
            if (self.number_of_each_lipid_per_leaflet[n_type] != 0) and self.make_files is True:
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
    
    def _translate_bottom_leaflet(self, lipid_bilayer, spacing_z):

        spacing_z *= self.unit_conversion
        lipid_bilayer['bottom_leaflet'].spin(np.pi, [1, 0, 0])
        top_z_min = np.amin(lipid_bilayer['top_leaflet'].xyz, axis=0)[2]
        bottom_z_max = np.amax(lipid_bilayer['bottom_leaflet'].xyz, axis=0)[2]
        space = top_z_min - bottom_z_max
        lipid_bilayer['bottom_leaflet'].translate([0, 0, (space - spacing_z)])
        if not self.cross_tilt:
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
        water_volume = .025 * np.power(self.unit_conversion, 3)
        # TODO: Support other solvents' volumes in a user-friendly way
        
        # Calculate box dimension
        box_volume = water_volume * self.solvent_per_layer
        x_min_box = min(self.pattern[:, 0])
        x_max_box = max(self.pattern[:, 0])
        y_min_box = min(self.pattern[:, 1])
        y_max_box = max(self.pattern[:, 1])
        
        box_area = (x_max_box - x_min_box) * (y_max_box - y_min_box)
        box_height = box_volume / box_area
        
        z_min_top = max(lipid_bilayer.xyz[:, 2])
        z_max_top = z_min_top + box_height
        
        z_max_bottom = min(lipid_bilayer.xyz[:, 2])
        z_min_bottom = z_max_bottom - box_height

        top_solvent_box = Box(mins=[x_min_box, y_min_box, z_min_top],
                                 maxs=[x_max_box, y_max_box, z_max_top])
        
        bottom_solvent_box = Box(mins=[x_min_box, y_min_box, z_min_bottom],
                                    maxs=[x_max_box, y_max_box, z_max_bottom])

        solvent_components.add(fill_region(compound=[self.solvent, self.solvent],
                                              n_compounds=[self.solvent_per_layer, self.solvent_per_layer],
                                              region=[top_solvent_box, bottom_solvent_box]))

        # Write solvent components to topology file
        if self.make_files is True:
            top_file.write("{:<10s}{:<10d}\n".format('SOL', self.solvent_per_layer * 2))

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

        for lpd in self.lipids[:-1]:
            self._number_of_each_lipid_per_layer.append(int(round(lpd[1] * self.lipids_per_leaflet)))
            estimate = abs(round(lpd[1] * self.lipids_per_leaflet))
            actual = (lpd[1] * self.lipids_per_leaflet)
            if (abs(estimate - actual) / actual) > 0.1:
                message = 'The Bilayer Builder was unable to create a bilayer with the exact fraction ' \
                          'specified for lipid {}. \nPlease check that the combination of bilayer size and' \
                          ' given composition fractions is logical'.format(lpd[0])
                raise MBuildError(message)
        self._number_of_each_lipid_per_layer.append(self.lipids_per_leaflet
                                                    - sum(self._number_of_each_lipid_per_layer))
        assert len(self._number_of_each_lipid_per_layer) == len(self.lipids)
        return self._number_of_each_lipid_per_layer

    def _sanitize_inputs(self, lipids, n_lipids_x, n_lipids_y, tilt_angle,
                         max_tail_randomization, z_spacing, cross_tilt, solvent, solvent_per_lipid,
                         itp_path, make_files, filename, unit_conversion):
        """Checks that the values passed by the user to the constructor are valid before instance attributes
        are created.
    
        Ensure that the user's lipid fractions add up to 1, 
        or raise a ValueError.
        Ensure that the user has not specified invalid bilayer dimensions,
        or raise a ValueError.
        Ensure that the user has not provided more lipids than available grid spaces,
        or raise a ValueError.
        """
        if not isinstance(lipids, (list, tuple, set)):
            raise TypeError('Parameter lipids must be an iterable: lipids parameter is of type {}'
                            .format(type(lipids)))

        # Check if lipids has more than one lipid element
        if isinstance(lipids[0], Compound):
            self.lipids = [lipids]
        else:
            self.lipids = lipids
        modellipid = (Compound(), 1.0, -0.1, 0)
        for i, lpd in enumerate(self.lipids):
            if len(lpd) != 4:
                raise ValueError('One or more lipid tuples are missing required elements')
            for j, element in enumerate(lpd):
                if not isinstance(element, type(modellipid[j])):
                    raise TypeError('Element {} of lipid {} is of invalid type: given type is {}, should be {}'
                                    .format(j, (i + 1), type(element), type(modellipid[j])))
        if sum([lpd[1] for lpd in self.lipids]) != 1.0:
            raise ValueError('Lipid fractions do not add up to 1.')
        if not (n_lipids_x > 0 and n_lipids_y > 0):
            raise ValueError('Bilayer dimensions must be greater than 0')
        if not isinstance(n_lipids_x, int) or not isinstance(n_lipids_y, int):
            raise TypeError('Both dimensions provided must be an integer: n_lipids_x is of type {}, '
                            'n_lipids_y is of type {}'.format(type(n_lipids_x), type(n_lipids_y)))
        if len(self.lipids) > (n_lipids_x * n_lipids_y):
            raise ValueError('Number of lipids provided exceeds number of available grid spaces')
        self.n_lipids_x = n_lipids_x
        self.n_lipids_y = n_lipids_y

        # Geometric attributes

        if not isinstance(tilt_angle, (float, int)):
            raise TypeError('Tilt angle must be a valid number.')
        self.tilt = tilt_angle * np.pi / 180
        if not isinstance(max_tail_randomization, (float, int)):
            raise TypeError('Z-axis randomization must be a valid number.')
        if max_tail_randomization >= 90:
            raise ValueError('The tail randomization maximum provided is too large. Please pick an angle < 90)')
        self.z_orientation = max_tail_randomization
        if not isinstance(z_spacing, float):
            raise TypeError('z_spacing parameter must be of type float: z_spacing provided is of type {}'
                            .format(type(z_spacing)))
        self.z_spacing = z_spacing
        if not isinstance(cross_tilt, bool):
            raise TypeError('cross_tilt parameter must be a valid boolean.')
        self.cross_tilt = cross_tilt

        # Solvent attributes
        if solvent is not None:
            if not isinstance(solvent, Compound):
                raise TypeError('Solvent provided must be a valid mb.Compound')
        else:
            from mbuild.lib.prototypes import H2O
            solvent = H2O()
        self.solvent = solvent
        if not isinstance(solvent_per_lipid, int):
            raise ValueError('solvent_per_lipid must be an integer')
        self.solvent_per_lipid = solvent_per_lipid

        # Path to .itp files and file attributes
        if not isinstance(make_files, bool):
            raise TypeError('make_files parameter must be a valid boolean')
        self.make_files = make_files
        if not isinstance(filename, str):
            raise TypeError('Filename must be a valid string')
        if not isinstance(itp_path, str):
            raise TypeError('Directory path to itp files must be a valid string')
        if make_files is True:
            if not os.path.exists(itp_path):
                raise IOError('The provided itp file path does not exist')
        self.itp_path = itp_path
        self.filename = filename
        if len(self.filename) == 0:
            components = [str(lip[0].name) + '_' + str(lip[1]) + '_' for lip in self.lipids]
            for component in components:
                self.filename += component
            self.filename += 'UA'
            # self.filename += 'AA' TODO: Make a command line flag
            # self.filename += 'CG' TODO: Same as above

        self.unit_conversion = unit_conversion

    def calculate_apl(self, area_per_lipid):
        """Calculate the proper APL based on lipid input"""
        multiplier = 0  # placeholder that will be used to calculate the necessary area per lipid

        for component in self.lipids:
            self.ref_atoms.append(component[3])
            if component[0] is ('DSPC' or 'DPPC' or 'DMPC' or 'ISIS'):
                multiplier += component[1]

        # Calculate the necessary area per lipid
        if area_per_lipid is None:
            area_per_lipid = uniform(0.2 + (0.3 * multiplier), 0.3 + (0.3 * multiplier))
        else:
            if not isinstance(area_per_lipid, (float, int)):
                raise TypeError('Area per lipid must be a valid number.')
            area_per_lipid = area_per_lipid
        return area_per_lipid

    def _set_grid_pattern(self, area_per_lipid):
        """Utilize an mBuild 2DGridPattern to create the scaffold of points that the lipids will be laid onto"""
        
        self.lipids_per_leaflet = self.n_lipids_x * self.n_lipids_y
        pattern = Grid2DPattern(self.n_lipids_x, self.n_lipids_y)
        pattern.scale(np.sqrt(area_per_lipid * self.lipids_per_leaflet) * self.unit_conversion)
        
        self.pattern = pattern

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Welcome to the Bilayer Builder!',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog='For more robust descriptions of variables and '
                                            'functionality of the Bilayer Builder, see tutorial_bilayer.ipynb, '
                                            'found in the tutorials directory.')
    
    parser.add_argument('n_lipids_x', type=int, help='The number of lipids in the x direction')
    parser.add_argument('n_lipids_y', type=int, help='The number of lipids in the y direction')
    parser.add_argument('-i', '--itp_path', type=str,
                        default="/home/loganguy/builds/setup/FF/gromos53a6/",
                        help='The absolute path on the local machine to the force field .itp files for the lipids in '
                             'the bilayer; necessary for proper construction of the topology file and for proper '
                             'GROMACS simulation')
    frac = parser.add_argument_group('Component Fractions', 'The fraction of the bilayer containing each component')
    frac.add_argument('--DSPC', type=float, default=0.0)
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
    geometry.add_argument('-z', '--spacing', type=float, default=-0.2,
                          help='The spacing in between the two bilayer leaflets')
    geometry.add_argument('-r', '--randomization', type=float, default=25,
                          help='The maximum degrees the lipid tails can be randomly rotated about the z-axis')
    geometry.add_argument('-m', '--mirror', action='store_true', default=False,
                          help='Makes the top and bottom leaflets mirror images of one another')
    geometry.add_argument('-c', '--cross_tilt', action='store_true', default=False,
                          help='Induces a cross-tilt in the bilayer')
    parser.add_argument('-s', '--solvate', type=int, default=20, help='Solvates the bilayer with a specified number '
                                                                      ' of solvent molecules per lipid')
    parser.add_argument('-u', '--units', type=float, default=1.0, help='Input any unit conversion that is necessary '
                                                                       'here (mBuild default units are in nanometers)')
    cmdline = parser.parse_args()

    from mbuild.lib.UA_molecules import DSPCUA, DMPCUA, DPPCUA, FFAUA, ISISUA, ALCUA

    lipids = [(DSPCUA(), cmdline.DSPC, 0.0, 0),
              (DPPCUA(), cmdline.DPPC, -0.3, 0),
              (DMPCUA(), cmdline.DMPC, -0.5, 0),
              (ALCUA(12), cmdline.alc12, -0.2, 13),
              (FFAUA(12, ester=False), cmdline.acd12, -0.2, 13),
              (ALCUA(14), cmdline.alc14, -0.2, 15),
              (FFAUA(14, ester=False), cmdline.acd14, -0.2, 15),
              (ALCUA(16), cmdline.alc16, -0.4, 17),
              (FFAUA(16, ester=False), cmdline.acd16, -0.4, 17),
              (ALCUA(18), cmdline.alc18, -0.4, 19),
              (FFAUA(18, ester=False), cmdline.acd18, -0.4, 19),
              (ALCUA(20), cmdline.alc20, -0.5, 21),
              (FFAUA(20, ester=False), cmdline.acd20, -0.5, 21),
              (ALCUA(22), cmdline.alc22, -0.5, 23),
              (FFAUA(22, ester=False), cmdline.acd22, -0.5, 23),
              (ALCUA(24), cmdline.alc24, -0.4, 25),
              (FFAUA(24, ester=False), cmdline.acd24, -0.4, 25),
              (ISISUA(), cmdline.ISIS, -0.4, 20)]

    # Remove all lipid tuples whose fractions are 0
    lipids_to_pop = []
    for index, lipid in enumerate(lipids):
        if lipid[1] == 0.0:
            lipids_to_pop.append(index)
    for index in sorted(lipids_to_pop, reverse=True):
        lipids.pop(index)

    # Set all parameter values
    n_x = cmdline.n_lipids_x
    n_y = cmdline.n_lipids_y
    path = cmdline.itp_path
    if cmdline.apl is None:
        apl = uniform(0.2 + (0.3 * (cmdline.DSPC + cmdline.DPPC + cmdline.DMPC + cmdline.ISIS)),
                      0.3 + (0.3 * (cmdline.DSPC + cmdline.DPPC + cmdline.DMPC + cmdline.ISIS)))
    else:
        apl = cmdline.apl
    tilt = cmdline.tilt_angle
    z_space = cmdline.spacing
    tail_randomization = cmdline.randomization
    mir = cmdline.mirror
    cross = cmdline.cross_tilt
    spl = cmdline.solvate
    units = cmdline.units

    bilayer = Bilayer(lipids, n_lipids_x=n_x, n_lipids_y=n_y, itp_path=path,
                      area_per_lipid=apl, tilt_angle=tilt, z_spacing=z_space,
                      max_tail_randomization=tail_randomization, mirror=mir, cross_tilt=cross,
                      solvent_per_lipid=spl, unit_conversion=units, make_files=True)
