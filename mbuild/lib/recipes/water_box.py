import numpy as np
import math as math

from mbuild import Compound, clone, force_overlap, load
from mbuild.exceptions import MBuildError
from collections.abc import Iterable

import mbuild.lib.molecules.water as water_models

__all__ = ["Water3SiteBox"]

def _flatten_list(c_list):
    """Flatten a list.

    Helper function to flatten a list that may be nested, e.g. [comp1, [comp2, comp3]].
    """
    if isinstance(c_list, Iterable) and not isinstance(c_list, str):
        for c in c_list:
            if isinstance(c, Iterable) and not isinstance(c, str):
                yield from _flatten_list(c)
            else:
                yield c
                
class Water3SiteBox(Compound):
    """Generate a box of 3-site water molecules.
    
    Efficiently create an mbuild Compound containing water at density ~1000 kg/m^3
    where local molecule orientations should exist in relatively low energy states.
    This loads in a configuration previously generated with packmol and relaxed with
    GROMACS via NVT simulation at 305K using tip3p model, simulated in a 4 nm^3 box.
    The code will duplicate/truncate the configuration as necessary to satisify the
    given box dimensions.
    
    Parameters
    ----------
    box : mb.Box
        The desired box to fill with water
    edge: float or list of floats, default=0.1 nm
        Specifies the gutter around the system to avoid overlaps at boundaries
    model: mb.Compound, optional, default=water_models.WaterTIP3P()
        The specified 3-site water model to be used. This uses the force overlap
        command to translate and orient the specified water model to the given coordinates.
        See mbuild/lib/molecules/water.py for available water models or extend the base model.
    mask: mb.Compound, optional, default=None
        Remove water molecules from the final configuration that overlap with the Compound
        specified by the mask. If the element field is set, the sum of the particle radii
        will be used.
    radii_dict = dict, optional, None
        User defined radii values based on the name field of a Compound.  This will supercede
    radii_overlap: float, optional, default=0.15 nm
        Default value if the radii_dict or element field are not defined for a particle.
    radii_scaling: float, optional, default=1.0
        Radii are multiplied by this factor during the masking routine. This allows the
        space between the masking particles and the water to be adjusted. This will apply to
        values in radii_dict, radii based on element, and radii_overlap.
        
        
    """
    def __init__(self, box, edge = 0.1, model = water_models.WaterTIP3P(), mask=None, radii_dict=None, radii_overlap = 0.15,  radii_scaling=1.0):

        super(Water3SiteBox, self).__init__()
        
        # if we do not define a dictionary, create an empty one
        if radii_dict is None:
            radii_dict = {}
        else:
            if not isinstance(radii_dict, dict):
                raise ValueError(f'radii_dict should be dictionary.')

        # check if we are given a list or single value
        if isinstance(edge, list):
            if len(edge) != 3:
                raise ValueError(f'edge should either be a single float or a list of length 3, not a list of {len(edge)}.')
            edges = np.array(edge)
        else:
            edges = np.array([edge,edge,edge])

        # If a model is specified, we will check to ensure that
        # the first particle in the compound corresponds to Oxygen.
        
        if model is not None:
            if isinstance(model, Compound) == False:
                raise MBuildError(f'Model must be a compound.')

            particles = [p for p in model.particles()]
            if particles[0].element.symbol != 'O':
                raise MBuildError('The first particle in model needs to correspond to oxygen.')
            if len(particles) !=3:
                raise MBuildError('The only works with 3-site models of water.')

 
        # check if mask is set
        if mask is not None:
            if not isinstance(mask, list):
                if not isinstance(mask, Compound):
                    raise MBuildError(f'Mask must be a Compound or a list of Compounds.')

            elif isinstance(mask, list):
                # in case we are specified a list of Compounds,
                # we will make sure it is a 1d list.
                mask = [e for e in _flatten_list(mask)]
                for entry in mask:
                    if not isinstance(entry, Compound):
                        raise MBuildError(f'Mask must be a Compound or a list of Compounds.')

                    
        # read in our propotype, a 4.0x4.0x4.0 nm box
        # our prototype was relaxed in GROMACs at 305 K, density 1000 kg/m^3 using tip3p
        aa_waters = load('water_proto.gro', relative_to_module=self.__module__,)

        # loop over each water in our configuration
        # add in the necessary bonds missing from the .gro file
        # rename particles/Compound according to the given water model
        for water in aa_waters.children:
           
            water.add_bond((water.children[0], water.children[1]))
            water.add_bond((water.children[0], water.children[2]))
                
            temp = clone(model)
            force_overlap(temp, temp, water, add_bond=False)
            water.name=model.name
            water.children[0].name = model.children[0].name
            water.children[1].name = model.children[1].name
            water.children[2].name = model.children[2].name
            water.xyz = temp.xyz
          
        # scaling parameters for the new box
        scale_Lx = math.ceil(box.Lx/aa_waters.box.Lx)
        scale_Ly = math.ceil(box.Ly/aa_waters.box.Ly)
        scale_Lz = math.ceil(box.Lz/aa_waters.box.Lz)
        
        water_system_list = []
        
        # we will create a list of particles for the mask
        # if specified now to save time later
        if mask is not None:
            if isinstance(mask, Compound):
                p_mask = [ p for p in mask.particles()]
            else:
                p_mask = []
                for entry in mask:
                    p_mask =  p_mask + [ p for p in entry.particles()]
                    
        # add water molecules to a list
        # note we add to a list first, as this is more efficient than calling
        # the Compound.add function repeatedly as the Compound size grows.
        for water in aa_waters.children:
            for i in range(0,scale_Lx):
                for j in range(0,scale_Ly):
                    for k in range(0,scale_Lz):
                        shift = np.array([i*aa_waters.box.Lx, j*aa_waters.box.Ly, k*aa_waters.box.Lz])
                        if all(water.pos+shift < (box.lengths-edges)):
                            if mask is not None:
                                particles = [p for p in water.particles()]
                                status = True
                                
                                # note this could be sped up using a cell list
                                # will have to wait until that PR is merged
                                for p1 in particles:
                                    for p2 in p_mask:
                                        dist= np.linalg.norm(p1.pos-p2.pos)
                                        
                                        if p1.name in radii_dict:
                                            c1 = radii_scaling*radii_dict[p1.name]
                                        elif p1.element is not None:
                                            c1 = radii_scaling*p1.element.radius_alvarez/10.0
                                        else:
                                            c1 = radii_scaling*radii_overlap
                                        
                                        if p2.name in radii_dict:
                                            c2 = radii_scaling*radii_dict[p2.name]
                                        elif p2.element is not None:
                                            c2 = radii_scaling*p2.element.radius_alvarez/10.0
                                        else:
                                            c2 = radii_scaling*radii_overlap
                                            
                                            
                                        cut_value = c1+c2
                                        if dist <= cut_value:
                                            status = False
                                if status:
                                    temp = clone(water)
                                    temp.translate(shift)
                                    water_system_list.append(temp)
                            else:

                                temp = clone(water)
                                temp.translate(shift)
                                water_system_list.append(temp)
        
            
        # add to the Compound and set box size
        self.add(water_system_list)
        self.box = box
