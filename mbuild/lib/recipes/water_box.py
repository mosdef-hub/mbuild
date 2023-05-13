import numpy as np
import math as math
from warnings import warn

from mbuild import Box, Compound, clone, force_overlap, load
from mbuild.exceptions import MBuildError

import mbuild.lib.molecules.water as water_models

__all__ = ["WaterBox"]

class WaterBox(Compound):
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
    edge: float or list of floats, default=0.1
        Specifies the gutter around the system to avoid overlaps at boundaries
    model: mb.Compound, optional, default=None
        The specified water model to be used in the resulting configuraiton.
        Uses the force overlap command to translate and orient the water model
        to the coordinates; if not specified, distances and angles correspond to tip3p.
        See mbuild/lib/molecules/water.py for available water models.
        
    """
    def __init__(self, box, edge = 0.1, model = None):

        super(WaterBox, self).__init__()
        
        # check if we are given a list or single value
        if isinstance(edge, list):
            assert(len(edge) == 3)
            edges = np.array(edge)
        else:
            edges = np.array([edge,edge,edge])

        # If a model is specified, we will check to ensure that
        # the first particle in the compound corresponds to Oxygen.
        
        if model is not None:
            assert isinstance(model, Compound)
            particles = [p for p in model.particles()]
            if 'O' not in particles[0].name:
                raise MBuildError('The first particle in model needs to correspond to oxygen')
                
        # read in our propotype, a 4.0x4.0x4.0 nm box
        # our prototype was relaxed in GROMACs at 305 K, density 1000 kg/m^3 using tip3p
        aa_waters = load('water_proto.gro')

        # loop over each water in our configuration
        # add in the necessary bonds missing from the .gro file
        # rename particles/Compound according to the given water model
        for water in aa_waters.children:
           
            water.add_bond((water.children[0], water.children[1]))
            water.add_bond((water.children[0], water.children[2]))
            if model is None:
                model = water_models.WaterTIP3P()
                
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

        # add water molecules to a list
        # note we add to a list first, as this is more efficient than calling
        # the Compound.add function repeatedly as the Compound size grows.
        for water in aa_waters.children:
            for i in range(0,scale_Lx):
                for j in range(0,scale_Ly):
                    for k in range(0,scale_Lz):
                        shift = np.array([i*aa_waters.box.Lx, j*aa_waters.box.Ly, k*aa_waters.box.Lz])

                        if all(water.pos+shift < (box.lengths-edges)):
                            temp = clone(water)
                            temp.translate(shift)
                                
                            water_system_list.append(temp)
        
            
        # add to the Compound and set box size
        self.add(water_system_list)
        self.box = box
