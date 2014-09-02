from copy import deepcopy
import numpy as np
from examples.pmpc.brush import Brush
from mbuild.coordinate_transform import rotate_around_z
from mbuild.file_formats.mol2file import write_mol2, load_mol2
from mbuild.plugins.system import FlatCompound

__author__ = 'sallai'

# create a compound and write it to file
brush1 = Brush()
write_mol2(brush1, filename="brush1.mol2")

# create another compound, rotate it and write it to file
brush2 = Brush()
rotate_around_z(brush2, np.pi/2)
write_mol2(brush2, filename="brush2.mol2")

# load brush2.mol2 into brush1, modifying the atom positions of brush1
sys = load_mol2("brush2.mol2")
sys.update_compound(brush1)
write_mol2(brush1, filename="modified_brush1.mol2")

# access the internals of brush2
print brush1.pmpc

for mpc in brush1.pmpc.monomer:
    print mpc

