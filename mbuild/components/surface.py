import sys
import os

from mbuild.coordinate_transform import *
from mbuild.compound import Compound
from mbuild.xyz import Xyz
from mbuild.port import Port
from mbuild.rules import RuleEngine
from mbuild.prototype import Prototype


class SurfaceRules(RuleEngine):
    """
    """
    def execute(self):
        self.add_bond("O", "Si", 1.4, 2.5, "o-si")


class Surface(Compound):
    """
    """
    def __init__(self, ctx={}):
        super(Surface, self).__init__(ctx=ctx)

        # Look for xyz file in same directory as this file.
        current_dir = os.path.dirname(os.path.realpath(sys.modules[__name__].__file__))
        xyz_path = os.path.join(current_dir, 'beta-cristobalite.xyz')

        s = Xyz(xyz_path)

        self.add(s, 'surface_xyz')
        self.bounds = [47.689, 41.3, 0.0]

        # Generate topology.
        r = SurfaceRules(self)
        r.execute()

        # We assume here that the surface is in the x-y plane.
        # We add ports pointing upwards in the z-direction
        # where there are undercoordinated oxygens on the surface.
        for atom in self.atoms():
            if atom.kind == "O" and len(atom.bonds) == 1:
                neighbors = self.getAtomsInRange(atom.pos, 3.0, maxItems=10)

                on_top = True
                for neighbor in neighbors:
                    if neighbor.pos[2] > atom.pos[2]:
                        on_top = False

                if on_top:
                    # This is an undercoordinated oxygen with just one bond.
                    p = Port()
                    p.transform(RotationAroundX(pi / 2))
                    p.transform(Translation(atom.pos))
                    p.transform(Translation([0.0, 0.0, 0.5]))
                    self.add(p, str(id(p)))
                    for b in atom.bonds:
                        b.kind = "si-o-top"

if __name__ == "__main__":
    m = Surface()
    Prototype('o-si', color='grey')
    from mbuild.plot import Plot
    Plot(m, bonds=True, verbose=True).show()
