from mbuild.prototype import Prototype
from mbuild.rules import RuleEngine

__author__ = 'sallai'
from mbuild.file_formats.xyz import *
from mbuild.port import *


class SurfaceRules(RuleEngine):

    def execute(self):
        self.add_bond("O", "Si", 1.4, 2.5, "o-si")

class Surface(Compound):

    def __init__(self, ctx={}):
        super(Surface, self).__init__(ctx=ctx)

        # s = Xyz('amorphous.xyz')
        s = Xyz('beta-cristobalite.xyz')

        self.add(s, 'surface_xyz')
        self.periodicity = [47.689, 41.3, 0.0]

        # create a rule engine
        r = SurfaceRules(self)
        r.execute()

        # we assume here that the surface is in the x-y plane
        # we add ports pointing upwards where there

        for atom in self.atoms():
            if atom.kind =="O" and len(atom.bonds) == 1:

                neighbors = self.getAtomsInRange(atom.pos, 3.0, maxItems=10)

                on_top = True
                for neighbor in neighbors:
                    if neighbor.pos[2] > atom.pos[2]:
                        on_top = False

                if on_top:
                    # this is an oxygen with just one bond
                    # we set up a port here
                    p = Port()
                    p.transform(RotationAroundX(pi/2))
                    p.transform(Translation(atom.pos))
                    p.transform(Translation([0,0,.5]))
                    self.add(p, str(id(p)))

                    for b in atom.bonds:
                        b.kind = "si-o-top"

if __name__ == "__main__":
    m = Surface()
    # print m
    #
    # for a in m.atoms():
    #     print a

    Prototype('o-si', color='grey')
    # Prototype('O', color='blue')
    Plot(m, bonds=True, verbose=True).show()
