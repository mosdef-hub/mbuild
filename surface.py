from mbuild.moleculemodel import MoleculeModel
from mbuild.rules import RuleEngine

__author__ = 'sallai'
from mbuild.compound import *
from mbuild.xyz import *
from mbuild.port import *


class SurfaceRules(RuleEngine):

    @classmethod
    def create(cls, model):
        re = super(SurfaceRules, cls).create(model)
        return re

    def execute(self):
        self.add_bond(O, Si, 1.4, 2.5, "o-si", (1, 1, 1))

class Surface(Compound):

    @classmethod
    def create(cls, ctx={}):
        s = Xyz.create('amorphous.xyz')

        # set up a molecule model
        mm = MoleculeModel.create(bounds=[47.689, 41.3, 0.0])
        mm.add([atom for label, atom in s.atoms()])

        # create a rule engine

        r = SurfaceRules.create(mm)
        r.execute()



        # we assume here that the surface is in the x-y plane
        # we add ports pointing upwards where there

        m = super(Surface, cls).create(ctx=ctx)
        for atom in mm.atoms:
            m.add(atom, str(id(atom)))
            if isinstance(atom, O) and len(atom.bonds) == 1:

                distances, indices = mm.atomKdtree.query(atom.pos, 10)

                on_top = True
                for index, distance in zip(indices, distances):
                    if distance > 3.0:
                        continue
                    neighbor = mm.atomsList[index]
                    if neighbor.pos[2] > atom.pos[2]:
                        on_top = False

                if on_top:
                    # this is an oxigen with just one bond
                    # we set up a port here
                    p = Port.create()
                    p.transform(RotationAroundX(pi/2))
                    p.transform(Translation(atom.pos))
                    p.transform(Translation([0,0,.5]))
                    m.add(p, str(id(p)))

                    for b in atom.bonds:
                        b.color = (0,0,0)


        mm.plot()


        return m

if __name__ == "__main__":
    m = Surface.create()
    print m

    for a in m.atoms():
        print a
    m.plot(labels=False, verbose=True)
