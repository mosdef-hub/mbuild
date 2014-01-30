import operator
from mbuild.compound import Compound
from mbuild.prototype import Prototype

__author__ = 'sallai'
from mayavi import mlab
import numpy as np
import webcolors

class Plot(object):

    def __init__(self, compound, verbose=False, labels=True, atoms=True, bonds=True, angles=True, dihedrals=True):
        assert(isinstance(compound, Compound))

        # display atoms
        if atoms:
            # sort atoms by type
            d = dict()

            for atom in compound.atoms():
                if atom.kind != 'G' or verbose:
                    if not atom.kind in d.keys():
                        d[atom.kind] = [atom]
                    else:
                        d[atom.kind].append(atom)

            # import pdb
            # pdb.set_trace()
            for (kind,atomList) in d.items():
                color = Prototype.getAttr(kind, "color", default="white")
                colorRGB=tuple(map(operator.div, webcolors.name_to_rgb(color), (256.0,256.0,256.0)))
                radius = Prototype.getAttr(kind, "radius", default=1.0)

                x = []
                y = []
                z = []
                r = []
                for atom in atomList:
                    x.append(atom.pos[0])
                    y.append(atom.pos[1])
                    z.append(atom.pos[2])
                    r.append(radius/5.0)


                fig = mlab.points3d(x,y,z,r,color=colorRGB, scale_factor=1.0, scale_mode='scalar')
                fig.glyph.glyph.clamping = False


        # display bonds
        if bonds:
            # sort bonds by type
            d=dict()
            for bond in compound.bonds:
                if not bond.kind in d.keys():
                    d[bond.kind] = [bond]
                else:
                    d[bond.kind].append(bond)

            for (kind, bondList) in d.items():
                color = Prototype.getAttr(kind, "color", default="white")
                colorRGB=tuple(map(operator.div, webcolors.name_to_rgb(color), (256.0,256.0,256.0)))

                x = []
                y = []
                z = []
                u = []
                v = []
                w = []

                # find max tube length
                epsilon = 0.2

                for bond in bondList:

                    pos1 = np.array(bond.atom1.pos)
                    pos2 = np.array(bond.atom2.pos)

                    v12 = pos2 - pos1 # vector from atom1 to atom2
                    d12 = np.linalg.norm(v12) # atom1-atom2 distance
                    p1 = pos1 + v12/d12 * epsilon
                    tube_length = d12-2*epsilon
                    v12 = v12/d12*tube_length

                    x.append(p1[0])
                    y.append(p1[1])
                    z.append(p1[2])
                    u.append(v12[0])
                    v.append(v12[1])
                    w.append(v12[2])

                src = mlab.pipeline.vector_scatter(x, y, z, u, v, w)
                fig = mlab.pipeline.vectors(src, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False

        if angles:
            for angle in compound.angles:
                epsilon = 0.3
                pos1 = np.array(angle.atom1.pos)
                pos2 = np.array(angle.atom2.pos)
                pos3 = np.array(angle.atom3.pos)
                v12 = pos2 - pos1 # vector from atom1 to atom2
                v23 = pos3 - pos2
                d12 = np.linalg.norm(v12) # atom1-atom2 distance
                d23 = np.linalg.norm(v23) # atom1-atom2 distance
                p1 = pos1 + v12/d12 * epsilon
                # p2 = pos1 + v12/d12 * (d12 - epsilon)
                p3 = pos3 - v23/d23 * epsilon
                mlab.plot3d([pos2[0], p1[0]],[pos2[1], p1[1]],[pos2[2], p1[2]],
                        tube_radius=0.20, color=angle.color, opacity=.5)

                mlab.plot3d([pos2[0], p3[0]],[pos2[1], p3[1]],[pos2[2], p3[2]],
                        tube_radius=0.20, color=angle.color, opacity=.5)

        if dihedrals:
            for d in compound.dihedrals:
                epsilon = 0.2
                pos1 = np.array(d.atom1.pos)
                pos2 = np.array(d.atom2.pos)
                pos3 = np.array(d.atom3.pos)
                pos4 = np.array(d.atom4.pos)
                v12 = pos2 - pos1 # vector from atom1 to atom2
                v34 = pos4 - pos3
                d12 = np.linalg.norm(v12) # atom1-atom2 distance
                d34 = np.linalg.norm(v34) # atom1-atom2 distance
                p1 = pos1 + v12/d12 * epsilon
                p4 = pos3 + v34/d34 * (d34 - epsilon)
                mlab.plot3d([pos2[0], p1[0]],[pos2[1], p1[1]],[pos2[2], p1[2]],
                        tube_radius=0.10, color=d.color)

                mlab.plot3d([pos2[0], pos3[0]],[pos2[1], pos3[1]],[pos2[2], pos3[2]],
                        tube_radius=0.10, color=d.color)

                mlab.plot3d([pos3[0], p4[0]],[pos3[1], p4[1]],[pos3[2], p4[2]],
                        tube_radius=0.10, color=d.color)

        self.mlab = mlab

    def show(self):
        mlab.show()