import operator
from mbuild.compound import Compound
from mbuild.prototype import Prototype

__author__ = 'sallai'
from mayavi import mlab
import numpy as np
import webcolors

class Plot(object):

    def __init__(self, compound, verbose=False, atoms=True, bonds=True, angles=True, dihedrals=True):
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

                x2 = []
                y2 = []
                z2 = []
                r = []
                for atom in atomList:
                    x2.append(atom.pos[0])
                    y2.append(atom.pos[1])
                    z2.append(atom.pos[2])
                    r.append(radius/5.0)


                fig = mlab.points3d(x2,y2,z2,r,color=colorRGB, scale_factor=1.0, scale_mode='scalar')
                fig.glyph.glyph.clamping = False


        # display bonds
        if bonds:
            # sort bonds by type
            d=dict()
            for item in compound.bonds:
                if not item.kind in d.keys():
                    d[item.kind] = [item]
                else:
                    d[item.kind].append(item)

            for (kind, itemList) in d.items():
                color = Prototype.getAttr(kind, "color", default="white")
                if isinstance(color, basestring):
                    print "color="+color
                    colorRGB=tuple(map(operator.div, webcolors.name_to_rgb(color), (256.0,256.0,256.0)))
                elif isinstance(color,tuple) and len(color) == 3:
                    colorRGB = color

                x2 = []
                y2 = []
                z2 = []
                u = []
                v = []
                w = []

                epsilon = 0.1

                for item in itemList:

                    pos1 = np.array(item.atom1.pos)
                    pos2 = np.array(item.atom2.pos)

                    v12 = pos2 - pos1 # vector from atom1 to atom2
                    d12 = np.linalg.norm(v12) # atom1-atom2 distance
                    p1 = pos1 + v12/d12 * epsilon
                    tube_length = d12-2*epsilon
                    v12 = v12/d12*tube_length

                    x2.append(p1[0])
                    y2.append(p1[1])
                    z2.append(p1[2])
                    u.append(v12[0])
                    v.append(v12[1])
                    w.append(v12[2])

                src = mlab.pipeline.vector_scatter(x2, y2, z2, u, v, w)
                fig = mlab.pipeline.vectors(src, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False

        if angles:
            # for angle in compound.angles:
            #     epsilon = 0.3
            #     pos1 = np.array(angle.atom1.pos)
            #     pos2 = np.array(angle.atom2.pos)
            #     pos3 = np.array(angle.atom3.pos)
            #     v12 = pos2 - pos1 # vector from atom1 to atom2
            #     v23 = pos3 - pos2
            #     d12 = np.linalg.norm(v12) # atom1-atom2 distance
            #     d23 = np.linalg.norm(v23) # atom1-atom2 distance
            #     p1 = pos1 + v12/d12 * epsilon
            #     # p2 = pos1 + v12/d12 * (d12 - epsilon)
            #     p3 = pos3 - v23/d23 * epsilon
            #     mlab.plot3d([pos2[0], p1[0]],[pos2[1], p1[1]],[pos2[2], p1[2]],
            #             tube_radius=0.20, color=angle.color, opacity=.5)
            #
            #     mlab.plot3d([pos2[0], p3[0]],[pos2[1], p3[1]],[pos2[2], p3[2]],
            #             tube_radius=0.20, color=angle.color, opacity=.5)

            # sort angles by kind
            d=dict()
            for item in compound.angles:
                if not item.kind in d.keys():
                    d[item.kind] = [item]
                else:
                    d[item.kind].append(item)

            for (kind, itemList) in d.items():
                color = Prototype.getAttr(kind, "color", default="white")
                if isinstance(color, basestring):
                    colorRGB=tuple(map(operator.div, webcolors.name_to_rgb(color), (256.0,256.0,256.0)))
                elif isinstance(color,tuple) and len(color) == 3:
                    colorRGB = color

                x2 = []
                y2 = []
                z2 = []
                u1 = []
                v1 = []
                w1 = []
                u3 = []
                v3 = []
                w3 = []


                epsilon = 0.1

                for item in itemList:

                    pos1 = np.array(item.atom1.pos)
                    pos2 = np.array(item.atom2.pos)
                    pos3 = np.array(item.atom3.pos)

                    v21 = pos1 - pos2 # vector from atom2 to atom1
                    d21 = np.linalg.norm(v21) # atom1-atom2 distance
                    v23 = pos3 - pos2 # vector from atom2 to atom3
                    d23 = np.linalg.norm(v23) # atom2-atom3 distance

                    tube_length21 = d21-epsilon
                    v21 = v21/d21*tube_length21

                    tube_length23 = d23-epsilon
                    v23 = v23/d23*tube_length23

                    x2.append(pos2[0])
                    y2.append(pos2[1])
                    z2.append(pos2[2])

                    u1.append(v21[0])
                    v1.append(v21[1])
                    w1.append(v21[2])

                    u3.append(v23[0])
                    v3.append(v23[1])
                    w3.append(v23[2])

                src21 = mlab.pipeline.vector_scatter(x2, y2, z2, u1, v1, w1)
                fig = mlab.pipeline.vectors(src21, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False
                src23 = mlab.pipeline.vector_scatter(x2, y2, z2, u3, v3, w3)
                fig = mlab.pipeline.vectors(src23, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False


        if dihedrals:

            # sort angles by kind
            d=dict()
            for item in compound.dihedrals:
                if not item.kind in d.keys():
                    d[item.kind] = [item]
                else:
                    d[item.kind].append(item)

            for (kind, itemList) in d.items():
                color = Prototype.getAttr(kind, "color", default="white")
                if isinstance(color, basestring):
                    colorRGB=tuple(map(operator.div, webcolors.name_to_rgb(color), (256.0,256.0,256.0)))
                elif isinstance(color,tuple) and len(color) == 3:
                    colorRGB = color

                x2 = []
                y2 = []
                z2 = []
                u1 = []
                v1 = []
                w1 = []
                u3 = []
                v3 = []
                w3 = []
                x4 = []
                y4 = []
                z4 = []
                u4 = []
                v4 = []
                w4 = []

                epsilon = 0.1

                for item in itemList:

                    pos1 = np.array(item.atom1.pos)
                    pos2 = np.array(item.atom2.pos)
                    pos3 = np.array(item.atom3.pos)
                    pos4 = np.array(item.atom4.pos)

                    v21 = pos1 - pos2 # vector from atom2 to atom1
                    d21 = np.linalg.norm(v21) # atom1-atom2 distance
                    tube_length21 = d21-epsilon
                    v21 = v21/d21*tube_length21

                    v23 = pos3 - pos2 # vector from atom2 to atom3
                    d23 = np.linalg.norm(v23) # atom2-atom3 distance

                    v43 = pos3 - pos4 # vector from atom4 to atom3
                    d43 = np.linalg.norm(v43) # atom4-atom3 distance
                    tube_length43 = d43-epsilon
                    v43 = v43/d43*tube_length43


                    x2.append(pos2[0])
                    y2.append(pos2[1])
                    z2.append(pos2[2])

                    u1.append(v21[0])
                    v1.append(v21[1])
                    w1.append(v21[2])

                    u3.append(v23[0])
                    v3.append(v23[1])
                    w3.append(v23[2])

                    x4.append(pos4[0])
                    y4.append(pos4[1])
                    z4.append(pos4[2])

                    u4.append(v43[0])
                    v4.append(v43[1])
                    w4.append(v43[2])


                src21 = mlab.pipeline.vector_scatter(x2, y2, z2, u1, v1, w1)
                fig = mlab.pipeline.vectors(src21, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False

                src23 = mlab.pipeline.vector_scatter(x2, y2, z2, u3, v3, w3)
                fig = mlab.pipeline.vectors(src23, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False

                src43 = mlab.pipeline.vector_scatter(x4, y4, z4, u4, v4, w4)
                fig = mlab.pipeline.vectors(src43, mode='cylinder', scale_mode='vector', scale_factor=1.0, color=colorRGB)
                fig.glyph.glyph.clamping = False


        self.mlab = mlab

    def show(self):
        mlab.show()