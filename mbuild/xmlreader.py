from treeview import TreeView

__author__ = 'sallai'
from lxml import etree
from mbuild.compound import *

class XmlReader(object):

    @classmethod
    def read(cls, f):

        root = etree.parse(f);

        r = XmlReader()
        r.root = root.getroot()
        r.cache = {}

        return r

    def build(self, elem):
        # assert(isinstance(elem,etree.ElementBase))

        if elem.tag == "compound":
            # build a compound

            # see if there's a ref in there
            ref = elem.get("ref")

            if ref:
                compound = deepcopy(self.cache[ref])
            else:
                compound = Compound.create(kind=elem.get("kind"))
                self.cache[compound.kind] = compound

            for c in elem.iterchildren():
                # assert(isinstance(c,etree.ElementBase))
                part = self.build(c)

                # part can be compound, atom, or transform
                if isinstance(part, Compound):
                    compound.add(part, c.get("label"))
                if isinstance(part, Atom):
                    compound.add(part, c.get("label"))
                if isinstance(part, CoordinateTransform):
                    compound.transform(part)

            return compound
        elif elem.tag == "translation":
            coords = elem.get("pos").split(",")
            assert(len(coords) == 3)
            pos = (float(coords[0]), float(coords[1]), float(coords[2]))
            return Translation(pos)

        elif elem.tag == "rotation":
            axis = elem.get("axis")
            angle = float(elem.get("angle"))*pi/180

            if axis == "x":
                return RotationAroundX(angle)
            elif axis == "y":
                return RotationAroundY(angle)
            elif axis == "z":
                return RotationAroundZ(angle)

            return None

        elif elem.tag == "transform":
            return None

        else:
            # build an atom
            coords = elem.get("pos").split(",")
            assert(len(coords) == 3)
            pos = (float(coords[0]), float(coords[1]), float(coords[2]))

            atom = Atom.create(kind=elem.tag, pos=pos)
            return atom


if __name__ == "__main__":
    xr = XmlReader.read("xml/ethane_full.xml")
    print etree.tostring(xr.root)
    c = xr.build(xr.root)
    TreeView(c).show()
