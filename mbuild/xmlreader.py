import imp
from mbuild.xyz import Xyz
from treeview import TreeView

__author__ = 'sallai'
from lxml import etree
from mbuild.compound import *
import os.path
import re

class XmlReader(object):

    @classmethod
    def read(cls, path, cwd=""):
        r = XmlReader()
        r.kind_to_compound = {}
        r.src_to_compound = {}
        r.items_by_path = {}

        r.path = os.path.join(cwd, path)

        r.cwd = os.path.dirname(r.path)

        root = etree.parse(r.path)
        r.root = root.getroot()

        r.build(r.root)

        return r.root_compound

    def build(self, elem, ancestor_compounds=[]):
        # assert(isinstance(elem,etree.ElementBase))

        if len(ancestor_compounds) == 0:
            root_compound=None
            parent_compound=None
        else:
            root_compound = ancestor_compounds[0]
            parent_compound = ancestor_compounds[-1]

        if elem.tag == "compound":
            # build a compound

            # see if there's a ref in there
            ref = elem.get("ref")
            src = elem.get("src")

            if ref:
                compound = deepcopy(self.kind_to_compound[ref])
            elif src:
                if src in self.src_to_compound:
                    compound = deepcopy(self.src_to_compound[src])
                else:

                    src_parts = src.split("?")
                    src = src_parts[0]

                    if len(src_parts) == 2:
                        src_params = src_parts[1]
                        regex = re.compile(r"\b(\w+)\s*=\s*([^=]*)(?=\+|&|$)")
                        src_params = dict(regex.findall(src_params))
                        print src_params

                    if src.endswith(".xml"):
                        compound = XmlReader.read(src, cwd=self.cwd)
                    elif src.endswith(".xyz"):
                        compound = Xyz.create(src, cwd=self.cwd)
                    elif src.endswith(".py"):
                        def load_from_file(filepath, expectedClass=None, baseClass=None):
                            class_inst = None

                            mod_name,file_ext = os.path.splitext(os.path.split(filepath)[-1])

                            if file_ext.lower() == '.py':
                                py_mod = imp.load_source(mod_name, filepath)

                            elif file_ext.lower() == '.pyc':
                                py_mod = imp.load_compiled(mod_name, filepath)

                            print py_mod


                            clsmembers = inspect.getmembers(py_mod, inspect.isclass)
                            print clsmembers
                            for label, cls in clsmembers:
                                print label, cls, issubclass(cls, Compound)

                                if not issubclass(cls, baseClass):
                                    continue
                                if not label.lower() == expectedClass.lower():
                                    continue

                                return cls

                        compoundClass = load_from_file(os.path.join(self.cwd,src), baseClass=Compound, expectedClass=os.path.splitext(os.path.split(src)[-1])[0])

                        print compoundClass
                        compound = compoundClass.create(**src_params)
                        # compound = compoundClass.create(*src_params, cwd=self.cwd )
                    else:
                        raise Exception, "don't know how to load " + src
                    self.src_to_compound[src] = compound
                    self.kind_to_compound[compound.kind] = compound
                    compound = deepcopy(compound)
                # self.elem_to_compound[elem] = compound
            else:
                compound = Compound.create(kind=elem.get("kind"))
                self.kind_to_compound[compound.kind] = compound
                # self.elem_to_compound[elem] = compound

            if parent_compound:
                parent_compound.add(compound, elem.get("label"))

            if not root_compound:
                root_compound = compound
                self.root_compound = compound

            for c in elem.iterchildren():
                # assert(isinstance(c,etree.ElementBase))
                part = self.build(c, ancestor_compounds + [compound])

                # # part can be compound, atom, or transform
                # if isinstance(part, Compound):
                #     compound.add(part, c.get("label"))
                # if isinstance(part, Atom):
                #     compound.add(part, c.get("label"))
                # if isinstance(part, CoordinateTransform):
                #     compound.transform(part)


        elif elem.tag == "alias":
            ancestor_compounds[-2].addAlias(parent_compound, elem.get("label") )

        elif elem.tag == "repeat":
            count = int(elem.get("count"))

            for i in range(0,count):
                for c in elem.iterchildren():
                # assert(isinstance(c,etree.ElementBase))
                    part = self.build(c, ancestor_compounds)

        elif elem.tag == "translation":
            coords = elem.get("pos").split(",")
            assert(len(coords) == 3)
            pos = (float(coords[0]), float(coords[1]), float(coords[2]))
            # return Translation(pos)
            parent_compound.transform(Translation(pos))

        elif elem.tag == "rotation":
            axis = elem.get("axis")
            angle = float(elem.get("angle"))*pi/180

            if axis == "x":
                 parent_compound.transform(RotationAroundX(angle))
            elif axis == "y":
                parent_compound.transform(RotationAroundY(angle))
            elif axis == "z":
                parent_compound.transform(RotationAroundZ(angle))

        elif elem.tag == "transform":
            equiv = list()
            for e in elem.iterchildren():
                assert(e.tag == "equivalence")
                my_cpath= e.get("my")
                target_cpath = e.get("target")

                my = root_compound
                for clabel in my_cpath.split("."):
                    if hasattr(my, clabel):
                        my = getattr(my, clabel)
                    else:
                        my = getattr(parent_compound, clabel)

                target = root_compound
                for clabel in target_cpath.split("."):
                    target = getattr(target, clabel)

                equiv.append((my, target))

            parent_compound.transform(Compound.createEquivalenceTransform(equiv))

        else:
            # build an atom
            if not elem.get("pos"):
                raise Exception("expecting an atom, but no pos attribute given (" + self.path + ":" + str(elem.sourceline) +")" )
            coords = elem.get("pos").split(",")
            assert(len(coords) == 3)
            pos = (float(coords[0]), float(coords[1]), float(coords[2]))

            atom = Atom.create(kind=elem.tag, pos=pos)
            parent_compound.add(atom, elem.get("label"))



if __name__ == "__main__":
    # compound = XmlReader.read("xml/methane.xml")
    # compound = XmlReader.read("xml/ethane.xml")
    # compound = XmlReader.read("xml/nalkane.xml")
    compound = XmlReader.read("xml/C_3.xml")
    # compound = XmlReader.read("xml/nalkane2.xml")
    TreeView(compound).show()
