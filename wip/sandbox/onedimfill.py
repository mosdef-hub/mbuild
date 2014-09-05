__author__ = 'sallai'

from mbuild.ndimmesh import *


class OneDimFill(NDimMesh):
    label_prefix = 'part_'

    @classmethod
    def create(cls, creator, n, left_port_name, right_port_name, ctx={}, kind=None):
        cls.creator = staticmethod(creator)
        cls.n = n
        cls.left_port_name = left_port_name
        cls.right_port_name = right_port_name
        mesh = super(OneDimFill, cls).create([cls.builder_0, cls.builder_horz], ctx=ctx, kind=kind)

        box = Box.create(lambda ctx: mesh, ctx)
        return box

    @classmethod
    def coord2label(cls, part):
        return cls.label_prefix + str(part.mesh_coord)

    @classmethod
    def builder_0(cls, mesh):
        # create the first part
        new_part = cls.creator(ctx=mesh.ctx)
        new_part.mesh_coord = 0
        mesh.add(new_part, cls.coord2label(new_part))
        return new_part

    @classmethod
    def builder_horz(cls, mesh, last_part):
        new_part = cls.creator(ctx=mesh.ctx)
        new_part.mesh_coord = last_part.mesh_coord + 1

        # bail out if new part already exists in the mesh
        if mesh.component(cls.coord2label(new_part)) is not None:
            return None

        # bail out if reached the bound
        if new_part.mesh_coord >= cls.n:
            return None

        # transform new part coordinates
        new_part.transform([(new_part.component(cls.left_port_name), last_part.component(cls.right_port_name))])

        # add it to the mesh
        mesh.add(new_part, cls.coord2label(new_part))

        return new_part

if __name__ == "__main__":

    from mbuild.examples.methane.methane import *
    from box import *

    methaneBox = Box.create(Methane.create,ctx={})

    methaneBoxLine = OneDimFill.create(lambda ctx: deepcopy(methaneBox), 5, kind="methane_line",
                          left_port_name='left_male_port',
                          right_port_name='right_female_port')

    methaneBoxPlane = OneDimFill.create(lambda ctx: deepcopy(methaneBoxLine), 5, kind="methane_plane",
                          left_port_name='top_male_port',
                          right_port_name='bottom_female_port')

    methaneBoxSpace = OneDimFill.create(lambda ctx: deepcopy(methaneBoxPlane), 5, kind="methane_space",
                          left_port_name='up_male_port',
                          right_port_name='down_female_port')

    from mbuild.treeview import *
    TreeView(methaneBoxSpace).show()
    methaneBoxSpace.plot3()
    # methaneBoxSpace.plot(labels=False, verbose=False)
    print methaneBoxSpace._components





