__author__ = 'sallai'

from tile import *
from mbuild.ndimmesh import *
from copy import copy, deepcopy

class TwoDimMesh(NDimMesh):
    # n = 3
    # m = 4
    label_prefix = 'part_'
    # left_port_name = 'left_male_port'
    # right_port_name = 'right_female_port'
    # top_port_name = 'top_male_port'
    # bottom_port_name = 'bottom_female_port'

    @classmethod
    def create(cls, creator, n, m, left_port_name, right_port_name, top_port_name, bottom_port_name, label=None):
        cls.creator = creator
        cls.m = m
        cls.n = n
        cls.left_port_name = left_port_name
        cls.right_port_name = right_port_name
        cls.top_port_name = top_port_name
        cls.bottom_port_name = bottom_port_name
        mesh = super(TwoDimMesh, cls).create([cls.builder_0, cls.builder_horz, cls.builder_vert], label)
        return mesh

    @classmethod
    def coord2label(cls, part):
        return cls.label_prefix + str(part.mesh_coords[0]) + '_' + str(part.mesh_coords[1])

    @classmethod
    def builder_0(cls, mesh):
        # create the first part
        new_part = cls.creator()
        new_part.mesh_coords = (0, 0)
        mesh.add(new_part, cls.coord2label(new_part))
        return new_part

    @classmethod
    def builder_horz(cls, mesh, last_part):
        new_part = cls.creator()
        new_part.mesh_coords = (last_part.mesh_coords[0] + 1, last_part.mesh_coords[1])

        # bail out if new part already exists in the mesh
        if mesh.component(cls.coord2label(new_part)) is not None:
            return None

        # bail out if reached the bound
        if new_part.mesh_coords[0] >= cls.n:
            return None

        # transform new part coordinates
        new_part.transform([(new_part.component(cls.left_port_name), last_part.component(cls.right_port_name))])

        # add it to the mesh
        mesh.add(new_part, cls.coord2label(new_part))

        return new_part

    @classmethod
    def builder_vert(cls, mesh, last_part):
        new_part = cls.creator()
        new_part.mesh_coords = (last_part.mesh_coords[0], last_part.mesh_coords[1] + 1)


        # bail out if new part already exists in the mesh
        if mesh.component(cls.coord2label(new_part)) is not None:
            return None

        # bail out if reached the bound
        if new_part.mesh_coords[1] >= cls.m:
            return None

        new_part.transform([(new_part.component(cls.top_port_name), last_part.component(cls.bottom_port_name))])
        mesh.add(new_part, cls.coord2label(new_part))
        return new_part


if __name__ == "__main__":
    seedTile = Tile.create()
    def tileCreator(tdm):
        return deepcopy(seedTile)
    m = TwoDimMesh.create(tileCreator, 7, 5,
                          left_port_name='left_male_port',
                          right_port_name='right_female_port',
                          top_port_name='top_male_port',
                          bottom_port_name='bottom_female_port')
    print m.atoms()
    m.plot(labels=False, verbose=False)
    print m._components





