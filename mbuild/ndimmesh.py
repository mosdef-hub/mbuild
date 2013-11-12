__author__ = 'sallai'

from tile import *


class NDimMesh(Compound):

    def build(self, builders):

        def part_factory(last_part):
            if last_part is None:
                return None

            for builder in builders[1:]:
                new_part = builder(self, last_part)
                part_factory(new_part)

        first_part = builders[0](self)
        part_factory(first_part)

    @classmethod
    def create(cls, builders, ctx={}, kind=None):

        mesh = super(NDimMesh, cls).create(ctx=ctx, kind=kind)

        mesh.build(builders)
        return mesh

if __name__ == "__main__":

    n = 13
    m = 4

    def builder_0(mesh):
        # create the first part
        new_part = Tile.create(ctx=mesh.ctx)
        new_part.mesh_coords = (0,0)
        new_part_label = 'part_'+str(new_part.mesh_coords[0])+'_'+str(new_part.mesh_coords[1])
        mesh.add(new_part, new_part_label)
        return new_part

    def builder_horz(mesh, last_part):
        new_part = Tile.create(ctx=mesh.ctx)
        new_part.mesh_coords = (last_part.mesh_coords[0] + 1, last_part.mesh_coords[1])
        new_part_label = 'part_'+str(new_part.mesh_coords[0])+'_'+str(new_part.mesh_coords[1])

        # bail out if new part already exists in the mesh
        if mesh.component(new_part_label) is not None:
            return None

        # bail out if reached the bound
        if new_part.mesh_coords[0] >= n:
            return None

        # transform new part coordinates
        new_part.transform([(new_part.component('left_male_port'), last_part.component('right_female_port'))])

        # add it to the mesh
        mesh.add(new_part, new_part_label)

        return new_part

    def builder_vert(mesh, last_part):
        new_part = Tile.create(ctx=mesh.ctx)
        new_part.mesh_coords = (last_part.mesh_coords[0], last_part.mesh_coords[1] + 1)
        new_part_label = 'part_'+str(new_part.mesh_coords[0])+'_'+str(new_part.mesh_coords[1])

        # bail out if new part already exists in the mesh
        if mesh.component(new_part_label) is not None:
            return None

        # bail out if reached the bound
        if new_part.mesh_coords[1] >= m:
            return None

        new_part.transform([(new_part.component('top_male_port'), last_part.component('bottom_female_port'))])
        mesh.add(new_part, new_part_label)
        return new_part


    m = NDimMesh.create([builder_0, builder_horz, builder_vert])
    print m.atoms()
    m.plot(labels=False, verbose=False)
    print m._components





