import numpy as np

from mbuild.compound import Compound
from mbuild.coordinate_transform import rotate_around_z, rotate_around_y, translate, equivalence_transform
from pegsilane import PegSilane
from mbuild.port import Port
from mbuild.tools.mask import sphere_mask


class Tnp(Compound):
    def __init__(self, n_chains, n_peg, monomer=None):
        super(Tnp, self).__init__()
        self.append_from_file('tnp.pdb')

        # recenter around center of mass
        com = self.center_of_mass()

        translate(self, -com)

        # Generate 65 points on the surface of a unit sphere.
        mask = sphere_mask(n_peg)
        radius = 100

        # Magnify the unit sphere by the provided radius.
        mask *= radius

        # Create particles and Ports at O atoms that are the closest to mask positions.
        pegs = []
        for i, pos in enumerate(mask):
            # find the closest atoms to pos
            neighbors = self.atoms_in_range(pos, 1000)
            # keep only the oxigenes
            o_neighbors = [atom for atom in neighbors if atom.kind == 'O']
            # choose the closest oxigen
            o_atom = o_neighbors[0]

            port = Port(anchor=o_atom)
            self.add(port, "port_{}".format(i))

            # Make the top of the port point toward the negative x axis.
            rotate_around_z(port, np.pi/2)
            # Raise up (or down) the top of the port in the z direction.
            rotate_around_y(port, -np.arcsin(pos[2]/radius))
            # Rotate the Port along the z axis.
            rotate_around_z(port, np.arctan2(pos[1], pos[0]))
            # Move the Port to the position of the oxygene atom on surface of the Sphere.
            translate(port, o_atom.pos)

            peg = PegSilane(n_chains)
            equivalence_transform(peg, peg.up, port)
            # don't add the peg yet to TNP, because this would mess up the nearest neighbor finding
            pegs.append(peg)

        # batch add the pegs
        for peg in pegs:
            self.add(peg, "peg[$]")

if __name__ == "__main__":

    tnp = Tnp(5, 25)
    tnp.visualize()
