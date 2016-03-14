import mbuild as mb
import networkx as nx
import numpy as np

from random import randint
from math import ceil

class SilicaInterface(mb.Compound):
    """ A recipe for creating an interface from bulk silica.

    Carves silica interface from bulk, adjusts to desired surface
    hydoxyl density by creating Si-O-Si bridges, and yields a 2:1
    Si:O ratio (excluding surface binding sites)

    Parameters
    ----------
    bulk_silica : mb.Compound
        Bulk silica from which to cleave an interface
    tile_x : int, optional, default=1
        Number of times to replicate bulk silica in x-direction
    tile_y : int, optional, default=1
        Number of times to replicate bulk silica in y-direction
    thickness : float, optional, default=1.0
        Thickness of the interface (in nm)
    oh_density : float, optional, default=5.0
        Desired density of reactive surface sites (sites/nm^2)

    """

    def __init__(self, bulk_silica, tile_x=1, tile_y=1, thickness=1.0, oh_density=5.0):
        super(SilicaInterface, self).__init__()

        O_buffer = 0.275

        self._cleave_interface(bulk_silica, tile_x, tile_y, thickness, O_buffer)
        self.generate_bonds(name_a='Si', name_b='O', dmin=0.0, dmax=0.20419)
        self._strip_stray_atoms()
        self._bridge_dangling_Os(oh_density, thickness)
        self._identify_surface_sites(thickness)
        self._adjust_stoichiometry(O_buffer)

    def _cleave_interface(self, bulk_silica, tile_x, tile_y, thickness, O_buffer):
        """ Carve interface from bulk silica, include a buffer of O's above and
            below the surface to ensure the interface is coated.
        """ 

        tile_z = int(ceil((thickness+(2*O_buffer))/bulk_silica.periodicity[2]))
        bulk = mb.TiledCompound(bulk_silica, n_tiles=(tile_x, tile_y, tile_z))

        interface = mb.Compound(periodicity=(bulk.periodicity[0],bulk.periodicity[1],
                                             0.0))
        for i, particle in enumerate(bulk.particles()):
            if particle.name == 'Si' and particle.pos[2] < (thickness+O_buffer) and particle.pos[2] > O_buffer:
                temp_particle = mb.Compound(name=particle.name, pos=particle.pos)
                interface.add(temp_particle, particle.name+"_{}".format(i))
            elif particle.name == 'O' and particle.pos[2] < (thickness+2.*(O_buffer)):
                temp_particle = mb.Compound(name=particle.name, pos=particle.pos)
                interface.add(temp_particle, particle.name+"_{}".format(i))
        self.add(interface)

    def _strip_stray_atoms(self):    
        """ Remove stray atoms and surface pieces """

        for i,atom in enumerate(list(self.particles())):
            if len(self.bond_graph.neighbors(atom)) == 0:
                self.remove(atom)

        for atom in list(self.particles()):
            if atom not in max(nx.connected_components(self.bond_graph), key=len):
                self.remove(atom)

    def _bridge_dangling_Os(self, oh_density, thickness):
        """ Create Si-O-Si bridges on the surface to yield the desired
            density of reactive surface sites
        """        

        area = self.periodicity[0] * self.periodicity[1]
        target = int(oh_density * area)

        dangling_Os = []
        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if atom.name == 'O' and atom.pos[2] > thickness:
                    dangling_Os.append(atom)

        n_bridges = int((len(dangling_Os) - target)/2.)

        print("Initial density of surface sites:\t{}".format(len(dangling_Os)/area))

        for i in range(n_bridges):
            bridged = False
            while not bridged:
                j = randint(0, len(dangling_Os)-1)
                O1 = dangling_Os[j]
                assert(len(self.bond_graph.neighbors(O1)) == 1)
                Si1 = self.bond_graph.neighbors(O1)[0]
                for O2 in dangling_Os:
                    assert(len(self.bond_graph.neighbors(O2)) == 1)
                    Si2 = self.bond_graph.neighbors(O2)[0]
                    if O2 == O1:
                        continue
                    if Si1 == Si2:
                        continue
                    if any(neigh in self.bond_graph.neighbors(Si2) for neigh in self.bond_graph.neighbors(Si1)):
                        continue
                    r = self.min_periodic_distance(Si1.pos,Si2.pos)
                    if r < 0.45:
                        bridged = True
                        self.add_bond((O1,Si2))
                        dangling_Os.remove(O1)
                        dangling_Os.remove(O2)
                        self.remove(O2)
                        break

        print("Final density of surface sites:\t{}".format(len(dangling_Os)/area))

    def _identify_surface_sites(self, thickness):
        """ Label surface sites and add ports above them """

        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if atom.name == 'O' and atom.pos[2] > thickness:
                    atom.name = 'OS'
                    port = mb.Port(anchor=atom)
                    mb.rotate_around_x(port, np.pi/2)
                    mb.translate(port, atom.pos + np.array([0.0, 0.0, 0.1]))
                    self.add(port, "port_{}".format(len(self.referenced_ports())))

    def _adjust_stoichiometry(self, O_buffer):
        """ Remove O's from the underside of the surface to yield a 2:1 Si:O ratio """

        O_buffer = 0.275
        num_O = len(np.array([atom for atom in self.particles() if atom.name == 'O']))
        num_Si = len(np.array([atom for atom in self.particles() if atom.name == 'Si']))
        n_deletions = num_O - 2*num_Si

        bottom_Os = []
        for atom in list(self.particles()):
            if len(self.bond_graph.neighbors(atom)) == 1:
                if atom.name == 'O' and atom.pos[2] < O_buffer:
                    bottom_Os.append(atom)

        for i in range(n_deletions):
            O1 = bottom_Os[randint(0, len(bottom_Os)-1)]
            bottom_Os.remove(O1)
            self.remove(O1) 

if __name__ == "__main__":
    from mbuild.lib.bulk_materials import AmorphousSilica
    silica_interface = mb.SilicaInterface(bulk_silica=AmorphousSilica(),thickness=1.2)
    silica_interface.save('silica_interface.mol2', show_ports=True)
