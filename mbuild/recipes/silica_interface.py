from __future__ import division

import math
import random

import mbuild as mb
import numpy as np


class SilicaInterface(mb.Compound):
    """ A recipe for creating an interface from bulk silica.

    Carves silica interface from bulk, adjusts to a reactive
    surface site density of 5.0 sites/nm^2 (agreeing with experimental
    results, see Zhuravlev 2000) by creating Si-O-Si bridges, and
    yields a 2:1 Si:O ratio (excluding the reactive surface sites).

    Parameters
    ----------
    bulk_silica : mb.Compound
        Bulk silica from which to cleave an interface
    tile_x : int, optional, default=1
        Number of times to replicate bulk silica in x-direction
    tile_y : int, optional, default=1
        Number of times to replicate bulk silica in y-direction
    thickness : float, optional, default=1.0
        Thickness of the slab to carve from the silica bulk. (in nm; not
        including oxygen layers on the top and bottom of the surface)

    References
    ----------
    .. [1] Hartkamp, R., Siboulet, B., Dufreche, J.-F., Boasne, B.
           "Ion-specific adsorption and electroosmosis in charged
           amorphous porous silica." (2015) Phys. Chem. Chem. Phys.
           17, 24683-24695
    .. [2] L.T. Zhuravlev, "The surface chemistry of amorphous silica.
           Zhuravlev model." (2000) Colloids Surf., A. 10, 1-38

    """

    def __init__(self, bulk_silica, tile_x=1, tile_y=1, thickness=1.0, seed=12345):
        super(SilicaInterface, self).__init__()

        random.seed(seed)
        self._oh_density = 5.0
        self._O_buffer = 0.275

        self._cleave_interface(bulk_silica, tile_x, tile_y, thickness)
        self.generate_bonds(name_a='Si', name_b='O', dmin=0.0, dmax=0.20419)
        self._strip_stray_atoms()
        self._bridge_dangling_Os(self._oh_density, thickness)
        self._identify_surface_sites(thickness)
        self._adjust_stoichiometry()

    def _cleave_interface(self, bulk_silica, tile_x, tile_y, thickness):
        """Carve interface from bulk silica.

        Also includes a buffer of O's above and below the surface to ensure the
        interface is coated.
        """
        O_buffer = self._O_buffer
        tile_z = int(math.ceil((thickness + 2*O_buffer) / bulk_silica.periodicity[2]))
        bulk = mb.TiledCompound(bulk_silica, n_tiles=(tile_x, tile_y, tile_z))

        interface = mb.Compound(periodicity=(bulk.periodicity[0],
                                             bulk.periodicity[1],
                                             0.0))
        for i, particle in enumerate(bulk.particles()):
            if ((particle.name == 'Si' and O_buffer < particle.pos[2] < (thickness + O_buffer)) or 
                    (particle.name == 'O' and particle.pos[2] < (thickness + 2*O_buffer))):
                interface_particle = mb.Compound(name=particle.name, pos=particle.pos)
                interface.add(interface_particle, particle.name + "_{}".format(i))
        self.add(interface)

    def _strip_stray_atoms(self):
        """Remove stray atoms and surface pieces. """
        components = self.bond_graph.connected_components()
        major_component = max(components, key=len)
        for atom in list(self.particles()):
            if atom not in major_component:
                self.remove(atom)

    def _bridge_dangling_Os(self, oh_density, thickness):
        """Form Si-O-Si bridges to yield desired density of reactive surface sites.

        References
        ----------
        .. [1] Hartkamp, R., Siboulet, B., Dufreche, J.-F., Boasne, B.
               "Ion-specific adsorption and electroosmosis in charged
               amorphous porous silica." (2015) Phys. Chem. Chem. Phys.
               17, 24683-24695
        """

        area = self.periodicity[0] * self.periodicity[1]
        target = int(oh_density * area)

        dangling_Os = [atom for atom in self.particles()
                       if atom.name == 'O' and
                       atom.pos[2] > thickness and
                       len(self.bond_graph.neighbors(atom)) == 1]

        n_bridges = int((len(dangling_Os) - target) / 2)

        for _ in range(n_bridges):
            bridged = False
            while not bridged:
                O1 = random.choice(dangling_Os)
                Si1 = self.bond_graph.neighbors(O1)[0]
                for O2 in dangling_Os:
                    if O2 == O1:
                        continue
                    Si2 = self.bond_graph.neighbors(O2)[0]
                    if Si1 == Si2:
                        continue
                    if any(neigh in self.bond_graph.neighbors(Si2)
                           for neigh in self.bond_graph.neighbors(Si1)):
                        continue
                    r = self.min_periodic_distance(Si1.pos, Si2.pos)
                    if r < 0.45:
                        bridged = True
                        self.add_bond((O1, Si2))
                        dangling_Os.remove(O1)
                        dangling_Os.remove(O2)
                        self.remove(O2)
                        break

    def _identify_surface_sites(self, thickness):
        """Label surface sites and add ports above them. """
        for atom in self.particles():
            if len(self.bond_graph.neighbors(atom)) == 1:
                if atom.name == 'O' and atom.pos[2] > thickness:
                    atom.name = 'OS'
                    port = mb.Port(anchor=atom)
                    port.spin(np.pi/2, [1, 0, 0])
                    port.translate(np.array([0.0, 0.0, 0.1]))
                    self.add(port, "port_{}".format(len(self.referenced_ports())))

    def _adjust_stoichiometry(self):
        """Remove O's from underside of surface to yield a 2:1 Si:O ratio. """
        num_O = len(list(self.particles_by_name('O')))
        num_Si = len(list(self.particles_by_name('Si')))
        n_deletions = num_O - 2*num_Si

        bottom_Os = [atom for atom in self.particles()
                     if atom.name == 'O' and
                        atom.pos[2] < self._O_buffer and
                        len(self.bond_graph.neighbors(atom)) == 1]

        for _ in range(n_deletions):
            O1 = random.choice(bottom_Os)
            bottom_Os.remove(O1)
            self.remove(O1)

if __name__ == "__main__":
    from mbuild.lib.bulk_materials import AmorphousSilica
    silica_interface = mb.SilicaInterface(bulk_silica=AmorphousSilica(), thickness=1.2)
    silica_interface.save('silica_interface.mol2', show_ports=True)
