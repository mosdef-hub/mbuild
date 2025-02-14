"""Amorphous silica surface."""

import mbuild as mb


class AmorphousSilicaSurface(mb.Compound):
    """Amorphous silica surface."""

    def __init__(self, surface_roughness=1.0):
        super(AmorphousSilicaSurface, self).__init__()

        if surface_roughness == 1.0:
            mb.load(
                "amorphous_silica_sr1.0.pdb",
                compound=self,
                relative_to_module=self.__module__,
            )
            self.periodicity = (True, True, False)
            self.box = mb.Box([5.4366, 4.7082, 1.0])
        else:
            raise ValueError(
                "Amorphous silica input file with surface "
                "roughness of {0:.1f} does not exist. If you have "
                "this structure, please submit a pull request to "
                "add it! ".format(surface_roughness)
            )
        count = 0
        for particle in list(self.particles()):
            if particle.name == "OB":
                count += 1
                port = mb.Port(anchor=particle, orientation=[0, 0, 1], separation=0.1)
                self.add(port, "port_{}".format(count))


if __name__ == "__main__":
    from mbuild.lib.recipes import TiledCompound

    single = AmorphousSilicaSurface()
    multiple = TiledCompound(single, n_tiles=(2, 1, 1), name="tiled")
    multiple.save("amorphous_silica_surface.mol2")
