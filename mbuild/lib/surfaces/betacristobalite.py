"""Beta-cristobalite surface."""
import itertools
import math

import numpy as np

import mbuild as mb


class Betacristobalite(mb.Compound):
    """The beta-cristobalite form of SiO2.

    Area per port specifies the density of attachment sites in nm^2.
    The crystal is expanded to yield an area per port of 0.25 nm^2, the
    typical density of alkane monolayers on SiO2 although these are actually
    grown on amorphous SiO2 in experiment.

    See http://www.wikiwand.com/en/Silicon_dioxide for more info on the various
    crystal forms.

    Note: Port sites are currently naively determined by placing them on all
    oxygens which are above 1.0 nm in the z-direction. This only holds true for
    the beta-cristobalite-expanded.mol2 file. If you add a new one, please
    modify the file or the method of determining port locations.
    """

    def __init__(self, dimensions=None):
        super(Betacristobalite, self).__init__(dimensions)

        betacristabolite = mb.load(
            "beta-cristobalite-expanded.mol2",
            relative_to_module=self.__module__,
            backend="gmso",
        )
        self.periodicity = (True, True, False)
        # 1.3200 taken from boundingbox length rounded to 4 decimal places
        ref_dims = [5.3888, 4.6669, 1.3200]
        if not dimensions:
            box = mb.Box(ref_dims)
        elif isinstance(dimensions, (list, tuple)) and len(dimensions) == 2:
            box = mb.Box([dimensions[0], dimensions[1], ref_dims[2]])
        elif isinstance(dimensions, (list, tuple)) and len(dimensions) == 3:
            box = mb.Box(dimensions)
        else:
            raise ValueError(
                "Unsupported value for dimension provided. "
                "Only accept `None` (default), list/tuple of 2 (x and y dimension) "
                "or list/tuple of 3(x, y, and z dimensions)"
            )

        self.box = box

        # Borrowed code from water_box recipe
        scale_Lx = math.ceil(box.Lx / ref_dims[0])
        scale_Ly = math.ceil(box.Ly / ref_dims[1])
        scale_Lz = math.ceil(box.Lz / ref_dims[2])

        silica_list = list()
        for particle in list(betacristabolite.particles()):
            for i, j, k in itertools.product(
                range(scale_Lx),
                range(scale_Ly),
                range(scale_Lz),
            ):
                shift = np.array(
                    [
                        i * ref_dims[0],
                        j * ref_dims[1],
                        k * ref_dims[2],
                    ]
                )
                if all(particle.pos + shift < box.lengths):
                    temp = mb.clone(particle)
                    temp.translate(shift)
                    silica_list.append(temp)

        self.add(silica_list)

        count = 0
        for particle in list(self.particles()):
            if particle.name.startswith("O") and particle.pos[2] > 1.0:
                count += 1
                port = mb.Port(
                    anchor=particle, orientation=[0, 0, 1], separation=0.1
                )
                self.add(port, "port_{}".format(count))
                particle.name = "O"  # Strip numbers required in .mol2 files.
            elif particle.name.startswith("Si"):
                particle.name = "Si"


if __name__ == "__main__":
    single = Betacristobalite()
    multiple = mb.recipes.TiledCompound(single, n_tiles=(2, 1, 1), name="tiled")
    multiple.save("betacristobalite.mol2", overwrite=True)
