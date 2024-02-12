"""A bulk structure of amorphous silica."""

import mbuild as mb


class AmorphousSilicaBulk(mb.Compound):
    """An amorphous silica box.

    density 2.2g/cm^3
    """

    def __init__(self):
        super(AmorphousSilicaBulk, self).__init__()

        mb.load(
            "amorphous_silica_bulk.pdb",
            compound=self,
            relative_to_module=self.__module__,
        )
        self.box = mb.Box([5, 5, 5])
        self.periodicity = (True, True, True)


if __name__ == "__main__":
    bulk = AmorphousSilicaBulk()
    bulk.save("bulk.mol2")
