import numpy as np
import mbuild as mb

class AmorphousSilica(mb.Compound):
    """ An amorphous silica box, 2.2g/cm^3"""

    def __init__(self):
        super(AmorphousSilica, self).__init__()

        mb.load('amorphous_silica.pdb', compound=self,
                relative_to_module=self.__module__)
        self.periodicity = np.array([5, 5, 5])

if __name__ == "__main__":
    bulk = AmorphousSilica()
    bulk.save('bulk.mol2')
