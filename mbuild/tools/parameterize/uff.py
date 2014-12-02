from atomtyper import (Element, NeighborCount, NeighborsAtLeast,
    NeighborsAtMost, NeighborsExactly, Whitelist, Blacklist, check_atom,
    InWhitelist)
from chemical_groups import benzene


# -------------- #
# House of rules #
# -------------- #


@Element('C')
@NeighborCount(4)
@NeighborsExactly('H', 4)
@Whitelist('C_3')
def uff_C_3(atom):
    """ """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist('H_')
def uff_H_(atom):
    """ """
    return True


if __name__ == "__main__":
    import pdb
    from mbuild.examples.methane.methane import Methane

    m = Methane()

    from atomtyper import find_atomtypes
    find_atomtypes(m, forcefield='UFF')

    for i, atom in enumerate(m.atoms):
        print "Atom kind={}, uff_type={}".format(
            atom.kind, atom.atomtype)