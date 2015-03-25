from mbuild.tools.parameterize.atomtyper import (
    Element, NeighborCount, NeighborsAtLeast, NeighborsExactly, Whitelist,
    Blacklist, check_atom, InWhitelist)
from mbuild.tools.parameterize.chemical_groups import benzene, dioxolane13


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

    import mbuild as mb
    from mbuild.tools.parameterize.atomtyper import find_atomtypes
    from mbuild.tools.parameterize.forcefield import prepare_atoms

    from mbuild.examples.methane.methane import Methane
    from mbuild.examples.ethane.ethane import Ethane

    m = Methane()
    # m = Ethane()

    traj = m.to_trajectory()
    prepare_atoms(traj.top)
    find_atomtypes(traj.top._atoms, forcefield='UFF')

    for atom in traj.top._atoms:
        print("Atom name={}, opls_type={}".format(atom.name, atom.atomtype))
