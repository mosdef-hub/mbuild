from atomtyper import (Element, NeighborCount, NeighborsAtLeast,
    NeighborsAtMost, NeighborsExactly, Whitelist, Blacklist, check_atom,
    InWhitelist)
from chemical_groups import benzene, dioxolane13


# -------------- #
# House of rules #
# -------------- #

@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 1)
@NeighborsExactly('H', 3)
@Whitelist(135)
def opls_135(atom):
    """alkane CH3 """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 2)
@NeighborsExactly('H', 2)
@Whitelist(136)
def opls_136(atom):
    """alkane CH2 """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 3)
@NeighborsExactly('H', 1)
@Whitelist(137)
def opls_137(atom):
    """alkane CH """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('H', 4)
@Whitelist(138)
def opls_138(atom):
    """alkane CH4 """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 4)
@Whitelist(139)
def opls_139(atom):
    """alkane C """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(140)
def opls_140(atom):
    """alkane H """
    return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 3)
@Whitelist(141)
def opls_141(atom):
    """alkene C (R2-C=) """
    return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 2)
@NeighborsExactly('H', 1)
@Whitelist(142)
def opls_142(atom):
    """alkene C (RH-C=) """
    return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 1)
@NeighborsExactly('H', 2)
@Whitelist(143)
def opls_143(atom):
    """alkene C (H2-C=) """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(144)
@Blacklist(140)
def opls_144(atom):
    """alkene H (H-C=) """
    # Make sure that the carbon is an alkene carbon.
    rule_ids = [141, 142, 143]
    return check_atom(atom.neighbors[0], rule_ids)


@Element('C')
@NeighborCount(3)
@NeighborsAtLeast('C', 2)
@Whitelist(145)
@Blacklist([141, 142])
def opls_145(atom):
    """Benzene C - 12 site JACS,112,4768-90. Use #145B for biphenyl """
    return benzene(atom)


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 3)
@Whitelist('145B')
@Blacklist([145])
def opls_145B(atom):
    """Biphenyl C1 """
    # Store for checking neighbors outside the first ring.
    ring_one = benzene(atom)
    if ring_one:
        for neighbor in atom.neighbors:
            if neighbor not in ring_one:
                if benzene(neighbor):
                    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(146)
@Blacklist([140, 144])
def opls_146(atom):
    """Benzene H - 12 site. """
    return check_atom(atom.neighbors[0], 145)


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 1)
@NeighborsExactly('H', 3)
@Whitelist(148)
@Blacklist(135)
def opls_148(atom):
    """C: CH3, toluene """
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            if check_atom(neighbor, 145):
                return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 2)
@NeighborsExactly('H', 2)
@Whitelist(149)
@Blacklist(136)
def opls_149(atom):
    """C: CH2, ethyl benzene """
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            if check_atom(neighbor, 145):
                return True


@Element('O')
@NeighborCount(2)
@NeighborsExactly('H', 1)
@Whitelist(154)
def opls_154(atom):
    """all-atom O: mono alcohols """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('O', 1)
@Whitelist(155)
def opls_155(atom):
    """all-atom H(O): mono alcohols, OP(=O)2 """
    return check_atom(atom.neighbors[0], 154)


@Element('C')
@NeighborCount(4)
@NeighborsExactly('O', 1)
@NeighborsExactly('C', 1)
@NeighborsExactly('H', 2)
@Whitelist(218)
def opls_218(atom):
    """C in CH2OH - benzyl alcohols """
    benzene_carbon = False
    alcohol_oxygen = False
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            benzene_carbon = check_atom(neighbor, 145)
        if neighbor.kind == 'O':
            alcohol_oxygen = check_atom(neighbor, 154)
    if benzene_carbon and alcohol_oxygen:
        return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 3)
@Whitelist(221)
@Blacklist([145, '145B'])
def opls_221(atom):
    """C(CH2OH)   - benzyl alcohols """
    if check_atom(atom, 145):  # Already identified as part of benzene.
        for neighbor in atom.neighbors:
            if check_atom(neighbor, 218):
                return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 1)
@NeighborsExactly('H', 1)
@NeighborsExactly('O', 1)
@Whitelist(232)
def opls_232(atom):
    """C: C=0 in benzaldehyde, acetophenone (CH) """
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            benzene_carbon = check_atom(neighbor, 145)
        if neighbor.kind == 'O':
            aldehyde_oxygen = check_atom(neighbor, 278)
    if benzene_carbon and aldehyde_oxygen:
        return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('C', 2)
@NeighborsExactly('O', 1)
@Whitelist(233)
@Blacklist(145)
def opls_233(atom):
    """C: C=0 in acetophenone (CMe) """
    for neighbor in atom.neighbors:
        if neighbor.kind == 'C':
            benzene_carbon = check_atom(neighbor, 145)
        if neighbor.kind == 'O':
            aldehyde_oxygen = check_atom(neighbor, 278)
    if benzene_carbon and aldehyde_oxygen:
        return True


@Element('C')
@InWhitelist(145)
@NeighborCount(3)
@NeighborsExactly('C', 2)
@NeighborsExactly('CL', 1)
@Whitelist(263)
@Blacklist(145)
def opls_263(atom):
    """C(Cl) chlorobenzene """
    return True


@Element('CL')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(264)
def opls_264(atom):
    """Cl chlorobenzene """
    return benzene(atom.neighbors[0])


@Element('O')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(278)
def opls_278(atom):
    """AA O: aldehyde """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly(277, 1)
@Whitelist(279)
@Blacklist(140)
def opls_279(atom):
    """AA H-alpha in aldehyde & formamide """
    # TODO: 232 needs to blacklist 277
    #return check_atom(atom.neighbors[0], [232, 277])
    return True


@Element('O')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly('233', 1)
@Whitelist(281)
@Blacklist(278)
def opls_281(atom):
    """AA O: ketone """
    return True


@Element('C')
@NeighborCount(3)
@NeighborsExactly('O', 2)
@NeighborsExactly('C', 1)
@Whitelist(465)
def opls_465(atom):
    """AA C: esters - for R on C=O use #280-#282 """
    return True


@Element('O')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly('465', 1)
@Whitelist(466)
@Blacklist(278)
def opls_466(atom):
    """AA =O: esters """
    return True


@Element('O')
@NeighborCount(2)
@NeighborsExactly('C', 2)
@NeighborsAtLeast('465', 1)
@Whitelist(467)
def opls_467(atom):
    """AA -OR: ester """
    return True


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly('490', 1)
@Whitelist(469)
@Blacklist(140)
def opls_469(atom):
    """methoxy Hs in esters """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 1)
@NeighborsExactly('O', 1)
@NeighborsExactly('467', 1)
@NeighborsExactly('H', 2)
@Whitelist(490)
def opls_490(atom):
    """C(H2OS) ethyl ester """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 3)
@NeighborsExactly('145', 1)
@NeighborsExactly('H', 1)
@Whitelist(515)
@Blacklist(137)
def opls_515(atom):
    """all-atom C: CH, isopropyl benzene """
    return True


@Element('C')
@InWhitelist('145')
@NeighborCount(3)
@NeighborsExactly('C', 3)
@NeighborsExactly('725', 1)
@Whitelist(724)
@Blacklist([145, 221])
def opls_724(atom):
    """C(CF3) trifluoromethylbenzene """
    return True


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 1)
@NeighborsExactly('145', 1)
@NeighborsExactly('F', 3)
@Whitelist(725)
def opls_725(atom):
    """CF3 trifluoromethylbenzene """
    return True


@Element('F')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly('725', 1)
@Whitelist(726)
def opls_726(atom):
    """F trifluoromethylbenzene """
    return True


@Element('N')
@NeighborCount(3)
@NeighborsExactly('O', 2)
@NeighborsExactly('C', 1)
@Whitelist(760)
def opls_760(atom):
    """N in nitro R-NO2 """
    return True


@Element('O')
@NeighborCount(1)
@NeighborsExactly('N', 1)
@Whitelist(761)
def opls_761(atom):
    """O in nitro R-NO2 """
    return True


@Element('C')
@InWhitelist('145')
@NeighborCount(3)
@NeighborsExactly('C', 2)
@NeighborsExactly('N', 1)
@NeighborsExactly('760', 1)
@Whitelist(768)
@Blacklist(145)
def opls_768(atom):
    """C(NO2) nitrobenzene """
    return True


@Element('O')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@Whitelist(771)
@Blacklist(278)
def opls_771(atom):
    """propylene carbonate O """
    if dioxolane13(atom.neighbors[0]):
        for neighbors_neighbor in atom.neighbors[0].neighbors:
            if neighbors_neighbor.kind != 'O':
                return False
        else:
            return True
    return False


@Element('C')
@NeighborCount(3)
@NeighborsExactly('O', 3)
@Whitelist(772)
def opls_772(atom):
    """propylene carbonate C=O """
    return dioxolane13(atom)


@Element('O')
@NeighborCount(2)
@NeighborsExactly('C', 2)
@Whitelist(773)
@Blacklist(467)
def opls_773(atom):
    """propylene carbonate OS """
    return dioxolane13(atom)


@Element('C')
@NeighborCount(4)
@NeighborsExactly('C', 1)
@NeighborsExactly('O', 1)
@NeighborsExactly('H', 2)
@Whitelist(774)
@Blacklist([218, 490])
def opls_774(atom):
    """propylene carbonate C in CH2 """
    return dioxolane13(atom)


@Element('H')
@NeighborCount(1)
@NeighborsExactly('C', 1)
@NeighborsExactly('774', 1)
@Whitelist(777)
@Blacklist(140)
def opls_777(atom):
    """propylene carbonate H in CH2 """
    return True


@Element('C')
@InWhitelist(145)
@NeighborCount(3)
@NeighborsExactly('C', 2)
@NeighborsExactly('N', 1)
@Whitelist(916)
@Blacklist(145)
def opls_916(atom):
    """C(NH2) aniline """
    return True


def get_opls_fn(name):
    """Get the full path to a file used to validate the OPLS-aa atomtyper.

    In the source distribution, these files are in ``opls_validation``.

    Args:
        name (str): Name of the file to load.
    """
    import os
    from pkg_resources import resource_filename
    fn = resource_filename('mbuild',
                           os.path.join('..', 'opls_validation', name))
    if not os.path.exists(fn):
        raise ValueError('Sorry! {} does not exists. If you just '
                         'added it, you\'ll have to re install'.format(fn))
    return fn


if __name__ == "__main__":
    import pdb
    from mbuild.compound import Compound
    from mbuild.examples.methane.methane import Methane
    #from mbuild.examples.ethane.ethane import Ethane

    m = Methane()
    # m = Ethane()
    # m = Compound.load(get_opls_fn('isopropane.pdb'))
    # m = Compound.load(get_opls_fn('cyclohexane.pdb'))
    # m = Compound.load(get_opls_fn('neopentane.pdb'))
    # m = Compound.load(get_opls_fn('benzene.pdb'))
    # m = Compound.load(get_opls_fn('1-propene.pdb'))
    # m = Compound.load(get_opls_fn('biphenyl.pdb'))

    # from mbuild.examples.alkane_monolayer.alkane_monolayer import AlkaneMonolayer
    # m = AlkaneMonolayer(chain_length=3)

    from atomtyper import find_atomtypes
    find_atomtypes(m, forcefield='OPLS-AA')

    for i, atom in enumerate(m.atoms):
        print "Atom kind={}, opls_type={}".format(
            atom.kind, atom.atomtype)
