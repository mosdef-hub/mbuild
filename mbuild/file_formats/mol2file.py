from mbuild.compound import Compound



def load_mol2(filename, part=None):
    """Load a TRIPOS mol2 file into a Compound.

    If no Compound is specified, this function will return a new Compound.

    Args:
        filename (str): Path to the mol2 file to be read.
        part (Compound, optional): Optionally read the mol2 data into an
            existing Compound.

    Returns:
        part (Compound): The Compound containing the mol2 file's data.

    """
    # coords = []
    # types = []
    # bonds = []
    # bond_types = []
    #
    # atom_list = list()
    # with open(filename, 'r') as mol2_file:
    #     data = dict((key, list(grp)) for key, grp in groupby(mol2_file, _parse_mol2_sections))
    #
    # for idx, atom in enumerate(data['@<TRIPOS>ATOM\n'][1:]):
    #     entries = atom.split()
    #     if len(entries) == 0:
    #         continue
    #     kind = entries[1]
    #     x = entries[2]
    #     y = entries[3]
    #     z = entries[4]
    #     # coords.append(float(x), float(y), float(z))
    #     coords.append(float(x))
    #     coords.append(float(y))
    #     coords.append(float(z))
    #     types.append(kind)
    #
    # # import pdb
    # # pdb.set_trace()
    # coords = np.reshape(coords, newshape=(len(coords)/3,3))
    # types = np.reshape(types, newshape=(len(types)))
    #
    # for bond in data['@<TRIPOS>BOND\n'][1:]:
    #     _, atom1_idx, atom2_idx, _ = bond.split()
    #     bonds.append(int(atom1_idx) - 1)
    #     bonds.append(int(atom2_idx) - 1)
    #
    # bonds = np.reshape(bonds, newshape=(len(bonds)/2, 2))
    #
    # sys = FlatCompound(coords=coords, types=types, bonds=bonds)

    from mbuild.trajectory import load
    t = load(filename)

    if part is None:
        return t
    elif isinstance(part, Compound):
        return t.to_compound(part=part)
    else:
        raise ValueError



def _parse_mol2_sections(x):
    """Helper functon for parsing a section in a mol2 file. """
    if x.startswith('@<TRIPOS>'):
        _parse_mol2_sections.key = x
    return _parse_mol2_sections.key
