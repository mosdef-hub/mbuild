import mbuild as mb
from mbuild.formats import compound_pb2

__all__ = ['write_pb3', 'read_pb3']

def write_pb3(cmpd, filename):
    """ Convert mb.Compound to Protobuff3 file

    Parameters
    ---------
    cmpd : mb.Compound
    filename : str
    """
    cmpd_to_proto = {}

    root_proto = compound_pb2.Compound()
    root_proto = _mb_to_proto(cmpd, root_proto)
    cmpd_to_proto[cmpd] = root_proto

    for sub_cmpd in cmpd.successors():
        parent_cmpd = sub_cmpd.parent
        sub_proto = cmpd_to_proto[parent_cmpd].children.add() 
        sub_proto = _mb_to_proto(sub_cmpd, sub_proto)
        cmpd_to_proto[sub_cmpd] = sub_proto

    _add_proto_bonds(cmpd, root_proto)

    with open(filename, 'wb') as f:
        f.write(root_proto.SerializeToString())

def read_pb3(filename):
    """ Convert a Protobuff3 file into mb.Compound

    Parameters
    ---------
    filename : str

    Returns
    ------
    root_compound : mb.Compound
    """
    root_proto = compound_pb2.Compound()
    with open(filename, 'rb') as f:
        root_proto.ParseFromString(f.read())

    proto_to_cmpd = {}
    root_compound = _proto_to_mb(root_proto)
    proto_to_cmpd[root_proto.id] = root_compound

    for sub_proto, parent_proto in _proto_successors(root_proto):
        if parent_proto.id not in proto_to_cmpd:
            parent_cmpd = _proto_to_mb(parent_proto)
            proto_to_cmpd[parent_proto.id] = parent_cmpd
        parent_cmpd = proto_to_cmpd[parent_proto.id]

        if sub_proto.id not in proto_to_cmpd:
            sub_cmpd = _proto_to_mb(sub_proto)
            proto_to_cmpd[sub_proto.id] = sub_cmpd
        sub_cmpd = proto_to_cmpd[sub_proto.id]

        parent_cmpd.add(sub_cmpd)

    _add_mb_bonds(root_proto, root_compound, proto_to_cmpd)
    return root_compound

def _mb_to_proto(cmpd, proto):
    """ Given mb.Compound, parse propertes into compound_pb2.Compound"""
    proto.name = cmpd.name
    proto.pos.x, proto.pos.y, proto.pos.z = cmpd.pos
    proto.charge = cmpd.charge
    proto.id = id(cmpd)
    proto.periodicity.x, proto.periodicity.y, proto.periodicity.z = cmpd.periodicity
   
    return proto

def _add_proto_bonds(cmpd, proto):
    """ Parse the mb.Compound bonds, add to the proto bonds

    Parameters
    ---------
    cmpd : mb.Compound
    proto : compound_pb2.Compound

    """
    for b in cmpd.bonds():
        proto_bond = proto.bonds.add()
        proto_bond.id1 = id(b[0])
        proto_bond.id2 = id(b[1])


def _proto_successors(proto):
    """ Recurisve method to look for a compound_pb2's children

    Parameters
    ---------
    proto : compound_pb2

    Notes
    -----
    Base Case: there are no children to the proto, just return 
    Recursion: First look at proto's children and return these children (sub_proto)
        Then make the recursive call to look at all the sub_proto's successors
    This is similar to mb.Compound().successors()
    Unlike mb.Compound(), we need to also keep track of parents in this recursion
    """
    if len(proto.children) == 0:
        return
    for sub_proto in proto.children:
        yield (sub_proto, proto)
        for sub_sub_proto, parent_proto in _proto_successors(sub_proto):
            yield (sub_sub_proto, parent_proto)

def _proto_to_mb(proto):
    """ Given compound_pb2.Compound, create mb.Compound 
    
    Parameters
    ----------
    proto: compound_pb2.Compound()
    
    """
    return  mb.Compound(name=proto.name,
                pos=[proto.pos.x, proto.pos.y, proto.pos.z],
                charge=proto.charge,
                periodicity=[proto.periodicity.x, proto.periodicity.y, 
                            proto.periodicity.z])

def _add_mb_bonds(proto, cmpd, proto_to_cmpd):
    """ Parse the compound_pb2.Compound bonds, add to mb.Compound

    Parameters
    ---------
    proto : compound_pb2.Compound
    cmpd : mb.Compound
    proto_to_cmpd : dict
        keys : compound_pb2.Compound.id
        value : mb.Compound
    """
    for bond in proto.bonds:
        cmpd.add_bond([
            proto_to_cmpd[bond.id1],
            proto_to_cmpd[bond.id2]
            ])


