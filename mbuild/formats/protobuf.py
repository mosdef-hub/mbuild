"""Protocol buffers (Protobuf) format.

A language-agnostic data serialization format developed by Google
https://developers.google.com/protocol-buffers
"""
import ele
import numpy as np
from google.protobuf.text_format import Merge, PrintMessage

from mbuild import Box, Compound
from mbuild.formats import compound_pb2

__all__ = ["write_pb2", "read_pb2"]


def write_pb2(cmpd, filename, binary=True):
    """Convert mb.Compound to Protobuf Message file.

    Parameters
    ----------
    cmpd : mb.Compound
    filename : str
    binary: bool, default True
        If True, will print a binary file
        If False, will print to a text file
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

    if binary:
        with open(filename, "wb") as f:
            f.write(root_proto.SerializeToString())
    else:
        with open(filename, "w") as f:
            PrintMessage(root_proto, f)


def read_pb2(filename, binary=True):
    """Convert a Protobuf Message file into mb.Compound.

    Parameters
    ----------
    filename : str
    binary: bool, default True
        If True, will print a binary file
        If False, will print to a text file

    Returns
    -------
    root_compound : mb.Compound
    """
    root_proto = compound_pb2.Compound()
    if binary:
        with open(filename, "rb") as f:
            root_proto.ParseFromString(f.read())
    else:
        with open(filename, "r") as f:
            Merge(f.read(), root_proto)

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
    """Given mb.Compound, parse propertes into compound_pb2.Compound."""
    proto.name = cmpd.name
    proto.pos.x, proto.pos.y, proto.pos.z = cmpd.pos
    proto.charge = cmpd.charge if cmpd.charge else 0.0
    proto.id = id(cmpd)
    (
        proto.periodicity.x,
        proto.periodicity.y,
        proto.periodicity.z,
    ) = cmpd.periodicity
    if cmpd.element:
        proto.element.name = cmpd.element.name
        proto.element.symbol = cmpd.element.symbol
        proto.element.atomic_number = cmpd.element.atomic_number
        proto.element.mass = cmpd.element.mass

    return proto


def _add_proto_bonds(cmpd, proto):
    """Parse the mb.Compound bonds, add to the proto bonds.

    Parameters
    ----------
    cmpd : mb.Compound
    proto : compound_pb2.Compound
    """
    for b in cmpd.bonds():
        proto_bond = proto.bonds.add()
        proto_bond.id1 = id(b[0])
        proto_bond.id2 = id(b[1])


def _proto_successors(proto):
    """Recurisve method to look for a compound_pb2's children.

    Parameters
    ----------
    proto : compound_pb2

    Notes
    -----
    Base Case: there are no children to the proto, just return
    Recursion: First look at proto's children and return these children
        (sub_proto) Then make the recursive call to look at all the sub_proto's
        successors. This is similar to mb.Compound().successors(). Unlike
        mb.Compound(), we need to also keep track of parents in this recursion
    """
    if len(proto.children) == 0:
        return
    for sub_proto in proto.children:
        yield (sub_proto, proto)
        for sub_sub_proto, parent_proto in _proto_successors(sub_proto):
            yield (sub_sub_proto, parent_proto)


def _proto_to_mb(proto):
    """Given compound_pb2.Compound, create mb.Compound.

    Parameters
    ----------
    proto: compound_pb2.Compound()
    """
    if proto.element.symbol == "":
        elem = None
    else:
        elem = ele.element_from_symbol(proto.element.symbol)
    lengths = [proto.periodicity.x, proto.periodicity.y, proto.periodicity.z]
    if np.all(np.array(lengths) != 0):
        box = Box(lengths)
    else:
        box = None
    return Compound(
        name=proto.name,
        pos=[proto.pos.x, proto.pos.y, proto.pos.z],
        charge=proto.charge,
        box=box,
        element=elem,
    )


def _add_mb_bonds(proto, cmpd, proto_to_cmpd):
    """Parse the compound_pb2.Compound bonds, add to mb.Compound.

    Parameters
    ----------
    proto : compound_pb2.Compound
    cmpd : mb.Compound
    proto_to_cmpd : dict
        keys : compound_pb2.Compound.id
        value : mb.Compound
    """
    for bond in proto.bonds:
        cmpd.add_bond([proto_to_cmpd[bond.id1], proto_to_cmpd[bond.id2]])
