"""JSON format."""
import json
from collections import OrderedDict

import ele

import mbuild as mb
from mbuild.exceptions import MBuildError


def compound_from_json(json_file):
    """Convert the given json file into a Compound.

    Given an input json file, this method scans for the particles, bond
    information as well as other hierarchical information regarding the
    compound and returns an mbuild.Compound.

    Parameters
    ----------
    json_file: path or str,
        Path to the json file

    Returns
    -------
    parent: mbuild.Compound,
        the compound equivelent of the json file

    Raises
    ------
    ValueError: This is raised when the JSON file cannot be parsed by python's
        json module
    MBuildError: This is raised on version incompatibility and missing JSON
        keys, when trying to convert the compound to JSON.
    """
    with open(json_file, "r") as cmpdfile:
        try:
            cmpd_dict_and_meta = json.load(cmpdfile)
        except ValueError as e:
            raise e
        try:
            _perform_sanity_check(cmpd_dict_and_meta)
        except MBuildError as e:
            raise e
        compound_dict = cmpd_dict_and_meta["Compound"]
        converted_dict = {}
        parent = _dict_to_mb(compound_dict)
        converted_dict[compound_dict["id"]] = parent
        for sub_compound, compound in _dict_successors(compound_dict):
            if compound["id"] not in converted_dict:
                parent_compound = _dict_to_mb(compound)
                converted_dict[compound["id"]] = parent_compound
            parent_compound = converted_dict[compound["id"]]

            if sub_compound["id"] not in converted_dict:
                sub_cmpd = _dict_to_mb(sub_compound)
                converted_dict[sub_compound["id"]] = sub_cmpd
            sub_cmpd = converted_dict[sub_compound["id"]]

            label_str = sub_compound["label"]
            label_list = compound.get("label_list", {})
            for key, vals in label_list.items():
                if not parent_compound.labels.get(key, None):
                    parent_compound.labels[key] = list()
                if sub_compound["id"] in vals:
                    parent_compound.labels[key].append(sub_cmpd)
            parent_compound.add(sub_cmpd, label=label_str)

        _add_ports(compound_dict, converted_dict)
        _add_bonds(compound_dict, parent, converted_dict)

        box = compound_dict.get("box")
        if box is not None:
            parent.box = mb.Box(lengths=box["lengths"], angles=box["angles"])
        return parent


def compound_to_json(cmpd, file_path, include_ports=False):
    """Convert the Compound into equivelent json representation.

    This method takes an mbuild.Compound and tries to save the hierarchical
    information of the Compound into a json file.

    Parameters
    ----------
    cmpd: mb.Compound,
    file_path: str,
        path to save the JSON file.
    include_ports: bool,
        whether to dump port information, default False
    """
    # Maintain a bookkeeping dict, to do the nesting of children correctly
    version = mb.__version__
    cmpd_info = {}
    compound_dict = _particle_info(cmpd, include_ports)
    cmpd_info[cmpd] = compound_dict

    # Iteratively collect all the information for the children/successors
    for sub_compound in cmpd.successors():
        if not sub_compound.port_particle:
            parent_compound = sub_compound.parent
            sub_compound_dict = _particle_info(sub_compound, include_ports)
            sub_compound_dict["parent_id"] = id(parent_compound)
            sub_compound_dict["is_port"] = False
            sub_compound_dict["label"] = None
            for key, val in sub_compound.parent.labels.items():
                if val == sub_compound:
                    sub_compound_dict["label"] = key
                if isinstance(val, list):
                    if not cmpd_info[sub_compound.parent].get(
                        "label_list", None
                    ):
                        cmpd_info[sub_compound.parent][
                            "label_list"
                        ] = OrderedDict()
                    cmpd_info[sub_compound.parent]["label_list"][key] = [
                        id(x) for x in val
                    ]

            if not cmpd_info[parent_compound].get("children", False):
                cmpd_info[parent_compound]["children"] = list()
            cmpd_info[parent_compound]["children"].append(sub_compound_dict)
            cmpd_info[sub_compound] = sub_compound_dict

    # Should this be nested as well? Not sure...
    compound_dict["bonds"] = _bond_info(cmpd)
    compound_dict["periodicity"] = cmpd.periodicity
    compound_json = OrderedDict()
    compound_json["mbuild-version"] = version
    compound_json["type"] = "Compound"
    compound_json["Compound"] = compound_dict
    with open(file_path, "w") as datafile:
        json.dump(compound_json, datafile, indent=2)


def _particle_info(cmpd, include_ports=False):
    """Return particle information, in a JSON serializable OrderedDict."""
    particle_dict = OrderedDict()
    particle_dict["id"] = id(cmpd)
    particle_dict["name"] = cmpd.name
    particle_dict["pos"] = cmpd.pos.tolist()
    particle_dict["charge"] = cmpd.charge
    particle_dict["element"] = cmpd.element
    particle_dict["box"] = _box_info(cmpd)

    if include_ports:
        particle_dict["ports"] = list()
        for port in cmpd.available_ports():
            port_info = OrderedDict()
            if port.anchor is not None:
                port_info["anchor"] = id(port.anchor)
            else:
                port_info["anchor"] = None
            port_info["label"] = None
            # Is this the most efficient way?
            for key, val in cmpd.labels.items():
                if (val == port) and val.port_particle:
                    port_info["label"] = key
            particle_dict["ports"].append(port_info)
    return particle_dict


def _bond_info(cmpd):
    """Given a compound, return the bond information."""
    bond_list = list()
    for bond in cmpd.bonds():
        bond_list.append((id(bond[0]), id(bond[1])))
    return bond_list


def _box_info(cmpd):
    """Given a compound, return the box information."""
    if cmpd.box:
        box_info = {"lengths": cmpd.box.lengths, "angles": cmpd.box.angles}
    else:
        box_info = None

    return box_info


def _dict_to_mb(compound_dict):
    """Given a dictionary, return the equivelent mb.Compound."""
    name = compound_dict.get("name", "Compound")
    pos = compound_dict.get("pos", [0.0, 0.0, 0.0])
    charge = compound_dict.get("charge", 0.0)
    periodicity = compound_dict.get("periodicity", (False, False, False))
    element = compound_dict.get("element", None)
    box = compound_dict.get("box", None)
    if isinstance(element, ele.element.Element):
        pass
    elif isinstance(element, list):
        atom_num = element[0]
        element = ele.element_from_atomic_number(atom_num)
    elif isinstance(element, str):
        pass
    else:
        pass

    if box is not None:
        box = mb.Box(lengths=box["lengths"], angles=box["angles"])

    this_particle = mb.Compound(
        name=name,
        pos=pos,
        charge=charge,
        periodicity=periodicity,
        element=element,
        box=box,
    )
    return this_particle


def _dict_successors(compound_dict):
    """Get all successors of a compound and its subcompounds recursively.

    Notes
    -----
    This implementation borrows concepts from protobuf.py's _proto_successors()
    """
    if not compound_dict.get("children", False):
        return
    else:
        for sub_compund in compound_dict["children"]:
            yield sub_compund, compound_dict
            for sub_sub_compound, parent_compound in _dict_successors(
                sub_compund
            ):
                yield (sub_sub_compound, parent_compound)


def _add_ports(compound_dict, converted_dict):
    """After adding all particles, this method will add ports if any exist."""
    for subcompound, compound in _dict_successors(compound_dict):
        ports = compound.get("ports", None)
        if ports:
            for port in ports:
                label_str = port["label"]
                port_to_add = mb.Port(anchor=converted_dict[port["anchor"]])
                converted_dict[compound["id"]].add(port_to_add, label_str)
            # Not necessary to add same port twice
            compound["ports"] = None
        ports = subcompound.get("ports", None)
        if ports:
            for port in ports:
                label_str = port["label"]
                port_to_add = mb.Port(anchor=converted_dict[port["anchor"]])
                converted_dict[subcompound["id"]].add(port_to_add, label_str)
            subcompound["ports"] = None


def _add_bonds(compound_dict, parent, converted_dict):
    """Add bonds from the json files to the compound."""
    for bond in compound_dict["bonds"]:
        parent.add_bond(
            particle_pair=(converted_dict[bond[0]], converted_dict[bond[1]])
        )


def _perform_sanity_check(json_dict):
    """Perform Sanity Check on the JSON File."""
    from warnings import warn

    warning_msg = (
        "This Json was written using {0}, current mbuild version is {1}."
    )
    this_version = mb.__version__
    json_mbuild_version = json_dict.get("mbuild-version", None)

    if not json_mbuild_version:
        raise MBuildError(
            "The uploaded JSON file doesn't isn't correctly formatted"
        )
    json_mb_type = json_dict.get("type", None)

    if (not json_mb_type) or (json_mb_type != "Compound"):
        raise MBuildError(
            "Error. Cannot convert JSON of type: {}".format(json_mb_type)
        )

    [major, minor, patch] = json_mbuild_version.split(".")
    [this_major, this_minor, this_patch] = this_version.split(".")
    if major != this_major:
        raise MBuildError(
            warning_msg.format(json_mbuild_version, this_version)
            + " Cannot Convert JSON to compound"
        )
    if minor != this_minor:
        warn(
            warning_msg.format(json_mbuild_version, this_version)
            + " Will Proceed."
        )
