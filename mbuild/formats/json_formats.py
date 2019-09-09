import json
from collections import OrderedDict


def compound_from_json(json_file):
    """
    Convert the given json file into a mb.Compound

    Given an input json file, this method scans for the particles, bond information
    as well as other hierarchical information regarding the compound and returns a
    mb.Compound.

    Parameters
    -----------
    json_file: (path, str) Path of the json file

    Returns
    -------
    parent: mb.Compound, the compound equivelent of the json file
    """
    with open(json_file, 'r') as cmpdfile:
        converted_dict = {}
        compound_dict = json.load(cmpdfile)
        parent = _dict_to_mb(compound_dict)
        converted_dict[compound_dict['id']] = parent
        for sub_compound, compound in _dict_successors(compound_dict):
            if compound['id'] not in converted_dict:
                parent_compound = _dict_to_mb(compound)
                converted_dict[compound['id']] = parent_compound
            parent_compound = converted_dict[compound['id']]

            if sub_compound['id'] not in converted_dict:
                sub_cmpd = _dict_to_mb(sub_compound)
                converted_dict[sub_compound['id']] = sub_cmpd
            sub_cmpd = converted_dict[sub_compound['id']]

            label_str = sub_compound['label']
            parent_compound.add(sub_cmpd, label=label_str)
        _add_bonds(compound_dict, parent, converted_dict)

        return parent


def compound_to_json(cmpd, file_path, include_ports=False):
    """Convert the mb.Compound into equivelent json representation

    This method takes in the mb.Compound and tries to save the hierarchical
    information of the mb.Compound into a json file.
    Parameters
    ----------
    cmpd: mb.Compound
    file_path: str, path to save the JSON file.
    include_ports: bool, whether to dump port information, default False

    Raises
    ------

    """
    # Maintain a bookkeeping dict, to do the nesting of children correctly
    cmpd_info = {}
    compound_dict = _particle_info(cmpd)
    cmpd_info[cmpd] = compound_dict

    # Iteratively collect all the information for the children/successors
    for sub_compound in cmpd.successors():
        if not sub_compound.port_particle:
            parent_compound = sub_compound.parent
            sub_compound_dict = _particle_info(sub_compound, include_ports)
            sub_compound_dict['parent_id'] = id(parent_compound)
            sub_compound_dict['is_port'] = False
            sub_compound_dict['label'] = None
            for key, val in cmpd.labels.items():
                if val == sub_compound:
                    sub_compound_dict['label'] = key
            if not cmpd_info[parent_compound].get('children', False):
                cmpd_info[parent_compound]['children'] = list()
            cmpd_info[parent_compound]['children'].append(sub_compound_dict)
            cmpd_info[sub_compound] = sub_compound_dict

    # Should this be nested as well? Not sure...
    compound_dict['bonds'] = _bond_info(cmpd)

    with open(file_path, 'w') as datafile:
        json.dump(compound_dict, datafile, indent=2)


def _particle_info(cmpd, include_ports=False):
    """Return information about a particle, in a JSON serializable OrderdDict"""
    particle_dict = OrderedDict()
    particle_dict['id'] = id(cmpd)
    particle_dict['name'] = cmpd.name
    particle_dict['pos'] = list(cmpd.pos)
    particle_dict['charge'] = cmpd.charge
    particle_dict['periodicity'] = list(cmpd.periodicity)

    if include_ports:
        particle_dict['ports'] = list()
        for port in cmpd.referenced_ports():
            port_info = OrderedDict()
            if port.anchor is not None:
                port_info['anchor'] = id(port.anchor)
            else:
                port_info['anchor'] = None
            port_info['used'] = port.used
            port_info['label'] = None
            # Is this the most efficient way?
            for key, val in cmpd.labels.items():
                if val == port:
                    port_info['label'] = key
            particle_dict['ports'].append(port_info)
    return particle_dict


def _bond_info(cmpd):
    """Given a compound, return the bond information"""
    bond_list = list()
    for bond in cmpd.bonds():
        bond_list.append((id(bond[0]), id(bond[1])))
    return bond_list


def _dict_to_mb(compound_dict):
    """Given a dictionary, return the equivelent mb.Compound."""
    import mbuild as mb
    name = compound_dict.get('name', "Compound")
    pos = compound_dict.get('pos', [0.0, 0.0, 0.0])
    charge = compound_dict.get('charge', 0.0)
    periodicity = compound_dict.get('periodicity', [0.0, 0.0, 0.0])
    this_particle = mb.Compound(name=name, pos=pos, charge=charge, periodicity=periodicity)
    ports = compound_dict.get('ports', None)
    if ports:
        for port in ports:
            label_str = port['label']
            port_to_add = mb.Port(anchor=this_particle)
            port_to_add.used = port['used']
            this_particle.add(port_to_add, label_str)
    return this_particle


def _dict_successors(compound_dict):
    """This is a recursive method to get all successors of a given compound and its subcompounds

    Notes
    -----
        This implementation burrows concept form protobuf.py's _proto_successors()
    """
    if not compound_dict.get('children', False):
        return
    else:
        for sub_compund in compound_dict['children']:
            yield sub_compund, compound_dict
            for sub_sub_compound, parent_compound in _dict_successors(sub_compund):
                yield (sub_sub_compound, parent_compound)


def _add_bonds(compound_dict, parent, converted_dict):
    """Add bonds from the json files to the compound"""
    for bond in compound_dict['bonds']:
        parent.add_bond(particle_pair=(converted_dict[bond[0]], converted_dict[bond[1]]))


