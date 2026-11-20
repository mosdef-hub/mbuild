"""Write mBuild Compounds to PDB files."""

from copy import deepcopy

import numpy as np

from mbuild.box import Box
from mbuild.exceptions import MBuildError

__all__ = ["write_pdb"]

# Residue names written as ATOM records: the standard amino acids and
# nucleotides. Every other residue is written as HETATM records.
STANDARD_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
    "TRP TYR VAL A C G U I DA DC DG DT DI".split()
)

# Bonds between consecutive standard residues of a chain that PDB readers
# infer from the record order: the peptide bond and the nucleic acid
# backbone link, as (atom of the earlier residue, atom of the later one).
_BACKBONE_LINKS = {("C", "N"), ("O3'", "P")}


def write_pdb(compound, filename, residues=None, box=None, include_ports=False):
    """Write a Compound to a PDB file.

    Parameters
    ----------
    compound : mb.Compound
        The Compound to write.
    filename : str
        Path of the file to write.
    residues : str or list of str, optional, default=None
        Names of the Compounds to write as residues, used for particles with
        no ``residue_number`` set on them or their ancestors.
    box : mb.Box, optional, default=None
        Box written to the CRYST1 record. If None, ``compound.box`` is used,
        or else the bounding box padded by 0.25 nm on each side.
    include_ports : bool, optional, default=False
        Write the particles of ports, as atoms named ``VS`` in a residue
        named ``PRT``.

    Notes
    -----
    Each particle is written in the residue given by its nearest ancestor
    (or itself) with a ``residue_number``; otherwise its nearest ancestor
    named in ``residues``; otherwise a single residue named ``RES``. Residues
    with no ``residue_number`` are numbered in the order they are written.
    The chain ID comes from the nearest ancestor with a ``chain_id``, and a
    TER record closes each chain.

    Residues with a standard amino acid or nucleotide name are written as
    ATOM records, and all others as HETATM records. CONECT records list every
    bond except the bonds within a standard residue and the backbone bonds
    between consecutive standard residues of a chain, which PDB readers take
    from the residue definitions (wwPDB format v3.3, section 10). A double
    or triple bond is listed two or three times, the convention RDKit and
    Open Babel use for bond orders. The formal charge of a particle is
    written in columns 79-80 when it is set and not zero.
    """
    if isinstance(residues, str):
        residues = [residues]
    residue_names = set(residues or ())

    if box is None:
        if compound.box is not None:
            box = deepcopy(compound.box)
        else:
            box = compound.get_boundingbox()
            box = Box(lengths=np.array(box.lengths) + 0.5, angles=box.angles)

    chains = _group_by_chain_and_residue(compound, residue_names, include_ports)

    lines = [_cryst1_line(box)]
    serials = {}
    position = {}
    serial = 0
    residue_count = 0
    for chain_id, chain_residues in chains.items():
        for index, (residue, particles) in enumerate(chain_residues.items()):
            residue_count += 1
            name, number, icode = _residue_fields(residue, residue_count)
            record = "ATOM  " if name in STANDARD_RESIDUES else "HETATM"
            position[chain_id, residue] = (index, name in STANDARD_RESIDUES)
            for particle in particles:
                serial += 1
                serials[particle] = serial
                lines.append(
                    _atom_line(record, serial, particle, name, chain_id, number, icode)
                )
        if chain_id:
            serial += 1
            lines.append(
                f"TER   {_encode_index(serial, 5)}      {name:<4.4s}"
                f"{chain_id}{_encode_index(number, 4)}{icode}"
            )

    residue_of = {
        particle: (chain_id, residue)
        for chain_id, chain_residues in chains.items()
        for residue, particles in chain_residues.items()
        for particle in particles
    }
    partners = {}
    for particle1, particle2, data in compound.bonds(return_bond_order=True):
        if particle1 not in serials or particle2 not in serials:
            continue
        if _implied_by_residues(particle1, particle2, residue_of, position):
            continue
        repeats = int(data["bond_order"]) if data["bond_order"] in (2, 3) else 1
        serial1, serial2 = serials[particle1], serials[particle2]
        partners.setdefault(serial1, []).extend([serial2] * repeats)
        partners.setdefault(serial2, []).extend([serial1] * repeats)
    for atom_serial in sorted(partners):
        bonded = sorted(partners[atom_serial])
        for start in range(0, len(bonded), 4):
            lines.append(
                "CONECT"
                + _encode_index(atom_serial, 5)
                + "".join(_encode_index(s, 5) for s in bonded[start : start + 4])
            )
    lines.append("END")

    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


def _group_by_chain_and_residue(compound, residue_names, include_ports):
    """Return {chain_id: {residue: [particles]}} in the order written."""
    chains = {}
    for particle in compound.particles(include_ports=include_ports):
        if particle.port_particle:
            residue, chain_id = "PRT", ""
        else:
            residue = _find_residue(particle, residue_names)
            chain_id = _find_chain_id(particle)
        chains.setdefault(chain_id, {}).setdefault(residue, []).append(particle)
    return chains


def _find_residue(particle, residue_names):
    """Return the Compound that is the residue of a particle, or "RES"."""
    lineage = [particle, *particle.ancestors()]
    for ancestor in lineage:
        if ancestor.residue_number is not None:
            return ancestor
    for ancestor in lineage:
        if ancestor.name in residue_names:
            return ancestor
    return "RES"


def _find_chain_id(particle):
    """Return the chain ID of the nearest ancestor that has one, or ""."""
    for ancestor in (particle, *particle.ancestors()):
        if ancestor.chain_id is not None:
            if len(ancestor.chain_id) != 1:
                raise MBuildError(
                    f"PDB chain IDs are one character, but {ancestor} has "
                    f"chain_id {ancestor.chain_id!r}."
                )
            return ancestor.chain_id
    return ""


def _residue_fields(residue, count):
    """Return the name, number and insertion code written for a residue."""
    if isinstance(residue, str):
        return residue, count, " "
    number = residue.residue_number
    return (
        residue.name,
        number if number is not None else count,
        (residue.insertion_code or " "),
    )


def _implied_by_residues(particle1, particle2, residue_of, position):
    """Return True if a PDB reader infers this bond from residue definitions."""
    key1, key2 = residue_of[particle1], residue_of[particle2]
    index1, standard1 = position[key1]
    index2, standard2 = position[key2]
    if not (standard1 and standard2):
        return False
    if key1 == key2:
        return True
    if key1[0] != key2[0] or abs(index1 - index2) != 1:
        return False
    earlier, later = (
        (particle1, particle2) if index1 < index2 else (particle2, particle1)
    )
    return (earlier.name, later.name) in _BACKBONE_LINKS


def _cryst1_line(box):
    """Return the CRYST1 record of a box, with lengths in Angstroms."""
    a, b, c = np.asarray(box.lengths) * 10.0
    alpha, beta, gamma = box.angles
    return (
        f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
    )


def _atom_line(record, serial, particle, resname, chain_id, resnum, icode):
    """Return the ATOM or HETATM record of one particle."""
    if particle.port_particle:
        name, symbol, charge = "VS", "", ""
    else:
        name = particle.name
        symbol = particle.element.symbol.upper() if particle.element else ""
        charge = _charge_field(particle.formal_charge)
    # Atom names of one-letter elements start in column 14 (wwPDB format v3.3)
    if len(name) >= 4 or len(symbol) == 2:
        name_field = f"{name[:4]:<4s}"
    else:
        name_field = f" {name:<3s}"
    x, y, z = particle.pos * 10.0
    return (
        f"{record}{_encode_index(serial, 5)} {name_field} {resname:<4.4s}"
        f"{chain_id or ' ':1s}{_encode_index(resnum, 4)}{icode}"
        f"   {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}"
        f"          {symbol:>2s}{charge:2s}"
    )


def _charge_field(formal_charge):
    """Return the formal charge as written in columns 79-80, e.g. "1-"."""
    if not formal_charge:
        return ""
    return f"{abs(formal_charge)}{'+' if formal_charge > 0 else '-'}"


def _encode_index(value, width):
    """Format an atom serial or residue number to fit its PDB columns.

    Values too large for the columns are written in hexadecimal, offset as
    OpenMM does (``openmm.app.PDBFile``), so 100000 is written as A0000.
    """
    if value < 10**width:
        return f"{value:{width}d}"
    limit = 10**width + 6 * 16 ** (width - 1)
    if value >= limit:
        raise MBuildError(f"{value} is too large for a {width}-column PDB field.")
    return f"{value - 10**width + 10 * 16 ** (width - 1):{width}X}"
