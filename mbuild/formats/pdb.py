"""Write mBuild Compounds to PDB files."""

from dataclasses import dataclass, field

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
        box = compound.box
    if box is None:
        bounding_box = compound.get_boundingbox()
        lengths = np.array(bounding_box.lengths) + 0.5
        box = Box(lengths=lengths, angles=bounding_box.angles)

    lines = [_cryst1_line(box)]
    serial_of = {}  # particle -> its atom serial number in the file
    residue_of = {}  # particle -> the _Residue it is written in
    serial = 0
    for chain in _chains_to_write(compound, residue_names, include_ports):
        for residue in chain:
            for particle in residue.particles:
                serial += 1
                serial_of[particle] = serial
                residue_of[particle] = residue
                lines.append(_atom_line(serial, particle, residue))
        last_residue = chain[-1]
        if last_residue.chain_id:
            serial += 1  # a TER record takes the next serial number
            lines.append(_ter_line(serial, last_residue))
    lines.extend(_conect_lines(compound, serial_of, residue_of))
    lines.append("END")

    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


@dataclass
class _Residue:
    """A residue as it is written to a PDB file."""

    name: str
    number: int
    insertion_code: str
    chain_id: str
    index_in_chain: int
    particles: list = field(default_factory=list)
    # Whether readers know the residue's bonds from its name; set from name
    is_standard: bool = field(init=False)

    def __post_init__(self):
        self.is_standard = self.name in STANDARD_RESIDUES

    @property
    def record_name(self):
        """The record name of the residue's atoms, ATOM or HETATM."""
        return "ATOM  " if self.is_standard else "HETATM"


def _chains_to_write(compound, residue_names, include_ports):
    """Group the particles into chains of residues, in the order written.

    Returns a list of chains, each a list of ``_Residue``.
    """
    # chain ID -> {residue Compound, or a name for a default residue -> particles}
    grouped = {}
    for particle in compound.particles(include_ports=include_ports):
        if particle.port_particle:
            owner, chain_id = "PRT", ""
        else:
            owner = _find_residue(particle, residue_names) or "RES"
            chain_id = _find_chain_id(particle)
        grouped.setdefault(chain_id, {}).setdefault(owner, []).append(particle)

    chains = []
    residue_count = 0
    for chain_id, residues in grouped.items():
        chain = []
        for index, (owner, particles) in enumerate(residues.items()):
            residue_count += 1
            if isinstance(owner, str):
                name, number, insertion_code = owner, residue_count, " "
            else:
                name = owner.name
                number = owner.residue_number
                if number is None:
                    number = residue_count
                insertion_code = owner.insertion_code or " "
            chain.append(
                _Residue(name, number, insertion_code, chain_id, index, particles)
            )
        chains.append(chain)
    return chains


def _find_residue(particle, residue_names):
    """Return the Compound a particle is written in as a residue, or None.

    That is the nearest ancestor (or the particle) with a ``residue_number``,
    or else the nearest one named in ``residue_names``.
    """
    lineage = [particle, *particle.ancestors()]
    for ancestor in lineage:
        if ancestor.residue_number is not None:
            return ancestor
    for ancestor in lineage:
        if ancestor.name in residue_names:
            return ancestor
    return None


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


def _needs_conect(particle1, particle2, residue_of):
    """Return False if PDB readers can infer this bond from residue names.

    Readers know the bonds within a standard residue, and the peptide or
    nucleic acid backbone bond between consecutive standard residues of a
    chain. Every other bond needs a CONECT record.
    """
    residue1 = residue_of[particle1]
    residue2 = residue_of[particle2]
    if not (residue1.is_standard and residue2.is_standard):
        return True
    if residue1 is residue2:
        return False
    if residue1.chain_id != residue2.chain_id:
        return True
    if abs(residue1.index_in_chain - residue2.index_in_chain) != 1:
        return True
    if residue1.index_in_chain < residue2.index_in_chain:
        earlier, later = particle1, particle2
    else:
        earlier, later = particle2, particle1
    return (earlier.name, later.name) not in _BACKBONE_LINKS


def _conect_repeats(bond_order):
    """Return how many times a bond is listed in CONECT records.

    A double bond is listed twice and a triple bond three times. Other bonds,
    including aromatic and unspecified ones, are listed once.
    """
    return int(bond_order) if bond_order in (2, 3) else 1


def _conect_lines(compound, serial_of, residue_of):
    """Return the CONECT records of the bonds PDB readers cannot infer."""
    bonded = {}  # atom serial -> serials it is bonded to, repeated by order
    for particle1, particle2, data in compound.bonds(return_bond_order=True):
        if particle1 not in serial_of or particle2 not in serial_of:
            continue
        if not _needs_conect(particle1, particle2, residue_of):
            continue
        repeats = _conect_repeats(data["bond_order"])
        serial1, serial2 = serial_of[particle1], serial_of[particle2]
        bonded.setdefault(serial1, []).extend([serial2] * repeats)
        bonded.setdefault(serial2, []).extend([serial1] * repeats)

    lines = []
    for serial in sorted(bonded):
        partners = sorted(bonded[serial])
        # One CONECT record has room for four bonded atoms
        for start in range(0, len(partners), 4):
            fields = "".join(_encode_index(p, 5) for p in partners[start : start + 4])
            lines.append(f"CONECT{_encode_index(serial, 5)}{fields}")
    return lines


def _cryst1_line(box):
    """Return the CRYST1 record of a box, with lengths in Angstroms."""
    a, b, c = np.asarray(box.lengths) * 10.0
    alpha, beta, gamma = box.angles
    return (
        f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
    )


def _atom_line(serial, particle, residue):
    """Return the ATOM or HETATM record of one particle.

    The columns, from the wwPDB format v3.3, are: 1-6 record name, 7-11 atom
    serial, 13-16 atom name, 18-21 residue name, 22 chain ID, 23-26 residue
    number, 27 insertion code, 31-54 x, y and z in Angstroms, 55-60
    occupancy, 61-66 temperature factor, 77-78 element and 79-80 charge.
    """
    if particle.port_particle:
        name, symbol, charge = "VS", "", ""
    else:
        name = particle.name
        symbol = particle.element.symbol.upper() if particle.element else ""
        charge = _charge_field(particle.formal_charge)
    x, y, z = particle.pos * 10.0
    return (
        f"{residue.record_name}{_encode_index(serial, 5)} "
        f"{_atom_name_field(name, symbol)} "
        f"{residue.name:<4.4s}{residue.chain_id or ' ':1s}"
        f"{_encode_index(residue.number, 4)}{residue.insertion_code}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
        f"{symbol:>2s}{charge:2s}"
    )


def _ter_line(serial, residue):
    """Return the TER record that closes a chain after ``residue``."""
    return (
        f"TER   {_encode_index(serial, 5)}      {residue.name:<4.4s}"
        f"{residue.chain_id}{_encode_index(residue.number, 4)}"
        f"{residue.insertion_code}"
    )


def _atom_name_field(name, symbol):
    """Return the atom name in columns 13-16.

    Names of atoms with one-letter elements start in column 14, so that
    the element symbol lines up with the two-letter ones (wwPDB format v3.3).
    """
    if len(name) >= 4 or len(symbol) == 2:
        return f"{name[:4]:<4s}"
    return f" {name:<3s}"


def _charge_field(formal_charge):
    """Return the formal charge as written in columns 79-80, e.g. "1-"."""
    if not formal_charge:
        return ""
    sign = "+" if formal_charge > 0 else "-"
    return f"{abs(formal_charge)}{sign}"


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
