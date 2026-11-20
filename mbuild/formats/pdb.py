"""Read and write mBuild Compounds as PDB files."""

import logging
from dataclasses import dataclass, field
from itertools import groupby
from operator import attrgetter

import numpy as np
from ele import element_from_symbol
from ele.exceptions import ElementError

import mbuild as mb
from mbuild.box import Box
from mbuild.exceptions import MBuildError

__all__ = ["read_pdb", "write_pdb"]

logger = logging.getLogger(__name__)

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


def read_pdb(filename, compound=None, coords_only=False, infer_hierarchy=True):
    """Read a PDB file into a Compound.

    Parameters
    ----------
    filename : str
        Path of the PDB file.
    compound : mb.Compound, optional, default=None
        Compound to add the structure to. A new Compound is made if None.
    coords_only : bool, optional, default=False
        Only set the positions of the particles of ``compound``, in file order.
    infer_hierarchy : bool, optional, default=True
        If True, each residue becomes a Compound. If False, atoms are added
        directly to their chain.

    Returns
    -------
    compound : mb.Compound

    Notes
    -----
    The hierarchy follows ``from_trajectory``: one Compound per chain,
    labeled ``chain[$]``, when the file has more than one chain, and one
    Compound per residue, named after the residue. A chain ends at a change of
    chain ID or a TER record. Residue numbers, insertion codes and chain IDs
    are stored in ``residue_number``, ``insertion_code`` and ``chain_id``.

    Elements are read from columns 77-78, or from the first letter of the atom
    name when those columns are blank, and formal charges from columns 79-80.
    Bonds are read from CONECT records; a bond listed two or three times is a
    double or triple bond. Bonds that are not listed, such as those within
    standard residues, are not added. Only the first model is read, and only
    the first alternate location of each residue is kept.
    """
    with open(filename) as f:
        chains, bond_orders, box = _read_records(f)
    chains = _keep_first_alt_locs(chains)

    if compound is not None and coords_only:
        atoms = [atom for chain in chains for atom in chain]
        particles = list(compound.particles())
        if len(atoms) != len(particles):
            raise ValueError(
                f"Number of atoms in {filename} ({len(atoms)}) does not "
                f"match {compound} ({len(particles)})."
            )
        for particle, atom in zip(particles, atoms):
            particle.pos = atom.position
        return compound

    if compound is None:
        compound = mb.Compound()
    particle_of = {}  # atom serial -> the Particle made for it
    chain_compounds = []
    for chain in chains:
        if len(chains) > 1:
            chain_compound = mb.Compound()
            chain_compounds.append(chain_compound)
        else:
            chain_compound = compound
        chain_compound.chain_id = chain[0].chain_id or None

        residue_compounds = []
        for _, residue_atoms in groupby(chain, key=attrgetter("residue_key")):
            residue_atoms = list(residue_atoms)
            particles = [atom.to_particle() for atom in residue_atoms]
            for atom, particle in zip(residue_atoms, particles):
                particle_of[atom.serial] = particle
            if infer_hierarchy:
                parent = residue_atoms[0].to_residue_compound()
                residue_compounds.append(parent)
            else:
                parent = chain_compound
            parent.add(particles, label=[f"{p.name}[$]" for p in particles])
        if residue_compounds:
            chain_compound.add(residue_compounds)
    if chain_compounds:
        compound.add(chain_compounds, label=["chain[$]"] * len(chain_compounds))

    for (serial1, serial2), bond_order in bond_orders.items():
        if serial1 in particle_of and serial2 in particle_of:
            compound.add_bond(
                (particle_of[serial1], particle_of[serial2]), bond_order=bond_order
            )
    if box is not None:
        compound.box = box
    return compound


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


def _bond_order_from_repeats(repeats):
    """Return the bond order of a bond listed ``repeats`` times.

    The inverse of ``_conect_repeats``.
    """
    return float(repeats) if repeats in (2, 3) else 1.0


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


@dataclass
class _AtomRecord:
    """The fields of one ATOM or HETATM line of a PDB file."""

    serial: int
    name: str
    alt_loc: str
    residue_name: str
    chain_id: str
    residue_number: int
    insertion_code: str
    position: np.ndarray  # nm
    element: object  # an ele Element, or None
    formal_charge: object  # an int, or None

    @classmethod
    def from_line(cls, line):
        """Parse a line, using the columns of ``_atom_line``."""
        return cls(
            serial=_decode_index(line[6:11], width=5),
            name=line[12:16].strip(),
            alt_loc=line[16].strip(),
            residue_name=line[17:21].strip(),
            chain_id=line[21].strip(),
            residue_number=_decode_index(line[22:26], width=4),
            insertion_code=line[26].strip(),
            position=np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
            / 10,
            element=_parse_element(line[76:78], atom_name=line[12:16]),
            formal_charge=_parse_charge(line[78:80]),
        )

    @property
    def residue_key(self):
        """The fields that together identify the atom's residue."""
        return (
            self.chain_id,
            self.residue_number,
            self.insertion_code,
            self.residue_name,
        )

    def to_particle(self):
        """Return a Particle with this atom's name, position and chemistry."""
        return mb.Particle(
            name=self.name,
            pos=self.position,
            element=self.element,
            formal_charge=self.formal_charge,
        )

    def to_residue_compound(self):
        """Return an empty Compound for this atom's residue."""
        residue = mb.Compound(name=self.residue_name)
        residue.residue_number = self.residue_number
        residue.insertion_code = self.insertion_code or None
        return residue


def _read_records(lines):
    """Read the atoms, bonds and box of the first model of a PDB file.

    Returns the atoms as a list of chains, each a list of ``_AtomRecord``,
    the bond orders keyed by pairs of atom serials, and the box (or None).
    A chain ends at a TER record or where the chain ID changes.
    """
    chains = [[]]
    times_listed = {}  # (atom serial, bonded atom serial) -> count
    box = None
    for line in lines:
        record_name = line[:6].strip()
        if record_name == "ENDMDL":
            break
        elif record_name == "CRYST1":
            box = _parse_cryst1(line)
        elif record_name in ("ATOM", "HETATM"):
            atom = _AtomRecord.from_line(line)
            current_chain = chains[-1]
            if current_chain and atom.chain_id != current_chain[-1].chain_id:
                chains.append([])
            chains[-1].append(atom)
        elif record_name == "TER":
            chains.append([])
        elif record_name == "CONECT":
            serial = _decode_index(line[6:11], width=5)
            # Up to four bonded atoms, in columns 12-31
            for start in (11, 16, 21, 26):
                columns = line[start : start + 5]
                if columns.strip():
                    key = (serial, _decode_index(columns, width=5))
                    times_listed[key] = times_listed.get(key, 0) + 1
    chains = [chain for chain in chains if chain]
    return chains, _bond_orders_from_conect(times_listed), box


def _bond_orders_from_conect(times_listed):
    """Return {(serial1, serial2): bond order} from CONECT listings.

    A bond is usually listed under both of its atoms, and sometimes under one
    only, so its order comes from the larger of the two counts.
    """
    repeats = {}
    for (serial1, serial2), count in times_listed.items():
        pair = (min(serial1, serial2), max(serial1, serial2))
        repeats[pair] = max(repeats.get(pair, 0), count)
    return {pair: _bond_order_from_repeats(n) for pair, n in repeats.items()}


def _keep_first_alt_locs(chains):
    """Drop the atoms of alternate locations other than the first of a residue.

    PDB files can hold several conformations of a residue, marked by an
    alternate location ID in column 17. The first ID seen in each residue is
    kept, and a warning says how many residues had others.
    """
    first_alt_loc = {}  # residue key -> the alternate location ID kept
    skipped = set()  # residue keys with another alternate location
    kept_chains = []
    for chain in chains:
        kept = []
        for atom in chain:
            if atom.alt_loc:
                key = atom.residue_key
                if first_alt_loc.setdefault(key, atom.alt_loc) != atom.alt_loc:
                    skipped.add(key)
                    continue
            kept.append(atom)
        kept_chains.append(kept)
    if skipped:
        logger.warning(
            f"Only the first alternate location was read for {len(skipped)} "
            f"residues, e.g. {sorted(skipped)[0]}."
        )
    return kept_chains


def _parse_cryst1(line):
    """Return the box of a CRYST1 record, with lengths in nm."""
    lengths = [float(line[6:15]), float(line[15:24]), float(line[24:33])]
    angles = [float(line[33:40]), float(line[40:47]), float(line[47:54])]
    return Box(lengths=np.array(lengths) / 10, angles=angles)


def _parse_element(symbol_field, atom_name):
    """Return the element in columns 77-78, or else the one the name starts with.

    Atom names begin with their element symbol, possibly after a digit, as
    in " CA " or "1HB ". Returns None when no element is found.
    """
    symbol = symbol_field.strip()
    if not symbol:
        symbol = next((c for c in atom_name if c.isalpha()), "")
    if not symbol:
        return None
    try:
        return element_from_symbol(symbol.capitalize())
    except ElementError:
        return None


def _parse_charge(columns):
    """Return the formal charge in columns 79-80, written like "1-" or "2+"."""
    text = columns.strip()
    if not text:
        return None
    magnitude, sign = int(text[:-1]), text[-1]
    return -magnitude if sign == "-" else magnitude


def _decode_index(columns, width):
    """Read an atom serial or residue number written by ``_encode_index``."""
    text = columns.strip()
    try:
        return int(text)
    except ValueError:
        return int(text, 16) - 10 * 16 ** (width - 1) + 10**width
