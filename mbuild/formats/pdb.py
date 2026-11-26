"""Write mBuild Compounds to PDB files."""

import logging
from collections import Counter
from dataclasses import dataclass, field
from itertools import groupby
from operator import attrgetter

import numpy as np

from mbuild.exceptions import MBuildError

__all__ = ["PDBResidue", "pdb_residues", "write_pdb"]

logger = logging.getLogger(__name__)

# Residue names written as ATOM records: the standard amino acids and
# nucleotides. Every other residue is written as HETATM records.
STANDARD_RESIDUES = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR "
    "TRP TYR VAL A C G U I DA DC DG DT DI".split()
)


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
        and if that is None too, no CRYST1 record is written.
    include_ports : bool, optional, default=False
        Write the particles of ports, as atoms named ``VS`` in a residue
        named ``PRT``.

    Notes
    -----
    The residues, residue numbers and atom names written are those returned
    by ``pdb_residues``, which describes them in detail. Each chain is closed
    by a TER record. Residues with a standard amino acid or nucleotide name
    are written as ATOM records, and all others as HETATM records.

    CONECT records list every bond, so readers do not have to infer bonds from
    residue names. A double or triple bond is listed two or three times, the
    convention RDKit and Open Babel use for bond orders. The formal charge of
    a particle is written in columns 79-80 when it is set and not zero.

    PDB fields have fixed widths. Atom serials above 99999 and residue
    numbers above 9999 are written in OpenMM's hexadecimal form; mmCIF is
    better suited to systems that large. Coordinates outside -999.999 to
    9999.999 Angstroms and formal charges beyond +/-9 raise an error.
    """
    if box is None:
        box = compound.box
    written = pdb_residues(compound, residues=residues, include_ports=include_ports)

    lines = [] if box is None else [_cryst1_line(box)]
    serial_of = {}  # particle -> its atom serial number in the file
    serial = 0
    for _, chain in groupby(written, key=attrgetter("chain_index")):
        chain = list(chain)
        for residue in chain:
            for particle in residue.particles:
                serial += 1
                serial_of[particle] = serial
                lines.append(_atom_line(serial, particle, residue))
        if chain[-1].chain_id:
            serial += 1  # a TER record takes the next serial number
            lines.append(_ter_line(serial, chain[-1]))
    lines.extend(_conect_lines(compound, serial_of))
    lines.append("END")

    largest_number = max((residue.number for residue in written), default=0)
    if serial > 99999 or largest_number > 9999:
        logger.warning(
            "Atom serials above 99999 or residue numbers above 9999 are "
            "written in OpenMM's hexadecimal form, which OpenMM, ParmEd and "
            "OpenFF read but tools such as RDKit and MDAnalysis do not. "
            "Consider mmCIF for systems this large."
        )

    with open(filename, "w") as f:
        f.write("\n".join(lines) + "\n")


def pdb_residues(compound, residues=None, include_ports=False):
    """Return the residues of a Compound as ``write_pdb`` writes them.

    Use this to find the residue, residue number and atom name that each
    particle is given in the written file. These can differ from the
    Compound, because repeated atom names are numbered and long names are cut
    to four characters.

    Parameters
    ----------
    compound : mb.Compound
        The Compound to describe.
    residues : str or list of str, optional, default=None
        Names of the Compounds to treat as residues, as in ``write_pdb``.
    include_ports : bool, optional, default=False
        Include the particles of ports, as in ``write_pdb``.

    Returns
    -------
    list of PDBResidue
        The residues in the order they are written.

    Notes
    -----
    Each particle belongs to the residue given by its nearest ancestor (or
    itself) with a ``residue_number``; otherwise its nearest ancestor named in
    ``residues``; otherwise a single residue named ``RES``. Each Compound with
    a ``chain_id`` is its own chain. Within a chain, a residue with no
    ``residue_number`` gets the number after the previous residue's, skipping
    numbers that other residues of the chain set explicitly.

    Residue and atom names longer than four characters are cut to four, with
    a warning. Four-character residue names are written in columns 18-21,
    which some readers (RDKit, gemmi, Biopython) do not read. Atom names
    repeated within a residue are numbered (for example ``H`` becomes ``H1``,
    ``H2``, ...), because readers such as OpenMM keep only one atom per name.
    If a numbered name would not fit in four characters, an error is raised;
    this happens when one residue holds many atoms of the same name, usually
    because the atoms were not grouped into residues with ``residues`` or
    ``residue_number``. The Compound itself is not changed.
    """
    if isinstance(residues, str):
        residues = [residues]
    residue_names = set(residues or ())

    # chain Compound (or None) -> {residue Compound, or a default name -> particles}
    grouped = {}
    for particle in compound.particles(include_ports=include_ports):
        if particle.port_particle:
            owner, chain = "PRT", None
        else:
            owner = _find_residue(particle, residue_names) or "RES"
            chain = _find_chain(particle)
        grouped.setdefault(chain, {}).setdefault(owner, []).append(particle)

    written = []
    for chain_index, (chain, owners) in enumerate(grouped.items()):
        chain_id = chain.chain_id if chain is not None else ""
        numbers = _residue_numbers(list(owners), chain_id)
        for (owner, particles), number in zip(owners.items(), numbers):
            name = owner if isinstance(owner, str) else owner.name
            written.append(
                PDBResidue(
                    name=name,
                    number=number,
                    insertion_code=_insertion_code(owner),
                    chain_id=chain_id,
                    chain_index=chain_index,
                    particles=particles,
                )
            )
    _warn_about_names(written)
    return written


@dataclass
class PDBResidue:
    """A residue as ``write_pdb`` writes it.

    Attributes
    ----------
    name : str
        The residue name, cut to four characters.
    number : int
        The residue number.
    insertion_code : str
        The insertion code, or "" if there is none.
    chain_id : str
        The chain ID, or "" if there is none.
    chain_index : int
        The position of the residue's chain in the file, counting from 0.
    particles : list of mb.Compound
        The particles of the residue, in the order they are written.
    atom_names : dict
        The atom name written for each particle, unique within the residue.
    original_name : str
        The residue name before it was cut to four characters.
    is_standard : bool
        Whether the name is a standard amino acid or nucleotide.
    """

    name: str
    number: int
    insertion_code: str
    chain_id: str
    chain_index: int
    particles: list = field(default_factory=list)
    atom_names: dict = field(init=False)
    original_name: str = field(init=False)
    is_standard: bool = field(init=False)

    def __post_init__(self):
        self.original_name = self.name
        self.name = self.name[:4]
        self.atom_names = _unique_atom_names(self.particles, self.name)
        self.is_standard = self.name in STANDARD_RESIDUES

    @property
    def record_name(self):
        """The record name of the residue's atoms, ATOM or HETATM."""
        return "ATOM  " if self.is_standard else "HETATM"


def _residue_numbers(owners, chain_id):
    """Return the residue number of each residue of one chain.

    A residue keeps a ``residue_number`` set on it. Any other residue gets the
    number after the previous residue's, skipping numbers set explicitly in
    the chain. Residues that share an explicit number and insertion code are
    reported, because readers would merge them.
    """
    explicit = [
        (owner.residue_number, _insertion_code(owner))
        for owner in owners
        if not isinstance(owner, str) and owner.residue_number is not None
    ]
    repeated = sorted(key for key, n in Counter(explicit).items() if n > 1)
    if repeated:
        logger.warning(
            f"Residues of chain {chain_id!r} share residue numbers "
            f"{[number for number, _ in repeated]}, so PDB readers may merge "
            "them into one residue."
        )
    taken = {number for number, _ in explicit}
    numbers = []
    previous = 0
    for owner in owners:
        number = None if isinstance(owner, str) else owner.residue_number
        if number is None:
            number = previous + 1
            while number in taken:
                number += 1
            taken.add(number)
        numbers.append(number)
        previous = number
    return numbers


def _insertion_code(owner):
    """Return the insertion code of a residue, or "" if it has none."""
    if isinstance(owner, str):
        return ""
    return owner.insertion_code or ""


def _unique_atom_names(particles, residue_name):
    """Return {particle: name}, numbering names repeated within a residue.

    Names are cut to four characters. A repeated name gets a number appended
    (``H`` -> ``H1``, ``H2``, ...), skipping names already used in the
    residue. Raises if a numbered name would be longer than four characters.
    """
    names = {}
    for particle in particles:
        name = "VS" if particle.port_particle else particle.name
        names[particle] = name[:4]
    counts = Counter(names.values())
    taken = set(names.values())
    next_number = {}
    for particle, name in names.items():
        if counts[name] == 1:
            continue
        number = next_number.get(name, 1)
        while True:
            candidate = f"{name}{number}"
            number += 1
            if candidate not in taken:
                break
        if len(candidate) > 4:
            raise MBuildError(
                f"Residue {residue_name!r} has too many atoms named {name!r} "
                "to give each a unique four-character PDB atom name. Group "
                "the atoms into residues with residues= or residue_number."
            )
        next_number[name] = number
        taken.add(candidate)
        names[particle] = candidate
    return names


def _warn_about_names(residues):
    """Log one warning each for long residue names and renamed atoms."""
    long_names = sorted({r.original_name for r in residues if len(r.original_name) > 3})
    if long_names:
        logger.warning(
            f"Residue names {long_names} are longer than three characters. "
            "Names are written in columns 18-21 and cut to four characters, "
            "and some readers (RDKit, gemmi, Biopython) read only three."
        )
    long_atom_names = sorted(
        {p.name for r in residues for p in r.particles if len(p.name) > 4}
    )
    if long_atom_names:
        logger.warning(
            f"Atom names {long_atom_names[:5]} are longer than four "
            "characters and were cut to four."
        )
    renamed = sum(
        any(
            name != ("VS" if p.port_particle else p.name)
            for p, name in r.atom_names.items()
        )
        for r in residues
    )
    if renamed:
        logger.warning(
            f"Atom names repeated within a residue were numbered in the file "
            f"(for example H -> H1, H2) in {renamed} residues, because PDB "
            "readers such as OpenMM keep only one atom per name."
        )


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


def _find_chain(particle):
    """Return the nearest ancestor (or the particle) with a chain_id, or None."""
    for ancestor in (particle, *particle.ancestors()):
        if ancestor.chain_id is not None:
            if len(ancestor.chain_id) != 1:
                raise MBuildError(
                    f"PDB chain IDs are one character, but {ancestor} has "
                    f"chain_id {ancestor.chain_id!r}."
                )
            return ancestor
    return None


def _conect_repeats(bond_order):
    """Return how many times a bond is listed in CONECT records.

    A double bond is listed twice and a triple bond three times. Other bonds,
    including aromatic and unspecified ones, are listed once.
    """
    return int(bond_order) if bond_order in (2, 3) else 1


def _conect_lines(compound, serial_of):
    """Return the CONECT records of every bond between written particles."""
    bonded = {}  # atom serial -> serials it is bonded to, repeated by order
    for particle1, particle2, data in compound.bonds(return_bond_order=True):
        if particle1 not in serial_of or particle2 not in serial_of:
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
    name = residue.atom_names[particle]
    if particle.port_particle:
        symbol, charge = "", ""
    else:
        symbol = particle.element.symbol.upper() if particle.element else ""
        charge = _charge_field(particle.formal_charge)
    x, y, z = particle.pos * 10.0
    if not all(-999.999 <= value <= 9999.999 for value in (x, y, z)):
        raise MBuildError(
            f"The coordinates of {particle} ({x:.3f}, {y:.3f}, {z:.3f} "
            "Angstroms) do not fit the PDB columns (-999.999 to 9999.999)."
        )
    return (
        f"{residue.record_name}{_encode_index(serial, 5)} "
        f"{_atom_name_field(name, symbol)} "
        f"{residue.name:<4.4s}{residue.chain_id or ' ':1s}"
        f"{_encode_index(residue.number, 4)}{residue.insertion_code or ' ':1s}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          "
        f"{symbol:>2s}{charge:2s}"
    )


def _ter_line(serial, residue):
    """Return the TER record that closes a chain after ``residue``."""
    return (
        f"TER   {_encode_index(serial, 5)}      {residue.name:<4.4s}"
        f"{residue.chain_id}{_encode_index(residue.number, 4)}"
        f"{residue.insertion_code or ' ':1s}"
    )


def _atom_name_field(name, symbol):
    """Return the atom name in columns 13-16.

    Names shorter than four characters that start with a letter and belong
    to a one-letter element start in column 14, so that the element symbol
    lines up with two-letter ones (wwPDB format v3.3). Others, such as "FE"
    or "1HB", start in column 13. This is the rule OpenMM uses.
    """
    if len(name) < 4 and name[:1].isalpha() and len(symbol) < 2:
        return f" {name:<3s}"
    return f"{name:<4s}"


def _charge_field(formal_charge):
    """Return the formal charge as written in columns 79-80, e.g. "1-"."""
    if not formal_charge:
        return ""
    if abs(formal_charge) > 9:
        raise MBuildError(
            f"A formal charge of {formal_charge} does not fit the two PDB "
            "charge columns."
        )
    sign = "+" if formal_charge > 0 else "-"
    return f"{abs(formal_charge)}{sign}"


def _encode_index(value, width):
    """Format an atom serial or residue number to fit its PDB columns.

    Values too large for the columns are written in hexadecimal, offset as
    OpenMM does (``openmm.app.PDBFile``), so 100000 is written as A0000.
    """
    if value <= -(10 ** (width - 1)):
        raise MBuildError(f"{value} is too small for a {width}-column PDB field.")
    if value < 10**width:
        return f"{value:{width}d}"
    limit = 10**width + 6 * 16 ** (width - 1)
    if value >= limit:
        raise MBuildError(f"{value} is too large for a {width}-column PDB field.")
    return f"{value - 10**width + 10 * 16 ** (width - 1):{width}X}"
