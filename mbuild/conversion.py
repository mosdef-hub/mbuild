"""Module for handling conversions in mBuild."""

import os
import sys
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from warnings import warn

import numpy as np
import parmed as pmd
from ele import (
    element_from_atomic_number,
    element_from_name,
    element_from_symbol,
)
from ele.exceptions import ElementError

import mbuild as mb
from mbuild.box import Box
from mbuild.exceptions import MBuildError
from mbuild.formats.gsdwriter import write_gsd
from mbuild.formats.hoomdxml import write_hoomdxml
from mbuild.formats.json_formats import compound_from_json, compound_to_json
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.formats.par_writer import write_par
from mbuild.formats.xyz import read_xyz, write_xyz
from mbuild.utils.io import (
    has_gmso,
    has_mdtraj,
    has_networkx,
    has_openbabel,
    has_rdkit,
    import_,
)


def load(
    filename_or_object,
    relative_to_module=None,
    compound=None,
    coords_only=False,
    rigid=False,
    smiles=False,
    infer_hierarchy=True,
    backend=None,
    ignore_box_warn=False,
    **kwargs,
):
    """Load a file or an existing topology into an mbuild Compound.

    Files are read using the predefined backend, unless otherwise specified by
    the user (through the `backend` flag).
    Supported backends include "pybel", "mdtraj", "parmed", "rdkit", and "internal".
    Please refer to http://mdtraj.org/1.8.0/load_functions.html for formats
    supported by MDTraj and https://parmed.github.io/ParmEd/html/readwrite.html
    for formats supported by ParmEd.

    Parameters
    ----------
    filename_or_object : str, mdtraj.Trajectory, parmed.Structure,
        mbuild.Compound, pybel.Molecule,
        Name of the file or topology from which to load atom and bond
        information.
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory, look for the file
        where this module is defined. This is typically used in Compound
        classes that will be instantiated from a different directory (such as
        the Compounds located in mbuild.lib).
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information into. New structure
        will be added to the existing compound as a sub compound.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    backend : str, optional, default=None
        Backend used to load structure from file or string. If not specified, a
        default backend (extension specific) will be used.
    smiles: bool, optional, default=False
        Use RDKit or OpenBabel to parse filename as a SMILES string or file
        containing a SMILES string. If this is set to True, `rdkit` is the
        default backend.
    infer_hierarchy : bool, optional, default=True
        If True, infer hierarchy from chains and residues
    ignore_box_warn : bool, optional, default=False
        If True, ignore warning if no box is present. Defaults to True when
        loading from SMILES
    ``**kwargs`` : keyword arguments
        Key word arguments passed to mdTraj, GMSO, RDKit, or pybel for loading.

    Returns
    -------
    compound : mb.Compound

    Notes
    -----
    If `smiles` is `True`, either `rdkit` (default) or `pybel` can be used, but
    RDkit is the only option of these that allows the user to specify a random
    number seed to reproducibly generate the same starting structure. This is
    NOT possible with `openbabel`, use `rdkit` if you need control over starting
    structure's position (recommended).
    """
    # First check if we are loading from an object
    if not isinstance(filename_or_object, str):
        return load_object(
            obj=filename_or_object,
            compound=compound,
            coords_only=coords_only,
            rigid=rigid,
            infer_hierarchy=infer_hierarchy,
            **kwargs,
        )
    # Second check if we are loading SMILES strings
    elif smiles:
        # Ignore the box info for SMILES (its never there)
        ignore_box_warn = True
        return load_smiles(
            smiles_or_filename=filename_or_object,
            compound=compound,
            infer_hierarchy=infer_hierarchy,
            ignore_box_warn=ignore_box_warn,
            backend=backend,
            **kwargs,
        )
    # Last, if none of the above, load from file
    else:
        return load_file(
            filename=filename_or_object,
            relative_to_module=relative_to_module,
            compound=compound,
            coords_only=coords_only,
            rigid=rigid,
            backend=backend,
            infer_hierarchy=infer_hierarchy,
            **kwargs,
        )


def load_object(
    obj,
    compound=None,
    coords_only=False,
    rigid=False,
    infer_hierarchy=True,
    **kwargs,
):
    """Load an obj into a mb.Compound.

    Functions to load on-disk obj to mb.Compound, supporting conversion from
    mb.Compound, pmd.Structure, md.Trajectory, and pybel.Molecule. If the
    compound flag is supplied, this will add the obj as a subcompound of that
    compound.

    Parameters
    ----------
    obj : mb.Compound, pmd.Structure, md.Trajectory, pybel.Molecule
        Molecular object to be converted.
    compound : mb.Compound, optional, default=None
        The host mbuild Compound
    coords_only : bool, optional, default=False
        Only load the coordinates into existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    infer_hierarchy : bool, optional, default=True
        If True, infer hiereachy from chains and residues
    **kwargs : keyword arguments
        Keyword arguments passed to mdTraj for loading

    Returns
    -------
    mb.Compound
    """
    # Create type_dict type -> loading method
    # Will need to add a gmso method soon
    type_dict = {pmd.Structure: from_parmed}
    if has_gmso:
        gmso = import_("gmso")
        type_dict.update({gmso.Topology: from_gmso})

    if has_openbabel:
        pybel = import_("pybel")
        type_dict.update({pybel.Molecule: from_pybel})

    if has_mdtraj:
        md = import_("mdtraj")
        type_dict.update({md.Trajectory: from_trajectory})

    # Check if the given object is an mb.Compound
    if isinstance(obj, mb.Compound):
        if not compound:
            warn("Given object is already an mb.Compound, doing nothing.")
            return obj
        else:
            warn("Given object is an mb.Compound, adding to the host compound.")
            compound.add(obj)
            return compound

    for type_ in type_dict:
        if isinstance(obj, type_):
            compound = type_dict[type_](
                obj,
                compound,
                coords_only=coords_only,
                infer_hierarchy=infer_hierarchy,
                **kwargs,
            )
            if rigid:
                compound.label_rigid_bodies()
            return compound

    # If nothing is return raise an error
    raise ValueError(f"Object of type {type(obj).__name__} is not supported")


def load_smiles(
    smiles_or_filename,
    compound=None,
    infer_hierarchy=True,
    ignore_box_warn=False,
    backend="rdkit",
    coords_only=False,
    **kwargs,
):
    """Load a SMILES string as an mBuild Compound.

    Loading SMILES string from a string, a list, or a file using RDKit by
    default. Must have rdkit or pybel packages installed.

    Parameters
    ----------
    smiles_or_filename : str
        SMILES string or file of SMILES string to load
    compound : mb.Compound
        The host mbuild Compound
    infer_hierarchy : bool, optional, default=True
    ignore_box_warn : bool, optional, default=False
        If True, ignore warning if no box is present.
    coords_only : bool, optional, default=False
        Only load the coordinates into a provided compound.
    backend : str, optional, default='rdkit'
        The smiles loading backend, either 'rdkit' or 'pybel'

    Returns
    -------
    compound : mb.Compound
    """
    # Initialize an mb.Compound if none is provided
    if not compound:
        compound = mb.Compound()

    test_path = Path(smiles_or_filename)

    # Will try to support list of smiles strings in the future
    if backend is None:
        backend = "rdkit"

    if backend == "rdkit":
        rdkit = import_("rdkit")
        from rdkit import Chem
        from rdkit.Chem import AllChem

        if test_path.exists():
            # assuming this is a smi file now
            mymol = Chem.SmilesMolSupplier(smiles_or_filename)
            if not mymol:
                raise ValueError(
                    "Provided smiles string or file was invalid. Refer to the "
                    "above RDKit error messages for additional information."
                )
            mol_list = [mol for mol in mymol]
            if len(mol_list) == 1:
                rdmol = mymol[0]
            else:
                rdmol = mymol[0]
                warn(
                    "More than one SMILES string in file, more than one SMILES "
                    f"string is not supported, using {Chem.MolToSmiles(rdmol)}"
                )
        else:
            rdmol = Chem.MolFromSmiles(smiles_or_filename)

        seed = kwargs.get("smiles_seed", 0)

        return from_rdkit(
            rdkit_mol=rdmol,
            compound=compound,
            coords_only=coords_only,
            smiles_seed=seed,
        )
    elif backend == "pybel":
        pybel = import_("pybel")
        # First we try treating filename_or_object as a SMILES string
        try:
            mymol = pybel.readstring("smi", smiles_or_filename)
        # Now we treat it as a filename
        except (OSError, IOError):
            # For now, we only support reading in a single smiles molecule,
            # but pybel returns a generator, so we get the first molecule
            # and warn the user if there is more

            mymol_generator = pybel.readfile("smi", smiles_or_filename)
            mymol_list = list(mymol_generator)
            if len(mymol_list) == 1:
                mymol = mymol_list[0]
            else:
                mymol = mymol_list[0]
                warn(
                    "More than one SMILES string in file, more than one SMILES "
                    f"string is not supported, using {mymol.write('smi')}"
                )
        mymol.make3D()
        return from_pybel(
            pybel_mol=mymol,
            compound=compound,
            infer_hierarchy=infer_hierarchy,
            ignore_box_warn=ignore_box_warn,
        )
    else:
        raise ValueError(
            "Expected SMILES loading backend 'rdkit' or 'pybel'. "
            f"Was provided: {backend}"
        )


def load_file(
    filename,
    relative_to_module=None,
    compound=None,
    coords_only=False,
    rigid=False,
    backend=None,
    infer_hierarchy=True,
    **kwargs,
):
    """Load from file into an mBuild Compound.

    Loading and converting a topology to mb.Compound from file. User can specify
    a prefered backend, or else it will be handled by a default backend based on
    the file extension.

    Parameters
    ----------
    filename : str
        Name of the file from which to load atom and bond information from
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory, look for the file
        where this module is defined. This is typically used in Compound classes
        that will be instantiated from a different directory (such as the
        Compounds located in mbuild.lib)
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information from. New structure
        will be added to the existing compound as a sub compound.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    backend : str, optional, default=None
        Backend used to load structure from file. If not specified, a default
        backend (extension specific) will be used.
    infer_hierarchy : bool, optional, default=True
        If True, infer hierarchy from chains and residues
    **kwargs : keyword arguments
        Keyword arguments passed to mdTraj for loading

    Returns
    -------
    compound : mb.Compound
    """
    # Initialize a compound if none is given
    if not compound:
        compound = mb.Compound()

    # Need to come up with a different dict structure
    default_backends = {
        ".json": "internal",
        ".xyz": "internal",
        ".sdf": "pybel",
        ".hoomdxml": "mdtraj",
        ".mol2": "mdtraj",
        ".pdb": "mdtraj",
    }

    # Handle mbuild *.py files containing a class that wraps a structure file
    # in its own folder. E.g., you build a system from ~/foo.py and it imports
    # from ~/bar/baz.py where baz.py loads ~/bar/baz.pdb.
    if relative_to_module:
        filename = str(
            Path(sys.modules[relative_to_module].__file__).parent / filename
        )
    extension = Path(filename).suffix

    if not backend:
        try:
            # Try matching backend based on extension
            backend = default_backends[extension]
        except KeyError:
            # Else use default backend
            backend = "mdtraj" if has_mdtraj else "parmed"

    # First check internal readers
    if backend == "internal":
        # Handle json format
        if extension == ".json":
            # This doesn't seem to handle the case when compound is given
            compound = compound_from_json(filename)
            return compound
        # Handle xyz file
        if extension == ".xyz" and not "top" in kwargs:
            if coords_only:
                tmp = read_xyz(filename)
                if tmp.n_particles != compound.n_particles:
                    raise ValueError(
                        "Number of atoms in {filename}"
                        "does not match {compound}".format(**locals())
                    )
                ref_and_compound = zip(
                    tmp._particles(include_ports=False),
                    compound.particles(include_ports=False),
                )
                for ref_particle, particle in ref_and_compound:
                    particle.pos = ref_particle.pos
            else:
                compound = read_xyz(filename, compound=compound)
        elif extension == ".xyz" and "top" in kwargs:
            backend = "mdtraj"

    # Then gmso reader
    if backend == "gmso":
        gmso = import_("gmso")

        top = gmso.Topology.load(filename=filename, **kwargs)
        compound = from_gmso(
            topology=top,
            compound=compound,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy,
        )

    # Then pybel reader
    elif backend == "pybel":
        pybel = import_("pybel")
        if extension == ".sdf":
            pybel_mol = pybel.readfile("sdf", filename)
            # pybel returns a generator, so we grab the first molecule of a
            # list of len 1.
            # Raise ValueError if there are more molecules
            pybel_mol = [i for i in pybel_mol]
            if len(pybel_mol) == 1:
                compound = from_pybel(
                    pybel_mol=pybel_mol[0],
                    compound=compound,
                    coords_only=coords_only,
                    infer_hierarchy=infer_hierarchy,
                )
            else:
                raise ValueError(
                    "More than one pybel molecule in file, more than one pybel "
                    "molecule is not supported"
                )

        # text file detected, assume contain smiles string
        elif extension == ".txt":
            warn(".txt file detected, loading as a SMILES string")
            # Fail-safe measure
            compound = load_smiles(filename, compound)

    # Then mdtraj reader
    elif backend == "mdtraj":
        md = import_("mdtraj")
        traj = md.load(filename, **kwargs)
        compound = from_trajectory(
            traj=traj,
            compound=compound,
            frame=-1,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy,
        )

    # Then parmed reader
    elif backend == "parmed":
        warn(
            "Using parmed reader. Bonds may be inferred from inter-particle "
            "distances and standard residue templates. Please check that the "
            "bonds in mb.Compound are accurate"
        )
        structure = pmd.load_file(filename, structure=True, **kwargs)
        compound = from_parmed(
            structure=structure,
            compound=compound,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy,
        )
    # Note: 'Input not supported' error will be handled
    # by the corresponding backend
    if rigid:
        compound.label_rigid_bodies()

    return compound


def from_parmed(
    structure, compound=None, coords_only=False, infer_hierarchy=True, **kwargs
):
    """Backend-specific loading function - parmed.

    Parameters
    ----------
    structure : pmd.Structure
        New structure to be loaded.
    compound : mb.Compound, optional, default=None
        Host mb.Compound that we are loading to.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates given by structure.
    infer_hierarchy : bool , optional, default=True
        If True, infer compound hierarchy from chain and reisdues

    Returns
    -------
    compound : mb.Compound
    """
    # Check coords_only option
    if compound and coords_only:
        if len(structure.atoms) != compound.n_particles:
            raise ValueError(
                f"Number of atoms in {structure} does not match {compound}"
                f"Structure: {len(structure.atoms)} atoms"
                f"Compound: {compound.n_particles} atoms"
            )
        if None in compound._particles(include_ports=False):
            raise ValueError("Some particles are None")

        for pmd_atom, particle in zip(
            structure.atoms, compound.particles(include_ports=False)
        ):
            particle.pos = (
                np.array([pmd_atom.xx, pmd_atom.xy, pmd_atom.xz]) / 10
            )
        return compound
    elif not compound and coords_only:
        raise MBuildError("coords_only=True but host compound is not provided")

    # Initialize a compound if none is provided
    if not compound:
        compound = mb.Compound()

    # Convert parmed structure to mbuild compound
    atom_mapping = dict()
    chain_id = None
    chains = defaultdict(list)

    # Build up chains dict
    # Could I change this to normal dict?
    for residue in structure.residues:
        chains[residue.chain].append(residue)

    # Build up compound
    chain_list = []
    for chain, residues in chains.items():
        if len(chain) > 1:
            chain_compound = mb.Compound()
            chain_list.append(chain_compound)
        else:
            chain_compound = compound
        res_list = []
        for residue in residues:
            if infer_hierarchy:
                residue_compound = mb.Compound(name=residue.name)
                parent_compound = residue_compound
                res_list.append(residue_compound)
            else:
                parent_compound = chain_compound
            atom_list = []
            atom_label_list = []
            for atom in residue.atoms:
                # Angstrom to nm
                pos = np.array([atom.xx, atom.xy, atom.xz]) / 10
                try:
                    element = element_from_atomic_number(atom.atomic_number)
                except ElementError:
                    element = None
                new_atom = mb.Particle(
                    name=str(atom.name), pos=pos, element=element
                )
                atom_list.append(new_atom)
                atom_label_list.append("{0}[$]".format(atom.name))
                atom_mapping[atom] = new_atom
            parent_compound.add(atom_list, label=atom_label_list)
        if infer_hierarchy:
            chain_compound.add(res_list)
    if len(chain) > 1:
        compound.add(chain_list)

    # Infer bonds information
    for bond in structure.bonds:
        atom1 = atom_mapping[bond.atom1]
        atom2 = atom_mapping[bond.atom2]
        compound.add_bond((atom1, atom2))

    # Convert box information
    if structure.box is not None:
        warn("All angles are assumed to be 90 degrees")
        compound.box = Box(structure.box[0:3] / 10)

    return compound


def from_trajectory(
    traj,
    compound=None,
    frame=-1,
    coords_only=False,
    infer_hierarchy=True,
    **kwargs,
):
    """Extract atoms and bonds from a md.Trajectory.

    Will create sub-compounds for every chain if there is more
    than one and sub-sub-compounds for every residue.

    Parameters
    ----------
    traj : mdtraj.Trajectory
        The trajectory to load.
    compound : mb.Compound, optional, default=None
        Host mb.Compound that we are loading to.
    frame : int, optional, default=-1 (last)
        The frame to take coordinates from.
    coords_only : bool, optional, default=False
        Only read coordinate information
    infer_hierarchy : bool, optional, default=True
        If True, infer compound hierarchy from chains and residues

    Returns
    -------
    compound : mb.Compound
    """
    # Check for coords_only option
    md = import_("mdtraj")
    if compound and coords_only:
        if traj.n_atoms != compound.n_particles:
            raise ValueError(
                "Number of atoms in {traj} does not match {compound}".format(
                    **locals()
                )
            )
        if None in compound._particles(include_ports=False):
            raise ValueError("Some particles are None")

        for mdtraj_atom, particle in zip(
            traj.topology.atoms, compound.particles(include_ports=False)
        ):
            particle.pos = traj.xyz[frame, mdtraj_atom.index]
        return compound

    elif coords_only and not compound:
        raise MBuildError("coords_only=True but host compound is not provided")

    # Initialize a compound if none is provided
    if not compound:
        compound = mb.Compound()

    atom_mapping = dict()
    # temporary lists to speed up add to the compound
    chains_list = []
    chains_list_label = []

    for chain in traj.topology.chains:
        if traj.topology.n_chains > 1:
            chain_compound = mb.Compound()
            chains_list.append(chain_compound)
            chains_list_label.append("chain[$]")
        else:
            chain_compound = compound

        res_list = []
        for res in chain.residues:
            if infer_hierarchy:
                res_compound = mb.Compound(name=res.name)
                parent_cmpd = res_compound
                res_list.append(res_compound)
            else:
                parent_cmpd = chain_compound

            atom_list = []
            atom_label_list = []
            for atom in res.atoms:
                try:
                    element = element_from_atomic_number(
                        atom.element.atomic_number
                    )
                except ElementError:
                    element = None
                new_atom = mb.Particle(
                    name=str(atom.name),
                    pos=traj.xyz[frame, atom.index],
                    element=element,
                )
                atom_list.append(new_atom)
                atom_label_list.append("{0}[$]".format(atom.name))
                atom_mapping[atom] = new_atom

            parent_cmpd.add(atom_list, label=atom_label_list)

        if infer_hierarchy:
            chain_compound.add(res_list)
    if traj.topology.n_chains > 1:
        compound.add(chains_list, label=chains_list_label)

    for mdtraj_atom1, mdtraj_atom2 in traj.topology.bonds:
        atom1 = atom_mapping[mdtraj_atom1]
        atom2 = atom_mapping[mdtraj_atom2]
        compound.add_bond((atom1, atom2))

    if np.any(traj.unitcell_lengths) and np.any(traj.unitcell_lengths[0]):
        compound.box = Box(traj.unitcell_lengths[0])

    return compound


def from_pybel(
    pybel_mol,
    compound=None,
    use_element=True,
    coords_only=False,
    infer_hierarchy=True,
    ignore_box_warn=False,
    **kwargs,
):
    """Create a Compound from a Pybel.Molecule.

    Parameters
    ----------
    pybel_mol : pybel.Molecule
        pybel Molecule that need to be converted.
    compound : mb.Compound, optional, default=None
        The host mbuild Compound.
    use_element : bool, optional, default=True
        If True, construct mb.Particle names based on the pybel Atom's element.
        If False, constructs mb.Particle names based on the pybel Atom's type.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates given by structure.
        Note: Not yet implemented, included only for parity with other
        conversion functions
    infer_hierarchy : bool, optional, default=False
        If True, infer hierarchy from residues
    ignore_box_warn: bool, optional, default=False
        If True, raise a warning if no unitcell detected for pybel molecule

    Return
    ------
    compound : mb.Compound
    """
    # Initialize a compound if none is provided
    if not compound:
        compound = mb.Compound()

    openbabel = import_("openbabel")
    compound.name = pybel_mol.title.split(".")[0]
    resindex_to_cmpd = {}

    if coords_only:
        raise Warning(
            "coords_only=True is not yet implemented for conversion from pybel"
        )

    # Iterating through pybel_mol for atom/residue information
    # This could just as easily be implemented by
    # an OBMolAtomIter from the openbabel library,
    # but this seemed more convenient at time of writing
    # pybel atoms are 1-indexed, coordinates in Angstrom
    for atom in pybel_mol.atoms:
        xyz = np.array(atom.coords) / 10
        try:
            element = element_from_atomic_number(atom.atomicnum)
        except ElementError:
            element = None
        if use_element:
            if element is None:
                warn(
                    "No element detected for atom at index "
                    f"{atom.idx} with number {atom.atomicnum}, type {atom.type}"
                )
                temp_name = atom.type
            else:
                temp_name = element.symbol
        else:
            temp_name = atom.type
        temp = mb.Particle(name=temp_name, pos=xyz, element=element)
        if infer_hierarchy and hasattr(atom, "residue"):
            # Is there a safer way to check for res?
            if atom.residue.idx not in resindex_to_cmpd:
                res_cmpd = mb.Compound(name=atom.residue.name)
                resindex_to_cmpd[atom.residue.idx] = res_cmpd
                compound.add(res_cmpd)
            resindex_to_cmpd[atom.residue.idx].add(temp)
        else:
            compound.add(temp)

    # Iterating through pybel_mol.OBMol for bond information
    # Bonds are 0-indexed, but the atoms are 1-indexed
    # Bond information doesn't appear stored in pybel_mol,
    # so we need to look into the OBMol object,
    # using an iterator from the openbabel library
    for bond in openbabel.OBMolBondIter(pybel_mol.OBMol):
        compound.add_bond(
            [
                compound[bond.GetBeginAtomIdx() - 1],
                compound[bond.GetEndAtomIdx() - 1],
            ]
        )

    if hasattr(pybel_mol, "unitcell"):
        box = Box(
            lengths=[
                pybel_mol.unitcell.GetA() / 10,
                pybel_mol.unitcell.GetB() / 10,
                pybel_mol.unitcell.GetC() / 10,
            ],
            angles=[
                pybel_mol.unitcell.GetAlpha(),
                pybel_mol.unitcell.GetBeta(),
                pybel_mol.unitcell.GetGamma(),
            ],
        )
        compound.box = box
    else:
        if not ignore_box_warn:
            warn(f"No unitcell detected for pybel.Molecule {pybel_mol}")
    return compound


def from_rdkit(rdkit_mol, compound=None, coords_only=False, smiles_seed=0):
    """Return an mbuild Compound based on a smiles string using RDKit.

    Parameters
    ----------
    rdkit_mol : rdkit.Chem.rdchem.Mol
        RDKit mol to generate an mBuild compound
    compound : mb.Compound, optional, default=None
        The host mbuild Compound.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates given by structure.
        Note: Not yet implemented, included only for parity with other
        conversion functions.
    smiles_seed : int, optional, default=0
        Random number seed for PRNG, set to -1 for non-deterministic behavior

    Returns
    -------
    mbuild.Compound

    Notes
    -----
    Option `coords_only` currently is not implemented, it is only provided to
        maintain parity with other conversion methods.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mymol = Chem.AddHs(rdkit_mol)
    if AllChem.EmbedMolecule(mymol, randomSeed=smiles_seed) != 0:
        raise MBuildError(
            f"RDKit was unable to generate 3D coordinates for {mymol}. Refer "
            "to the RDKit error messages for possible fixes. You can also "
            "install openbabel and use the backend='pybel' instead"
        )
    AllChem.UFFOptimizeMolecule(mymol)
    single_mol = mymol.GetConformer(0)
    # convert from Angstroms to nanometers
    xyz = single_mol.GetPositions() / 10

    if compound is None:
        comp = mb.Compound()
    else:
        comp = compound

    part_list = []
    for i, atom in enumerate(mymol.GetAtoms()):
        part = mb.Particle(
            name=atom.GetSymbol(),
            element=element_from_atomic_number(atom.GetAtomicNum()),
            pos=xyz[i],
        )
        part_list.append(part)

    comp.add(part_list)
    bond_order_dict = {
        Chem.BondType.SINGLE: "single",
        Chem.BondType.DOUBLE: "double",
        Chem.BondType.TRIPLE: "triple",
        Chem.BondType.AROMATIC: "aromatic",
        Chem.BondType.UNSPECIFIED: "unspecified",
    }

    for bond in mymol.GetBonds():
        comp.add_bond(
            [comp[bond.GetBeginAtomIdx()], comp[bond.GetEndAtomIdx()]],
            bond_order=bond_order_dict[bond.GetBondType()],
        )

    return comp


def from_gmso(
    topology, compound=None, coords_only=False, infer_hierarchy=True, **kwargs
):
    """Convert a GMSO Topology to mBuild Compound.

    Parameter
    ---------
    topology : gmso.Topology
        The GMSO Topology to be converted.
    compound : mb.Compound, optional, default=None
        Host mb.Compound that we are loading to.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates given by Topology.
    infer_hierarchy : bool, optional, default=True
        If True, infer compound hierarchy from Topology residue, to be implemented.
        Pending new GMSO release.

    Returns
    -------
    compound : mb.Compound
    """
    import unyt as u
    from gmso.external.convert_mbuild import to_mbuild

    if compound and coords_only:
        if topology.n_sites != compound.n_particles:
            raise ValueError(
                f"Number of sites in {topology} does not match {compound}"
                f"Topology: {topology.n_sites} sites"
                f"Compound: {compound.n_particles} particles"
            )

        if None in compound._particles(include_ports=False):
            raise ValueError("Some particles are None")

        for site, particle in zip(
            topology.sites, compound.particles(include_ports=False)
        ):
            particle.pos = np.array(site.position.to(u.nm).value)
        return compound
    elif not compound and coords_only:
        raise MBuildError("coords_only=True but host compound is not provided")

    # Convert gmso Topology to mbuild Compound
    if not compound:
        return to_mbuild(
            topology,
            infer_hierarchy=infer_hierarchy,
            **kwargs,
        )
    else:
        compound.add(
            to_mbuild(
                topology,
                infer_hierarchy=infer_hierarchy,
                **kwargs,
            )
        )
    return compound


def save(
    compound,
    filename,
    show_ports=False,
    forcefield_name=None,
    forcefield_files=None,
    forcefield_debug=False,
    box=None,
    overwrite=False,
    residues=None,
    combining_rule=None,
    foyer_kwargs=None,
    parmed_kwargs=None,
    **kwargs,
):
    """Save the Compound to a file.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that to be saved.
    filename : str
        Filesystem path in which to save the trajectory. The extension or prefix
        will be parsed and control the format. Supported extensions are:
        'hoomdxml', 'gsd', 'gro', 'top', 'lammps', 'lmp', 'mcf', 'xyz', 'pdb',
        'sdf', 'mol2', 'psf'. See parmed/structure.py for more information on
        savers.
    show_ports : bool, optional, default=False
        Save ports contained within the compound.
    forcefield_files : str, optional, default=None
        Apply a forcefield to the output file using a forcefield provided by the
        `foyer` package.
    forcefield_name : str, optional, default=None
        Apply a named forcefield to the output file using the `foyer` package,
        e.g. 'oplsaa'. Forcefields listed here:
        https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields
    forcefield_debug : bool, optional, default=False
        Choose level of verbosity when applying a forcefield through `foyer`.
        Specifically, when missing atom types in the forcefield xml file,
        determine if the warning is condensed or verbose.
    box : mb.Box, optional, default=compound.boundingbox (with buffer)
        Box information to be written to the output file. If 'None', a bounding
        box is used with 0.25nm buffers at each face to avoid overlapping atoms.
    overwrite : bool, optional, default=False
        Overwrite if the filename already exists
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by checking
        against Compound.name.
    combining_rule : str or None, optional, default=None
        Specify the combining rule for nonbonded interactions. Only relevant
        when the `foyer` package is used to apply a forcefield. Valid options
        are 'lorentz' and 'geometric' and None, specifying Lorentz-Berthelot and
        geometric combining rules respectively. If None is provided, the combining_rule
        specified from the forcefield will be used.
    foyer_kwargs : dict, optional, default=None
        Keyword arguments to provide to `foyer.Forcefield.apply`.
    parmed_kwargs : dict, optional, default=None
        Keyword arguments to provide to :meth:`mbuild.Compound.to_parmed`
    **kwargs
        Depending on the file extension these will be passed to either
        `write_gsd`, `write_hoomdxml`, `write_lammpsdata`, `write_mcf`, or
        `parmed.Structure.save`.
        See https://parmed.github.io/ParmEd/html/structobj/parmed.structure.
        Structure.html#parmed.structure.Structure.save

    Other Parameters
    ----------------
    ref_distance : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats for
        converting distance values to reduced units.
    ref_energy : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats for
        converting energy values to reduced units.
    ref_mass : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats for
        converting mass values to reduced units.
    atom_style: str, default='full'
        Defines the style of atoms to be saved in a LAMMPS data file. The
        following atom styles are currently supported: 'full', 'atomic',
        'charge', 'molecular'. See http://lammps.sandia.gov/doc/atom_style.html
        for more information on atom styles.
    unit_style: str, default='real'
        Defines to unit style to be save in a LAMMPS data file.  Defaults to
        'real' units. Current styles are supported: 'real', 'lj'. See
        https://lammps.sandia.gov/doc/99/units.html for more information on
        unit styles

    Notes
    -----
    When saving the compound as a json, only the following arguments are used:
        - filename
        - show_ports

    See Also
    --------
    formats.gsdwriter.write_gsd : Write to GSD format
    formats.hoomdxml.write_hoomdxml : Write to Hoomd XML format
    formats.xyzwriter.write_xyz : Write to XYZ format
    formats.lammpsdata.write_lammpsdata : Write to LAMMPS data format
    formats.cassandramcf.write_mcf : Write to Cassandra MCF format
    formats.json_formats.compound_to_json : Write to a json file
    """
    extension = os.path.splitext(filename)[-1]

    if extension == ".json":
        compound_to_json(compound, file_path=filename, include_ports=show_ports)
        return

    # Savers supported by mbuild.formats
    savers = {
        ".hoomdxml": write_hoomdxml,
        ".gsd": write_gsd,
        ".xyz": write_xyz,
        ".lammps": write_lammpsdata,
        ".lmp": write_lammpsdata,
        ".par": write_par,
    }
    if has_networkx:
        from mbuild.formats.cassandramcf import write_mcf

        savers.update({".mcf": write_mcf})

    try:
        saver = savers[extension]
    except KeyError:
        saver = None

    if os.path.exists(filename) and not overwrite:
        raise IOError("{0} exists; not overwriting".format(filename))

    if not parmed_kwargs:
        parmed_kwargs = {}
    structure = compound.to_parmed(
        box=box,
        residues=residues,
        show_ports=show_ports,
        **parmed_kwargs,
    )
    # Apply a force field with foyer if specified
    if forcefield_name or forcefield_files:
        foyer = import_("foyer")
        ff = foyer.Forcefield(
            forcefield_files=forcefield_files,
            name=forcefield_name,
            debug=forcefield_debug,
        )

        if not foyer_kwargs:
            foyer_kwargs = {}
        structure = ff.apply(structure, **foyer_kwargs)

        if combining_rule is None:
            combining_rule = structure.combining_rule

        if structure.combining_rule != combining_rule:
            warn(
                f"Overwriting forcefield-specified combining rule ({combining_rule})"
                f"to new combining rule ({combining_rule})."
                "This can cause inconsistent between the 1-4 pair interactions,"
                "calculated in foyer, and the new combining rule."
                "Consider directly changing the metadata of the Forcefield."
            )

        structure.combining_rule = combining_rule
        if structure.__dict__.get("defaults"):
            structure.defaults.comb_rule = (
                2 if combining_rule == "lorentz" else 3
            )

    total_charge = sum([atom.charge for atom in structure])
    if round(total_charge, 4) != 0.0:
        warn(f"System is not charge neutral. Total charge is {total_charge}.")

    # Provide a warning if rigid_ids are not sequential from 0
    if compound.contains_rigid:
        unique_rigid_ids = sorted(
            set([p.rigid_id for p in compound.rigid_particles()])
        )
        if max(unique_rigid_ids) != len(unique_rigid_ids) - 1:
            warn("Unique rigid body IDs are not sequential starting from zero.")

    if saver:  # mBuild supported saver.
        if extension in [".gsd", ".hoomdxml"]:
            kwargs["rigid_bodies"] = [p.rigid_id for p in compound.particles()]
        saver(filename=filename, structure=structure, **kwargs)

    elif extension == ".sdf":
        pybel = import_("pybel")
        new_compound = mb.Compound()
        # Convert pmd.Structure to mb.Compound
        new_compound.from_parmed(structure)
        # Convert mb.Compound to pybel molecule
        pybel_molecule = new_compound.to_pybel()
        # Write out pybel molecule to SDF file
        output_sdf = pybel.Outputfile("sdf", filename, overwrite=overwrite)
        output_sdf.write(pybel_molecule)
        output_sdf.close()

    else:  # ParmEd supported saver.
        structure.save(filename, overwrite=overwrite, **kwargs)


def catalog_bondgraph_type(compound, bond_graph=None):
    """Identify type of subgraph found at this stage of the compound.

    Parameters
    ----------
    compound : obj
        An instance of :class:`mbuild.compound.Compound`

    Returns
    -------
    str:
        "particle_graph" if compound is at the particle level,
        "one_graph" if compound is a single molecule piece,
        "multiple_graphs" if compound has multiple molecules
    """
    if not compound.children:
        return "particle_graph"
    elif bond_graph:
        # at a subgraph level
        multiple_connectionsBool = (
            len(bond_graph.subgraph(compound).connected_components()) == 1
            and len(bond_graph.subgraph(compound).connected_components()[0])
            == compound.n_particles
        )
    elif compound.bond_graph:
        # check at the top level
        multiple_connectionsBool = (
            len(compound.bond_graph.connected_components()) == 1
            and len(compound.bond_graph.connected_components()[0])
            == compound.n_particles
        )
    else:
        msg = f"`bond_graph` argument was not passed, but compound {compound} has no bond_graph attribute."
        raise ValueError(msg)

    if multiple_connectionsBool:
        return "one_graph"
    else:
        return "multiple_graphs"


def pull_residues(
    compound, segment_level=0, include_base_level=False, bond_graph=None
):
    """Pull residues from a Compound object.

    Search class instance for completed compounds based on the number of
    particles and bonds. If for example and peptide chain with three
    individual peptides that were connected and hydrated with thirty water
    molecules, the list will contain 31 residues. However, if the water and
    the individual peptides should be labeled as individual residues, set
    ``segment_level==1`` to receive a list with 33 residues. Depending on
    the method used to assemble the peptides, this procedure may continue to
    set ``segment_level=2`` and breakdown the peptide into functional groups.
    Setting ``include_base_level==True`` will allow this procedure to
    function with coarse-grained systems, while the default behavior will
    end a search at the level before atoms are obtained.

    Parameters
    ----------
    compound : obj
        An instance of :class:`mbuild.compound.Compound`
    segment_level : int, optional, default=0
        Level of full residue architecture to be identified
    include_base_level : bool, optional, default=False
        Set whether a search should continue if the list of children are single particles.
    bond_graph: obj
        An instance of :class:`mbuild.BondGraph`.
        The top level bondgraph that contains the compound to get residues from.

    Returns
    -------
    residuesList : list
        List of residue ids (str).
    """
    residuesList = []

    if (
        not bond_graph
    ):  # generate the bond graph is a top level bondgraph was not passed,
        # useful for recursion
        bond_graph = compound.bond_graph
    compound_graphtype = catalog_bondgraph_type(compound, bond_graph=bond_graph)

    # Checks segment_level and graphtype for adding particles to the residuesList
    if (
        segment_level == 0 and compound_graphtype == "multiple_graphs"
    ):  # All want to write out the children of this state
        residuesList.extend(list(map(id, compound.children)))
    elif segment_level == 0 and compound_graphtype == "one_graph":
        # At top level and a single molecule is here
        residuesList.append(id(compound))
    elif (
        compound_graphtype == "particle_graph"
    ):  # Currently at the particle level
        if include_base_level:
            # only consider adding particles if specified
            residuesList.append(id(compound))
        else:
            residuesList.append(id(compound.parent))

    # Checks for recursion
    if segment_level > 0 and compound_graphtype == "one_graph":
        # start reducing segment_level once you hit single molecules
        segment_level -= 1
        for i, child in enumerate(compound.children):
            residuesList.extend(
                pull_residues(
                    child,
                    segment_level=segment_level,
                    include_base_level=include_base_level,
                    bond_graph=bond_graph,
                )
            )
    elif segment_level > 0 and compound_graphtype == "multiple_graphs":
        # Check the next tier until you hit molecules
        for i, child in enumerate(compound.children):
            residuesList.extend(
                pull_residues(
                    child,
                    segment_level=segment_level,
                    include_base_level=include_base_level,
                    bond_graph=bond_graph,
                )
            )
    elif segment_level < 0:
        raise ValueError("`segment_level` must be greater than zero.")

    return residuesList


def to_parmed(
    compound,
    box=None,
    title="",
    residues=None,
    show_ports=False,
    infer_residues=False,
    infer_residues_kwargs={},
):
    """Create a Parmed Structure from a Compound.

    Parameters
    ----------
    compound : mb.Compound
        mbuild Compound that need to be converted.
    box : mb.Box, optional, default=None
        Box information to be used when converting to a `Structure`. If 'None'
        and the box attribute is set, the box is used with 0.25nm buffers at
        each face to avoid overlapping atoms. Otherwise the boundingbox is used
        with the same 0.25nm buffers.
    title : str, optional, default=compound.name
        Title/name of the ParmEd Structure
    residues : str of list of str, optional, default=None
        Labels of residues in the Compound. Residues are assigned by checking
        against Compound.name.
    show_ports : boolean, optional, default=False
        Include all port atoms when converting to a `Structure`.
    infer_residues : bool, optional, default=False
        Attempt to assign residues based on the number of bonds and particles in
        an object. This option is not used if `residues == None`
    infer_residues_kwargs : dict, optional, default={}
        Keyword arguments for :func:`mbuild.conversion.pull_residues`

    Returns
    -------
    parmed.structure.Structure
        ParmEd Structure object converted from compound

    See Also
    --------
    parmed.structure.Structure : Details on the ParmEd Structure object
    """
    structure = pmd.Structure()
    structure.title = title if title else compound.name
    atom_mapping = {}  # For creating bonds below
    guessed_elements = set()

    # Attempt to grab residue names based on names of children for the first
    # level of hierarchy without a box definition
    flag_res_str = True
    if not residues and infer_residues:
        residues = pull_residues(compound, **infer_residues_kwargs)
        flag_res_str = False

    if isinstance(residues, str):
        residues = [residues]
    if isinstance(residues, (list, set)):
        residues = tuple(residues)

    default_residue = pmd.Residue("RES")
    port_residue = pmd.Residue("PRT")
    compound_residue_map = dict()
    atom_residue_map = dict()

    # Loop through particles and add initialize ParmEd atoms
    for atom in compound.particles(include_ports=show_ports):
        if atom.port_particle:
            current_residue = port_residue
            atom_residue_map[atom] = current_residue

            if current_residue not in structure.residues:
                structure.residues.append(current_residue)

            pmd_atom = pmd.Atom(atomic_number=0, name="VS", mass=0, charge=0)
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10  # Angstroms

        else:
            tmp_check = atom.name if flag_res_str else id(atom)
            if residues and tmp_check in residues:
                current_residue = pmd.Residue(atom.name)
                atom_residue_map[atom] = current_residue
                compound_residue_map[atom] = current_residue
            elif residues:
                for parent in atom.ancestors():
                    tmp_check = parent.name if flag_res_str else id(parent)
                    if residues and tmp_check in residues:
                        if parent not in compound_residue_map:
                            current_residue = pmd.Residue(
                                parent.name
                                if parent.name
                                else default_residue.name
                            )
                            compound_residue_map[parent] = current_residue
                        atom_residue_map[atom] = current_residue
                        break
                else:  # Did not find specified residues in ancestors.
                    current_residue = default_residue
                    atom_residue_map[atom] = current_residue
            else:
                current_residue = default_residue
                atom_residue_map[atom] = current_residue

            if current_residue not in structure.residues:
                structure.residues.append(current_residue)

            # If we have an element attribute assigned this is easy
            if atom.element is not None:
                atomic_number = atom.element.atomic_number
                mass = atom.element.mass
            # Else we try to infer from the name
            else:
                element = _infer_element_from_compound(atom, guessed_elements)
                if element is not None:
                    atomic_number = element.atomic_number
                    mass = element.mass
                else:
                    atomic_number = 0
                    mass = 0.0

            pmd_atom = pmd.Atom(
                atomic_number=atomic_number,
                name=atom.name,
                mass=mass,
                charge=atom.charge,
            )
            # nm to Angstroms
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10.0

        residue = atom_residue_map[atom]
        structure.add_atom(pmd_atom, resname=residue.name, resnum=residue.idx)

        atom_mapping[atom] = pmd_atom

    # "Claim" all of the items it contains and subsequently index all its items
    structure.residues.claim()

    # Create and add bonds to ParmEd Structure
    for atom1, atom2 in compound.bonds():
        bond = pmd.Bond(atom_mapping[atom1], atom_mapping[atom2])
        structure.bonds.append(bond)

    # If a box is not explicitly provided:
    # (1) Grab from compound.box
    # (2) Grab from compound.get_boundingbox()
    if box is None:
        if compound.box is not None:
            box = deepcopy(compound.box)
        else:
            box = compound.get_boundingbox()
            # Pad by an extra 0.5 nm (0.25 on each side) from bounding box
            box = Box(lengths=np.array(box.lengths) + 0.5, angles=box.angles)

    box_vector = np.empty(6)
    box_vector[3:6] = box.angles
    for dim in range(3):
        box_vector[dim] = box.lengths[dim] * 10
    structure.box = box_vector

    return structure


def to_trajectory(
    compound, show_ports=False, chains=None, residues=None, box=None
):
    """Convert to an md.Trajectory and flatten the compound.

    Parameters
    ----------
    show_ports : bool, optional, default=False
        Include all port atoms when converting to trajectory.
    chains : mb.Compound or list of mb.Compound
        Chain types to add to the topology
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by checking
        against Compound.name.
    box : mb.Box, optional, default=compound.boundingbox (with buffer)
        Box information to be used when converting to a `Trajectory`.
        If 'None', a bounding box is used with a 0.5nm buffer in each
        dimension to avoid overlapping atoms.

    Returns
    -------
    trajectory : md.Trajectory

    See Also
    --------
    _to_topology
    """
    md = import_("mdtraj")
    atom_list = [particle for particle in compound.particles(show_ports)]

    top = _to_topology(compound, atom_list, chains, residues)

    # Coordinates.
    xyz = np.ndarray(shape=(1, top.n_atoms, 3), dtype="float")
    for idx, atom in enumerate(atom_list):
        xyz[0, idx] = atom.pos

    if box is None:
        box = compound.box

    # Unitcell information.
    if box is None:
        unitcell_lengths = np.array(compound.get_boundingbox().lengths) + 0.5
        unitcell_angles = [90.0, 90.0, 90.0]
    else:
        unitcell_lengths = box.lengths
        unitcell_angles = box.angles

    return md.Trajectory(
        xyz,
        top,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles,
    )


def _to_topology(compound, atom_list, chains=None, residues=None):
    """Create a mdtraj.Topology from a Compound.

    Helper function for to_trajectory.

    Parameters
    ----------
    atom_list : list of mb.Compound
        Atoms to include in the topology
    chains : mb.Compound or list of mb.Compound
        Chain types to add to the topology
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by checking
        against Compound.name.

    Returns
    -------
    top : mdtraj.Topology

    See Also
    --------
    mdtraj.Topology : Details on the mdtraj Topology object
    """
    md = import_("mdtraj")
    from mdtraj.core.element import get_by_symbol
    from mdtraj.core.topology import Topology

    if isinstance(chains, str):
        chains = [chains]
    if isinstance(chains, (list, set)):
        chains = tuple(chains)

    if isinstance(residues, str):
        residues = [residues]
    if isinstance(residues, (list, set)):
        residues = tuple(residues)
    top = Topology()
    atom_mapping = {}

    default_chain = top.add_chain()
    default_residue = top.add_residue("RES", default_chain)

    compound_residue_map = dict()
    atom_residue_map = dict()
    compound_chain_map = dict()
    atom_chain_map = dict()

    for atom in atom_list:
        # Chains
        if chains:
            if atom.name in chains:
                current_chain = top.add_chain()
                compound_chain_map[atom] = current_chain
            else:
                for parent in atom.ancestors():
                    if chains and parent.name in chains:
                        if parent not in compound_chain_map:
                            current_chain = top.add_chain()
                            compound_chain_map[parent] = current_chain
                            current_residue = top.add_residue(
                                "RES", current_chain
                            )
                        break
                else:
                    current_chain = default_chain
        else:
            current_chain = default_chain
        atom_chain_map[atom] = current_chain

        # Residues
        if residues:
            if atom.name in residues:
                current_residue = top.add_residue(atom.name, current_chain)
                compound_residue_map[atom] = current_residue
            else:
                for parent in atom.ancestors():
                    if residues and parent.name in residues:
                        if parent not in compound_residue_map:
                            current_residue = top.add_residue(
                                parent.name, current_chain
                            )
                            compound_residue_map[parent] = current_residue
                        break
                else:
                    current_residue = default_residue
        else:
            if chains:
                try:  # Grab the default residue from the custom chain.
                    current_residue = next(current_chain.residues)
                except StopIteration:  # Add the residue to the current chain
                    current_residue = top.add_residue("RES", current_chain)
            else:  # Grab the default chain's default residue
                current_residue = default_residue
        atom_residue_map[atom] = current_residue

        # Add the actual atoms
        if atom.element is not None:
            try:
                elem = get_by_symbol(atom.element.symbol)
            except KeyError:
                elem = get_by_symbol("VS")
        else:
            try:
                elem = get_by_symbol(atom.name)
            except KeyError:
                elem = get_by_symbol("VS")

        at = top.add_atom(atom.name, elem, atom_residue_map[atom])
        at.charge = atom.charge
        atom_mapping[atom] = at

    # Remove empty default residues.
    chains_to_remove = [chain for chain in top.chains if chain.n_atoms == 0]
    residues_to_remove = [res for res in top.residues if res.n_atoms == 0]
    for chain in chains_to_remove:
        top._chains.remove(chain)
    for res in residues_to_remove:
        for chain in top.chains:
            try:
                chain._residues.remove(res)
            except ValueError:  # Already gone.
                pass

    for atom1, atom2 in compound.bonds():
        # Ensure that both atoms are part of the compound. This becomes an
        # issue if you try to convert a sub-compound to a topology which is
        # bonded to a different subcompound.
        if all(a in atom_mapping.keys() for a in [atom1, atom2]):
            top.add_bond(atom_mapping[atom1], atom_mapping[atom2])
    return top


def to_pybel(
    compound,
    box=None,
    title="",
    residues=None,
    show_ports=False,
    infer_residues=False,
):
    """Create a pybel.Molecule from a Compound.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that need to be converted.
    box : mb.Box, optional, default=None
    title : str, optional, default=compound.name
        Title/name of the ParmEd Structure
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by checking
        against Compound.name.
    show_ports : boolean, optional, default=False
        Include all port atoms when converting to a `Structure`.
    infer_residues : bool, optional, default=False
        Attempt to assign residues based on names of children

    Returns
    -------
    pybelmol : pybel.Molecule

    Notes
    -----
    Most of the mb.Compound is first converted to openbabel.OBMol and then pybel
    creates a pybel.Molecule from the OBMol. Bond orders are assumed to be 1
    OBMol atom indexing starts at 1, with spatial dimension Angstrom
    """
    openbabel = import_("openbabel")
    pybel = import_("pybel")

    mol = openbabel.OBMol()
    particle_to_atom_index = {}
    guessed_elements = set()

    if not residues and infer_residues:
        residues = list(set([child.name for child in compound.children]))
    if isinstance(residues, str):
        residues = [residues]
    if isinstance(residues, (list, set)):
        residues = tuple(residues)

    compound_residue_map = dict()
    atom_residue_map = dict()

    for i, part in enumerate(compound.particles(include_ports=show_ports)):
        if residues and part.name in residues:
            current_residue = mol.NewResidue()
            current_residue.SetName(part.name)
            atom_residue_map[part] = current_residue
            compound_residue_map[part] = current_residue
        elif residues:
            for parent in part.ancestors():
                if residues and parent.name in residues:
                    if parent not in compound_residue_map:
                        current_residue = mol.NewResidue()
                        current_residue.SetName(parent.name)
                        compound_residue_map[parent] = current_residue
                    atom_residue_map[part] = current_residue
                    break
            else:  # Did not find specified residues in ancestors.
                current_residue = mol.NewResidue()
                current_residue.SetName("RES")
                atom_residue_map[part] = current_residue
        else:
            current_residue = mol.NewResidue()
            current_residue.SetName("RES")
            atom_residue_map[part] = current_residue

        temp = mol.NewAtom()
        residue = atom_residue_map[part]
        temp.SetResidue(residue)
        if part.port_particle:
            temp.SetAtomicNum(0)
        else:
            if part.element is not None:
                temp.SetAtomicNum(part.element.atomic_number)
            else:
                element = _infer_element_from_compound(part, guessed_elements)
                if element is not None:
                    temp.SetAtomicNum(element.atomic_number)
                else:
                    temp.SetAtomicNum(0)

        temp.SetVector(*(part.xyz[0] * 10))
        particle_to_atom_index[part] = i

    ucell = openbabel.OBUnitCell()
    if box is None:
        box = compound.get_boundingbox()
    (a, b, c) = box.lengths
    a *= 10
    b *= 10
    c *= 10
    alpha, beta, gamma = np.radians(box.angles)
    ucell.SetData(a, b, c, alpha, beta, gamma)
    mol.CloneData(ucell)

    for bond in compound.bonds():
        bond_order = 1
        mol.AddBond(
            particle_to_atom_index[bond[0]] + 1,
            particle_to_atom_index[bond[1]] + 1,
            bond_order,
        )

    pybelmol = pybel.Molecule(mol)
    pybelmol.title = title if title else compound.name

    return pybelmol


def to_rdkit(compound):
    """Create an RDKit RWMol from an mBuild Compound.

    Parameters
    ----------
    compound : mbuild.Compound; required
        The compound to convert to a Chem.RWmol

    Returns
    -------
    rdkit.Chem.RWmol
    """
    rdkit = import_("rdkit")
    from rdkit import Chem
    from rdkit.Chem import AllChem

    for particle in compound.particles():
        if particle.element is None:
            try:
                particle._element = element_from_symbol(particle.name)
            except ElementError:
                try:
                    particle._element = element_from_name(particle.name)
                except ElementError:
                    raise MBuildError(
                        f"No element assigned to {particle};"
                        "element could not be"
                        f"inferred from particle name {particle.name}."
                        " Cannot perform an energy minimization."
                    )

    temp_mol = Chem.RWMol()
    p_dict = {particle: i for i, particle in enumerate(compound.particles())}

    bond_order_dict = {
        "single": Chem.BondType.SINGLE,
        "double": Chem.BondType.DOUBLE,
        "triple": Chem.BondType.TRIPLE,
        "aromatic": Chem.BondType.AROMATIC,
        "unspecified": Chem.BondType.UNSPECIFIED,
        "default": Chem.BondType.SINGLE,
    }

    for particle in compound.particles():
        temp_atom = Chem.Atom(particle.element.atomic_number)

        # this next line is necessary to prevent rdkit from adding hydrogens
        # this will also set the label to be the element with particle index
        temp_atom.SetProp(
            "atomLabel", f"{temp_atom.GetSymbol()}:{p_dict[particle]}"
        )

        temp_mol.AddAtom(temp_atom)

    for bond in compound.bonds(return_bond_order=True):
        bond_indices = (p_dict[bond[0]], p_dict[bond[1]])
        temp_mol.AddBond(*bond_indices)
        rdkit_bond = temp_mol.GetBondBetweenAtoms(*bond_indices)
        rdkit_bond.SetBondType(bond_order_dict[bond[2]["bond_order"]])

    return temp_mol


def to_smiles(compound, backend="pybel"):
    """Create a SMILES string from an mbuild compound.

    Parameters
    ----------
    compound : mb.Compound.
        The mbuild compound to be converted.
    backend : str, optional, default="pybel"
        Backend used to do the conversion.

    Return
    ------
    smiles_string : str
    """
    if backend == "pybel":
        mol = to_pybel(compound)

        warn(
            "The bond orders will be guessed using pybel"
            "OBMol.PerceviedBondOrders()"
        )
        mol.OBMol.PerceiveBondOrders()
        smiles_string = mol.write("smi").replace("\t", " ").split(" ")[0]

        return smiles_string
    else:
        raise NotImplementedError(f"Backend {backend} not implemented.")


def to_networkx(compound, names_only=False):
    """Create a NetworkX graph representing the hierarchy of a Compound.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that need to be converted.
    names_only : bool, optional, default=False
        Store only the names of the compounds in the graph, appended with their
        IDs, for distinction even if they have the same name. When set to False,
        the default behavior, the nodes are the compounds themselves.

    Returns
    -------
    G : networkx.DiGraph

    Notes
    -----
    This digraph is not the bondgraph of the compound.

    See Also
    --------
    mbuild.bond_graph
    """
    nx = import_("networkx")

    nodes = list()
    edges = list()
    if names_only:
        nodes.append(compound.name + "_" + str(id(compound)))
    else:
        nodes.append(compound)
    nodes, edges = _iterate_children(
        compound, nodes, edges, names_only=names_only
    )

    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    return graph


def _iterate_children(compound, nodes, edges, names_only=False):
    """Create nodes and edges that connect parents and their children.

    Helper function for to_networkx
    """
    if not compound.children:
        return nodes, edges
    for child in compound.children:
        if names_only:
            unique_name = child.name + "_" + str(id(child))
            unique_name_parent = (
                child.parent.name + "_" + str((id(child.parent)))
            )
            nodes.append(unique_name)
            edges.append([unique_name_parent, unique_name])
        else:
            nodes.append(child)
            edges.append([child.parent, child])
        nodes, edges = _iterate_children(
            child, nodes, edges, names_only=names_only
        )
    return nodes, edges


def to_gmso(compound, box=None, **kwargs):
    """Create a GMSO Topology from a mBuild Compound.

    Parameters
    ----------
    compound : mb.Compound
        The mb.Compound to be converted.
    box : mb.Box, optional, default=None
        The mb.Box to be converted, if different that compound.box

    Returns
    -------
    topology : gmso.Topology
        The converted gmso Topology
    """
    from gmso.external.convert_mbuild import from_mbuild

    return from_mbuild(compound, box=None, **kwargs)


def to_intermol(compound, molecule_types=None):  # pragma: no cover
    """Create an InterMol system from a Compound.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that need to be converted.
    molecule_types : list or tuple of subclasses of Compound

    Returns
    -------
    intermol_system : intermol.system.System
    """
    import simtk.unit as u
    from intermol.atom import Atom as InterMolAtom
    from intermol.molecule import Molecule
    from intermol.system import System

    if isinstance(molecule_types, list):
        molecule_types = tuple(molecule_types)
    elif molecule_types is None:
        molecule_types = (compound.name,)
    intermol_system = System()

    last_molecule_compound = None
    for atom_index, atom in enumerate(compound.particles()):
        for parent in atom.ancestors():
            # Don't want inheritance via isinstance().
            if parent.name in molecule_types:
                # Check if we have encountered this molecule type before.
                if parent.name not in intermol_system.molecule_types:
                    _add_intermol_molecule_type(intermol_system, parent)
                if parent != last_molecule_compound:
                    last_molecule_compound = parent
                    last_molecule = Molecule(name=parent.name)
                    intermol_system.add_molecule(last_molecule)
                break
        else:
            # Should never happen if molecule_types only contains
            # type(compound)
            raise ValueError(
                f"Found an atom {atom} that is not part of any of the "
                f"specified molecule types {molecule_types}"
            )

        # Add the actual intermol atoms.
        intermol_atom = InterMolAtom(
            atom_index + 1, name=atom.name, residue_index=1, residue_name="RES"
        )
        intermol_atom.position = atom.pos * u.nanometers
        last_molecule.add_atom(intermol_atom)
    return intermol_system


def _add_intermol_molecule_type(intermol_system, parent):  # pragma: no cover
    """Create a molecule type for the parent and add bonds.

    This method takes an intermol system and adds a
    parent compound, including its particles and bonds, to it.
    """
    from intermol.forces.abstract_bond_type import (
        AbstractBondType as InterMolBond,
    )
    from intermol.moleculetype import MoleculeType

    molecule_type = MoleculeType(name=parent.name)
    intermol_system.add_molecule_type(molecule_type)

    for index, parent_atom in enumerate(parent.particles()):
        parent_atom.index = index + 1

    for atom1, atom2 in parent.bonds():
        intermol_bond = InterMolBond(atom1.index, atom2.index)
        molecule_type.bond_forces.add(intermol_bond)


def _infer_element_from_compound(compound, guessed_elements):
    """Infer the element from the compound name.

    Parameters
    ----------
    compound : mbuild.Compound
        the compound to infer the element for
    guessed_elements : list
        a list of the already-guessed-elements

    Returns
    -------
    element : ele.Element or None
    """
    try:
        element = element_from_symbol(compound.name)
    except ElementError:
        try:
            element = element_from_name(compound.name)
            warn_msg = (
                f"No element attribute associated with '{compound}'; "
                f"Guessing it is element '{element}'"
            )
        except ElementError:
            element = None
            warn_msg = (
                f"No element attribute associated with '{compound}'; "
                "and no matching elements found based upon the "
                "compound name. Setting atomic number to zero."
            )
        if compound.name not in guessed_elements:
            warn(warn_msg)
            guessed_elements.add(compound.name)

    return element
