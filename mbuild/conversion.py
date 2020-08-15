import os
import sys
from warnings import warn
from pathlib import Path
from collections import defaultdict
import numpy as np

import parmed as pmd
from parmed.periodic_table import AtomicNum, element_by_name, Mass, Element

import mbuild as mb
from mbuild.box import Box
from mbuild.exceptions import MBuildError
from mbuild.formats.xyz import read_xyz, write_xyz
from mbuild.formats.json_formats import compound_to_json, compound_from_json
from mbuild.formats.hoomdxml import write_hoomdxml
from mbuild.formats.lammpsdata import write_lammpsdata
from mbuild.formats.gsdwriter import write_gsd
from mbuild.formats.par_writer import write_par
from mbuild.utils.io import import_, has_networkx, has_openbabel, has_mdtraj


def load(filename_or_object,
         relative_to_module=None,
         compound=None,
         coords_only=False,
         rigid=False,
         smiles=False,
         infer_hierarchy=True,
         backend=None,
         ignore_box_warn=False,
         **kwargs):
    """Load a file or an existing topology into an mbuild compound.

    Files are read using the MDTraj package unless the `use_parmed` argument is
    specified as True. Please refer to http://mdtraj.org/1.8.0/load_functions.html
    for formats supported by MDTraj and https://parmed.github.io/ParmEd/html/
    readwrite.html for formats supported by ParmEd.

    Parameters
    ----------
    filename_or_object : str, mdtraj.Trajectory, parmed.Structure, mbuild.Compound,
            pybel.Molecule
        Name of the file or topology from which to load atom and bond information.
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory, look for the file
        where this module is defined. This is typically used in Compound
        classes that will be instantiated from a different directory
        (such as the Compounds located in mbuild.lib).
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information into.
        New structure will be added to the existing compound as a sub compound.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    backend : str, optional, default=None
        Backend used to load structure from file. If not specified, a default
        backend (extension specific) will be used.
    smiles: bool, optional, default=False
        Use Open Babel to parse filename as a SMILES string
        or file containing a SMILES string.
    infer_hierarchy : bool, optional, default=True
        If True, infer hierarchy from chains and residues
    ignore_box_warn : bool, optional, default=False
        If True, ignore warning if no box is present.
    **kwargs : keyword arguments
        Key word arguments passed to mdTraj for loading.

    Returns
    -------
    compound : mb.Compound
    """
    # First check if we are loading from an object
    if not isinstance(filename_or_object, str):
        return load_object(
            obj=filename_or_object,
            compound=compound,
            coords_only=coords_only,
            rigid=rigid,
            infer_hierarchy=infer_hierarchy,
            **kwargs
        )
    # Second check if we are loading SMILES strings
    elif smiles:
        return load_smiles(
            smiles_or_filename=filename_or_object,
            compound=compound,
            infer_hierarchy=infer_hierarchy,
            ignore_box_warn=ignore_box_warn
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
            **kwargs
        )


def load_object(obj,
                compound=None,
                coords_only=False,
                rigid=False,
                infer_hierarchy=True,
                **kwargs):
    """Helper function to load an obj into a mb.Compound

    Functions to load on-disk obj to mb.Compound, supporting conversion
    from mb.Compound, pmd.Structure, md.Trajectory, and pybel.Molecule. If the
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
    if has_openbabel:
        pybel = import_('pybel')
        type_dict.update({pybel.Molecule: from_pybel})

    if has_mdtraj:
        md = import_('mdtraj')
        type_dict.update({md.Trajectory: from_trajectory})

    # Check if the given object is an mb.Compound
    if isinstance(obj, mb.Compound):
        if not compound:
            warn('Given object is an mb.Compound, \
                  do nothing and return the object.')
            return obj
        else:
            warn('Given object is an mb.Compound, \
                  adding object to the host compound.')
            compound.add(obj)
            return compound

    for type_ in type_dict:
        if isinstance(obj, type_):
            compound = type_dict[type_](
                obj,
                compound,
                coords_only=coords_only,
                infer_hierarchy=infer_hierarchy,
                **kwargs
            )
            if rigid:
                compound.label_rigid_bodies()
            return compound

    # If nothing is return raise an error
    raise ValueError(f'Object of type {type(obj).__name__} is not supported')


def load_smiles(smiles_or_filename,
                compound=None,
                infer_hierarchy=True,
                ignore_box_warn=False):
    """Helper function to load a SMILES string

    Loading SMILES string from a string, a list, or a file using pybel.
    Must have pybel packages installed.

    Parameters
    ----------
    smiles_or_filename : str
        SMILES string or file of SMILES string to load
    compound : mb.Compound
        The host mbuild Compound
    infer_hierarchy : bool, optional, default=True
    ignore_box_warn : bool, optional, default=False
        If True, ignore warning if no box is present.

    Returns
    -------
    compound : mb.Compound

    """
    # Will try to support list of smiles strings in the future
    pybel = import_('pybel')

    # Initialize an mb.Compound if none is provided
    if not compound:
        compound = mb.Compound()

    # First we try treating filename_or_object as a SMILES string
    try:
        mymol = pybel.readstring("smi", smiles_or_filename)
    # Now we treat it as a filename
    except(OSError, IOError):
        # For now, we only support reading in a single smiles molecule,
        # but pybel returns a generator, so we get the first molecule
        # and warn the user if there is more

        mymol_generator = pybel.readfile("smi", smiles_or_filename)
        mymol_list = list(mymol_generator)
        if len(mymol_list) == 1:
            mymol = mymol_list[0]
        else:
            mymol = mymol_list[0]
            warn("More than one SMILES string in file, more than one SMILES "
                 "string is not supported, using {}".format(mymol.write("smi")))
    mymol.make3D()

    return from_pybel(
        pybel_mol=mymol,
        compound=compound,
        infer_hierarchy=infer_hierarchy,
        ignore_box_warn=ignore_box_warn
    )


def load_file(filename,
              relative_to_module=None,
              compound=None,
              coords_only=False,
              rigid=False,
              backend=None,
              infer_hierarchy=True,
              **kwargs):
    """ Helper function to load from files

    Loading and converting a topology to mb.Compound from file. User can specify
    a prefered backend, or else it will be handled by a default backend based on
    the file extension.

    Parameters
    ----------
    filename : str
        Name of the file from which to load atom and bond information from
    relative_to_module : str, optional, default=None
        Instead of looking in the current working directory,
        look for the file where this module is defined. This
        is typically used in Compound classes that will be
        instantiated from a different directory (such as the
        Compounds located in mbuid.lib
    compound : mb.Compound, optional, default=None
        Existing compound to load atom and bond information from.
        New structure will be added to the existing compound
        as a sub compound.
    coords_only : bool, optional, default=False
        Only load the coordinates into an existing compound.
    rigid : bool, optional, default=False
        Treat the compound as a rigid body
    backend : str, optional, default=None
        Backend used to load structure from file. if not specified, a default
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
        '.json': 'internal',
        '.xyz': 'internal',
        '.sdf': 'pybel',
        '.hoomdxml': 'mdtraj',
        '.mol2': 'mdtraj',
        '.pdb': 'mdtraj'
    }

    # Handle mbuild *.py files containing a class that wraps a structure file
    # in its own folder. E.g., you build a system from ~/foo.py and it imports
    # from ~/bar/baz.py where baz.py loads ~/bar/baz.pdb.
    if relative_to_module:
        filename = str(Path(sys.modules[relative_to_module].__file__).parent / filename)
    extension = Path(filename).suffix

    if not backend:
        try:
            # Try matching backend based on extension
            backend = default_backends[extension]
        except KeyError:
            # Else use default backend
            backend = 'mdtraj' if has_mdtraj else 'parmed'

    # First check internal readers
    if backend == 'internal':
        # Handle json format
        if extension == '.json':
            # This doesn't seem to handle the case when compound is given
            compound = compound_from_json(filename)
            return compound
        # Handle xyz file
        if extension == '.xyz' and not 'top' in kwargs:
            if coords_only:
                tmp = read_xyz(filename)
                if tmp.n_particles != compound.n_particles:
                    raise ValueError(
                        'Number of atoms in {filename}'
                        'does not match {compound}'.format(**locals())
                    )
                ref_and_compound = zip(
                    tmp._particles(include_ports=False),
                    compound.particles(include_ports=False)
                )
                for ref_particle, particle in ref_and_compound:
                    particle.pos = ref_particle.pos
            else:
                compound = read_xyz(filename, compound=compound)
        elif extension == '.xyz' and 'top' in kwargs:
            backend = 'mdtraj'

    # Then pybel reader
    elif backend == 'pybel':
        pybel = import_('pybel')
        if extension == '.sdf':
            pybel_mol = pybel.readfile('sdf', filename)
            # pybel returns a generator, so we grab the first molecule of a
            # list of len 1.
            # Raise ValueError if there are more molecules
            pybel_mol = [i for i in pybel_mol]
            if len(pybel_mol) == 1:
                compound = from_pybel(
                    pybel_mol=pybel_mol[0],
                    compound=compound,
                    coords_only=coords_only,
                    infer_hierarchy=infer_hierarchy
                )
            else:
                raise ValueError('More than one pybel molecule in file, '
                                 'more than one pybel molecule is not supported')

        # text file detected, asssume contain smiles string
        elif extension == '.txt':
            warn('.txt file detected, loading as a SMILES string')
            # Fail-safe measure
            compound = load_smiles(filename, compound)

    # Then mdtraj reader
    elif backend == 'mdtraj':
        md = import_('mdtraj')
        traj = md.load(filename, **kwargs)
        compound = from_trajectory(
            traj=traj,
            compound=compound,
            frame=-1,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy
        )

    # Then parmed reader
    elif backend == 'parmed':
        warn('Using parmed reader. Bonds may be inferred '
             'from inter-particle distances and standard '
             'residue templates. Please check that the bonds '
             'in mb.Compound are accurate')
        structure = pmd.load_file(filename, structure=True, **kwargs)
        compound = from_parmed(
            structure=structure,
            compound=compound,
            coords_only=coords_only,
            infer_hierarchy=infer_hierarchy
        )
    # Note: 'Input not supported' error will be handled
    # by the corresponding backend
    if rigid:
        compound.label_rigid_bodies()

    return compound


def from_parmed(structure,
                compound=None,
                coords_only=False,
                infer_hierarchy=True):
    """Backend-specific loading function - parmed

    Parameters
    ----------
    structure : pmd.Structure
        New structure to be loaded.
    compound : mb.Compound, optional, default=None
        Host mb.Compound that we are loading to.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates
        given by structure.
    infer_hierarchy : bool , optional, default=True
        If True, infer compound hierarchy from chain
        and reisdues

    Returns
    -------
    compound : mb.Compound

    """
    # Check coords_only option
    if compound and coords_only:
        if len(structure.atoms) != compound.n_particles:
            raise ValueError(
                'Number of atoms in {structure} does not '
                '{compound}'.formats(**locals())
            )
        atoms_particles = zip(
            structure.atoms,
            compound.particles(include_ports=False)
        )
        if None in compound._particles(include_ports=False):
            raise ValueError('Some particles are None')

        for pmd_atom, particle in atoms_particles:
            particle.pos = np.array([pmd_atom.xx,
                                     pmd_atom.xy,
                                     pmd_atom.xz]) / 10
        return compound
    elif not compound and coords_only:
        raise MBuildError(
            'coords_only=True but host compound is not provided'
        )

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
    for chain, residues in chains.items():
        if len(chain) > 1:
            chain_compound = mb.Compound()
            compound.add(chain_compound, chain_id)
        else:
            chain_compound = compound
        for residue in residues:
            if infer_hierarchy:
                residue_compound = mb.Compound(name=residue.name)
                chain_compound.add(residue_compound)
                parent_compound = residue_compound
            else:
                parent_compound = chain_compound
            for atom in residue.atoms:
                pos = np.array([atom.xx,
                                atom.xy,
                                atom.xz]) / 10
                new_atom = mb.Particle(name=str(atom.name), pos=pos)
                parent_compound.add(new_atom, label='{0}[$]'.format(atom.name))
                atom_mapping[atom] = new_atom

    # Infer bonds information
    for bond in structure.bonds:
        atom1 = atom_mapping[bond.atom1]
        atom2 = atom_mapping[bond.atom2]
        compound.add_bond((atom1, atom2))

    # Convert box information
    if structure.box is not None:
        warn('All angles are assumed to be 90 degrees')
        compound.periodicity = structure.box[0:3] / 10
    else:
        warn('No box information detected, periodicity is set to [0, 0, 0]')
        compound.periodicity = np.array([0., 0., 0.])

    return compound


def from_trajectory(traj,
                    compound=None,
                    frame=-1,
                    coords_only=False,
                    infer_hierarchy=True):
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
    md = import_('mdtraj')
    if compound and coords_only:
        if traj.n_atoms != compound.n_particles:
            raise ValueError(
                'Number of atoms in {traj} does not match {compound}'.format(**locals())
            )
        atoms_particles = zip(
            traj.topology.atoms,
            compound.particles(include_ports=False)
        )

        if None in compound._particles(include_ports=False):
            raise ValueError('Some particles are None')

        for mdtraj_atom, particle in atoms_particles:
            particle.pos = traj.xyz[frame, mdtraj_atom.index]
        return compound

    elif coords_only and not compound:
        raise MBuildError(
            'coords_only=True but host compound is not provided'
        )

    # Initialize a compound if none is provided
    if not compound:
        compound = mb.Compound()

    atom_mapping = dict()
    for chain in traj.topology.chains:
        if traj.topology.n_chains > 1:
            chain_compound = mb.Compound()
            compound.add(chain_compound, 'chain[$]')
        else:
            chain_compound = compound
        for res in chain.residues:
            if infer_hierarchy:
                res_compound = mb.Compound(name=res.name)
                chain_compound.add(res_compound)
                parent_cmpd = res_compound
            else:
                parent_cmpd = chain_compound
            for atom in res.atoms:
                new_atom = mb.Particle(
                    name=str(atom.name),
                    pos=traj.xyz[frame, atom.index]
                )
                parent_cmpd.add(
                    new_atom,
                    label='{0}[$]'.format(atom.name)
                )
                atom_mapping[atom] = new_atom

    for mdtraj_atom1, mdtraj_atom2 in traj.topology.bonds:
        atom1 = atom_mapping[mdtraj_atom1]
        atom2 = atom_mapping[mdtraj_atom2]
        compound.add_bond((atom1, atom2))

    if np.any(traj.unitcell_lengths) and np.any(traj.unitcell_lengths[0]):
        compound.periodicity = traj.unitcell_lengths[0]
    else:
        compound.periodicity = np.array([0., 0., 0.])

    return compound


def from_pybel(pybel_mol,
               compound=None,
               use_element=True,
               coords_only=False,
               infer_hierarchy=True,
               ignore_box_warn=False):
    """Create a Compound from a Pybel.Molecule

    Parameters
    ---------
    pybel_mol : pybel.Molecule
        pybel Molecule that need to be converted.
    compound : mb.Compound, optional, default=None
        The host mbuild Compound.
    use_element : bool, optional, default=True
        If True, construct mb Particles based on the pybel Atom's element.
        If False, construcs mb Particles based on the pybel Atom's type.
    coords_only : bool, optional, default=False
        Set preexisting atoms in compound to coordinates given
        by structure. Note: Not yet implemented, included only
        for parity with other conversion functions
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
    compound.name = pybel_mol.title.split('.')[0]
    resindex_to_cmpd = {}

    if coords_only:
        raise Warning(
            'coords_only=True is not yet implemented '
            'for conversion from pybel'
        )

    # Iterating through pybel_mol for atom/residue information
    # This could just as easily be implemented by
    # an OBMolAtomIter from the openbabel library,
    # but this seemed more convenient at time of writing
    # pybel atoms are 1-indexed, coordinates in Angstrom
    for atom in pybel_mol.atoms:
        xyz = np.array(atom.coords)/10
        if use_element:
            try:
                temp_name = Element[atom.atomicnum]
            except KeyError:
                warn(
                    "No element detected for atom at index "
                    "{} with number {}, type {}".format(
                        atom.idx,
                        atom.atomicnum,
                        atom.type
                    )
                )
                temp_name = atom.type
        else:
            temp_name = atom.type
        temp = mb.Particle(name=temp_name, pos=xyz)
        if infer_hierarchy and hasattr(atom, 'residue'):
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
        compound.add_bond([compound[bond.GetBeginAtomIdx()-1],
                        compound[bond.GetEndAtomIdx()-1]])

    if hasattr(pybel_mol, 'unitcell'):
        box = Box(lengths=[pybel_mol.unitcell.GetA()/10,
                            pybel_mol.unitcell.GetB()/10,
                            pybel_mol.unitcell.GetC()/10],
                    angles=[pybel_mol.unitcell.GetAlpha(),
                            pybel_mol.unitcell.GetBeta(),
                            pybel_mol.unitcell.GetGamma()])
        compound.periodicity = box.lengths
    else:
        if not ignore_box_warn:
            warn("No unitcell detected for pybel.Molecule {}".format(pybel_mol))
#       TODO: Decide how to gather PBC information from openbabel. Options may
#             include storing it in .periodicity or writing a separate function
#             that returns the box.
    return compound


def save(compound,
         filename,
         show_ports=False,
         forcefield_name=None,
         forcefield_files=None,
         forcefield_debug=False,
         box=None,
         overwrite=False,
         residues=None,
         combining_rule='lorentz',
         foyer_kwargs=None,
         **kwargs):
    """Save the Compound to a file.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that to be saved.
    filename : str
        Filesystem path in which to save the trajectory. The extension or
        prefix will be parsed and control the format. Supported
        extensions are: 'hoomdxml', 'gsd', 'gro', 'top',
        'lammps', 'lmp', 'mcf'.
    show_ports : bool, optional, default=False
        Save ports contained within the compound.
    forcefield_files : str, optional, default=None
        Apply a forcefield to the output file using a forcefield provided
        by the `foyer` package.
    forcefield_name : str, optional, default=None
        Apply a named forcefield to the output file using the `foyer`
        package, e.g. 'oplsaa'. Forcefields listed here:
        https://github.com/mosdef-hub/foyer/tree/master/foyer/forcefields
    forcefield_debug : bool, optional, default=False
        Choose level of verbosity when applying a forcefield through `foyer`.
        Specifically, when missing atom types in the forcefield xml file,
        determine if the warning is condensed or verbose.
    box : mb.Box, optional, default=compound.boundingbox (with buffer)
        Box information to be written to the output file. If 'None', a
        bounding box is used with 0.25nm buffers at each face to avoid
        overlapping atoms.
    overwrite : bool, optional, default=False
        Overwrite if the filename already exists
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by
        checking against Compound.name.
    combining_rule : str, optional, default='lorentz'
        Specify the combining rule for nonbonded interactions. Only relevant
        when the `foyer` package is used to apply a forcefield. Valid
        options are 'lorentz' and 'geometric', specifying Lorentz-Berthelot
        and geometric combining rules respectively.
    foyer_kwargs : dict, optional, default=None
        Keyword arguments to provide to `foyer.Forcefield.apply`.
    **kwargs
        Depending on the file extension these will be passed to either
        `write_gsd`, `write_hoomdxml`, `write_lammpsdata`,
        `write_mcf`, or `parmed.Structure.save`.
        See https://parmed.github.io/ParmEd/html/structobj/parmed.structure.Structure.html#parmed.structure.Structure.save


    Other Parameters
    ----------------
    ref_distance : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats
        for converting distance values to reduced units.
    ref_energy : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats
        for converting energy values to reduced units.
    ref_mass : float, optional, default=1.0
        Normalization factor used when saving to .gsd and .hoomdxml formats
        for converting mass values to reduced units.
    atom_style: str, default='full'
        Defines the style of atoms to be saved in a LAMMPS data file. The following atom
        styles are currently supported: 'full', 'atomic', 'charge', 'molecular'
        see http://lammps.sandia.gov/doc/atom_style.html for more
        information on atom styles.
    unit_style: str, default='real'
        Defines to unit style to be save in a LAMMPS data file.  Defaults to 'real' units.
        Current styles are supported: 'real', 'lj'
        see https://lammps.sandia.gov/doc/99/units.html for more information
        on unit styles

    Notes
    ------
    When saving the compound as a json, only the following arguments are used:
        - filename
        - show_ports

    See Also
    --------
    formats.gsdwrite.write_gsd : Write to GSD format
    formats.hoomdxml.write_hoomdxml : Write to Hoomd XML format
    formats.xyzwriter.write_xyz : Write to XYZ format
    formats.lammpsdata.write_lammpsdata : Write to LAMMPS data format
    formats.cassandramcf.write_mcf : Write to Cassandra MCF format
    formats.json_formats.compound_to_json : Write to a json file

    """
    extension = os.path.splitext(filename)[-1]

    if extension == '.json':
        compound_to_json(compound,
                         file_path=filename,
                         include_ports=show_ports)
        return

    # Savers supported by mbuild.formats
    savers = {
        '.hoomdxml': write_hoomdxml,
        '.gsd': write_gsd,
        '.xyz': write_xyz,
        '.lammps': write_lammpsdata,
        '.lmp': write_lammpsdata,
        '.par': write_par
    }
    if has_networkx:
        from mbuild.formats.cassandramcf import write_mcf
        savers.update({'.mcf': write_mcf})

    try:
        saver = savers[extension]
    except KeyError:
        saver = None

    if os.path.exists(filename) and not overwrite:
        raise IOError('{0} exists; not overwriting'.format(filename))

    structure = compound.to_parmed(box=box, residues=residues,
                               show_ports=show_ports)
    # Apply a force field with foyer if specified
    if forcefield_name or forcefield_files:
        foyer = import_('foyer')
        ff = foyer.Forcefield(
            forcefield_files=forcefield_files,
            name=forcefield_name,
            debug=forcefield_debug
        )

        if not foyer_kwargs:
            foyer_kwargs = {}
        structure = ff.apply(structure, **foyer_kwargs)
        structure.combining_rule = combining_rule

    total_charge = sum([atom.charge for atom in structure])
    if round(total_charge, 4) != 0.0:
        warn('System is not charge neutral. Total charge is {}.'
             ''.format(total_charge))

    # Provide a warning if rigid_ids are not sequential from 0
    if compound.contains_rigid:
        unique_rigid_ids = sorted(set([
            p.rigid_id for p in compound.rigid_particles()]))
        if max(unique_rigid_ids) != len(unique_rigid_ids) - 1:
            warn("Unique rigid body IDs are not sequential starting from zero.")

    if saver:  # mBuild supported saver.
        if extension in ['.gsd', '.hoomdxml']:
            kwargs['rigid_bodies'] = [
                    p.rigid_id for p in compound.particles()]
        saver(filename=filename, structure=structure, **kwargs)

    elif extension == '.sdf':
        pybel = import_('pybel')
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


def to_parmed(compound,
              box=None,
              title='',
              residues=None,
              show_ports=False,
              infer_residues=False):
    """ Create a Parmed Structure from a Compound.

    Parameters
    ----------
    compound : mb.Compound
        mbuild Compound that need to be converted.
    box : mb.Box, optional, default=compound.boundingbox (with buffer)
        Box information to be used when converting to a `Structure`.
        If 'None', a bounding box is used with 0.25nm buffers at
        each face to avoid overlapping atoms, unless `compound.periodicity`
        is not None, in which case those values are used for the
        box lengths.
    title : str, optional, default=compound.name
        Title/name of the ParmEd Structure
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by
        checking against Compound.name.
    show_ports : boolean, optional, default=False
        Include all port atoms when converting to a `Structure`.
    infer_residues : bool, optional, default=False
        Attempt to assign residues based on names of children.

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

    # Attempt to grab residue names based on names of children
    if not residues and infer_residues:
        residues = list(set([child.name for child in compound.children]))

    if isinstance(residues, str):
        residues = [residues]
    if isinstance(residues, (list, set)):
        residues = tuple(residues)

    default_residue = pmd.Residue('RES')
    port_residue = pmd.Residue('PRT')
    compound_residue_map = dict()
    atom_residue_map = dict()

    # Loop through particles and add initialize ParmEd atoms
    for atom in compound.particles(include_ports=show_ports):
        if atom.port_particle:
            current_residue = port_residue
            atom_residue_map[atom] = current_residue

            if current_residue not in structure.residues:
                structure.residues.append(current_residue)

            pmd_atom = pmd.Atom(atomic_number=0, name='VS',
                                mass=0, charge=0)
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10  # Angstroms

        else:
            if residues and atom.name in residues:
                current_residue = pmd.Residue(atom.name)
                atom_residue_map[atom] = current_residue
                compound_residue_map[atom] = current_residue
            elif residues:
                for parent in atom.ancestors():
                    if residues and parent.name in residues:
                        if parent not in compound_residue_map:
                            current_residue = pmd.Residue(parent.name)
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

            atomic_number = None
            name = ''.join(char for char in atom.name if not char.isdigit())
            try:
                atomic_number = AtomicNum[atom.name.capitalize()]
            except KeyError:
                element = element_by_name(atom.name.capitalize())
                if name not in guessed_elements:
                    warn(
                        'Guessing that "{}" is element: "{}"'.format(
                            atom, element))
                    guessed_elements.add(name)
            else:
                element = atom.name.capitalize()

            atomic_number = atomic_number or AtomicNum[element]
            mass = Mass[element]
            pmd_atom = pmd.Atom(
                atomic_number=atomic_number,
                name=atom.name,
                mass=mass,
                charge=atom.charge
            )
            pmd_atom.xx, pmd_atom.xy, pmd_atom.xz = atom.pos * 10  # Angstroms

        residue = atom_residue_map[atom]
        structure.add_atom(
            pmd_atom,
            resname=residue.name,
            resnum=residue.idx
        )

        atom_mapping[atom] = pmd_atom

    # "Claim" all of the items it contains and subsequently index all of its items
    structure.residues.claim()

    # Create and add bonds to ParmEd Structure
    for atom1, atom2 in compound.bonds():
        bond = pmd.Bond(atom_mapping[atom1], atom_mapping[atom2])
        structure.bonds.append(bond)
    # pad box with .25nm buffers
    if box is None:
        box = compound.boundingbox
        box_vec_max = box.maxs.tolist()
        box_vec_min = box.mins.tolist()
        for dim, val in enumerate(compound.periodicity):
            if val:
                box_vec_max[dim] = val
                box_vec_min[dim] = 0.0
            if not val:
                box_vec_max[dim] += 0.25
                box_vec_min[dim] -= 0.25
        box = Box(mins=box_vec_min, maxs=box_vec_max)

    box_vector = np.empty(6)
    if box.angles is not None:
        box_vector[3:6] = box.angles
    else:
        box_vector[3] = box_vector[4] = box_vector[5] = 90.0
    for dim in range(3):
        box_vector[dim] = box.lengths[dim] * 10
    structure.box = box_vector
    return structure


def to_trajectory(compound,
                  show_ports=False,
                  chains=None,
                  residues=None,
                  box=None):
    """Convert to an md.Trajectory and flatten the compound.

    Parameters
    ----------
    show_ports : bool, optional, default=False
        Include all port atoms when converting to trajectory.
    chains : mb.Compound or list of mb.Compound
        Chain types to add to the topology
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by
        checking against Compound.name.
    box : mb.Box, optional, default=compound.boundingbox (with buffer)
        Box information to be used when converting to a `Trajectory`.
        If 'None', a bounding box is used with a 0.5nm buffer in each
        dimension. to avoid overlapping atoms, unless `compound.periodicity`
        is not None, in which case those values are used for the
        box lengths.

    Returns
    -------
    trajectory : md.Trajectory

    See also
    --------
    _to_topology

    """
    md = import_('mdtraj')
    atom_list = [particle for particle in compound.particles(show_ports)]

    top = _to_topology(compound, atom_list, chains, residues)

    # Coordinates.
    xyz = np.ndarray(shape=(1, top.n_atoms, 3), dtype='float')
    for idx, atom in enumerate(atom_list):
        xyz[0, idx] = atom.pos

    # Unitcell information.
    unitcell_angles = [90.0, 90.0, 90.0]
    if box is None:
        unitcell_lengths = np.empty(3)
        for dim, val in enumerate(compound.periodicity):
            if val:
                unitcell_lengths[dim] = val
            else:
                unitcell_lengths[dim] = compound.boundingbox.lengths[dim] + 0.5
    else:
        unitcell_lengths = box.lengths
        unitcell_angles = box.angles

    return md.Trajectory(
        xyz,
        top,
        unitcell_lengths=unitcell_lengths,
        unitcell_angles=unitcell_angles
    )


def _to_topology(compound,
                 atom_list,
                 chains=None,
                 residues=None):
    """Create a mdtraj.Topology from a Compound.

    Helper function for to_trajectory.
    Parameters
    ----------
    atom_list : list of mb.Compound
        Atoms to include in the topology
    chains : mb.Compound or list of mb.Compound
        Chain types to add to the topology
    residues : str of list of str
        Labels of residues in the Compound. Residues are assigned by
        checking against Compound.name.

    Returns
    -------
    top : mdtraj.Topology

    See Also
    --------
    mdtraj.Topology : Details on the mdtraj Topology object

    """
    md = import_('mdtraj')
    from mdtraj.core.topology import Topology
    from mdtraj.core.element import get_by_symbol

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
    default_residue = top.add_residue('RES', default_chain)

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
                                'RES', current_chain)
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
                                parent.name, current_chain)
                            compound_residue_map[parent] = current_residue
                        break
                else:
                    current_residue = default_residue
        else:
            if chains:
                try:  # Grab the default residue from the custom chain.
                    current_residue = next(current_chain.residues)
                except StopIteration:  # Add the residue to the current chain
                    current_residue = top.add_residue('RES', current_chain)
            else:  # Grab the default chain's default residue
                current_residue = default_residue
        atom_residue_map[atom] = current_residue

        # Add the actual atoms
        try:
            elem = get_by_symbol(atom.name)
        except KeyError:
            elem = get_by_symbol("VS")
        at = top.add_atom(atom.name, elem, atom_residue_map[atom])
        at.charge = atom.charge
        atom_mapping[atom] = at

    # Remove empty default residues.
    chains_to_remove = [
        chain for chain in top.chains if chain.n_atoms == 0]
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


def to_pybel(compound,
             box=None,
             title='',
             residues=None,
             show_ports=False,
             infer_residues=False):
    """Create a pybel.Molecule from a Compound

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that need to be converted.
    box : mb.Box, optional, default=None
    title : str, optional, default=compound.name
        Title/name of the ParmEd Structure
    residues : str of list of str
        Labels of residues in the Compound. Residues are
        assigned by checking against Compound.name.
    show_ports : boolean, optional, default=False
        Include all port atoms when converting to a `Structure`.
    infer_residues : bool, optional, default=False
        Attempt to assign residues based on names of children

    Returns
    -------
    pybelmol : pybel.Molecule

    Notes
    -----
    Most of the mb.Compound is first converted to openbabel.OBMol
    And then pybel creates a pybel.Molecule from the OBMol
    Bond orders are assumed to be 1
    OBMol atom indexing starts at 1, with spatial dimension Angstrom
    """

    openbabel = import_('openbabel')
    pybel = import_('pybel')

    mol = openbabel.OBMol()
    particle_to_atom_index = {}

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
            try:
                temp.SetAtomicNum(AtomicNum[part.name.capitalize()])
            except KeyError:
                warn("Could not infer atomic number from "
                     "{}, setting to 0".format(part.name))
                temp.SetAtomicNum(0)

        temp.SetVector(*(part.xyz[0]*10))
        particle_to_atom_index[part] = i

    ucell = openbabel.OBUnitCell()
    if box is None:
        box = compound.boundingbox
    a, b, c = 10.0 * box.lengths
    alpha, beta, gamma = np.radians(box.angles)

    cosa = np.cos(alpha)
    cosb = np.cos(beta)
    sinb = np.sin(beta)
    cosg = np.cos(gamma)
    sing = np.sin(gamma)
    mat_coef_y = (cosa - cosb * cosg) / sing
    mat_coef_z = np.power(sinb, 2, dtype=float) - \
                np.power(mat_coef_y, 2, dtype=float)

    if mat_coef_z > 0.:
        mat_coef_z = np.sqrt(mat_coef_z)
    else:
        raise Warning(
            'Non-positive z-vector. Angles {} '
            'do not generate a box with the z-vector in the '
            'positive z direction'.format(box.angles)
        )

    box_vec = [[1, 0, 0], [cosg, sing, 0], [cosb, mat_coef_y, mat_coef_z]]
    box_vec = np.asarray(box_vec)
    box_mat = (np.array([a, b, c]) * box_vec.T).T
    first_vector = openbabel.vector3(*box_mat[0])
    second_vector = openbabel.vector3(*box_mat[1])
    third_vector = openbabel.vector3(*box_mat[2])
    ucell.SetData(first_vector, second_vector, third_vector)
    mol.CloneData(ucell)

    for bond in compound.bonds():
        bond_order = 1
        mol.AddBond(
            particle_to_atom_index[bond[0]]+1,
            particle_to_atom_index[bond[1]]+1,
            bond_order
        )

    pybelmol = pybel.Molecule(mol)
    pybelmol.title = title if title else compound.name

    return pybelmol


def to_networkx(compound, names_only=False):
    """Create a NetworkX graph representing the hierarchy of a Compound.

    Parameters
    ----------
    compound : mb.Compound
        The mbuild Compound that need to be converted.
    names_only : bool, optional, default=False
        Store only the names of the compounds in the graph,
        appended with their IDs, for distinction even if they
        have the same name. When set to False, the default
        behavior, the nodes are the compounds themselves.

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
    nx = import_('networkx')

    nodes = list()
    edges = list()
    if names_only:
        nodes.append(compound.name + '_' + str(id(compound)))
    else:
        nodes.append(compound)
    nodes, edges = _iterate_children(compound, nodes,
                                edges, names_only=names_only)

    graph = nx.DiGraph()
    graph.add_nodes_from(nodes)
    graph.add_edges_from(edges)
    return graph


def _iterate_children(compound, nodes, edges, names_only=False):
    """ Create nodes and edges that connect parents and their corresponding children

    Helper function for to_networkx
    """
    if not compound.children:
        return nodes, edges
    for child in compound.children:
        if names_only:
            unique_name = child.name + '_' + str(id(child))
            unique_name_parent = child.parent.name + '_' + str((id(child.parent)))
            nodes.append(unique_name)
            edges.append([unique_name_parent, unique_name])
        else:
            nodes.append(child)
            edges.append([child.parent, child])
        nodes, edges = _iterate_children(child, nodes, edges, names_only=names_only)
    return nodes, edges


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
    from intermol.atom import Atom as InterMolAtom
    from intermol.molecule import Molecule
    from intermol.system import System
    import simtk.unit as u

    if isinstance(molecule_types, list):
        molecule_types = tuple(molecule_types)
    elif molecule_types is None:
        molecule_types = (type(compound),)
    intermol_system = System()

    last_molecule_compound = None
    for atom_index, atom in enumerate(compound.particles()):
        for parent in atom.ancestors():
            # Don't want inheritance via isinstance().
            if type(parent) in molecule_types:
                # Check if we have encountered this molecule type before.
                if parent.name not in intermol_system.molecule_types:
                    compound._add_intermol_molecule_type(
                        intermol_system, parent)
                if parent != last_molecule_compound:
                    last_molecule_compound = parent
                    last_molecule = Molecule(name=parent.name)
                    intermol_system.add_molecule(last_molecule)
                break
        else:
            # Should never happen if molecule_types only contains
            # type(compound)
            raise ValueError(
                'Found an atom {} that is not part of any of '
                'the specified molecule types {}'.format(atom, molecule_types)
            )

        # Add the actual intermol atoms.
        intermol_atom = InterMolAtom(atom_index + 1, name=atom.name,
                                     residue_index=1, residue_name='RES')
        intermol_atom.position = atom.pos * u.nanometers
        last_molecule.add_atom(intermol_atom)
    return intermol_system


def _add_intermol_molecule_type(intermol_system, parent):  # pragma: no cover
    """Create a molecule type for the parent and add bonds.

    This method takes an intermol system and adds a
    parent compound, including its particles and bonds, to it.
    """
    from intermol.moleculetype import MoleculeType
    from intermol.forces.bond import Bond as InterMolBond

    molecule_type = MoleculeType(name=parent.name)
    intermol_system.add_molecule_type(molecule_type)

    for index, parent_atom in enumerate(parent.particles()):
        parent_atom.index = index + 1

    for atom1, atom2 in parent.bonds():
        intermol_bond = InterMolBond(atom1.index, atom2.index)
        molecule_type.bonds.add(intermol_bond)
