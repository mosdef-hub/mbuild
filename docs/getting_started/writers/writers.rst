----------------------------------------------
File Writers
----------------------------------------------

The mBuild library also supports simulation engine-specific file writers.  These writers create a complete set of simulation writers to input files or a partial set of file writers, where the other required files are generated via another means.

mBuild utilizes ParmEd to write ``Compound`` information to a variety of file
formats (e.g. PDB, MOL2, GRO.  The full list of formats supported by ParmEd
can be found at the `ParmEd website <http://parmed.github.io/ParmEd/html/readwrite.html>`_).
Additionally, mBuild features several internal writers for file formats not yet
supported by ParmEd. Information on these internal writers can be found below.

By default, many mBuild functions will only write coordinate and bond information to these files,
i.e. no angles or dihedrals, and no atom typing is performed (atom names are used
as atom types). However, force fields can be applied to Compounds by passing force
field XML files (used by the `Foyer package <https://github.com/mosdef-hub/foyer>`_)
to the ``Compound.save`` function if Foyer is installed. If a force field is applied to a
Compound, the mBuild internal writers will also write angle and dihedral information
to the file in addition to labelling atoms by the atom types specified by the force
field.  The CHARMM-style GOMC writers (supported through the `MoSDeF-GOMC extension <https://github.com/GOMC-WSU/MoSDeF-GOMC>`_) are the exception to this default rule since
they need a force field to build the files, as these files depend on the force field parameters (Example: charge and MW in the PSF files).

The simulation engine writers that use mBuild or are currently contained in the mBuild library:


* `Cassandra <https://cassandra.nd.edu/>`_
* `GROMACS <https://www.gromacs.org/>`_
* `HOOMD-blue <http://glotzerlab.engin.umich.edu/hoomd-blue//>`_
* `Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS) <https://lammps.sandia.gov/>`_

Support for `GPU Optimized Monte Carlo (GOMC) <http://gomc.eng.wayne.edu/>`_  is also available through the `MoSDeF-GOMC library <https://github.com/GOMC-WSU/MoSDeF-GOMC>`_

.. toctree::

	cassandra_file_writers
   	HOOMD_blue_file_writers
	LAMMPS_file_writers
