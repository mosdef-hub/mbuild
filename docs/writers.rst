=======
Writers
=======

mBuild utilizes ParmEd to write ``Compound`` information to a variety of file
formats (e.g. PDB, MOL2, GRO.  The full list of formats supported by ParmEd
can be found at http://parmed.github.io/ParmEd/html/readwrite.html).
Additionally, mBuild features several internal writers for file formats not yet
supported by ParmEd. Information on these internal writers can be found below.

By default, mBuild will only write coordinate and bond information to these files,
i.e. no angles or dihedrals, and no atom typing is performed (atom names are used
as atom types). However, force fields can be applied to Compounds by passing force
field XML files (used by the Foyer package - https://github.com/mosdef-hub/foyer) 
to the ``save`` function if Foyer is installed. If a force field is applied to a 
Compound, the mBuild internal writers will also write angle and dihedral information 
to the file in addition to labelling atoms by the atom types specified by the force 
field.

GSD (General Simulation Data)
-----------------------------
Default data file format for HOOMD v2
+++++++++++++++++++++++++++++++++++++

.. automodule:: mbuild.formats.gsdwriter
    :members:

HOOMD XML
---------
Default data file format for HOOMD v1
+++++++++++++++++++++++++++++++++++++

.. automodule:: mbuild.formats.hoomdxml
    :members:

LAMMPS data
-----------
Default data file format for LAMMPS
+++++++++++++++++++++++++++++++++++

.. automodule:: mbuild.formats.lammpsdata
    :members:
