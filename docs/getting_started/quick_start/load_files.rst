.. _QuickStart_Load_files:

Load files
========================


mol2 files
------------------------

Create an ``mbuild.Compound`` (i.e., the "pentane" variable) by loading a molecule from a `mol2 <http://chemyang.ccnu.edu.cn/ccb/server/AIMMS/mol2.pdf>`_ file.

Import the required mbuild packages.

.. code:: ipython3

    import mbuild as mb


Load the "pentane.mol2" file from its directory.

.. code:: ipython3

    pentane = mb.load("path_to_mol2_file/pentane.mol2")


CIF files
------------------------

Build an ``mbuild.Compound`` (i.e., the "ETV_triclinic" variable) by loading a `Crystallographic Information File (CIF) <https://www.iucr.org/resources/cif>`_ file and selecting the number of cell units to populate in the x, y, and z-dimensions.


Import the required mbuild packages.

.. code:: ipython3

    import mbuild as mb
    from mbuild.lattice import load_cif



The `CIF <https://www.iucr.org/resources/cif>`_ file is loaded using the ``load_cif`` function. Next, three (3) cell units shall be built for all the x, y, and z-dimensions with the populate function.  Finally, the `CIF <https://www.iucr.org/resources/cif>`_'s residues are named 'ETV'.

.. code:: ipython3

    lattice_cif_ETV_triclinic = load_cif("path_to_cif_file/ETV_triclinic.cif")
    ETV_triclinic = lattice_cif_ETV_triclinic.populate(x=3, y=3, z=3)
    ETV_triclinic.name = 'ETV'


Other file types
------------------------
mBuild also supports :ref:`loading_data` or files via hoomd_snapshot, GSD, SMILES strings, and ParmEd structures.
