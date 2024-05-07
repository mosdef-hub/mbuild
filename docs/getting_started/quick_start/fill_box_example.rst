Fill Box
========


All-Atom (AA) Hexane and Ethanol System
---------------------------------------

.. note::
    `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with ``mBuild`` in
    the following example to demonstrate how the `MoSDeF <https://mosdef.org>`_
    libraries can be used together to generate a simulation box.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct an all-atom (AA) hexane and ethanol using the OPLS-AA force field (FF),
which is shipped as a standard `foyer <https://foyer.mosdef.org/en/stable/>`_ FF.
The hexane and ethanol molecules will be created using `smiles strings <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_.
The hexane and ethanol residues will be named `"HEX"` and `"ETO"`, respectively.
Lastly, the hexane and ethanol molecule's configuration will be energy minimized, properly reorienting the molecule to the specified FF, which is sometimes needed for some simulation engines to ensure the initial configuration energy is not too high.

.. note::
    The energy minimize step requires the `foyer <https://foyer.mosdef.org/en/stable/>`_ package.

.. code:: ipython3

    hexane = mb.load('CCCCCC', smiles=True)
    hexane.name = 'HEX'
    hexane.energy_minimize(forcefield='oplsaa', steps=10**4)


    ethanol = mb.load('CCO', smiles=True)
    ethanol.name = 'ETO'
    ethanol.energy_minimize(forcefield='oplsaa', steps=10**4)

The liquid box is built to a density of 680 kg/m\ :sup:`3`, with a 50/50 mol ratio of hexane and ethanol, and will be in an orthogonal box measuring 5.0 nm in the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound= [hexane, ethanol],
                          density=680,
                          compound_ratio=[0.5, 0.5],
                          box=[5.0, 5.0, 5.0])


United Atom (UA) Methane System
-------------------------------

.. note::
    `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with ``mBuild`` in
    the following example to demonstrate how the `MoSDeF <https://mosdef.org>`_ libraries
    integrate to generate a simulation box.  A subset of the `TraPPE-United Atom <http://trappe.oit.umn.edu>`_
    force field (FF) comes standard with the `foyer <https://foyer.mosdef.org/en/stable/>`_ software package.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct a pseudo-monatomic molecule (united atom (UA) methane), for use with the
`TraPPE <http://trappe.oit.umn.edu>`_ FF.  The UA methane, bead type `"_CH4"`, will be built as a child (``mbuild.Compound.children``), so the parent (``mbuild.Compound``) will
allow a user-selected residue name (``mbuild.Compound.name``). If the methane is built using ``methane = mb.Compound(name="_CH4")``, then the user must keep the residue name `"_CH4"` or `foyer <https://foyer.mosdef.org/en/stable/>`_ will not recognize the bead type when using the standard TraPPE force field XML file.

.. code:: ipython3

    methane = mb.Compound(name="MET")
    methane_child_bead = mb.Compound(name="_CH4")
    methane.add(methane_child_bead, inherit_periodicity=False)

.. note::
    The ``inherit_periodicity`` flag is an optional boolean (default=True), which replaces
    the periodicity of self with the periodicity of the Compound being added.

The orthogonal liquid box contains 1230 methane molecules and measures 4.5 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound=methane,
                          n_compounds=1230,
                          box=[4.5, 4.5, 4.5]
                          )
