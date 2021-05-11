Fill Box Example
========================


Example of an All-Atom (AA) Hexane and Ethanol System
------------------------

**Note**: `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with **mBuild** in the following example to demonstrate how all the `MoSDeF <https://mosdef.org>`_ software integrates to generate a simulation box.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct an all-atom (AA) hexane and ethanol using the OPLS-AA force field (FF),
which is comes as a standard `foyer <https://foyer.mosdef.org/en/stable/>`_ force field (FF).
The hexane and ethanol molecules will be created using `smiles strings <https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html>`_.
The hexane and ethanol residues will be named HEX and ETO, respectively.
Lastly, the hexane and ethanol molecule's configuration will be energy minimized, properly reorienting the molecule to the specified FF, which is sometimes needed for some simulation engines to ensure the initial configuration energy is not too high.

**Note**: The energy minimize step requires the `foyer <https://foyer.mosdef.org/en/stable/>`_ package.

.. code:: ipython3

    hexane = mb.load('CCCCCC', smiles=True)
    hexane.name = 'HEX'
    hexane.energy_minimize(forcefield='oplsaa' , steps=10**4)


    ethanol = mb.load('CCO', smiles=True)
    ethanol.name = 'ETO'
    ethanol.energy_minimize(forcefield='oplsaa' , steps=10**4)


The liquid box is built to a density of 680 kg/m^3, with a 50/50 mol ratio of hexane and ethanol,
and will be an orthogonal box measuring 5.0 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound= [hexane, ethanol],
                          density=680,
                          compound_ratio=[0.5, 0.5],
                          box=[5.0, 5.0, 5.0])
