Example of an All-Atom (AA) Hexane and Ethanol System
========================

**Note**: `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with **mBuild** in the following example to demonstrate how all the `MoSDeF <https://mosdef.org>`_ software integrates to generate a simulation box.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct an all-atom (AA) hexane and ethanol using the OPLS-AA force field (FF),
which is comes as a standard `foyer <https://foyer.mosdef.org/en/stable/>`_ force field (FF).
The hexane and ethanol residues will be named HEX and ETO, respectively.  
Lastly, the hexane and ethanol molecule's configuration will be energy minimized. 

**Note**: The energy minimize step requires the `foyer <https://foyer.mosdef.org/en/stable/>`_ package. 

.. code:: ipython3

    hexane = mb.load('CCCCCC', smiles=True, name='HEX')
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
