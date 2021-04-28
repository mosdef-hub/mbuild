Example of an All-Atom (AA) Hexane and Ethanol System
========================

**Note**: `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with **mBuild** in the following example to demonstrate how all the **MoSDeF** software integrates to generate all the required to execute simulations.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct an all-atom (AA) hexane and ethanol using the OPLS-AA force field (FF),
which is comes as a standard `foyer <https://foyer.mosdef.org/en/stable/>`_ force field (FF).
The hexane and ethanol residues will be named HEX and ETO, respectively.  
Lastly, the hexane and ethanol molecule's configuration will be energy minimized. 

**Note**: The energy minimize step does also require the `foyer <https://foyer.mosdef.org/en/stable/>`_ package. 

.. code:: ipython3

    hexane = mb.load('CCCCCC', smiles=True, name='HEX')
    hexane.name = 'HEX'
    hexane.energy_minimize(forcefield='oplsaa' , steps=10**4)


    ethanol = mb.load('CCO', smiles=True)
    ethanol.name = 'ETO'
    ethanol.energy_minimize(forcefield='oplsaa' , steps=10**4)


Build the liquid box designed for 300 K using and isobaric-isothermal (NPT) simulation.
The liquid box is built to a density of 680 kg/m^3, with a 50/50 mol ratio of hexane and ethanol, 
and will be an orthogonal box measuring 5.0 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound= [hexane, ethanol],
                          density=680,
                          compound_ratio=[0.5, 0.5],
                          box=[5.0, 5.0, 5.0])



Writing all the GOMC files
----------------------------

**Note**: This portion of the code is specific to building the files for GPU Optimized Monte Carlo 
(`GOMC <http://gomc.eng.wayne.edu>`_), and requires the `foyer <https://foyer.mosdef.org/en/stable/>`_ package.

First, import the required GOMC packages.

.. code:: ipython3

    import mbuild.formats.charmm_writer as mf_charmm
    import mbuild.formats.gomc_conf_writer as gomc_control


The mbuild.Charmm object is built.  Then, the force field/parameter file, 
PSF, and PDB files are output, in that order, from the mbuild.Charmm object.


.. code:: ipython3

    charmm = mf_charmm.Charmm(box_liq,
                              'NPT_n_hexane_ethanol_liq',
                              ff_filename ="NPT_n_hexane_ethanol_FF",
                              forcefield_selection='oplsaa',
                              residues=[hexane.name, ethanol.name],
                              )

    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()


**Note**: This portion of the code is specific to building the files for GPU Optimized Monte Carlo 
(`GOMC <http://gomc.eng.wayne.edu>`_).


The `GOMC <http://gomc.eng.wayne.edu>`_ control file is written from the mbuild.Charmm object using the NPT ensemble 
with a pressure of 10 bar.  
The first five required variables are the mbuild.Charmm object, desired GOMC control file name, ensemble type, number of steps, and the temperature (K).  
The `GOMC <http://gomc.eng.wayne.edu>`_ input variables (input_variables_dict) are optional, except for the grand canonical ensemble (GCMC), which requires the input for the chemical potential or fugacity. 

**Note**: Most input_variables_dict keys are the same as the GOMC Manual commands or may 
have "_box_0" or "_box_1" added at the end of the `GOMC <http://gomc.eng.wayne.edu>`_ 
Manual naming convention. 


.. code:: ipython3

    gomc_control.write_gomc_control_file(charmm, 'in_NPT.conf', 'NPT', 100, 300,
                                         input_variables_dict={"Pressure": 10}
                                         )
