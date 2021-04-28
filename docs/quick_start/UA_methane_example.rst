Example System Containing United Atom (UA) Methane
========================

**Note**: `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with **mBuild** in the following example to demonstrate how all the **MoSDeF** software integrates to generate all the required to execute simulations.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct a pseudo-monatomic molecule (united atom (UA) methane) using the 
TraPPE force field (FF).  The UA methane, bead type "_CH4", will be built as a child (mbuild.Compound.children), so the parent (mbuild.Compound) will 
allow a user-selected residue name (mbuild.Compound.name). Otherwise, if the 
parent methane bead is not "_CH4", the standard 
`foyer <https://foyer.mosdef.org/en/stable/>`_  TraPPE FF will not 
recognize the bead.

.. code:: ipython3

    methane = mb.Compound(name="MET")
    methane_child_bead = mb.Compound(name="_CH4")
    methane.add(methane_child_bead, inherit_periodicity=False)


Build the liquid and vapor boxes for methane at 149.14 K, which are designed to calculate
the liquid-vapor phase equilibrium using the Gibbs ensemble Monte Carlo (GEMC) simulation.
The orthogonal liquid box contains 1230 methane molecules and measures 4.5 nm in all the x, y, and z-dimensions.
The orthogonal vapor box contains 302 methane molecules and measures 8.0 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound=methane,
                          n_compounds=1230,
                          box=[4.5, 4.5, 4.5]
                          )

    box_vap = mb.fill_box(compound=methane,
                          n_compounds=302,
                          box=[8.0, 8.0, 8.0]
                          )


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
                              'GEMC_NVT_methane_liq_box_0',
                              structure_box_1=box_vap,
                              filename_box_1='GEMC_NVT_methane_vap_box_1',
                              ff_filename="GEMC_NVT_methane_FF",
                              forcefield_selection='trappe-ua',
                              residues=[methane.name],
                              bead_to_atom_name_dict={'_CH4':'C'},
                              )

    charmm.write_inp()

    charmm.write_psf()

    charmm.write_pdb()


**Note**: This portion of the code is specific to building the files for GPU Optimized Monte Carlo 
(`GOMC <http://gomc.eng.wayne.edu>`_).

The `GOMC <http://gomc.eng.wayne.edu>`_ control file is written from the mbuild.Charmm object using the GEMC-NVT ensemble.  
The first five required variables are the mbuild.Charmm object, desired GOMC control file name, ensemble type, number of steps, and the temperature (K).  
The `GOMC <http://gomc.eng.wayne.edu>`_ input variables (input_variables_dict) are optional, except for the grand canonical ensemble (GCMC), which requires the input for the chemical potential or fugacity. 

**Note**: Most input_variables_dict keys are the same as the GOMC Manual commands or may 
have "_box_0" or "_box_1" added at the end of the `GOMC <http://gomc.eng.wayne.edu>`_ 
Manual naming convention. 


.. code:: ipython3

    gomc_control.write_gomc_control_file(charmm, 'in_GEMC_NVT.conf', 'GEMC_NVT', 100, 149.14,
                                         input_variables_dict={"Potential" : 'SWITCH',
                                                               "Rswitch" : 10 ,
                                                               "Rcut": 12,
                                                               "RcutLow": 1,
                                                               }
                                         )

