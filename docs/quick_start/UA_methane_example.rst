Example System Containing United Atom (UA) Methane
========================

**Note**: `foyer <https://foyer.mosdef.org/en/stable/>`_ is used in conjunction with **mBuild** in the following example to demonstrate how all the `MoSDeF <https://mosdef.org>`_ software integrates to generate a simulation box.  The `TraPPE <http://trappe.oit.umn.edu>`_ force field (FF) comes standard with the `foyer <https://foyer.mosdef.org/en/stable/>`_ software package.

Import the required mbuild package.

.. code:: ipython3

    import mbuild as mb


Construct a pseudo-monatomic molecule (united atom (UA) methane), for use with the
`TraPPE <http://trappe.oit.umn.edu>`_ force field (FF).  The UA methane, bead type "_CH4", will be built as a child (mbuild.Compound.children), so the parent (mbuild.Compound) will 
allow a user-selected residue name (mbuild.Compound.name). If the methane is built using methane = mb.Compound(name="_CH4"), then the user must keep the residue name "_CH4" or `foyer <https://foyer.mosdef.org/en/stable/>`_ will not recognize the bead type.

.. code:: ipython3

    methane = mb.Compound(name="MET")
    methane_child_bead = mb.Compound(name="_CH4")
    methane.add(methane_child_bead, inherit_periodicity=False)


The orthogonal liquid box contains 1230 methane molecules and measures 4.5 nm in all the x, y, and z-dimensions.

.. code:: ipython3

    box_liq = mb.fill_box(compound=methane,
                          n_compounds=1230,
                          box=[4.5, 4.5, 4.5]
                          )
