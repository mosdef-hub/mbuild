Ethane: Reading from files, Ports and coordinate transforms
-----------------------------------------------------------

**Note**: mBuild expects all distance units to be in nanometers.

In this example, we’ll cover reading molecular components from files,
introduce the concept of ``Ports`` and start using some coordinate
transforms.

First, we need to import the mbuild package:

.. code:: ipython3

    import mbuild as mb

As you probably noticed while creating your methane molecule in the last
tutorial, manually adding ``Particles`` and ``Bonds`` to a ``Compound``
is a bit cumbersome. The easiest way to create small, reusable
components, such as methyls, amines or monomers, is to hand draw them
using software like `Avogadro <https://avogadro.cc/>`__ and
export them as either a .pdb or .mol2 file (the file should contain
connectivity information).

Let’s start by reading a methyl group from a ``.pdb`` file:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.visualize()


Now let’s use our first coordinate transform to center the methyl at its
carbon atom:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.translate(-ch3[0].pos)  # Move carbon to origin.

Now we have a methyl group loaded up and centered. In order to connect
``Compounds`` in mBuild, we make use of a special type of ``Compound``:
the ``Port``. A ``Port`` is a ``Compound`` with two sets of four “ghost”
``Particles`` that assist in bond creation. In addition, ``Ports`` have an ``anchor`` attribute which
typically points to a particle that the ``Port`` should be associated
with. In our methyl group, the ``Port`` should be anchored to the carbon
atom so that we can now form bonds to this carbon:

.. code:: ipython3

    import mbuild as mb

    ch3 = mb.load('ch3.pdb')
    ch3.translate(-ch3[0].pos)  # Move carbon to origin.

    port = mb.Port(anchor=ch3[0])
    ch3.add(port, label='up')

    # Place the port at approximately half a C-C bond length.
    ch3['up'].translate([0, -0.07, 0])

By default, ``Ports`` are never output from the mBuild structure.
However, it can be useful to look at a molecule with the ``Ports`` to
check your work as you go:

.. code:: ipython3

    ch3.visualize(show_ports=True)

Now we wrap the methyl group into a python class, so that we can reuse
it as a component to build more complex molecules later.

.. code:: ipython3

    import mbuild as mb

    class CH3(mb.Compound):
        def __init__(self):
            super(CH3, self).__init__()

            mb.load('ch3.pdb', compound=self)
            self.translate(-self[0].pos)  # Move carbon to origin.

            port = mb.Port(anchor=self[0])
            self.add(port, label='up')
            # Place the port at approximately half a C-C bond length.
            self['up'].translate([0, -0.07, 0])

When two ``Ports`` are connected, they are forced to overlap in space
and their parent ``Compounds`` are rotated and translated by the same
amount.

**Note:** If we tried to connect two of our methyls right now using only
one set of four ghost particles, not only would the ``Ports`` overlap
perfectly, but the carbons and hydrogens would also perfectly overlap -
the 4 ghost atoms in the ``Port`` are arranged identically with respect
to the other atoms. For example, if a ``Port`` and its direction is
indicated by “<-”, forcing the port in <-CH3 to overlap with <-CH3 would
just look like <-CH3 (perfectly overlapping atoms).

To solve this problem, every port contains a second set of 4 ghost atoms
pointing in the opposite direction. When two ``Compounds`` are
connected, the port that places the anchor atoms the farthest away from
each other is chosen automatically to prevent this overlap scenario.

When <->CH3 and <->CH3 are forced to overlap, the CH3<->CH3 is
automatically chosen.

Now the fun part: stick ’em together to create an ethane:

.. code:: ipython3

    ethane = mb.Compound()

    ethane.add(CH3(), label="methyl_1")
    ethane.add(CH3(), label="methyl_2")
    mb.force_overlap(move_this=ethane['methyl_1'],
                             from_positions=ethane['methyl_1']['up'],
                             to_positions=ethane['methyl_2']['up'])

Above, the ``force_overlap()`` function takes a ``Compound`` and then
rotates and translates it such that two other ``Compounds`` overlap.
Typically, as in this case, those two other ``Compounds`` are ``Ports``
- in our case, ``methyl1['up']`` and ``methyl2['up']``.

.. code:: ipython3

    ethane.visualize()

.. code:: ipython3

    ethane.visualize(show_ports=True)

Similarly, if we want to make ethane a reusable component, we need to
wrap it into a python class.

.. code:: ipython3

    import mbuild as mb

    class Ethane(mb.Compound):
        def __init__(self):
            super(Ethane, self).__init__()

            self.add(CH3(), label="methyl_1")
            self.add(CH3(), label="methyl_2")
            mb.force_overlap(move_this=self['methyl_1'],
                             from_positions=self['methyl_1']['up'],
                             to_positions=self['methyl_2']['up'])

.. code:: ipython3

    ethane = Ethane()
    ethane.visualize()

.. code:: ipython3

    # Save to .mol2
    ethane.save('ethane.mol2', overwrite=True)
