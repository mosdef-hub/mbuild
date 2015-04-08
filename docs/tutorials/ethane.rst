Ethane: Reading from files, Ports and coordinate transforms
-----------------------------------------------------------

In this example, we'll cover reading molecular components from files, introduce
the concept of ``Ports`` and start using some coordinate transforms.

As you probably noticed while creating your methane mocule in the last tutorial,
manually adding ``Atoms`` and ``Bonds`` to a ``Compound`` is a bit cumbersome.
The easiest way to create small, reusable components, such as methyls, amines or
monomers, is to hand draw them using software like `Avogadro <http://avogadro.cc/wiki/Main_Page>`_
and export them as either a .pdb or .mol2 file (the file should contain
connectivity information).

Let's start by reading a methyl group from a .pdb file::

    import mbuild as mb

    class CH3(mb.Compound):
        """A methyl group. """
        def __init__(self):
            super(CH3, self).__init__(self)

            mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)

.. note:: The ``relative_to_module`` argument allows you to look for the file
          in directories relative to the current module, rather than where your
          main script is being executed. This is useful because we have a
          collection of components in ``mbuild/components`` where we store
          files like the one that we're creating along with the corresponding
          .pdb files.

Now let's use our first coordinate transform to center the methyl at its carbon
atom::

   import mbuild as mb


    class CH3(mb.Compound):

        def __init__(self):
            super(CH3, self).__init__()

            mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
            mb.translate(self, -self.C[0])  # Move carbon to origin.

Here, when referring to the ``Atom`` object in a mathematical expression, we
operate directly on its coordinates. This functionality is currently implemented
for addition and subtraction only.

So now we have a methyl group loaded up and centered. In order to connect
``Compounds`` in mBuild, we make use of a special type of ``Compound``: the ``Port``.
A ``Port`` is a ``Compound`` with two sets of four "ghost" ``Atoms`` (of kind
'G' by convention). In addition ``Ports`` have an ``anchor`` attribute which
typically points to an ``Atom`` that the ``Port`` should be associated with. In
our methyl group, the ``Port`` should be anchored to the carbon atom so that we
can now form ``Bonds`` to this carbon::


    import mbuild as mb


    class CH3(mb.Compound):

        def __init__(self):
            super(CH3, self).__init__()

            mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
            mb.translate(self, -self.C[0])  # Move carbon to origin.

            port = mb.Port(anchor=self.C[0])
            self.add(port, label='up')

By default, ``Ports`` are created to point in the y-direction. Typically,
``Ports`` should point in the direction of the ``Bond`` we want to create.
Since the methyl group is already oriented towards the y-direction, let's simply
move the port a tiny bit away from the carbon::

    import mbuild as mb


    class CH3(mb.Compound):

        def __init__(self):
            super(CH3, self).__init__()

            mb.load('ch3.pdb', compound=self, relative_to_module=self.__module__)
            mb.translate(self, -self.C[0])  # Move carbon to origin.

            port = mb.Port(anchor=self.C[0])
            self.add(port, label='up')
            mb.translate(self.up, [0, -0.07, 0])

By default, ``Ports`` are never output from the mBuild structure. However,
it can be useful to look at a molecule with the ``Ports`` to check your work as
you go. The default VMD settings don't show ports very well so let's change
the representation style to small vdW spheres:

.. image:: ../images/methyl_port.png
    :align: center
    :scale: 50%
    :alt: Methyl group with one set of the four ghost atoms contained in the Port.

When two ``Ports`` are connected, they are forced to overlap in space and their
parent ``Compounds`` are rotated and translated by the same amount.
If we tried to connect two of our methyls right now using only one set of four
ghost atoms, not only would the ``Ports`` overlap perfectly, but the carbons and
hydrogens would also perfectly overlap. This is why every port contains a second
set of 4 ghost atoms pointing in the opposite direction. When two ``Compounds`` are
connected, the port that places the anchor atoms the farthest away from each other
is chosen automatically to prevent this overlap scenario.

.. note:: By convention, we try to label ``Ports`` successively as 'up', 'down', 'left',
          'right', 'front', 'back' which should roughly correspond to the relative
          This is a bit tricky to enforce because the system is so flexible so
          use your best judgement and try to be consistent! The more components
          we collect in our library with the same labeling conventions, the
          easier it becomes to build ever more complex structures.

Now the fun part: stick 'em together to create an ethane::

    import mbuild as mb


    class Ethane(mb.Compound):
        """An ethane molecule. """
        def __init__(self):
            """Connect two methyl groups to form an ethane. """
            super(Ethane, self).__init__()

            self.add(Ch3(), "methyl1")
            self.add(Ch3(), "methyl2")
            equivalence_transform(self.methyl1, self.methyl1.up, self.methyl2.up)

    if __name__ == '__main__':
        ethane = Ethane()
        ethane.visualize(show_ports=True)



.. image:: ../images/ethane.png
    :align: center
    :scale: 50%
    :alt: Ethane with all Ports shown.

The ``equivalence_transform()`` function takes a ``Compound`` and then rotates
and translates it such that two other ``Compounds`` overlap. Typically, as in
this case, those two other ``Compounds`` are ``Ports`` - in our case, ``methyl1.up``
and ``methyl2.up``.

