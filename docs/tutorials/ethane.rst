
Ethane: Reading from files, Ports and coordinate transforms
-----------------------------------------------------------

In this example, we'll cover reading molecular components from files,
introduce the concept of ``Ports`` and start using some coordinate
transforms.

As you probably noticed while creating your methane mocule in the last
tutorial, manually adding ``Atoms`` and ``Bonds`` to a ``Compound`` is a
bit cumbersome. The easiest way to create small, reusable components,
such as methyls, amines or monomers, is to hand draw them using software
like `Avogadro <http://avogadro.cc/wiki/Main_Page>`__ and export them as
either a .pdb or .mol2 file (the file should contain connectivity
information).

Let's start by reading a methyl group from a ``.pdb`` file:

.. code:: python

    import mbuild as mb
    
    class CH3(mb.Compound):
        def __init__(self):
            super(CH3, self).__init__(self)
    
            mb.load('ch3.pdb', compound=self)


Now let's use our first coordinate transform to center the methyl at its
carbon atom:

.. code:: python

    import mbuild as mb
    
    class CH3(mb.Compound):
        def __init__(self):
            super(CH3, self).__init__()
    
            mb.load('ch3.pdb', compound=self)
            mb.translate(self, -self.C[0])  # Move carbon to origin.

Here, when referring to the ``Atom`` object in a mathematical
expression, we operate directly on its coordinates. This functionality
is currently implemented for addition and subtraction only.

So now we have a methyl group loaded up and centered. In order to
connect ``Compounds`` in mBuild, we make use of a special type of
``Compound``: the ``Port``. A ``Port`` is a ``Compound`` with two sets
of four "ghost" ``Atoms`` (of kind ``'G'`` by convention). In addition
``Ports`` have an ``anchor`` attribute which typically points to an
``Atom`` that the ``Port`` should be associated with. In our methyl
group, the ``Port`` should be anchored to the carbon atom so that we can
now form ``Bonds`` to this carbon:

.. code:: python

    import mbuild as mb
    
    class CH3(mb.Compound):
        def __init__(self):
            super(CH3, self).__init__()
    
            mb.load('ch3.pdb', compound=self)
            mb.translate(self, -self.C[0])  # Move carbon to origin.
    
            port = mb.Port(anchor=self.C[0])
            self.add(port, label='up')
            mb.translate(self.up, [0, -0.07, 0])

By default, ``Ports`` are never output from the mBuild structure.
However, it can be useful to look at a molecule with the ``Ports`` to
check your work as you go:

.. code:: python

    CH3().visualize(show_ports=True)


.. raw:: html

   <div id="widgetview"></div>
   <script>loadWidget('/assets/json/example-trajectoryview-data.json', '#widgetview', 400, 400);</script>

.. parsed-literal::

    <IPython.core.display.Javascript object>



.. parsed-literal::

    <IPython.core.display.Javascript object>


When two ``Ports`` are connected, they are forced to overlap in space
and their parent ``Compounds`` are rotated and translated by the same
amount.

**Note:** If we tried to connect two of our methyls right now using only
one set of four ghost atoms, not only would the ``Ports`` overlap
perfectly, but the carbons and hydrogens would also perfectly overlap.
This is why every port contains a second set of 4 ghost atoms pointing
in the opposite direction. When two ``Compounds`` are connected, the
port that places the anchor atoms the farthest away from each other is
chosen automatically to prevent this overlap scenario. By convention, we
try to label ``Ports`` successively as 'up', 'down', 'left', 'right',
'front', 'back' which should roughly correspond to the relative. This is
a bit tricky to enforce because the system is so flexible so use your
best judgement and try to be consistent! The more components we collect
in our library with the same labeling conventions, the easier it becomes
to build ever more complex structures.

Now the fun part: stick 'em together to create an ethane:

.. code:: python

    import mbuild as mb
    
    class Ethane(mb.Compound):
        def __init__(self):
            super(Ethane, self).__init__()
    
            self.add(CH3(), "methyl1")
            self.add(CH3(), "methyl2")
            mb.equivalence_transform(self.methyl1, self.methyl1.up, self.methyl2.up)

.. code:: python

    Ethane().visualize(show_ports=True)



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. parsed-literal::

    <IPython.core.display.Javascript object>


Above, the ``equivalence_transform()`` function takes a ``Compound`` and
then rotates and translates it such that two other ``Compounds``
overlap. Typically, as in this case, those two other ``Compounds`` are
``Ports`` - in our case, ``methyl1.up`` and ``methyl2.up``.

