Methane: Compounds and bonds
----------------------------

**Note**: mBuild expects all distance units to be in nanometers.

The primary building block in mBuild is a ``Compound``. Anything you
construct will inherit from this class. Let’s start with some basic
imports and initialization:

.. code:: ipython

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()

Any ``Compound`` can contain other ``Compounds`` which can be added
using its ``add()`` method. ``Compounds`` at the bottom of such a
hierarchy are referred to as ``Particles``. Note however, that this is
purely semantic in mBuild to help clearly designate the bottom of a
hierarchy.

.. code:: ipython

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')

By default a created ``Compound/Particle`` will be placed at ``0, 0, 0``
as indicated by its ``pos`` attribute. The ``Particle`` objects
contained in a ``Compound``, the bottoms of the hierarchy, can be
referenced via the ``particles`` method which returns a generator of all
``Particle`` objects contained below the ``Compound`` in the hierarchy.

**Note:** All positions in mBuild are stored in nanometers.

Any part added to a ``Compound`` can be given an optional, descriptive
string label. If the label ends with the characters ``[$]``, a list will
be created in the labels. Any subsequent parts added to the ``Compound``
with the same label prefix will be appended to the list. In the example
above, we’ve labeled the hydrogen as ``HC[$]``. So this first part, with
the label prefix ``HC``, is now referenceable via ``self['HC'][0]``. The
next part added with the label ``HC[$]`` will be referenceable via
``self['HC'][1]``.

Now let’s use these styles of referencing to connect the carbon to the
hydrogen. Note that for typical use cases, you will almost never have to
explicitly define a bond when using mBuild - this is just to show you
what’s going on under the hood:

.. code:: ipython

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')

            self.add_bond((self[0], self['HC'][0]))

As you can see, the carbon is placed in the zero index of ``self``. The
hydrogen could be referenced via ``self[1]`` but since we gave it a
fancy label, it’s also referenceable via ``self['HC'][0]``.

Alright now that we’ve got the basics, let’s finish building our
``Methane`` and take a look at it:

.. code:: ipython

    import mbuild as mb

    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Particle(name='C')
            self.add(carbon, label='C[$]')

            hydrogen = mb.Particle(name='H', pos=[0.1, 0, -0.07])
            self.add(hydrogen, label='HC[$]')

            self.add_bond((self[0], self['HC'][0]))

            self.add(mb.Particle(name='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0, 0.1, 0.07]), label='HC[$]')
            self.add(mb.Particle(name='H', pos=[0, -0.1, 0.07]), label='HC[$]')

            self.add_bond((self[0], self['HC'][1]))
            self.add_bond((self[0], self['HC'][2]))
            self.add_bond((self[0], self['HC'][3]))

.. code:: ipython

    methane = Methane()
    methane.visualize()

.. code:: ipython

    # Save to .mol2
    methane.save('methane.mol2',overwrite=True)
