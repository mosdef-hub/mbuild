
Methane: Atoms, Bonds and Compounds
-----------------------------------

The primary building block in mBuild is a ``Compound``. Anything you
construct will inherit from this class. Let's start with some basic
imports and initialization:

.. code:: python

    import mbuild as mb
    
    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()

The other parts used in building molecules are ``Atoms`` and ``Bonds``.
Either of them can be added to a ``Compound`` using its ``add()``
method. Let's add a carbon and a hydrogen atom to our ``Methane``:

.. code:: python

    import mbuild as mb
    
    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Atom(name='C')
            self.add(carbon)
    
            hydrogen = mb.Atom(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')

By default a created ``Atom`` will be placed at ``0, 0, 0`` as indicated
by its ``pos`` attribute. The ``Atom`` objects contained in a
``Compound`` can be referenced via the ``atoms`` attribute which returns
a list of all ``Atoms`` including those in any sub-compounds.

**Note:** All positions in mBuild are stored in nanometers.

Any part added to a ``Compound`` can be given an optional, descriptive
string label. If the label ends with the characters ``[$]``, a list will
be created in the labels. Any subsequent parts added to the ``Compound``
with the same label prefix will be appended to the list. In the example
above, we've labeled the hydrogen as ``HC[$]``. So this first part, with
the label prefix ``HC``, is now referenceable via ``self.HC[0]``. The
next part added with the label ``HC[$]`` will be referenceable via
``self.HC[1]``.

Now let's use these styles of referencing to connect the carbon to the
hydrogen. Note that for typical use cases, you will almost never have to
explicitly define a bond when using mBuild - this is just to show you
what's going on under the hood:

.. code:: python

    import mbuild as mb
    
    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Atom(name='C')
            self.add(carbon)
    
            hydrogen = mb.Atom(name='H', pos=[0.11, 0, 0])
            self.add(hydrogen, label='HC[$]')
    
            ch_bond = mb.Bond(self.atoms[0], self.HC[0])
            self.add(ch_bond)

As you can see, the carbon is placed in the zero index of the ``atoms``
attribute. The hydrogen could be referenced via ``self.atoms[1]`` but
since we gave it a fancy label, it's also referenceable via ``HC[0]``.

Like ``Atoms``, ``Bonds`` also have a descriptive ``kind`` attribute. By
default, ``kind`` is set to
``'{0}-{1}'.format(atom1.name, atom2.name)``.

Alright now that we've got the basics, let's finish building our
``Methane`` and take a look at it:

.. code:: python

    import mbuild as mb
    
    class Methane(mb.Compound):
        def __init__(self):
            super(Methane, self).__init__()
            carbon = mb.Atom(name='C')
            self.add(carbon)
    
            hydrogen = mb.Atom(name='H', pos=[0.1, 0, -0.07])
            self.add(hydrogen, label='HC[$]')
            
            ch_bond = mb.Bond(self.atoms[0], self.HC[0])
            self.add(ch_bond)
    
            self.add(mb.Atom(name='H', pos=[-0.1, 0, -0.07]), label='HC[$]')
            self.add(mb.Atom(name='H', pos=[0, 0.1, 0.07]), label='HC[$]')
            self.add(mb.Atom(name='H', pos=[0, -0.1, 0.07]), label='HC[$]')
    
            self.add(mb.Bond(self.atoms[0], self.HC[1]))
            self.add(mb.Bond(self.atoms[0], self.HC[2]))
            self.add(mb.Bond(self.atoms[0], self.HC[3]))

.. code:: python

    methane = Methane()
    methane.visualize()



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. parsed-literal::

    <IPython.core.display.Javascript object>


