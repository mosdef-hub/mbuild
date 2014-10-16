Methane: Atoms, Bonds and Compounds
-----------------------------------
In this example, we'll cover the core classes that make up an mBuild structure
by building a methane molecule.

The primary building block in mBuild is a `Compound`. Anything you construct
will inherit from this class. Let's start with some basic imports and
initialization::

    $ from mbuild.compound import Compound
    $
    $
    $ class Methane(Compound):
    $    def __init__(self):
    $        super(Methane, self).__init__()

.. note:: The use of the `super()` method is required here to resolve
          `Compound`'s `multiple inheritance <http://stackoverflow.com/questions/3277367/how-does-pythons-super-work-with-multiple-inheritance>`_
          from the `MBase`, `PartMixin` and `HasPartsMixin` classes.

The other classes used in building molecules are `Atoms` and `Bonds`. Either of
them can be added to a `Compound` using its `add()` method. Let's add a carbon
and a hydrogen atom to our `Methane`::

    $ from mbuild.compound import Compound
    $ from mbuild.atom import Atom
    $
    $
    $ class Methane(Compound):
    $    def __init__(self):
    $        super(Methane, self).__init__()
    $        carbon = Atom(kind='C')
    $        self.add(carbon)
    $
    $        hydrogen = Atom(kind='H', pos=[0.15, 0, 0])
    $        self.add(hydrogen, label='hc[$]')

By default a created `Atom` will be placed at `0, 0, 0` as indicated by its
`pos` attribute. All positions in mBuild are stored in nanometers.

Explain labeling

Now let's connect the carbon to the hydrogen::

    $ from mbuild.compound import Compound
    $ from mbuild.atom import Atom
    $ from mbuild.bond import Bond
    $
    $
    $ class Methane(Compound):
    $    def __init__(self):
    $        super(Methane, self).__init__()
    $        carbon = Atom(kind='C')
    $        self.add(carbon)
    $
    $        hydrogen = Atom(kind='H', pos=[0.15, 0, 0])
    $        self.add(hydrogen, label='hc[$]')
    $
    $        ch_bond = Bond(self.atom[0], self.hc[0])
    $        self.add(ch_bond)

Explainn labeling here


Like `Atoms`, `Bonds` also have a descriptive `kind` attribute. By default,
`kind` is set to `'{0}-{1}'.format(atom1.kind, atom2.kind)`.

Alright now that we've got the basics, let's finish building our `Methane` and
take a look at it!::


    $ from mbuild.compound import Compound
    $ from mbuild.atom import Atom
    $ from mbuild.bond import Bond
    $
    $
    $ class Methane(Compound):
    $     def __init__(self):
    $         super(Methane, self).__init__()
    $         carbon = Atom(kind='C')
    $         self.add(carbon)
    $
    $        hydrogen = Atom(kind='H', pos=[0.15, 0, 0])
    $        self.add(hydrogen, label='hc[$]')
    $
    $        ch_bond = Bond(self.atoms[0], self.hc[0])
    $        self.add(ch_bond)
    $
    $         self.add(Atom(kind='H', pos=[0, 0.15, 0]), 'hc[$]')
    $         self.add(Bond(self.atoms[0], self.hc[1]))
    $         self.add(Atom(kind='H', pos=[-0.15, 0, 0]), 'hc[$]')
    $         self.add(Bond(self.atoms[0], self.hc[2]))
    $         self.add(Atom(kind='H', pos=[0, -0.15, 0]), 'hc[$]')
    $         self.add(Bond(self.atoms[0], self.hc[3]))
    $
    $ if __name__ == "__main__":
    $     methane = Methane()
    $     methane.visualize()

.. note:: The `visualize()` method currently invokes a very primative call to
          VMD from the command-line. If it fails for you but you do have VMD
          installed, the method works by writing an intermediate output file
          named `visualize_Methane.pdb` which you can load yourself. We are
          currently working on creating a more robust and useful interface VMD
          but any help would be appreciated (see issue #32).

