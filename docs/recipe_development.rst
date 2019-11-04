==============
Recipe Development
==============

mBuild, owing to its object-oriented-design, encourages users to subclass ``Compound`` to disseminate reproducible scripts for building their systems of interest.
Sharing a script to reproduce the initial configuration and force field parameters (if applicable) of a system is a great step towards increasing the ability of others to reproduce another's work.

However, there might be some cases where those subclasses and scripts can be generalized to support even more systems than initially realized.
These subclasses would be a valuable resource for many and might justify development of a new Python package to distribute this work to the community.


mBuild has been developed with this in mind, in the form of a plug-in recipe system.
Detailed below are the specifications of this system, how to convert your own Python project into an mBuild-discoverable plug-in, and an example.

Entry Points
--------

The basis of the recipe system in mBuild relies on `setuptools.entry_points <https://packaging.python.org/guides/creating-and-discovering-plugins/#using-package-metadata>`_.
This allows other packages to register themselves with the ``entry_point`` group we defined in mBuild, so they are accessible through the ``mbuild.recipes.my_recipe_foo`` location.
In this case, ``my_recipe_foo`` would be the recipe you want to share.
You can install many mBuild recipes, and only have to import mBuild to have access to them all!
It is not necessary to use calls like this: ``from my_foo_project import my_recipe_foo``, the call ``import mbuild`` discovers all recipes that fit the ``entry_point`` group specification and will then make them available under ``mbuild.recipes``. 

Registering a Recipe
________


In the case that a user already has a Python project set up with a structure similar to the layout below:

::

    mbuild_fcc
    ├── LICENSE
    ├── README.md
    ├── mbuild_fcc
    │   ├── mbuild_fcc.py
    │   └── tests
    │       ├── __init__.py
    │       └── test_fcc.py
    └── setup.py


The two important files for the user to convert their ``mBuild`` recipe to a discoverable plug-in are ``setup.py`` and ``mbuild_fcc.py``.

To begin, lets first inspect the ``mbuild_fcc.py`` file, a shortened snippet is below.

::

    import mbuild


    class FCC(mbuild.Compound):
        """Create a mBuild Compound with a repeating unit of the FCC unit cell.

        ... (shortened for viewability)

        """

        def __init__(self, lattice_spacing=None, compound_to_add=None, x=1, y=1, z=1):
            super(FCC, self).__init__()

            # ... (shortened for viewability)

    if __name__ == "__main__":
        au_fcc_lattice = FCC(lattice_spacing=0.40782,
                             compound_to_add=mbuild.Compound(name="Au"),
                             x=5, y=5, z=1)
        print(au_fcc_lattice)



There are two very important lines in this file that we need to know when developing a plug-in for mBuild.

The first is the import statement ``import mbuild``.
We must make sure that mbuild is installed since we are inheriting from ``mbuild.Compound``. When you decide to distribute your plug-in,
the dependencies must be listed.

The second is the name of the plug-in itself. It is in general good practice to name it the name of your ``class``.
In this case, the name of the plug-in would be ``FCC``.

The last step is to edit the ``setup.py`` file such that the plug-in can be registered under the entry_point group ``mbuild.plugins``.

::
    from setuptools import setup

    setup(
        name="mbuild_fcc",
        install_requires="mbuild",
        entry_points={
            "mbuild.plugins":[
                "FCC = mbuild_fcc.mbuild_fcc:FCC"
            ]
        },
        py_modules=["mbuild_fcc"],
    )

This is once again, a very minimal setup file, a more thoroughly tested and developed package will have more information contained within.

The important section is the ``entry_points`` argument. Here we define the entry_point group we want to plug in to ``"mbuild.plugins"``.
Finally, we tell python what name we want to have our plug-in be discoverable by, in this case we call it ``FCC`` as denoted by the name before the assignment operator ``FCC =``.
Next, we pass the "path" from this ``setup.py`` file to the python file that contains the class we want to register as a plug-in: ``mbuild_fcc.mbuild_fcc``.
Then, we provide the name of the class within that python file we want to make discoverable ``:FCC``.


Putting it all together
________

Finally, we have ``FCC = mbuild_fcc.mbuild_fcc:FCC``.

Make sure you have mBuild installed, then install your plug-in project in the same location as the ``setup.py`` file.

``pip install -e .``

Note that this command will install this project from source in an editable format.


Trying it Out
_________

To test that you set up your plug-in correctly, try importing mBuild!

``import mbuild``

If you received no error messages, your recipe should be discoverable!

``help(mbuild.recipes.FCC)``

