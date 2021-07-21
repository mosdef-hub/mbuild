==================
Recipe Development
==================

There may be cases where your ``Compounds`` and/or building scripts can be generalized to support a broad range of systems.
Such objects would be a valuable resource for many researchers, and might justify development of a Python package that could be distributed to the community.


``mBuild`` has been developed with this in mind, in the form of a plug-in system.
Detailed below are the specifications of this system, how to convert an existing Python project into an mBuild-discoverable plug-in, and an example.

Entry Points
------------

The basis of the plug-in system in mBuild is the `setuptools.entry_points package <https://packaging.python.org/guides/creating-and-discovering-plugins/#using-package-metadata>`_.
This allows other packages to register themselves with the ``entry_point`` group we defined in mBuild, so they are accessible through the ``mbuild.recipes`` location.
Imagine you have a class named ``my_foo`` that inherits from ``mb.Compound``.
It is currently inside of a project ``my_project`` and is accessed via a direct import, i.e. ``from my_project import my_foo``.
You can register this class as an entry point associated with ``mbuild.recipes``.
It will then be accessible from inside mBuild as a plug-in via ``mbuild.recipes.my_foo`` and a direct import will be unncessary.
The call ``import mbuild`` discovers all plug-ins that fit the ``entry_point`` group specification and makes them available under ``mbuild.recipes``.

Registering a Recipe
____________________

Here we consider the case that a user already has a Python project set up with a structure similar to the layout below.

This project can be found `here <https://github.com/justinGilmer/mbuild-fcc>`_.

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


The two important files for the user to convert their ``mBuild`` plug-in to a discoverable plug-in are ``setup.py`` and ``mbuild_fcc.py``.

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



There are two notable lines in this file that we need to focus on when developing this as a plug-in for mBuild.

The first is the import statement ``import mbuild``.
We must make sure that mbuild is installed since we are inheriting from ``mbuild.Compound``. When you decide to distribute your plug-in,
the dependencies must be listed.

The second is to select the name of the plug-in itself.
It is considered good practice to name it the name of your ``class``.
In this case, we will name the plug-in ``FCC``.

The last step is to edit the ``setup.py`` file such that the plug-in can be registered under the entry_point group ``mbuild.plugins``.

::

    from setuptools import setup

    setup(
        ...
        entry_points={ "mbuild.plugins":[ "FCC = mbuild_fcc.mbuild_fcc:FCC"]},
        ...
    )

The important section is the ``entry_points`` argument. Here we define the entry_point group we want to register with: ``"mbuild.plugins"``.
Finally, we tell Python what name to use when accessing this plug-in.
Earlier, we decided to call it ``FCC``.
This is denoted here by the name before the assignment operator ``FCC =``.
Next, we pass the location of the file with our plug-in: ``mbuild_fcc.mbuild_fcc`` as if we were located at the ``setup.py`` file.
Then, we provide the name of the class within that Python file we want to make discoverable ``:FCC``.

Since the ``setup.py`` file is located in the top folder of the python project, the first ``mbuild_fcc`` is the name of the folder, and the second is the name of the python file. The colon (``:``) is used when accessing the class that is in the python file itself.


Putting it all together
_______________________

Finally, we have ``FCC = mbuild_fcc.mbuild_fcc:FCC``.

To test this feature, you should clone the ``mbuild-fcc`` project listed above.

``git clone https://github.com/justinGilmer/mbuild-fcc``


Make sure you have mBuild installed, then run the command below after changing into the ``mbuild-fcc`` directory.

``cd mbuild-fcc``

``pip install -e .``

Note that this command will install this example from source in an editable format.


Trying it Out
_____________

To test that you set up your plug-in correctly, try importing mBuild:

``import mbuild``

If you do not receive error messages, your plug-in should be discoverable!

``help(mbuild.recipes.FCC)``
`
