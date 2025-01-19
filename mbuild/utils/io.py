# ruff: noqa: F401
"""Module for working with external libraries.

Portions of this code are adapted from MDTraj and are released under the
following license.

##############################################################################
# MDTraj is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 2.1
# of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with MDTraj. If not, see <http://www.gnu.org/licenses/>.
##############################################################################
"""

import importlib
import inspect
import os
import sys
import textwrap
import warnings
from unittest import SkipTest

import importlib_resources as resources


class DelayImportError(ImportError, SkipTest):
    """Error to allow better import handling."""

    pass


MESSAGES = dict()
MESSAGES["gsd"] = """
The code at {filename}:{line_number} requires the "gsd" package

gsd can be installed with conda using:

# conda install -c conda-forge gsd
"""

MESSAGES["nglview"] = """
The code at {filename}:{line_number} requires the "nglview" package

nglview can be installed using:

# conda install -c conda-forge nglview

or

# pip install nglview
"""

MESSAGES["py3Dmol"] = """
The code at {filename}:{line_number} requires the "py3Dmol" package

py3Dmol can be installed using:

# conda install -c conda-forge py3Dmol

or

# pip install py3Dmol
"""

MESSAGES["rdkit"] = """
The code at {filename}:{line_number} requires the "rdkit" package

rdkit can be installed with conda using:

# conda install -c conda-forge rdkit

or from source following instructions at:

https://www.rdkit.org/docs/Install.html#installation-from-source
"""

MESSAGES["openbabel"] = """
The code at {filename}:{line_number} requires the "openbabel" package

openbabel can be installed with conda using:

# conda install -c conda-forge openbabel

or from source following instructions at:

# http://openbabel.org/docs/current/UseTheLibrary/PythonInstall.html
"""

MESSAGES["pybel"] = MESSAGES["openbabel"]

MESSAGES["mdtraj"] = """
The code at {filename}:{line_number} requires the "mdtraj" package
mdtraj can be installed using:
# conda install -c conda-forge mdtraj
or
# pip install mdtraj
"""

MESSAGES["foyer"] = """
The code at {filename}:{line_number} requires the "foyer" package

foyer can be installed using:

# conda install -c conda-forge foyer

or

# pip install foyer
"""

MESSAGES["garnett"] = """
The code at {filename}:{line_number} requires the "garnett" package

garnett can be installed with conda using:

# conda install -c conda-forge garnett
"""

MESSAGES["pycifrw"] = """
The code at {filename}:{line_number} requires the "pycifrw" package

pycifrw can be installed with conda using:

# conda install -c conda-forge pycifrw
"""

MESSAGES["freud"] = """
The code at {filename}:{line_number} requires the "freud" package

freud can be installed with conda using:

# conda install -c conda-forge freud
"""


def import_(module):
    """Import a module and issue a nice message to stderr if it isn't installed.

    Parameters
    ----------
    module : str
        The module you'd like to import, as a string

    Returns
    -------
    module : {module, object}
        The module object

    Examples
    --------
    >>> # the following two lines are equivalent. the difference is that the
    >>> # second will check for an ImportError and print you a very nice
    >>> # user-facing message about what's wrong (where you can install the
    >>> # module from, etc) if the import fails
    >>> import tables
    >>> tables = import_('tables')

    Notes
    -----
    The pybel/openbabel block is meant to resolve compatibility between
    openbabel 2.x and 3.0.  There may be other breaking changes but the change
    in importing them is the major one we are aware of. For details, see
    https://open-babel.readthedocs.io/en/latest/UseTheLibrary/migration.
    html#python-module
    """
    if module == "pybel":
        try:
            return importlib.import_module("openbabel.pybel")
        except ModuleNotFoundError:
            pass
        try:
            pybel = importlib.import_module("pybel")
            msg = (
                "openbabel 2.0 detected and will be dropped in a future "
                "release. Consider upgrading to 3.x."
            )
            warnings.warn(msg, DeprecationWarning)
            return pybel
        except ModuleNotFoundError:
            pass
    if module == "openbabel":
        try:
            return importlib.import_module("openbabel.openbabel")
        except ModuleNotFoundError:
            pass
        try:
            openbabel = importlib.import_module("openbabel")
            msg = (
                "openbabel 2.0 detected and will be dropped in a future "
                "release. Consider upgrading to 3.x."
            )
            warnings.warn(msg, DeprecationWarning)
            return openbabel
        except ModuleNotFoundError:
            pass
    try:
        return importlib.import_module(module)
    except ImportError:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = (
                "The code at {filename}:{line_number} requires the "
                f"{module} package"
            )
            raise ImportError(f"No module named {module}")

        (
            frame,
            filename,
            line_number,
            function_name,
            lines,
            index,
        ) = inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(filename=os.path.basename(filename), line_number=line_number)
        m = textwrap.dedent(m)

        bar = (
            "\033[91m"
            + "#" * max(len(line) for line in m.split(os.linesep))
            + "\033[0m"
        )

        print("", file=sys.stderr)
        print(bar, file=sys.stderr)
        print(m, file=sys.stderr)
        print(bar, file=sys.stderr)
        raise DelayImportError(m)


try:
    import hoomd

    has_hoomd = True
    del hoomd
except ImportError:
    has_hoomd = False

try:
    import intermol

    has_intermol = True
    del intermol
except ImportError:
    has_intermol = False

try:
    import gsd

    has_gsd = True
    del gsd
except ImportError:
    has_gsd = False

try:
    from openbabel import openbabel

    has_openbabel = True
    del openbabel
except ImportError:
    has_openbabel = False

try:
    import mdtraj

    has_mdtraj = True
    del mdtraj
except ImportError:
    has_mdtraj = False

try:
    import foyer

    has_foyer = True
    del foyer
except ImportError:
    has_foyer = False

try:
    import gmso

    has_gmso = True
except ImportError:
    has_gmso = False

try:
    import networkx

    has_networkx = True
    del networkx
except ImportError:
    has_networkx = False

try:
    import hoomd

    has_hoomd = True
    del hoomd
except ImportError:
    has_hoomd = False

try:
    import py3Dmol

    has_py3Dmol = True
    del py3Dmol
except ImportError:
    has_py3Dmol = False

try:
    import garnett

    has_garnett = True
    del garnett
except ImportError:
    has_garnett = False

try:
    import CifFile

    has_pycifrw = True
    del CifFile
except ImportError:
    has_pycifrw = False

try:
    import rdkit

    has_rdkit = True
    del rdkit
except ImportError:
    has_rdkit = False

try:
    import freud

    has_freud = True
    del freud
except ImportError:
    has_freud = False


def get_fn(name):
    """Get the full path to one of the reference files shipped for utils.

    In the source distribution, these files are in ``mbuild/utils/reference``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    name : str
        Name of the file to load (with respect to the reference/ folder).
    """
    fn = resources.files("mbuild").joinpath("utils", "reference", name)
    # fn = resource_filename("mbuild", os.path.join("utils", "reference", name))
    if not os.path.exists(fn):
        raise IOError(f"Sorry! {fn} does not exists.")
    return str(fn)


def run_from_ipython():
    """Get whether python is being run interactively."""
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
