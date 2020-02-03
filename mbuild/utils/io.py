# Portions of this code are adapted from MDTraj and are released under the
# following license.

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
import inspect
import importlib
import os
from pkg_resources import resource_filename
import sys
import textwrap
import warnings
from unittest import SkipTest


class DelayImportError(ImportError, SkipTest):
    pass


MESSAGES = dict()
MESSAGES['gsd'] = '''
The code at {filename}:{line_number} requires the "gsd" package

gsd can be installed with conda using:

# conda install -c conda-forge gsd
'''

MESSAGES['nglview'] = '''
The code at {filename}:{line_number} requires the "nglview" package

nglview can be installed using:

# conda install -c conda-forge nglview

or

# pip install nglview
'''

MESSAGES['openbabel'] = '''
The code at {filename}:{line_number} requires the "openbabel" package

openbabel can be installed with conda using:

# conda install -c conda-forge openbabel

or from source following instructions at:

# http://openbabel.org/docs/current/UseTheLibrary/PythonInstall.html
'''

MESSAGES['pybel'] = MESSAGES['openbabel']

MESSAGES['foyer'] = '''
The code at {filename}:{line_number} requires the "foyer" package

foyer can be installed using:

# conda install -c mosdef foyer

or

# pip install foyer
'''

MESSAGES['protobuf'] = ''' 
The code at {filename}:{line_number} requires the "protobuf" package

protobuf can be installed using:

# conda install -c anaconda protobuf

or 

# pip install protobuf
'''

def import_(module):
    """Import a module, and issue a nice message to stderr if the module isn't installed.

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
    https://open-babel.readthedocs.io/en/latest/UseTheLibrary/migration.html#python-module
    """
    if module == 'pybel':
        try:
            return importlib.import_module('openbabel.pybel')
        except ModuleNotFoundError:
            pass
        try:
            pybel = importlib.import_module('pybel')
            msg = ('openbabel 2.0 detected and will be dropped in a future '
                   'release. Consider upgrading to 3.x.')
            warnings.warn(msg, DeprecationWarning)
            return pybel
        except ModuleNotFoundError:
            pass
    try:
        return importlib.import_module(module)
    except ImportError as e:
        try:
            message = MESSAGES[module]
        except KeyError:
            message = 'The code at {filename}:{line_number} requires the ' + module + ' package'
            e = ImportError('No module named %s' % module)

        frame, filename, line_number, function_name, lines, index = \
            inspect.getouterframes(inspect.currentframe())[1]

        m = message.format(filename=os.path.basename(filename), line_number=line_number)
        m = textwrap.dedent(m)

        bar = '\033[91m' + '#' * max(len(line) for line in m.split(os.linesep)) + '\033[0m'

        print('', file=sys.stderr)
        print(bar, file=sys.stderr)
        print(m, file=sys.stderr)
        print(bar, file=sys.stderr)
        raise DelayImportError(m)


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
    import openbabel
    has_openbabel = True
    del openbabel
except ImportError:
    has_openbabel = False

try:
    import foyer
    has_foyer = True
    del foyer
except ImportError:
    has_foyer = False

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
    has_hoomd  = False

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
    fn = resource_filename('mbuild', os.path.join('utils', 'reference', name))
    if not os.path.exists(fn):
        raise IOError('Sorry! {} does not exists.'.format(fn))
    return fn


def run_from_ipython():
    try:
        __IPYTHON__
        return True
    except NameError:
        return False
