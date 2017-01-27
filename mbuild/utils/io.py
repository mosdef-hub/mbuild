import os
from pkg_resources import resource_filename

try:
    import intermol
    has_intermol = True
except ImportError:
    has_intermol = False

try:
    import foyer
    has_foyer = True
except ImportError:
    has_foyer = False

try:
    import gsd
    has_gsd = True
except ImportError:
    has_gsd = False


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




