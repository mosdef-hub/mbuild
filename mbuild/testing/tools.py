import os
from pkg_resources import resource_filename


def get_fn(name):
    """Get the full path to one of the reference files shipped for testing

    In the source distribution, these files are in ``mbuild/testing/reference``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Args:
        name (str): Name of the file to load (with respect to the reference/ folder).

    """
    fn = resource_filename('mbuild', os.path.join('testing', 'reference', name))

    if not os.path.exists(fn):
        raise ValueError('Sorry! %s does not exists. If you just '
                         'added it, you\'ll have to re install' % fn)

    return fn