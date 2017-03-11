"""Execute notebooks as tests, reporting an error if cell throws an exception.

Adapted from MDTraj https://github.com/mdtraj/mdtraj/pull/1130 which adapted
from https://gist.github.com/minrk/2620876.
"""

import glob
import os
import sys

from jupyter_client import KernelManager
import nbformat
import pytest


EXAMPLE_NOTEBOOKS = [f for f in glob.glob('mbuild/examples/*/*.ipynb')]


@pytest.mark.skipif(sys.platform in ['win32'],
                    reason="Not testing examples on Appveyor")
@pytest.mark.parametrize("filepath", EXAMPLE_NOTEBOOKS)
def test_examples(filepath):
    check_one_notebook(filepath)


def check_one_notebook(filepath):
    folder, filename = os.path.split(filepath)
    os.chdir(folder)
    with open(filename) as f:
        nb = nbformat.reads(f.read(), nbformat.NO_CONVERT)
    run_notebook(nb)
    os.chdir('../../../')


def run_notebook(nb):
    km = KernelManager()
    km.start_kernel(stderr=open(os.devnull, 'w'))
    kc = km.client()
    kc.start_channels()
    shell = kc.shell_channel
    # simple ping:
    kc.execute("pass")
    shell.get_msg()

    failures = 0
    for cell in nb.cells:
        if cell.cell_type != 'code':
            continue
        kc.execute(cell.source)
        # wait for finish, maximum 20s
        reply = shell.get_msg(timeout=20)['content']
        if reply['status'] == 'error':
            failures += 1
            print("\nFAILURE:")
            print('\n'.join(reply['traceback']))
            print()

    kc.stop_channels()
    km.shutdown_kernel()
    del km
    if failures > 0:
        raise Exception()
