import difflib

import numpy as np
import pytest

from mbuild.tests.base_test import BaseTest
from mbuild.utils.io import get_fn, import_
from mbuild.utils.validation import assert_port_exists
from mbuild.utils.jsutils import nglview_custom_tooltip


class TestUtils(BaseTest):

    def test_assert_port_exists(self, ch2):
        assert_port_exists('up', ch2)
        with pytest.raises(ValueError):
            assert_port_exists('dog', ch2)

    def test_structure_reproducibility(self, alkane_monolayer):
        filename = 'monolayer-tmp.pdb'
        alkane_monolayer.save(filename)
        with open(get_fn('monolayer.pdb')) as file1:
            with open('monolayer-tmp.pdb') as file2:
                diff = difflib.ndiff(file1.readlines(), file2.readlines())
        changes = [l for l in diff if l.startswith('+ ') or l.startswith('- ')]
        assert not changes

    def test_fn(self):
        get_fn('benzene.mol2')

        with pytest.raises((IOError, OSError)):
            get_fn('garbage_file_name.foo')

    def test_import(self):
        assert np == import_('numpy')

        with pytest.raises(ImportError):
            import_('garbagepackagename')

    def test_js_utils(self):
        nglview = import_('nglview')
        with pytest.raises(TypeError):
            nglview_custom_tooltip(object())
        test_widget = nglview.NGLWidget()
        nglview_custom_tooltip(test_widget)
        assert hasattr(test_widget, 'stage')
        assert isinstance(test_widget._ngl_msg_archive, list)
        assert len(test_widget._ngl_msg_archive) == 1
        assert isinstance(test_widget._ngl_msg_archive[0], dict)
        message_dict = test_widget._ngl_msg_archive[0]
        assert message_dict['target'] == 'Widget'
        assert message_dict['type'] == 'call_method'
        assert message_dict['methodName'] == 'executeCode'
        assert message_dict['args'] == [
                    """
                    this.stage.mouseControls.add('hoverPick', (stage, pickingProxy) => {
                        let tooltip = this.stage.tooltip;
                        if(pickingProxy && pickingProxy.atom && !pickingProxy.bond){
                            let atom = pickingProxy.atom;
                            tooltip.innerText = "ATOM: " + atom.qualifiedName() + ", Index: " + atom.index;
                        }
                    });
                 """
                ]

