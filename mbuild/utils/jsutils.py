"""
These are the set of utility methods which are used to communicate with underlying 'js'
libraries by the various notebook visualization libraries used by mbuild.
"""
from .io import import_


def nglview_custom_tooltip(widget):
    """Change visualization tooltip by adding a mouseControl signal

    This method takes in a nglview.NGLWidget and changes the default hover
    behaviour of the widget to add the atom index when it is hovered over
    the atom.
    Parameters:
    ----------
    widget: nglview.NGLWidget, the ipython widget view.
    Returns:
    --------
    None
    Raises:
    ------
    TypeError: If widget is not of type nglview.NGLWidget
    """
    nglview = import_('nglview')
    if not isinstance(widget, nglview.NGLWidget):
        raise TypeError("The argument widget can only be of type nglview.NGLWidget not {}".format(type(widget)))
    tooltip_js = """
                    this.stage.mouseControls.add('hoverPick', (stage, pickingProxy) => {
                        let tooltip = this.stage.tooltip;
                        if(pickingProxy && pickingProxy.atom && !pickingProxy.bond){
                            let atom = pickingProxy.atom;
                            tooltip.innerText = "ATOM: " + atom.qualifiedName() + ", Index: " + atom.index;
                        }
                    });
                 """
    widget._js(tooltip_js)
