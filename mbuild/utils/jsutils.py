"""Utilities for communicating with Javascript (js) libraries.

These are the set of utility methods which are used to communicate with
underlying 'js' libraries by the various notebook visualization libraries used
by mBuild.
"""

from .io import import_


def overwrite_nglview_default(widget):
    """Change the default visualization in nglview.

    This method takes in a nglview.NGLWidget and changes the default hover
    behaviour of the widget to add the atom index when it is hovered over the
    atom. It also overwrites the click signal from the stage to include extra
    information(atom index) in the text display, whenever an atom or bond is
    clicked.

    Parameters
    ----------
    widget: nglview.NGLWidget,
        the ipython widget view.

    Raises
    ------
    TypeError: If widget is not of type nglview.NGLWidget
    """
    nglview = import_("nglview")
    if not isinstance(widget, nglview.NGLWidget):
        raise TypeError(
            "The argument widget can only be of type nglview.NGLWidget not "
            f"{type(widget)}."
        )
    tooltip_js = """
                    this.stage.mouseControls.add('hoverPick', (stage, pickingProxy) => {
                        let tooltip = this.stage.tooltip;
                        if(pickingProxy && pickingProxy.atom && !pickingProxy.bond){
                            let atom = pickingProxy.atom;
                            tooltip.innerText = "ATOM: " + atom.qualifiedName() + ", Index: " + atom.index;
                        }
                    });
                 """

    infotext_js = """
                    this.stage.signals.clicked.removeAll();
                    this.stage.signals.clicked.add((pickingProxy) => {
                            if(pickingProxy){
                               let pickingText = null;
                               this.model.set('picked', {});
                               this.touch();
                               let currentPick = {};
                               if(pickingProxy.atom){
                                    currentPick.atom1 = pickingProxy.atom.toObject();
                                    currentPick.atom1.name = pickingProxy.atom.qualifiedName();
                                    pickingText = "Atom: " + currentPick.atom1.name + ", Index: "
                                                  + pickingProxy.atom.index;
                               }
                               else if(pickingProxy.bond){
                                    currentPick.bond = pickingProxy.bond.toObject();
                                    currentPick.atom1 = pickingProxy.bond.atom1.toObject();
                                    currentPick.atom1.name = pickingProxy.bond.atom1.qualifiedName();
                                    currentPick.atom2 = pickingProxy.bond.atom2.toObject();
                                    currentPick.atom2.name = pickingProxy.bond.atom2.qualifiedName();
                                    pickingText = "Bond: " + currentPick.atom1.name +
                                                    `(${pickingProxy.bond.atom1.index})` +
                                                    " - " + currentPick.atom2.name    +
                                                    `(${pickingProxy.bond.atom2.index})`;
                               }

                               if(pickingProxy.instance){
                                    currentPick.instance = pickingProxy.instance;
                               }
                               var nComponents = this.stage.compList.length;
                               for(let i = 0; i < nComponents; i++){
                                    let comp = this.stage.compList[i];
                                    if(comp.uuid == pickingProxy.component.uuid){
                                        currentPick.component = i;
                                    }
                               }
                               this.model.set('picked', currentPick);
                               this.touch();
                               this.$pickingInfo.text(pickingText);
                            }
                    });
                """
    widget._js(tooltip_js)
    widget._js(infotext_js)
