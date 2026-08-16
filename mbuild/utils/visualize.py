"""Methods for visualizing mBuild Compound and Path instances."""

import os
import tempfile
from copy import deepcopy

import networkx as nx
import numpy as np

from mbuild.compound import Compound, clone
from mbuild.path.formats import to_mol3000
from mbuild.utils.io import import_, run_from_ipython
from mbuild.utils.jsutils import overwrite_nglview_default


def visualize_path(path, radius=0.1, hide_periodic_bonds=False):
    """Visualize in 3D space using py3Dmol of the Path as a Compound.

    Parameters
    ----------
    radius : float, default 0.1
        Radius for sphere and stick representation
    hide_periodic_bonds : bool, default False
        If ``True`` bonds crossing periodic boundaries aren't shown.
    """
    py3Dmol = import_("py3Dmol")

    G = deepcopy(path.bond_graph)
    if hide_periodic_bonds:
        max_pos = np.max(path.coordinates, axis=0)
        min_pos = np.min(path.coordinates, axis=0)
        half_box_l = np.max(max_pos - min_pos) / 2
        remove_edges = []
        for n1, n2 in G.edges():
            if np.linalg.norm(path.coordinates[n1] - path.coordinates[n2]) > half_box_l:
                remove_edges.append((int(n1), int(n2)))
        print(f"Hiding {len(remove_edges)} periodic edges")
        G.remove_edges_from(remove_edges)

    view = py3Dmol.view(width=600, height=600)
    # view.addModel(mol2_string, "mol2", keepH=True)

    # Get unique bead names
    unique_names = list(dict.fromkeys(node for node in path.beads))

    # Color palette
    colors = [
        "#e6194b",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#46f0f0",
        "#f032e6",
        "#bcf60c",
        "#fabebe",
        "#008080",
        "#e6beff",
        "#9a6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000075",
        "#808080",
        "#ffffff",
        "#000000",
    ]

    data = to_mol3000(path=path, G=G)
    view = py3Dmol.view(data=data)

    for i, name in enumerate(unique_names):
        color = colors[i % len(colors)]
        indices = [int(idx) for idx, bead in enumerate(path.beads) if bead == name]
        view.addStyle(
            {"index": indices},  # Select all atoms with this bead name
            {
                "sphere": {"color": color, "radius": radius, "scale": 0.5},
                "stick": {"radius": radius / 4, "color": "grey"},
            },
        )
    view.setBackgroundColor("white")
    view.zoomTo()

    return view


def visualize_compound(
    compound,
    show_ports=False,
    backend="py3dmol",
    color_scheme={},
    bead_size=0.3,
    periodic_bond_opacity=False,
):  # pragma: no cover
    """Visualize the Compound using py3dmol (default) or nglview.

    Allows for visualization of a Compound within a Jupyter Notebook.

    Parameters
    ----------
    compound : mb.Compound
        The compound to show.
    show_ports : bool, optional, default=False
        Visualize Ports in addition to Particles
    backend : str, optional, default='py3dmol'
        Specify the backend package to visualize compounds
        Currently supported: py3dmol, nglview
    color_scheme : dict, optional
        Specify coloring for non-elemental particles
        keys are strings of the particle names
        values are strings of the colors
        i.e. {'_CGBEAD': 'blue'}
    bead_size : float, Optional, default=0.3
        Size of beads in visualization
    periodic_bond_opacity : bool, float, Optional, default=False
        Specify as a float from 0 to 1 to set the bond opacity
        for bonds that cross periodic boundaries.
    """
    viz_pkg = {
        "nglview": _visualize_nglview,
        "py3dmol": _visualize_py3dmol,
    }
    if run_from_ipython():
        if backend.lower() in viz_pkg:
            if backend.lower() == "nglview":
                return viz_pkg[backend.lower()](compound, show_ports=show_ports)
            else:
                return viz_pkg[backend.lower()](
                    compound,
                    show_ports=show_ports,
                    color_scheme=color_scheme,
                    bead_size=bead_size,
                    periodic_bond_opacity=periodic_bond_opacity,
                )
        else:
            raise RuntimeError(
                f"Unsupported visualization backend ({backend}). "
                "Currently supported backends include nglview and py3dmol"
            )

    else:
        raise RuntimeError("Visualization is only supported in Jupyter Notebooks.")


def _visualize_py3dmol(
    compound,
    show_ports=False,
    color_scheme={},
    bead_size=0.3,
    periodic_bond_opacity=False,
):
    """Visualize the Compound using py3Dmol.

    Allows for visualization of a Compound within a Jupyter Notebook.

    Parameters
    ----------
    show_ports : bool, optional, default=False
        Visualize Ports in addition to Particles
    color_scheme : dict, optional
        Specify coloring for non-elemental particles
        keys are strings of the particle names
        values are strings of the colors
        i.e. {'_CGBEAD': 'blue'}
    bead_size : float, Optional, default=0.3
        Size of beads in visualization
    periodic_bond_opacity : bool, float, Optional, default=False
        Specify as a float from 0 to 1 to set the bond opacity
        for bonds that cross periodic boundaries.

    Returns
    -------
    view : py3Dmol.view
    """
    py3Dmol = import_("py3Dmol")

    cloned = clone(compound)
    for edge in cloned.bond_graph.edges(data=True):
        if edge[2]["bond_order"] == 0.0:
            edge[2]["bond_order"] = 1.0

    modified_color_scheme = {}
    for name, color in color_scheme.items():
        # Py3dmol does some element string conversions,
        # first character is as-is, rest of the characters are lowercase
        new_name = name[0] + name[1:].lower()
        modified_color_scheme[new_name] = color
        modified_color_scheme[name] = color

    for particle in cloned.particles():
        if not particle.name:
            particle.name = "UNK"
    tmp_dir = tempfile.mkdtemp()
    # bin bonds into periodic and aperiodic bonds
    if isinstance(periodic_bond_opacity, float):
        # save into two mol2 files, one with periodic bonds and one without
        periodic_bonds, aperiodic_bonds = cloned._classify_periodic_bonds()
        periodicGraph = nx.subgraph_view(
            cloned.bond_graph,
            filter_edge=lambda n1, n2: (
                (n1, n2) in periodic_bonds or (n2, n1) in periodic_bonds
            ),
        )
        aperiodicGraph = nx.subgraph_view(
            cloned.bond_graph,
            filter_edge=lambda n1, n2: (
                (n1, n2) in aperiodic_bonds or (n2, n1) in aperiodic_bonds
            ),
        )
        cpd1 = Compound.from_bondgraph(periodicGraph)
        cpd2 = Compound.from_bondgraph(aperiodicGraph)
        cpd1.save(
            os.path.join(tmp_dir, "periodic.mol2"),
            include_ports=show_ports,
        )
        cpd2.save(
            os.path.join(tmp_dir, "aperiodic.mol2"),
            include_ports=show_ports,
        )
        view = py3Dmol.view()
        with open(os.path.join(tmp_dir, "periodic.mol2"), "r") as f:
            view.addModel(f.read(), "mol2", keepH=True)
        with open(os.path.join(tmp_dir, "aperiodic.mol2"), "r") as f:
            view.addModel(f.read(), "mol2", keepH=True)

        view.setStyle(
            {"model": 0},
            {
                "stick": {
                    "radius": bead_size * 0.3,
                    "color": "grey",
                    "opacity": periodic_bond_opacity,
                },
                "sphere": {
                    "scale": bead_size,
                    "colorscheme": modified_color_scheme,
                },
            },
        )
        view.setStyle(
            {"model": 1},
            {
                "stick": {"radius": bead_size * 0.6, "color": "grey"},
                "sphere": {
                    "scale": bead_size,
                    "colorscheme": modified_color_scheme,
                },
            },
        )
        view.zoomTo()

    else:
        cloned.save(
            os.path.join(tmp_dir, "tmp.mol2"),
            include_ports=show_ports,
            overwrite=True,
        )

        view = py3Dmol.view()
        with open(os.path.join(tmp_dir, "tmp.mol2"), "r") as f:
            view.addModel(f.read(), "mol2", keepH=True)

        view.setStyle(
            {
                "stick": {"radius": bead_size * 0.6, "color": "grey"},
                "sphere": {
                    "scale": bead_size,
                    "colorscheme": modified_color_scheme,
                },
            }
        )
        view.zoomTo()

    return view


def _visualize_nglview(compound, show_ports=False):
    """Visualize the Compound using nglview.

    Allows for visualization of a Compound within a Jupyter Notebook.

    Parameters
    ----------
    show_ports : bool, optional, default=False
        Visualize Ports in addition to Particles
    """
    nglview = import_("nglview")
    mdtraj = import_("mdtraj")  # noqa: F841
    from mdtraj.geometry.sasa import _ATOMIC_RADII

    def remove_digits(x):
        return "".join(i for i in x if not i.isdigit() or i == "_")

    for particle in compound.particles():
        particle.name = remove_digits(particle.name).upper()
        if not particle.name:
            particle.name = "UNK"
    tmp_dir = tempfile.mkdtemp()
    compound.save(
        os.path.join(tmp_dir, "tmp.mol2"),
        include_ports=show_ports,
        overwrite=True,
    )
    widget = nglview.show_file(os.path.join(tmp_dir, "tmp.mol2"))
    widget.clear()
    widget.add_ball_and_stick(cylinderOnly=True)
    elements = set([particle.name for particle in compound.particles()])
    scale = 50.0
    for element in elements:
        try:
            widget.add_ball_and_stick(
                f"_{element.upper()}",
                aspect_ratio=_ATOMIC_RADII[element.title()] ** 1.5 * scale,
            )
        except KeyError:
            ids = [
                str(i)
                for i, particle in enumerate(compound.particles())
                if particle.name == element
            ]
            widget.add_ball_and_stick(
                f"@{','.join(ids)}",
                aspect_ratio=0.17**1.5 * scale,
                color="grey",
            )
    if show_ports:
        widget.add_ball_and_stick("_VS", aspect_ratio=1.0, color="#991f00")
    overwrite_nglview_default(widget)
    return widget
