"""Methods for visualizing mBuild Compound and Path instances."""

from copy import deepcopy

from mbuild.utils.io import import_
import numpy as np


def visualize_path(path, radius=0.1, hide_periodic_bonds=False):
    """Visualize in 3D space using py3Dmol of the Path as a Compound.

    Parameters
    ----------
    radius : float, default 0.06
        Radius for sphere and stick representation
    """
    py3Dmol = import_("py3Dmol")

    G = deepcopy(path.bond_graph)
    if hide_periodic_bonds:
        max_pos = np.max(path.coordinates, axis=0)
        min_pos = np.min(path.coordinates, axis=0)
        half_box_l = np.max(max_pos - min_pos) / 2
        remove_edges = []
        for n1, n2 in G.edges():
            if (
                np.linalg.norm(path.coordinates[n1] - path.coordinates[n2])
                > half_box_l
            ):
                remove_edges.append((int(n1), int(n2)))
        print(f"Hiding {len(remove_edges)} periodic edges")
        G.remove_edges_from(remove_edges)

    view = py3Dmol.view(width=600, height=600)
    # view.addModel(mol2_string, "mol2", keepH=True)

    # Get unique bead names
    unique_names = list(dict.fromkeys(G.nodes[node]["name"] for node in G.nodes()))

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

    data = path.to_mol3000(G)
    view = py3Dmol.view(data=data)
    for i, name in enumerate(unique_names):
        color = colors[i % len(colors)]
        view.addStyle(
            {"elem": name.strip("_")},  # Select all atoms with this name
            {
                "sphere": {"color": color, "radius": radius, "scale": 0.5},
                "stick": {"radius": radius / 4, "color": "grey"},
            },
        )
    view.setBackgroundColor("white")
    view.zoomTo()
    scale_factor = max(1, 5 - int(np.log10(len(path.coordinates))))
    view.zoom(scale_factor)  # helps zoom on smaller systems

    return view
