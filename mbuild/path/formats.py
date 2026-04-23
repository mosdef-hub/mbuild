def to_mol2(path):
    """Convert NetworkX graph with xyz attribute to MOL2 format including bonds"""
    G = path.bond_graph

    mol2_lines = []

    # Get unique names and create mapping
    unique_names = list(set(G.nodes[node]["name"] for node in G.nodes()))
    node_to_atom = {
        node: atom_num for atom_num, node in enumerate(G.nodes(), start=1)
    }

    n_atoms = len(G.nodes())
    n_bonds = len(G.edges())

    # @<TRIPOS>MOLECULE section
    mol2_lines.append("@<TRIPOS>MOLECULE")
    mol2_lines.append("NETWORKX_GRAPH")
    mol2_lines.append(f"{n_atoms} {n_bonds} 0 0 0")
    mol2_lines.append("SMALL")
    mol2_lines.append("NO_CHARGES")
    mol2_lines.append("")

    # @<TRIPOS>ATOM section
    mol2_lines.append("@<TRIPOS>ATOM")
    for atom_num, node in enumerate(G.nodes(), start=1):
        x, y, z = path.coordinates[node]
        name = path.beads[node]
        subst_id = unique_names.index(name) + 1
        atom_type = "CG.A"  # Generic CG atom type
        charge = 0.0

        mol2_lines.append(
            f"{atom_num:7d} {name:4s} {x:10.4f} {y:10.4f} {z:10.4f} "
            f"{atom_type:5s} {subst_id:5d} {name:4s} {charge:10.4f}"
        )

    mol2_lines.append("")

    # @<TRIPOS>BOND section
    mol2_lines.append("@<TRIPOS>BOND")
    for bond_num, (node1, node2) in enumerate(G.edges(), start=1):
        atom1 = node_to_atom[node1]
        atom2 = node_to_atom[node2]
        bond_type = "1"  # Single bond

        mol2_lines.append(f"{bond_num:6d} {atom1:5d} {atom2:5d} {bond_type:>4s}")

    return "\n".join(mol2_lines)

def to_mol(path):
    """
    Convert mBuild Path to SDF/MOL format

    Parameters:
    -----------
    path : mbuild.Path
    mol_name : str
        Name of the molecule
    atom_types : dict
        Mapping of atom indices to atom type symbols (e.g., {0: 'A', 5: 'B'})
        If None, all atoms will be 'A'
    """
    lines = []

    # Header
    lines.append("PATH GRAPH\n")
    lines.append("     RDKit          3D\n")
    lines.append("\n")

    # Counts line: natoms nbonds nlist 3D chiral stext nrxn nreac nproduct v2000
    n_atoms = len(path.coordinates)
    n_bonds = len(path.bond_graph.edges())
    lines.append(f"{n_atoms:4d} {n_bonds:4d}  0  0  0  0  0  0  0  0999 V3000\n")

    # Atom block
    for i, (coord, bead_name) in enumerate(zip(path.coordinates, path.beads)):
        lines.append(
            f"{coord[0]:10.4f}{coord[1]:10.4f}{coord[2]:10.4f} {bead_name.strip('_'):<3s} 0  0  0  0  0  0  0  0  0  0  0  0\n"
        )

    # Bond block
    for edge in path.bond_graph.edges():
        # atom1, atom2, bond_type (1=single), stereo
        lines.append(f"{edge[0] + 1:4d} {edge[1] + 1:4d}  1  0\n")

    # End
    lines.append("M  END\n")
    lines.append("$$\n")

    return "".join(lines)

def to_mol3000(path, G=None):
    """
    Convert mBuild Path to SDF/MOL V3000 format

    Parameters:
    -----------
    G : nx.Graph, default None
        Bondgraph to use for visualization.
    """
    if G is None:
        G = path.bond_graph
    lines = []

    # Header block (3 lines)
    lines.append("PATH GRAPH\n")
    lines.append("     RDKit          3D\n")
    lines.append("\n")

    # Counts line for V3000
    n_atoms = len(path.coordinates)
    n_bonds = len(G.edges())
    lines.append("  0  0  0     0  0            999 V3000\n")

    # Begin CTAB
    lines.append("M  V30 BEGIN CTAB\n")
    lines.append(f"M  V30 COUNTS {n_atoms} {n_bonds} 0 0 0\n")

    # Atom block
    lines.append("M  V30 BEGIN ATOM\n")
    for i, (coord, bead_name) in enumerate(
        zip(path.coordinates, path.beads), start=1
    ):
        atom_type = bead_name.strip("_")
        lines.append(
            f"M  V30 {i} {atom_type} {coord[0]:.4f} {coord[1]:.4f} {coord[2]:.4f} 0\n"
        )
    lines.append("M  V30 END ATOM\n")

    # Bond block
    lines.append("M  V30 BEGIN BOND\n")
    for bond_idx, edge in enumerate(G.edges(), start=1):
        # V3000: bond_index bond_type atom1 atom2
        lines.append(f"M  V30 {bond_idx} 1 {edge[0] + 1} {edge[1] + 1}\n")
    lines.append("M  V30 END BOND\n")

    # End CTAB
    lines.append("M  V30 END CTAB\n")

    # End of record
    lines.append("M  END\n")

    return "".join(lines)
