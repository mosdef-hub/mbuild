"""CHARMM Par format."""

import warnings

__all__ = ["write_par"]


def write_par(structure, filename):
    """Write CHARMM Par file given a parametrized structure.

    Notes
    -----
    Follows format according to
    https://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-unix-html/
    node25.html
    Furthermore, ParmEd should support writing CHARMM par, rtf, str files
    by converting the parmed.Structure into parmed.CharmmParameterSet

    Parmed stores rmin/2 in "rmin"
    """
    # ATOMS
    with open(filename, "w") as f:
        f.write("ATOMS\n")
        unique_atoms = set()
        for atom in structure.atoms:
            unique_atoms.add((atom.atom_type.name, atom.atom_type.mass))
        for atom in unique_atoms:
            f.write("MASS -1 {:8s} {:8.4f}\n".format(atom[0], atom[1]))

        f.write("\nBONDS\n")
        unique_bonds = set()
        for bond in structure.bonds:
            unique_bonds.add(
                (
                    bond.atom1.atom_type.name,
                    bond.atom2.atom_type.name,
                    bond.type,
                )
            )

        for bond in unique_bonds:
            f.write(
                "{:8s} {:8s} {:.5f} {:.5f}\n".format(
                    bond[0], bond[1], bond[2].k, bond[2].req
                )
            )

        f.write("\nANGLES\n")
        unique_angles = set()
        unique_ubs = set()
        for angle in structure.angles:
            associated_ub = False
            for ub in structure.urey_bradleys:
                if ((angle.atom1, angle.atom3) == (ub.atom1, ub.atom2)) or (
                    angle.atom3,
                    angle.atom1,
                ) == (ub.atom1, ub.atom2):
                    unique_ubs.add(
                        (
                            angle.atom1.atom_type.name,
                            angle.atom2.atom_type.name,
                            angle.atom3.atom_type.name,
                            angle.type,
                            ub.type,
                        )
                    )
                    associated_ub = True

            if not associated_ub:
                unique_angles.add(
                    (
                        angle.atom1.atom_type.name,
                        angle.atom2.atom_type.name,
                        angle.atom3.atom_type.name,
                        angle.type,
                    )
                )

        for ub in unique_ubs:
            f.write(
                "{:8s} {:8s} {:8s} {:.5f} {:.5f} {:.5f} {:.5f}\n".format(
                    ub[0],
                    ub[1],
                    ub[2],
                    ub[3].k,
                    ub[3].theteq,
                    ub[4].k,
                    ub[4].req,
                )
            )
        for angle in unique_angles:
            f.write(
                "{:8s} {:8s} {:8s} {:.5f} {:.5f}\n".format(
                    angle[0], angle[1], angle[2], angle[3].k, angle[3].theteq
                )
            )

        # These dihedrals need to be PeriodicTorsion Style (Charmm style)
        if len(structure.rb_torsions) > 0:
            warnings.warn("RB Torsions detected, but unsupported in par writer")
        f.write("\nDIHEDRALS\n")
        unique_dihedrals = set()
        scnb = set()
        for dihedral in structure.dihedrals:
            if not dihedral.improper:
                unique_dihedrals.add(
                    (
                        dihedral.atom1.atom_type.name,
                        dihedral.atom2.atom_type.name,
                        dihedral.atom3.atom_type.name,
                        dihedral.atom4.atom_type.name,
                        dihedral.type,
                    )
                )
                scnb.add(dihedral.type.scnb)
            else:
                msg = (
                    "AMBER-style improper detected between "
                    + "{} {} {} {}".format(
                        dihedral.atom1,
                        dihedral.atom2,
                        dihedral.atom3,
                        dihedral.atom4,
                    )
                    + ", but unsupported in par writer"
                )
                warnings.warn(msg)

        for dihedral in unique_dihedrals:
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n".format(
                    dihedral[0],
                    dihedral[1],
                    dihedral[2],
                    dihedral[3],
                    dihedral[4].phi_k,
                    dihedral[4].per,
                    dihedral[4].phase,
                )
            )

        f.write("\nIMPROPER\n")
        unique_impropers = set()
        for improper in structure.impropers:
            unique_impropers.add(
                (
                    improper.atom1.atom_type.name,
                    improper.atom2.atom_type.name,
                    improper.atom3.atom_type.name,
                    improper.atom4.atom_type.name,
                    improper.type,
                )
            )
        for improper in unique_impropers:
            f.write(
                "{:8s} {:8s} {:8s} {:8s} {:.5f} {:5d} {:.5f}\n".format(
                    improper[2],
                    improper[0],
                    improper[1],
                    improper[3],
                    improper[4].psi_k,
                    0,
                    improper[4].psi_eq,
                )
            )

        sc_nb = [a for a in scnb]
        if len(sc_nb) > 1:
            warnings.warn(
                "Multiple 1-4 LJ scalings were detected, "
                "defaulting to first LJ scaling detected, {}".format(sc_nb[0])
            )
            sc_nb = sc_nb[0]
        elif len(sc_nb) == 1:
            sc_nb = sc_nb[0]
        elif len(sc_nb) == 0:
            warnings.warn("No 1-4 LJ scaling was detected, defaulting 1")
            sc_nb = 1.0

        f.write("\nNONBONDED\n")
        unique_atypes = set()
        for atom in structure.atoms:
            unique_atypes.add(atom.atom_type)
        for atype in unique_atypes:
            # atype, 0.0, epsilon, rmin/2, 0.0, epsilon(1-4), rmin/2 (1-4)
            f.write(
                "{:8s} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f} {:8.3f}\n".format(
                    atype.name,
                    0.0,
                    -1 * atype.epsilon,
                    atype.rmin,
                    0.0,
                    -1 * sc_nb * atype.epsilon,
                    atype.rmin,
                )
            )

        if structure.has_NBFIX():
            warnings.warn("NBFixes detected but unsupported in par writer")

        f.write("\nEND")
