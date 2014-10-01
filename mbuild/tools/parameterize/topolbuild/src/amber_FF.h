/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage for amber type forcefields
*/

/* amber type force field parameters */
          /* the atom type parameters */

extern char **ffatom_name;               /* name of an Amber atom type */
extern float *ffatom_mass;               /* mass of an Amber atom type (may be unified atom
                                            model mass and not necessarily element mass) */
extern float *ffatom_pol;                /* a polarizability factor */

          /* end the atom type parameters */

          /* the bond parameters */

extern char **bond_name1;               /* name of an Amber atom type */
extern char **bond_name2;               /* name of an Amber atom type */
extern float *bond_force;               /* Amber FF bond force */
extern float *bond_length;              /* Amber FF bond length */
extern float *bond_beta;                /* second force constant for some types of bonds */
extern short *FF_bond_type;             /* bond force constant type */
                                        /* defined as: */
                                        /*     1     bond kJ/(mol nm^2) */
                                        /*     2     G96  kJ/(mol nm^4) */
                                        /*     3     morse D kJ/mol and beta 1/nm */
                                        /*     4     C2 kJ/(mol nm^2) and C3 kJ/(mol nm^3) */
                                        /* Gromacs type 5, connection, not used */
                                        /* Gromacs type 6, harmonic potential reserved for restraints */
                                        /*     7     FENE kJ/(mol nm^2) */

          /* end the bond parameters */

          /* the angle parameters */

extern char **angle_name1;              /* name of an Amber atom type */
extern char **angle_name2;              /* name of an Amber atom type */
extern char **angle_name3;              /* name of an Amber atom type */
extern float *angle_force;              /* Amber FF angle force */
extern float *angle_angle;              /* Amber FF angle itself */
extern float **angle_quartic;           /* quartic angle force constants */
extern short *FF_angle_type;            /* angle force constant type */
                                        /* defined as: */
                                        /*     1     angle kJ/(mol rad^2) */
                                        /*     2     G96   kJ/mol */
                                        /*     6     quartic C0 kJ/mol, C1 kJ/(mol rad) */
                                        /*           C2 kJ/(mol rad^2), C3 kJ/(mol rad^3) */
                                        /*           C4 kJ/(mol rad^4) */

          /* end the angle parameters */

          /* the dihedral angle (torsion) parameters */

extern char **tors_name1;               /* name of an Amber atom type */
extern char **tors_name2;               /* name of an Amber atom type */
extern char **tors_name3;               /* name of an Amber atom type */
extern char **tors_name4;               /* name of an Amber atom type */
extern int *tors_mult;                  /* Amber FF force multiplier factor */
extern float *tors_force;               /* Amber FF torsion force */
extern float *tors_phase;               /* Amber FF torsion angle */
extern float *tors_term;                /* an Amber FF torsion term */
extern float **tors_RB;                 /* allow for RB constants */
extern int *tors_more;                  /* flag more torsion entries for atom type */
extern short *FF_tors_type;             /* type of torsion term */
                                        /* defined as: */
                                        /*     1     proper kJ/mol */
                                        /*     3     RB 6 constants all kJ/mol */

          /* end the dihedral angle (torsion) parameters */

          /* the inproper dihedral angle (torsion) parameters */

extern char **impro_name1;              /* name of an Amber atom type */
extern char **impro_name2;              /* name of an Amber atom type */
extern char **impro_name3;              /* name of an Amber atom type */
extern char **impro_name4;              /* name of an Amber atom type */
extern int *impro_mult;                 /* Amber FF force multiplier factor (set to 1 in code) */
extern int *impro_xcount;               /* Amber FF count of X's in improper */
extern float *impro_force;              /* Amber FF improper torsion force */
extern float *impro_phase;              /* Amber FF improper torsion angle */
extern float *impro_term;               /* an Amber FF improper torsion term */

          /* end the inproper dihedral angle (torsion) parameters */

          /* the van der Waals parameters */

extern char **vdw_atom_name;            /* name of an Amber atom type */
extern float *vdw_atom_radius;          /* Amber FF van der Waals radius (may be unified atom
                                           model radius and not necessarily element radius) */
extern float *vdw_atom_pot;             /* Amber FF van der Waals potential */

          /* end the van der Waals parameters */

          /* the equivalent van der Waals lists */

extern char ***equiv_name;             /* list of names of Amber atom types with the first name also
                                          appearing in the van der Waals parameters  */
extern int *equiv_count;               /* number of equivalent Amber atom types given */

          /* end the equivalent van der Waals lists */

/* end amber type force field parameters */
