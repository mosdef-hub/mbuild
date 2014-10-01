/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage for molecule from mol2 file and for the ancillary derived
        atom types, forcefield, and measurements data
*/

#define MAXNAME 128

/* information directly from .mol2 file */
                    /* header information */

extern char mol_name[MAXNAME];     /* molecule name record */
extern char model_type[MAXNAME];   /* mol2 Molecule model type record */
extern char charge_type[MAXNAME];  /* mol2 Molecule charge method */
extern float adj_chg;              /* charge adjustment factor */
extern int use_charge;             /* charge usage flag */

                    /* end header information */

                    /* atoms information */

extern char **atom_name;           /* list of atom names in mol2 atoms record */
extern char **atom_type;           /* list of atom types in mol2 atoms record */
extern char **atom_residue;        /* residue in which atom appears in mol2 atoms record */
extern int *atom_segno;            /* residue or chain id in mol2 atoms record */
extern float *atom_charge;         /* atom charge from mol2 atoms record */

                    /* end atoms information */

                    /* coordinates */

extern float *atom_x;              /* atom x coordinate from mol2 atoms record */
extern float *atom_y;              /* atom y coordinate from mol2 atoms record */
extern float *atom_z;              /* atom z coordinate from mol2 atoms record */

                    /* end coordinates */

                    /* bond information */

extern char **bond_type;           /* type of bond from mol2 bonds record */
extern int *bond_i;                /* first atom in bond from mol2 bonds record */
extern int *bond_j;                /* second atom in bond from mol2 bonds record */
extern int *bond_count;            /* number of atoms bound to this atom from mol2 bonds */
extern int **connect;              /* table of atom connections based on mol2 bonds record */
extern int **bonds_list;           /* table of which bonds each atom participates in */

                    /* end bond information */

                    /* restraints information */

extern int rstr_dist;              /* number of simple distance restraints */
extern int rstr_range;             /* number of distance range restraints */
extern int rstr_angle;             /* number of angle restraints */
extern int rstr_torsn;             /* number of torsion restraints */
extern int **restr;                /* type and atoms [i][0] is type, & [i][1 to 4] are atom no. */
extern float **frest;              /* restraint values and force constants */

                    /* end restraints information */
/* end information directly from .mol2 file */

/* derived molecule information */

extern int *atom_sat;
extern int *atom_ewd;
extern float *atom_mass;           /* atom mass from table search */
extern int *atom_atno;             /* atomic number from table search */
extern int *how_bonded;            /* translation of bond type */
extern float *meas_len;            /* measured bond length based on coordinates */
extern char **atom_symbl;          /* element symbol in proper case */
extern char **a_ff_type;           /* generalized amber force field type */
extern int **angle_tab;            /* table of atom numbers in an angle */
extern float *meas_ang;            /* measured angles based on coordinates */
extern int **torsion_tab;          /* table of atom numbers in a torsion or improper torsion */
extern float *meas_tors;           /* measured torsions based on coordinates */

                     /* improper locus information */

extern int improper_num;           /* number marked improper by simple search */
extern int *improper;              /* flag for marked improper atoms */
extern int **improperid;           /* atoms connected to each marked improper */
extern float *meas_impro;          /* measured improper torsion based on coordinates */

                    /* end improper locus information */

                    /* residue atom ordering information */

extern int *atom_order;            /* ordering of the atoms for output.  Solves
                                      problem of all residue atoms required to be
                                      together in gromacs topologies */
extern int *back_order;            /* reverse of ordering of the atoms for output. */
extern int true_atoms;             /* true count of atoms */

                    /* end residue atom ordering information */
/* end derived molecule information */

/* information from force field files */
               /* read from selected amber force field */

extern float *a_pol;               /* atom polarizability from selected amber */

extern float *a_bond_force;        /* bond force from selecter amber */
extern float *a_bond_length;       /* bond length from selecter amber */
extern float *a_bond_beta;         /* for morse and cubic bonds */
extern short *a_bond_type;         /* type of bond force constant */

extern float *a_angle_force;       /* angle force from selected amber */
extern float *a_angle_angle;       /* angle length from selected amber */
extern float **a_angle_quartic;    /* dimensioned 5 by # of atoms to allow for quartics */
extern int *a_meth;                /* angle assignment method */
extern short *a_angle_type;        /* type of angle force constant */
extern char *opposite;             /* opposite or adjacent flag */

extern int tors_imp;               /* count torsion that should be improper */
extern int regs_wanted;            /* number of wanted torsions of type 0 or 1 */
extern int regs_all;               /* total number of torsions of type 0 or 1 */
extern int rbs_wanted;             /* number of wanted RB torsions */
extern int rbs_all;                /* total number of RB torsions */
extern int fourier_wanted;         /* number of wanted fourier coefficient torsions */
extern int fourier_all;            /* total number of fourier coefficient torsions */
extern int table_wanted;           /* number of wanted table entry torsions */
extern int table_all;              /* total number of table entry torsions */

extern int *a_tors_mult;           /* torsion multiplier factor from selected amber */
extern float *a_tors_phase;        /* torsion phase from selected amber */
extern float *a_tors_force;        /* torsion force from selected amber */
extern float *a_tors_term;         /* the torsion term from selected amber */
extern float **a_tors_set;         /* dimensioned 6 by # of atoms to allow for RB's */
extern short *a_tors_type;         /* type of dihedral */

extern int *a_impr_mult;           /* improper multiplier factor from selected amber */
extern float *a_impr_force;        /* improper force from selected amber */
extern float *a_impr_phase;        /* improper phase from selected amber */
extern float *a_impr_term;         /* the improper term from selected amber */

extern float *a_vdw_atom_radius;   /* van der Waals radius from selected amber */
extern float *a_vdw_atom_pot;      /* van der Waals potential from selected amber */

/*         Define constraint storage elements */

extern int *constrain_i;           /* first element of a constraint */
extern int *constrain_j;           /* second element of a constraint */
extern float *dist_constr;         /* constraint distance */
extern short *constr_type;         /* type of constraint */
extern int no_constraints;         /* count of constraints in molecule */

/* end constraint storage elements */

               /* end read from selected amber force field */
/* end information from force field files */
