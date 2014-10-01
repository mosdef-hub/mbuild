/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
       Defining declarations for external variables
       Created solely for increased portability based on
       a suggestion in Harbison & Steele (1984) C: A
       Reference Manual, Prentice Hall, Inc., Englewood
       Cliffs, NJ    ISBN 0-13-110008-4
            initializations.c
*/

#include <stdio.h>

#define MAXNAME 128
#define MAXCHAR 256

                    /* main program options storage */

FILE *stdlog = NULL;                /* define log file pointer */

char line[(2 * MAXCHAR)];           /* the input line */
char mess[MAXCHAR];                 /* for error messages */
char *in_command;                   /* command line */

char *name[11];                     /* tokens pointers */

char *use_define[6];                /* topology defines pointers */

int lstart_at = 0;
int res_start_at = 0;
int rstart_at = 0;
int lstart2_at = 0;
int use_amber = 0;
int center_it = 0;
short charge_set = 0;
short no_oppose = 0;

char ambername[MAXCHAR];
char mycorr_name[MAXCHAR];

                    /* end main program options storage */

/*      Define storage for molecule from mol2 file and for the ancillary derived
        atom types, forcefield, and measurements data
*/

/* information directly from .mol2 file */
                    /* header information */

char mol_name[MAXNAME];             /* molecule name record */
char model_type[MAXNAME];           /* mol2 Molecule model type record */
char charge_type[MAXNAME];          /* mol2 Molecule charge method */
float adj_chg = 1.0;                /* charge adjustment factor */
int use_charge = 0;                 /* charge usage flag */

                    /* end header information */

                    /* atoms information */

char **atom_name = NULL;            /* list of atom names in mol2 atoms record */
char **atom_type = NULL;            /* list of atom types in mol2 atoms record */
char **atom_residue = NULL;         /* residue in which atom appears in mol2 atoms record */
int *atom_segno = NULL;             /* residue or chain id in mol2 atoms record */
float *atom_charge = NULL;          /* atom charge from mol2 atoms record */

                    /* end atoms information */

                    /* coordinates */

float *atom_x = NULL;               /* atom x coordinate from mol2 atoms record */
float *atom_y = NULL;               /* atom y coordinate from mol2 atoms record */
float *atom_z = NULL;               /* atom z coordinate from mol2 atoms record */

                    /* end coordinates */

                    /* bond information */

char **bond_type = NULL;            /* type of bond from mol2 bonds record */
int *bond_i = NULL;                 /* first atom in bond from mol2 bonds record */
int *bond_j = NULL;                 /* second atom in bond from mol2 bonds record */
int *bond_count = NULL;             /* number of atoms bound to this atom from mol2 bonds */
int **connect = NULL;               /* table of atom connections based on mol2 bonds record */
int **bonds_list = NULL;            /* table of which bonds each atom participates in */

                    /* end bond information */

                    /* restraints information */

int rstr_dist = 0;                  /* number of simple distance restraints */
int rstr_range = 0;                 /* number of distance range restraints */
int rstr_angle = 0;                 /* number of angle restraints */
int rstr_torsn = 0;                 /* number of torsion restraints */
int **restr = NULL;                 /* type and atoms [i][0] is type, [i][1 to 4] are atom no. */
float **frest = NULL;               /* restraint values and force constants */

                    /* end restraints information */
/* end information directly from .mol2 file */

/* derived molecule information */

int *atom_sat = NULL;
int *atom_ewd = NULL;
float *atom_mass = NULL;            /* atom mass from table search */
int *atom_atno = NULL;              /* atomic number from table search */
int *how_bonded = NULL;             /* translation of bond type */
float *meas_len = NULL;             /* measured bond length based on coordinates */
char **atom_symbl = NULL;           /* element symbol in proper case */
char **a_ff_type = NULL;            /* generalized amber force field type */
int **angle_tab = NULL;             /* table of atom numbers in an angle */
float *meas_ang = NULL;             /* measured angles based on coordinates */
int **torsion_tab = NULL;           /* table of atom numbers in a torsion or improper torsion */
float *meas_tors = NULL;            /* measured torsions based on coordinates */

                     /* improper locus information */

int improper_num = 0;               /* number marked improper by simple search */
int *improper = NULL;               /* flag for marked improper atoms */
int **improperid = NULL;            /* atoms connected to each marked improper */
float *meas_impro = NULL;           /* measured improper torsion based on coordinates */

                    /* end improper locus information */

                    /* residue atom ordering information */

int *atom_order = NULL;             /* ordering of the atoms for output.  Solves
                                       problem of all residue atoms required to be
                                       together in gromacs topologies */
int *back_order = NULL;             /* reverse of ordering of the atoms for output. */
int true_atoms = 0;                 /* true count of atoms */

                    /* end residue atom ordering information */
/* end derived molecule information */

/* information from force field files */
               /* read from selected amber force field */

float *a_pol = NULL;                /* atom polarizability from selected amber */

float *a_bond_force = NULL;         /* bond force from selected amber */
float *a_bond_length = NULL;        /* bond length from selected amber */
float *a_bond_beta = NULL;          /* for morse and cubic bonds */
short *a_bond_type = NULL;          /* type of bond force constant */

float *a_angle_force = NULL;        /* angle force from selected amber */
float *a_angle_angle = NULL;        /* angle length from selected amber */
float **a_angle_quartic = NULL;     /* dimensioned 5 by # of atoms to allow for quartics */
int *a_meth = NULL;                 /* angle assignment method */
short *a_angle_type = NULL;         /* type of angle force constant */
char *opposite = NULL;              /* opposite or adjacent flag */

int tors_imp = 0;                   /* count torsion that should be improper */
int *a_tors_mult = NULL;            /* torsion multiplier factor from selected amber */
float *a_tors_phase = NULL;         /* torsion phase from selected amber */
float *a_tors_force = NULL;         /* torsion force from selected amber */
float *a_tors_term = NULL;          /* the torsion term from selected amber */
float **a_tors_set = NULL;          /* dimensioned 6 by # of atoms to allow for RB's */
short *a_tors_type = NULL;          /* type of dihedral */
short *wanted_tors = NULL;          /* For Gromacs, do I want this torsion? */

int *a_impr_mult = NULL;            /* improper multiplier factor from selected amber */
float *a_impr_force = NULL;         /* improper force from selected amber */
float *a_impr_phase = NULL;         /* improper phase from selected amber */
float *a_impr_term = NULL;          /* the improper term from selected amber */

float *a_vdw_atom_radius = NULL;    /* van der Waals radius from selected amber */
float *a_vdw_atom_pot = NULL;       /* van der Waals potential from selected amber */

               /* end read from selected amber force field */
/* end information from force field files */

/*****************************************************************************/

/*         Define storage for amber type forcefields */

/* amber type force field parameters */
          /* the atom type parameters */

char **ffatom_name = NULL;               /* name of an Amber atom type */
float *ffatom_mass = NULL;               /* mass of an Amber atom type (may be unified atom
                                            model mass and not necessarily element mass) */
float *ffatom_pol = NULL;                /* a polarizability factor */

          /* end the atom type parameters */

          /* the bond parameters */

char **bond_name1 = NULL;               /* name of an Amber atom type */
char **bond_name2 = NULL;               /* name of an Amber atom type */
float *bond_force = NULL;               /* Amber FF bond force */
float *bond_length = NULL;              /* Amber FF bond length */
float *bond_beta = NULL;                /* second force constant for some types of bonds */
short *FF_bond_type = NULL;             /* bond force constant type */
                                        /* defined as: */
                                        /*     1     bond kJ/(mol nm^2) */
                                        /*     2     G96  kJ/(mol nm^4) */
                                        /*     3     morse D kJ/mol and beta 1/nm */
                                        /*     4     C2 kJ/(mol nm^2) and C3 kJ/(mol nm^3) */
                                        /* Gromacs type 5, connection, not used */
                                        /* Gromacs type 6, harmonic potential for restraints */
                                        /*     7     FENE kJ/(mol nm^2) */

          /* end the bond parameters */

          /* the angle parameters */

char **angle_name1 = NULL;              /* name of an Amber atom type */
char **angle_name2 = NULL;              /* name of an Amber atom type */
char **angle_name3 = NULL;              /* name of an Amber atom type */
float *angle_force = NULL;              /* Amber FF angle force */
float *angle_angle = NULL;              /* Amber FF angle itself */
float **angle_quartic = NULL;           /* quartic angle force constants */
short *FF_angle_type = NULL;            /* angle force constant type */
                                        /* defined as: */
                                        /*     1     angle kJ/(mol rad^2) */
                                        /*     2     G96   kJ/mol */
                                        /*     6     quartic C0 kJ/mol, C1 kJ/(mol rad) */
                                        /*           C2 kJ/(mol rad^2), C3 kJ/(mol rad^3) */
                                        /*           C4 kJ/(mol rad^4) */

          /* end the angle parameters */

          /* the dihedral angle (torsion) parameters */

char **tors_name1 = NULL;               /* name of an Amber atom type */
char **tors_name2 = NULL;               /* name of an Amber atom type */
char **tors_name3 = NULL;               /* name of an Amber atom type */
char **tors_name4 = NULL;               /* name of an Amber atom type */
int *tors_mult = NULL;                  /* Amber FF force multiplier factor */
float *tors_force = NULL;               /* Amber FF torsion force */
float *tors_phase = NULL;               /* Amber FF torsion angle */
float *tors_term = NULL;                /* an Amber FF torsion term */
float **tors_RB = NULL;                 /* allow for RB constants */
int *tors_more = NULL;                  /* flag more torsion entries for atom type */
short *FF_tors_type = NULL;             /* type of torsion term */
                                        /* defined as: */
                                        /*     1     proper kJ/mol */
                                        /*     3     RB 6 constants all kJ/mol */

          /* end the dihedral angle (torsion) parameters */

          /* the inproper dihedral angle (torsion) parameters */

char **impro_name1 = NULL;              /* name of an Amber atom type */
char **impro_name2 = NULL;              /* name of an Amber atom type */
char **impro_name3 = NULL;              /* name of an Amber atom type */
char **impro_name4 = NULL;              /* name of an Amber atom type */
int *impro_mult = NULL;                 /* Amber FF force multiplier factor (set to 1 in code) */
int *impro_xcount = NULL;               /* Amber FF count of X's in improper */
float *impro_force = NULL;              /* Amber FF improper torsion force */
float *impro_phase = NULL;              /* Amber FF improper torsion angle */
float *impro_term = NULL;               /* an Amber FF improper torsion term */
                                        /* all improper force constants are in units of kJ/(mol rad^2) */

          /* end the inproper dihedral angle (torsion) parameters */

          /* the van der Waals parameters */

char **vdw_atom_name = NULL;            /* name of an Amber atom type */
float *vdw_atom_radius = NULL;          /* Amber FF van der Waals radius (may be unified atom
                                           model radius and not necessarily element radius) */
float *vdw_atom_pot = NULL;             /* Amber FF van der Waals potential */

          /* end the van der Waals parameters */

          /* the equivalent van der Waals lists */

char ***equiv_name = NULL;             /* list of names of Amber atom types with the first name
                                          also appearing in the van der Waals parameters  */
int *equiv_count = NULL;               /* number of equivalent Amber atom types given */

          /* end the equivalent van der Waals lists */

/* end amber type force field parameters */

/*****************************************************************************/

/*         Define storage and functions for atom types and atom type correlations */

/* atom type descriptors */

int num_wilds = 0;                     /* number of wild atom entries */
char **wild_name = NULL;               /* wild atom names array */
char ***wild_elements = NULL;          /* list of wild atom element assignments by name */
int *wild_count = NULL;                /* number of elements in each wild atom name */
int **wild_atno = NULL;                /* atomic number of wild elements */

char **at_type_name = NULL;            /* atom type name */
int *at_type_residue = NULL;           /* atom type residue usage flag */
char **residue_name = NULL;            /* residue name storage if needed */
int *at_type_no = NULL;                /* atom type atomic number */
int *at_type_attached = NULL;          /* atom type attached atoms number */
int *at_type_attached_H = NULL;        /* atom type attached hydrogens number */
int *at_type_ewd = NULL;               /* for H, no. ewd attached to atom attached to H */
char **at_type_prop = NULL;            /* atom type properties */
char **at_type_env = NULL;             /* atom type chemical environments */
char **at_env_bonds = NULL;            /* terminator to atom type line */

/* end atom type descriptors */

/* atom correlation data */

char **corr_name = NULL;               /* name of the atom to correlate */
int *index_improper = NULL;            /* correlation impropers index */
char ***corr_to = NULL;                /* atom names correlated to */
int *corr_num = NULL;                  /* number of correlations here */

/* end atom correlation data */

/*****************************************************************************/

/*      Define storage for judging atom types on the basis
        of bonding and of chemical environment
*/

int scan_num = 0;
int initial = 0;

/* bond information arrays */

int *sb = NULL;
int *SB = NULL;
int *db = NULL;
int *DB = NULL;
int *tb = NULL;
int *TB = NULL;
int *AB = NULL;
int *DL = NULL;

/* end bond information arrays */

/* chemical environment arrays */

char **apch = NULL;
char **env_bonds = NULL;
char **env_atom1 = NULL;
char **env_atom2 = NULL;
char **atom_chem = NULL;
char **envblockstr = NULL;
char ***env_atom_name = NULL;
char ***env_ap = NULL;
char ***env_strname = NULL;

int *env_len = NULL;
int *chem_indx = NULL;
int **env_con = NULL;
int **env_index = NULL;

int *selchain = NULL;
int *selindx = NULL;
int *sel_chain_indx = NULL;

int *schain_id = NULL;
int *schain_num = NULL;
int **schain_atom = NULL;

int env_bond_num = 0;

/* end chemical environment arrays */

/*****************************************************************************/

/* main chain derivation information */

int main_num = 0;                /* actual number of main chain atoms */
int *main_chain = NULL;          /* flag for main chain atoms */

/* end main chain derivation information */

int *sel_chain = NULL;           /* local chain selection flags */
int *sel_indx = NULL;

char **new_names = NULL;

/* end main chain derivation information */

/*****************************************************************************/

/*        Define storage for ring structure information */

int num_rngs = 0;                       /* count of number of rings */
int *ring_num = NULL;                   /* general number of a ring */
int **ring_atomno = NULL;               /* atoms in ring */
int *arom_num = NULL;                   /* general number of aromatic ring */
int **arom_atomno = NULL;               /* atoms in aromatic ring */
int **ar_set = NULL;                    /* set for aromatic ring */
short *ring_type = NULL;                /* set to type of ring */

/* end ring structure information */


/* temporary ring structure storage */

int *t_ring_num = NULL;
int **t_ring_atomno = NULL;
short *t_ring_type = NULL;

/*       end storage for ring structure information */

/*****************************************************************************/

/*         Define storage for atomic correlation data */

/* correlated names arrays */
          /* for gaff to gaff */

char **gaff_name_corr1 = NULL;           /* first correlation name */
char **gaff_name_corr2 = NULL;           /* second correlation name */
char **gaff_name_corr3 = NULL;           /* third correlation name */
int *gaff_corr_num = NULL;               /* number of correlated names */

          /* end for gaff to gaff */
/* end correlated names arrays */

/*****************************************************************************/

/*         Define auxillary storage for multiple torsions entries */

/* cases where there are multiple torsion entries for a bond type */

int *multors_cnt = NULL;                 /* count of multiple entries for each atom */
float ***multiples = NULL;               /* an atom's multiple torsion entries */

/* end cases where there are multiple torsion entries for a bond type */

/*****************************************************************************/

/*         Define additional storage elements for Tripos type forcefields */

/* Tripos type force field additional parameters */

char **trp_bnd_type = NULL;
char **tors_bond = NULL;

/* end Tripos type force field additional parameters */

/*****************************************************************************/

/*         Define special storage elements for Gromacs type force fields */

/* Gromacs type force field additional parameters */

int num_macro = 0;                       /* number of macro names loaded */
int *macro_nmbr = NULL;                  /* number of atom types in each macro name */
char **macro_name = NULL;                /* macro names */
char ***macro_types = NULL;              /* atom types in each macro name */
char ***gmx_tors_resid = NULL;           /* list of residues to which a torsion applies */
int *gmx_tors_num = NULL;                /* number of residues listed for this torsion */
int constrain = 0;                       /* number of constraint types in force field */
char **constrain_name1 = NULL;           /* first constraint element name */
char **constrain_name2 = NULL;           /* second constraint element name */
float *constrain_dist = NULL;            /* constraint distance */
short *type_constr = NULL;               /* type of constraint */

/* end Gromacs type force field additional parameters */

/*****************************************************************************/

/*         Define constraint storage elements */

int *constrain_i = NULL;                 /* first element of a constraint */
int *constrain_j = NULL;                 /* second element of a constraint */
float *dist_constr = NULL;               /* constraint distance */
short *constr_type = NULL;               /* type of constraint */
int no_constraints = 0;                  /* counter of constraints in molecule */

/* end constraint storage elements */

/*****************************************************************************/

/*         Define counters for torsion types */

int regs_wanted = 0;                     /* number of wanted torsions of type 0 or 1 */
int regs_all = 0;                        /* total number of torsions of type 0 or 1 */
int rbs_wanted = 0;                      /* number of wanted RB torsions */
int rbs_all = 0;                         /* total number of RB torsions */
int fourier_wanted = 0;                  /* number of wanted fourier coefficient torsions */
int fourier_all = 0;                     /* total number of fourier coefficient torsions */
int table_wanted = 0;                    /* number of wanted table entry torsions */
int table_all = 0;                       /* total number of table entry torsions */

/* end counters for torsion types */

/*****************************************************************************/

/*       end defining declarations for external variables */
