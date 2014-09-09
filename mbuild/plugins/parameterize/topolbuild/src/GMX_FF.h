/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

        Define further storage for Gromacs type forcefields
        GMX_FF.h
*/

/* Gromacs type force field parameters */

extern int num_macro;                       /* number of macro names loaded */
extern int *macro_nmbr;                     /* number of atom types in each macro name */
extern char **macro_name;                   /* macro names */
extern char ***macro_types;                 /* atom types in each macro name */

/* torsion residues to handle special cases in gromacs force fields
   where the same chemical formula gets different dihedral settings
   in different conformaitons (mostly alpha vs. beta sugars)
*/
extern char ***gmx_tors_resid;           /* list of residues to which a torsion applies */
extern int *gmx_tors_num;                /* number of residues listed for this torsion */
extern short *wanted_tors;               /* For Gromacs, do I want this torsion? */

/* handle constraints in gromacs force fields */
extern int constrain;                    /* number of constraint types in force field */
extern char **constrain_name1;           /* first constraint element name */
extern char **constrain_name2;           /* second constraint element name */
extern float *constrain_dist;            /* constraint distance */
extern short * type_constr;              /* type of constraint */

/* end Gromacs type force field additional parameters */
