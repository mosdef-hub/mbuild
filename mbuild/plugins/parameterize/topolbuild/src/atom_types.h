/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage and functions for atom types and atom type correlations
        atom_types.h
*/

#define MAXCHAR 256

extern char line[2*MAXCHAR];               /* the input line */
extern char mess[MAXCHAR];                 /* for error messages */

extern char *name[11];                     /* tokens pointers */

/* atom type descriptors */

extern int num_wilds;                      /* number of wild atom entries */
extern char **wild_name;                   /* wild atom names array */
extern char ***wild_elements;              /* list of wild atom element assignments by name */
extern int *wild_count;                    /* number of elements in each wild atom name */
extern int **wild_atno;                    /* atomic number of wild elements */

extern char **at_type_name;                /* atom type name */
extern int *at_type_residue;               /* atom type residue usage flag */
extern char **residue_name;                /* residue name storage if needed */
extern int *at_type_no;                    /* atom type atomic number */
extern int *at_type_attached;              /* atom type attached atoms number */
extern int *at_type_attached_H;            /* atom type attached hydrogens number */
extern int *at_type_ewd;                   /* for H, no. ewd attached to atom attached to H */
extern char **at_type_prop;                /* atom type properties */
extern char **at_type_env;                 /* atom type chemical environments */
extern char **at_env_bonds;                /* terminator to atom type line */

/* end atom type descriptors */

/* atom correlation data */

extern char **corr_name;                   /* name of the atom to correlate */
extern int *index_improper;                /* correlation impropers index */
extern char ***corr_to;                    /* atom names correlated to */
extern int *corr_num;                      /* number of correlations here */

/* end atom correlation data */

extern void read_at_types(int *at_type_cnt, int *wild_cnt, int max_type,
                          int max_wild, int just_count, char *filename);
extern void read_at_corr(int *corr_cnt, int max_corr, int just_count,
                         char *filename);
extern void release_types(int at_type_cnt, int wild_cnt);
extern void release_corr(int corr_count);
extern int str_balanced(char cin, char cout, char *str);
