/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage for atomic correlation data and
        set the similar names lists for an assigned type
        based on atom correlation tables header
        set_similars.h
*/

extern void similars_alloc(int num_atoms);
extern void release_similars(int num_atoms);
extern void set_similars(int num_atoms, int corr_cnt);

/* correlated names arrays */
          /* for gaff to gaff */

extern char **gaff_name_corr1;           /* first correlation name */
extern char **gaff_name_corr2;           /* second correlation name */
extern char **gaff_name_corr3;           /* third correlation name */
extern int *gaff_corr_num;               /* number of correlated names */

          /* end for gaff to gaff */
/* end correlated names arrays */
