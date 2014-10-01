/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage for ring structure information
*/

extern int num_rngs;                    /* count of number of rings */
extern int *ring_num;                   /* general number of a ring */
extern int **ring_atomno;               /* atoms in ring */
extern int *arom_num;                   /* general number of aromatic ring */
extern int **arom_atomno;               /* atoms in aromatic ring */
extern int **ar_set;                    /* set for aromatic ring */
extern short *ring_type;                /* set to type of ring */

/* end ring structure information */

