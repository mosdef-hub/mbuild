/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

         Define auxillary storage for multiple torsions entries
         multors.h
*/

/* cases where there are multiple torsion entries for an atom type */

extern int *multors_cnt;                 /* count of multiple entries for each atom */
extern float ***multiples;               /* an atom's multiple torsion entries */

/* end cases where there are multiple torsion entries for an atom type */
