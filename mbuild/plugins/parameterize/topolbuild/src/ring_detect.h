/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        detect ring structures header   ring_detect.h
*/

#define MAXCHAR 256

extern void cycle_detect(int *rings, int selected, int st_num, int *arom,
                         int *arom_ind, int max_rng);
extern void cycle_clean(int *ringers);
extern void cycle_props(int rings);
extern void aromaticity(int no_atms, int no_bnds, int no_rings);
extern void detect_ring(int atom_num, int bond_num, int max_rng,
                        int *num_rings);

extern char mess[MAXCHAR];

/* temporary ring structure storage */

extern int *t_ring_num;
extern int **t_ring_atomno;
extern short *t_ring_type;
