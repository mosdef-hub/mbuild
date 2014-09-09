/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        determine the molecule's main chain
        determine impropers
        and rename based on mainchain    header    mainchain.h
*/

extern void set_mainchain(int num_atoms);
extern void get_mainchain(int num_atoms, int sel_num, int starter);
extern void renameit(int num_atoms, int num_rings);
extern void loc_improper(int atom_num, int do_ua);

/* main chain derivation information */

extern int main_num;                    /* actual number of main chain atoms */
extern int *main_chain;                 /* flag for main chain atoms */

/* end main chain derivation information */

extern int *sel_chain;                  /* local chain selection flags */
extern int *sel_indx;

extern char **new_names;
