/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

        Read and use gromacs force field data header gromacs_use.h
*/

#include "utilities.h"
#include "block_memory.h"
#include "use_amber.h"
#include "compare_FF.h"
#include "rings.h"
#include "similars.h"
#include "mol2.h"
#include "amber_FF.h"
#include "GMX_FF.h"
#include "mainchain.h"
#include "multors.h"
#include "read_amber.h"
#include "param_tables.h"

extern short charge_set;
extern short no_oppose;

extern void do_gmx_FFs(char *ante_name, char *lev1_name, int num_atoms,
                       int num_bonds, char *amberff, int angle_max,
                       int tors_max, int ua_model, char *residue,
                       int by_corr);
extern void united_atoms(int num_atoms, char *resid, int nmol);
extern void read_ffgmx1(char *ffname, int *nmmol, int *nmbond, int *nmangle,
                        int *nmtors, int *nmimpro, int *nmvdw);
extern void set_gromacs_torsFF(int num_tors, int tors_cntFF);
extern int comp_gmx_tors(int tors_atom1, int tors_atom2, int tors_atom3,
                         int tors_atom4, int tors_cntFF);
extern void special_tors(int numb, FILE *ffptr, int i, char *myfiles);
extern void set_constraints(int num_cons);
extern void set_charge(int num_atoms, int num_mol);
