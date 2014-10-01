/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        judge the bond type header   judge_bond.h
*/

#include "block_memory.h"
#include "mol2.h"
#include "rings.h"

#define MAXCHAR 256

extern void judge_bond(int num_bonds, int num_atoms, int num_rings);
extern void is_carbon(int i, int index[]);
extern void is_nitrogen(int i, int index[]);
extern void is_P_or_S(int i, int index[], int k);
extern void is_S_or_O(int i, int index[], int k);
