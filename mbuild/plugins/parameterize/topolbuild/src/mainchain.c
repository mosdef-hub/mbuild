/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions determine the molecule's main chain     mainchain.c

   NOTICE: This is a derivative work based on study of the routines
   found in the antechamber 1.27 program parmchk, version 1.0,
   dated October, 2001 by Junmei Wang, Department of Pharmaceutical
   Chemistry, School of Pharmacy, University of California,
   San Francisco, CA  94143

   Portions of this work simplify storage allocation, and clarify
   decision trees compared to the work cited above.
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "mainchain.h"
#include "block_memory.h"
#include "mol2.h"

/* set up to find main chain atoms
   parameters are:
        num_atoms             number of atoms in molecule
*/
void set_mainchain(int num_atoms)
{
     int i, j, k, l, sel_num;

     main_num = 0;
     main_chain = (int *)allocator((num_atoms + 10), sizeof(int), FARGS);
     sel_chain = (int *)allocator((num_atoms + 10), sizeof(int), FARGS);
     sel_indx = (int *)allocator((num_atoms + 10), sizeof(int), FARGS);

     for(i = 0; i < num_atoms; i++)
          main_chain[i] = -1;

     for(i = 0; i < num_atoms; i++) {
          if(!atom_atno[i] ||
             (atom_atno[i] == 1))
                    continue;
          sel_num = 0;
          for(j = 0; j < num_atoms; j++) {
               sel_chain[j] = -1;
               sel_indx[j] = -1;
          }                             /* end for(j = 0; j < num_atoms, j++) */

          get_mainchain(num_atoms, sel_num, i);
     }                                  /* end for(i = 0; i < num_atoms; i++) */

     free_me(sel_chain, FARGS);
     free_me(sel_indx, FARGS);
     return;
}

/* recursive routine to locate the main chain atoms
   parameters are:
        num_atoms             number of atoms in molecule
        sel_num               current selection
        starter               current starting point
*/
void get_mainchain(int num_atoms, int sel_num, int starter)
{
     int i, j, k, starter1;
     int visiter;

     visiter = 0;
     sel_indx[starter] = sel_num;
     sel_chain[sel_num++] = starter;

     if(sel_num > main_num) {
          for(j = 0; j < sel_num; j++)
               main_chain[j] = sel_chain[j];
          main_num = sel_num;
     }                                  /* end if(sel_num > main_num) */

     if(sel_num > num_atoms) {
          if(main_num > num_atoms) main_num = num_atoms;
          return;
     }                                  /* end if(sel_num > num_atoms) */

     for(i = 0; i < bond_count[starter]; i++) {
          if(!atom_atno[(connect[starter][i])] ||
             (atom_atno[(connect[starter][i])] == 1))
                    continue;
          starter1 = connect[starter][i];
          if(sel_indx[starter1] != -1) continue;

          for(j = 0; j < sel_num; j++)
               if(starter1 == sel_chain[j]) {
                    visiter = 1;
                    k = j;
                    break;
               }                        /* end if(starter1 == sel_chain[j]) */

          if(visiter) {
               visiter = 0;
               sel_num = sel_indx[(sel_chain[k])] + 2;
               for(j = sel_num; j < num_atoms; j++)
                    sel_chain[j] = -1;
               continue;
          }                             /* end if(visiter) */

          get_mainchain(num_atoms, sel_num, starter1);
     }                                  /* end for(i = 0; i < bond_count[starter]; i++) */

     return;
}
