/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to set the similar names lists for an assigned type
        based on atom correlation tables.  set_similars.c

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
#include "block_memory.h"
#include "similars.h"
#include "mol2.h"
#include "atom_types.h"

/* allocate the similar names storage
   parameter is:
        num_atoms              number of atoms in molecule
*/
void similars_alloc(int num_atoms)
{
     gaff_name_corr1 = (char **)mat_alloc(num_atoms, 4, FARGS);
     gaff_name_corr2 = (char **)mat_alloc(num_atoms, 4, FARGS);
     gaff_name_corr3 = (char **)mat_alloc(num_atoms, 4, FARGS);
     gaff_corr_num = (int *)allocator(num_atoms, sizeof(int), FARGS);

     return;
}

/* free the similar names storage
   parameter is:
        num_atoms              number of atoms in molecule
*/
void release_similars(int num_atoms)
{
     int i;

     for(i = 0; i < num_atoms; i++) {
          free_me(gaff_name_corr1[i], FARGS);
          free_me(gaff_name_corr2[i], FARGS);
          free_me(gaff_name_corr3[i], FARGS);
     }
     free_me(gaff_name_corr1, FARGS);
     free_me(gaff_name_corr2, FARGS);
     free_me(gaff_name_corr3, FARGS);
     free_me(gaff_corr_num, FARGS);

     return;
}

/* set the similar names in the molecule based on atomic correlation data
   parameter is:
        num_atoms              number of atoms in molecule
        corr_cnt               count of atomic correlation entries
*/
void set_similars(int num_atoms, int corr_cnt)
{
     int i, j, ibrk;

     for(i = 0; i < num_atoms; i++) {
          ibrk = 0;
          for(j = 0; j < corr_cnt; j++) {
               if(strcmp(a_ff_type[i], corr_name[j])) continue;
               gaff_corr_num[i] = corr_num[j];
               improper[i] = index_improper[j];
               ibrk = 1;
               if(!corr_num[j]) break;
               strcpy(gaff_name_corr1[i], &corr_to[j][0][0]);
               if(corr_num[j] == 1) break;
               strcpy(gaff_name_corr2[i], &corr_to[j][1][0]);
               if(corr_num[j] == 2) break;
               strcpy(gaff_name_corr3[i], &corr_to[j][2][0]);
               break;
          }                                       /* end for(j = 0; j < corr_cnt; j++) */
          if(!ibrk) gaff_corr_num[i] = -1;        /* flag failures and move on */
     }                                            /* end for(i = 0; i < num_atoms; i++) */

     return;
}
