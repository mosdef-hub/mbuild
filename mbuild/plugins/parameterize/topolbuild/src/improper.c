/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Function to locate impropers         improper.c

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
#include "mainchain.h"
#include "mol2.h"
#include "rings.h"

#define MAXCHAR 256
#define MAXNAME 128

extern char mess[MAXCHAR];

/* locate impropers
   parameters are:
        atom_num          number of atoms
        do_ua             united atoms model flag
*/
void loc_improper(int atom_num, int do_ua)
{
     int i, j, indx, tmpint;

     improper_num = 0;

     for(i = 0; i < atom_num; i++) {
          improper[i] = 0;
          if((atom_atno[i] == 6) && (bond_count[i] == 3))
               improper[i] = 1;
/* special needs for cases of tetrahedral C with single H in united atoms models */
          if(do_ua)
               if((atom_atno[i] == 6) && (bond_count[i] == 4)) {
                    indx = 0;
                    for(j = 0; j < 4; j++) {
                         if(connect[i][j] < 0) {
                              indx = 0;
                              break;
                         }
                         if(atom_atno[(connect[i][j])] == 1)
                              indx++;
                    }
                    if(indx == 1)
                         improper[i] = 1;
               }
          if((atom_atno[i] == 7) && (bond_count[i] == 3) &&
             ((ar_set[i][1] >= 1) || (ar_set[i][2] >= 1) ||
              (ar_set[i][3] >= 1))) improper[i] = 1;

          if(!(strcmp(&atom_type[i][0], "C.2"))) {		       
               indx = 0;
               for(j = 0; j < 3; j++) {
                    if((tmpint = connect[i][j]) == -1) break;
                    if(atom_atno[tmpint] == 7) {
                         indx = 1;
                         break;
                    }
               }
               if(indx) improper[i] = 1;
          }

          if((atom_atno[i] == 7) && (bond_count[i] == 3) &&
             (!strcmp(&atom_type[i][0], "N.pl3") ||
              !strcmp(&atom_type[i][0], "N.am")))
                    improper[i] = 1;

          if((atom_atno[i] == 7) && !(ar_set[i][1]) && !(ar_set[i][2]) &&
              !(ar_set[i][3]) && (bond_count[i] == 3)) {
                    indx = 0;
                    for(j = 0; j < 3; j++) {
                         if((tmpint = connect[i][j]) == -1) break;
                         if((ar_set[tmpint][1] >= 1) || (ar_set[tmpint][2] >= 1) ||
                            (ar_set[tmpint][3] >= 1)) {
                                   indx = 1;
                                   break;
                         }
                    }
                    if(indx) improper[i] = 1;
          }


          if(improper[i]) {
               if((connect[i][0] < 0) || (connect[i][1] < 0) ||
                  (connect[i][2] < 0)) continue;
               improper_num++;
          }
     }

/* make a table of atom numbers in each improper candidate */

     if(improper_num) {                              /* skip if no candidates */
          improperid = (int **)allocator(improper_num, sizeof(int*), FARGS);

          for(i = 0; i < improper_num; i++)
               improperid[i] = (int *)allocator(4, sizeof(int), FARGS);

          j = 0;
          for(i = 0; i < atom_num; i++) {               
               if(improper[i]) {
                    if((connect[i][0] < 0) || (connect[i][1] < 0) ||
                       (connect[i][2] < 0)) continue;                  /* should never take this branch */
/* special needs for cases of tetrahedral C with single H in united atoms models */
                    if(do_ua)
                         if((atom_atno[i] == 6) && (bond_count[i] == 4)) {
                              if(connect[i][3] < 0) continue;          /* should never take this branch */
                              tmpint = 0;
                              for(indx = 0; indx < 4; indx++) {
                                   if(atom_atno[(connect[i][indx])] == 1)
                                        continue;
                                   improperid[j][tmpint] = connect[i][indx];
                                   tmpint++;
                                   if(tmpint == 3) break;
                              }
                              improperid[j][3] = improperid[j][2];     /* adjust positions */
                              improperid[j][2] = i;
                              j++;
                              if(j > improper_num ) {
                                   sprintf(mess, "Improper locus counts mismatch %d > %d\n",
                                           j, improper_num);
                                   my_fatal(FARGS, mess);
                                   exit(1);
                              }
                              continue;                                /* skip the rest */
                         }
/* end special case in united atoms models */
                    improperid[j][2] = i;
                    improperid[j][0] = connect[i][0];
                    improperid[j][1] = connect[i][1];
                    improperid[j][3] = connect[i][2];
                    j++;
               }
               if(j > improper_num ) {
                    sprintf(mess, "Improper locus counts mismatch %d > %d\n",
                             j, improper_num);
                    my_fatal(FARGS, mess);
                    exit(1);
               }
          }
     }

     return;
}
