/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to judge the bond type

   NOTICE: This is a derivative work based on study of the finalize
   function found in the antechamber 1.27 program bondtype, version
   1.0, dated October, 2001 by Junmei Wang, Department of Pharmaceutical
   Chemistry, School of Pharmacy, University of California,
   San Francisco, CA  94143

   Portions of this work correct errors found in the code of the
   above cited work.  Other changes simplify storage allocation,
   and clarify decision trees.
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "judge_bond.h"

/* basic bond type judgment
   parameters are:
        num_bonds            number of bonds
        num_atoms            number of atoms
        num_rings            number of rings
*/
void judge_bond(int num_bonds, int num_atoms, int num_rings)
{
     int i, j, k, ibond, jbond, flag, flag0, flag1, flag2;
     int num, ichg, jchg, kchg, index[11];
     int *atom_junct;

/* set junctions index */

     atom_junct = (int *)allocator(num_atoms, sizeof(int), FARGS);

     for(i = 0; i < num_atoms; i++) {
          for(j = 0; j < 11; j++) index[j] = 0;       /* reset index for each atom */

          switch(atom_atno[i]) {
               case 6:  is_carbon(i, index);
                        break;

               case 7:  is_nitrogen(i, index);
                        break;

               case 8:  is_S_or_O(i, index, 8);
                        break;

               case 15: is_P_or_S(i, index, 2);
                        break;

               case 16: is_P_or_S(i, index, 3);
                        is_S_or_O(i, index, 9);
                        break;

               default: break;
          }

          atom_junct[i] = 0;
          if((index[1] >= 2) || (index[2] >= 2)) atom_junct[i] = 1;
          if((index[3] >= 2) || (index[4] >= 2)) atom_junct[i] = 1;
     }                           /* end for(i = 0; i < atom_num; i++) */

/* correct aromatic bond related connection types */

     for(i = 0; i < num_bonds; i++) {
          if(how_bonded[i] != 10) continue;             /* only for aromatic */
          ibond = bond_i[i];
          jbond = bond_j[i];
          flag0 = 0;
          for(j = 0; j < num_bonds; j++) {
               if(i == j) continue;
               if((bond_i[j] == ibond) ||
                  (bond_j[j] == ibond) ||
                  (bond_i[j] == jbond) ||
                  (bond_j[j] == jbond))
                         if((how_bonded[j] == 2) ||
                            (how_bonded[j] == 8)) {
                                   flag0 = 1;
                                   break;
                         }
          }                    /* end for(j = 0; j < num_bonds; j++) */

          if(flag0) how_bonded[i] = 7;
          else how_bonded[i] = 8;
     }                         /* end for(i = 0; i < num_atoms; i++) */

/* now determine bond types */

     for(i = 0; i < num_bonds; i++) {
          ichg = 0;
          jchg = 0;
          ibond = bond_i[i];
          jbond = bond_j[i];
          if((atom_charge[ibond] < -0.5) ||
             (atom_charge[ibond] > 0.5))
                    ichg = 1;
          if((atom_charge[jbond] < -0.5) ||
             (atom_charge[jbond] > 0.5))
                    jchg = 1;

/* part 1 */
          if(((ar_set[ibond][1] > 0)  &&
              (ar_set[jbond][1] > 0)) ||
             ((ar_set[ibond][2] > 0)  &&
              (ar_set[jbond][2] > 0))) {
                    for(j = 0; j < num_rings; j++) {
                         if(ring_num[j] >= 7) continue;
                         if(ring_num[j] <= 4) continue;
                         flag = 0;
                         for(k = 0; k < ring_num[j]; k++) {
                              if((ar_set[(ring_atomno[j][k])][1] <= 0) &&
                                 (ar_set[(ring_atomno[j][k])][2] <= 0)) {
                                      flag = 1;
                                      break;
                              }      /* end if((ar_set[(ring_atomno[j][k])][1] etc. */
                         }           /* end for(k = 0; k < ring_num[j]; k++) */
                         if(flag == 1) continue;
                         flag1 = 0;
                         flag2 = 0;
                         for(k = 0; k < ring_num[j]; k++) {
                              if(ring_atomno[j][k] == ibond) flag1 = 1;
                              if(ring_atomno[j][k] == jbond) flag2 = 1;
                         }     /* end for(k = 0; k < ring_num[j]; k++) */
                         if(flag1 && flag2) {
                              if(how_bonded[i] == 1)
                                   how_bonded[i] = 7;
                              if(how_bonded[i] == 2)
                                   how_bonded[i] = 8;
                              break;
                         }          /* end if(flag1 && flag2) */
                    }          /* end for(j = 0; j < num_rings; j++) */
          }                    /* end if((ar_set[ibond][1] > 0) etc. */

/* part 2 */
          if(how_bonded[i] == 2) {
               if(atom_junct[ibond]) {
                    if(((jchg + bond_count[jbond]) == 1) &&
                       ((atom_atno[jbond] == 8) ||
                        (atom_atno[jbond] == 16))) {
                              how_bonded[i] = 9;
                              continue;
                    }
               }               /* end if(atom_junct[ibond]) */
               if(atom_junct[jbond]) {
                    if(((ichg + bond_count[ibond]) == 1) &&
                       ((atom_atno[ibond] == 8) ||
                        (atom_atno[ibond] == 16))) {
                              how_bonded[i] = 9;
                              continue;
                    }
               }               /* end if(atom_junct[jbond]) */
          }                    /* end if(how_bonded[i] == 2) */

/* part 3 */
          if((((ichg + bond_count[ibond]) == 2) ||
              ((ichg + bond_count[ibond]) == 3)) &&
               (atom_atno[ibond] == 7) &&
             (((jchg + bond_count[jbond]) == 1) &&
              ((atom_atno[jbond] == 8) ||
               (atom_atno[jbond] == 16)))) {
                    for(j = 0; j < 5; j++) {
                         kchg = 0;
                         if(connect[ibond][j] < 0) break;
                         if(connect[ibond][j] == jbond) continue;
                         if((atom_charge[(connect[ibond][j])] < -0.5) ||
                            (atom_charge[(connect[ibond][j])] > 0.5))
                                   kchg = 1;
                         if(((kchg + bond_count[(connect[ibond][j])]) == 1) &&
                            ((atom_atno[(connect[ibond][j])] == 8) ||
                            (atom_atno[(connect[ibond][j])] == 16))) {
                                 how_bonded[i] = 6;
                                 break;
                         }      /* end if(((kchg + bond_count[(connect[ibond][j])]) == 1) etc. */
                    }           /* end for(j = 0; j < 5; j++) */
                    if(how_bonded[i] == 6)
                         continue;
          }                     /* end if((((ichg + bond_count[ibond]) == 2) etc. */

          if((((jchg + bond_count[jbond]) == 2) ||
              ((jchg + bond_count[jbond]) == 3)) &&
               (atom_atno[jbond] == 7) &&
             (((ichg + bond_count[ibond]) == 1) &&
              ((atom_atno[ibond] == 8) ||
               (atom_atno[ibond] == 16)))) {
                    for(j = 0; j < 5; j++) {
                         kchg = 0;
                         if(connect[jbond][j] < 0) break;
                         if((atom_charge[(connect[jbond][j])] < -0.5) ||
                            (atom_charge[(connect[jbond][j])] > 0.5))
                                   kchg = 1;
                         if(connect[jbond][j] == ibond) continue;
                         if(((kchg + bond_count[(connect[jbond][j])]) == 1) &&
                            ((atom_atno[(connect[jbond][j])] == 8) ||
                            (atom_atno[(connect[jbond][j])] == 16))){
                                 how_bonded[i] = 6;
                                 break;
                         }      /* end if(((kchg + bond_count[(connect[jbond][j])]) == 1) etc. */
                    }           /* end for(j = 0; j < 5; j++) */
                    if(how_bonded[i] == 6) continue;
          }                     /* end if((((jchg + bond_count[jbond]) == 2) etc. */

/* part 4 */
          if(how_bonded[i] == 1) {
               if((((atom_atno[ibond] == 8) ||
                    (atom_atno[ibond] == 16)) &&
                  ((ichg + bond_count[ibond]) == 1)) ||
                  (((atom_atno[jbond] == 8) ||
                    (atom_atno[jbond] == 16)) &&
                  ((jchg + bond_count[jbond]) == 1))) {
                       how_bonded[i] =9;
                       continue;
               }
          }
     }                          /* end for(i = 0; i < num_bonds; i++) */

/* part 5 */
     for(i = 0; i < num_atoms; i++) {
          ichg = 0;
          if((atom_charge[i] < -0.5) ||
             (atom_charge[i] > 0.5))
                    ichg = 1;
          if(((ichg + bond_count[i]) == 3) &&
              (atom_atno[i] == 16)) {
                    num = 0;
                    for(j = 0; j <= 5; j++)
                         index[j] = 0;
                    for(j = 0; j <= 2; j++) {
                         ibond = connect[i][j];
                         kchg = 0;
                         if((atom_charge[ibond] < -0.5) ||
                            (atom_charge[ibond] > 0.5))
                                   kchg = 1;
                         if((bond_count[ibond] == 1) &&
                            ((atom_atno[ibond] == 8) ||
                             (atom_atno[ibond] == 16))) {
                                   num++;
                                   if(!kchg) index[j] = 1;
                         }
                    }          /* end for(j = 0; j <= 5; j++) */

                    if(num == 2) {
                         for(j = 0; j < num_bonds; j++) {
                              ibond = bond_i[j];
                              jbond = bond_j[j];
                              for(k = 0; k <=2; k++)
                                   if(index[k] &&
                                     (((ibond == i) &&
                                       (jbond == connect[i][k])) ||
                                      ((jbond == i) &&
                                       (ibond == connect[i][k]))))
                                                how_bonded[j] = 9;
                         }     /* end for(j = 0; j < num_bonds; j++) */
                    }          /* end if(num == 2) */
          }                    /* end if(((ichg + bond_count[i]) == 3) etc. */

          if(((ichg + bond_count[i]) == 4) &&
              (atom_atno[i] == 15)) {
                    num = 0;
                    for(j = 0; j <= 5; j++)
                         index[j] = 0;
                    for(j = 0; j <= 3; j++) {
                         ibond = connect[i][j];
                         kchg = 0;
                         if((atom_charge[ibond] < -0.5) ||
                            (atom_charge[ibond] > 0.5))
                                   kchg = 1;
                         if((bond_count[ibond] == 1) &&
                            ((atom_atno[ibond] == 8) ||
                             (atom_atno[ibond] == 16))) {
                                   num++;
                                   if(!kchg) index[j] = 1;
                         }
                    }          /* end for(j = 0; j <= 5; j++) */

                    if(num >= 2) {
                         for(j = 0; j < num_bonds; j++) {
                              ibond = bond_i[j];
                              jbond = bond_j[j];
                              for(k = 0; k <=3; k++)
                                   if(index[k] &&
                                      (((ibond == i) &&
                                        (jbond == connect[i][k])) ||
                                       ((jbond == i) &&
                                        (ibond == connect[i][k]))))
                                                how_bonded[j] = 9;
                         }     /* end for(j = 0; j < num_bonds; j++) */
                    }          /* end if(num >= 2) */
          }                    /* end if(((ichg + bond_count[i]) == 4) etc. */

     }                         /* end for(i = 0; i < num_atoms; i++) */

     free_me(atom_junct, FARGS);         /* not needed after bond types set */

     return;
}

/* function to set bond type indices for carbon atoms
   parameters are:
        i                                atom number
        index                            array of indices to set
*/
void is_carbon(int i, int index[])
{
     int j;

     if(bond_count[i] == 3) {
          for(j = 0; j < 3; j++) {
               if((atom_atno[(connect[i][j])] == 8) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[1]++;
               if((atom_atno[(connect[i][j])] == 16) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[1]++;
          }
     }

     if(bond_count[i] == 1) {
          if((atom_atno[(connect[i][0])] == 7) &&
             (bond_count[(connect[i][0])] == 2))
                    index[10] = 1; 
     }

     return;
}

/* function to set bond type indices for nitrogen atoms
   parameters are:
        i                                atom number
        index                            array of indices to set
*/
void is_nitrogen(int i, int index[])
{
     int j;

     if(bond_count[i] == 1) {
          if((atom_atno[(connect[i][0])] == 7) &&
             (bond_count[(connect[i][0])] == 2))
                    index[4]++;
     }

     if(bond_count[i] == 2) {
          if(((atom_atno[(connect[i][0])] == 7) ||
              (atom_atno[(connect[i][0])] == 6)) &&
              (bond_count[(connect[i][0])] == 2))
                    index[5]++;
          if(((atom_atno[(connect[i][1])] == 7) ||
              (atom_atno[(connect[i][1])] == 6)) &&
              (bond_count[(connect[i][1])] == 2))
                    index[5]++;
     }

     if(bond_count[i] == 3) {
          for(j = 0; j < 3; j++) {
               if((atom_atno[(connect[i][j])] == 8) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[6]++;
               if((atom_atno[(connect[i][j])] == 16) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[6]++;
          }
     }

     if((bond_count[i] == 3) && (index[6] < 2)) {
          for(j = 0; j < 3; j++) {
               if((atom_atno[(connect[i][j])] == 8) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[7] = 1;
               if((atom_atno[(connect[i][j])] == 16) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[7] = 1;
          }
     }

     return;
}

/* function to set bond type indices for phosphorous or sulfur atoms
   parameters are:
        i                                atom number
        index                            array of indices to set
        k                                index type to set
*/
void is_P_or_S(int i, int index[], int k)
{
     int j;

     if(bond_count[i] == 4) {
          for(j = 0; j < 4; j++) {
               if((atom_atno[(connect[i][j])] == 8) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[k]++;
               if((atom_atno[(connect[i][j])] == 16) &&
                  (bond_count[(connect[i][j])] == 1))
                         index[k]++;
          }
     }

     return;
}

/* function to set bond type indices for sulfur or oxygen atoms
   parameters are:
        i                                atom number
        index                            array of indices to set
        k                                index type to set
*/
void is_S_or_O(int i, int index[], int k)
{
     int j, con[3];

     if(bond_count[i] == 1) {
          if((atom_atno[(connect[i][0])] == 7) &&
             (bond_count[(connect[i][0])] == 3)) {
                    for(j = 0; j < 3; j++)
                         con[j] = connect[(connect[i][0])][j];
                    for(j = 0; j < 3; j++) {
                         if((atom_atno[(con[j])] == 8) &&
                            (bond_count[(con[j])] == 1))
                                   index[k]++;
                         if((atom_atno[(con[j])] == 16) &&
                            (bond_count[(con[j])] == 1))
                                   index[k]++;
                    }
                    if(index[k] <= 1)
                         index[k] = 1;
                    else
                         index[k] = 0;
          }
     }

     return;
}
