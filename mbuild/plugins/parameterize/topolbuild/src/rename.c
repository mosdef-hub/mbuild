/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
   Function to rename atoms according to a rule set
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

/* rename the atoms in this molecule according to a rule set
   parameters are:
        num_atoms               number of atoms in molecule
        num_rings               number of rings in molecule
*/
void renameit(int num_atoms, int num_rings)
{
     int i, j, k, sat_count, has_arom, has_hetro, indx_H;
     int atom_count, rptr, count2, start, start2;
     int unbound, has_prime, prime_count;
     int *indx_chain;
     char number[5];
     char *alpha_desig[] = { "a", "b", "c", "d" };

     new_names = (char **)mat_alloc(num_atoms, 12, FARGS);

/* rule 0: must have more than two main chain atoms */

     if(main_num < 3) {
          for(i = 0; i < num_atoms; i++)
               strcpy(new_names[i], atom_name[i]);
          return;
     }

     for(i = 0; i < num_atoms; i++) {
          strcpy(new_names[i], "OLD");
          strcat(new_names[i], atom_name[i]);
     }

     indx_chain = (int *)allocator(num_atoms, sizeof(int), FARGS);

     atom_count = 1;
     unbound = 1;
     has_prime = 0;
     prime_count = 1;

     for(i = 0; i < num_atoms; i++)
          indx_chain[i] = -1;

/* rule 1:  first atom must be in main chain and must have two or more atoms bound to
   it that are not H */

     for(i = 0; i < main_num; i++) {
          if(main_chain[i] < 0) continue;
          sat_count = 0;

          for(j = 0; j < bond_count[(main_chain[i])]; j++)
               if(!atom_atno[(connect[(main_chain[i])][j])] ||
                  (atom_atno[(connect[(main_chain[i])][j])] == 1))
                         sat_count++;

          if(bond_count[(main_chain[i])] <= (sat_count + 1))
               continue;                  /* if only one atom not H bound to this, skip */ 

          indx_chain[(main_chain[i])] = main_chain[i];
     }                                    /* end for(i = 0; i < main_num; i++) */

/* rule 2:  if there are rings in main chain, first atom must be in a ring and not at
            junction of rings */

     if(num_rings) {
          for(i = 0; i < num_atoms; i++) {
               if(indx_chain[i] < 0)
                    continue;             /* skip atoms not covered by rule 1 */
               if(arom_atomno[i][0] == 1)
                    continue;             /* retain atoms in rings and not at junctions */
               indx_chain[i] = -1;        /* do not consider atoms not in rings */
          }                               /* end for(i = 0; i < num_atoms; i++) */

/* rule 3:  if any of the rings are aromatic, first atom must be in an aromatic ring */

          has_arom = 0;
          for(i = 0; i < num_rings; i++) {
               for(j = 0; j < ring_num[i]; j++) {
                    if(indx_chain[(ring_atomno[i][j])] < 0)
                         continue;       /* skip atoms not covered by rules 1 and 2 */
                    if(ar_set[(ring_atomno[i][j])][1] >= 1) {
                         has_arom = 1;
                         break;
                    }
               }                         /* end for(j = 0; j < ring_num[i]; j++) */
               if(has_arom) break;
          }                              /* end for(i = 0; i < num_atoms; i++) */

          if(has_arom)
               for(i = 0; i < num_atoms; i++) {
                    if(indx_chain[i] < 0)
                         continue;       /* skip atoms not covered by rules 1, and 2 */
                    if(ar_set[i][1])
                         continue;       /* skip aromatic atoms */
                    indx_chain[i] = -1;  /* do not consider atoms not in aromatic rings */
               }                         /* end for(i = 0; i < num_atoms; i++) */

/* rule 4:  if any of the above are hetero-atoms, first atom must be a hetero-atom */

          has_hetro = 0;
          for(i = 0; i < num_atoms; i++) {
               if(indx_chain[i] < 0)
                    continue;            /* skip atoms not covered by previous rules */
               if(atom_atno[i] != 6) {
                    has_hetro = 1;
                    break;
               }                         /* end if(atom_atno[i] != 6) */
          }                              /* end for(i = 0; i < num_atoms; i++) */

          if(has_hetro)
               for(i = 0; i < num_atoms; i++) {
                    if(indx_chain[i] < 0)
                         continue;       /* skip atoms not covered by previous rules */
                    if(atom_atno[i] != 6)
                         continue;       /* skip aromatic hetero-atoms */
                    indx_chain[i] = -1;  /* do not consider non-aromatic hetero-atoms */
               }                         /* end for(i = 0; i < num_atoms; i++) */
     }                                   /* end if(num_rings) */

     for(i = 0; i < main_num; i++) {     /* assign new name to first atom */
          if(indx_chain[(main_chain[i])] < 0)
               continue;                 /* skip atoms not covered by previous rules */
          strcpy(new_names[(main_chain[i])], atom_symbl[(main_chain[i])]);
          sprintf(number, "%d", atom_count);
          atom_count++;
          strcat(new_names[(main_chain[i])], number);
          indx_chain[(main_chain[i])] = main_chain[i];       /* mark assignment made */
          k = i;
          break;
     }                                   /* end for(i = 0; i < main_num; i++) */

     start = -1;
     if(num_rings) {
          for(i = 0; i < num_rings; i++) {
               for(j = 0; j < ring_num[i]; j++)
                    if(ring_atomno[i][j] == main_chain[k]) {
                         start = j;
                         start2 = i;
                         break;
                    }
               if(start > -1) break;
          }

          for(j = 0; j < ring_num[start2]; j++) {
               if(j == start) continue;
               strcpy(new_names[(ring_atomno[start2][j])],
                      atom_symbl[(ring_atomno[start2][j])]);
               sprintf(number, "%d", atom_count);
               atom_count++;
               strcat(new_names[(ring_atomno[start2][j])], number);
               indx_chain[(ring_atomno[start2][j])] = ring_atomno[start2][j];
          }
     }

     for(i = 0; i < main_num; i++) {     /* work way down main chain assigning new names */
          if(strncmp(new_names[(main_chain[i])], "OLD", 3))
               continue;                 /* skip if already set name */

          sat_count = 0;
          for(j = 0; j < bond_count[(main_chain[i])]; j++)
               if(!atom_atno[(connect[(main_chain[i])][j])] ||
                  (atom_atno[(connect[(main_chain[i])][j])] == 1))
                         sat_count++;

          if(bond_count[(main_chain[i])] <= (sat_count + 1))
               continue;                 /* if only one atom not H bound to this, skip */ 

          if(!has_prime &&
             arom_atomno[(main_chain[i])][0] &&
             (num_rings > 1))            /* find any rings that should be labeled as primed */
                    for(j = 1; j < num_rings; j++) {
                         for(k = 0; k < ring_num[j]; k++) {
                              if(arom_atomno[(ring_atomno[j][k])][0] != 1) {
                                   has_prime = 0;
                                   break;
                              }          /* end if(arom_atomno[(ring_atomno[j][k])][0] != 1) */
                              if(ring_atomno[j][k] == main_chain[i])
                                   has_prime = 1;
                         }               /* end for(k = 0; k < ring_num[j]; k++) */
                         if(has_prime) break;
                    }                    /* end for(j = 1; j < num_rings; j++) */

          if(has_prime) {                /* label primed rings as primed */
               if(!arom_atomno[(main_chain[i])][0])
                    has_prime = 0;  /* end the primed ring labeling at the first atom after it */
               strcpy(new_names[(main_chain[i])], atom_symbl[(main_chain[i])]);
               sprintf(number, "%d", prime_count);
               prime_count++;
               strcat(new_names[(main_chain[i])], number);
               strcat(new_names[(main_chain[i])], "p");
               indx_chain[(main_chain[i])] = main_chain[i];    /* mark assignment made */
               continue;
          }

          strcpy(new_names[(main_chain[i])], atom_symbl[(main_chain[i])]);
          sprintf(number, "%d", atom_count);
          atom_count++;
          strcat(new_names[(main_chain[i])], number);
          indx_chain[(main_chain[i])] = main_chain[i];         /* mark assignment made */
     }                                   /* end for(i = j; i < main_num; i++) */

     rptr = 1;
     while(rptr) {
          rptr = 0;
          for(i = 0; i < num_atoms; i++) {
               if(indx_chain[i] != -1) continue;

               if(bond_count[i] == 0) {
                    strcpy(new_names[i], atom_symbl[i]);
                    if(unbound > 1) {
                         sprintf(number, "%d", unbound);
                         strcat(new_names[i], number);
                    }
                    unbound++;
                    indx_chain[i] = i;
                    continue;
               }

               if(!atom_atno[i] || (atom_atno[i] == 1)) {     /* special naming rules for H */
                    j = connect[i][0];
                    if(indx_chain[j] < 0) {
                         rptr = 1;
                         continue;       /* have to repeat this section if I haven't assigned */
                    }                    /* new name to atom to which this H is bound */

                    count2 = strcspn(new_names[j], "0123456789");
                    indx_H = 0;
                    for(k = 0; k < bond_count[j]; k++)
                         if(!atom_atno[(connect[j][k])] ||
                            (atom_atno[(connect[j][k])] == 1))
                                   indx_H++;

                    if(indx_H == 1) {    /* label single H with H followed by number of attached */
                         strcpy(new_names[i], "H");
                         strcat(new_names[i], &new_names[j][count2]);
                         indx_chain[i] = i;
                         continue;       /* once I've assigned a name move on to next atom */
                    }                    /* end if(indx_H == 1) */

                    indx_H = 0;
                    for(k = 0; k < bond_count[j]; k++)
                         if(!atom_atno[(connect[j][k])] ||
                            (atom_atno[(connect[j][k])] == 1)) {
                                   strcpy(new_names[(connect[j][k])],
                                          "H");
                                   strcat(new_names[(connect[j][k])],
                                          &new_names[j][count2]);
                                   strcat(new_names[(connect[j][k])],
                                          alpha_desig[indx_H]);
                                   indx_H++;
                                   indx_chain[(connect[j][k])] = connect[j][k];
                         }               /* end if(!atom_atno[(connect[j][k])] etc. */
                    continue;            /* once I've assigned a name move on to next atom */
               }                         /* end if(!atom_atno[i] || (atom_atno[i] == 1)) */

               if(bond_count[i] > 0) {   /* label atom with only one bond to other than H */
                    sat_count = 0;       /* with number of non-H atom to which it is attached */
                    start = -1;          /* label attached H to this atom appropriately */
                    for(j = 0; j < bond_count[i]; j++)
                         if(!atom_atno[(connect[i][j])] ||
                            (atom_atno[(connect[i][j])] == 1))
                                   sat_count++;
                         else start = j;

                    if((start != -1) &&
                       (sat_count != 0) &&
                       (atom_atno[i] != atom_atno[(connect[i][start])]) &&
                       (bond_count[i] <= (sat_count + 1))) {
                              j = connect[i][start];
                              if(indx_chain[j] < 0) {
                                   rptr = 1;
                                   continue;  /* have to repeat this section if I haven't assigned */
                              }               /* new name to parent atom */
                              count2 = strcspn(new_names[j], "0123456789");
                              strcpy(new_names[i], atom_symbl[i]);
                              strcat(new_names[i], &new_names[j][count2]);
                              indx_chain[i] = i;
                              indx_H = 0;
                              for(k = 0; k < bond_count[j]; k++)
                                   if(atom_atno[(connect[j][k])] == 1)
                                        indx_H++;
                              if(indx_H)
                                   for(k = 0; k < bond_count[i]; k++)
                                        if(atom_atno[(connect[i][k])] == 1) {
                                             strcpy(new_names[(connect[i][k])],
                                                    "H");
                                             strcat(new_names[(connect[i][k])],
                                                    &new_names[j][count2]);
                                              strcat(new_names[(connect[i][k])],
                                                     alpha_desig[(indx_H - 1)]);
                                              indx_H++;
                                              indx_chain[(connect[i][k])] = connect[i][k];
                                        }     /* end if(!atom_atno[(connect[i][k])] etc. */
                              
                         continue;
                    }                    /* end if(bond_count[i] <= (sat_count + 1)) */
               }                         /* end if(bond_count[i] > 0) */

               if(bond_count[i] > 1) {   /* for main chain loops, label non-main chain */
                    start = 0;           /* atom in loop between two main chain atoms */
                    start2 = -1;         /* with higher atom number of those two */
                    for(j = 0; j < bond_count[i]; j++) {
                         if(!atom_atno[(connect[i][j])] ||
                            (atom_atno[(connect[i][j])] == 1)) continue;
                         for(k = 0; k < main_num; k++)
                              if(connect[i][j] == main_chain[k]) {
                                   start++;
                                   if(atom_atno[(connect[i][j])] == atom_atno[i])
                                        continue;
                                   if(start2 < main_chain[k])
                                        start2 = main_chain[k];
                              }          /* end if(connect[i][j] == main_chain[k]) */
                    }                    /* end for(j = 0; j < bond_count[i]; j++) */
                    if(start2 > -1) {
                         if((start == 2) &&
                            (atom_atno[start2] != atom_atno[i])) {
                                   strcpy(new_names[i], atom_symbl[i]);
                                   count2 = strcspn(new_names[start2], "0123456789");
                                   strcat(new_names[i], &new_names[start2][count2]);
                                   indx_chain[i] = i;
                                   continue;
                              }          /* end if((start == 2)  etc. */
                    }                    /* end if(start2 > -1) */
               }                         /* end if(bond_count[i] > 1) */

               strcpy(new_names[i], atom_symbl[i]);
               sprintf(number, "%d", atom_count);
               atom_count++;
               strcat(new_names[i], number);
               indx_chain[i] = i;
          }                              /* end for(i = 0; i < num_atoms; i++) */
     }                                   /* end while(rptr) */

     free_me(indx_chain, FARGS);
     return;
}
