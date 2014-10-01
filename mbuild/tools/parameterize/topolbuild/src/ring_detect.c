/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to detect ring structures   ring_detect.c

   NOTICE: This is a derivative work based on study of the routines
   found in the antechamber 1.27 function file rings.c,
   by Junmei Wang, Department of Pharmaceutical Chemistry, School
   of Pharmacy, University of California, San Francisco, CA  94143

   Portions of this work simplify storage allocation, and clarify
   decision trees compared to the work cited above.
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "ring_detect.h"
#include "mol2.h"
#include "rings.h"

#define MAXNAME 128

/* detect rings in structure
   parameters are:
        atom_num                   number of atoms in molecule
        bond_num                   number of bonds in molecule
        max_rng                    maximum possible number of rings in molecule
        num_rings                  return count of actual number of rings
*/
void detect_ring(int atom_num, int bond_num, int max_rng, int *num_rings)
{
     int *sel_arom_ind;
     int *sel_arom;
     int i, j, cnt_ring, nmbr;

     arom_num = (int *)allocator(atom_num, sizeof(int), FARGS);
     arom_atomno = (int **)allocator(atom_num, sizeof(int*), FARGS);

     for(i = 0; i < atom_num; i++) {
          arom_num[i] = 1;
          atom_sat[i] = 0;
          atom_ewd[i] = 0;
          arom_atomno[i] = (int *)allocator(12, sizeof(int), FARGS);
          for(j = 0; j <= 10; j++)
               arom_atomno[i][j] = 0;
     }

     ar_set = (int **)allocator(atom_num, sizeof(int*), FARGS);
     for(i = 0; i < atom_num; i++){
          ar_set[i] = (int *)allocator(6, sizeof(int), FARGS);
          for(j = 0; j < 6; j++ ) ar_set[i][j] = 0;
     }

     if(!max_rng) {
          *num_rings = 0;
          return;
     }
     sel_arom_ind = (int *)allocator(atom_num, sizeof(int), FARGS);
     sel_arom = (int *)allocator(atom_num, sizeof(int), FARGS);
     t_ring_num = (int *)allocator(max_rng, sizeof(int), FARGS);
     t_ring_atomno = (int **)allocator(max_rng, sizeof(int*), FARGS);
     t_ring_type = (short *)allocator(max_rng, sizeof(short), FARGS);

     for(i = 0; i < max_rng; i++) {
          t_ring_num[i] = 0;
          t_ring_type[i] = 0;
          t_ring_atomno[i] = (int *)allocator(12, sizeof(int), FARGS);
          for(j = 0; j < 10; j++)
               t_ring_atomno[i][j] = -1;
     }

     cnt_ring = 0;
     *num_rings = 0;
     for(i = 0; i < atom_num; i++) {
          if((atom_atno[i] != 6)  &&
             (atom_atno[i] != 7)  &&
             (atom_atno[i] != 8)  &&
             (atom_atno[i] != 16) &&
             (atom_atno[i] != 15))
                    continue;
          if((atom_atno[i] == 6) &&
             (bond_count[i] <= 2))
                    continue;
          if((atom_atno[i] == 8) &&
             (bond_count[i] == 1))
                    continue;
          if((atom_atno[i] == 16) &&
             (bond_count[i] == 1))
                    continue;

          for(j = 0; j < atom_num; j++)
               sel_arom_ind[j] = -1;
          for(j = 0; j < 10; j++)
               sel_arom[j] = -1;

          cycle_detect(&cnt_ring, 0, i, sel_arom, sel_arom_ind, max_rng);

          if(cnt_ring >= max_rng) {
               sprintf(mess,
                  "In ring detection, no. detected rings, %d, too large.\n",
                  cnt_ring);
               my_fatal(FARGS, mess);
               exit(1);
          }
          cycle_clean(&cnt_ring);
     }

     cycle_clean(&cnt_ring);
     cycle_props(cnt_ring);
     aromaticity(atom_num, bond_num, cnt_ring);

     free_me(sel_arom, FARGS);
     free_me(sel_arom_ind, FARGS);

     if(cnt_ring) {
          ring_num = (int *)allocator(cnt_ring, sizeof(int), FARGS);
          ring_atomno = (int **)allocator(cnt_ring, sizeof(int*), FARGS);
          ring_type = (short *)allocator(cnt_ring, sizeof(short), FARGS);
          for(i = 0; i < cnt_ring; i++) {
               ring_atomno[i] = (int *)allocator(12, sizeof(int), FARGS);
               ring_num[i] = t_ring_num[i];
               ring_type[i] = t_ring_type[i];
               for(j = 0; j < 10; j++)
                    ring_atomno[i][j] = t_ring_atomno[i][j];
          }
     }

     free_me(t_ring_num, FARGS);
     free_me(t_ring_type, FARGS);
     for(i = 0; i < max_rng; i++)
          free_me(t_ring_atomno[i], FARGS);
     free_me(t_ring_atomno, FARGS);

     *num_rings = cnt_ring;
     return;
}

/* recursive determination of aromatic atoms
   parameters are:
        rings                  number of rings
        selected               current position in proposed cycle
        st_num                 current atom number to work on
        arom                   array of atom numbers proposed for cycle
        arom_ind               marker that atom is proposed for cycle 
        max_rng                maximum number of rings allowed
*/
void cycle_detect(int *rings, int selected, int st_num, int *arom,
                  int *arom_ind, int max_rng)
{
     int i, j, k, m;
     int start;
     int breaker = 0;

     start = -1;
     arom[selected++] = st_num;
     arom_ind[st_num] = 1;

     for(i = 0; i < 4; i++) {
          if(connect[st_num][i] < 0) return;
          start = connect[st_num][i];
          if((atom_atno[start] != 6)  &&
             (atom_atno[start] != 7)  &&
             (atom_atno[start] != 8)  &&
             (atom_atno[start] != 16) &&
             (atom_atno[start] != 15))
                    continue;
          if((atom_atno[start] == 6) &&
             (bond_count[start] <= 2))
                    continue;
          if((atom_atno[start] == 8) &&
             (bond_count[start] == 1))
                    continue;
          if((atom_atno[start] == 16) &&
             (bond_count[start] == 1))
                    continue;

          for(j = 0; j < selected; j++) {
               if(arom[j] == start) {
                    breaker = 1;
                    break;
               }
          }                              /* end for(j = 0; j < selected; j++) */
          if(breaker) {
               breaker = 0;
               continue;
          }
          if(selected > 10) return;     /* have already been here */

          for(j = 2; j < 10; j++) {     /* deal with 3 to 10 member ring */
               if(selected == j) {
                    for(k = 0; k < 4; k++) {
                         if(connect[(arom[0])][k] < 0) break;
                         if(connect[(arom[0])][k] == start) {
                              for(m = 0; m < j; m++)
                                   t_ring_atomno[(*rings)][m] = arom[m];
                              t_ring_atomno[(*rings)][j] = start;
                              t_ring_num[(*rings)] = j + 1;
                              (*rings)++;
                              if((*rings) >= max_rng) {
                                   sprintf(mess,
                                      "In ring detection, no. detected rings, %d, too large.\n",
                                           (*rings));
                                   my_fatal(FARGS, mess);
                                   exit(1);
                              }
                              break;
                         }          /* end if(connect[(arom[0])][k]) == start) */
                    }               /* end for(k = 0; k < 4; k++) */
               }                    /* end if(selected == j) */
          }                         /* end for(j = 2; j < 10; j++) */
          if(start == -1) return;   /* did not find anything */

          cycle_detect(rings, selected, start, arom, arom_ind, max_rng);
     }                              /* end for(i = 0; i < 4; i++) */

     return;
}

/* order ring atoms and remove duplicated rings
   parameter is
        ringers                  number of rings
*/
void cycle_clean(int *ringers)
{
     int i, j, k;
     int *back_num;
     int **back_atoms;
     int backing_no = 0;
     int tmpint, indx, rings;

     rings = (*ringers);

     if(!rings) return;      /* don't bother if no rings */

     back_num = (int *)allocator(rings, sizeof(int), FARGS);
     back_atoms = (int **)allocator(rings, sizeof(int*), FARGS);
     for(i = 0; i < rings; i++)
          back_atoms[i] = (int *)allocator(12, sizeof(int), FARGS);

     for(i = 0; i < rings; i++) {              /* put ring atom numbers in order */
          for(j = 0; j < t_ring_num[i]; j++) {
               for(k = (j + 1); k < t_ring_num[i]; k++) {
                    if(t_ring_atomno[i][j] > t_ring_atomno[i][k]) {
                         tmpint = t_ring_atomno[i][k];
                         t_ring_atomno[i][k] = t_ring_atomno[i][j];
                         t_ring_atomno[i][j] = tmpint;
                    }     /* end if(t_ring_atomno[i][j] > t_ring_atomno[i][k]) */
               }          /* end for(k = (j + 1); k < t_ring_num[i]; k++) */
          }               /* end for(j = 0; j < t_ring_num[i]; j++) */
     }                    /* end for(i = 0; i < rings; i++) */

     for(i = 0; i < rings; i++) {
          for(j = (i + 1); j < rings; j++) {
               if((t_ring_num[i] == t_ring_num[j]) &&
                  (t_ring_num[i] != 0)) {
                         indx = 1;
                         for(k = 0; k < t_ring_num[i]; k++) {
                              if(t_ring_atomno[i][k] != t_ring_atomno[j][k]) {
                                   indx = 0;
                                   break;
                              }   /* end if(t_ring_atomno[i][k] != t_ring_atomno[j][k]) */
                         }        /* end for(k = 0; k < t_ring_num[i]; k++) */
                         if(indx)
                              t_ring_num[j] = 0;
               }               /* end if((t_ring_num[i] == t_ring_num[j]) &&
                                         !(t_ring_num[i])) */
          }                    /* end for(j = (i + 1); j < rings; j++) */
     }                         /* end for(i = 0; i < rings; i++) */

     for(i = 0; i < rings; i++) {
          for(j = 0; j < t_ring_num[i]; j++) {
               indx = 0;
               for(k = 0; k < t_ring_num[i]; k++) {
                    if(connect[(t_ring_atomno[i][j])][0] == t_ring_atomno[i][k])
                         indx++;
                    if(connect[(t_ring_atomno[i][j])][1] == t_ring_atomno[i][k])
                         indx++;
                    if(connect[(t_ring_atomno[i][j])][2] == t_ring_atomno[i][k])
                         indx++;
                    if(connect[(t_ring_atomno[i][j])][3] == t_ring_atomno[i][k])
                         indx++;
                    if(connect[(t_ring_atomno[i][j])][4] == t_ring_atomno[i][k])
                         indx++;
                    if(connect[(t_ring_atomno[i][j])][5] == t_ring_atomno[i][k])
                         indx++;
               }          /* end for(k = 0; k < t_ring_num[i]; k++) */

               if(indx == 3) {
                    t_ring_num[i] = 0;
                    break;

               }
          }               /* end for(j = 0; j < t_ring_num[i]; j++) */
     }                    /* end for(i = 0; i < rings; i++) */

     for(i = 0; i < rings; i++) {
          if(t_ring_num[i]) {
               for(j = 0; j < t_ring_num[i]; j++)
                    back_atoms[backing_no][j] = t_ring_atomno[i][j];
               back_num[backing_no] = t_ring_num[i];
               backing_no++;
          }               /* end if(t_ring_num[i]) */
     }                    /* end for(i = 0; i < rings; i++) */

     for(i = 0; i < backing_no; i++) {
          t_ring_num[i] = back_num[i];
          for(j = 0; j < t_ring_num[i]; j++)
               t_ring_atomno[i][j] = back_atoms[i][j];
     }                    /* end for(i = 0; i < backing_no; i++) */

     free_me(back_num, FARGS);
     for(i = 0; i < rings; i++)
          free_me(back_atoms[i], FARGS);
     free_me(back_atoms, FARGS);

     *ringers = backing_no;                 /* correct the count of rings */

     return;
}

/* set counters for aromaticity checks
   parameter is:
        rings                           count of number of rings
*/
void cycle_props(int rings)
{
     int i, j, tmpint;

     if(!rings) return;

     for(i = 0; i < rings; i++) {
          for(j = 0; j < t_ring_num[i]; j++) {
               tmpint = t_ring_atomno[i][j];
               arom_atomno[tmpint][0]++;
               arom_atomno[tmpint][(t_ring_num[i])]++;
          }          /* end for(j = 0; j < t_ring_num[i]; j++) */
     }               /* end for(i = 0; i < rings; i++) */

     return;
}

/* check aromaticities of rings found
   parameters are:
        no_atms          number of atoms in molecule
        no_bnds          number of bonds in molecule
        no_rings         number of rings in molecule
*/
void aromaticity(int no_atms, int no_bnds, int no_rings)
{
     int i, j, k, tmpint, tmpint1, indx, indx0;
     int *init_arom;
     int cont_ind = 0;

     init_arom = (int *)allocator(no_atms, sizeof(int), FARGS);

     for(i = 0; i < no_atms; i++) {
          init_arom[i] = 0;
          ar_set[i][1] = 0;
          ar_set[i][2] = 0;
          ar_set[i][3] = 0;
          ar_set[i][4] = 0;
          ar_set[i][5] = 0;
          atom_sat[i] = -1;

          switch(atom_atno[i]) {
               case 6:   atom_ewd[i] = 0;
                         if(bond_count[i] == 3)
	                      init_arom[i] = 2;
                         if(bond_count[i] == 4) {
                              init_arom[i] = -2;
                              atom_sat[i] = 1;
                         }
                         break;

               case 7:   atom_ewd[i] = 1;
                         if(bond_count[i] <= 3)
                              init_arom[i] = 2;
                         if(bond_count[i] >= 3)
                              atom_sat[i] = 1;
                         break;

               case 8:   atom_ewd[i] = 1;
                         if(bond_count[i] == 2) {
                              init_arom[i] = 1;
                              atom_sat[i] = 1;
                         }
                         break;

               case 15:  atom_ewd[i] = 0;
                         if(bond_count[i] == 2)
                              init_arom[i] = 2;
                         if(bond_count[i] > 2)
                              init_arom[i] = 1;
                         if(bond_count[i] >= 3)
                              atom_sat[i] = 1;
                         break;

               case 16:  atom_ewd[i] = 1;
                         if(bond_count[i] >= 2)
                              init_arom[i] = 1;
                         if(bond_count[i] >= 3)
                              atom_sat[i] = 1;
                         break;

               case 9:
               case 17:
               case 35:
               case 53:  atom_ewd[i] = 1;
                         atom_sat[i] = 1;
                         break;

                case 1:  atom_ewd[i] = 0;
                         atom_sat[i] = 1;
                         break;

               default:  atom_ewd[i] = 0;
                         break;
          }          /* end switch(atom_atno[i]) */
     }               /* end for(i = 0; i < no_atms; i++) */

     for(i = 0; i < no_rings; i++) {
          tmpint = 0;
          for(j = 0; j < t_ring_num[i]; j++)
               tmpint += init_arom[(t_ring_atomno[i][j])];

/* pure aliphatic rings */
          if(tmpint == ( -2 * t_ring_num[i])) {
               for(j = 0; j < t_ring_num[i]; j++)
                    ar_set[(t_ring_atomno[i][j])][5]++;
               t_ring_type[i] = 5;
               continue;
          }          /* end if(tmpint == ( -2 * t_ring_num[i])) */

/* for ring containing sp3 carbon */
          for(k = 0; k < t_ring_num[i]; k++) {
               if(init_arom[(t_ring_atomno[i][k])] < 0) {
                    for(j = 0; j < t_ring_num[i]; j++)
                         ar_set[(t_ring_atomno[i][j])][4]++;
                    t_ring_type[i] = 4;
                    cont_ind = 1;
                    break;
               }     /* end if(init_arom[(t_ring_atomno[i][k])] < 0) */
          }          /* end for(k = 0; k < t_ring_num[i]; k++) */
          if(cont_ind == 1) {
               cont_ind = 0;
               continue;
          }          /* end if(cont_ind == 1) */

/* for planar ring and outside double bonds */
          if((tmpint >= t_ring_num[i]) &&
             (tmpint <= ( 2 * t_ring_num[i]))) {
                    for(j = 0; j  < no_bnds; j++) {
                         indx = 0;
                         for(k = 0; k < t_ring_num[i]; k++) {
                              if(bond_i[j] == t_ring_atomno[i][k])
                                   if((arom_atomno[(bond_j[j])][0]) == 0)
                                        indx++;
                              if(bond_j[j] == t_ring_atomno[i][k])
                                   if((arom_atomno[(bond_i[j])][0]) == 0)
                                        indx++;
                         }     /* end for(k = 0; k < t_ring_num[i]; k++) */
                         if((indx == 1) &&
                            ((how_bonded[j] == 2) ||
                             (how_bonded[j] == 8))) {
                                   for(k = 0; k < t_ring_num[i]; k++)
                                        ar_set[(t_ring_atomno[i][k])][3]++;
                                   t_ring_type[i] = 3;
                                   cont_ind = 1;
                                   break;
                         }     /* end if((indx ==1) etc.  */
                    }          /* end for(j = 0; j  < no_bnds; j++) */
          }               /* end if((tmpint >= t_ring_num[i]) &&
                                    (tmpint <= ( -2 * t_ring_num[i]))) */
          if(cont_ind == 1) {
               cont_ind = 0;
               continue;
          }          /* end if(cont_ind == 1) */

/* Pure Aromatic Rings */
          if((tmpint == 12) &&
             (t_ring_num[i] == 6)) {
                    indx = 0;
                    for(j = 0; j < t_ring_num[i]; j++) {
                         if((atom_atno[(t_ring_atomno[i][j])] == 7) ||
                            (atom_atno[(t_ring_atomno[i][j])] == 15)) {
                                   tmpint1 = t_ring_atomno[i][j];
                                   indx0 = 0;
                                   for( k = 0; k < no_bnds; k++) {
                                        if((bond_i[k] == tmpint1) &&
                                           ((how_bonded[k] == 8) ||
                                            (how_bonded[k] == 2) ||
                                            (how_bonded[k] == 10)))
                                                  indx0 = 1;
                                        if((bond_j[k] == tmpint1) &&
                                           ((how_bonded[k] == 8) ||
                                            (how_bonded[k] == 2) ||
                                            (how_bonded[k] == 10)))
                                                  indx0 = 1;
                                   }     /* end for( k = 0; k < no_bnds; k++) */
                                   if(!indx0) indx = 1;
                         }     /* end if((atom_atno[(t_ring_atomno[i][j])] == 7) ||
                                         (atom_atno[(t_ring_atomno[i][j])] == 15) */
                    }          /* end for(j = 0; j < t_ring_num[i]; j++) */

                    if(!indx) {
                         for(j = 0; j < t_ring_num[i]; j++)
                              ar_set[(t_ring_atomno[i][j])][1]++;
                         t_ring_type[i] = 1;
                         continue;
                    }          /* end if(!indx) */
          }                    /* end if((tmpint == 12) &&
                                         (t_ring_num[i] == 6) */

/* for other planar rings */
          if(tmpint >= (t_ring_num[i] + 3)) {
               for(j = 0; j < t_ring_num[i]; j++)
                    ar_set[(t_ring_atomno[i][j])][2]++;
               t_ring_type[i] = 2;
               continue;
          }                    /* end if(tmpint >= (t_ring_num[i] + 3) */

/* for the other rings */
          for(j = 0; j < t_ring_num[i]; j++)
               ar_set[(t_ring_atomno[i][j])][4]++;
          t_ring_type[i] = 4;
          continue;            /* Not strictly necessary, but placed here for
                                  similarity in structures */

     }                         /* end for(i = 0; i < no_rings; i++) */

     for(i = 0; i < no_atms; i++) {
          if(ar_set[i][1] > 0) {
               arom_num[i] = 0;
               continue;
          }                    /* indicator 1 set. */
          if(ar_set[i][2] > 0) {
               arom_num[i] = 0;
               continue;
          }                    /* indicator 2 set. */
          if(ar_set[i][3] > 0) {
               arom_num[i] = 0;
               continue;
          }                    /* indicator 3 set. */
          if(ar_set[i][4] > 0) {
               arom_num[i] = 0;
               continue;
          }                    /* indicator 4 set. */
          if(ar_set[i][5] > 0) {
               arom_num[i] = 0;
               continue;
          }                    /* indicator 5 set. */

     }                         /* end for(i = 0; i < no_atms; i++) */

     free_me(init_arom, FARGS);
     return;
}
