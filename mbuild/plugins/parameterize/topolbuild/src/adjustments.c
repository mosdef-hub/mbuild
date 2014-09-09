/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to adjust general amber forcefield atom type assignments

   NOTICE: This is a derivative work based on study of routines found in
   the antechamber 1.27 program atomtype, version 1.0, dated October, 2001
   by Junmei Wang, Department of Pharmaceutical Chemistry, School of Pharmacy,
   University of California, San Francisco, CA  94143
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "adjustments.h"
#include "atom_types.h"
#include "mol2.h"
#include "rings.h"

#define MAXCHAR 256

short *adjindxa = NULL;
short *adjindxb = NULL;


/* adjust gaff types cc, ce, cg, pc, pe, nc, and ne to be
   gaff types cd, cf, ch, pd, pf, nd, or nf with some bonding
   situations
   parameters are:
        num_atoms                    number of atoms
        num_bonds                    number of bonds
        assigned_type                the types as assigned by type judging
*/
void adjustments(int num_atoms, int num_bonds, char **assigned_type)
{
     int i, indx, num, flag, bonda, bondb;

     i = num_atoms + 10;
     adjindxa = (short *)allocator(i, sizeof(short), FARGS);
     adjindxb = (short *)allocator(i, sizeof(short), FARGS);

     for(i = 0; i < num_atoms; i++) {
          adjindxa[i] = 0;
          adjindxb[i] = 0;
     }

     indx = 0;
     num = 0;

     for(i = 0; i < num_atoms; i++)
          if(!strncmp(assigned_type[i], "cc", 2) ||
             !strncmp(assigned_type[i], "ce", 2) ||
             !strncmp(assigned_type[i], "cg", 2) ||
             !strncmp(assigned_type[i], "pc", 2) ||
             !strncmp(assigned_type[i], "pe", 2) ||
             !strncmp(assigned_type[i], "nc", 2) ||
             !strncmp(assigned_type[i], "ne", 2)) {
                  adjindxb[i] = 1;
                  if(!indx) {
                       adjindxa[i] = 1;
                       indx = 1;
                  }
                  num++;
          }

     if(!num) {
          free_me(adjindxa, FARGS);
          free_me(adjindxb, FARGS);
          return;
     }

     num--;
     while(num > 0) {
          num--;
          flag = 0;
          for(i = 0; i < num_bonds; i++) {
               bonda = bond_i[i];
               bondb = bond_j[i];
               if((adjindxb[bonda] + adjindxb[bondb]) != 2) continue;
               if(!flag && !adjindxa[bonda] && !adjindxa[bondb])
                    adjindxa[bonda] = 1;
               if(!adjindxa[bonda] && adjindxa[bondb]) {
                    flag = 1;
                    if((how_bonded[i] == 1) ||
                       (how_bonded[i] == 7)) adjindxa[bonda] = adjindxa[bondb];
                    if((how_bonded[i] == 2) ||
                       (how_bonded[i] == 8) ||
                       (how_bonded[i] == 3)) adjindxa[bonda] = -adjindxa[bondb];
               }                    /* end if(!adjindxa[bonda] && adjindxa[bondb]) */

               if(!adjindxa[bondb] && adjindxa[bonda]) {
                    flag = 1;
                    if((how_bonded[i] == 1) ||
                       (how_bonded[i] == 7))
                            adjindxa[bondb] = adjindxa[bonda];
                    if((how_bonded[i] == 2) ||
                       (how_bonded[i] == 8) ||
                       (how_bonded[i] == 3))
                            adjindxa[bondb] = -adjindxa[bonda];
               }                    /* end if(!adjindxa[bondb] && adjindxa[bonda]) */
          }                         /* end for(i = 0; i < num_bonds; i++) */
     }                              /* end while(num > 0) */

     for(i = 0; i < num_atoms; i++)
          if(adjindxa[i] == -1) {
               if(!strncmp(assigned_type[i], "cc", 2))
                    assigned_type[i][1] = 'd';
               if(!strncmp(assigned_type[i], "ce", 2))
                    assigned_type[i][1] = 'f';
               if(!strncmp(assigned_type[i], "cg", 2))
                    assigned_type[i][1] = 'h';
               if(!strncmp(assigned_type[i], "pc", 2))
                    assigned_type[i][1] = 'd';
               if(!strncmp(assigned_type[i], "pe", 2))
                    assigned_type[i][1] = 'f';
               if(!strncmp(assigned_type[i], "nc", 2))
                    assigned_type[i][1] = 'd';
               if(!strncmp(assigned_type[i], "ne", 2))
                    assigned_type[i][1] = 'f';
          }                          /* end if(adjindxa[i] == -1) */

     free_me(adjindxa, FARGS);
     free_me(adjindxb, FARGS);

     return;
}

/* adjust gaff type cp to be cq in some bonding situations
   parameters are:
        num_atoms                    number of atoms
        num_bonds                    number of bonds
        assigned_type                the types as assigned by type judging
*/
void adjust_type_cp(int num_atoms, int num_bonds, char **assigned_type)
{
     int i, indx, num, flag, bonda, bondb;

     i = num_atoms + 10;
     adjindxa = (short *)allocator(i, sizeof(short), FARGS);
     adjindxb = (short *)allocator(i, sizeof(short), FARGS);

     for(i = 0; i < num_atoms; i++) {
          adjindxa[i] = 0;
          adjindxb[i] = 0;
     }

     indx = 0;
     num = 0;

     for(i = 0; i < num_atoms; i++)
          if(!strncmp(assigned_type[i], "cp", 2)) {
               adjindxb[i] = 1;
               if(!indx) {
                    adjindxa[i] = 1;
                    indx = 1;
               }
               num++;
          }

     if(!num) {
          free_me(adjindxa, FARGS);
          free_me(adjindxb, FARGS);
          return;
     }

     num--;

     while(num > 0) {
          num--;
          for(i = 0; i < num_bonds; i++) {
               bonda = bond_i[i];
               bondb = bond_j[i];
               if((adjindxb[bonda] + adjindxb[bondb]) != 2) continue;
               if(!adjindxa[bonda] && adjindxa[bondb]) {
                    if(how_bonded[i] == 1)
                         adjindxa[bonda] = adjindxa[bondb];
                    else
                         adjindxa[bonda] = -adjindxa[bondb];
               }                    /* end if(!adjindxa[bonda] && adjindxa[bondb]) */

               if(!adjindxa[bondb] && adjindxa[bonda]) {
                    if(how_bonded[i] == 1)
                         adjindxa[bonda] = adjindxa[bondb];
                    else
                         adjindxa[bonda] = -adjindxa[bondb];
               }                    /* end if(!adjindxa[bondb] && adjindxa[bonda]) */
          }                         /* end for(i = 0; i < num_bonds; i++) */
     }                              /* end while(num > 0) */

     for(i = 0; i < num_atoms; i++)
          if((adjindxa[i] == -1) && !strncmp(assigned_type[i], "cp", 2))
               assigned_type[i][1] = 'q';

     free_me(adjindxa, FARGS);
     free_me(adjindxb, FARGS);

     return;
}

/* check for gaff typing errors for types cc, ce, cg, pc, pe, nc, ne,
   cd, cf, ch, pd, pf, nd, nf, cp, and cq
   parameters are:
        num_atoms                    number of atoms
        num_bonds                    number of bonds
        assigned_type                the types as assigned by type judging
*/
void check_type_errors(int num_atoms, int num_bonds, char **assigned_type)
{
     int bondi, bondj, i;

     i = num_atoms + 10;
     adjindxa = (short *)allocator(i, sizeof(short), FARGS);

     for(i = 0; i < num_atoms; i++) {
          adjindxa[i] = 0;
          if(!strncmp(assigned_type[i], "cc", 2) ||
             !strncmp(assigned_type[i], "ce", 2) ||
             !strncmp(assigned_type[i], "nc", 2) ||
             !strncmp(assigned_type[i], "ne", 2) ||
             !strncmp(assigned_type[i], "pc", 2) ||
             !strncmp(assigned_type[i], "pe", 2) ||
             !strncmp(assigned_type[i], "cg", 2))  adjindxa[i] = 1;

          if(!strncmp(assigned_type[i], "cd", 2) ||
             !strncmp(assigned_type[i], "cf", 2) ||
             !strncmp(assigned_type[i], "nd", 2) ||
             !strncmp(assigned_type[i], "nf", 2) ||
             !strncmp(assigned_type[i], "pd", 2) ||
             !strncmp(assigned_type[i], "pf", 2) ||
             !strncmp(assigned_type[i], "ch", 2))  adjindxa[i] = -1;
     }

     for(i = 0; i < num_bonds; i++) {
          bondi = bond_i[i];
          bondj = bond_j[i];
          if((adjindxa[bondi] == 1) && (adjindxa[bondj] == 1) ||
             (adjindxa[bondi] == -1) && (adjindxa[bondj] == -1))
                  if((how_bonded[i] != 1) && (how_bonded[i] != 7))
                       warn_me(bondi, bondj, assigned_type);

          if(adjindxa[bondi]*adjindxa[bondj] == -1)
               if((how_bonded[i] == 1) || (how_bonded[i] == 7))
                    warn_me(bondi, bondj, assigned_type);

          if((!strncmp(assigned_type[bondi], "cp", 2) &&
              !strncmp(assigned_type[bondj], "cp", 2)) ||
             (!strncmp(assigned_type[bondi], "cq", 2) &&
              !strncmp(assigned_type[bondj], "cq", 2)))
                  if(how_bonded[i] != 1)
                       warn_me(bondi, bondj, assigned_type);

          if((!strncmp(assigned_type[bondi], "cp", 2) &&
              !strncmp(assigned_type[bondj], "cq", 2)) ||
             (!strncmp(assigned_type[bondi], "cq", 2) &&
              !strncmp(assigned_type[bondj], "cp", 2)))
                  if(how_bonded[i] == 1)
                       warn_me(bondi, bondj, assigned_type);
     }

     free_me(adjindxa, FARGS);

     return;
}

/* print warning of type errors
   parameters are:
        id1                          id of first atom in bad pair
        id2                          id of second atom in bad pair
        assigned_type                the types as assigned by type judging
*/
void warn_me(int id1, int id2, char **assigned_type)
{
     printf("\nWARNING: atom type of %5s (%s) and %5s (%s) may be wrong\n",
            atom_name[id1], assigned_type[id1], atom_name[id2], assigned_type[id2]);

     return;
}

/* Set proper ion types
   parameters are:
        num_atoms                    number of atoms
        assigned_type                the types as assigned by type judging
*/
void adjust_ion_types(int num_atoms, char **assigned_type)
{
     int i, j;
     char dummy[4];

     for(i = 0; i < num_atoms; i++) {
          if(bond_count[i]) continue;                       /* ions are not bonded */
          if(strlen(assigned_type[i]) > 2) continue;        /* only adjust 2 character types */
          j = fabs(atom_charge[i]);
          if((j > 4) || (j == 0)) continue;
          switch(j) {
               case 1:  if(atom_charge[i] < 0.0) dummy[0] = '-';
                        else dummy[0] = '+';
                        dummy[1] = '\0';
                        break;

               case 2:  dummy[0] = '2';
                        if(atom_charge[i] < 0.0) dummy[1] = '-';
                        else dummy[1] = '+';
                        dummy[2] = '\0';
                        break;

               case 3:  dummy[0] = '3';
                        if(atom_charge[i] < 0.0) dummy[1] = '-';
                        else dummy[1] = '+';
                        dummy[2] = '\0';
                        break;

               case 4:  dummy[0] = '4';
                        if(atom_charge[i] < 0.0) dummy[1] = '-';
                        else dummy[1] = '+';
                        dummy[2] = '\0';
                        break;
          }
          strcat(assigned_type[i], dummy);
     }

     return;
}
