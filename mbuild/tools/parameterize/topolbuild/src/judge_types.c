/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to judge the atom types and to allocate storage
        for the molecule's atom type judgments 

   NOTICE: This is a derivative work based on study of routines
   found in the antechamber 1.27 program atomtype, version 1.0,
   dated October, 2001 by Junmei Wang, Department of Pharmaceutical
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
#include "judge_types.h"

/* allocate atom type bond information and chemical environment
   arrays
   parameters are:
        atom_num               number of atoms
        num_bonds              number of bonds
*/
void init_atom_types(int atom_num, int num_bonds)
{
     int i, j;

/* bond information arrays */

     j = atom_num + 10;
     sb = (int *)allocator(j, sizeof(int), FARGS);
     SB = (int *)allocator(j, sizeof(int), FARGS);
     db = (int *)allocator(j, sizeof(int), FARGS);
     DB = (int *)allocator(j, sizeof(int), FARGS);
     tb = (int *)allocator(j, sizeof(int), FARGS);
     TB = (int *)allocator(j, sizeof(int), FARGS);
     AB = (int *)allocator(j, sizeof(int), FARGS);
     DL = (int *)allocator(j, sizeof(int), FARGS);

/* end bond information arrays */

/* chemical environment arrays */

     apch = mat_alloc(MAXLEVELS, MAXENVSTR, FARGS);
     env_bonds = mat_alloc(MAXENVBOND, 10, FARGS);
     env_atom1 = mat_alloc(MAXENVBOND, 10, FARGS);
     env_atom2 = mat_alloc(MAXENVBOND, 10, FARGS);
     atom_chem = mat_alloc(j, 10, FARGS);
     envblockstr = mat_alloc(MAXLEVELS, MAXENVSTR, FARGS);
     env_atom_name = tabl_alloc(MAXCHAIN, MAXLEVELS, 10, FARGS);
     env_ap = tabl_alloc(MAXCHAIN, MAXLEVELS, MAXENVSTR, FARGS);
     env_strname = tabl_alloc(MAXCHAIN, MAXLEVELS, 10, FARGS);

     chem_indx = (int *)allocator(MAXCHAIN, sizeof(int), FARGS);
     env_len = (int *)allocator(MAXCHAIN, sizeof(int), FARGS);
     env_con = (int **)allocator(MAXCHAIN, sizeof(int*), FARGS);
     for(i = 0; i < MAXCHAIN; i++)
          env_con[i] = (int *)allocator(MAXLEVELS, sizeof(int), FARGS);

     env_index = (int **)allocator(MAXCHAIN, sizeof(int*), FARGS);
     for(i = 0; i < MAXCHAIN; i++)
          env_index[i] = (int *)allocator(MAXLEVELS, sizeof(int), FARGS);

     selchain = (int *)allocator(j, sizeof(int), FARGS);
     selindx = (int *)allocator(j, sizeof(int), FARGS);
     sel_chain_indx = (int *)allocator(MAXCHAIN, sizeof(int), FARGS);
     schain_id = (int *)allocator(MAXSCAN, sizeof(int), FARGS);
     schain_num = (int *)allocator(MAXSCAN, sizeof(int), FARGS);
     schain_atom = (int **)allocator(MAXSCAN, sizeof(int*), FARGS);
     for(i = 0; i < MAXSCAN; i++)
          schain_atom[i] = (int *)allocator(j, sizeof(int), FARGS);

/* end chemical environment arrays */

     return;
}

/* free atom type bond information and chemical environment
   arrays
   parameter is:
        atom_num               number of atoms
*/
void free_atom_types(int atom_num)
{
     int i, j;

/* bond information arrays */

     free_me(sb, FARGS);
     free_me(SB, FARGS);
     free_me(db, FARGS);
     free_me(DB, FARGS);
     free_me(tb, FARGS);
     free_me(TB, FARGS);
     free_me(AB, FARGS);
     free_me(DL, FARGS);

/* end bond information arrays */

/* chemical environment arrays */

     for(i = 0; i < MAXLEVELS; i++) {
          free_me(apch[i], FARGS);
          free_me(envblockstr[i], FARGS);
     }
     free_me(apch, FARGS);
     free_me(envblockstr, FARGS);

     for(i = 0; i < MAXENVBOND; i++) {
          free_me(env_bonds[i], FARGS);
          free_me(env_atom1[i], FARGS);
          free_me(env_atom2[i], FARGS);
     }
     free_me(env_bonds, FARGS);
     free_me(env_atom1, FARGS);
     free_me(env_atom2, FARGS);

     for(i = 0; i < (atom_num + 10); i++)
          free_me(atom_chem[i], FARGS);
     free_me(atom_chem, FARGS);

     for(i = 0; i < MAXCHAIN; i++) {
          for(j = 0; j < MAXLEVELS; j++) {
               free_me(env_atom_name[i][j], FARGS);
               free_me(env_ap[i][j], FARGS);
               free_me(env_strname[i][j], FARGS);
          }
          free_me(env_atom_name[i], FARGS);
          free_me(env_ap[i], FARGS);
          free_me(env_strname[i], FARGS);
          free_me(env_con[i], FARGS);
          free_me(env_index[i], FARGS);
     }
     free_me(env_atom_name, FARGS);
     free_me(env_ap, FARGS);
     free_me(env_strname, FARGS);
     free_me(env_con, FARGS);
     free_me(env_index, FARGS);

     free_me(selchain, FARGS);
     free_me(selindx, FARGS);
     free_me(sel_chain_indx, FARGS);
     free_me(schain_id, FARGS);
     free_me(schain_num, FARGS);

     for(i = 0; i < MAXSCAN; i++)
          free_me(schain_atom[i], FARGS);
     free_me(schain_atom, FARGS);

/* end chemical environment arrays */

     return;
}

/* basic judgment of atom types
   parameters are:
        num_types          number of atom type definitions available
        num_wild           number of wild atom definitions available
        atom_num           number of atoms in molecule
        num_bonds          number of bonds in molecule
        assigned_type      array into which type assignments are placed
*/
void judge_atom_types(int num_types, int num_wild, int atom_num, int num_bonds,
                      char **assigned_type)
{
     int i, j, k, indx_brk, resindex, temp, temp2, bond_index;
     char resid[6], *c, *d;
     char mymask = 0x07;
     float temflt;

     initial = 0;

     for(i = 0; i < atom_num; i++) {
          assigned_type[i][0] = NULLCHAR;
          for(j = 0; j < num_types; j++) {
               resindex = 0;
               indx_brk = 0;
               switch(at_type_residue[j]) {
                    case 4:  strncpy(resid, atom_residue[i], 5);
                             resid[5] = NULLCHAR;
                             if((c = strpbrk(resid, "0123456789")) != NULL) {
                                  k = strlen(resid);
                                  if((c != &resid[0]) &&
                                     ((d = strpbrk(&resid[k-1], "0123456789")) != NULL))
                                            *c = NULLCHAR;
                             }

                             if(!strcmp(resid, residue_name[j]))
                                  resindex = 1;
                             break;

                    case -1: resindex = 1;
                             break;

                    case 0:  strcpy(assigned_type[i], at_type_name[j]);
                             if(strlen(assigned_type[i]) == 1) {
                                  assigned_type[i][1] = ' ';
                                  assigned_type[i][2] = NULLCHAR;
                             }
                             indx_brk = 1;
                             break;

                    case 1:  resindex = check_AA(i);
                             break;

                    case 2:  resindex = check_NA(i);
                             break;

                    case 3:  resindex = check_BIO(i);
                             break;

                    default: break; 
               }

               if(indx_brk) break;
               if(!resindex) continue;

               if(!at_type_no[j]){
                    if(atom_atno[i] > 0 )
                         continue;             /* not a dummy or lone pair type */
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }
               if(at_type_no[j] > 0)
                    if(at_type_no[j] != atom_atno[i])
                         continue;

/* charged names definitions check */
               if(use_charge &&
                  ((k = strlen(at_type_name[j]) - 1) > 1)) {
                         if(at_type_name[j][k] == '+') {
                              temflt = 1.0;
                              if(isdigit(at_type_name[j][(k - 1)]))
                                   temflt = (float)(at_type_name[j][(k - 1)] & mymask);
                                   if(atom_charge[i] != temflt)
                                        continue;
                         }
                         else
                              if((at_type_name[j][(k - 1)] == '+') &&
                                 (isdigit(at_type_name[j][k]))) {
                                        temflt = (float)(at_type_name[j][k] & mymask);
                                        if(atom_charge[i] != temflt)
                                             continue;
                              }
               }
/* end charged names definitions check */

               if(at_type_attached[j] == 999) {
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }
               if(at_type_attached[j] > -1) {              /* check number of attached atoms */
                    bond_index = bond_count[i];
                    if(strcmp(atom_type[i], "O.co2") &&
                       strcmp(atom_type[i], "O.2")   &&
                       !(atom_atno[i] == 6)          &&
                       !((atom_atno[i] == 7)         &&
                         (bond_count[i] > 3))) {
                              if((bond_index > 0)        &&
                                 (atom_atno[i] != 1)     &&
                                 ((atom_charge[i] < 0.8) ||
                                  (atom_charge[i] > -0.8))) {       /* correct for charges */
                                        if((atom_charge[i]/adj_chg) > 0.501) bond_index--;
                                        if((atom_charge[i]/adj_chg) < -0.501) bond_index++;
                              }
                    }

                    for(k = 0; k < bond_count[i]; k++)     /* correct for LP */
                         if(!atom_atno[(connect[i][k])])
                              bond_index--;

                    if(at_type_attached[j] != bond_index)
                         continue;
               }

               if(at_type_attached_H[j] == 999) {
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }

               if(at_type_attached_H[j] > -1) {            /* check number of attached H */
                    temp2 = 0;
                    for(k = 0; k < 6; k++) {
                         if((temp = connect[i][k]) < 0) break;
                         if(atom_atno[temp] == 1)
                              temp2++;
                    }
                    if(strcmp(atom_type[i], "O.co2") &&
                       strcmp(atom_type[i], "O.2")   &&
                       !(atom_atno[i] == 6)          &&
                       !((atom_atno[i] == 7)         &&
                         (bond_count[i] > 3))) {
                              if((bond_count[i] > 0) &&
                                 (atom_atno[i] != 1)     &&
                                 ((atom_charge[i] < 0.8) ||
                                  (atom_charge[i] > -0.8))) {         /* correct for charges */
                                        if((atom_charge[i]/adj_chg) > 0.501) temp2--;
                                        if((atom_charge[i]/adj_chg) < -0.501) temp2++;
                              }                            /* Somewhere loose in the solvent */
                    }                                      /* There's an H+ */
                    if(at_type_attached_H[j] != temp2)
                         continue;
               }               /* end if(at_type_attached_H[j] > -1) */

               if(at_type_ewd[j] == 999) {
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }

               if(at_type_ewd[j] > -1) {         /* for H, check number of ewd attached */
                    temp = connect[i][0];        /* to attached atom */
                    temp2 = 0;
                    for(k = 0; k < 6; k++) {
                         if(connect[temp][k] < 0) break;
                         if(atom_ewd[(connect[temp][k])])
                              temp2++;
                    }
                    if(at_type_ewd[j] != temp2) continue;
               }               /* end if(at_type_ewd[j] > -1) */

               if(at_type_prop[j][0] == NULLCHAR) {
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }

               if(strcmp(at_type_prop[j], "*"))            /* check atom properties */
                    if(!prop_check(i, -1, num_bonds, at_type_prop[j]))
                         continue;

               if(at_type_env[j][0] == NULLCHAR) {
                    strcpy(assigned_type[i], at_type_name[j]);
                    if(strlen(assigned_type[i]) == 1) {
                         assigned_type[i][1] = ' ';
                         assigned_type[i][2] = NULLCHAR;
                    }
                    break;
               }

               if(strcmp(at_type_env[j], "*")) {              /* check chemical */
                    bond_index = 1;                           /* environment string */
                    if(at_env_bonds[j][0] == NULLCHAR)
                         bond_index = 0;
                    initial = i;
                    if(chem_environ(initial, at_type_env[j], at_env_bonds[j],
                                    bond_index, atom_num, num_bonds)) {
                         strcpy(assigned_type[i], at_type_name[j]);
                         if(strlen(assigned_type[i]) == 1) {
                              assigned_type[i][1] = ' ';
                              assigned_type[i][2] = NULLCHAR;
                         }
                         break;
                    } else continue;                     /* not strictly necessary */
               }
          }                    /* end for(j = 0; j < num_types; j++) */
          if(assigned_type[i][0] == NULLCHAR) {               /* unassigned? */
               strncpy(assigned_type[i], atom_symbl[i], 2);   /* set type as symbol */
               assigned_type[i][2] = NULLCHAR;                /* with ?? appended */
               strcat(assigned_type[i], "??");
          }                    /* end if(assigned_type[i][0] == NULLCHAR) */
     }                         /* end for(i = 0; i < atom_num; i++) */

     return;
}

/* check residue name against amino acid names and return 1 if match or 0 if not
   parameter is:
        number            number of the atom to check against
*/
int check_AA(int number)
{
     static char amino_acids[26][4] = {"ALA", "GLY", "SER", "THR", "LEU",
                                       "ILE", "VAL", "ASN", "GLN", "ARG",
                                       "HID", "HIE", "HIP", "TRP", "PHE",
                                       "TYR", "GLU", "ASP", "LYS", "PRO",
                                       "CYS", "CYX", "MET", "ASH", "GLH",
                                       "LYH" };
#define NUM_AAS 26

     int i, k;
     char resid[6], *c, *d;

     strncpy(resid, atom_residue[number], 5);
     resid[5] = NULLCHAR;
     if((c = strpbrk(resid, "0123456789")) != NULL) {
          k = strlen(resid);
          if((c != &resid[0]) &&
             ((d = strpbrk(&resid[k-1], "0123456789")) != NULL))
                    *c = NULLCHAR;
     }

     for(i = 0; i < NUM_AAS; i++) {
          if(!strcmp(resid, amino_acids[i]))
               return(1);
     }

     return(0);
}

/* check residue name against nucleic acid names and return 1 if match or 0 if not
   parameter is:
        number            number of the atom to check against
*/
int check_NA(int number)
{
     static char nucleic_acids[8][5] = {"DADE", "RADE", "DTRY", "RURA", "DGUA",
                                        "RGUA", "DCYT", "RCYT" };

#define NUM_NAS 8

     int i, k;
     char resid[6], *c, *d;

     strncpy(resid, atom_residue[number], 5);
     resid[5] = NULLCHAR;
     if((c = strpbrk(resid, "0123456789")) != NULL) {
          k = strlen(resid);
          if((c != &resid[0]) &&
             ((d = strpbrk(&resid[k-1], "0123456789")) != NULL))
                    *c = NULLCHAR;
     }

     for(i = 0; i < NUM_NAS; i++) {
          if(!strcmp(resid, nucleic_acids[i]))
               return(1);
     }

     return(0);
}

/* check residue name against both amino acid names and nucleic acid names and
   return 1 if match or 0 if not
   parameter is:
        number            number of the atom to check against
*/
int check_BIO(int number)
{
     if((check_AA(number)) || (check_NA(number))) return(1);

     return(0);
}

/* set bond information arrays
   parameters are:
        atom_num               number of atoms
        num_bonds              number of bonds
*/
void bond_info(int atom_num, int num_bonds)
{
     int i;

     for(i = 0; i < atom_num; i++) {
          sb[i] = 0;
          SB[i] = 0;
          db[i] = 0;
          DB[i] = 0;
          tb[i] = 0;
          TB[i] = 0;
          AB[i] = 0;
          DL[i] = 0;
     }

     for(i = 0; i < num_bonds; i++) {
          if(!atom_atno[(bond_i[i])] ||
             !atom_atno[(bond_j[i])])
                    continue;                          /* don't count bonds to lone pairs */
          switch(how_bonded[i]) {
               case  1:  sb[(bond_i[i])]++;
                         sb[(bond_j[i])]++;
                         SB[(bond_i[i])]++;
                         SB[(bond_j[i])]++;
                         break;

               case  2:  db[(bond_i[i])]++;
                         db[(bond_j[i])]++;
                         DB[(bond_i[i])]++;
                         DB[(bond_j[i])]++;
                         break;

               case  3:  tb[(bond_i[i])]++;
                         tb[(bond_j[i])]++;
                         TB[(bond_i[i])]++;
                         TB[(bond_j[i])]++;
                         break;

               case  7:  sb[(bond_i[i])]++;
                         sb[(bond_j[i])]++;
                         AB[(bond_i[i])]++;
                         AB[(bond_j[i])]++;
                         break;

               case  8:  db[(bond_i[i])]++;
                         db[(bond_j[i])]++;
                         AB[(bond_i[i])]++;
                         AB[(bond_j[i])]++;
                         break;

               case  9:  sb[(bond_i[i])]++;
                         sb[(bond_j[i])]++;
                         SB[(bond_i[i])]++;
                         SB[(bond_j[i])]++;
                         DL[(bond_i[i])]++;
                         DL[(bond_j[i])]++;
                         break;

               case 10:  AB[(bond_i[i])]++;
                         AB[(bond_j[i])]++;
                         break;

               default:  break;
          }                    /* end switch(how_bonded[i]) */
     }                         /* end for(i = 0; i < num_bonds; i++) */

/* correct bonding information for charges on some atoms */

     for(i = 0; i < atom_num; i++) {
          if(!bond_count[i]) continue;              /* skip ions because they don't have bonds */
          if(!strcmp(atom_type[i], "O.co2") ||
             !strcmp(atom_type[i], "O.2")   ||
             (atom_atno[i] == 1)            ||
             (atom_atno[i] == 6)            ||
             ((atom_atno[i] == 7)           &&
              (bond_count[i] > 3)))
                    continue;                       /* skip carbonyls, charged but no H+ */
                                                    /* skip N with 4 bonds, H+ and charge, */
                                                    /* but has special type */
          if((atom_charge[i] > 0.8) ||
             (atom_charge[i] < -0.8)) continue;   /* skip large charges */
          if((atom_charge[i]/adj_chg) > 0.501) {
               if(sb[i]) sb[i]--;              /* positive charge subtracts one if have */
               if(SB[i]) SB[i]--;              /* single bond */
               continue;
          }
          if((atom_charge[i]/adj_chg) < -0.501) {  /* Somewhere floating in solvent,
                                                       there's an H+ */
               sb[i]++;                            /* negative charge adds one */
               SB[i]++;                            /* single bond */
               continue;
          }
     }

     return;
}

/* check on atom properties
   parameters are:
        num                atom number in molecule to check
        mark               environment property check parameter
        num_bonds          number of bonds in the molecule
        prop_str           atom properties string from atom type definition
*/
int prop_check(int num, int mark, int num_bonds, char *prop_str)
{
     int i, j, index1, index2, featnum, nmbr;
     int isw, my_rg;
     char tmpstr[5], tmp1[5], tmp2[5];
     char *feat_str[NUM_FEAT] = { "RG3", "RG4", "RG5", "RG6", "RG7",
                                  "RG8", "RG9", "RG10", "RG", "NR",
                                  "AR1", "AR2", "AR3", "AR4", "AR5",
                                  "SB",  "sb",  "DB",  "db",  "TB",
                                  "tb",  "DL",  "AB",  "PL",  "MI" };
     int featlen[NUM_FEAT] = { 3, 3, 3, 3, 3, 3, 3, 4, 2, 2,
                               3, 3, 3, 3, 3, 2, 2, 2, 2, 2,
                               2, 2, 2, 2, 2 };

     nmbr = 0;
     index2 = 0;

     for(i = 0; i < strlen(prop_str); i++) {
          if(prop_str[i] =='[') continue;
          if((prop_str[i] == '.') ||
             (prop_str[i] == ',') ||
             (prop_str[i] == ']')) {
                  index1 = 1;
                  if((prop_str[i] == '.') &&
                     (index2 == 0))
                            index2 = 1;
                  if((prop_str[i] == '.') &&
                     (index2 == 2))
                            continue;
                  if((prop_str[i] == '.') &&
                     (index2 == -1))
                            index2 = 1;
                  if((prop_str[i] == ',') &&
                     (index2 == 2)) {
                            index2 = 0;
                            nmbr = 0;
                            continue;
                  }
                  if((prop_str[i] == ',') &&
                     (index2 == 1))
                            return(0);
                  if(prop_str[i] == ']')
                       index2 = 0;
                  tmpstr[nmbr] = NULLCHAR;
                  featnum = -1;
                  if(isdigit(tmpstr[0])) {
                       tmp1[0] = tmpstr[0];
                       tmp1[1] = NULLCHAR;
                       featnum = atoi(tmp1);
                       for(j = 0; j < nmbr; j++)
                            tmpstr[j] = tmpstr[(j + 1)];
                  }

                  for(j = strlen(tmpstr); j < 5; j++)
                       tmpstr[j] = NULLCHAR;
                  nmbr = 0;

                  isw = -1;
                  for(j = 0; j < NUM_FEAT; j++)
                       if(!strncmp(tmpstr, feat_str[j], featlen[j])) {
                            isw = j;
                            break;
                       }

                  switch(isw) {
/* "RG3"  */
                       case 0:  check_feat(featnum, arom_atomno[num][3],
                                           &index1);
                                break;
/* "RG4"  */
                       case 1:  check_feat(featnum, arom_atomno[num][4],
                                           &index1);
                                break;
/* "RG5"  */
                       case 2:  check_feat(featnum, arom_atomno[num][5],
                                           &index1);
                                break;
/* "RG6"  */
                       case 3:  check_feat(featnum, arom_atomno[num][6],
                                           &index1);
                                break;
/* "RG7"  */
                       case 4:  check_feat(featnum, arom_atomno[num][7],
                                           &index1);
                                break;
/* "RG8"  */
                       case 5:  check_feat(featnum, arom_atomno[num][8],
                                           &index1);
                                break;
/* "RG9"  */
                       case 6:  check_feat(featnum, arom_atomno[num][9],
                                           &index1);
                                break;
/* "RG10" */
                       case 7:  check_feat(featnum, arom_atomno[num][10],
                                           &index1);
                                break;
/* "RG"   */
                       case 8:  check_feat(featnum, arom_atomno[num][0],
                                           &index1);
                                break;
/* "NR" */
                       case 9:  check_feat(featnum, arom_num[num],
                                           &index1);
                                break;
/* "AR1" */
                       case 10: check_feat(featnum, ar_set[num][1],
                                           &index1);
                                break;
/* "AR2" */
                       case 11: check_feat(featnum, ar_set[num][2],
                                           &index1);
                                break;
/* "AR3" */
                       case 12: check_feat(featnum, ar_set[num][3],
                                           &index1);
                                break;
/* "AR4" */
                       case 13: check_feat(featnum, ar_set[num][4],
                                           &index1);
                                break;
/* "AR5" */
                       case 14: check_feat(featnum, ar_set[num][5],
                                           &index1);
                                break;
/* "SB" */
                       case 15: bond_check(tmpstr, featnum, SB[num], &index1,
                                           mark, num_bonds, num, -1);
                                break;

/* "sb" */
                       case 16: bond_check(tmpstr, featnum, sb[num], &index1,
                                           mark, num_bonds, num, 1);
                                break;

/* "DB" */
                       case 17: bond_check(tmpstr, featnum, DB[num], &index1,
                                           mark, num_bonds, num, -2);
                                break;

/* "db" */
                       case 18: bond_check(tmpstr, featnum, db[num], &index1,
                                           mark, num_bonds, num, 2);
                                break;

/* "TB" */
                       case 19: bond_check(tmpstr, featnum, TB[num], &index1,
                                           mark, num_bonds, num, -3);
                                break;

/* "tb" */
                       case 20: bond_check(tmpstr, featnum, tb[num], &index1,
                                           mark, num_bonds, num, 3);
                                break;

/* "DL" */
                       case 21: bond_check(tmpstr, featnum, DL[num], &index1,
                                           mark, num_bonds, num, 9);
                                break;

/* "AB" */
                       case 22: bond_check(tmpstr, featnum, AB[num], &index1,
                                           mark, num_bonds, num, 4);
                                break;
/* "PL" */
                       case 23: if(!use_charge ||
                                   ((atom_atno[num] !=  6) &&
                                    (atom_atno[num] !=  7) &&
                                    (atom_atno[num] !=  8) &&
                                    (atom_atno[num] != 15) &&
                                    (atom_atno[num] != 16))) {
                                          index1 = 0;
                                          break;
                                }
                                if(atom_charge[num] < 0.85)
                                     index1 = 0;
                                break;
/* "MI" */
                       case 24: if(!use_charge ||
                                   ((atom_atno[num] !=  6) &&
                                    (atom_atno[num] !=  7) &&
                                    (atom_atno[num] !=  8) &&
                                    (atom_atno[num] != 15) &&
                                    (atom_atno[num] != 16))) {
                                          index1 = 0;
                                          break;
                                }
                                if(atom_charge[num] > -0.85)
                                     index1 = 0;
                       default: break;
                  }                      /* end switch(isw) */

                  if((index1 == 0) &&
                     (index2 == 0))
                            return(0);
                  if((index1 != 0) &&
                     (index2 == 1))
                            index2 = 2;
                  if((index1 == 0) &&
                     (index2 == 1))
                            index2 = -1;
                  continue;
          }                              /* end if(prop_str[i] == '.')  etc. */
          if(index2 == 2)
               continue;
          tmpstr[nmbr++] = prop_str[i];
          if(nmbr > 5) {
               my_fatal(FARGS,
                        "Property string element too long.\n");
               exit(1);
          }
     }                              /* end for(i = 0; i < strlen(prop_str); i++) */

     return(index1);
}

/* generalized check for RG and AR properties
   parameters are:
        featindx           flag to determine if further chexking needed
        status             the actual property of the atom in question
        setindx            return setting from the check
*/
void check_feat(int featindx, int status, int *setindx)
{
     if(featindx == -1) {
          if(!status)
               *setindx = 0;
     }
     else {
          if(featindx != status)
               *setindx = 0;
     }
     return;
}

/* generalized bond features check
   parameters are:
        thestr             a bond type property to check
        featindx           flag to determine if further checking needed
        status             the actual property of the atom in question
        setindx            return setting from the check
        mark               atom number for must(/must not) have this bond to 
        num_bonds          number of bonds in the molecule
        num                atom number in molecule to pass to further checking
        bndtyp             bond type must(/must not) have to atom mark
*/
void bond_check(char *thestr, int featindx, int status, int *setindx,
                int mark, int num_bonds, int num, int bndtyp)
{
     int indexa;

     indexa = 0;
     if(strlen(thestr) > 2)
          if(thestr[2] == '\'')
               indexa = 1;
     if(strlen(thestr) > 3)
          if(indexa && (thestr[3] == '\''))
               indexa = 2;
     if(mark < 0)
          indexa = 0;
     if(!indexa) {
          if(featindx == -1) {
               if(!status)
                    *setindx = 0;
          }
          else {
               if(featindx != status)
                    *setindx = 0;
          }
     }
     if((indexa == 1) && (mark != -1) &&
        !junctbond(num, mark, num_bonds, bndtyp))
             *setindx = 0;
     if((indexa == 2) && (mark != -1) &&
        junctbond(num, mark, num_bonds, bndtyp))
             *setindx = 0;

     return;
}

/* bond type further check
   parameters are:
        num                atom number in molecule to check
        mark               atom number for must(/must not) have this bond to
        num_bonds          number of bonds in the molecule
        bndtyp             a bond type to check against
*/
int junctbond(int num, int mark, int num_bonds, int bndtyp)
{
     int i, j;

     j = 0;
     for(i = 0; i < num_bonds; i++) {
          if(((bond_i[i] == num) && (bond_j[i] == mark)) ||
             ((bond_j[i] == num) && (bond_i[i] == mark))) {
                  switch(bndtyp) {
                       case  1: if((how_bonded[i] == 1) ||
                                   (how_bonded[i] == 7) ||
                                   (how_bonded[i] == 9) ||
                                   (how_bonded[i] == 10)) j = 1;
                                break;

                       case -1: if((how_bonded[i] == 1) ||
                                   (how_bonded[i] == 9)) j = 1;
                                break;

                       case  2: if((how_bonded[i] == 2) ||
                                   (how_bonded[i] == 8) ||
                                   (how_bonded[i] == 9) ||
                                   (how_bonded[i] == 10)) j = 1;
                                break;

                       case -2: if((how_bonded[i] == 2) ||
                                   (how_bonded[i] == 9)) j = 1;
                                break;

                       case -3:
                       case  3: if(how_bonded[i] == 3) j = 1;
                                break;

                       case  4: if((how_bonded[i] == 7) ||
                                   (how_bonded[i] == 8) ||
                                   (how_bonded[i] == 10)) j = 1;
                                break;

                       case  9: if(how_bonded[i] == 9) j = 1;
                                break;

                       default: break;
                  }            /* end switch(bndtyp) */

                  break;
          }                    /* end if(((bond_i[i] == num)  etc. */
     }                         /* end for(i = 0; i < num_bonds; i++) */

     return(j);
}

/* check chemical environments
   parameters are:
        num                atom number in molecule to check
        env_str            chemical environment from atom type definitions
        bonds_str          environment bonds string from atom type
                           definitions
        bndindx            flag for existence of valid bonds string
        atom_num           number of atoms in the molecule
        num_bonds          number of bonds in the molecule
*/
int chem_environ(int num, char *env_str, char *bonds_str, int bndindx,
                 int atom_num, int num_bonds)
{
     int i, j, level, t_aps, env_name, apcntr, env_num, chain, indx0;
     char t_envstr[MAXENVSTR], t_apch[MAXENVSTR], t_ccs[MAXENVSTR];
     char t_atom[MAXLEVELS][10], tmpc[10];
     int t_index[MAXLEVELS], t_con[MAXLEVELS];
     int env_cntr, env_indx, my_index, maxchain, success;

     if(strlen(env_str) >= (4 * MAXENVSTR)) {
          my_fatal(FARGS,
                   "Environment string length exceeds maximum.\n");
          exit(1);
     }

     if(strlen(bonds_str) >= MAXENVSTR) {
          my_fatal(FARGS,
                   "Bonds string length exceeds maximum.\n");
          exit(1);
     }

/* check that chemical environment string has balanced sets of brackets */
     balanced_str('(',')', env_str);
     balanced_str('[',']', env_str);
     balanced_str('<','>', env_str);

     indx0 = 0;
     level = 0;
     chain = 0;
     t_aps = 0;
     env_name = 0;
     for(i = 0; i < MAXCHAIN; i++)
          for(j = 0; j < MAXLEVELS; j++) {
               env_con[i][j] = 0;
               env_atom_name[i][j][0] = NULLCHAR;
               env_ap[i][j][0] = NULLCHAR;
          }

     for(j = 0; j < MAXLEVELS; j++) {
          t_con[j] = 0;
          t_index[j] = 0;
          apch[j][0] = NULLCHAR;
          envblockstr[j][0] = NULLCHAR;
          t_atom[j][0] = NULLCHAR;
     }

/* parse the environment string */
     for(i = 0; i < strlen(env_str); i++) {
          if((!t_aps && !env_name) &&
             (isalpha(env_str[i]) != 0) &&
             (isalpha(env_str[(i + 1)]) != 0)) continue;
          if(env_str[i] == '(') {
               level++;
               if(level >= MAXLEVELS) {
                    my_fatal(FARGS,
                             "Number of levels exceeds maximum.\n");
                    exit(1);
               }
          }

          if(env_str[i] == ')') {
               level--;
               if(level < 0) {
                    my_fatal(FARGS,
                       "Too many \")\" preceeding their balancing \"(\".\n");
                    exit(1);
               }
          }

          if(!t_aps && (env_str[i] == '[')) {
               t_aps = 1;
               t_apch[0] = '[';
               t_apch[1] = NULLCHAR;
               apcntr = 1;
               continue;
          }

          if(t_aps && (env_str[i] == ']')) {
               t_index[level] = 1;
               t_apch[apcntr++] = ']';
               t_apch[apcntr] = NULLCHAR;
               strcpy(apch[level], t_apch);
               t_aps = 0;
               continue;
          }

          if(t_aps) {
               t_apch[apcntr++] = env_str[i];
               if(apcntr > (MAXENVSTR - 3)) {
                    my_fatal(FARGS, 
                             "Property string length exceeds maximum.\n");
                    exit(1);
               }
               continue;
          }

          if(!env_name && (env_str[i] == '<')) {
               env_name = 1;
               t_envstr[0] = NULLCHAR;
               env_num = 0;
               continue;
          }

          if(env_name && (env_str[i] == '>')) {
               t_envstr[env_num] = NULLCHAR;
               strcpy(envblockstr[level], t_envstr);
               env_name = 0;
               continue;
          }

          if(env_name) {
               t_envstr[env_num++] = env_str[i];
               if(env_num > (MAXENVSTR - 3)) {
                    my_fatal(FARGS, 
                             "Environment block length exceeds maximum.\n");
                    exit(1);
               }
               continue;
          }

          if((env_str[i] == ',') && (env_str[(i - 1)] != ')' )) {
               if(level >= MAXLEVELS) {
                    my_fatal(FARGS, "Number of levels exceeds maximum.\n");
                    exit(1);
               }

               for(j = 0; j < level; j++) {     /* level is always 1 or greater here */
                    strcpy(&env_atom_name[chain][j][0], t_atom[(j + 1)]);
                    env_con[chain][j] = t_con[(j + 1)];
                    env_index[chain][j] = t_index[(j + 1)];
                    strcpy(&env_ap[chain][j][0], apch[(j + 1)]);
                    strcpy(&env_strname[chain][j][0], envblockstr[(j + 1)]);
               }
               env_len[chain] = level;
               chain++;
               if(chain >= MAXCHAIN) {
                    my_fatal(FARGS, "Number of chains exceeds maximum.\n");
                    exit(1);
               }
          }

          if((env_str[i] == ')') &&
             (env_str[(i - 1)] != ')' )) {
               for(j = 0; j < (level + 1); j++) {   /* level had 1 subtracted above */
                    strcpy(&env_atom_name[chain][j][0], t_atom[(j + 1)]);
                    env_con[chain][j] = t_con[(j + 1)];
                    env_index[chain][j] = t_index[(j + 1)];
                    strcpy(&env_ap[chain][j][0], apch[(j + 1)]);
                    strcpy(&env_strname[chain][j][0], envblockstr[(j + 1)]);
               }
               env_len[chain] = level + 1;
               chain++;
               if(chain >= MAXCHAIN) {
                    my_fatal(FARGS, "Number of chains exceeds maximum.\n");
                    exit(1);
               }
          }

          if(isalpha(env_str[i]) &&
             isalpha(env_str[(i + 1)]))
                    continue;

          if(isalpha(env_str[i])) {
               indx0 = 1;
               if(isalpha(env_str[i - 1])) {
                    t_atom[level][0] = env_str[i - 1];
                    t_atom[level][1] = env_str[i];
                    t_atom[level][2] = NULLCHAR;
               }
               else {
                    t_atom[level][0] = env_str[i];
                    t_atom[level][1] = NULLCHAR;
               }
               envblockstr[level][0] = NULLCHAR;
               apch[level][0] = NULLCHAR;
               t_index[level] = 0;
          }

          if(isdigit(env_str[i])) {
               tmpc[0] = env_str[i];
               tmpc[1] = NULLCHAR;
               t_con[level] = atoi(tmpc);
          }
          else if(indx0) {
                    t_con[level] = 0;
                    indx0 = 0;
               }
     }                   /* end for(i = 0; i < strlen(env_str); i++) */

     env_bond_num = 0;
     if(bndindx) {
          env_cntr = 0;
          env_indx = 0;
          t_ccs[0] = NULLCHAR;
          for(i = 0; i < strlen(bonds_str); i++) {
               if((bonds_str[i] == ',') || (i == strlen(bonds_str))) {
                    t_ccs[env_cntr] = NULLCHAR;
                    strcpy(env_bonds[env_bond_num], t_ccs);
                    env_cntr = 0;
                    env_indx = 0;
                    t_ccs[0] = NULLCHAR;
                    env_bond_num++;
                    if(env_bond_num >= MAXENVBOND) {
                         my_fatal(FARGS,
                            "Number of connected atom pairs exceeds maximum\n");
                         exit(1);
                    }
                    continue;
               }

               if(bonds_str[i] == ':') {
                    t_ccs[env_cntr] = NULLCHAR;
                    if(!env_indx)
                         strcpy(env_atom1[env_bond_num], t_ccs);
                    if(env_indx)
                         strcpy(env_atom2[env_bond_num], t_ccs);
                    env_cntr = 0;
                    t_ccs[0] = NULLCHAR;
                    env_indx++;
                    continue;
               }

               t_ccs[env_cntr++] = bonds_str[i];
               if(env_cntr > (MAXENVSTR - 3)) {
                    my_fatal(FARGS, 
                       "Environment bond string length exceeds maximum.\n");
                    exit(1);
               }
          }              /* end for(i = 0; i < strlen(bonds_str); i++) */
     }                   /* end if(bndindx) */

     maxchain = -1;
     scan_num = 0;
     for(i = 0; i < chain; i++) {
          chem_indx[i] = 0;
          if(env_len[i] > maxchain)
               maxchain = env_len[i];
     }

     for(i = 0; i < atom_num; i++) {
          selindx[i] = -1;
          selchain[i] = -1;
          atom_chem[i][0] = NULLCHAR;
     }

     enviromatch(0, num, chain, maxchain, num_bonds);

     my_index = 1;
     for(i = 0; i < chain; i++)
          if(!chem_indx[i]) {
               my_index = 0;
               break;
          }

     if(my_index) {
          for(i = 0; i < chain; i++)
               sel_chain_indx[i] = -1;
          success = 0;
          dccheck(0, &success, chain, atom_num, num_bonds);
          my_index = success;
     }

     return(my_index);
}

/* check string for balanced characters
   parameters are
        cin               opening character that must be balanced
        cout              closing character that balances opening
                          character
        str               string to check for balanced opening
                          and closing
*/
void balanced_str(char cin, char cout, char *str)
{
     if(str_balanced(cin, cout, str)) {
          sprintf(mess, " %c and %c do not match for %s\n", cin, cout, str);
          my_fatal(FARGS, mess);
          exit(1);
     }

     return;
}

/* check string for balanced characters
   parameters are
        cin               opening character that must be balanced
        cout              closing character that balances opening character
        str               striing to check for balanced opening and closing
*/
int str_balanced(char cin, char cout, char *str)
{
     int i, countin, countout;

     countin = 0;
     countout = 0;

     for(i = 0; i < strlen(str); i++) {
          if(str[i] == cin) countin++;
          if(str[i] == cout) countout++;
          if(countout > countin)
               return(1);
     }

     i = countin - countout;

     return(i);
}

/* environmental matching
   parameters are:
        selection          position within parsed environment
        startnum           atom number in molecule to check
        chain              current chain in parsed environment
        maximum            longest environment sequence
        num_bonds          number of bonds in the molecule
*/
void enviromatch(int selection, int startnum, int chain,
                 int maximum, int num_bonds)
{
     int i, j, k, m, start, indx, indx1, timetobreak;
     int bond_index, ijk;

     start = -1;
     selchain[selection] = startnum;
     selection++;
     selindx[startnum] = 1;

     for(k = 0; k < chain; k++) {
          indx = 1;
          if((selection - 1) == env_len[k]) {
               for(j = 1; j < selection; j++) {
                    indx1 = 1;
                    bond_index = bond_count[(selchain[j])];
                    for(ijk = 0; ijk < bond_count[(selchain[j])]; ijk++)     /* correct for LP */
                         if(!atom_atno[(connect[(selchain[j])][ijk])])
                              bond_index--;
                    if(env_con[k][(j - 1)] &&
                       (bond_index != env_con[k][(j - 1)])) {
                            indx = 0;
                            break;
                    }

                    if(env_atom_name[k][(j-1)][0] == NULLCHAR) {
                         indx = 0;       /* if the atom name is NULLCHAR,
                         break;          /* then wildmatch would return 0 */
                    }
                    if(!strncmp(env_atom_name[k][(j-1)] ,"EW", 2))
                         if(atom_ewd[(selchain[j])] != 1) {
                              indx = 0;
                              break;
                         }

                    if(atom_name[(selchain[j])][0] != env_atom_name[k][(j-1)][0])
                         indx1 = 0;
                    if(indx1 && strlen(env_atom_name[k][(j-1)]) >= 2)
                         if(atom_name[(selchain[j])][1] != env_atom_name[k][(j-1)][1])
                              indx1 = 0;

                    if(!indx1)
                         if(!wildmatch(env_atom_name[k][(j-1)], selchain[j])) {
                              indx = 0;
                              break;
                         }

                    if(env_index[k][(j-1)]) {
                         if(strlen(env_ap[k][(j-1)])) {
                              if(j != 1)
                                   if(!prop_check(selchain[j], selchain[(j-1)],
                                             num_bonds, env_ap[k][(j-1)])) {
                                        indx = 0;
                                        break;
                                   }
                              if(j == 1)
                                   if(!prop_check(selchain[j], initial,
                                                  num_bonds, env_ap[k][(j-1)])) {
                                        indx = 0;
                                        break;
                                   }
                         }         /* end if(strlen(env_ap[k][(j-1)])) */
                    }              /* end if(env_index[k][(j-1)]) */
               }                   /* end for(j = 1; j < selection; j++) */

               if(indx) {
                    chem_indx[k]++;
                    schain_num[scan_num] = selection - 1;
                    schain_id[scan_num] = k;
                    for(m = 0; m < (selection - 1); m++)
                         schain_atom[scan_num][m] = selchain[(m + 1)];
                    scan_num++;
                    if(scan_num >= MAXSCAN) {
                         my_fatal(FARGS, "Maximum number of chains exceeded\n");
                         exit(1);
                    }
               }                   /* end if(indx) */
          }                        /* end if((selection - 1) == env_len[k]) */
     }                             /* end for(k = 0; k < chain; k++) */

     timetobreak = 0;
     for(i = 0; i < 6; i++) {
          if(connect[startnum][i] == -1) return;
          start = connect[startnum][i];
          for(j = 0; j < selection; j++)
               if(selchain[j] == start) {
                    timetobreak = 1;
                    break;
               }
          if(timetobreak) {
               timetobreak = 0;
               continue;
          }
          if(selection > maximum) return;
          enviromatch(selection, start, chain, maximum, num_bonds);
     }                             /* end for(i = 0; i < 6; i++) */

     return;
}

/* wild atom matching
   parameters are:
        my_atom               name in chemical environment to match against
                              wild atom name
        my_id                 atom to check against
*/
int wildmatch(char *my_atom, int my_id)
{
     int i, j, indxr, tobrk, k, bondings, bond_index, ijk;
     char atomid[5];

     if(!num_wilds) return(0);

     indxr = 0;
     tobrk = 0;
     for(i = 0; i < num_wilds; i++) {
          if(!strncmp(my_atom, wild_name[i], 2)) {
               for(j = 0; j < wild_count[i]; j ++) {
                    bondings = 0;
                    if((k = strcspn(&wild_elements[i][j][0], "0123456789")) <
                       strlen(&wild_elements[i][j][0])) {
                            strcpy(atomid, &wild_elements[i][j][k]);
                            bondings = atoi(atomid);
                    }

                    bond_index = bond_count[my_id];
                    for(ijk = 0; ijk < bond_count[my_id]; ijk++)     /* correct for LP */
                         if(!atom_atno[(connect[my_id][ijk])])
                              bond_index--;
                    if((atom_atno[my_id] == wild_atno[i][j]) &&
                       ((bondings  == 0) ||
                        (bond_index == bondings))) {
                            indxr = 1;
                            tobrk = 1;
                            break;
                    }
               }               /* end for(j = 0; j < wild_count[i]; j ++) */
          }                    /* end if(!strncmp(my_atom, wild_name, 2)) */
          if(tobrk) break;
     }                         /* end for(i = 0; i < num_wilds; i++) */

     return(indxr);
}

/* check of chains
   parameters are:
        startpoint         current point in chain of check
        success            flag successful return
        chain              current environment chain to check against
        num_atoms          number of atoms in the molecule
        num_bonds          number of bonds in the molecule
*/
void dccheck(int startpoint, int *success, int chain,  int num_atoms,
             int num_bonds)
{
     int i;

     if(*success) return;

     startpoint++;
     for(i = 0; i < scan_num; i++) {
          if(schain_id[i] != (startpoint - 1)) continue;
          sel_chain_indx[(startpoint - 1)] = i;

          if(startpoint == chain)
               if((*success = check_chain(chain, num_atoms, num_bonds)))
                    return;

          if(startpoint < chain)
               dccheck(startpoint, success, chain, num_atoms,
                       num_bonds);
     }

     return;
}

/* check chain
   parameters are:
        chain              current environment chain to check against
        num_atoms          number of atoms in the molecule
        num_bonds          number of bonds in the molecule
*/
int check_chain(int chain, int num_atoms, int num_bonds)
{
     int i, j, k, the_min, flag, tmpry1, tmpry2;

     for(i = 0; i < chain; i++)
          for(j = (i + 1); j < chain; j++) {
               tmpry1 = sel_chain_indx[i];
               tmpry2 = sel_chain_indx[j];
               if(tmpry1 == tmpry2) return(0);
               the_min = schain_num[tmpry2];
               if(schain_num[tmpry1] < schain_num[tmpry2])
                    the_min = schain_num[tmpry1];

               flag = 0;
               for(k = 0; k < the_min; k++)
                    if(schain_atom[tmpry1][k] != schain_atom[tmpry2][k]) {
                         flag = k + 1;
                         break;
                    }

               if(!flag) return(0);

               for(k = (flag - 1); k < the_min; k++)
                    if(schain_atom[tmpry1][k] == schain_atom[tmpry2][k]) {
                         if(env_strname[i][k][0] == NULLCHAR)
                              return(0);

                         if(strcmp(env_strname[i][k], env_strname[j][k]))
                              return(0);
                    }
          }                    /* end for(j = (i + 1); j < chain; j++) */

     if(env_bond_num) {
          for(i = 0; i < num_atoms; i++)
               atom_chem[i][0] = NULLCHAR;
          for(i = 0; i < chain; i++)
               for(j = 0; j < scan_num; j++)
                    if(sel_chain_indx[i] == j)
                         for(k = 0; k < env_len[i]; k++)
                              strcpy(atom_chem[(schain_atom[j][k])],
                                     env_strname[i][k]);

          strcpy(atom_chem[initial], "sa");
          for(i = 0; i < env_bond_num; i++) {
               tmpry1 = -1;
               tmpry2 = -1;

               for(j = 0; j < num_atoms; j++) {
                    if(!strcmp(atom_chem[j], env_atom1[i]))
                         tmpry1 = j;
                    if(!strcmp(atom_chem[j], env_atom2[i]))
                         tmpry2 = j;
               }               /* end for(j = 0; j < num_atoms; j++) */

               if((tmpry1 != -1) && (tmpry1 == tmpry2))
                    return(0);
               if(!check_bonds(tmpry1, tmpry2, env_bonds[i], num_bonds))
                    return(0);
          }                    /* end for(i = 0; i < env_bond_num; i++) */
     }                         /* end if(env_bond_num) */

     return(1);
}

/* check environment bonds settings.  returns 1 for a match
   and 0 on failure
   parameters are:
        ida                first atom number in bond
        idb                second atom number in bond
        bondstr            string of the bond type to check
        num_bonds          number of bonds in the molecule
*/
int check_bonds(int ida, int idb, char *bondstr, int num_bonds)
{
     int i, index;

     index = 0;

     for(i = 0; i < num_bonds; i++)
          if(((bond_i[i] == ida) && (bond_j[i] == idb)) ||
             ((bond_j[i] == ida) && (bond_i[i] == idb))) {
                  if(!strcmp(bondstr, "any") || !strcmp(bondstr, "ANY") ||
                     !strcmp(bondstr, "Any")) {
                          index = 1;
                          break;
                  }

                  switch(how_bonded[i]) {
                       case 1:  if(!strcmp(bondstr, "sb")) {
                                     index = 1;
                                     break;
                                }

                       case 2:  if(!strcmp(bondstr, "db")) {
                                     index = 1;
                                     break;
                                }

                       case 3:  if(!strcmp(bondstr, "tb")) {
                                     index = 1;
                                     break;
                                }

                       case 7:  if(!strcmp(bondstr, "sb") ||
                                   !strcmp(bondstr, "SB") ||
                                   !strcmp(bondstr, "AB") ||
                                   !strcmp(bondstr, "ab")) {
                                        index = 1;
                                        break;
                                }

                       case 8:  if(!strcmp(bondstr, "db") ||
                                   !strcmp(bondstr, "DB") ||
                                   !strcmp(bondstr, "AB") ||
                                   !strcmp(bondstr, "ab")) {
                                        index = 1;
                                        break;
                                }

                       case 9:  if(!strcmp(bondstr, "sb") ||
                                   !strcmp(bondstr, "SB") ||
                                   !strcmp(bondstr, "db") ||
                                   !strcmp(bondstr, "DB")) {
                                        index = 1;
                                        break;
                                }

                       case 10: if(!strcmp(bondstr, "AB") ||
                                   !strcmp(bondstr, "ab")) {
                                        index = 1;
                                        break;
                                }

                       default: break;
                  }                 /* end switch(how_bonded[i]) */

                  if(index) break;
          }                         /* end if(((bond_i[i] == ida) etc. */

     return(index);
}
