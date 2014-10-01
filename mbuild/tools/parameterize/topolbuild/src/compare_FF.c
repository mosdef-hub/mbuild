/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        force field comparison routines compare_FF.c

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
#include "compare_FF.h"

/* comparison of bond to force field
   parameters are:
        atom types of each atom in bond:                  bond_atom1
                                                          bond_atom2
        number of bond parameter sets in force field:     bind_cntFF
        position of bond in molecule bonds table:         where
*/
int compare_bonds(char *bond_atom1, char *bond_atom2, int bind_cntFF,
                  int where)
{
     int j;

     for(j = 0; j < bind_cntFF; j++) {
          if(!((!atom_cmp(bond_atom1, bond_name1[j])   &&
                !atom_cmp(bond_atom2, bond_name2[j]))  ||
               (!atom_cmp(bond_atom1, bond_name2[j])   &&
                !atom_cmp(bond_atom2, bond_name1[j])))) continue;

          a_bond_type[where] = FF_bond_type[j];
          switch(a_bond_type[where]) {
               case 5:  a_bond_force[where]  = 0.0;
                        a_bond_length[where] = 0.0;
                        return(j + 1);

               case 7:
               case 8:
               case 9:
               case 6:
               case 2:
               case 1:  a_bond_force[where]  = bond_force[j];
                        a_bond_length[where] = bond_length[j];
                        return(j + 1);

               case 3:
               case 4:  a_bond_force[where]  = bond_force[j];
                        a_bond_length[where] = bond_length[j];
                        a_bond_beta[where] = bond_beta[j];
                        return(j + 1);

               default: sprintf(mess,
                           "Bond type error, force field bond number %d.\n",
                                j);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
     }

     return(0);
}

/* comparison of angle to force field
   parameters are:
        atom types of each atom in angle:                 angle_atom1
                                                          angle_atom2
                                                          angle_atom3
        number of angle parameter sets in force field:    angle_cntFF
        position of angle in molecule angles table:       where
        angle determination method setting                method
*/
int angles_comp(char *angle_atom1, char *angle_atom2, char *angle_atom3,
                int angle_cntFF, int where, int method)
{
     int j, k, l;

     for(j = 0; j < angle_cntFF; j++) {
          if(atom_cmp(angle_atom2, angle_name2[j]))
               continue;                   /* central atom must always match */

          if(!((!atom_cmp(angle_atom1, angle_name1[j])  &&
                !atom_cmp(angle_atom3, angle_name3[j])) ||
               (!atom_cmp(angle_atom1, angle_name3[j])  &&
                !atom_cmp(angle_atom3, angle_name1[j])))) continue;

          a_angle_type[where] = FF_angle_type[j];
          a_meth[where] = method + 1;              /* signal found in table */
          switch(a_angle_type[where]) {
               case 1:
               case 2:  a_angle_force[where] = angle_force[j];
                        a_angle_angle[where] = angle_angle[j];
                        return(j + 1);

               case 3:  for(k = 0; k < 2; k++)
                             a_angle_quartic[k][where] = angle_quartic[k][j];
                        return(j + 1);

               case 4:  for(k = 0; k < 4; k++)
                             a_angle_quartic[k][where] = angle_quartic[k][j];
                        return(j + 1);

               case 5:  a_angle_angle[where] = angle_angle[j];
                        for(k = 0; k < 2; k++)
                             a_angle_quartic[k][where] = angle_quartic[k][j];
                        return(j + 1);

               case 6:  a_angle_angle[where] = angle_angle[j];
                        for(k = 0; k < 5; k++)
                             a_angle_quartic[k][where] = angle_quartic[k][j];
                        return(j + 1);

               case 8:  a_angle_quartic[0][where] = angle_quartic[0][j];
                        a_angle_quartic[1][where] = angle_quartic[1][j];
                        return(j + 1);

               default: sprintf(mess,
                           "Angle type error, force field angle number %d.\n",
                                j);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
     }

     return(0);
}

/* comparison of angle to force field
   parameters are:
        atom types of each atom in angle:                 angle_atom1
                                                          angle_atom2
                                                          angle_atom3
        number of angle parameter sets in force field:    angle_cntFF
        angle force returned:                             force
        angle value returned:                             angle
        method of determination flag returned:            method
*/
int compare_angles(char *angle_atom1, char *angle_atom2, char *angle_atom3,
                   int angle_cntFF, float *force, float *angle, int *method)
{
     int j, k;

     k = *method;
     for(j = 0; j < angle_cntFF; j++) {
          if(atom_cmp(angle_atom2, angle_name2[j]))
               continue;              /* central atom must always match */

          if(!((!atom_cmp(angle_atom1, angle_name1[j])  &&
                !atom_cmp(angle_atom3, angle_name3[j])) ||
               (!atom_cmp(angle_atom1, angle_name3[j])  &&
                !atom_cmp(angle_atom3, angle_name1[j])))) continue;

          if(FF_angle_type[j] != 1) continue;
          *force = angle_force[j];
          *angle = angle_angle[j];
          *method = k + 1;                    /* signal found in table */
          return(1);
     }

     return(0);
}

/* empirical method to calculate angle phase and force 
   parameters are:
        atom index numbers of angle atoms:                ang1, ang2, ang3
        computed angle force returned:                    force
        computed angle value returned:                    angle
        molecule bond lengths to use already set:         bnd_length
        similar atoms index pointer:                      sim_atom_indx
        similar atoms list pointers:                      corr_name1
                                                          corr_name2
                                                          corr_name3
        number of angle parameter sets in force field:    angle_cntFF
        assigned atom types to use:                       assigned_type
*/
void calc_angle(int ang1, int ang2, int ang3, float *force, float *angle,
                float *bnd_length, int *method, int *sim_atom_indx,
                char **corr_name1, char **corr_name2, char **corr_name3,
                int angle_cntFF, char **assigned_type)
{
     int i, k, l, brkr;
     double length1, length2, my_angle, my_force, param[4];
     float angle1, angle2, use_force;

     brkr = 0;
     while(!brkr) {              /* device to allow me to skip things without */
                                 /* a really nasty nesting of if's */
          if((brkr = compare_angles(assigned_type[ang1], assigned_type[ang2],
                                    assigned_type[ang1], angle_cntFF,
                                    &use_force, &angle1, &k))) continue;

          if(corr_name1 == NULL) break;      /* flag lack of alternate names */

          if(sim_atom_indx[ang1] > 1)
               if((brkr = compare_angles(corr_name2[ang1], assigned_type[ang2],
                                         corr_name2[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if(sim_atom_indx[ang2] > 1)
               if((brkr = compare_angles(assigned_type[ang1], corr_name2[ang2],
                                         assigned_type[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if((sim_atom_indx[ang1] > 1) && (sim_atom_indx[ang2] > 1))
               if((brkr = compare_angles(corr_name2[ang1], corr_name2[ang2],
                                         corr_name2[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if(sim_atom_indx[ang1] > 2)
               if((brkr = compare_angles(corr_name3[ang1], assigned_type[ang2],
                                         corr_name3[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if(sim_atom_indx[ang2] > 2)
               if((brkr = compare_angles(assigned_type[ang1], corr_name3[ang2],
                                         assigned_type[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if((sim_atom_indx[ang1] > 2) && (sim_atom_indx[ang2] > 2))
               if((brkr = compare_angles(corr_name3[ang1], corr_name3[ang2],
                                         corr_name3[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if(sim_atom_indx[ang1] > 0)
               if((brkr = compare_angles(corr_name1[ang1], assigned_type[ang2],
                                         corr_name1[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if(sim_atom_indx[ang2] > 0)
               if((brkr = compare_angles(assigned_type[ang1], corr_name1[ang2],
                                         assigned_type[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          if((sim_atom_indx[ang1] > 0) && (sim_atom_indx[ang2] > 1))
               if((brkr = compare_angles(corr_name1[ang1], corr_name1[ang2],
                                         corr_name1[ang1], angle_cntFF,
                                         &use_force, &angle1, &k))) continue;

          break;                     /* failure, exit loop */
     }

     if(!brkr) return;               /* failed to calculate angle parameters */

     brkr = 0;
     while(!brkr) {                  /* device to allow me to skip things without */
                                     /* a really nasty nesting of if's */
          if((brkr = compare_angles(assigned_type[ang3], assigned_type[ang2],
                                    assigned_type[ang3], angle_cntFF,
                                    &use_force, &angle2, &k))) continue;

          if(corr_name1 == NULL) break;    /* flag lack of alternate names */

          if(sim_atom_indx[ang3] > 1)
               if((brkr = compare_angles(corr_name2[ang3], assigned_type[ang2],
                                         corr_name2[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if(sim_atom_indx[ang2] > 1)
               if((brkr = compare_angles(assigned_type[ang3], corr_name2[ang2],
                                         assigned_type[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if((sim_atom_indx[ang3] > 1) && (sim_atom_indx[ang2] > 1))
               if((brkr = compare_angles(corr_name2[ang3], corr_name2[ang2],
                                         corr_name2[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if(sim_atom_indx[ang3] > 2)
               if((brkr = compare_angles(corr_name3[ang3], assigned_type[ang2],
                                         corr_name3[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if(sim_atom_indx[ang2] > 2)
               if((brkr = compare_angles(assigned_type[ang3], corr_name3[ang2],
                                         assigned_type[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if((sim_atom_indx[ang3] > 2) && (sim_atom_indx[ang2] > 2))
               if((brkr = compare_angles(corr_name3[ang3], corr_name3[ang2],
                                         corr_name3[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if(sim_atom_indx[ang3] > 0)
               if((brkr = compare_angles(corr_name1[ang3], assigned_type[ang2],
                                         corr_name1[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if(sim_atom_indx[ang2] > 0)
               if((brkr = compare_angles(assigned_type[ang3], corr_name1[ang2],
                                         assigned_type[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          if((sim_atom_indx[ang3] > 0) && (sim_atom_indx[ang2] > 1))
               if((brkr = compare_angles(corr_name1[ang3], corr_name1[ang2],
                                         corr_name1[ang3], angle_cntFF,
                                         &use_force, &angle2, &k))) continue;

          break;                 /* failure, exit loop */
     }

     if(!brkr) return;           /* failed to calculate angle parameters */

     my_angle = 0.5 * (angle1 + angle2);

     brkr = 0;                  /* retrieve lengths and forces of both bonds */
     k = 0;
     l = 0;
     for(i = 0; i < bond_count[ang2]; i++) {
          if(!k)
               if((bond_i[(bonds_list[ang2][i])] == ang1)  ||
                  (bond_j[(bonds_list[ang2][i])] == ang1)) {
                         length1 = bnd_length[(bonds_list[ang2][i])];
                         brkr++;
                         k = 1;
               }

          if(!l)
               if((bond_i[(bonds_list[ang2][i])] == ang3)  ||
                  (bond_j[(bonds_list[ang2][i])] == ang3)) {
                         length2 = bnd_length[(bonds_list[ang2][i])];
                         brkr++;
                         l = 1;
          }
          if(brkr == 2) break;
     }

     if(brkr != 2)
          return;            /* something wrong. failed to find the bonds! */
     if((length1 > 9999.0) || (length2 > 9999.0))
          return;                     /* something not set */

     param[1] = (length1 - length2) * (length1 - length2);
     param[1] = param[1] / ((length1 + length2) * (length1 + length2));

     switch(atom_atno[ang2]) {
          case  6: param[2] = 1.339;
                   break;

          case  7: param[2] = 1.300;
                   break;

          case  8: param[2] = 1.249;
                   break;

          case 15: param[2] = 1.448;
                   break;

          case 16: param[2] = 0.906;
                   break;

          default: param[2] = 0.0;
                   break;
     }

     param[3] = force_param(atom_atno[ang1]);
     param[4] = force_param(atom_atno[ang3]);

     my_force = 143.9 * param[3] * param[2] * param[4] *
                     exp(-2.0 * param[1]) / (length1 + length2);
     my_force /= sqrt(my_angle * 3.1415926 / 180.0);
     *force = my_force;
     *angle = my_angle;
     *method = 9;                          /* signal empirical calculation */

     return;
}

/* based on atomic number, return a parameter to compute an angle force
   input parameter is:
        atomic number         atomic
*/
double force_param(int atomic)
{
     switch(atomic) {
          case  1: return(0.784);
                   break;

          case  6: return(1.183);
                   break;

          case  7: return(1.212);
                   break;

          case  8: return(1.219);
                   break;

          case  9: return(1.166);
                   break;

          case 15: return(1.280);
                   break;

          case 16: return(1.620);
                   break;

          case 17: return(1.272);
                   break;

          case 35: return(1.378);
                   break;

          case 53: return(1.398);
                   break;

          default: return(0.0);
                   break;
     }

     return(0.0);                           /* should never get here */
}

/* comparison of torsion to force field
   parameters are:
        atom types of each atom in torsion:               tors_atom1
                                                          tors_atom2
                                                          tors_atom3
                                                          tors_atom4
        number of torsion parameter sets in force field:  tors_cntFF
        torsion multiplier returned:                      mult
        torsion force returned:                           force
        torsion angle returned:                           phase
        torsion term returned:                            term
*/
int compare_tors(char *tors_atom1, char *tors_atom2, char *tors_atom3,
                 char *tors_atom4, int tors_cntFF, int *mult,
                 float *force, float *phase, float *term)
{
     int j;

/* First try for match to a definite specification */
     for(j = 0; j < tors_cntFF; j++) {
          if((!atom_cmp(tors_atom2, tors_name2[j]) &&
              !atom_cmp(tors_atom3, tors_name3[j])))          /* two center must match */
                    if((!atom_cmp(tors_atom1,
                                  tors_name1[j])) &&          /* 1 goes to 2 and 4 goes to 3 */
                       (!atom_cmp(tors_atom4,
                                  tors_name4[j]))) {          /* order here is important */
                              *mult  = tors_mult[j];
                              *force = tors_force[j];
                              *phase = tors_phase[j];
                              *term  = tors_term[j];
                              return(j + 1);
                    }

          if((!atom_cmp(tors_atom2, tors_name3[j]) &&
              !atom_cmp(tors_atom3, tors_name2[j])))          /* two center must match */
                    if((!atom_cmp(tors_atom1,
                                  tors_name4[j])) &&          /* 4 goes to 2 and 1 goes to 3 */
                       (!atom_cmp(tors_atom4,
                                  tors_name1[j]))) {          /* order here is important */
                              *mult  = tors_mult[j];
                              *force = tors_force[j];
                              *phase = tors_phase[j];
                              *term  = tors_term[j];
                              return(j + 1);
                    }
     }

/* Didn't work, so try for those that have wild cards */
     for(j = 0; j < tors_cntFF; j++) {
          if((!atom_cmp(tors_atom2, tors_name2[j]) &&
              !atom_cmp(tors_atom3, tors_name3[j])))          /* two center must match */
                    if((!atom_cmp(tors_atom1,
                                  tors_name1[j]) ||           /* 1 goes to 2 and 4 goes to 3 */
                        (tors_name1[j][0] == 'X')) &&
                       (!atom_cmp(tors_atom4, tors_name4[j]) ||
                        (tors_name4[j][0] == 'X'))) {         /* order here is important */
                              *mult  = tors_mult[j];
                              *force = tors_force[j];
                              *phase = tors_phase[j];
                              *term  = tors_term[j];
                              return(j + 1);
                    }

          if((!atom_cmp(tors_atom2, tors_name3[j]) &&
              !atom_cmp(tors_atom3, tors_name2[j])))          /* two center must match */
                    if((!atom_cmp(tors_atom1,
                                  tors_name4[j]) ||           /* 4 goes to 2 and 1 goes to 3 */
                        (tors_name4[j][0] == 'X')) &&
                       (!atom_cmp(tors_atom4,
                                  tors_name1[j]) ||
                        (tors_name1[j][0] == 'X'))) {         /* order here is important */
                              *mult  = tors_mult[j];
                              *force = tors_force[j];
                              *phase = tors_phase[j];
                              *term  = tors_term[j];
                              return(j + 1);
                    }
     }

     return(0);                                            /* nothing matched */
}

/* comparison of lexigraphically ordered improper to force field
   parameters are:
        atom types of each atom in improper:               impro_atom1
                                                           impro_atom2
                                                           impro_atom3
                                                           impro_atom4
        number of improper parameter sets in force field:  impro_cntFF
        improper multiplier returned:                      mult
        improper force returned:                           force
        improper angle returned:                           phase
        improper term returned:                            term
*/
int compare_improp(char *impro_atom1, char *impro_atom2, char *impro_atom3,
                   char *impro_atom4, int impro_cntFF, int *mult,
                   float *force, float *phase, float *term)
{
     int j, ibrk;
     char tmp1[8], tmp2[8], tmp3[8], tmp4[8];

     ibrk = -1;
     for(j = 0; j < impro_cntFF; j++) {
          strcpy(tmp1, impro_name1[j]);
          strcpy(tmp2, impro_name2[j]);
          strcpy(tmp3, impro_name3[j]);
          strcpy(tmp4, impro_name4[j]);

          if(impro_xcount[j] > 0)
               if(impro_name3[j][0] == 'X')
                    strcpy(tmp3, impro_atom3);  /* handle the case of 'X' in the key atom */

          if(atom_cmp(impro_atom3, tmp3)) continue;

          if((!atom_cmp(impro_atom1, tmp1) || (tmp1[0] == 'X')) &&              /* 1 2 4 */
             (!atom_cmp(impro_atom2, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom4, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom1) etc. */

          if((!atom_cmp(impro_atom2, tmp1) || (tmp1[0] == 'X')) &&              /* 2 1 4 */
             (!atom_cmp(impro_atom1, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom4, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom2) etc. */

          if((!atom_cmp(impro_atom2, tmp1) || (tmp1[0] == 'X')) &&              /* 2 4 1 */
             (!atom_cmp(impro_atom4, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom1, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom2) etc. */

          if((!atom_cmp(impro_atom4, tmp1) || (tmp1[0] == 'X')) &&              /* 4 1 2 */
             (!atom_cmp(impro_atom1, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom2, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom4) etc. */

          if((!atom_cmp(impro_atom4, tmp1) || (tmp1[0] == 'X')) &&              /* 4 2 1 */
             (!atom_cmp(impro_atom2, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom1, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom4) etc. */

          if((!atom_cmp(impro_atom1, tmp1) || (tmp1[0] == 'X')) &&              /* 1 4 2 */
             (!atom_cmp(impro_atom4, tmp2) || (tmp2[0] == 'X')) &&
             (!atom_cmp(impro_atom2, tmp4) || (tmp4[0] == 'X'))) {
                    ibrk = j;
                    break;
          }                                 /* end if((!atom_cmp(tmp1, impro_atom1) etc. */
     }

     if(ibrk < 0) return(0);                            /* none found */

     *mult  = impro_mult[j];
     *force = impro_force[j];
     *phase = impro_phase[j];
     *term  = impro_term[j];
     return(1);
}

/* comparisons against macro_atoms
   parameters are:
        my_atom               the atom of potential interest
        at_type               the atom to compare against

   external data are:
        num_macro             number of macro names available to match
        macro_name            array of macro names to test
        macro_nmbr            array of number of atom names under each
                              macro name
        macro_types           table of atom names under each macro name
*/
int atom_cmp(char *my_atom, char *at_type)
{
     int i, j;

     if(num_macro)                                      /* any macro names? */
          for(i = 0; i < num_macro; i++)
               if(!strcmp(at_type, macro_name[i])) {    /* match a macro name? */

/* check for a universal macro name such as all c's or all n's */
                    if((macro_nmbr[i] == 1) &&
                       (macro_types[i][0][1] == '#'))    /* a universal macro name? */
                            if(my_atom[0] == macro_types[i][0][0])
                                 return(0);         /* return a universal types match */

/* otherwise, search for macro name atom type match */
                    for(j = 0; j < macro_nmbr[i]; j++)
                         if(!strcmp(my_atom, macro_types[i][j]))
                              return(0);

/* Atom type was a macro name, but this atom is not in that set of names */
                         return(1);
               }

/* either did not have any macro names, or atom type failed to match a
   macro name compare atom names directly
*/

     j = strcmp(my_atom, at_type);

     return(j);
}
