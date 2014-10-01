/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to use the gaff force field data to set parameters use_gaff.c

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
#include "use_gaff.h"
#include "compare_FF.h"
#include "similars.h"
#include "mol2.h"
#include "amber_FF.h"
#include "multors.h"
#include "block_memory.h"

/* set atom parameters for gaff from force field read in.
   parameters are:
        number of atoms in molecule                  num_atoms
        number of atom types in force field          num_mol
*/
void set_gaff_atomFF(int num_atoms, int num_mol)
{
     int i, j, gotit;

     for(i = 0; i < num_atoms; i++) {
          gotit = 0;
          for(j = 0; j < num_mol; j++) {
               if(strcmp(a_ff_type[i],
                         ffatom_name[j])) continue;      /* this atom type? */
               gotit = 1;
               a_pol[i] = ffatom_pol[j];
               break;
          }                             /* end for(j = 0; j < num_mol; j++) */

          if(gotit) continue;                     /* skip rest if have assignment */
          if(gaff_corr_num[i] < 1) continue;      /* skip rest if not assigned and do not have */
                                                  /* any alternative names */
          for(j = 0; j < num_mol; j++) {          /* try first alternative name */
               if(strcmp(gaff_name_corr1[i],
                         ffatom_name[j])) continue;      /* anything yet? */
               gotit = 1;
               a_pol[i] = ffatom_pol[j];
               break;
          }                             /* end for(j = 0; j < num_mol; j++)  again */
     }                                  /* end for(i = 0; i < num_atoms; i++) for gaff */

     return;
}

/* set van der Waals atom parameters for gaff from force field read in.
   parameters are:
        number of atoms in molecule:                       num_atoms
        number of van der Walls types in force field:      num_vdw
*/
void set_gaff_vdwFF(int num_atoms, int num_vdw)
{
     int i, j, gotit;

     for(i = 0; i < num_atoms; i++) {
          gotit = 0;

          for(j = 0; j < num_vdw; j++) {
               if(strcmp(a_ff_type[i],
                         vdw_atom_name[j])) continue;
               gotit = 1;
               a_vdw_atom_radius[i] = vdw_atom_radius[j];
               a_vdw_atom_pot[i]    = vdw_atom_pot[j];
               break;
          }                             /* end for(j = 0; j < num_vdw; j++) */

          if(gotit) continue;                     /* skip rest if have assignment */
          if(gaff_corr_num[i] < 1) continue;      /* skip rest if not assigned and do not have */
                                                  /* any alternative names */
          for(j = 0; j < num_vdw; j++) {          /* try first alternative name */
               if(strcmp(gaff_name_corr1[i],
                         vdw_atom_name[j])) continue;      /* anything yet? */
               gotit = 1;
               a_vdw_atom_radius[i] = vdw_atom_radius[j];
               a_vdw_atom_pot[i]    = vdw_atom_pot[j];
               break;
          }                             /* end for(j = 0; j < num_vdw; j++) again */

          if(gotit) continue;                     /* skip rest if have assignment */
          if(gaff_corr_num[i] < 2) continue;      /* skip rest if not assigned and do not have */
                                                  /* any alternative names */
          for(j = 0; j < num_vdw; j++) {          /* try first alternative name */
               if(strcmp(gaff_name_corr2[i],
                         vdw_atom_name[j])) continue;      /* anything yet? */
               gotit = 1;
               a_vdw_atom_radius[i] = vdw_atom_radius[j];
               a_vdw_atom_pot[i]    = vdw_atom_pot[j];
               break;
          }                             /* end for(j = 0; j < num_vdw; j++) second */

          if(gotit) continue;                     /* skip rest if have assignment */
          if(gaff_corr_num[i] < 3) continue;      /* skip rest if not assigned and do not have */
                                                  /* any alternative names */
          for(j = 0; j < num_vdw; j++) {          /* try first alternative name */
               if(strcmp(gaff_name_corr3[i],
                         vdw_atom_name[j])) continue;      /* anything yet? */
               gotit = 1;
               a_vdw_atom_radius[i] = vdw_atom_radius[j];
               a_vdw_atom_pot[i]    = vdw_atom_pot[j];
               break;
          }                             /* end for(j = 0; j < num_vdw; j++) third */
     }
     return;
}

/* set bond parameters for gaff from force field read in.
   parameters are:
        number of bonds in molecule:              num_bonds
        number of bond types in force field:      bind_cntFF
*/
void set_gaff_bondFF(int num_bonds, int bind_cntFF)
{
     int i, bnd1, bnd2;

     for(i = 0; i < num_bonds; i++) {
          bnd1 = bond_i[i];
          bnd2 = bond_j[i];

          if(compare_bonds(a_ff_type[bnd1],
                           a_ff_type[bnd2], bind_cntFF, i))
                                continue;                       /* these atom types? */
/* if have alternate name on first atom, try it */
          if(gaff_corr_num[bnd1] > 1)
               if(compare_bonds(gaff_name_corr2[bnd1],
                                a_ff_type[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */
/* if have alternate name on second atom, try it */
          if(gaff_corr_num[bnd2] > 1)
               if(compare_bonds(a_ff_type[bnd1],
                                gaff_name_corr2[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */

          if((gaff_corr_num[bnd1] > 1) &&
             (gaff_corr_num[bnd2] > 1))
                    if(compare_bonds(gaff_name_corr2[bnd1],
                                     gaff_name_corr2[bnd2], bind_cntFF, i))
                                          continue;                       /* these atom types? */
/* if have another alternate name on first atom, try it */
          if(gaff_corr_num[bnd1] > 2)
               if(compare_bonds(gaff_name_corr3[bnd1],
                                a_ff_type[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */
/* if have another alternate name on second atom, try it */
          if(gaff_corr_num[bnd2] > 2)
               if(compare_bonds(a_ff_type[bnd1],
                                gaff_name_corr3[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */

          if((gaff_corr_num[bnd1] > 1) &&
             (gaff_corr_num[bnd2] > 2))
                    if(compare_bonds(gaff_name_corr2[bnd1],
                                     gaff_name_corr3[bnd2], bind_cntFF, i))
                                          continue;                       /* these atom types? */

          if((gaff_corr_num[bnd1] > 2) &&
             (gaff_corr_num[bnd2] > 1))
                    if(compare_bonds(gaff_name_corr3[bnd1],
                                     gaff_name_corr2[bnd2], bind_cntFF, i))
                                          continue;                       /* these atom types? */

          if((gaff_corr_num[bnd1] > 2) &&
             (gaff_corr_num[bnd2] > 2))
                    if(compare_bonds(gaff_name_corr3[bnd1],
                                     gaff_name_corr3[bnd2], bind_cntFF, i))
                                          continue;                       /* these atom types? */
/* if have third alternate name on first atom, try it */
          if(gaff_corr_num[bnd1] > 0)
               if(compare_bonds(gaff_name_corr1[bnd1],
                                a_ff_type[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */
/* if have third alternate name on second atom, try it */
          if(gaff_corr_num[bnd2] > 0)
               if(compare_bonds(a_ff_type[bnd1],
                                gaff_name_corr1[bnd2], bind_cntFF, i))
                                     continue;                       /* these atom types? */

          if((gaff_corr_num[bnd1] > 0) &&
             (gaff_corr_num[bnd2] > 0))
                    if(compare_bonds(gaff_name_corr1[bnd1],
                                     gaff_name_corr1[bnd2], bind_cntFF, i))
                                          continue;                       /* these atom types? */
     }                                  /* end for(i = 0; i < num_bonds; i++) */

     return;
}

/* set angle parameters for gaff from force field read in.
   parameters are:
        number of angles in molecule:              num_angles
        number of angle types in force field:      angle_cntFF
*/
void set_gaff_angleFF(int num_angles, int angle_cntFF)
{
     int i, ang1, ang2, ang3;

     for(i = 0; i < num_angles; i++) {
          ang1 = angle_tab[i][0];
          ang2 = angle_tab[i][1];
          ang3 = angle_tab[i][2];
          a_meth[i] = 0;
          if(compare_angles(a_ff_type[ang1], a_ff_type[ang2],
                      a_ff_type[ang3], angle_cntFF, &a_angle_force[i],
                      &a_angle_angle[i], &a_meth[i])) {
                           a_angle_type[i] = 1;
                           if(opposite[i] == 'y')
                                a_angle_angle[i] *= 2.0;
                           continue;
          }
/* if have alternate name on first atom, try it */
          a_meth[i] = 10010;
          if(gaff_corr_num[ang1] > 1)
               if(compare_angles(gaff_name_corr2[ang1],
                           a_ff_type[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have alternate name on second atom, try it */
          a_meth[i] = 10100;
          if(gaff_corr_num[ang2] > 1)
               if(compare_angles(a_ff_type[ang1],
                           gaff_name_corr2[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
 /* if have alternate name on third atom, try it */
          a_meth[i] = 11000;
          if(gaff_corr_num[ang3] > 1)
               if(compare_angles(a_ff_type[ang1],
                           a_ff_type[ang2], gaff_name_corr2[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have 3 alternate names, try them */
          a_meth[i] = 11110;
          if((gaff_corr_num[ang1] > 1) &&
             (gaff_corr_num[ang2] > 1) &&
             (gaff_corr_num[ang3] > 1))
                    if(compare_angles(gaff_name_corr2[ang1],
                                gaff_name_corr2[ang2], gaff_name_corr2[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have alternate name on first atom, try it */
          a_meth[i] = 20020;
          if(gaff_corr_num[ang1] > 2)
               if(compare_angles(gaff_name_corr3[ang1],
                           a_ff_type[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
 /* if have alternate name on second atom, try it */
          a_meth[i] = 20200;
          if(gaff_corr_num[ang2] > 2)
               if(compare_angles(a_ff_type[ang1],
                           gaff_name_corr3[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
 /* if have alternate name on third atom, try it */
          a_meth[i] = 22000;
          if(gaff_corr_num[ang3] > 2)
               if(compare_angles(a_ff_type[ang1],
                           a_ff_type[ang2], gaff_name_corr3[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have 3 alternate names, try them */
          a_meth[i] = 22210;
          if((gaff_corr_num[ang1] > 1) &&
             (gaff_corr_num[ang2] > 2) &&
             (gaff_corr_num[ang3] > 2))
                    if(compare_angles(gaff_name_corr2[ang1],
                                gaff_name_corr3[ang2], gaff_name_corr3[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 22120;
          if((gaff_corr_num[ang1] > 2) &&
             (gaff_corr_num[ang2] > 1) &&
             (gaff_corr_num[ang3] > 2))
                    if(compare_angles(gaff_name_corr3[ang1],
                                gaff_name_corr2[ang2], gaff_name_corr3[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 21220;
          if((gaff_corr_num[ang1] > 2) &&
             (gaff_corr_num[ang2] > 2) &&
             (gaff_corr_num[ang3] > 1))
                    if(compare_angles(gaff_name_corr3[ang1],
                                gaff_name_corr3[ang2], gaff_name_corr2[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 22110;
          if((gaff_corr_num[ang1] > 1) &&
             (gaff_corr_num[ang2] > 1) &&
             (gaff_corr_num[ang3] > 2))
                    if(compare_angles(gaff_name_corr2[ang1],
                                gaff_name_corr2[ang2], gaff_name_corr3[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 21210;
          if((gaff_corr_num[ang1] > 1) &&
             (gaff_corr_num[ang2] > 2) &&
             (gaff_corr_num[ang3] > 1))
                    if(compare_angles(gaff_name_corr2[ang1],
                                gaff_name_corr3[ang2], gaff_name_corr2[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 21120;
          if((gaff_corr_num[ang1] > 2) &&
             (gaff_corr_num[ang2] > 1) &&
             (gaff_corr_num[ang3] > 1))
                    if(compare_angles(gaff_name_corr3[ang1],
                                gaff_name_corr2[ang2], gaff_name_corr2[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have 3 alternate names, try them */
          a_meth[i] = 22220;
          if((gaff_corr_num[ang1] > 2) &&
             (gaff_corr_num[ang2] > 2) &&
             (gaff_corr_num[ang3] > 2))
                    if(compare_angles(gaff_name_corr3[ang1],
                                gaff_name_corr3[ang2], gaff_name_corr3[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }
/* if have alternate name on first atom, try it */
          a_meth[i] = 10;
          if(gaff_corr_num[ang1] > 0)
               if(compare_angles(gaff_name_corr1[ang1],
                           a_ff_type[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have alternate name on second atom, try it */
          a_meth[i] = 100;
          if(gaff_corr_num[ang2] > 0)
               if(compare_angles(a_ff_type[ang1],
                           gaff_name_corr1[ang2], a_ff_type[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
 /* if have alternate name on third atom, try it */
          a_meth[i] = 1000;
          if(gaff_corr_num[ang3] > 0)
               if(compare_angles(a_ff_type[ang1],
                           a_ff_type[ang2], gaff_name_corr1[ang3],
                           angle_cntFF, &a_angle_force[i],
                           &a_angle_angle[i], &a_meth[i])) {
                                a_angle_type[i] = 1;
                                if(opposite[i] == 'y')
                                     a_angle_angle[i] *= 2.0;
                                continue;
               }
/* if have 3 alternate names, try them */
          a_meth[i] = 1110;
          if((gaff_corr_num[ang1] > 0) &&
             (gaff_corr_num[ang2] > 0) &&
             (gaff_corr_num[ang3] > 0))
                    if(compare_angles(gaff_name_corr1[ang1],
                                gaff_name_corr1[ang2], gaff_name_corr1[ang3],
                                angle_cntFF, &a_angle_force[i],
                                &a_angle_angle[i], &a_meth[i])) {
                                     a_angle_type[i] = 1;
                                     if(opposite[i] == 'y')
                                          a_angle_angle[i] *= 2.0;
                                     continue;
                    }

/* if got here, failed to find an appropriate entry.  Try to compute from other
   entries, emperically */

          a_meth[i] = 0;
          calc_angle(ang1, ang2, ang3, &a_angle_force[i], &a_angle_angle[i],
                     a_bond_length, &a_meth[i], gaff_corr_num, gaff_name_corr1,
                     gaff_name_corr2, gaff_name_corr3, angle_cntFF, a_ff_type);
          if(a_meth[i]) {
               a_angle_type[i] = 1;
               if(opposite[i] == 'y')
                    a_angle_angle[i] *= 2.0;
          }
    }

     return;
}

/* set torsion parameters for gaff from force field read in.
   parameters are:
        number of torsions in molecule:              num_tors
        number of torsion types in force field:      tors_cntFF
*/
void set_gaff_torsFF(int num_tors, int tors_cntFF)
{
     int i, j, tors1, tors2, tors3, tors4;

     for(i = 0; i < num_tors; i++) {
          tors1 = torsion_tab[i][0];
          tors2 = torsion_tab[i][1];
          tors3 = torsion_tab[i][2];
          tors4 = torsion_tab[i][3];

          if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                               a_ff_type[tors3], a_ff_type[tors4],
                               tors_cntFF, &a_tors_mult[i],
                               &a_tors_force[i], &a_tors_phase[i],
                               &a_tors_term[i]))) {
                                    a_tors_type[i] = 1;
                                    check_more_tors(j, i);
                                    continue;
          }

          if(gaff_corr_num[tors1] > 1)
               if((j = compare_tors(gaff_name_corr2[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors2] > 1)
               if((j = compare_tors(a_ff_type[tors1], gaff_name_corr2[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors3] > 1)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    gaff_name_corr2[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors4] > 1)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], gaff_name_corr2[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if(gaff_corr_num[tors1] > 2)
               if((j = compare_tors(gaff_name_corr3[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors2] > 2)
               if((j = compare_tors(a_ff_type[tors1], gaff_name_corr3[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors3] > 2)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    gaff_name_corr3[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if(gaff_corr_num[tors4] > 2)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], gaff_name_corr3[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                         a_tors_type[i] = 1;
                                         check_more_tors(j, i);
                                         continue;
               }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 2))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr3[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 2) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr3[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 1) &&
             (gaff_corr_num[tors2] > 2) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr2[tors1], gaff_name_corr3[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 2) &&
             (gaff_corr_num[tors2] > 1) &&
             (gaff_corr_num[tors3] > 1) &&
             (gaff_corr_num[tors4] > 1))
                    if((j = compare_tors(gaff_name_corr3[tors1], gaff_name_corr2[tors2],
                                         gaff_name_corr2[tors3], gaff_name_corr2[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if(gaff_corr_num[tors1] > 0)
               if((j = compare_tors(gaff_name_corr1[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if(gaff_corr_num[tors2] > 0)
               if((j = compare_tors(a_ff_type[tors1], gaff_name_corr1[tors2],
                                    a_ff_type[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if(gaff_corr_num[tors3] > 0)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    gaff_name_corr1[tors3], a_ff_type[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if(gaff_corr_num[tors4] > 0)
               if((j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                                    a_ff_type[tors3], gaff_name_corr1[tors4],
                                    tors_cntFF, &a_tors_mult[i],
                                    &a_tors_force[i], &a_tors_phase[i],
                                    &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }

          if((gaff_corr_num[tors1] > 0) &&
             (gaff_corr_num[tors2] > 0) &&
             (gaff_corr_num[tors3] > 0) &&
             (gaff_corr_num[tors4] > 0))
                    if((j = compare_tors(gaff_name_corr1[tors1], gaff_name_corr1[tors2],
                                         gaff_name_corr1[tors3], gaff_name_corr1[tors4],
                                         tors_cntFF, &a_tors_mult[i],
                                         &a_tors_force[i], &a_tors_phase[i],
                                         &a_tors_term[i]))) {
                                              a_tors_type[i] = 1;
                                              check_more_tors(j, i);
                                              continue;
                    }
     }

    return;
}

/* if this has multiple torsions, add links and set them as well
   parameter:
        where               points to the next torsion entry which is the next
                            component of a multiple torsion and these facts are
                            used below.
        what                location in atoms record of the torsions
*/
void check_more_tors(int where, int what)
{
     int k;

     if(tors_more[(where - 1)]) {
          multors_cnt[what] = tors_more[(where - 1)];     /* set count */
          multiples[what] = (float **)allocator(multors_cnt[what], sizeof(float **),
                                                FARGS);
          for(k = 0; k < multors_cnt[what]; k++) {
               multiples[what][k] = (float *)allocator(4, sizeof(float), FARGS);
               multiples[what][k][0] = (float)tors_mult[where];
               multiples[what][k][1] = tors_force[where];
               multiples[what][k][2] = tors_phase[where];
               multiples[what][k][3] = tors_term[where];
               where++;
          }
     }

     return;
}

/* set improper parameters for gaff from force field read in.
   parameters are:
        number of improper in molecule:              num_impr
        number of improper types in force field:     impr_cntFF
*/
void set_gaff_improFF(int num_impr, int impr_cntFF)
{
     int i, impr1, impr2, impr3, impr4;

     if(num_impr)
          for(i = 0; i < num_impr; i++) {
               impr1 = improperid[i][0];
               impr2 = improperid[i][1];
               impr3 = improperid[i][2];
               impr4 = improperid[i][3];

               if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                 a_ff_type[impr3], a_ff_type[impr4],
                                 impr_cntFF, &a_impr_mult[i],
                                 &a_impr_force[i], &a_impr_phase[i],
                                 &a_impr_term[i])) continue;

               if(gaff_corr_num[impr3] > 0)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      gaff_name_corr1[impr3], a_ff_type[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if(gaff_corr_num[impr3] > 1)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      gaff_name_corr2[impr3], a_ff_type[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if(gaff_corr_num[impr3] > 2)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      gaff_name_corr3[impr3], a_ff_type[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if(gaff_corr_num[impr4] > 0)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      a_ff_type[impr3], gaff_name_corr1[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if(gaff_corr_num[impr4] > 1)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      a_ff_type[impr3], gaff_name_corr2[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if(gaff_corr_num[impr4] > 2)
                    if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                      a_ff_type[impr3], gaff_name_corr3[impr4],
                                      impr_cntFF, &a_impr_mult[i],
                                      &a_impr_force[i], &a_impr_phase[i],
                                      &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 0) &&
                  (gaff_corr_num[impr4] > 0))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr1[impr3], gaff_name_corr1[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 1) &&
                  (gaff_corr_num[impr4] > 0))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr2[impr3], gaff_name_corr1[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 0) &&
                  (gaff_corr_num[impr4] > 1))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr1[impr3], gaff_name_corr2[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 1) &&
                  (gaff_corr_num[impr4] > 1))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr2[impr3], gaff_name_corr2[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;


               if((gaff_corr_num[impr3] > 2) &&
                  (gaff_corr_num[impr4] > 0))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr3[impr3], gaff_name_corr1[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 2) &&
                  (gaff_corr_num[impr4] > 1))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr3[impr3], gaff_name_corr2[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 0) &&
                  (gaff_corr_num[impr4] > 2))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr1[impr3], gaff_name_corr3[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 1) &&
                  (gaff_corr_num[impr4] > 2))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr2[impr3], gaff_name_corr3[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

               if((gaff_corr_num[impr3] > 2) &&
                  (gaff_corr_num[impr4] > 2))
                         if(compare_improp(a_ff_type[impr1], a_ff_type[impr2],
                                           gaff_name_corr3[impr3], gaff_name_corr3[impr4],
                                           impr_cntFF, &a_impr_mult[i],
                                           &a_impr_force[i], &a_impr_phase[i],
                                           &a_impr_term[i])) continue;

/* given the X allocations, there is hardly any point in altering on the first two
   actually, there wasn't much point in the others either */

/* if I get here, I use default values */

               a_impr_mult[i] = 101;
               a_impr_force[i] = 1.1;
               a_impr_term[i] = 2.0;
               a_impr_phase[i] = 180.0;
          }

     return;
}
