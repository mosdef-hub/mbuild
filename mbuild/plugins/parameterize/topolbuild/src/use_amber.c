/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to use the amber force field data to set parameters use_amber.c

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
#include "use_amber.h"
#include "compare_FF.h"
#include "similars.h"
#include "mol2.h"
#include "amber_FF.h"
#include "multors.h"
#include "block_memory.h"

extern short no_oppose;

/* set atom parameters for amber from force field read in.
   parameters are:
        number of atoms in molecule                  num_atoms
        number of atom types in force field          num_mol
*/
void set_amber_atomFF(int num_atoms, int num_mol)
{
     int i, j;

     for(i = 0; i < num_atoms; i++) 
          for(j = 0; j < num_mol; j++) {
               if(strcmp(a_ff_type[i],
                         ffatom_name[j])) continue;      /* this atom type? */
               a_pol[i] = ffatom_pol[j];
               break;
          }                                       /* end for(j = 0; j < num_mol; j++) */
     return;
}

/* set van der Waals atom parameters for amber from force field read in.
   parameters are:
        number of atoms in molecule:                       num_atoms
        number of van der Walls types in force field:      num_vdw
*/
void set_amber_vdwFF(int num_atoms, int num_vdw)
{
     int i, j;

     for(i = 0; i < num_atoms; i++)
          for(j = 0; j < num_vdw; j++) {
               if(strcmp(a_ff_type[i],
                         vdw_atom_name[j])) continue;
               a_vdw_atom_radius[i] = vdw_atom_radius[j];
               a_vdw_atom_pot[i]    = vdw_atom_pot[j];
               break;
          }                                       /* end for(j = 0; j < num_vdw; j++) */
     return;
}

/* set bond parameters for amber from force field read in.
   parameters are:
        number of bonds in molecule:              num_bonds
        number of bond types in force field:      bind_cntFF
*/
void set_amber_bondFF(int num_bonds, int bind_cntFF)
{
     int i, j, bnd1, bnd2;

     for(i = 0; i < num_bonds; i++) {
          bnd1 = bond_i[i];
          bnd2 = bond_j[i];

          compare_bonds(a_ff_type[bnd1], a_ff_type[bnd2],
                        bind_cntFF, i);
     }
     return;
}

/* set angle parameters for amber from force field read in.
   parameters are:
        number of angles in molecule:              num_angles
        number of angle types in force field:      angle_cntFF
*/
void set_amber_angleFF(int num_angles, int angle_cntFF)
{
     int i, j, ang1, ang2, ang3;

     for(i = 0; i < num_angles; i++) {
          if((opposite[i] == 'y') &&
             (no_oppose != 0)) {
                    a_angle_type[i] = -1;
                    continue;
          }
          ang1 = angle_tab[i][0];
          ang2 = angle_tab[i][1];
          ang3 = angle_tab[i][2];

          a_meth[i] = 0;

          angles_comp(a_ff_type[ang1], a_ff_type[ang2],
                      a_ff_type[ang3], angle_cntFF, i, 0);

          if(!a_meth[i]) {
               calc_angle(ang1, ang2, ang3, &a_angle_force[i],
                          &a_angle_angle[i], a_bond_length,
                          &a_meth[i], NULL, NULL, NULL, NULL,
                          angle_cntFF, a_ff_type);
               if(a_meth[i]) a_angle_type[i] = 1;
          }
          if(opposite[i] == 'y') {
               if(a_angle_type[i] == 1)
                    a_angle_angle[i] *= 2.0;
               else
                    a_angle_type[i] = -1;
          }
     }
     return;
}

/* set torsion parameters for amber from force field read in.
   parameters are:
        number of torsions in molecule:              num_tors
        number of torsion types in force field:      tors_cntFF
*/
void set_amber_torsFF(int num_tors, int tors_cntFF)
{
     int i, tors1, tors2, tors3, tors4, j, k;

     for(i = 0; i < num_tors; i++) {
          tors1 = torsion_tab[i][0];
          tors2 = torsion_tab[i][1];
          tors3 = torsion_tab[i][2];
          tors4 = torsion_tab[i][3];

          j = compare_tors(a_ff_type[tors1], a_ff_type[tors2],
                           a_ff_type[tors3], a_ff_type[tors4],
                           tors_cntFF, &a_tors_mult[i],
                           &a_tors_force[i], &a_tors_phase[i],
                           &a_tors_term[i]);

          if(j)
               a_tors_type[i] = 1;
          if((j == 0) || (j == (tors_cntFF - 1)))
               continue;               /* if return 0 or limit, skip tests below */

/* if this has multiple torsions, add links and set them as well
   note:  j points to the next torsion entry which is the next
   component of a multiple torsion and these facts are used below.
*/
          if(tors_more[(j - 1)]) {
               multors_cnt[i] = tors_more[(j - 1)];     /* set count */
               multiples[i] = (float **)allocator(multors_cnt[i], sizeof(float **), FARGS);
               for(k = 0; k < multors_cnt[i]; k++) {
                    multiples[i][k] = (float *)allocator(4, sizeof(float), FARGS);
                    multiples[i][k][0] = (float)tors_mult[j];
                    multiples[i][k][1] = tors_force[j];
                    multiples[i][k][2] = tors_phase[j];
                    multiples[i][k][3] = tors_term[j];
                    j++;
               }
          }
     }

    return;
}

/* set improper parameters for amber from force field read in.
   parameters are:
        number of improper in molecule:              num_impr
        number of improper types in force field:     impr_cntFF
*/
void set_amber_improFF(int num_impr, int impr_cntFF)
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
                                 &a_impr_term[i]))
                    continue;

/* if I get here, I use default values */

               a_impr_mult[i] = 101;
               a_impr_force[i] = 1.1;
               a_impr_term[i] = 2.0;
               a_impr_phase[i] = 180.0;
          }

     return;
}


