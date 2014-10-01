/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to use Tripos force field data to set parameters use_tripos.c
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "similars.h"
#include "mol2.h"
#include "amber_FF.h"
#include "Tripos_FF.h"

#define MAXCHAR 256

extern char mess[MAXCHAR];

/* set bond parameters for Tripos from force field read in.
   parameters are:
        number of bonds in molecule:              num_bonds
        number of bond types in force field:      bind_cntFF
*/
void set_tripos_bondFF(int num_bonds, int bind_cntFF)
{
     int i, j, bnd1, bnd2;

     for(i = 0; i < num_bonds; i++) {
          bnd1 = bond_i[i];
          bnd2 = bond_j[i];
          for(j = 0; j < bind_cntFF; j++) {
               if(strcmp(bond_type[i], trp_bnd_type[j]))
                    continue;                                 /* bond types must match */
               if((compare_types(bnd1, bond_name1[j])   &&
                   compare_types(bnd2, bond_name2[j]))  ||
                  (compare_types(bnd1, bond_name2[j])   &&
                   compare_types(bnd2, bond_name1[j]))) {     /* atom names must match */
                         a_bond_type[i] = 1;
                         a_bond_force[i]  = bond_force[j];
                         a_bond_length[i] = bond_length[j];
                         break;
               }
          }
     }

     return;
}

/* set angle parameters for tripos from force field read in.
   parameters are:
        number of angles in molecule:              num_angles
        number of angle types in force field:      angle_cntFF
*/
void set_tripos_angleFF(int num_angles, int angle_cntFF)
{
     int i, j, ang1, ang2, ang3;

     for(i = 0; i < num_angles; i++) {
          ang1 = angle_tab[i][0];
          ang2 = angle_tab[i][1];
          ang3 = angle_tab[i][2];

          a_meth[i] = 0;
          for(j = 0; j < angle_cntFF; j++) {
               if(strcmp(a_ff_type[ang2], angle_name2[j]))
                    continue;          /* central atom must always match */

/* try to match all atoms */
               if((compare_types(ang1, angle_name1[j])  &&
                   compare_types(ang3, angle_name3[j])) ||
                  (compare_types(ang1, angle_name3[j])  &&
                   compare_types(ang3, angle_name1[j]))) {
                         a_angle_type[i] = 1;
                         a_angle_force[i] = angle_force[j];
                         a_angle_angle[i] = angle_angle[j];
                         a_meth[i] = 1;                       /* signal found in table */
                         if(opposite[i] == 'y')
                              a_angle_angle[i] *= 2.0;
                         break;
               }                                              /* end if((compare_types    etc. */
          }
     }

     return;
}

/* set torsion parameters for tripos from force field read in.
   parameters are:
        number of torsions in molecule:              num_tors
        number of torsion types in force field:      tors_cntFF
*/
void set_tripos_torsFF(int num_tors, int tors_cntFF)
{
     int i, j, k, mybond, tors1, tors2, tors3, tors4;

     for(i = 0; i < num_tors; i++) {
          tors1 = torsion_tab[i][0];
          tors2 = torsion_tab[i][1];
          tors3 = torsion_tab[i][2];
          tors4 = torsion_tab[i][3];

          mybond = -1;

          for(j = 0; j < bond_count[tors2]; j++) {            /* find the central bond number */
               for(k = 0; k < bond_count[tors3]; k++)
                    if(bonds_list[tors2][j] == bonds_list[tors3][k]) {
                         mybond = bonds_list[tors2][j];
                         break;
                    }
               if(mybond > -1) break;
          }

          if(mybond < 0) {
               sprintf(mess,
                  "Central bond of torsion %d-%d-%d-%d not in bonds list\n",
                       tors1, tors2, tors3, tors4);
               my_fatal(FARGS, mess);
               exit(1);
          }

          for(j = 0; j < tors_cntFF; j++) {
               if(strcmp(bond_type[mybond], tors_bond[j]))
                    continue;                                 /* central bond must match */
               if((!strcmp(a_ff_type[tors2], tors_name2[j]) &&
                   !strcmp(a_ff_type[tors3], tors_name3[j]))) {      /* two center must match */
                         if(compare_types(tors1, tors_name1[j])  &&  /* 1 to 2 and 4 to 3 */
                            compare_types(tors4, tors_name4[j])) {   /* order here is important */
                                   a_tors_type[i] = 1;
                                   a_tors_force[i] = tors_force[j];
                                   a_tors_phase[i] = 0.0;
                                   if(tors_term[j] < 0.0)
                                        a_tors_phase[i] = 180.0;
                                   a_tors_term[i] = tors_term[j];
                                   a_tors_mult[i] = 1;
                                   break;
                         }
               }

               if((!strcmp(a_ff_type[tors2], tors_name3[j]) &&
                   !strcmp(a_ff_type[tors3], tors_name2[j]))) {      /* two center must match */
                         if(compare_types(tors1, tors_name4[j])  &&  /* 4 to 2 and 1 to 3 */
                            compare_types(tors4, tors_name1[j])) {   /* order here is important */
                                   a_tors_type[i] = 1;
                                   a_tors_force[i] = tors_force[j];
                                   a_tors_phase[i] = 0.0;
                                   if(tors_term[j] < 0.0)
                                        a_tors_phase[i] = 180.0;
                                   a_tors_term[i] = tors_term[j];
                                   a_tors_mult[i] = 1;
                                   break;
                         }
               }
          }
     }

     return;
}

/* set improper parameters for Tripos from force field read in.
   parameters are:
        number of improper in molecule:                       num_impr
        number of out of plane bend types in force field:     impr_cntFF
*/
void set_tripos_improFF(int num_impr, int impr_cntFF)
{
     int i, j, impr3;

     if(!num_impr)
          return;
     for(i = 0; i < num_impr; i++) {
          impr3 = improperid[i][2];

          for(j = 0; j < impr_cntFF; j++) {
               if(impro_force[j] < 0.0)
                    continue;
               if(strcmp(a_ff_type[impr3], impro_name3[j]))
                    continue;
/* correct improper force to account for change from pyramid height out
   of plane used in Sybyl to cos improper out of plane.  Correction
   factor tends to be in the neighborhood of 39.3 */
               a_impr_force[i] = impro_force[j]/39.3;
               a_impr_mult[i] = 1;
               a_impr_term[i] = 2.0;
               a_impr_phase[i] = 180.0;
          }
     }

     return;
}

/* compare Tripos types and account for the Tripos wild cards.
   return 1 for a match and 0 for failure to match.
   parameters are:
        id_atom                   number in molecule of atom to compare against
        name_of_type              the Tripos type to compare against
*/
int compare_types(int id_atom, char *name_of_type)
{
     int i;
     int halogens[4] = { 9, 17, 35, 53 };
     int hetatoms[4] = { 7,  8, 15, 16 };

     if(!strcmp(a_ff_type[id_atom], name_of_type))
          return(1);                                      /* if it matches, it matches */

     if(!strcmp(name_of_type, "X") ||
        !strcmp(name_of_type, "ANY"))                     /* wild card always matches */
               return(1);

     if(!strcmp(name_of_type, "HAL"))                     /* halogens check */
          for(i = 0; i < 4; i++) 
               if(atom_atno[id_atom] == halogens[i])
                    return(1);

     if(!strcmp(name_of_type, "HEV") &&
        (atom_atno[id_atom] > 1))                         /* atomic number > 1 check */
               return(1);

     if(!strcmp(name_of_type, "HET"))                     /* heteroatoms check */
          for(i = 0; i < 4; i++)
               if(atom_atno[id_atom] == hetatoms[i])
                    return(1);

     return(0);                                           /* did not match */
}
