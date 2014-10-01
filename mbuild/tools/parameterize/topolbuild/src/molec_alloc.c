/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Function to allocate remainder of needed memory for molecule's
        parameters            molec_alloc.c
*/
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "mol2.h"
#include "readmol2.h"
#include "multors.h"

#define NULLCHAR 0

/* allocate force field storage for the molecule
   parameters are:
        num_atoms          number of atoms in the molecule
        num_bonds          number of bonds in the molecule
        num_angles         number of angles between bonds in the molecule
        num_tors           number of torsions around bonds in the molecule
        num_impr           number of improper torsions in the molecule
*/
void molec_alloc(int num_atoms, int num_bonds, int num_angles, int num_tors,
                 int num_impr)
{
     int i;

     a_pol = (float *)allocator(num_atoms, sizeof(float), FARGS);
     a_bond_force = (float *)allocator(num_bonds, sizeof(float), FARGS);
     a_bond_length = (float *)allocator(num_bonds, sizeof(float), FARGS);
     a_bond_beta = (float *)allocator(num_bonds, sizeof(float), FARGS);
     a_bond_type = (short *)allocator(num_bonds, sizeof(short), FARGS);
     a_angle_quartic = (float **)allocator(5, sizeof(float*), FARGS);
     for(i = 0; i < 5; i++)
          a_angle_quartic[i] = (float *)allocator(num_angles,
                                                  sizeof(float), FARGS);
     a_angle_force = a_angle_quartic[0];
     a_angle_angle = (float *)allocator(num_angles, sizeof(float), FARGS);
     a_meth = (int *)allocator(num_angles, sizeof(int), FARGS);
     a_angle_type = (short *)allocator(num_angles,  sizeof(short), FARGS);
     a_tors_set = (float **)allocator(6, sizeof(float*), FARGS);
     for(i = 0; i < 6; i++)
          a_tors_set[i] = (float *)allocator(num_tors, sizeof(float),
                                             FARGS);
     a_tors_mult = (int *)a_tors_set[0];
     a_tors_force = a_tors_set[1];
     a_tors_phase = a_tors_set[2];
     a_tors_term = a_tors_set[3];
     multors_cnt = (int *)allocator(num_tors, sizeof(int), FARGS);
     multiples = (float ***)allocator(num_tors, sizeof(float **), FARGS);
     a_tors_type = (short *)allocator(num_tors, sizeof(short), FARGS);
     if(num_impr) {
          a_impr_mult = (int *)allocator(num_impr, sizeof(int), FARGS);
          a_impr_force = (float *)allocator(num_impr, sizeof(float), FARGS);
          a_impr_phase = (float *)allocator(num_impr, sizeof(float), FARGS);
          a_impr_term = (float *)allocator(num_impr, sizeof(float), FARGS);
     }
     else {
          a_impr_mult = NULL;
          a_impr_force = NULL;
          a_impr_phase = NULL;
          a_impr_term = NULL;
     }
     a_vdw_atom_radius = (float *)allocator(num_atoms, sizeof(float),
                                            FARGS);
     a_vdw_atom_pot = (float *)allocator(num_atoms, sizeof(float), FARGS);

/* set everything to mark unused (don't use 0 here except on method flags) */

     for(i = 0; i < num_atoms; i++){
          a_pol[i] = 99999.0;
          a_vdw_atom_radius[i] = 99999.0;
          a_vdw_atom_pot[i] =  99999.0;
     }

     for(i = 0; i < num_bonds; i++){
          a_bond_force[i] = -1.0;
          a_bond_type[i] = 0;
     }

     for(i = 0; i < num_angles; i++){
          a_angle_force[i] = -1.0;
          a_angle_type[i] = 0;
          a_meth[i] = 0;
     }

     for(i = 0; i < num_tors; i++){
          a_tors_mult[i] = 9999;
          a_tors_force[i] = -1.0;
          a_tors_phase[i] = 99999.0;
          a_tors_term[i] = 99999.0;
          multors_cnt[i] = 0;
          multiples[i] = NULL;
          a_tors_type[i] = 0;
     }

     if(num_impr)
          for(i = 0; i < num_impr; i++){
               a_impr_mult[i] = 9999;
               a_impr_force[i] = -1.0;
               a_impr_phase[i] = 99999.0;
               a_impr_term[i] = 99999.0;
          }
     return;
}

/* Clean up 0 residue numbers and set order of atoms according to  
   residue. to residue.  Solves problem of all residue atoms required
   to be together in gromacs topologies and of zero mass things.
   Much easier to set up a forward and reverse table look-up
   system than to try sort atoms, bonds, and connection tables.
   parameters are:
        num_atoms                number of atoms in molecule
        resid                    residue name if set by -r option.
                                 NULLCHAR otherwise.
*/
void atom_residue_order(int num_atoms, char *resid)
{
     int i, j, k, nex, count, maxsegno;

     maxsegno = 1;
     for(i = 0; i < num_atoms; i++)
          if(atom_segno[i] > maxsegno) maxsegno = atom_segno[i];

     maxsegno++;

     for(i = 0; i < num_atoms; i++)
          if(!atom_segno[i]) atom_segno[i] = maxsegno;

     for(i = 0; i < num_atoms; i++) {
          atom_order[i] = -1;
          back_order[i] = -1;
     }

     count = 0;
     j = 0;
     if(resid[0] != NULLCHAR) {
          while(j < num_atoms) {
               for(i = 0; i < num_atoms; i++) {
                    if(back_order[i] > -1)
                         continue;                 /* already set? */
                    if(atom_mass[i] == 0.0) {
                         atom_order[j] = -1;
                         back_order[i] = -1;
                         j++;
                         if(bond_count[i])
                              atom_charge[(connect[i][0])] += atom_charge[i];
                         continue;
                    }
                    atom_order[j] = i;
                    back_order[i] = count;
                    j++;
                    count++;
/* set attached H's next to parent */
                    if((bond_count[i] > 0) && (atom_atno[i] > 1)) {
                         for(nex = 0; nex < bond_count[i]; nex++) {
                              if(atom_mass[(connect[i][nex])] == 0.0)
                                   continue;
                              if(atom_atno[(connect[i][nex])] == 1) {
                                   atom_order[j] = connect[i][nex];
                                   back_order[(connect[i][nex])] = count;
                                   j++;
                                   count++;
                              }       /* end if(atom_atno[(connect[i][nex])] == 1) */
                         }            /* end for(nex = 0; nex < bond_count[i]; nex++) */
                    }                 /* end if((bond_count[i] > 0) && (atom_atno[i] > 1)) */
               }                      /* end for(i = 0; i < num_atoms; i++) */
          }                           /* end while(j < num_atoms) */

          true_atoms = count;
          return;
     }                               /* end if(resid[0] != NULLCHAR) */

     j = 0;
     k = 0;

     while(j < num_atoms) {
          k++;
          for(i = 0; i < num_atoms; i++) {
               if(back_order[i] > -1)
                    continue;                      /* already set? */
               if(k == atom_segno[i]) {
                    if(atom_mass[i] == 0.0) {
                         atom_order[j] = -1;
                         back_order[i] = -1;
                         j++;
                         if(bond_count[i])
                              atom_charge[(connect[i][0])] += atom_charge[i];
                         continue;
                    }
                    atom_order[j] = i;
                    back_order[i] = count;
                    j++;
                    count++;
/* set attached H's next to parent */
                    if((bond_count[i] > 0) && (atom_atno[i] > 1)) {
                         for(nex = 0; nex < bond_count[i]; nex++) {
                              if(atom_mass[(connect[i][nex])] == 0.0)
                                   continue;
                              if(atom_atno[(connect[i][nex])] == 1) {
                                   atom_order[j] = connect[i][nex];
                                   back_order[(connect[i][nex])] = count;
                                   j++;
                                   count++;
                              }       /* end if(atom_atno[(connect[i][nex])] == 1) */
                         }            /* end for(nex = 0; nex < bond_count[i]; nex++) */
                    }                 /* end if((bond_count[i] > 0) && (atom_atno[i] > 1)) */
               }                      /* end if(k == atom_segno[i]) */
          }                           /* end for(i = 0; i < num_atoms; i++) */
     }                                /* while(j < num_atoms) */

     true_atoms = count;
     return;
}
