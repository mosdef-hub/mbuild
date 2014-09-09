/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Function for angle and torsion tables set up param_tables.c
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "param_tables.h"
#include "mol2.h"
#include "rings.h"
#include "multors.h"
#include "GMX_FF.h"
#include "mainchain.h"

#define MAXCHAR 256
#define MAXNAME 128
#define NULLCHAR 0

/* determine the atoms that participate in angles  
   parameters are:
        just_count                    flag to just count the number of table entries
        max_ang                       size of table to allocate for angles
        num_atoms                     number of atoms in molecule
        num_bonds                     number of bonds in molecule
        num_angs                      return of number of angles found.
*/
void find_angles(int just_count, int max_ang, int num_atoms, int num_bonds,
                 int *num_angs)
{
     int i, j, k, bond1, bond2, ang_count, l, m;
     int *used_bonds;
     float measured[8], mytors;

     if(just_count) {
          angle_tab = (int **)allocator(max_ang, sizeof(int*), FARGS);
          opposite = (char *)allocator(max_ang, sizeof(char), FARGS);
          for(i = 0; i < max_ang; i++) {
               angle_tab[i] = (int *)allocator(3, sizeof(int), FARGS);
               opposite[i] = NULLCHAR;
          }
     }

     used_bonds = (int *)allocator(num_bonds, sizeof(int), FARGS);
     for(i = 0; i < num_bonds; i++) used_bonds[i] = 0;

     ang_count = 0;

     for(i = 0; i < num_bonds; i++) {
          bond1 = bond_i[i];
          bond2 = bond_j[i];

          if(bond_count[bond1] > 1) {
               for(j = 0; j < bond_count[bond1]; j++) {
                    if(bonds_list[bond1][j] == i) continue;
                    if(used_bonds[(bonds_list[bond1][j])]) continue;
                    if(just_count) {
                         if(ang_count >= max_ang) {
                              my_fatal(FARGS,
                                       "Found more angles than maximum allowed\n");
                              exit(1);
                         }
                         k = bond_i[(bonds_list[bond1][j])];
                         if((k == bond1) || (k == bond2))
                              k = bond_j[(bonds_list[bond1][j])];
                         angle_tab[ang_count][1] = bond1;
                         angle_tab[ang_count][2] = bond2;
                         angle_tab[ang_count][0] = k;
                    }                   /* end if(just_count) */
                    ang_count++;
               }                        /* end for(j = 0; j < bond_count[bond1]; j++) */
          }                             /* end if(bond_count[bond1] > 1) */

          if(bond_count[bond2] > 1) {
               for(j = 0; j < bond_count[bond2]; j++) {
                    if(bonds_list[bond2][j] == i) continue;
                    if(used_bonds[(bonds_list[bond2][j])]) continue;
                    if(just_count) {
                         if(ang_count >= max_ang) {
                              my_fatal(FARGS,
                                       "Found more angles than maximum allowed\n");
                              exit(1);
                         }
                         k = bond_i[(bonds_list[bond2][j])];
                         if((k == bond1) || (k == bond2))
                              k = bond_j[(bonds_list[bond2][j])];
                         angle_tab[ang_count][0] = bond1;
                         angle_tab[ang_count][1] = bond2;
                         angle_tab[ang_count][2] = k;
                    }                   /* end if(just_count) */
                    ang_count++;
               }                        /* end for(j = 0; j < bond_count[bond2]; j++) */
          }                             /* end if(bond_count[bond2] > 1) */
          used_bonds[i] = 1;
     }                                  /* end for(i = 0; i < num_bonds; i++) */

     if(just_count) {                   /* set planar opposition flags */
          for(i = 0; i < ang_count; i++) {
               bond1 = angle_tab[i][1];
               if((bond_count[bond1] < 4) ||
                  (bond_count[bond1] > 8) || 
                  (atom_atno[bond1] == 6) ||
                  (atom_atno[bond1] == 7))
                         continue;      /* need 4 bonds to have opposites */
               bond2 = angle_tab[i][0];
               m = angle_tab[i][2];
               if((atom_atno[bond2] != 7) ||
                  (atom_atno[m] != 7))
                         continue;      /* and need bonded to N's */
               if((bond_count[bond2] != 3) ||
                  (bond_count[m] != 3))
                         continue;      /* N's with 3 bonds to make rings */
               l = 0;
               for(j = 0; j < 3; j++)
                    if((atom_atno[(connect[bond2][j])] < 6) ||
                       (atom_atno[(connect[m][j])] < 6))
                              l = 1;
               if(l) continue;          /* N but never in a ring */
               measured[0] = angle_comp(bond2, bond1, m);
               if(measured[0] < 0.0)
                    measured[0] += 360.0;
               l = 1;
               for(j = 0; j < bond_count[bond1]; j++) {
                    k = connect[bond1][j];
                    if((k == bond2) ||
                       (k == m))
                              continue;
                    if(atom_atno[k] != 7)
                         continue;      /* bond not to N */
/* check for all approximately in same plane so that the 
   dihedral angle is either around 0 or around 180 deg. */
                    mytors = dihed(bond2,bond1,m,k);
                    if(mytors > 180.0)
                         mytors -= 360.0;
                    if((fabs(mytors) > 30.0) &&
                       (fabs(mytors) < 150.0))
                              continue;     /* Out of range */

                    measured[l] = angle_comp(bond2, bond1, k);
                    if(measured[l] < 0.0)
                         measured[l] += 360.0;
                    l++;
               }                        /* end for(j = 0; j < 4; j++) */

               if(l < 2) continue;         /* broke on not all bonded to N's */
               k = 0;
               for(j = 1; j < l; j++)
                    if(measured[0] > (1.2 * measured[j]))
                         k++;
               if(k < 2) continue;      /* Atom not intervening in space */
               opposite[i] = 'y';
          }                             /* end for(i = 0; i < ang_count; i++) */
     }

     free_me(used_bonds, FARGS);

     *num_angs = ang_count;
     return;
}

/* use the list of bonds as the bonds for which torsions are desired 
   parameters are:
        just_count                    flag to just count the number of table entries
        max_tors                      size of table to allocate for torsions
        num_atoms                     number of atoms in molecule
        num_bonds                     number of bonds in molecule
        num_tors                      return of number of torsions found.
*/
void find_torsions(int just_count, int max_tors, int num_atoms, int num_bonds,
                   int *num_tors)
{
     int i, j, k, bond1, bond2, tors_count, start, end;

     if(just_count) {
          torsion_tab = (int **)allocator(max_tors, sizeof(int*), FARGS);
          for(i = 0; i < max_tors; i++)
               torsion_tab[i] = (int *)allocator(4, sizeof(int), FARGS);
     }

     tors_count = 0;

     for(i = 0; i < num_bonds; i++) {
          bond1 = bond_i[i];
          bond2 = bond_j[i];

          if((bond_count[bond1] == 1) ||
             (bond_count[bond2] == 1))
                    continue;    /* need four atoms for torsion */

          for(j = 0; j < bond_count[bond1]; j++) {
               if(bonds_list[bond1][j] == i) continue;
               start = bond_i[(bonds_list[bond1][j])];
                    if(start == bond1)
                         start = bond_j[(bonds_list[bond1][j])];
               for(k = 0; k < bond_count[bond2]; k++) {
                    if(bonds_list[bond2][k] == i) continue;
                    end = bond_i[(bonds_list[bond2][k])];
                    if(end == bond2)
                         end = bond_j[(bonds_list[bond2][k])];
                    if(just_count) {
                         if(tors_count >= max_tors) {
                              my_fatal(FARGS,
                                 "Found more torsions than maximum allowed\n");
                         }
                         torsion_tab[tors_count][0] = start;
                         torsion_tab[tors_count][1] = bond1;
                         torsion_tab[tors_count][2] = bond2;
                         torsion_tab[tors_count][3] = end;
                    }                   /* end if(just_count) */
                    tors_count++;
               }                        /* end for(k = 0; k < bond_count[bond2]; k++) */
          }                             /* end for(j = 0; j < bond_count[bond1]; j++) */
     }                                  /* end for(i = 0; i < num_bonds; i++) */

     *num_tors = tors_count;

     return;
}

/* compute angles from coordinate data.
   parameter is number of angles found
*/
void measure_angles(int num_angs)
{
     int i;

     meas_ang = (float *)allocator(num_angs, sizeof(float), FARGS);

     for(i = 0; i < num_angs; i++) {
          meas_ang[i] = angle_comp(angle_tab[i][0], angle_tab[i][1], angle_tab[i][2]);
     }

     return;
}

/* compute an angle
   parameters are
        i, j, k          the indices of the atoms of the angle
*/
double angle_comp(int i, int j, int k)
{
     double delx1, delx2, dely1, dely2, delz1, delz2, denom;
     double my_cos, ang_meas;

     delx1 = atom_x[i] - atom_x[j];
     delx2 = atom_x[k] - atom_x[j];
     dely1 = atom_y[i] - atom_y[j];
     dely2 = atom_y[k] - atom_y[j];
     delz1 = atom_z[i] - atom_z[j];
     delz2 = atom_z[k] - atom_z[j];
     denom = ((delx1 * delx1) + (dely1 * dely1) + (delz1 * delz1));
     denom = denom * ((delx2 * delx2) + (dely2 * dely2) + (delz2 * delz2));
     denom = sqrt(denom);
     my_cos = ((delx1 * delx2) + (dely1 * dely2) + (delz1 * delz2)) / denom;
     ang_meas = acos(my_cos) * (180.0/M_PI);

     return(ang_meas);
}

/* compute torsions from coordinate data.
   parameter is number of torsions found
*/
void measure_torsions(int num_tors)
{
     int i;
     double my_tor;

     meas_tors = (float *)allocator(num_tors, sizeof(float), FARGS);

     for(i = 0; i < num_tors; i++) {
          my_tor = dihed(torsion_tab[i][0], torsion_tab[i][1],
                        torsion_tab[i][2], torsion_tab[i][3]);
          if(my_tor < 0.0) my_tor += 360.0;
          meas_tors[i] = my_tor;
     }

     return;
}

/* compute bond lengths from coordinates
   parameter is number of bonds in molecule
*/
void measure_bondlen(int num_bonds)
{
     int i;
     double delx, dely, delz, my_dist;

     meas_len = (float *)allocator(num_bonds, sizeof(float), FARGS);

     for(i = 0; i < num_bonds; i++) {
          delx = atom_x[(bond_i[i])] - atom_x[(bond_j[i])];
          dely = atom_y[(bond_i[i])] - atom_y[(bond_j[i])];
          delz = atom_z[(bond_i[i])] - atom_z[(bond_j[i])];
          my_dist = sqrt((delx * delx) + (dely * dely) + (delz * delz));
          meas_len[i] = my_dist;
     }

     return;
}

/* compute improper torsions from coordinate data.
   parameter is number of impropers found
*/
void measure_improp(int num_improp)
{
     int i;
     double my_tor;

     if(!num_improp) {
          meas_impro = NULL;
          return;
     }

     meas_impro = (float *)allocator(num_improp, sizeof(float), FARGS);

     for(i = 0; i <  num_improp; i++) {
          my_tor = dihed(improperid[i][2], improperid[i][0],
                         improperid[i][1], improperid[i][3]);
          if(my_tor < 0.0) my_tor += 360.0;
          my_tor = ((double)((int)((my_tor / 5.0) + 0.5)) * 5.0);
          if(my_tor == 360.0) my_tor = 0.0;
          meas_impro[i] = my_tor;
     }

     return;
}

/* compute a dihedral
   parameters are:
        i, j, k, l         the indices of the atoms of the dihedral
*/

double dihed(int i, int j, int k, int l)
{
     double delx1, delx2, delx3, dely1, dely2, dely3, delz1, delz2, delz3;
     double prodax, proday, prodaz, prodbx, prodby, prodbz, prodcx, prodcy, prodcz;
     double norma, normb, normc, mysin, mycos, denom1, denom2, asign, mytor;

     delx1 = atom_x[i] - atom_x[j];
     dely1 = atom_y[i] - atom_y[j];
     delz1 = atom_z[i] - atom_z[j];
     delx2 = atom_x[j] - atom_x[k];
     dely2 = atom_y[j] - atom_y[k];
     delz2 = atom_z[j] - atom_z[k];
     delx3 = atom_x[k] - atom_x[l];
     dely3 = atom_y[k] - atom_y[l];
     delz3 = atom_z[k] - atom_z[l];

     prodax = (dely1 * delz2) - (delz1 * dely2);
     proday = (delz1 * delx2) - (delx1 * delz2);
     prodaz = (delx1 * dely2) - (dely1 * delx2);
     norma = (prodax * prodax) + (proday * proday) + (prodaz * prodaz);

     prodbx = (dely2 * delz3) - (delz2 * dely3);
     prodby = (delz2 * delx3) - (delx2 * delz3);
     prodbz = (delx2 * dely3) - (dely2 * delx3);
     normb = (prodbx * prodbx) + (prodby * prodby) + (prodbz * prodbz);

     prodcx = (dely2 * prodaz) - (delz2 * proday);
     prodcy = (delz2 * prodax) - (delx2 * prodaz);
     prodcz = (delx2 * proday) - (dely2 * prodax);
     normc = (prodcx * prodcx) + (prodcy * prodcy) + (prodcz * prodcz);

     denom1 = sqrt((norma * normb));
     denom2 = sqrt((normb * normc));
     if(denom1 < 1.0e-10) denom1 = 1.0e-10;
     if(denom2 < 1.0e-10) denom2 = 1.0e-10;

     mycos = ((prodax * prodbx) + (proday * prodby) + (prodaz * prodbz))/denom1;
     mysin = ((prodcx * prodbx) + (prodcy * prodby) + (prodcz * prodbz))/denom2;

     asign = -1.0;
     if(mysin < 0.0) asign = 1.0;
     mytor = asign * acos(mycos) * (180.0/M_PI);

     return(mytor);
}

/* fix measured torsions to fit the multiplicities from the force field better.
   Also fix measured improper values.
   parameters are:
        count                 number of torsions
        num_impr              total number of improper dihedrals to check
*/
void fix_tor(int count, int num_impr)
{
     int i;

     for(i = 0; i < count; i++) {
          if((a_tors_type[i] != 1) &&
             (a_tors_type[i] != 9))
                    continue;
          if(abs(((int)a_tors_term[i])) < 3) {    /* multiplicities of 1 or 2 are */
                                                  /* set to 180.0 or 0.0 */
               if((meas_tors[i] <= 185.0) &&
                  (meas_tors[i] >= 175.0))        /* but only if in range of 180.0, */
                                                  /* 360.0 or 0.0 */
                         meas_tors[i] = 180.0;    /* don't alter outliers.  They're real */
               if(meas_tors[i] >= 355.0)
                    meas_tors[i] = 0.0;
               if((meas_tors[i] >= -5.0) &&
                  (meas_tors[i] <= 5.0))
                         meas_tors[i] = 0.0;
               continue;                                 /* done here */
          }
          if(abs(((int)a_tors_term[i])) == 3.0) {   /* multiplicities of 3 are set to */
                                                    /* 0.0, 120.0 or 240.0 */
               if(meas_tors[i] >= 355.0)
                    meas_tors[i] = 0.0;
               if((meas_tors[i] >= -5.0) &&
                  (meas_tors[i] <= 5.0))     /* but only if in range of 0.0, 120.0, or 240.0 */
                         meas_tors[i] = 0.0;          /* don't alter outliers.  They're real */
               if((meas_tors[i] >= 115.0) &&
                  (meas_tors[i] <= 125.0))
                         meas_tors[i] = 120.0;
               if((meas_tors[i] >= 235.0) &&
                  (meas_tors[i] <= 245.0))
                         meas_tors[i] = 240.0;
          }
     }

     if(num_impr)
          for(i = 0; i < num_impr; i++){
               if((meas_impro[i] <= 185.0)  &&
                  (meas_impro[i] >= 175.0))
                         meas_impro[i] = 180.0;
/* scale to between -180 and +180 */
               if(meas_impro[i] > 180.0)
                    meas_impro[i] -= 360.0;
               if(meas_impro[i] < -180.0)
                    meas_impro[i] += 360.0;
               if((meas_impro[i] >= -5.0) &&
                  (meas_impro[i] <= 5.0))
                          meas_impro[i] = 0.0;
          }                             /* end for(i = 0; i < num_impr; i++) */

     return;
}

/* test for flat ring torsion within a single ring
   and set appropriate angular values
   parameter is
        i          torsion to test
*/
int single_flat_ring(int i)
{
     int j, k, l;
     j = 0;
     if(a_tors_type[i] != 1)
          return(0);               /* only for type 1 torsion */
     if(multors_cnt[i])
          return(0);               /* only if single torsion present */

/* All 4 must be in a single planar ring */
     if((ar_set[(torsion_tab[i][0])][1] >= 1) &&
        (ar_set[(torsion_tab[i][1])][1] >= 1) &&
        (ar_set[(torsion_tab[i][2])][1] >= 1) &&
        (ar_set[(torsion_tab[i][3])][1] >= 1)) j++;
     if((ar_set[(torsion_tab[i][0])][2] >= 1) &&
        (ar_set[(torsion_tab[i][1])][2] >= 1) &&
        (ar_set[(torsion_tab[i][2])][2] >= 1) &&
        (ar_set[(torsion_tab[i][3])][2] >= 1)) j++;
     if((ar_set[(torsion_tab[i][0])][3] >= 1) &&
        (ar_set[(torsion_tab[i][1])][3] >= 1) &&
        (ar_set[(torsion_tab[i][2])][3] >= 1) &&
        (ar_set[(torsion_tab[i][3])][3] >= 1)) j++;

     if(!j)
          return(0);               /* all 4 atoms in planar ring? */

     for(j = 0; j < num_rngs; j++) {
          l = 0;
          if(ring_num[j] < 4)
               continue;              /* must be cyclobutane or larger */
          if((ring_type[j] > 3) || (ring_type[j] < 1))
               continue;              /* only test if planar ring type */
          for(k = 0; k < ring_num[j]; k++) {
               if((ring_atomno[j][k] == torsion_tab[i][0]) ||
                  (ring_atomno[j][k] == torsion_tab[i][1]) ||
                  (ring_atomno[j][k] == torsion_tab[i][2]) ||
                  (ring_atomno[j][k] == torsion_tab[i][3]))
                         l++;
          }                        /* end for(k = 0; k < ring_num[j]; k++) */
          if(l == 4) {             /* all 4 in the same ring? */
               if((meas_tors[i] <= 185.0)  &&
                  (meas_tors[i] >= 175.0))
                         meas_tors[i] = 180.0;
/* scale to between -180 and +180 */
               if(meas_tors[i] >  180.0)
                    meas_tors[i] -= 360.0;
               if(meas_tors[i] < -180.0)
                    meas_tors[i] += 360.0;
               if((meas_tors[i] >= -5.0) &&
                  (meas_tors[i] <= 5.0))
                         meas_tors[i] = 0.0;
               return(1);
          }                        /* end if(l == 4) */
     }                             /* end for(j = 0; j < num_rngs; j++) */

     return(0);
}

/* test for flatness torsion across two rings and set
   angular values
   parameter is
        i          torsion to test
*/
int cross_flat_ring(int i)
{
     int j, k, l;
     j = 0;
     if(a_tors_type[i] != 1)
          return(0);               /* only for type 1 torsion */
     if(multors_cnt[i])
          return(0);               /* only if single torsion present */

/* all 4 must be in planar rings */
     if((ar_set[(torsion_tab[i][0])][1] >= 1) ||
        (ar_set[(torsion_tab[i][0])][2] >= 1) ||
        (ar_set[(torsion_tab[i][0])][3] >= 1)) j++;
     if((ar_set[(torsion_tab[i][1])][1] >= 1) ||
        (ar_set[(torsion_tab[i][1])][2] >= 1) ||
        (ar_set[(torsion_tab[i][1])][3] >= 1)) j++;
     if((ar_set[(torsion_tab[i][2])][1] >= 1) ||
        (ar_set[(torsion_tab[i][2])][2] >= 1) ||
        (ar_set[(torsion_tab[i][2])][3] >= 1)) j++;
     if((ar_set[(torsion_tab[i][3])][1] >= 1) ||
        (ar_set[(torsion_tab[i][3])][2] >= 1) ||
        (ar_set[(torsion_tab[i][3])][3] >= 1)) j++;

     if(j != 4)
          return(0);
     l = ar_set[(torsion_tab[i][1])][1] +
         ar_set[(torsion_tab[i][1])][2] +
         ar_set[(torsion_tab[i][1])][3];
     k = ar_set[(torsion_tab[i][2])][1] +
         ar_set[(torsion_tab[i][2])][2] +
         ar_set[(torsion_tab[i][2])][3];

     if((k != 2) && (l != 2))
          return(0);              /* must have at least one in two planar rings */
     for(j = 0; j < num_rngs; j++) {   /* must not have both ends in same ring */
          if(ring_num[j] < 4)
               continue;              /* must be cyclobutane or larger */
          if((ring_type[j] > 3) || (ring_type[j] < 1))
               continue;              /* only test if planar ring type */
          l = 0;
          for(k = 0; k < ring_num[j]; k++) {
               if((ring_atomno[j][k] == torsion_tab[i][0]) ||
                  (ring_atomno[j][k] == torsion_tab[i][3]))
                         l++;
          }                        /* end for(k = 0; k < ring_num[j]; k++) */
          if(l == 2)
               return(0);
     }

     a_tors_phase[i] = 180.0;       /* correct the phase setting for two planar rings */
     if((meas_tors[i] <= 185.0)  &&
        (meas_tors[i] >= 175.0))
               meas_tors[i] = 180.0;
/* scale to between -180 and +180 */
     if(meas_tors[i] >  180.0)
          meas_tors[i] -= 360.0;
     if(meas_tors[i] < -180.0)
          meas_tors[i] += 360.0;
     if((meas_tors[i] >= -5.0) &&
        (meas_tors[i] <=  5.0))
               meas_tors[i] = 0.0;
     return(1);
}

void amber_flat_rings(int num_tors)
{
     int i, j;
     if(!num_rngs) return;
     tors_imp = 0;
     for(i = 0; i < num_tors; i++)
          j = single_flat_ring(i);

/* Case of a torsion that bridges two planar rings.
   Note:  In this case, the correct dihedral angle is 180.
*/
     for(i = 0; i < num_tors; i++)
          j = cross_flat_ring(i);

     return;
}

/* Set dihedrals that are for flat rings to type 2 to be improper
   parameters are:
        num_tors              total number of dihedrals to check
*/
void correct_dihed_types(int num_tors)
{
     int i;

     if(!num_rngs) return;
     tors_imp = 0;
     for(i = 0; i < num_tors; i++) {
          if(single_flat_ring(i)) {
               a_tors_type[i] = 2;
               tors_imp++;
          }
     }                                  /* end for(i = 0; i < num_tors; i++) */

/* Case of a torsion that bridges two planar rings.
   Note:  In this case, the correct dihedral angle is 180.
*/
     for(i = 0; i < num_tors; i++) {
          if(cross_flat_ring(i)) {
               a_tors_type[i] = 2;
               tors_imp++;
          }
     }

     return;
}

/* code to set wanted dihedrals flags for gromacs force fields
   parameter is:
        tors_max                   number of torsions in molecule
        purge_level                level at which to remove
*/
void set_wanted_tors(int tors_max, int purge_level)
{
     int i, j, k, l;

/* first assign storage for flags for needed torsions */

     wanted_tors = (short *)allocator(tors_max, sizeof(short), FARGS);

     for(i = 0; i < tors_max; i++)
          wanted_tors[i] = 0;             /* initialize as not wanted */

/* flag torsions that are absolutely needed */
/* First flag those for which the central bond is main chain atoms in rings. */
     for(i = 0; i < tors_max; i++) {
          wanted_tors[i] = 0;             /* initialize as not wanted */
          if((atom_atno[(torsion_tab[i][0])] < 2 )  ||
             (atom_atno[(torsion_tab[i][3])] < 2 ))
                    continue;             /* No H, lp, or du entries needed */
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          if((arom_num[(torsion_tab[i][0])] == 0) &&
             (arom_num[(torsion_tab[i][1])] == 0) &&
             (arom_num[(torsion_tab[i][2])] == 0) &&
             (arom_num[(torsion_tab[i][3])] == 0)) {  /* only those in rings */
                    k = 0;
                    for(j = 0; j < main_num; j++)
                         if((main_chain[j] == torsion_tab[i][1])  ||
                            (main_chain[j] == torsion_tab[i][2])) k++;
                    if(k == 2)
                         wanted_tors[i] = 1;    /* only if both are in the main chain */
          }
     }                                  /* end for(i = 0; i < tors_max; i++) */

     if(purge_level > 4) return;

/* Next flag others for which all atoms are in the same ring */
     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i])
               continue;           /* already wanted.  More checks not needed */
          if((atom_atno[(torsion_tab[i][0])] < 2 )  ||
             (atom_atno[(torsion_tab[i][3])] < 2 ))
                    continue;               /* No H, lp, or du entries needed */
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          if(!num_rngs)
               continue;      /* ring checks not needed if don't have any rings */
          if((arom_num[(torsion_tab[i][0])] > 0) ||
             (arom_num[(torsion_tab[i][1])] > 0) ||
             (arom_num[(torsion_tab[i][2])] > 0) ||
             (arom_num[(torsion_tab[i][3])] > 0))
                     continue;              /* skip if dihedral not in any ring */
          for(j = 0; j < num_rngs; j++) {
               if(!ring_num[j]) continue;
               l = 0;
               for(k = 0; k < ring_num[j]; k++)
                    if((ring_atomno[j][k] == torsion_tab[i][0]) ||
                       (ring_atomno[j][k] == torsion_tab[i][1]) ||
                       (ring_atomno[j][k] == torsion_tab[i][2]) ||
                       (ring_atomno[j][k] == torsion_tab[i][3]))
                              l++;
               if(l == 4) {
                    wanted_tors[i] = 2;
                    break;
               }
          }                             /* end for(j = 0; j < num_rngs; j++) */
     }                                  /* end for(i = 0; i < tors_max; i++) */

     if(purge_level > 3) return;

/* next flag all for which all are main chain atoms */
     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i])
               continue;           /* already wanted.  More checks not needed */
          if((atom_atno[(torsion_tab[i][0])] < 2 )  ||
             (atom_atno[(torsion_tab[i][3])] < 2 ))
                    continue;               /* No H, lp, or du entries needed */
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          l = 0;
          for(j = 0; j < main_num; j++)
               if((main_chain[j] == torsion_tab[i][0])  ||
                  (main_chain[j] == torsion_tab[i][1])  ||
                  (main_chain[j] == torsion_tab[i][2])  ||
                  (main_chain[j] == torsion_tab[i][3])) l++;

          if(l != 4) continue;
          k = 0;
          for(j = 0; j < tors_max; j++) {
               if(j == i) continue;
               if(((torsion_tab[i][1] == torsion_tab[j][1])   &&
                   (torsion_tab[i][2] == torsion_tab[j][2]))  ||
                  ((torsion_tab[i][2] == torsion_tab[j][1])   &&
                   (torsion_tab[i][1] == torsion_tab[j][2])))
                         if(wanted_tors[j]) {
                              k = 1;
                              break;
                         }              /* end if(wanted_tors[j]) */
          }                             /* end for(j = 0; j < tors_max; j++) */
          if(!k)
               wanted_tors[i] = 3;      /* only if all in the main chain */
     }                                  /* end for(i = 0; i < tors_max; i++) */

     if(purge_level > 2) return;

/* Next be sure every central bond of every acceptable torsion has at least one
   wanted torsion set */
/* Do this first without anything that ends in an H to ensure that most times the
   wanted torsion uses heavy atoms at both ends rather than the H's. */
     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i]) continue;
          if((atom_atno[(torsion_tab[i][0])] < 2 )  ||
             (atom_atno[(torsion_tab[i][3])] < 2 ))
                    continue;               /* No H, lp, or du entries needed */
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          k = 0;
          for(j = 0; j < tors_max; j++) {
               if(j == i) continue;
               if(((torsion_tab[i][1] == torsion_tab[j][1])   &&
                   (torsion_tab[i][2] == torsion_tab[j][2]))  ||
                  ((torsion_tab[i][2] == torsion_tab[j][1])   &&
                   (torsion_tab[i][1] == torsion_tab[j][2])))
                         if(wanted_tors[j]) {
                              k = 1;
                              break;
                         }              /* end if(wanted_tors[j]) */
          }                             /* end for(j = 0; j < tors_max; j++) */
          
          if(!k)
               wanted_tors[i] = 4;      /* At least one lacks an improper */
     }                                  /* end for(i = 0; i < tors_max; i++) */

/* Repeat allowing H on one end and not the other to ensure preference for
   wanted torsions with maximum number of heavy atoms. */ 
     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i]) continue;
          if((atom_atno[(torsion_tab[i][0])] < 2 )  &&
             (atom_atno[(torsion_tab[i][3])] < 2 ))
                    continue;               /* No H, lp, or du entries needed */
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          k = 0;
          for(j = 0; j < tors_max; j++) {
               if(j == i) continue;
               if(((torsion_tab[i][1] == torsion_tab[j][1])   &&
                   (torsion_tab[i][2] == torsion_tab[j][2]))  ||
                  ((torsion_tab[i][2] == torsion_tab[j][1])   &&
                   (torsion_tab[i][1] == torsion_tab[j][2])))
                         if(wanted_tors[j]) {
                              k = 1;
                              break;
                         }              /* end if(wanted_tors[j]) */
          }                             /* end for(j = 0; j < tors_max; j++) */
          
          if(!k)
               wanted_tors[i] = 4;      /* At least one lacks an improper */
     }                                  /* end for(i = 0; i < tors_max; i++) */

     if(purge_level > 1) return;

/* Last include those that have an H on an end. */
     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i]) continue;
          if((atom_atno[(torsion_tab[i][0])] != 1) &&
             (atom_atno[(torsion_tab[i][3])] != 1))
                    continue;
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          k = 0;
          for(j = 0; j < tors_max; j++) {
               if(j == i) continue;
               if(((torsion_tab[i][1] == torsion_tab[j][1])   &&
                   (torsion_tab[i][2] == torsion_tab[j][2]))  ||
                  ((torsion_tab[i][2] == torsion_tab[j][1])   &&
                   (torsion_tab[i][1] == torsion_tab[j][2])))
                         if(wanted_tors[j]) {
                              k = 1;
                              break;
                         }              /* end if(wanted_tors[j]) */
          }                             /* end for(j = 0; j < tors_max; j++) */
          
          if(!k) wanted_tors[i] = 5;
     }                                  /* end for(i = 0; i < tors_max; i++) */

     if(purge_level) return;

     for(i = 0; i < tors_max; i++) {
          if(wanted_tors[i]) continue;
          if((back_order[(torsion_tab[i][0])] < 0)  ||
             (back_order[(torsion_tab[i][1])] < 0)  ||
             (back_order[(torsion_tab[i][2])] < 0)  ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          wanted_tors[i] = 6;
     }

     return;
}

/* count torsions
   parameter is:
        num_tors              total number of dihedrals to check
*/
void count_torsions(int num_tors)
{
     int i;

     regs_wanted = 0;
     regs_all = 0;
     rbs_wanted = 0;
     fourier_wanted = 0;
     fourier_all = 0;
     table_wanted = 0;
     table_all = 0;

     for(i = 0; i < num_tors; i++) {
          switch(a_tors_type[i]) {
               case 0:
               case 9:
               case 1:  if(wanted_tors[i])
                             regs_wanted++;
                        regs_all++;
                        break;

               case 2:  break;

               case 3:  if(wanted_tors[i])
                             rbs_wanted++;
                        rbs_all++;
                        break;

               case 5:  if(wanted_tors[i])
                             fourier_wanted++;
                        fourier_all++;
                        break;

               case 8:  if(wanted_tors[i])
                             table_wanted++;
                        table_all++;
                        break;

               default: fprintf(stdlog,
                                "Torsion type error, torsion number %d.\n",
                                i);
                        break;
          }                       /* end switch(a_tors_type[i]) */

     }                            /* end for(i = 0; i < num_tors; i++) */

     return;
}

/* finalize determination of planar opposition angles
   Are the 4 atoms bound to the central atom each in a planar ring?
   Is the central atom not listed as in a ring?
   parameter is:
        num_angles            total number of angles
*/
void planar_angles(int num_angles)
{
     int i, j, center, itest, l;

     for(i = 0; i < num_angles; i++) {
          if(opposite[i] != 'y')
               continue;
          center = angle_tab[i][1];
          if(arom_atomno[center][0] > 0) {
               opposite[i] = NULLCHAR;     /* Atom in ring */
               continue;
          }
          j = angle_tab[i][0];
          l = angle_tab[i][2];

/* both ends of angle in planar rings? */
          if((ar_set[j][1] == 0) &&
             (ar_set[j][2] == 0) &&
             (ar_set[j][3] == 0)) {
                    opposite[i] = NULLCHAR;
                    continue;           /* No! */
          }
          if((ar_set[l][1] == 0) &&
             (ar_set[l][2] == 0) &&
             (ar_set[l][3] == 0)) {
                    opposite[i] = NULLCHAR;
                    continue;           /* No! */
          }

          l = 0;

/* And bonded to two more planar rings? */
          for(j = 0; j < bond_count[center]; j++) {
               itest = connect[center][j];
               if(atom_atno[itest] != 7)
                    continue;
               if((ar_set[itest][1] != 0) ||
                  (ar_set[itest][2] != 0) ||
                  (ar_set[itest][3] != 0))
                         l++;
          }                       /* end for(j = 0; j < bond_count[center]; j++) */
          if(l < 4) {
               opposite[i] = NULLCHAR;     /* Not 4 planar N's */
               continue;
          }
     }                            /* end for(i = 0; i < num_angles; i++) */

     return;
}

