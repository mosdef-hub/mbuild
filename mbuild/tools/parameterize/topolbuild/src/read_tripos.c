/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to read tripos forcefield TAFF files read_tripos.c
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "amber_FF.h"
#include "Tripos_FF.h"

#define MAXCHAR 256
#define MAXNAME 128

extern char line[2*MAXCHAR];               /* input line */
extern char mess[MAXCHAR];                 /* for error messages */

/* read a tripos type forcefield parm file
   parameters are:
        dir                     directory with Tripos forcefield files to read
        just_count              flag to just count entries
        num_bind                number of bonds in forcefield
        num_angle               number of angles in forcefield
        num_tors                number of torsions in forcefield
        num_impro               number of impropers in forcefield
        num_vdw                 number of van der Waals parameters in forcefield
*/
void read_tripos(char *dir, int just_count, int *num_bind, int *num_angle,
                 int *num_tors, int *num_impro, int *num_vdw)
{
     int incase, i, j, k, oldnum_vdw;
     int bind_max, angle_max, tors_max, impro_max, vdw_max;
     int bnds, angls, mytors, impros, vdws;
     char filename[MAXCHAR];
     FILE *parmfp;

     bind_max = *num_bind;
     angle_max = *num_angle;
     tors_max = *num_tors;
     impro_max = *num_impro;
     vdw_max = *num_vdw;

     if(just_count) {
          i = *num_bind;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field bond array sizes\n");
               exit(1);
          }
          bond_name1 = (char **)mat_alloc(i, 7, FARGS);
          bond_name2 = (char **)mat_alloc(i, 7, FARGS);
          bond_force = (float *)allocator(i, sizeof(float), FARGS);
          bond_length = (float *)allocator(i, sizeof(float), FARGS);
          trp_bnd_type = (char **)mat_alloc(i, 6, FARGS);

          i = *num_angle;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field angle array sizes\n");
               exit(1);
          }
          angle_name1 = (char **)mat_alloc(i, 7, FARGS);
          angle_name2 = (char **)mat_alloc(i, 7, FARGS);
          angle_name3 = (char **)mat_alloc(i, 7, FARGS);
          angle_force = (float *)allocator(i, sizeof(float), FARGS);
          angle_angle = (float *)allocator(i, sizeof(float), FARGS);

          i = *num_tors;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field torsion array sizes\n");
               exit(1);
          }
          tors_name1 = (char **)mat_alloc(i, 7, FARGS);
          tors_name2 = (char **)mat_alloc(i, 7, FARGS);
          tors_name3 = (char **)mat_alloc(i, 7, FARGS);
          tors_name4 = (char **)mat_alloc(i, 7, FARGS);
          tors_force = (float *)allocator(i, sizeof(float), FARGS);
          tors_term = (float *)allocator(i, sizeof(float), FARGS);
          tors_bond = (char **)mat_alloc(i, 6, FARGS);

          i = *num_impro;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field improper array sizes\n");
               exit(1);
          }
          impro_name3 = (char **)mat_alloc(i, 7, FARGS);
          impro_force = (float *)allocator(i, sizeof(float), FARGS);

          i = *num_vdw;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field van der Waals array sizes\n");
               exit(1);
          }
          vdw_atom_name = (char **)mat_alloc(i, 7, FARGS);
          vdw_atom_radius = (float *)allocator(i, sizeof(float), FARGS);
          vdw_atom_pot = (float *)allocator(i, sizeof(float), FARGS);
     }

     bnds = 0;
     angls = 0;
     mytors = 0;
     impros = 0;
     vdws = 0;

/* read bond forcefield data */
     strcpy(filename, dir);
     strcat(filename, "/TAFF_BOND_STRETCH");

     if((parmfp = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          if(just_count) {
               if(bnds > bind_max) {
                    my_fatal(FARGS,
                      "Number of Tripos FF bonds records exceeds maximum\n");
                    exit(1);
               }
               bond_length[bnds] = 99999.0;
               bond_force[bnds] = -1.0;
               sscanf(line, "%s%s%s%f%f", bond_name1[bnds], bond_name2[bnds],
                      trp_bnd_type[bnds], &bond_length[bnds], &bond_force[bnds]);
          }                     /* end if(just_count) */
          bnds++;
     }                          /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     *num_bind = bnds;

/* read angles forcefield data */
     strcpy(filename, dir);
     strcat(filename, "/TAFF_ANGLE_BEND");

     if((parmfp = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          if(just_count) {
               if(angls > angle_max) {
                    my_fatal(FARGS,
                      "Number of Tripos FF angle records exceeds maximum\n");
                    exit(1);
               }
               angle_force[angls] = -1.0;
               angle_angle[angls] = 99999.0;
               sscanf(line, "%s%s%s%f%f", angle_name1[angls], angle_name2[angls],
                      angle_name3[angls], &angle_angle[angls], &angle_force[angls]);
               if(angle_force[angls] < 9999.0)     /* convert from deg^-2 to rad^-2 */
                    angle_force[angls] = angle_force[angls] * 3282.80635;
               if(!strcmp(angle_name1[angls], "*") ||
                  !strcmp(angle_name1[angls], "ANY"))
                         strcpy(angle_name1[angls], "X");
               if(!strcmp(angle_name2[angls], "*") ||
                  !strcmp(angle_name2[angls], "ANY"))
                         strcpy(angle_name2[angls], "X");
               if(!strcmp(angle_name3[angls], "*") ||
                  !strcmp(angle_name3[angls], "ANY"))
                         strcpy(angle_name3[angls], "X");
          }                     /* end if(just_count) */
          angls++;
     }                          /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     *num_angle = angls;

/* read torsions forcefield data */
     strcpy(filename, dir);
     strcat(filename, "/TAFF_TORS");

     if((parmfp = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          if(just_count) {
               if(mytors > tors_max) {
                    my_fatal(FARGS,
                      "Number of Tripos FF torsion records exceeds maximum\n");
                    exit(1);
               }
               tors_term[mytors] = 99999.0;
               tors_force[mytors] = -1.0;
               sscanf(line, "%s%s%s%s%s%f%f", tors_name1[mytors],
                      tors_name2[mytors], tors_name3[mytors],
                      tors_name4[mytors], tors_bond[mytors],
                      &tors_force[mytors], &tors_term[mytors]);
               if(!strcmp(tors_name1[mytors], "*") ||
                  !strcmp(tors_name1[mytors], "ANY"))
                         strcpy(tors_name1[mytors], "X");
               if(!strcmp(tors_name2[mytors], "*") ||
                  !strcmp(tors_name2[mytors], "ANY"))
                         strcpy(tors_name2[mytors], "X");
               if(!strcmp(tors_name3[mytors], "*") ||
                  !strcmp(tors_name3[mytors], "ANY"))
                         strcpy(tors_name3[mytors], "X");
               if(!strcmp(tors_name4[mytors], "*") ||
                  !strcmp(tors_name4[mytors], "ANY"))
                         strcpy(tors_name4[mytors], "X");
          }                     /* end if(just_count) */
          mytors++;
     }                          /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     *num_tors = mytors;

/* read van der Waals forcefield data */
     strcpy(filename, dir);
     strcat(filename, "/TAFF_VDW");

     if((parmfp = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          if(just_count) {
               if(vdws >= vdw_max) {
                    my_fatal(FARGS,
                      "Number of Tripos FF van der Waals records exceeds maximum\n");
                    exit(1);
               }
               vdw_atom_radius[vdws] = 99999.0;
               vdw_atom_pot[vdws] = 99999.0;
               sscanf(line, "%s%f%f", vdw_atom_name[vdws], &vdw_atom_radius[vdws],
                      &vdw_atom_pot[vdws]);
          }                     /* end if(just_count) */
          vdws++;
     }                          /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     *num_vdw = vdws;

/* read out of plane bending forcefield data */
     strcpy(filename, dir);
     strcat(filename, "/TAFF_OOP_BEND");

     if((parmfp = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          if(just_count) {
               if(impros >= impro_max) {
                    my_fatal(FARGS,
                      "Number of Tripos FF out of plane bending records exceeds maximum\n");
                    exit(1);
               }
               impro_force[impros] = -1.0;
               sscanf(line, "%s%f", impro_name3[impros], &impro_force[impros]);
          }                     /* end if(just_count) */
          impros++;
     }                          /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     *num_impro = impros;

     return;
}

/* free Tripos forcefield storage
   parameters are:
        num_bind           number of bonds in forcefield
        num_angle          number of angles in forcefield
        num_tors           number of torsions in forcefield
        num_impro          number of impropers in forcefield
        num_vdw            number of van der Waals parameters in forcefield
*/

void free_triposFF(int num_bind, int num_angle, int num_tors, int num_impro,
                   int num_vdw)
{
     int i;

     for(i = 0; i < num_bind; i++) {
          free_me(bond_name1[i], FARGS);
          free_me(bond_name2[i], FARGS);
          free_me(trp_bnd_type[i], FARGS);
     }
     free_me(bond_name1, FARGS);
     free_me(bond_name2, FARGS);
     free_me(trp_bnd_type, FARGS);
     free_me(bond_force, FARGS);
     free_me(bond_length, FARGS);

     for(i = 0; i < num_angle; i++) {
          free_me(angle_name1[i], FARGS);
          free_me(angle_name2[i], FARGS);
          free_me(angle_name3[i], FARGS);
     }
     free_me(angle_name1, FARGS);
     free_me(angle_name2, FARGS);
     free_me(angle_name3, FARGS);
     free_me(angle_force, FARGS);
     free_me(angle_angle, FARGS);

     for(i = 0; i < num_tors; i++) {
          free_me(tors_name1[i], FARGS);
          free_me(tors_name2[i], FARGS);
          free_me(tors_name3[i], FARGS);
          free_me(tors_name4[i], FARGS);
          free_me(tors_bond[i], FARGS);
     }
     free_me(tors_name1, FARGS);
     free_me(tors_name2, FARGS);
     free_me(tors_name3, FARGS);
     free_me(tors_name4, FARGS);
     free_me(tors_force, FARGS);
     free_me(tors_bond, FARGS);
     free_me(tors_term, FARGS);

     if(num_impro) {
          for(i = 0; i < num_impro; i++)
               free_me(impro_name3[i], FARGS);

          free_me(impro_name3, FARGS);
          free_me(impro_force, FARGS);
     }

     for(i = 0; i < num_vdw; i++)
          free_me(vdw_atom_name[i], FARGS);
     free_me(vdw_atom_name, FARGS);
     free_me(vdw_atom_radius, FARGS);
     free_me(vdw_atom_pot, FARGS);

     return;
}

