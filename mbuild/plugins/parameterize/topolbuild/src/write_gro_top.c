/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Function to write forcefield include files for a topology

   NOTICE: Portions of this are derivative work based on study of the
   routines found GROMACS version 3.3.1 written by David van der Spoel,
   Erik Lindahl, Berk Hess, and others and copyright (c) 1991-2000,
   University of Groningen, The Netherlands.  Copyright (c) 2001-2004,
   The GROMACS development team
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "gromacs_FF.h"
#include "mol2.h"

#define MAXCHAR 256
#define NULLCHAR 0

extern char mess[MAXCHAR];
extern char *use_define[6];
extern char *in_command;
extern char vers_no[];

extern int lstart2_at;

/* write Gromacs topology include files
   parameters are:
        commonname                common name string for these files
        assigned_types            list of assigned types for atoms
        numatom                   number of atoms
        vdw_rad                   usually a van derWaals radius
                                  always parameter for calculating either sigma or a Cg term
        vdw_pot                   usually a van derWaals potential
                                  always parameter for calculating either sigma or a Cg term
        amber                     flag for calculation from amber derived parameters
        altnb                     either NULL or a file name for a nonbonded parameters file
*/
void write_ff(char *commonname, char **assigned_types, int numatom,
              float *vdw_rad, float *vdw_pot, int amber, char *altnb)
{
     char filename[MAXCHAR], file2[MAXCHAR], file3[MAXCHAR];
     char file4[MAXCHAR];
     char **unique;
     int *indx_unq;
     int i, j, num_unq, is_unq, count;
     float zero = 0.0;
     double myvdwr, sigma, epsilon;
     double conv = 4.18400000;

     FILE *nb_itp, *itp, *atp, *posre;

     strcpy(filename, "ff");
     strcat(filename, commonname);
     strcat(filename, ".itp");

     if(altnb == NULL) {
          strcpy(file2, "ff");
          strcat(file2, commonname);
          strcat(file2, "nb.itp");
     } else strcpy(file2, altnb);

     strcpy(file4, "posre");
     strcat(file4, commonname);
     strcat(file4, ".itp");

     if((itp = fopen(filename, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     fprintf(itp, "; topolbuild %s\n", vers_no);
     fprintf(itp, "; Command line:\n;     %s\n;\n",
             in_command);
     fprintf(itp, "#define %s\n", use_define[0]);
     fprintf(itp, "#define %s\n", use_define[1]);
     fprintf(itp, "#define _FF_USER\n\n[ defaults ]\n");
     fprintf(itp, ";nbfunc     comb-rule      gen-pairs     fudgeLJ      fudgeQQ\n");
     fprintf(itp, "%s\n\n", use_define[5]);
     fprintf(itp, "#include \"%s\"\n", file2);
     if(altnb == NULL) 
          fprintf(itp, "#include \"ffusernb.itp\"\n");

     fclose(itp);

     if((posre = fopen(file4, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", file4);
          my_fatal(FARGS, mess);
          exit(1);
     }

     fprintf(posre, "; topolbuild %s\n", vers_no);
     fprintf(posre, "; Command line:\n;     %s\n;\n",
             in_command);
     fprintf(posre, "; In this topology include file, you will find position restraint\n");
     fprintf(posre, "; entries for all the heavy atoms in your original mol2 file.\n");
     fprintf(posre, "; This means that protons are not restrained.\n\n");
     fprintf(posre, "[ position_restraints ]\n");
     fprintf(posre, "; atom  type      fx      fy      fz\n");

     count = -1;
     for(i = 0; i < numatom; i++) {
          if(atom_order[i] < 0) continue;
          count++;
          if(atom_atno[(atom_order[i])] < 2) continue;
          fprintf(posre, "%6d     1  1000  1000  1000\n", (count + lstart2_at));
     }

     fclose(posre);

     if(altnb != NULL) return;

     if((nb_itp = fopen(file2, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", file2);
          my_fatal(FARGS, mess);
          exit(1);
     }

     unique = (char**)mat_alloc(numatom, 12, FARGS);
     indx_unq = (int *)allocator(numatom, sizeof(int), FARGS);

     num_unq = 0;
     unique[0][0] = NULLCHAR;
     for(i = 0; i < numatom; i++) {
          if(back_order[i] < 0) continue;
          if((use_define[4] != NULL) &&
             (isdigit(assigned_types[(atom_order[i])][0])))
                    strcpy(unique[0], use_define[4]);
          strcat(unique[0], assigned_types[i]);              /* seed unique with first type */
          indx_unq[0] = i;
          num_unq = 1;
          break;
     }

     if(!num_unq) {
          my_fatal(FARGS, "No atoms with nonzero mass in molecule\n");
          exit(1);
     }

     for(i = 1; i < numatom; i++) {
          if(back_order[i] < 0) continue;
          is_unq = 1;
          for(j = 0; j < num_unq; j++)
               if(!strcmp(unique[j], assigned_types[i])) {
                    is_unq = 0;
                    break;
               }                                        /* end if(!strcmp(unique[j], assigned_types[i]) */

          if(is_unq) {
               unique[num_unq][0] = NULLCHAR;
               if((use_define[4] != NULL) &&
                  (isdigit(assigned_types[(atom_order[i])][0])))
                         strcpy(unique[num_unq], use_define[4]);

               strcat(unique[num_unq], assigned_types[i]);
               indx_unq[num_unq] = i;
               num_unq++;
          }                                             /* end if(is_unq) */
     }                                                  /* end for(i = 1; i < numatom; i++) */

     fprintf(nb_itp, "; topolbuild %s\n", vers_no);
     fprintf(nb_itp, "; Command line:\n;     %s\n;\n",
             in_command);
     fprintf(nb_itp, "[ atomtypes ]\n");
     fprintf(nb_itp, "; name        mass      charge    ptype               sigma             epsilon\n");

     for(i = 0; i < num_unq; i++) {
          myvdwr = vdw_rad[(indx_unq[i])];
          epsilon = vdw_pot[(indx_unq[i])];

          if(amber) myvdwr = myvdwr * 0.1;
          if(amber) epsilon = epsilon * conv;

          if(amber) sigma = (myvdwr / 1.12246) * 2.0;
          else sigma = vdw_rad[(indx_unq[i])];

          if(!bond_count[(indx_unq[i])])
               fprintf(nb_itp, "; ");

          fprintf(nb_itp, "%-9s%10.5f%12.5f%9s", unique[i],
                  atom_mass[(indx_unq[i])], zero, "A");

          if(vdw_rad[(indx_unq[i])] > 9999.0)
               fprintf(nb_itp,
                       "                                        ; Needs attention\n");
          else
               fprintf(nb_itp, "%20.6e%20.6e ; %8.2f%8.4f\n",
                       sigma, epsilon, vdw_rad[(indx_unq[i])],
                       vdw_pot[(indx_unq[i])]);
     }

     fclose(nb_itp);

     free_me(indx_unq, FARGS);
     for(i = 0; i < numatom; i++)
          free_me(unique[i], FARGS);
     free_me(unique, FARGS);

     return;
}
