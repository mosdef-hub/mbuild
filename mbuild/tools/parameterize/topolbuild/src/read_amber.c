/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to read amber forcefield parm files read_amber.c

   NOTICE: This is a derivative work based on study of the routines
   found in the antechamber 1.27 programs parmchk and prepgen,
   version 1.0, dated October, 2001 by Junmei Wang, Department
   of Pharmaceutical Chemistry, School of Pharmacy, University
   of California, San Francisco, CA  94143

   Portions of this work simplify storage allocation, and clarify
   decision trees compared to the work cited above.
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include "block_memory.h"
#include "read_amber.h"
#include "amber_FF.h"

#define MAXNAME 128

extern char line[2*MAXCHAR];               /* input line */
extern char mess[MAXCHAR];                 /* for error messages */

int set_same_tors = 0;

/* read an amber type forcefield parm file
   parameters are:
        parmfile           file to read
        just_count         flag to just count entries
        num_mol            number of atoms in forcefield
        num_bind           number of bonds in forcefield
        num_angle          number of angles in forcefield
        num_tors           number of torsions in forcefield
        num_impro          number of impropers in forcefield
        num_vdw            number of van der Waals parameters in forcefield
        num_equiv          number of van der Waals equivalency records
        max_equiv          maximum length of van der Waals equivalency records
*/
void read_amber(char *parmfile, int just_count, int *num_mol, int *num_bind,
                int *num_angle, int *num_tors, int *num_impro, int *num_vdw,
                int *num_equiv, int *max_equiv)
{
     int incase, i, j, k, oldnum_vdw;
     int mol_max, bind_max, angle_max, tors_max, impro_max, vdw_max;
     int equiv_max, eqlen_max;
     int mols, bnds, angls, mytors, impros, vdws, equivs, lengths;
     FILE *parmfp;

     mol_max = *num_mol;
     bind_max = *num_bind;
     angle_max = *num_angle;
     tors_max = *num_tors;
     impro_max = *num_impro;
     vdw_max = *num_vdw;
     equiv_max = *num_equiv;
     eqlen_max = *max_equiv;

     if(just_count) {
          i = *num_mol;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field atom array sizes\n");
               exit(1);
          }
          ffatom_name = (char **)mat_alloc(i, 5, FARGS);
          ffatom_mass = (float *)allocator(i, sizeof(float), FARGS);
          ffatom_pol = (float *)allocator(i, sizeof(float), FARGS);

          i = *num_bind;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field bond array sizes\n");
               exit(1);
          }
          bond_name1 = (char **)mat_alloc(i, 5, FARGS);
          bond_name2 = (char **)mat_alloc(i, 5, FARGS);
          bond_force = (float *)allocator(i, sizeof(float), FARGS);
          bond_length = (float *)allocator(i, sizeof(float), FARGS);
          FF_bond_type = (short *)allocator(i, sizeof(short), FARGS);
          for(j = 0; j < i; j++)
               FF_bond_type[j] = 1;            /* always type 1 for gaff/amber */

          i = *num_angle;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field angle array sizes\n");
               exit(1);
          }
          angle_name1 = (char **)mat_alloc(i, 5, FARGS);
          angle_name2 = (char **)mat_alloc(i, 5, FARGS);
          angle_name3 = (char **)mat_alloc(i, 5, FARGS);
          angle_force = (float *)allocator(i, sizeof(float), FARGS);
          angle_angle = (float *)allocator(i, sizeof(float), FARGS);
          FF_angle_type = (short *)allocator(i, sizeof(short), FARGS);
          for(j = 0; j < i; j++)
               FF_angle_type[j] = 1;           /* always type 1 for gaff/amber */

          i = *num_tors;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field torsion array sizes\n");
               exit(1);
          }
          tors_name1 = (char **)mat_alloc(i, 5, FARGS);
          tors_name2 = (char **)mat_alloc(i, 5, FARGS);
          tors_name3 = (char **)mat_alloc(i, 5, FARGS);
          tors_name4 = (char **)mat_alloc(i, 5, FARGS);
          tors_mult = (int *)allocator(i, sizeof(float), FARGS);
          tors_force = (float *)allocator(i, sizeof(float), FARGS);
          tors_phase = (float *)allocator(i, sizeof(float), FARGS);
          tors_term = (float *)allocator(i, sizeof(float), FARGS);
          tors_more = (int *)allocator(i, sizeof(int), FARGS);
          FF_tors_type = (short *)allocator(i, sizeof(short), FARGS);
          for(j = 0; j < i; j++)
               FF_tors_type[j] = 1;            /* always type 1 for gaff/amber */

          i = *num_impro;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field improper array sizes\n");
               exit(1);
          }
          impro_name1 = (char **)mat_alloc(i, 5, FARGS);
          impro_name2 = (char **)mat_alloc(i, 5, FARGS);
          impro_name3 = (char **)mat_alloc(i, 5, FARGS);
          impro_name4 = (char **)mat_alloc(i, 5, FARGS);
          impro_mult = (int *)allocator(i, sizeof(float), FARGS);
          impro_xcount = (int *)allocator(i, sizeof(int), FARGS);
          impro_force = (float *)allocator(i, sizeof(float), FARGS);
          impro_phase = (float *)allocator(i, sizeof(float), FARGS);
          impro_term = (float *)allocator(i, sizeof(float), FARGS);

          i = *num_vdw;
          if(i < 1) {
               my_fatal(FARGS,
                 "Less than 1 specified for force field van der Waals array sizes\n");
               exit(1);
          }
          vdw_atom_name = (char **)mat_alloc(i, 5, FARGS);
          vdw_atom_radius = (float *)allocator(i, sizeof(float), FARGS);
          vdw_atom_pot = (float *)allocator(i, sizeof(float), FARGS);

          i = *num_equiv;
          if(i > 0) {                 /* allow for 0 equivalent names entries */
               equiv_name = (char ***)tabl_alloc(i, eqlen_max, 5, FARGS);
               equiv_count = (int *)allocator(i, sizeof(float), FARGS);
          }
     }

     incase = 0;                        /* switch flag */
     mols = 0;
     bnds = 0;
     angls = 0;
     mytors = 0;
     impros = 0;
     vdws = 0;
     equivs = 0;
     lengths = 0;
     set_same_tors = 0;

     if((parmfp = fopen(parmfile, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", parmfile);
          my_fatal(FARGS, mess);
          exit(1);
     }

     while(fgets(line, MAXCHAR, parmfp) != NULL) {
          switch(incase) {
               case 0:  incase = 1;
                        break;            /* first line is skipped */

               case 1:  if(spaceline(line)) {
                             incase = 2;
                             break;
                        }
                        read_molparms(just_count, mol_max, &mols);
                        break;

               case 2:  incase = 3;
                        break;            /* a line is skipped */

               case 3:  if(spaceline(line)) {
                             incase = 4;
                             break;
                        }
                        read_bindparms(just_count, bind_max, &bnds);
                        break;

               case 4:  if(spaceline(line)) {
                             incase = 5;
                             break;
                        }
                        read_angleparms(just_count, angle_max, &angls);
                        break;

               case 5:  if(spaceline(line)) {
                             incase = 6;
                             break;
                        }
                        read_torsparms(just_count, tors_max, &mytors);
                        break;

               case 6:  if(spaceline(line)) {
                             incase = 7;
                             break;
                        }
                        read_improperparms(just_count, impro_max, &impros);
                        break;

               case 7: if(spaceline(line)) {
                             incase = 8;
                             break;
                        }
                        break;

               case 8:  if(spaceline(line)) {
                             incase = 9;
                             break;
                        }
                        read_vdwequivs(just_count, equiv_max, eqlen_max,
                                       &equivs, &lengths);
                        break;

               case 9:  if(!strncmp(line, "MOD4", 4)) {
                             incase = 10;
                             break;
                        }

               case 10: if(spaceline(line)) {
                            fclose(parmfp);           /* normal exit */
                            if(equivs) {
                                 oldnum_vdw = vdws;

                                 for(i = 0; i < equivs; i++)
                                      for(k = 0; k < oldnum_vdw; k++)
                                           if(just_count)
                                                if(!strcmp(vdw_atom_name[k],
                                                           equiv_name[i][0])) {
                                                     for(j = 1; j < equiv_count[i]; j++) {
                                                          if(vdws >= vdw_max) {
                                                               my_fatal(FARGS,
                  "With equivalents, number of Amber FF van der Waals records exceeds maximum\n");
                                                               exit(1);
                                                          }
                                                          if(strlen(equiv_name[i][j]) < 1)
                                                               continue;
                                                          strcpy(vdw_atom_name[vdws],
                                                                        equiv_name[i][j]);
                                                          vdw_atom_radius[vdws] =
                                                                        vdw_atom_radius[k];
                                                          vdw_atom_pot[vdws] = 
                                                                           vdw_atom_pot[k];
                                                          vdws++;

                                                     } /* end for(j = 1;   etc. */
                                                }      /* end if(!strcmp(vdw_atom_name[k],
                                                                         equiv_name[i][0])) */
                            }                          /* end if(num_equiv) */

                            *num_mol = mols;
                            *num_bind = bnds;
                            *num_angle = angls;
                            *num_tors = mytors;
                            *num_impro = impros;
                            *num_vdw = vdws;
                            *num_equiv = equivs;
                            *max_equiv = lengths;
                            return;                    /* DONE! we read all we need */
                            break;
                       }
                       read_vdwparms(just_count, vdw_max, &vdws);
                       break;

               default: break;
          }                    /* end switch(incase) */
     }                         /* end while(fgets(line, MAXCHAR, parmfp) != NULL) */

     fclose(parmfp);                /* we have an error exit */
     sprintf(mess,
             "Premature termination of Amber FF file %s\n",
             parmfile);
     my_fatal(FARGS, mess);
     exit(1);
}

/* read atomic parameters
   parameters are:
        just_count         flag to just count entries
        mol_max            maximum number of atoms to read
        num_mol            count of atoms in forcefield to this point
*/
void read_molparms(int just_count, int mol_max, int *num_mol)
{
     int cntr;

     cntr = *num_mol;
     if(just_count) {
          if(cntr >= mol_max) {
               my_fatal(FARGS,
                 "Number of Amber FF atom records exceeds maximum\n");
               exit(1);
          }
          strncpy(&ffatom_name[cntr][0], &line[0], 2);
          ffatom_name[cntr][2] = NULLCHAR;
          ffatom_mass[cntr] = -1.0;
          ffatom_pol[cntr] = -1.0;
          sscanf(&line[3], "%f %f", &ffatom_mass[cntr],
                 &ffatom_pol[cntr]);
     }

     cntr++;
     *num_mol = cntr;
     return;
}

/* read bond parameters
   parameters are:
        just_count         flag to just count entries
        bind_max           maximum number of bond types to read
        num_bind           count of bond types in forcefield to this point
*/
void read_bindparms(int just_count, int bind_max, int *num_bind)
{
     int cntr;

     cntr = *num_bind;
     if(just_count) {
          if(cntr >= bind_max) {
               my_fatal(FARGS,
                 "Number of Amber FF bonds records exceeds maximum\n");
               exit(1);
          }
          strncpy(&bond_name1[cntr][0], &line[0], 2);
          strncpy(&bond_name2[cntr][0], &line[3], 2);
          bond_name1[cntr][2] = NULLCHAR;
          bond_name2[cntr][2] = NULLCHAR;
          bond_force[cntr] = -1.0;
          bond_length[cntr] = 99999.0;
          sscanf(&line[5], "%f%f", &bond_force[cntr],
                                   &bond_length[cntr]);
     }

     cntr++;
     *num_bind = cntr;
     return;
}

/* read angle parameters
   parameters are:
        just_count         flag to just count entries
        angle_max          maximum number of angle types to read
        num_angle          count of angle types in forcefield to this point
*/
void read_angleparms(int just_count, int angle_max, int *num_angle)
{
     int cntr;

     cntr = *num_angle;
     if(just_count) {
          if(cntr >= angle_max) {
               my_fatal(FARGS,
                 "Number of Amber FF angles records exceeds maximum\n");
               exit(1);
          }
          strncpy(&angle_name1[cntr][0], &line[0], 2);
          strncpy(&angle_name2[cntr][0], &line[3], 2);
          strncpy(&angle_name3[cntr][0], &line[6], 2);
          angle_name1[cntr][2] = NULLCHAR;
          angle_name2[cntr][2] = NULLCHAR;
          angle_name3[cntr][2] = NULLCHAR;
          angle_force[cntr] = -1.0;
          angle_angle[cntr] = 99999.0;
          sscanf(&line[8], "%f%f", &angle_force[cntr],
                                   &angle_angle[cntr]);
     }

     cntr++;
     *num_angle = cntr;
     return;
}

/* read torsion parameters
   parameters are:
        just_count         flag to just count entries
        tors_max           maximum number of torsion types to read
        num_tors           count of torsion types in forcefield to this point
*/
void read_torsparms(int just_count, int tors_max, int *num_tors)
{
     int cntr, i;

     cntr = *num_tors;
     i = 0;
     if(just_count) {
          if(cntr >= tors_max) {
               my_fatal(FARGS,
                 "Number of Amber FF torsions records exceeds maximum\n");
               exit(1);
          }

/* Deal with cases where the following lines are further torsion entries for
   the same atom type by setting the atom type to that of the prior line, and
   counting the number of entries.  This situation is found with the glycam
   parameters.
*/
          if(cntr)
               if(!strncmp(line, "           ", 11)) {
                    strcpy(tors_name1[cntr], tors_name1[(cntr - 1)]);
                    strcpy(tors_name2[cntr], tors_name2[(cntr - 1)]);
                    strcpy(tors_name3[cntr], tors_name3[(cntr - 1)]);
                    strcpy(tors_name4[cntr], tors_name4[(cntr - 1)]);
                    i = 1;             /* flag that the name is set */
                    flag_more_tors(cntr);
               }

          if(!i) {                     /* read the name if it is not set */
               strncpy(&tors_name1[cntr][0], &line[0], 2);
               strncpy(&tors_name2[cntr][0], &line[3], 2);
               strncpy(&tors_name3[cntr][0], &line[6], 2);
               strncpy(&tors_name4[cntr][0], &line[9], 2);
               tors_name1[cntr][2] = NULLCHAR;
               tors_name2[cntr][2] = NULLCHAR;
               tors_name3[cntr][2] = NULLCHAR;
               tors_name4[cntr][2] = NULLCHAR;
               if(set_same_tors) {
                    if(!strcmp(tors_name1[cntr], tors_name1[(cntr - 1)]) &&
                       !strcmp(tors_name2[cntr], tors_name2[(cntr - 1)]) &&
                       !strcmp(tors_name3[cntr], tors_name3[(cntr - 1)]) &&
                       !strcmp(tors_name4[cntr], tors_name4[(cntr - 1)])) 
                              flag_more_tors(cntr);
                    else set_same_tors = 0;
               }
          }
          tors_mult[cntr] = 9999;
          tors_force[cntr] = -1.0;
          tors_phase[cntr] = 99999.0;
          tors_term[cntr] = 99999.0;
          tors_more[cntr] = NULLCHAR;
          sscanf(&line[11], "%d%f%f%f", &tors_mult[cntr],
                                        &tors_force[cntr],
                                        &tors_phase[cntr],
                                        &tors_term[cntr]);

          if(tors_term[cntr] < 0.0) set_same_tors = 1;
     }

     cntr++;
     *num_tors = cntr;
     return;
}

/* flag that there are more torsions in the list
   parameter is
       cntr           count of torsion types in forcefield to this point
*/
void flag_more_tors(int cntr)
{
     int i;

     tors_more[(cntr - 1)]++;
     if((cntr - 2) > -1)
          for(i = (cntr - 2); i > -1; i--) {
               if(!tors_more[i]) return;
               if(tors_more[i] > 3) {
                    my_fatal(FARGS,
                      "Too many torsion entries for same atom type\n");
                    exit(1);
               }
               tors_more[i]++;
          }

     return;
}

/* read improper parameters
   parameters are:
        just_count              flag to just count entries
        impro_max               maximum number of impropers types to read
        num_impro               count of impropers types in forcefield to this point
*/
void read_improperparms(int just_count, int impro_max, int *num_impro)
{
     int cntr;

     cntr = *num_impro;
     int tmpnum;
     char tmpry[5];

     tmpnum = 0;
     if(just_count) {
          if(cntr >= impro_max) {
               my_fatal(FARGS,
                 "Number of Amber FF impropers records exceeds maximum\n");
               exit(1);
          }
          if(line[0] == 'X') tmpnum++;
          if(line[3] == 'X') tmpnum++;
          if(line[6] == 'X') tmpnum++;
          if(line[9] == 'X') tmpnum++;
          strncpy(&impro_name1[cntr][0], &line[0], 2);
          strncpy(&impro_name2[cntr][0], &line[3], 2);
          strncpy(&impro_name3[cntr][0], &line[6], 2);
          strncpy(&impro_name4[cntr][0], &line[9], 2);
          impro_name1[cntr][2] = NULLCHAR;
          impro_name2[cntr][2] = NULLCHAR;
          impro_name3[cntr][2] = NULLCHAR;
          impro_name4[cntr][2] = NULLCHAR;

          if(strcmp(&impro_name1[cntr][0], &impro_name2[cntr][0]) > 0) {
               strcpy(tmpry, &impro_name2[cntr][0]);
               strcpy(&impro_name2[cntr][0], &impro_name1[cntr][0]);
               strcpy(&impro_name1[cntr][0], tmpry);
          }
          if(strcmp(&impro_name1[cntr][0], &impro_name4[cntr][0]) > 0) {
               strcpy(tmpry, &impro_name4[cntr][0]);
               strcpy(&impro_name4[cntr][0], &impro_name1[cntr][0]);
               strcpy(&impro_name1[cntr][0], tmpry);
          }
          if(strcmp(&impro_name2[cntr][0], &impro_name4[cntr][0]) > 0) {
               strcpy(tmpry, &impro_name4[cntr][0]);
               strcpy(&impro_name4[cntr][0], &impro_name2[cntr][0]);
               strcpy(&impro_name2[cntr][0], tmpry);
          }

          impro_force[cntr] = -1.0;
          impro_phase[cntr] = 99999.0;
          impro_term[cntr] = 99999.0;
          sscanf(&line[11], "%f%f%f", &impro_force[cntr],
                                      &impro_phase[cntr],
                                      &impro_term[cntr]);
          impro_mult[cntr] = 1;
          impro_xcount[cntr] = tmpnum;
     }

     cntr++;
     *num_impro = cntr;
     return;
}

/* read van der Waals parameters
   parameters are:
        just_count         flag to just count entries
        vdw_max            maximum number of van der Waals types to read
        num_vdw            count of van der Waals types in forcefield to
                           this point
*/
void read_vdwparms(int just_count, int vdw_max, int *num_vdw)
{
     int cntr;

     cntr = *num_vdw;
     if(just_count) {
          if(cntr >= vdw_max) {
               my_fatal(FARGS,
                 "Number of Amber FF van der Waals records exceeds maximum\n");
               exit(1);
          }
          vdw_atom_radius[cntr] = 99999.0;
          vdw_atom_pot[cntr] = 99999.0;
          sscanf(line, "%s%f%f", vdw_atom_name[cntr], 
                                 &vdw_atom_radius[cntr],
                                 &vdw_atom_pot[cntr]);

          if(strlen(vdw_atom_name[cntr]) < 2) {
               vdw_atom_name[cntr][1] = ' ';
               vdw_atom_name[cntr][2] = NULLCHAR;
          }
     }

     cntr++;
     *num_vdw = cntr;
     return;
}

/* read van der Waals equivalency records
   parameters are:
        just_count         flag to just count entries
        equiv_max          maximum number of van der Waals equivalency records
                           to read
        eqlen_max          maximum length of any van der Waals equivalency records
                           in forcefield
        num_equiv          count of van der Waals equivalency records in forcefield
                           to this point
        max_equiv          maximum length of van der Waals equivalency records in
                           forcefield to this point
*/
void read_vdwequivs(int just_count, int equiv_max, int eqlen_max, int *num_equiv,
                    int *max_equiv)
{
     char tmpry[5];
     int i, eqlen, point, cntr;

     cntr = *num_equiv;
     eqlen = 0;
     point = 0;

     if(just_count)
          if(cntr >= equiv_max) {
               my_fatal(FARGS,
                 "Number of Amber FF van der Waals equivalence records exceeds maximum\n");
               exit(1);
          }

     while((point < strlen(line)) &&
           sscanf(&line[point], "%2s", tmpry)) {      /* did we read a name? */
               if(strlen(tmpry) == 1) {
                    tmpry[1] = ' ';
                    tmpry[2] = NULLCHAR;
               }

          if(just_count)
               strcpy(&equiv_name[cntr][eqlen][0],
                      tmpry);     /* save the name */
          if(strlen(tmpry)) eqlen++;                       /* count it */
          for(i = point; i < strlen(line); i++)        /* work way past name */
               if(!strncmp(&line[i], tmpry, 2)) {      /* name we just read must be
                                                          here somewhere */
                    point = i + 2;                     /* found it. pass it and read on */
                    break;                                     
               }
               else point++;
          if(spaceline(&line[point])) break;           /* rest of line blank? */
     }

/* no more names to read.  prepare to exit */

     if(eqlen < 2) return;                             /* if did not get any names, exit
                                                          without updating counts */

     if(just_count) equiv_count[cntr] = eqlen;      /* store names gotten count */

     if(eqlen > *max_equiv) *max_equiv = eqlen;        /* update maximum length */
     cntr++;
     *num_equiv = cntr;

     return;
}

/* frees everything except the equivalent van der Waals lists, which
   ought to have been freed immediately after their use.
   parameters are:
        num_mol            number of atoms in forcefield
        num_bind           number of bonds in forcefield
        num_angle          number of angles in forcefield
        num_tors           number of torsions in forcefield
        num_impro          number of impropers in forcefield
        num_vdw            number of van der Waals parameters in
                           forcefield
*/

void free_amberFF(int num_mol, int num_bind, int num_angle,
                  int num_tors, int num_impro, int num_vdw)
{
     int i;
     for(i = 0; i < num_mol; i++)
          free_me(ffatom_name[i], FARGS);
     free_me(ffatom_name, FARGS);
     free_me(ffatom_mass, FARGS);
     free_me(ffatom_pol, FARGS);

     for(i = 0; i < num_bind; i++) {
          free_me(bond_name1[i], FARGS);
          free_me(bond_name2[i], FARGS);
     }
     free_me(bond_name1, FARGS);
     free_me(bond_name2, FARGS);
     free_me(bond_force, FARGS);
     free_me(bond_length, FARGS);
     if(bond_beta != NULL)
          free_me(bond_beta, FARGS);
     free_me(FF_bond_type, FARGS);

     for(i = 0; i < num_angle; i++) {
          free_me(angle_name1[i], FARGS);
          free_me(angle_name2[i], FARGS);
          free_me(angle_name3[i], FARGS);
     }
     free_me(angle_name1, FARGS);
     free_me(angle_name2, FARGS);
     free_me(angle_name3, FARGS);
     if(angle_quartic != NULL) {
          for(i = 0; i < 5; i++)
               free_me(angle_quartic[i], FARGS);
          free_me(angle_quartic, FARGS);
          angle_force = NULL;
     }
     else free_me(angle_force, FARGS);
     free_me(angle_angle, FARGS);
     free_me(FF_angle_type, FARGS);

     for(i = 0; i < num_tors; i++) {
          free_me(tors_name1[i], FARGS);
          free_me(tors_name2[i], FARGS);
          free_me(tors_name3[i], FARGS);
          free_me(tors_name4[i], FARGS);
     }
     free_me(tors_name1, FARGS);
     free_me(tors_name2, FARGS);
     free_me(tors_name3, FARGS);
     free_me(tors_name4, FARGS);
     if(tors_RB != NULL) {
          for(i = 0; i < 6; i++)
               free_me(tors_RB[i], FARGS);
          free_me(tors_RB, FARGS);
          tors_mult = NULL;
          tors_force = NULL;
          tors_phase = NULL;
          tors_term = NULL;
          tors_more = NULL;
     }
     else {
          free_me(tors_mult, FARGS);
          free_me(tors_force, FARGS);
          free_me(tors_phase, FARGS);
          free_me(tors_term, FARGS);
          free_me(tors_more, FARGS);
     }
     free_me(FF_tors_type, FARGS);

     if(num_impro) {
          for(i = 0; i < num_impro; i++) {
               free_me(impro_name1[i], FARGS);
               free_me(impro_name2[i], FARGS);
               free_me(impro_name3[i], FARGS);
               free_me(impro_name4[i], FARGS);
          }
          free_me(impro_name1, FARGS);
          free_me(impro_name2, FARGS);
          free_me(impro_name3, FARGS);
          free_me(impro_name4, FARGS);
          free_me(impro_mult, FARGS);
          free_me(impro_xcount, FARGS);
          free_me(impro_force, FARGS);
          free_me(impro_phase, FARGS);
          free_me(impro_term, FARGS);
     }

     for(i = 0; i < num_vdw; i++)
          free_me(vdw_atom_name[i], FARGS);
     free_me(vdw_atom_name, FARGS);
     free_me(vdw_atom_radius, FARGS);
     free_me(vdw_atom_pot, FARGS);

     return;
}

/* free van der Waals equivalents records lists
   parameters are:
        num_equiv          number of van der Waals equivalency records
        max_equiv          maximum length of van der Waals equivalency records
*/
void free_equivs(int num_equiv, int max_equiv)
{
     int i, j;

     if(!num_equiv) return;               /* allow for 0 equivalent names entries */

     for(i = 0; i < num_equiv; i++) {     /* don't need the equivalents data any more */
          for(j = 0; j < max_equiv; j++)
               free_me(equiv_name[i][j], FARGS);
          free_me(equiv_name[i], FARGS);
     }
     free_me(equiv_name, FARGS);
     free_me(equiv_count, FARGS);

     return;
}
