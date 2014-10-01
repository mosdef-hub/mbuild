/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

        Functions to read and use gromacs force field data
        gromacs_use.c
        Portions of this work are copied from the files read_amber.c,
        and the do_an_amber function in the file topolbuild.c
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "gromacs_use.h"
#include "utilities.h"

extern char mess[MAXCHAR];
extern char ambername[MAXCHAR];
extern char mycorr_name[MAXCHAR];
extern char line[(2 * MAXCHAR)];

/* code to use a gromacs atoms force field, including united atoms model
   options
   parameters are:
        ante_name        directory atom type and parameter data files
        lev1_name        name of level 1 atom type file
        num_atoms        number of atoms in molecule
        num_bonds        number of bonds in molecule
        amberff          specific name of force field to use
        angle_max        number of angles in molecule
        tors_max         number of torsions in molecule
        ua_model         flag for use of united atoms model
        residue          residue name or NULLCHAR if do not have
                         common residue.
        by_corr          flag for later use
*/

void do_gmx_FFs(char *ante_name, char *lev1_name, int num_atoms,
                int num_bonds, char *amberff, int angle_max,
                int tors_max, int ua_model, char *residue, int by_corr)
{
     int num, num2, count1, count2, i;
     int FFmol, FFbond, FFangle, FFtors, FFimpro, FFvdw;

/* get a first level gromacs atom types assignment  Needed in order to
   properly sort out the various parameter differences for the same
   standard gromacs atom type in similar situations, but in different
   amino acid residues. Use a modified gaff or amber atom types
   definition table for this step.
*/

     strcpy(ambername, ante_name);
     strcat(ambername, lev1_name);

     read_at_types(&num, &num2, (4 * MAXCHAR), MAXCHAR,
                   0, ambername);                   /* first count lines */

     fprintf(stdlog,
             "\nIn %s there are %d WILDATOM lines and %d ATD lines\n",
             ambername, num2, num);
     count1 = num;
     count2 = num2;
     read_at_types(&num, &num2, count1, count2, 1,
                   ambername);                      /* now read them */

     init_atom_types(num_atoms, num_bonds);
     bond_info(num_atoms, num_bonds);

     judge_atom_types(count1, count2, num_atoms, num_bonds, a_ff_type);
     adjustments(num_atoms, num_bonds, a_ff_type);
     adjust_type_cp(num_atoms, num_bonds, a_ff_type);
     check_type_errors(num_atoms, num_bonds, a_ff_type);

     free_atom_types(num_atoms);
     release_types(count1, count2);

/* now have the first level of gromacs type assignments and can work
   on the force field entries and assignments
*/

     strcpy(ambername, ante_name);
     strcat(ambername, "/ff");
     strcat(ambername, amberff);
     strcat(ambername, ".dat");

/* This is my own forcefield format and has the counts of each type of
   before entry as part of the file.  Therefore I don't have to count
   anything before allocating memory
*/

     read_ffgmx1(ambername, &FFmol, &FFbond, &FFangle, &FFtors,
                 &FFimpro, &FFvdw);

     fprintf(stdlog,
             "\n\tFile has\n\t\t%d atom records\n\t\t%d bond records\n",
             FFmol, FFbond);
     fprintf(stdlog, "\t\t%d angle records\n", FFangle);
     fprintf(stdlog,
             "\t\t%d torsion records\n\t\t%d impropers records\n",
             FFtors, FFimpro);
     fprintf(stdlog, "\t\t%d van der Waals records\n", FFvdw);

/* I arranged my format so that I can use the same matching routines that
   I use with amber force field parameters
*/

     no_constraints = 0;
     set_amber_bondFF(num_bonds, FFbond);
     if(constrain)
          set_constraints(angle_max);
     set_amber_angleFF(angle_max, FFangle);
     set_gromacs_torsFF(tors_max, FFtors);
     if(improper_num)
          set_amber_improFF(improper_num, FFimpro);

     strcpy(mycorr_name, ante_name);
     strcat(mycorr_name, "/ATOMTYPE_");
     strcat(mycorr_name, amberff);
     strcat(mycorr_name, ".DEF");

     read_at_types(&num, &num2, (4 * MAXCHAR), MAXCHAR,
                   0, mycorr_name);            /* first count lines */

     fprintf(stdlog,
             "\nIn %s there are %d WILDATOM lines and %d ATD lines\n",
             mycorr_name, num2, num);
     count1 = num;
     count2 = num2;
     read_at_types(&num, &num2, count1, count2, 1,
                   mycorr_name);               /* now read them */

     init_atom_types(num_atoms, num_bonds);
     bond_info(num_atoms, num_bonds);

     judge_atom_types(count1, count2, num_atoms, num_bonds, a_ff_type);

/* not a gaff type atom type definition file.  gaff atom type
   corrections won't work here and aren't done.
*/
     free_atom_types(num_atoms);
     release_types(count1, count2);

     if(charge_set) {
          set_charge(num_atoms, FFmol);     /* set force field based charges */
          use_charge = 1;
     }

     if(FFvdw > 2)
          set_amber_vdwFF(num_atoms, FFvdw);
     set_amber_atomFF(num_atoms, FFmol);

     if(ua_model)                            /* convert to united atoms model */
          united_atoms(num_atoms, residue, FFmol);   /* if called for */

     free_me(gmx_tors_resid, FARGS);
     free_me(gmx_tors_num, FARGS);
     free_amberFF(FFmol, FFbond, FFangle, FFtors, FFimpro, FFvdw);

     if(constrain) {
          free_me(constrain_name1, FARGS);
          free_me(constrain_name2, FARGS);
          free_me(constrain_dist, FARGS);
          free_me(type_constr, FARGS);
     }

     free_me(macro_types, FARGS);
     free_me(macro_name, FARGS);
     free_me(macro_nmbr, FARGS);

     return;
}

/* convert C to united atoms model
   parameter are:
        num_atoms          number of atoms
        resid              residue name or NULLCHAR if do not have
                           common residue.
*/
void united_atoms(int num_atoms, char *resid, int nmol)
{
     int i, j, k, l, any, count;
     char mixer[6], *c, *d;

     any = 0;
     for(i = 0; i < num_atoms; i++) {
          if(atom_atno[i] != 6) continue;      /* Only looking for C */
          if(!bond_count[i]) continue;         /* Only C with bonds to
                                                  something else */
/* check if mass includes attached H's.  If not, can skip. */
          k = 0;
          for(j = 0; j < nmol; j++) {
               if(!strcmp(a_ff_type[i], ffatom_name[j])) {
                    if(ffatom_mass[j] < 13.0) {
                         k = 1;                /* Mass says keep H's */
                         break;
                    }
                    else break;                /* Mass says remove H's */
               }      /* end if(!strcmp(a_ff_type[i], ffatom_name[j])) */
          }           /* end for(j = 0; j < nmol; j++) */
          if(k) continue;

          for(j = 0; j < bond_count[i]; j++) {
               if(atom_atno[(connect[i][j])] != 1)
                    continue;   /* if not H, skip it */
/* add charge of the H to the C */
               atom_charge[i] += atom_charge[(connect[i][j])];
/* add mass of the H to the C */
               atom_mass[i] += atom_mass[(connect[i][j])];
/* remove mass from the H */
               atom_mass[(connect[i][j])] = 0.0;
               any = 1;               /* yes, I did one */
          }                           /* end for(j = 0; j < bond_count[i]; j++) */
     }                                /* end for(i = 0; i < num_atoms; i++) */

     if(!any) return;                 /* Only need to do more if changes made */

     for(i = 0; i < num_atoms; i++) {
          atom_order[i] = -1;
          back_order[i] = -1;
     }

     count = 0;
     if(resid[0] != NULLCHAR) {
          for(i = 0; i < num_atoms; i++) {
               if(atom_mass[i] == 0.0) {
                    atom_order[i] = -1;
                    back_order[i] = -1;
                    continue;
               }
               atom_order[i] = i;
               back_order[i] = count;
               count++;
          }                /* end for(i = 0; i < num_atoms; i++) */

          true_atoms = count;

          return;
     }                     /* end if(resid[0] != NULLCHAR) */

     j = 0;
     k = 0;

     while(j < num_atoms) {
          k++;
          for(i = 0; i < num_atoms; i++) {
               if(back_order[i] > -1)
                    continue;                 /* already set? */
               if(k == atom_segno[i]) {
                    if(atom_mass[i] == 0.0) {
                         atom_order[j] = -1;
                         back_order[i] = -1;
                         j++;
                         continue;
                    }
                    atom_order[j] = i;
                    back_order[i] = count;
                    j++;
                    count++;
/* set attached H's next to parent */
                    if((bond_count[i] > 0) && (atom_atno[i] > 1))
                         for(l = 0; l < bond_count[i]; l++) {
                              if(atom_mass[(connect[i][l])] == 0.0)
                                   continue;
                              if(atom_atno[(connect[i][l])] == 1) {
                                   atom_order[j] = connect[i][l];
                                   back_order[(connect[i][l])] = count;
                                   j++;
                                   count++;
                              }  /* end if(atom_atno[(connect[i][l])] == 1) */
                         }       /* end for(l = 0; l < bond_count[i]; l++) */
                                 /* and end if((bond_count[i] > 0) &&
                                    (atom_atno[i] > 1)) */
               }                 /* end if(k == atom_segno[i]) */
          }                      /* end for(i = 0; i < num_atoms; i++) */
     }                           /* while(j < num_atoms) */

     true_atoms = count;

     return;
}

/* read a gromacs intermediate atom types force field
   parameters are:
        ffname            name of force field file to read
        nmmol             number of molecule atom types read
        nmbond            number of bond parameter lines read
        nmangle           number of angle parameter lines read
        nmtors            number of torsion parameter lines read
        nmimpro           number of improper parameter lines read
        nmvdw             number of van der Waals parameter entries
*/
void read_ffgmx1(char *ffname, int *nmmol, int *nmbond, int *nmangle,
                 int *nmtors, int *nmimpro, int *nmvdw)
{
     FILE *ffptr;
     char myfiles[MAXCHAR], my_type[8], *c, *d;
     int i, j, k, countin, vdwtype;
     double dum1, dum2;
     int set_same_tors = 0;

     strcpy(myfiles, ffname);

     if((ffptr = fopen(myfiles, "r")) == NULL) {
          sprintf(mess,
                  "Cannot open atomic parameters file %s\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

/* get macro names parameters */
     countin = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "equivs", 6)) {
               sscanf(&line[7], "%d", &countin);
               break;
          }
     }

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks macro names count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     macro_name = (char **)mat_alloc(countin, 4, FARGS);
     macro_nmbr = (int *)allocator(countin, sizeof(int), FARGS);
     macro_types = (char ***)allocator(countin, sizeof(char**), FARGS);
     num_macro = countin;
     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          strncpy(macro_name[i], line, 3);  /* all macro names are */
          macro_name[i][3] = NULLCHAR;      /* 3 characters */
          sscanf(&line[5], "%d", &macro_nmbr[i]);
          macro_types[i] = (char **)mat_alloc(macro_nmbr[i], 6, FARGS);

/* list of atom types without comments or blanks follows */
          for(j = 0; j < macro_nmbr[i]; j++) {
               if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
                    sprintf(mess,
                            "Premature end to force field file %s\n",
                            myfiles);
                    my_fatal(FARGS, mess);
                    exit(1);
               }
               c = get_name(line, macro_types[i][j], 0);
          }          /* end for(j = 0; j < macro_nmbr[i]; j++) */
          i++;
     }               /* end while(i < countin) */

/* get atoms parameters including van der Waals */
     countin = 0;
     vdwtype = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "atoms", 5)) {
               sscanf(&line[6], "%d %d", &countin, &vdwtype);
               break;
          }
     }

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks atom types count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     ffatom_name = (char **)mat_alloc(countin, 6, FARGS);
     ffatom_mass = (float *)allocator(countin, sizeof(float), FARGS);
     ffatom_pol = (float *)allocator(countin, sizeof(float), FARGS);
     vdw_atom_name = (char **)mat_alloc(countin, 6, FARGS);
     vdw_atom_radius = (float *)allocator(countin, sizeof(float), FARGS);
     vdw_atom_pot = (float *)allocator(countin, sizeof(float), FARGS);
     *nmmol = countin;
     *nmvdw = countin;

     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          c = get_name(line, ffatom_name[i], 0);
          strcpy(vdw_atom_name[i], ffatom_name[i]);
          vdw_atom_radius[i] = 99999.0;
          vdw_atom_pot[i] = 99999.0;
          ffatom_mass[i] = 99999.0;
          ffatom_pol[i] = 99999.0;
          sscanf(c, "%f %f %lf %lf", &ffatom_mass[i], &ffatom_pol[i],
                 &dum1, &dum2);
          if((vdwtype == 6) && (dum2 > 0.0)) {
               vdw_atom_radius[i] = CubeRoot((fast_InvSqrt((dum1/dum2))));
               vdw_atom_pot[i] = (dum1 * dum1)/4/dum2;
          }
          else {
               vdw_atom_radius[i] = dum1;
               vdw_atom_pot[i] = dum2;
          }
          i++;
     }

/* get bonds parameters */
    countin = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "bonds", 5)) {
               sscanf(&line[6], "%d", &countin);
               break;
          }
     }

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks bonds count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     bond_name1 = (char **)mat_alloc(countin, 6, FARGS);
     bond_name2 = (char **)mat_alloc(countin, 6, FARGS);
     bond_force = (float *)allocator(countin, sizeof(float), FARGS);
     bond_beta = (float *)allocator(countin, sizeof(float), FARGS);
     bond_length = (float *)allocator(countin, sizeof(float), FARGS);
     FF_bond_type = (short *)allocator(countin, sizeof(short), FARGS);
     *nmbond = countin;

     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          bond_force[i] = -1.0;
          bond_length[i] = 99999.0;
          FF_bond_type[i] = 0;
          c = get_name(line, bond_name1[i], 0);
          c = get_name(c, bond_name2[i], 0);
          c = get_name(c, my_type, 0);
          if((d = strchr(my_type, ' ')) != NULL)
               *d = NULLCHAR;
          FF_bond_type[i] = make_int(my_type,
                                     "Illegal characters in bond type %s\n",
                                     FARGS);
          switch(FF_bond_type[i]) {
               case 5:  break;

               case 1:
               case 2:
               case 6:
               case 8:
               case 9:
               case 7:  sscanf(c, "%f %f", &bond_length[i],
                               &bond_force[i]);
                        break;

               case 3:
               case 4:  sscanf(c, "%f %f %f", &bond_length[i],
                               &bond_force[i], &bond_beta[i]);
                        break;

               default: sprintf(mess,
                          "Illegal bond type in parameters, force field file %s\n",
                          myfiles);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
          i++;
     }

/* get angles parameters or constraints parameters */
     countin = 0;
     constrain = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "angles", 6)) {
               sscanf(&line[7], "%d", &countin);
               break;
          }          /* end if(!strncmp(line, "angles", 6)) */
          if(!strncmp(line, "constraints", 11)) {
               sscanf(&line[12], "%d", &constrain);
               break;
          }          /* end if(!strncmp(line, "constraints", 11)) */
     }               /* end while(get_a_line(line, MAXCHAR, ffptr) != NULL) */

     if(constrain) {
          constrain_name1 = (char **)mat_alloc(constrain, 6, FARGS);
          constrain_name2 = (char **)mat_alloc(constrain, 6, FARGS);
          constrain_dist = (float *)allocator(constrain, sizeof(float), FARGS);
          type_constr = (short *)allocator(constrain, sizeof(short), FARGS);
          i = 0;
/* use of while instead of for to permit comments and blank lines */
          while(i < constrain) {
               if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
                    sprintf(mess,
                            "Premature end to force field file %s\n",
                            myfiles);
                    my_fatal(FARGS, mess);
                    exit(1);
               }     /* end if(get_a_line(line, MAXCHAR, ffptr) == NULL) */
               constrain_dist[i] = 99999.0;
               type_constr[i] = 0;
               c = get_name(line, constrain_name1[i], 0);
               c = get_name(c, constrain_name2[i], 0);
               c = get_name(c, my_type, 0);
               if((d = strchr(my_type, ' ')) != NULL)
                    *d = NULLCHAR;
               type_constr[i] = make_int(my_type,
                                     "Illegal characters in constraint type %s\n",
                                         FARGS);
               sscanf(c, "%f", &constrain_dist[i]);
               i++;
          }          /* end while(i < constrain) */

/* now need to get angles since there were constraint types */
          while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
               if(!strncmp(line, "angles", 6)) {
                    sscanf(&line[7], "%d", &countin);
                    break;
               }     /* end if(!strncmp(line, "angles", 6)) */
          }          /* end while(get_a_line(line, MAXCHAR, ffptr) != NULL) */
     }               /* end if(constrain) */

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks angles count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     angle_name1 = (char **)mat_alloc(countin, 6, FARGS);
     angle_name2 = (char **)mat_alloc(countin, 6, FARGS);
     angle_name3 = (char **)mat_alloc(countin, 6, FARGS);
     angle_quartic = (float **)allocator(5, sizeof(float *), FARGS);
     for(i = 0; i < 5; i++)
          angle_quartic[i] = (float *)allocator(countin,
                                                sizeof(float),
                                                FARGS);
     angle_force = angle_quartic[0];
     angle_angle = (float *)allocator(countin, sizeof(float),
                                      FARGS);
     FF_angle_type = (short *)allocator(countin, sizeof(short),
                                        FARGS);
     *nmangle = countin;

     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          angle_force[i] = -1.0;
          angle_angle[i] = 99999.0;
          FF_angle_type[i] = 0;
          c = get_name(line, angle_name1[i], 0);
          c = get_name(c, angle_name2[i], 0);
          c = get_name(c, angle_name3[i], 0);
          c = get_name(c, my_type, 0);
          if((d = strchr(my_type, ' ')) != NULL)
               *d = NULLCHAR;
          FF_angle_type[i] = make_int(my_type,
                                 "Illegal characters in angle type %s\n",
                                      FARGS);
          switch(FF_angle_type[i]) {
               case 1:
               case 2:  sscanf(c, "%f %f", &angle_angle[i],
                               &angle_force[i]);
                        break;

               case 3:  sscanf(c, "%f %f %f", &angle_quartic[0][i],
                               &angle_quartic[1][i],
                               &angle_quartic[2][i]);
                        break;

               case 4:  sscanf(c, "%f %f %f %f", &angle_quartic[0][i],
                               &angle_quartic[1][i], &angle_quartic[2][i],
                               &angle_quartic[3][i]);
                        break;

               case 5:  sscanf(c, "%f %f %f %f", &angle_angle[i],
                               &angle_quartic[0][i], &angle_quartic[1][i],
                               &angle_quartic[2][i]);
                        break;

               case 6:  sscanf(c, "%f %f %f %f %f %f", &angle_angle[i],
                               &angle_quartic[0][i], &angle_quartic[1][i],
                               &angle_quartic[2][i], &angle_quartic[3][i],
                               &angle_quartic[4][i]);
                        break;

               case 8:  sscanf(c, "%f %f", &angle_quartic[0][i],
                               &angle_quartic[1][i]);
                        break;

               default: sprintf(mess,
                          "Illegal angle type in parameters, force field file %s\n",
                                 myfiles);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
          i++;
     }

/* get regular dihedrals parameters */
     countin = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "torsions", 8)) {
               sscanf(&line[9], "%d", &countin);
               break;
          }
     }

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks torsions count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     tors_name1 = (char **)mat_alloc(countin, 6, FARGS);
     tors_name2 = (char **)mat_alloc(countin, 6, FARGS);
     tors_name3 = (char **)mat_alloc(countin, 6, FARGS);
     tors_name4 = (char **)mat_alloc(countin, 6, FARGS);
     tors_RB = (float **)allocator(6, sizeof(float *), FARGS);
     for(i = 0; i < 6; i++)
          tors_RB[i] = (float *)allocator(countin, sizeof(float),
                                          FARGS);
     tors_mult = (int *)tors_RB[0];
     tors_force = tors_RB[1];
     tors_phase = tors_RB[2];
     tors_term = tors_RB[3];
     tors_more = (int *)tors_RB[4];
     FF_tors_type = (short *)allocator(countin, sizeof(short),
                                       FARGS);
     gmx_tors_resid = (char ***)allocator(countin, sizeof(char **),
                                          FARGS);
     gmx_tors_num = (int *)allocator(countin, sizeof(int), FARGS);
     *nmtors = countin;
     set_same_tors = 0;

     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          tors_more[i] = 0;
          gmx_tors_resid[i] = NULL;
          gmx_tors_num[i] = 0;
          tors_mult[i] = 9999;
          tors_force[i] = -1.0;
          tors_phase[i] = 99999.0;
          tors_term[i] = 1.0;
          c = get_name(line, tors_name1[i], 0);
          c = get_name(c, tors_name2[i], 0);
          c = get_name(c, tors_name3[i], 0);
          c = get_name(c, tors_name4[i], 0);
          c = get_name(c, my_type, 0);
          if((d = strchr(my_type, ' ')) != NULL)
               *d = NULLCHAR;
          FF_tors_type[i] = make_int(my_type,
                                 "Illegal characters in torsion type %s\n",
                                     FARGS);
          switch(FF_tors_type[i]) {
               case 9:
               case 1:  sscanf(c, "%f %f %d %d", &tors_phase[i],
                               &tors_force[i], &tors_mult[i],
                               &gmx_tors_num[i]); 
                        tors_term[i] = (float)abs(tors_mult[i]);
                        if(gmx_tors_num[i])
                             special_tors(gmx_tors_num[i], ffptr, i, myfiles);

                        while(set_same_tors) {
                             if(gmx_tors_num[i] != gmx_tors_num[(i - 1)]) {
                                  set_same_tors = 0;
                                  break;
                             }
                             if(!gmx_tors_num[i])
                                  break;
                             for(j = 0; j < gmx_tors_num[i]; j++) 
                                  if(strcmp(gmx_tors_resid[i][j],
                                            gmx_tors_resid[(i - 1)][j])) {
                                       set_same_tors = 0;
                                       break;
                                  }
                             break;      /* needs a break when using
                                               a while as an if */
                        }                /* end while(set_same_tors) */
                        if(set_same_tors) {
                             if(!strcmp(tors_name1[i],
                                        tors_name1[(i - 1)]) &&
                                !strcmp(tors_name2[i],
                                        tors_name2[(i - 1)]) &&
                                !strcmp(tors_name3[i],
                                        tors_name3[(i - 1)]) &&
                                !strcmp(tors_name4[i],
                                        tors_name4[(i - 1)])) 
                                     flag_more_tors(i);
                             else set_same_tors = 0;
                        }

                        if(tors_mult[i] < 0) set_same_tors = 1;
                        else set_same_tors = 0;
                        break;

               case 3:  sscanf(c, "%f %f %f %f %f %f %d", &tors_RB[0][i],
                               &tors_RB[1][i], &tors_RB[2][i],
                               &tors_RB[3][i], &tors_RB[4][i],
                               &tors_RB[5][i], &gmx_tors_num[i]);
                        set_same_tors = 0;
                        if(gmx_tors_num[i])
                             special_tors(gmx_tors_num[i], ffptr, i, myfiles);
                        break;

               case 5:  sscanf(c, "%f %f %f %f %d", &tors_RB[0][i],
                               &tors_RB[1][i], &tors_RB[2][i],
                               &tors_RB[3][i], &gmx_tors_num[i]);
                        set_same_tors = 0;
                        if(gmx_tors_num[i])
                             special_tors(gmx_tors_num[i], ffptr, i, myfiles);
                        break;

               case 8:  sscanf(c, "%f %f %d", &tors_RB[0][i], &tors_RB[1][i],
                               &gmx_tors_num[i]);
                        set_same_tors = 0;
                        if(gmx_tors_num[i])
                             special_tors(gmx_tors_num[i], ffptr, i, myfiles);
                        break;

               default: sprintf(mess,
                          "Illegal torsion type in parameters, force field file %s\n",
                                myfiles);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
          i++;
     }

/* get improper dihedrals parameters */
     countin = 0;
     while(get_a_line(line, MAXCHAR, ffptr) != NULL) {
          if(!strncmp(line, "impropers", 9)) {
               sscanf(&line[10], "%d", &countin);
               break;
          }
     }

     if(!countin) {
          sprintf(mess,
                  "Force field file %s lacks impropers count data\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     impro_name1 = (char **)mat_alloc(countin, 6, FARGS);
     impro_name2 = (char **)mat_alloc(countin, 6, FARGS);
     impro_name3 = (char **)mat_alloc(countin, 6, FARGS);
     impro_name4 = (char **)mat_alloc(countin, 6, FARGS);
     impro_mult = (int *)allocator(countin, sizeof(float), FARGS);
     impro_force = (float *)allocator(countin, sizeof(float), FARGS);
     impro_phase = (float *)allocator(countin, sizeof(float), FARGS);
     impro_term = (float *)allocator(countin, sizeof(float), FARGS);
     impro_xcount = (int *)allocator(countin, sizeof(int), FARGS);
     *nmimpro = countin;

     i = 0;
/* use of while instead of for to permit comments and blank lines */
     while(i < countin) {
          if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
               sprintf(mess,
                       "Premature end to force field file %s\n",
                       myfiles);
               my_fatal(FARGS, mess);
               exit(1);
          }
          impro_force[i] = -1.0;
          impro_phase[i] = 99999.0;
          impro_term[i] = 1.0;
          c = get_name(line, impro_name3[i], 0);
          c = get_name(c, impro_name1[i], 0);
          c = get_name(c, impro_name2[i], 0);
          c = get_name(c, impro_name4[i], 0);
          impro_xcount[i] = 0;
          if(impro_name1[i][0] == 'X') impro_xcount[i]++;
          if(impro_name2[i][0] == 'X') impro_xcount[i]++;
          if(impro_name3[i][0] == 'X') impro_xcount[i]++;
          if(impro_name4[i][0] == 'X') impro_xcount[i]++;
          sscanf(c, "%f %f", &impro_phase[i], &impro_force[i]);
          impro_mult[i] = 1;
          if(impro_force[i] < 0.0 ) impro_mult[i] = 9999;
          i++;
     }

     fclose(ffptr);

     return;
}

/* Set the names array for a torsion that is only applicable to certain
   residues
   parameters are:
        numb               number of residue names to read from the
                           next line of the file.
        ffptr              file pointer for force field file to read
        i                  number of the current force field record
        myfiles            name of force field file being read
*/
void special_tors(int numb, FILE *ffptr, int i, char *myfiles)
{
     int j;
     char *c;

     if(numb > 10) {
          sprintf(mess,
                  "Bad torsions residue count in force field file %s\n",
                  myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }

     if(get_a_line(line, MAXCHAR, ffptr) == NULL) {
          sprintf(mess,
                  "Premature end to force field file %s\n",
                   myfiles);
          my_fatal(FARGS, mess);
          exit(1);
     }
     gmx_tors_resid[i] = (char **)mat_alloc(numb, 6, FARGS);
     c = &line[0];
     for(j = 0; j < numb; j++)
           c = get_name(c, &gmx_tors_resid[i][j][0], 0);

     return;
}

/* Set torsion parameters from gromacs force field read in.
   Includes checks for residue specific torsion entries.
   parameters are:
        number of torsions in molecule:              num_tors
        number of torsion types in force field:      tors_cntFF
*/
void set_gromacs_torsFF(int num_tors, int tors_cntFF)
{
     int i, tors1, tors2, tors3, tors4, j, k;

     for(i = 0; i < num_tors; i++) {
          tors1 = torsion_tab[i][0];
          tors2 = torsion_tab[i][1];
          tors3 = torsion_tab[i][2];
          tors4 = torsion_tab[i][3];

          j = comp_gmx_tors(tors1, tors2, tors3, tors4,
                            tors_cntFF);

          if(j) a_tors_type[i] = FF_tors_type[(j - 1)];
          else a_tors_type[i] = 0;

          switch(a_tors_type[i]) {
               case 0: break;          /* torsion type not set. */

               case 9:
               case 1: a_tors_force[i] = tors_force[(j - 1)];
                       a_tors_mult[i] = tors_mult[(j - 1)];
                       a_tors_phase[i] = tors_phase[(j - 1)];
                       a_tors_term[i] = tors_term[(j - 1)];
                       if(j >= (tors_cntFF - 1))
                             break;
/* if this has multiple torsions, add links and set them as well
   note:  j points to the next torsion entry which is the next
   component of a multiple torsion and these facts are used below.
*/
                       if(tors_more[(j - 1)]) {
                            multors_cnt[i] = tors_more[(j - 1)];     /* set count */
                            multiples[i] = (float **)allocator(multors_cnt[i],
                                                               sizeof(float **),
                                                               FARGS);
                            for(k = 0; k < multors_cnt[i]; k++) {
                                 multiples[i][k] = (float *)allocator(4,
                                                                      sizeof(float),
                                                                      FARGS);
                                 multiples[i][k][0] = (float)tors_mult[j];
                                 multiples[i][k][1] = tors_force[j];
                                 multiples[i][k][2] = tors_phase[j];
                                 multiples[i][k][3] = tors_term[j];
                                 j++;
                            }
                       }
                       break;

               case 3: for(k = 0; k < 6; k++)
                            a_tors_set[k][i] = tors_RB[k][(j - 1)];
                       break;

               case 5: for(k = 0; k < 4; k++)
                            a_tors_set[k][i] = tors_RB[k][(j - 1)];
                       break;

               case 8: a_tors_set[0][i] = tors_RB[0][(j - 1)];
                       a_tors_set[1][i] = tors_RB[1][(j - 1)];
                       break;

               default: my_fatal(FARGS, "Illegal torsion type returned\n");
                        exit(1);
          }

     }

    return;
}

/* comparison of torsion to force field
   parameters are:
        atom types of each atom in torsion:               tors_atom1
                                                          tors_atom2
                                                          tors_atom3
                                                          tors_atom4
        number of torsion parameter sets in force field:  tors_cntFF
*/
int comp_gmx_tors(int tors_atom1, int tors_atom2, int tors_atom3,
                  int tors_atom4, int tors_cntFF)
{
     int j, k, l, no;

/* first try without wild cards */
     for(j = 0; j < tors_cntFF; j++) {
          no = 1;                            /* flag for special residue requirements */
/* check for dihedrals with special residue requirements.
   Uses the molecule's original residue names and not
   changed residue name.
*/
          while(gmx_tors_num[j]) {        /* use while instead of if for escapes below */
               if(strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom2]) ||
                  strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom3]) ||
                  strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom4]))
                    break;           /* all atoms must be in the same original residue */
               l = 0;
               for(k = 0; k < gmx_tors_num[j]; k++)
                    if(!strcmp(atom_residue[tors_atom1],
                               gmx_tors_resid[j][k])) {
                         l = 1;
                         break;
                    }
/* if failed to find a needed match.  skip this force field entry */
               if(!l) no = 0;
               break;             /* really the while was an if, so need a break here */
          }
          if(!no) continue;              /* special residue requirements not met */

          if((!atom_cmp(a_ff_type[tors_atom2], tors_name2[j]) &&
              !atom_cmp(a_ff_type[tors_atom3], tors_name3[j])))    /* two center must match */
                    if((!atom_cmp(a_ff_type[tors_atom1],
                                   tors_name1[j])) &&        /* 1 goes to 2 and 4 goes to 3 */
                       (!atom_cmp(a_ff_type[tors_atom4],
                                  tors_name4[j])))          /* order here is important */
                              return(j + 1);

          if((!atom_cmp(a_ff_type[tors_atom2], tors_name3[j]) &&
              !atom_cmp(a_ff_type[tors_atom3], tors_name2[j])))   /* two center must match */
                    if((!atom_cmp(a_ff_type[tors_atom1],
                                  tors_name4[j])) &&       /* 4 goes to 2 and 1 goes to 3 */
                       (!atom_cmp(a_ff_type[tors_atom4],
                                  tors_name1[j])))             /* order here is important */
                              return(j + 1);
     }

/* then try with wild cards */
     for(j = 0; j < tors_cntFF; j++) {
          no = 1;                            /* flag for special residue requirements */
/* check for dihedrals with special residue requirements.
   Uses the molecule's original residue names and not
   changed residue name.
*/
          while(gmx_tors_num[j]) {        /* use while instead of if for escapes below */
               if(strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom2]) ||
                  strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom3]) ||
                  strcmp(atom_residue[tors_atom1],
                         atom_residue[tors_atom4]))
                    break;           /* all atoms must be in the same original residue */
               l = 0;
               for(k = 0; k < gmx_tors_num[j]; k++)
                    if(!strcmp(atom_residue[tors_atom1],
                               gmx_tors_resid[j][k])) {
                         l = 1;
                         break;
                    }
/* if failed to find a needed match.  skip this force field entry */
               if(!l) no = 0;
               break;             /* really the while was an if, so need a break here */
          }
          if(!no) continue;              /* special residue requirements not met */

          if((!atom_cmp(a_ff_type[tors_atom2], tors_name2[j]) &&
              !atom_cmp(a_ff_type[tors_atom3], tors_name3[j])))    /* two center must match */
                    if(((!atom_cmp(a_ff_type[tors_atom1], tors_name1[j])) ||
                        (tors_name1[j][0] == 'X')) &&      /* 1 goes to 2 and 4 goes to 3 */
                       ((!atom_cmp(a_ff_type[tors_atom4], tors_name4[j])) ||
                        (tors_name4[j][0] == 'X')))        /* order here is important */
                              return(j + 1);

          if((!atom_cmp(a_ff_type[tors_atom2], tors_name3[j]) &&
              !atom_cmp(a_ff_type[tors_atom3], tors_name2[j])))   /* two center must match */
                    if(((!atom_cmp(a_ff_type[tors_atom1], tors_name4[j])) ||
                        (tors_name4[j][0] == 'X')) &&     /* 4 goes to 2 and 1 goes to 3 */
                       ((!atom_cmp(a_ff_type[tors_atom4], tors_name1[j])) ||
                        (tors_name1[j][0] == 'X')))       /* order here is important */
                              return(j + 1);
     }

     return(0);                                            /* nothing matched */
}

/* set constraints data
   parameter:
        number of angles in molecule          num_cons
*/
void set_constraints(int num_cons)
{
     int i, j;

     constrain_i = (int *)allocator(num_cons, sizeof(int), FARGS);
     constrain_j = (int *)allocator(num_cons, sizeof(int), FARGS);
     dist_constr = (float *)allocator(num_cons, sizeof(float), FARGS);
     constr_type = (short *)allocator(num_cons, sizeof(short), FARGS);
     no_constraints = 0;

     for(i = 0; i < num_cons; i++) {
          constrain_i[i] = -1;
          constrain_j[i] = -1;
          dist_constr[i] = 99999.0;
          constr_type[i] = 0;

          for(j = 0; j < constrain; j++) {
               if(((!atom_cmp(a_ff_type[(angle_tab[i][0])],
                              constrain_name1[j])   &&
                    !atom_cmp(a_ff_type[(angle_tab[i][2])],
                              constrain_name2[j]))) ||
                  ((!atom_cmp(a_ff_type[(angle_tab[i][2])],
                              constrain_name1[j])   &&
                    !atom_cmp(a_ff_type[(angle_tab[i][0])],
                              constrain_name2[j])))) {
                                   constrain_i[i] = angle_tab[i][0];
                                   constrain_j[i] = angle_tab[i][2];
                                   dist_constr[i] = constrain_dist[j];
                                   constr_type[i] = type_constr[j];
                                   no_constraints++;
                                   break;
               }
          }               /* end for(j = 0; j < constraints; j++) */
     }                    /* end for(i = 0; i < num_cons; i++) */

     if(!no_constraints) {
          free_me(constrain_i, FARGS);
          free_me(constrain_j, FARGS);
          free_me(dist_constr, FARGS);
          free_me(constr_type, FARGS);
     }

     return;
}

/* set atom charge from force field read in.
   parameters are:
        number of atoms in molecule                  num_atoms
        number of atom types in force field          num_mol
*/
void set_charge(int num_atoms, int num_mol)
{
     int i, j;

     for(i = 0; i < num_atoms; i++) {
          for(j = 0; j < num_mol; j++) {
               if(strcmp(a_ff_type[i],
                         ffatom_name[j])) continue;      /* this atom type? */
               atom_charge[i] = ffatom_pol[j];
               break;
          }                             /* end for(j = 0; j < num_mol; j++) */
     }                                  /* end for(i = 0; i < num_atoms; i++) */

     return;
}
