/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to write a simple gromacs format .gro coordinates file
        and to write a gromacs format .top topolgy file

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
#include "gromacs_FF.h"
#include "block_memory.h"
#include "mol2.h"
#include "multors.h"
#include "GMX_FF.h"

#define MAXCHAR 256
#define NULLCHAR 0

extern int lstart_at;
extern int res_start_at;
extern int rstart_at;
extern int lstart2_at;
extern char *use_define[6];
extern int use_amber;
extern int center_it;

int *cgnr_grp = NULL;

extern char mess[MAXCHAR];
extern char *in_command;
extern char vers_no[];

/* this is a very simplistic write of a simple gromacs .gro format coordinates file
   that computes the box.
   parameters are:
        filename                 name of .gro file
        atomnum                  number of atoms
        names                    list of atom names
        residue                  the optional residue name for all the atoms
*/
void write_gro (char *filename, int atomnum, char **names, char *residue)
{
     FILE *fpout;
     int i, j, k, count;
     float xmin, ymin, zmin, xmax, ymax, zmax;
     float box_x, box_y, box_z;
     double xavg, yavg, zavg, sig_mass;
     char dummy[6], *c, *d;

     xavg = 0.0;
     yavg = 0.0;
     zavg = 0.0;

     if(center_it) {
          sig_mass = 0.0;
          for(i = 0; i < atomnum; i++) {
               if(atom_order[i] < 0) continue;
               xavg += atom_x[(atom_order[i])] *
                       (double)atom_mass[(atom_order[i])];
               yavg += atom_y[(atom_order[i])] *
                       (double)atom_mass[(atom_order[i])];
               zavg += atom_z[(atom_order[i])] *
                       (double)atom_mass[(atom_order[i])];
               sig_mass += atom_mass[(atom_order[i])];
          }

          xavg = xavg/sig_mass;
          yavg = yavg/sig_mass;
          zavg = zavg/sig_mass;
     }

     if((fpout = fopen(filename, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     j = 1;
     dummy[0] = NULLCHAR;
     fprintf(fpout, "%s\n", mol_name);
     fprintf(fpout, "%5d\n", true_atoms);

     xmin = atom_x[0] * 0.1;
     xmax = atom_x[0] * 0.1;
     ymin = atom_y[0] * 0.1;
     ymax = atom_y[0] * 0.1;
     zmin = atom_z[0] * 0.1;
     zmax = atom_z[0] * 0.1;

     if(residue[0] != NULLCHAR) {
          strncpy(dummy, residue, 5);
          dummy[5] = NULLCHAR;
     }

     count = 0;
     for(i = 0; i < atomnum; i++) {
          if(atom_order[i] < 0) continue;
          if(residue[0] == NULLCHAR) {
               j = atom_segno[(atom_order[i])];
               strncpy(dummy, atom_residue[(atom_order[i])], 5);
               dummy[5] = NULLCHAR;
               if((c = strpbrk(dummy, "0123456789")) != NULL) {
                    k = strlen(dummy);
                    if((c != &dummy[0]) &&
                       ((d = strpbrk(&dummy[k-1], "0123456789")) != NULL))
                              *c = NULLCHAR;
               }
          }

          fprintf(fpout, "%5d%-5.5s%5.5s%5d%8.3f%8.3f%8.3f\n",
                  (j + res_start_at), dummy, names[(atom_order[i])],
                  ((count + lstart_at)%100000),
                  ((atom_x[(atom_order[i])] - (float)xavg) * 0.1),
                  ((atom_y[(atom_order[i])] - (float)yavg) * 0.1),
                  ((atom_z[(atom_order[i])] - (float)zavg) * 0.1));
          count++;

          if((atom_x[(atom_order[i])] * 0.1) > xmax)
               xmax = atom_x[(atom_order[i])] * 0.1;
          if((atom_y[(atom_order[i])] * 0.1) > ymax)
               ymax = atom_y[(atom_order[i])] * 0.1;
          if((atom_z[(atom_order[i])] * 0.1) > zmax)
               zmax = atom_z[(atom_order[i])] * 0.1;
          if((atom_x[(atom_order[i])] * 0.1) < xmin)
               xmin = atom_x[(atom_order[i])] * 0.1;
          if((atom_y[(atom_order[i])] * 0.1) < ymin)
               ymin = atom_y[(atom_order[i])] * 0.1;
          if((atom_z[(atom_order[i])] * 0.1) < zmin)
               zmin = atom_z[(atom_order[i])] * 0.1;
     }

     box_x = 3.0*(xmax - xmin);
     box_y = 3.0*(ymax - ymin);
     box_z = 3.0*(zmax - zmin);

     fprintf(fpout, "%10.5f %10.5f %10.5f\n", box_x, box_y, box_z);

     fclose(fpout);

     return;
}

/* write the itp file that includes all else
   parameters are:
        filename                     the itp file name for all output
        residue                      the residue name for all the atoms
        numatoms                     number of atoms
        numbonds                     number of bonds
        numpairs                     number of pairs
        numangles                    number of angles
        numdihed                     number of dihedrals
        numimpro                     number of impropers
        assigned_types               list of atom types
        names                        list of atom names
        usage                        flag to use measured or force field
                                     derived lengths and angles
        amber                        flag for data in kcal/mol and angstroms
*/
void write_top(char *filename, char *residue, int numatoms, int numbonds,
               int numpairs, int numangles, int numdihed, int numimpro,
               char **assigned_types, char **names, int usage, int amber)
{
     int tot_rstr;
     char *c;
     double conv = 4.18400000;
     FILE *fptop;

     tot_rstr = rstr_dist + rstr_range + rstr_angle + rstr_torsn;

     if((fptop = fopen(filename, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     fprintf(fptop, ";\n; Topology from .mol2 file\n");
     fprintf(fptop, "; topolbuild %s\n", vers_no);
     fprintf(fptop, "; Command line:\n;     %s\n;\n",
             in_command);
     if(!res_start_at) {
          fprintf(fptop, "; The force field files to be included\n");
          strcpy(mess, "\"ff");
          strcat(mess, filename);
          c = strrchr(mess, '.');
          *c = NULLCHAR;
          strcat(mess, ".itp\"");
          fprintf(fptop, "#include %s\n\n", mess);
     }

     write_atoms(fptop, residue, numatoms, assigned_types, names);
     write_bonds(fptop, numbonds, names, amber, conv, usage,
                 tot_rstr);
     if(no_constraints > 0)
          write_constr(fptop, numangles, names, amber, usage);
     write_pairs(fptop, numpairs, names);
     write_angles(fptop, numangles, names, usage, amber, conv);
     if((numdihed > 0) && (numdihed != tors_imp))
          write_dihed(fptop, numdihed, names, usage, amber, conv);
     if((numimpro > 0) || (tors_imp > 0))
          write_impropers(fptop, numimpro, numdihed, names, usage,
                          amber, conv);
     write_restraints(fptop, tot_rstr, names);

     if(!res_start_at) {
/* write position restraints inclusion lines unless writing an addition */
          strcpy(mess, "\"posre");
          strcat(mess, filename);
          c = strrchr(mess, '.');
          *c = NULLCHAR;
          strcat(mess, ".itp\"");
          fprintf(fptop, "\n; Include Position restraint file\n");
          fprintf(fptop,
  "; WARNING: Position restraints and distance restraints ought not be done together\n");
          fprintf(fptop, "#ifdef POSRES\n");
          fprintf(fptop, "#include %s\n", mess);
          fprintf(fptop, "#endif\n");

/* write water topology inclusion lines unless writing an addition */
          fprintf(fptop, "\n; Include water topology\n");
          fprintf(fptop, "#include \"%s\"\n", use_define[2]);

          fprintf(fptop, "\n#ifdef POSRES_WATER\n");
          fprintf(fptop, "; Position restraint for each water oxygen\n");
          fprintf(fptop, "[ position_restraints ]\n");
          fprintf(fptop, ";  i funct       fcx        fcy        fcz\n");
          fprintf(fptop, "   1    1       1000       1000       1000\n#endif\n");

/* write ions inclusion lines unless writing an addition */
          fprintf(fptop, "\n; Include generic topology for ions\n");
          fprintf(fptop, "#include \"%s\"\n", use_define[3]);

/* write system information lines unless writing an addition */
          fprintf(fptop, "\n [ system ]\n; title from mol2 input\n");
          fprintf(fptop, "%s\n\n", mol_name);

/* write molecules information lines unless writing an addition */
          fprintf(fptop, " [ molecules ]\n; molecule name    nr.\n");
          fprintf(fptop, "%8s           1\n", mol_name);
     }

     fclose(fptop);
     return;
}

/* write the atoms records
   parameters are:
        fptop                        file pointer to write to
        residue                      the residue name for all the atoms
        numatoms                     number of atoms
        assigned_types               list of atom types
        names                        list of atom names
*/
void write_atoms(FILE *fptop, char *residue, int numatoms, char **assigned_types,
                 char **names)
{
     char *c, *d;
     char dummy[6], thistype[12];
     int i, j, k, count, sum_pos, sum_neg;
     double sigma_chg;

     cgnr_assgn(numatoms, residue);

     dummy[0] = NULLCHAR;

     if(residue[0] != NULLCHAR) {
          strncpy(dummy, residue, 5);
          dummy[5] = NULLCHAR;
     }

     fprintf(fptop, " [ moleculetype ]\n");
     fprintf(fptop, "; name  nrexcl\n%8s   3\n\n", mol_name);

/* atoms */
     fprintf(fptop, " [ atoms ]\n");
     fprintf(fptop,
        ";  nr    type   resnr   residu   atom   cgnr    charge      mass\n");

     sum_neg = 0;
     sum_pos = 0;
     j = 1;

     count = 0;
     for(i = 0; i < numatoms; i++) {
          if(atom_order[i] < 0) continue;
          if(residue[0] == NULLCHAR) {
               j = atom_segno[(atom_order[i])];
               strncpy(dummy, atom_residue[(atom_order[i])], 5);
               dummy[5] = NULLCHAR;
               if((c = strpbrk(dummy, "0123456789")) != NULL) {
                    k = strlen(dummy);
                    if((c != &dummy[0]) &&
                       ((d = strpbrk(&dummy[k-1], "0123456789")) != NULL))
                              *c = NULLCHAR;
               }
          }

          if(atom_charge[(atom_order[i])] < 0.0)
               sum_neg += (int)(atom_charge[(atom_order[i])] * 1e5 - 0.5);
          else
               sum_pos += (int)(atom_charge[(atom_order[i])] * 1e5 + 0.5);
          sigma_chg = ((double)(sum_pos + sum_neg))/1e5;
          thistype[0] = NULLCHAR;
          if((use_define[4] != NULL) &&
             (isdigit(assigned_types[(atom_order[i])][0])))
                    strcpy(thistype, use_define[4]);

          if((strlen(thistype) + strlen(assigned_types[(atom_order[i])])) > 11) {
               sprintf(mess, "Name %s%s longer than 11 characters.",
                       thistype, assigned_types[(atom_order[i])]);
               my_fatal(FARGS, mess);
               exit(1);
          }

          strcat(thistype, assigned_types[(atom_order[i])]);
          fprintf(fptop, "%5d %11s %5d  %8s%8s%5d%10.5f%10.5f       ; %11.7f\n",
                  (count + lstart2_at), thistype, (j + rstart_at), dummy,
                  names[(atom_order[i])], cgnr_grp[(atom_order[i])],
                  atom_charge[(atom_order[i])], atom_mass[(atom_order[i])],
                  sigma_chg);
          count++;
     }

     fprintf(fptop, "; total molecule charge = %11.7f\n", sigma_chg);

     return;
}

/* write bond records and simple harmonic bond restraints
   parameters are:
        fptop                        file pointer to write to
        numbonds                     number of bonds
        names                        list of atom names
        amber                        flag for data in kcal/mol and angstroms
        conv                         conversion factor
        usage                        flag to use measured or force field
                                     derived lengths and angles
        tot_rstr                     total number of restraints
*/
void write_bonds(FILE *fptop, int numbonds, char **names, int amber,
                 double conv, int usage, int tot_rstr)
{
     int i, j, idone;
     double bfc, blen;

/* bonds */
     fprintf(fptop, "\n [ bonds ]\n");
     fprintf(fptop, ";   ai  aj   funct      b0          kb\n");

     for(i = 0; i < numbonds; i++) {
          if((back_order[(bond_i[i])] < 0) ||
             (back_order[(bond_j[i])] < 0))
                    continue;
          if(usage) blen = (double)meas_len[i];
          else blen = a_bond_length[i];
          if(amber) blen = blen * 0.1;
          bfc = (double)a_bond_force[i];
          if(amber) bfc = bfc * conv * 200.0;

          fprintf(fptop, "  %6d%6d", (back_order[(bond_i[i])] + lstart2_at),
                  (back_order[(bond_j[i])] + lstart2_at));
          switch(a_bond_type[i]) {
               case 0:  fprintf(fptop,
                           "                                    ;%6s-%6s\n", 
                                names[(bond_i[i])], names[(bond_j[i])]);
                        break;

               case 5:  fprintf(fptop,
                           "   5                                ;%6s-%6s\n", 
                                names[(bond_i[i])], names[(bond_j[i])]);
                        break;


               case 7:
               case 6:
               case 1:
               case 2:  fprintf(fptop, "   %1d%12.5f%12.0f.       ;%6s-%6s\n",
                                a_bond_type[i], blen, bfc, names[(bond_i[i])],
                                names[(bond_j[i])]);
                        break;

               case 3:
               case 4:  fprintf(fptop, "   %1d%12.5%12.0f.%12.0f.  ;%6s-%6s\n",
                                a_bond_type[i], blen, bfc, a_bond_beta[i],
                                names[(bond_i[i])], names[(bond_j[i])]);
                        break;

               case 8:
               case 9:  fprintf(fptop, "   %1d%12d%12.0f.       ;%6s-%6s\n",
                                a_bond_type[i], (int)a_bond_length[i], bfc,
                                names[(bond_i[i])], names[(bond_j[i])]);
                        break;


               default: sprintf(mess,
                           "Bond type error, bond number %d.\n", i);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
     }

/* treat simple distance restraints as bonds with function type 6, harmonic potential */
     if(rstr_dist) {
          idone = 1;
          j = rstr_dist;
          for(i = 0; i < tot_rstr; i++) {
               if(restr[i][0] != 1) continue;
               if((back_order[(restr[i][1] - 1)] < 0) ||
                  (back_order[(restr[i][2] - 1)] < 0)) {
                         j--;
                         if(!j) break;
                         continue;
               }
               if(idone) {
                    fprintf(fptop, "#ifdef DISTREST\n");
                    idone = 0;
               }
               blen = (double)(frest[i][0]) * 0.1;            /* these are always amber type */
               bfc = (double)(frest[i][1]) * conv * 200.0;
               fprintf(fptop, "  %6d%6d%s%12.5f%12.0f%s       ;%6s-%6s\n",
                       (back_order[(restr[i][1] - 1)] + lstart2_at),
                       (back_order[(restr[i][2] - 1)] + lstart2_at),
                       "   6", blen, bfc, names[(restr[i][1] - 1)],
                       names[(restr[i][2] - 1)]);
               j--;
               if(!j) break;
          }
          if(!idone)
               fprintf(fptop, "#endif\n");
     }

     return;
}

/* write constraints records
   parameters are:
        fptop                        file pointer to write to
        numangles                    number of angles
        names                        list of atom names
        amber                        flag for data in kcal/mol and angstroms
        usage                        flag to use measured or force field
                                     derived lengths and angles
*/
void write_constr(FILE *fptop, int numangles, char **names, int amber,
                  int usage)
{
     int i, prnt;
     float conlen;

     prnt = 0;
     for(i = 0; i < numangles; i++) {
          if((constrain_i[i] < 0) ||
             (constrain_j[i] < 0) ||
             (constr_type[i] < 1) ||
             (constr_type[i] > 2))
                    continue;
          if((back_order[(constrain_i[i])] < 0) ||
             (back_order[(constrain_j[i])] < 0))
                    continue;
          if(!prnt) {
/* constraints */
               fprintf(fptop, "\n [ constraints ]\n");
               prnt = 1;
          }
          conlen = dist_constr[i];
          if(amber)
               conlen = conlen * 0.1;
          if(usage)
               conlen = get_dist(constrain_i[i], constrain_j[i]) * 0.1;

          fprintf(fptop, "%8d%8d   %1d%12.5f       ;%6s  %6s\n",
                  (back_order[(constrain_i[i])] + lstart2_at),
                  (back_order[(constrain_j[i])] + lstart2_at),
                  constr_type[i], conlen, names[(constrain_i[i])],
                  names[(constrain_j[i])]);
     }

     return;
}

/* write pairs records
   parameters are:
        fptop                        file pointer to write to
        numpairs                     number of pairs
        names                        list of atom names
*/
void write_pairs(FILE *fptop, int numpairs, char **names)
{
     int i, j, k, not_pair, itry;

/* pairs */
     fprintf(fptop, "\n [ pairs ]\n");

     for(i = 0; i < numpairs; i++) {
          if((back_order[(torsion_tab[i][0])] < 0) ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          if((a_tors_type[(torsion_tab[i][0])] == 3) &&
             (use_amber != 5))
                    continue;    /* Not pair if RB dihedral and FF not oplsaa. */
          not_pair = 0;          /* must not have a 1 or 2 bond connection elsewhere */
          if((bond_count[(torsion_tab[i][0])] > 1) &&
             (bond_count[(torsion_tab[i][3])] > 1)) {
                    for(j = 0; j < (bond_count[(torsion_tab[i][0])]); j++) {
                         itry = connect[(torsion_tab[i][0])][j];
                         if(bond_count[itry] > 1) {
                              if(itry == torsion_tab[i][3]) {  /* direct bond elsewhere? */
                                   not_pair = 1;
                                   break;
                              }       /* end if(itry == torsion_tab[i][3])*/

/* only 1 intervening atom? */
                              for(k = 0; k < bond_count[itry]; k++)
                                   if(connect[itry][k] == torsion_tab[i][3]) {
                                        not_pair = 1;
                                        break;
                                   }  /* end if(connect[itry][k]    etc. */
                              if(not_pair) break;
                         }            /* end if(bond_count[itry] > 1) */
                    }                 /* end for(j = 0; j <    etc. */
          }                           /* end if((bond_count[(torsion_tab[i][0])] > 1) etc. */

          if(i)                       /* eliminate duplicate pairs */
               for(j = 0; j < i; j++)
                    if(((torsion_tab[i][0] == torsion_tab[j][0]) &&
                        (torsion_tab[i][3] == torsion_tab[j][3])) || 
                       ((torsion_tab[i][0] == torsion_tab[j][3]) &&
                        (torsion_tab[i][3] == torsion_tab[j][0]))) {
                              not_pair = 1;
                              break;
                    }

          if(not_pair) continue;
          fprintf(fptop, "%8d%8d  1       ;%6s-%6s\n",
                  (back_order[(torsion_tab[i][0])] + lstart2_at),
                  (back_order[(torsion_tab[i][3])] + lstart2_at),
                  names[(torsion_tab[i][0])], names[(torsion_tab[i][3])]);
     }

     return;
}

/* write the angle parameters records
   parameters are:
        fptop                        file pointer to write to
        numangles                    number of angles
        names                        list of atom names
        usage                        flag to use measured or force field
                                     derived lengths and angles
        amber                        flag for data in kcal/mol and angstroms
        conv                         conversion factor
*/
void write_angles(FILE *fptop, int numangles, char **names, int usage,
                  int amber, double conv)
{
     int i, j;
     double afc;
     float theta, r1, r2, r3;

/* angles */
     fprintf(fptop, "\n[ angles ]\n");
     fprintf(fptop, "; ai  aj  ak  funct      th0         cth\n");

     for(i = 0; i < numangles; i++) {
          if((back_order[(angle_tab[i][0])] < 0) ||
             (back_order[(angle_tab[i][1])] < 0) ||
             (back_order[(angle_tab[i][2])] < 0))
                    continue;
          if(usage) theta = meas_ang[i];
          else theta = a_angle_angle[i];
          afc = (double)a_angle_force[i];
          if(amber) afc = afc * conv * 2.0;

          if(a_angle_type[i] > -1)
               fprintf(fptop, "%6d%6d%6d",
                       (back_order[(angle_tab[i][0])] + lstart2_at),
                       (back_order[(angle_tab[i][1])] + lstart2_at),
                       (back_order[(angle_tab[i][2])] + lstart2_at));
          switch(a_angle_type[i]) {
               case -1: break;

               case 0:  fprintf(fptop,
                           "                                 ;%6s-%6s-%6s\n",
                                names[(angle_tab[i][0])], names[(angle_tab[i][1])],
                                names[(angle_tab[i][2])]);
                        break;

               case 1:
               case 2:  fprintf(fptop, "   %1d%12.3f%12.4f     ;%6s-%6s-%6s\n",
                                a_angle_type[i], theta, afc,
                                names[(angle_tab[i][0])], names[(angle_tab[i][1])],
                                names[(angle_tab[i][2])]);
                        break;

               case 3:  if(usage) {
                             r1 = get_dist(angle_tab[i][0], angle_tab[i][1]) * 0.1;
                             r2 = get_dist(angle_tab[i][1], angle_tab[i][2]) * 0.1;
                        }
                        else {
                             r1 = a_angle_quartic[0][i];
                             r2 = a_angle_quartic[1][i];
                        }
                        fprintf(fptop, "   3%12.3f%12.3f%12.4f     ;%6s-%6s-%6s\n",
                                r1, r2, a_angle_quartic[2][i],
                                names[(angle_tab[i][0])], names[(angle_tab[i][1])],
                                names[(angle_tab[i][2])]);
                        break;

               case 4:  if(usage) {
                             r1 = get_dist(angle_tab[i][0], angle_tab[i][1]) * 0.1;
                             r2 = get_dist(angle_tab[i][1], angle_tab[i][2]) * 0.1;
                             r3 = get_dist(angle_tab[i][0], angle_tab[i][2]) * 0.1;
                        }
                        else {
                             r1 = a_angle_quartic[0][i];
                             r2 = a_angle_quartic[1][i];
                             r3 = a_angle_quartic[2][i];
                        }
                        fprintf(fptop,
                           "   4%12.3f%12.3f%12.3f%12.4f     ;%6s-%6s-%6s\n",
                                r1, r2, r3, a_angle_quartic[3][i],
                                names[(angle_tab[i][0])], names[(angle_tab[i][1])],
                                names[(angle_tab[i][2])]);
                        break;

               case 5:  if(usage) r1 = get_dist(angle_tab[i][0], angle_tab[i][2]) * 0.1;
                        else r1 = a_angle_quartic[1][i];
                        fprintf(fptop, "   5%12.3f%12.4f%12.3f%12.4f     ;%6s-%6s-%6s\n",
                                theta, a_angle_quartic[0][i], r1,
                                a_angle_quartic[2][i], names[(angle_tab[i][0])],
                                names[(angle_tab[i][1])], names[(angle_tab[i][2])]);
                        break;

               case 6:  fprintf(fptop,
                           "   6%12.3f%12.4f%12.4f%12.4f%12.4f%12.4f ;%6s-%6s-%6s\n",
                                theta, a_angle_quartic[0][i], a_angle_quartic[1][i],
                                a_angle_quartic[2][i], a_angle_quartic[3][i],
                                a_angle_quartic[4][i], names[(angle_tab[i][0])],
                                names[(angle_tab[i][1])], names[(angle_tab[i][2])]);
                        break;

               case 8:  fprintf(fptop, "   8%12d%12.4f     ;%6s-%6s-%6s\n",
                                (int)a_angle_quartic[0][i], a_angle_quartic[1][i],
                                names[(angle_tab[i][0])], names[(angle_tab[i][1])],
                                names[(angle_tab[i][2])]);
                        break;

               default: sprintf(mess, "Angle type error, angle number %d.\n", i);
                        my_fatal(FARGS, mess);
                        exit(1);
          }
     }

     return;
}

/* write dihedral force records of simple form
   parameters are:
        fptop                        file pointer to write to
        numdihed                     number of dihedrals
        names                        list of atom names
        usage                        flag to use measured or force field
                                     derived lengths and angles
        amber                        flag for data in kcal/mol and angstroms
        conv                         conversion factor
*/
void write_dihed(FILE *fptop, int numdihed, char **names, int usage,
                 int amber, double conv)
{
     int i, j, k, dihmul;
     double dfc, amult;
     float aphi;

/* dihedrals */
     if(regs_wanted) {
          k = regs_wanted;
          fprintf(fptop, "\n[ dihedrals ]\n");
          fprintf(fptop,
                  "; ai  aj   ak  al  funct    phi0          cp      mult\n");

          for(i = 0; i < numdihed; i++) {
               if(wanted_tors != NULL)
                    if(!wanted_tors[i])
                         continue;     /* only give wanted torsions */
               if((back_order[(torsion_tab[i][0])] < 0) ||
                  (back_order[(torsion_tab[i][1])] < 0) ||
                  (back_order[(torsion_tab[i][2])] < 0) ||
                  (back_order[(torsion_tab[i][3])] < 0))
                         continue;
               if((a_tors_type[i] > 1) &&
                  (a_tors_type[i] != 9))
                         continue;
               dihmul = abs(((int)a_tors_term[i]));
               amult = (double)(abs(a_tors_mult[i]));
               aphi = a_tors_phase[i];
               if(dihmul > 10) dihmul = 1;
               dfc = (double)a_tors_force[i];
               if(amber) dfc = (dfc * conv) / amult;

               fprintf(fptop, "%6d%6d%6d%6d", (back_order[(torsion_tab[i][0])] + lstart2_at),
                       (back_order[(torsion_tab[i][1])] + lstart2_at),
                       (back_order[(torsion_tab[i][2])] + lstart2_at),
                       (back_order[(torsion_tab[i][3])] + lstart2_at));
               switch(a_tors_type[i]) {
                    case 0:  fprintf(fptop,
                 "                                        ; dih %6s-%6s-%6s-%6s\n",
                                names[(torsion_tab[i][0])], names[(torsion_tab[i][1])],
                                names[(torsion_tab[i][2])], names[(torsion_tab[i][3])]);
                             break;

                    case 9:
                    case 1:  fprintf(fptop,
                 "   1%12.3f%12.3f%10d  ; dih %6s-%6s-%6s-%6s\n",
                                     aphi, dfc, dihmul, names[(torsion_tab[i][0])],
                                     names[(torsion_tab[i][1])], names[(torsion_tab[i][2])],
                                     names[(torsion_tab[i][3])]);
/* handle cases with multiple torsion entries */
                             if(multors_cnt[i])
                                  for(j = 0; j < multors_cnt[i]; j++) {
                                       fprintf(fptop, "%6d%6d%6d%6d%s",
                                          (back_order[(torsion_tab[i][0])] + lstart2_at),
                                          (back_order[(torsion_tab[i][1])] + lstart2_at),
                                          (back_order[(torsion_tab[i][2])] + lstart2_at),
                                          (back_order[(torsion_tab[i][3])] + lstart2_at),
                                          "   1");
                                       dihmul = abs(((int)multiples[i][j][3]));
                                       if(dihmul > 10) dihmul = 1;
                                       amult = (double)(abs(((int)multiples[i][j][0])));
                                       aphi = multiples[i][j][2];
                                       dfc = (double)multiples[i][j][1];
                                       if(amber) dfc = (dfc * conv) / amult;
                                       fprintf(fptop,
                                          "%12.3f%12.3f%10d  ; dih %6s-%6s-%6s-%6s\n",
                                          aphi, dfc, dihmul,
                                          names[(torsion_tab[i][0])],
                                          names[(torsion_tab[i][1])],
                                          names[(torsion_tab[i][2])],
                                          names[(torsion_tab[i][3])]);
                                  }         /* end for(j = 0; j < multors_cnt[i]; j++) */
                    default: break;
               }             /* end switch(a_tors_type[i]) */
               k--;
               if(!k) break;
          }                  /* end for(i = 0; i < numdihed; i++) */
     }                       /* end if(regs_wanted) */

     if(rbs_wanted) {
          k = rbs_wanted;
          fprintf(fptop, "\n[ dihedrals ]\n");
          fprintf(fptop,
             ";  ai    aj    ak    al funct            c0            c1            c2");
          fprintf(fptop,
                  "            c3            c4            c5\n");
          for(i = 0; i < numdihed; i++) {
               if(wanted_tors != NULL)
                    if(!wanted_tors[i])
                         continue;     /* only give wanted torsions */
               if((back_order[(torsion_tab[i][0])] < 0) ||
                  (back_order[(torsion_tab[i][1])] < 0) ||
                  (back_order[(torsion_tab[i][2])] < 0) ||
                  (back_order[(torsion_tab[i][3])] < 0))
                         continue;
               if(a_tors_type[i] != 3)
                    continue;
               fprintf(fptop, "%6d%6d%6d%6d", (back_order[(torsion_tab[i][0])] + lstart2_at),
                       (back_order[(torsion_tab[i][1])] + lstart2_at),
                       (back_order[(torsion_tab[i][2])] + lstart2_at),
                       (back_order[(torsion_tab[i][3])] + lstart2_at));
               fprintf(fptop,
                  "   3%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f ; dih %6s-%6s-%6s-%6s\n",
                        a_tors_set[0][i], a_tors_set[1][i], a_tors_set[2][i],
                        a_tors_set[3][i], a_tors_set[4][i], a_tors_set[5][i],
                        names[(torsion_tab[i][0])], names[(torsion_tab[i][1])],
                        names[(torsion_tab[i][2])], names[(torsion_tab[i][3])]);
               k--;
               if(!k) break;
          }                  /* end for(i = 0; i < numdihed; i++) */
     }                       /* end if(rbs_wanted) */

     if(fourier_wanted) {
          k = fourier_wanted;
          fprintf(fptop, "\n[ dihedrals ]\n");
          fprintf(fptop,
                  ";  ai    aj    ak    al funct            c0            c1");
          fprintf(fptop,
                  "            c2            c3\n");

          for(i = 0; i < numdihed; i++) {
               if(wanted_tors != NULL)
                    if(!wanted_tors[i])
                         continue;     /* only give wanted torsions */
               if((back_order[(torsion_tab[i][0])] < 0) ||
                  (back_order[(torsion_tab[i][1])] < 0) ||
                  (back_order[(torsion_tab[i][2])] < 0) ||
                  (back_order[(torsion_tab[i][3])] < 0))
                         continue;
               if(a_tors_type[i] != 5)
                    continue;
               fprintf(fptop, "%6d%6d%6d%6d", (back_order[(torsion_tab[i][0])] + lstart2_at),
                       (back_order[(torsion_tab[i][1])] + lstart2_at),
                       (back_order[(torsion_tab[i][2])] + lstart2_at),
                       (back_order[(torsion_tab[i][3])] + lstart2_at));
               fprintf(fptop,
                  "   5%10.5f%10.5f%10.5f%10.5f ; dih %6s-%6s-%6s-%6s\n",
                       a_tors_set[0][i], a_tors_set[1][i], a_tors_set[2][i],
                       a_tors_set[3][i], names[(torsion_tab[i][0])],
                       names[(torsion_tab[i][1])], names[(torsion_tab[i][2])],
                       names[(torsion_tab[i][3])]);
               k--;
               if(!k) break;
          }                  /* end for(i = 0; i < numdihed; i++) */
     }                       /* end if(fourier_wanted) */

     if(table_wanted) {
          k = table_wanted;
          fprintf(fptop, "\n[ dihedrals ]\n");
          fprintf(fptop,
                  ";  ai    aj    ak    al funct            n0            n1");
          for(i = 0; i < numdihed; i++) {
               if(wanted_tors != NULL)
                    if(!wanted_tors[i])
                         continue;     /* only give wanted torsions */
               if((back_order[(torsion_tab[i][0])] < 0) ||
                  (back_order[(torsion_tab[i][1])] < 0) ||
                  (back_order[(torsion_tab[i][2])] < 0) ||
                  (back_order[(torsion_tab[i][3])] < 0))
                         continue;
               if(a_tors_type[i] != 8)
                    continue;
               fprintf(fptop, "%6d%6d%6d%6d", (back_order[(torsion_tab[i][0])] + lstart2_at),
                       (back_order[(torsion_tab[i][1])] + lstart2_at),
                       (back_order[(torsion_tab[i][2])] + lstart2_at),
                       (back_order[(torsion_tab[i][3])] + lstart2_at));
               fprintf(fptop, "   8%12d%12.4f     ; dih %6s-%6s-%6s-%6s\n",
                       (int)a_tors_set[0][i], a_tors_set[1][i],
                       names[(torsion_tab[i][0])], names[(torsion_tab[i][1])],
                       names[(torsion_tab[i][2])], names[(torsion_tab[i][3])]);
               k--;
               if(!k) break;
          }                  /* end for(i = 0; i < numdihed; i++) */
     }                       /* end if(table_wanted) */

     return;
}

/* write improper force records of simple form
   parameters are:
        fptop                        file pointer to write to
        numimpro                     number of impropers
        numdihed                     number of dihedrals to search for type 2
        names                        list of atom names
        usage                        flag to use measured or force field
                                     derived lengths and angles
        amber                        flag for data in kcal/mol and angstroms
        conv                         conversion factor
*/
void write_impropers(FILE *fptop, int numimpro, int numdihed, char **names,
                     int usage, int amber, double conv)
{
     int i, impmul;
     double ifc;
     float iphi;
     char atwo[]   = { "   2" };
     char aone[]   = { "   1" };
     char foursp[] = { "    " };
     char *ifunct;
     char *imultp;

/* amber type force fields just append improper dihedrals to the others
   as type 1 dihedrals with multiplicity of 1.  Only have standard
   improper dihedrals.  Ring flatness dihedrals are left with all
   other dihedrals.
*/
     if(amber) {
          if(numimpro)
               for(i = 0; i < numimpro; i++) {
                    if((back_order[(improperid[i][0])] < 0) ||
                       (back_order[(improperid[i][1])] < 0) ||
                       (back_order[(improperid[i][2])] < 0) ||
                       (back_order[(improperid[i][3])] < 0))
                              continue;
                    if(usage) iphi = meas_impro[i];
                    else iphi = a_impr_phase[i];
                    ifc = (double)a_impr_force[i] * conv;
                    impmul = abs(((int)a_impr_term[i]));
                    if(impmul > 10) impmul = 1;
/* change from amber order */
                    fprintf(fptop, "%6d%6d%6d%6d",
                            (back_order[(improperid[i][2])] + lstart2_at),
                            (back_order[(improperid[i][0])] + lstart2_at),
                            (back_order[(improperid[i][1])] + lstart2_at),
                            (back_order[(improperid[i][3])] + lstart2_at));
                    if(abs(a_impr_mult[i]) > 999) 
                         fprintf(fptop,
                            "   1                                    ; imp %6s-%6s-%6s-%6s\n",
                                 names[(improperid[i][2])], names[(improperid[i][0])],
                                 names[(improperid[i][1])], names[(improperid[i][3])]);
                    else fprintf(fptop,
                            "   1%12.3f%12.3f%10d  ; imp %6s-%6s-%6s-%6s\n",
                                 iphi, ifc, impmul, names[(improperid[i][2])],
                                 names[(improperid[i][0])], names[(improperid[i][1])],
                                 names[(improperid[i][3])]);
               }

          return;
     }

/* non-amber impropers, just like dihedrals.
   gmx impropers use function 2 and no multiplicities.
   oplsaa impropers use function 1 and multiplicity 2
*/
     switch(use_amber) {
          case 5:  fprintf(fptop, "\n[ dihedrals ]\n");
                   fprintf(fptop,
                  "; ai  aj   ak  al  funct    phi0          cp      mult\n");

                   ifunct = &aone[0];
                   imultp = &atwo[0];
                   break;
          case 4:  fprintf(fptop, "\n[ dihedrals ]\n");
                   fprintf(fptop,
                  ";  ai    aj    ak    al funct           c0           c1\n");

                   ifunct = &atwo[0];
                   imultp = &foursp[0];
          default: break;
     }
     for(i = 0; i < numdihed; i++) {
          if((back_order[(torsion_tab[i][0])] < 0) ||
             (back_order[(torsion_tab[i][1])] < 0) ||
             (back_order[(torsion_tab[i][2])] < 0) ||
             (back_order[(torsion_tab[i][3])] < 0))
                    continue;
          if(a_tors_type[i] != 2)
               continue;
          if(usage)
               iphi = meas_tors[i];
          else
               iphi = a_tors_phase[i];
          ifc = (double)a_tors_force[i];
          fprintf(fptop, "%6d%6d%6d%6d", (back_order[(torsion_tab[i][0])] + lstart2_at),
                  (back_order[(torsion_tab[i][1])] + lstart2_at),
                  (back_order[(torsion_tab[i][2])] + lstart2_at),
                  (back_order[(torsion_tab[i][3])] + lstart2_at));
          fprintf(fptop,
                       "%s%12.3f%12.3f%s        ; imp %6s-%6s-%6s-%6s\n",
                       ifunct, iphi, ifc, imultp, names[(torsion_tab[i][0])],
                       names[(torsion_tab[i][1])], names[(torsion_tab[i][2])],
                       names[(torsion_tab[i][3])]);
     }

     if(numimpro)
          for(i = 0; i < numimpro; i++) {
               if((back_order[(improperid[i][0])] < 0) ||
                  (back_order[(improperid[i][1])] < 0) ||
                  (back_order[(improperid[i][2])] < 0) ||
                  (back_order[(improperid[i][3])] < 0))
                         continue;
               if(usage)
                    iphi = meas_impro[i];
               else
                    iphi = a_impr_phase[i];
               ifc = (double)a_impr_force[i];
/* change from amber order */
               fprintf(fptop, "%6d%6d%6d%6d", (back_order[(improperid[i][2])] + lstart2_at),
                       (back_order[(improperid[i][0])] + lstart2_at),
                       (back_order[(improperid[i][1])] + lstart2_at),
                       (back_order[(improperid[i][3])] + lstart2_at));
               if(abs(a_impr_mult[i]) > 999) 
                    fprintf(fptop,
                            "%s                                    ; imp %6s-%6s-%6s-%6s\n",
                            ifunct, names[(improperid[i][2])], names[(improperid[i][0])],
                            names[(improperid[i][1])], names[(improperid[i][3])]);
               else
                    fprintf(fptop,
                            "%s%12.3f%12.3f%s        ; imp %6s-%6s-%6s-%6s\n",
                            ifunct, iphi, ifc, imultp, names[(improperid[i][2])],
                            names[(improperid[i][0])], names[(improperid[i][1])],
                            names[(improperid[i][3])]);
          }

     return;
}

/* write restraint records
   parameters are:
        fptop                        file pointer to write to
        tot_rstr                     total number of restraints
        names                        list of atom names
*/
void write_restraints(FILE *fptop, int tot_rstr, char **names)
{
     int i, j, rrstr1, rrstr2, idone;

/* write distance range restraints information */
     if(rstr_range) {
          idone = 1;
          j = rstr_range;
          for(i = 0; i < tot_rstr; i++) {
               if(restr[i][0] != 2) continue;
               if((back_order[(restr[i][1] - 1)] < 0) ||
                  (back_order[(restr[i][2] - 1)] < 0)) {
                         j--;
                         if(!j) break;
                         continue;
               }
               if(idone) {
                    fprintf(fptop, "\n#ifdef DISTRANGE\n");
                    fprintf(fptop, "[ distance_restraints ]\n");
                    fprintf(fptop,
                       ";   ai  aj   type  index  type'  low  up1  up2  fac\n");
                    idone = 0;
               }
               rrstr1 = back_order[(restr[i][1] - 1)] + lstart2_at;
               rrstr2 = back_order[(restr[i][2] - 1)] + lstart2_at;
               fprintf(fptop,
                  "  %6d%6d%s%6d%s%12.5f%12.5f%12.5f%s       ;%6s-%6s\n",
                       rrstr1, rrstr2,  "   1", i, "   1",
                       (frest[i][0] * 0.1), (frest[i][1] * 0.1),
                       ((frest[i][1] * 0.1) + 1.0), "   1.0",
                       names[(restr[i][1] - 1)], names[(restr[i][2] - 1)]);
               j--;
               if(!j) break;
          }

          if(!idone)
               fprintf(fptop, "#endif\n");
     }

/* write torsion restraints information */
     if(rstr_torsn) {
          idone = 1;
          j = rstr_torsn;
          for(i = 0; i < tot_rstr; i++) {
               if(restr[i][0] != 4) continue;
               if((back_order[(restr[i][1] - 1)] < 0) ||
                  (back_order[(restr[i][2] - 1)] < 0) ||
                  (back_order[(restr[i][3] - 1)] < 0) ||
                  (back_order[(restr[i][4] - 1)] < 0)) {
                         j--;
                         if(!j) break;
                         continue;
               }
               if(idone) {
                    fprintf(fptop, "\n#ifdef DIHEDREST\n");
                    fprintf(fptop, "[ dihedral_restraints ]\n");
                    fprintf(fptop,
                       ";   ai  aj  ak  al   type  label  phi  dphi  kfac  power\n");
                    idone = 0;
               }
               fprintf(fptop,
                  " %6d%6d%6d%6d%s%6d%12.3f%s%s%s       ;%6s-%6s-%6s-%6s\n",
                       (back_order[(restr[i][1] - 1)] + lstart2_at),
                       (back_order[(restr[i][2] - 1)] + lstart2_at),
                       (back_order[(restr[i][3] - 1)] + lstart2_at),
                       (back_order[(restr[i][4] - 1)] + lstart2_at), "   1", i,
                       frest[i][0], "         0.0",  "   1", "   2",
                       names[(restr[i][1] - 1)], names[(restr[i][2] - 1)],
                       names[(restr[i][3] - 1)], names[(restr[i][4] - 1)]);
               j--;
               if(!j) break;
          }
          if(!idone)
               fprintf(fptop, "#endif\n");
     }

     return;
}

/* assign charge group numbers
   parameters are:
        numatoms                     number of atoms
        residue                      the residue name for all the atoms
*/
void cgnr_assgn(int num_atoms, char *residue)
{
     int i, j, k, jend, jatom, grp_num, new_grp, redo, assgnd;
     char curr[6];

     cgnr_grp = (int *)allocator(num_atoms, sizeof(int), FARGS);
     for(i = 0; i < num_atoms; i++)
          cgnr_grp[i] = -1;

     curr[0] = NULLCHAR;
     grp_num = 0;
     new_grp = 0;
     redo = 1;

     if(residue[0] != NULLCHAR)
          grp_num = 1;

     for(i = 0; i < num_atoms; i++) {
          if(atom_order[i] < 0) continue;
          if(cgnr_grp[(atom_order[i])] > -1)
               continue;
/* each residue starts a new charge group */
          if(residue[0] == NULLCHAR)
               if(strncmp(atom_residue[(atom_order[i])], curr, 5)) {
                    if(!new_grp) grp_num++;
                    strncpy(curr, atom_residue[(atom_order[i])], 5);
                    curr[5] = NULLCHAR;
                    new_grp = 1;
               }
          if((atom_atno[(atom_order[i])] != 1)  &&
             (atom_atno[(atom_order[i])] != 7)  &&
             (atom_atno[(atom_order[i])] != 8)) {
                    cgnr_grp[(atom_order[i])] = grp_num;
                    if(!bond_count[(atom_order[i])]) continue;
                    jend = bond_count[(atom_order[i])];
                    jatom = atom_order[i];
                    for(j = 0; j < jend; j++) {
                         if(atom_order[(connect[jatom][j])] < 0)
                              continue;
                         if(cgnr_grp[(connect[jatom][j])] > -1)
                              continue;
                         if(residue[0] == NULLCHAR)
                              if(strncmp(atom_residue[(connect[jatom][j])],
                                         curr, 5))
                                   continue;
                         if((atom_atno[(connect[jatom][j])] != 1) &&
                            (atom_atno[(connect[jatom][j])] != 7) &&
                            (atom_atno[(connect[jatom][j])] != 8))
                                   continue;
                         cgnr_grp[(connect[jatom][j])] = grp_num;
                    }
                    grp_num++;
                    new_grp = 1;
          }
     }

     for(k = 0; k < num_atoms; k++) {
          redo = 0;
          for(i = 0; i < num_atoms; i++) {
               if(cgnr_grp[i] > -1) continue;
               if(!bond_count[i]) {
                    cgnr_grp[i] = grp_num;
                    grp_num++;
                    continue;
               }
               assgnd = 0;
               for(j = 0; j < bond_count[i]; j++) {
                    if(cgnr_grp[(connect[i][j])] < 0)
                         continue;
/* each residue starts a new charge group */
                    if(residue[0] == NULLCHAR)
                         if(strncmp(atom_residue[i],
                                    atom_residue[(connect[i][j])], 5))
                              continue;
                    cgnr_grp[i] = cgnr_grp[(connect[i][j])];
                    assgnd = 1;
                    break;
               }
               if(!assgnd) redo = 1;
          }

          if(!redo) break;
     }

     if(redo)
          for(i = 0; i < num_atoms; i++) {
               if(cgnr_grp[i] > -1) continue;
               cgnr_grp[i] = grp_num;
               grp_num++;
          }

     return;
}

float get_dist(int point1, int point2)
{
     double delx, dely, delz;
     float my_dist;

          delx = atom_x[point1] - atom_x[point2];
          dely = atom_y[point1] - atom_y[point2];
          delz = atom_z[point1] - atom_z[point2];
          my_dist = sqrt((delx * delx) + (dely * dely) + (delz * delz));

     return(my_dist);
}
