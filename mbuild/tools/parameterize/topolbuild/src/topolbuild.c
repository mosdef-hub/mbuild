/* Topology Builder
   Bruce D. Ray
   IUPUI Physics Dept.
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

     Version 1.3
             Main program with functions to print usage and version
             information and to handle amber / gaff processing calls.
             Produces gromacs topology and coordinates files from
             a Tripos Sybyl mol2 file, and the user's choice of
             amber / gaff, glycam, Tripos, gromacs, or oplsaa force
             field

   NOTICE: Portions of this are derivative work based on study of:
      1. The routines found GROMACS version 3.3.1 written by David
         van der Spoel, Erik Lindahl, Berk Hess, and others and
         copyright (c) 1991-2000, University of Groningen, The Netherlands.
         Copyright (c) 2001-2004, The GROMACS development team

      2. The routines found in antechamber 1.27 particularly in:
         a. program atomtype, version 1.0, dated October, 2001 by
            Junmei Wang, Department of Pharmaceutical Chemistry,
            School of Pharmacy, University of California, San
            Francisco, CA  94143
         b. programs parmchk and prepgen, version 1.0, dated
            October, 2001 by Junmei Wang, Department of
            Pharmaceutical Chemistry, School of Pharmacy,
            University of California, San Francisco, CA  94143
         c. function file rings.c, by Junmei Wang, Department of
            Pharmaceutical Chemistry, School of Pharmacy, University
            of California, San Francisco, CA  94143
         d. program charmgen, version 2.0, dated June, 2004 by Victor
            E. Bazterra, Center for High Performance Computing, University
            of Utah, Salt Lake City, UT  84114
         e. program bondtype, version 1.0, dated October, 2001 by
            Junmei Wang, Department of Pharmaceutical Chemistry,
            School of Pharmacy, University of California,
            San Francisco, CA  94143
         f. function file mol2.c  by Junmei Wang, Department of
            Pharmaceutical Chemistry, School of Pharmacy, University
            of California, San Francisco, CA  94143

      3. Routines found in AMBCONV by Filip Ryjacek, 2002 as distributed
         at the GROMACS web site.

      4. The routines found in Dock 6.1 file dockmol.cpp version dated
         December 2006, Molecular Design Institute, University of
         California, San Francisco, CA  94143
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "compare_FF.h"
#include "judge_bond.h"
#include "judge_types.h"
#include "mainchain.h"
#include "param_tables.h"
#include "readmol2.h"
#include "ring_detect.h"
#include "read_amber.h"
#include "use_amber.h"
#include "use_gaff.h"
#include "adjustments.h"
#include "similars.h"
#include "Tripos_FF.h"
#include "gromacs_FF.h"
#include "gromacs_use.h"
#include "multors.h"

#define MAXCHAR 256
#define MAXNAME 128
#define NULLCHAR 0
#define MAXATOM 2000
#define NUMFORCEFIELDS 6
#define NUMOPTS 17

void usage(char *name);
void do_an_amber(char *ante_name, int use_amber, char *leap_name,
                 int num_atoms, int num_bonds, char *tp_amber,
                 int angle_max, int tors_max);
void version(char *name);
void regs_log(int i_count);
void rbs_log(int i_count);
void fourier_log(int i_count);
void table_log(int i_count);

extern int lstart_at;
extern int res_start_at;
extern int rstart_at;
extern int lstart2_at;
extern int use_amber;
extern int center_it;

extern char mess[MAXCHAR];
extern char ambername[MAXCHAR];
extern char mycorr_name[MAXCHAR];

extern char *use_define[6];
extern char *in_command;
char vers_no[] = { "version 1.3" };

/* force field define lines for topology */
static char *amber_def[]  = { "_FF_AMBER", "_FF_AMBER" };
static char *gaff_def[]   = { "_FF_GAFF", "_FF_AMBER" };
static char *gromos_def[] = { "_FF_GROMOS96", "_FF_GROMOS" };
static char *tripos_def[] = { "_FF_TRIPOS", "_FF_TRIPOS" };
static char *glycam_def[] = { "_FF_AMBER", "_FF_AMBER94" };
static char *oplsaa_def[] = { "_FF_OPLS", "_FF_OPLSAA" };
       
/* end force field define lines for topology */

/* water and ion include lines for topology */
static char gaff_spce[] = { "gaff_spce.itp" };
static char norm_spce[] = { "spce.itp" };
static char gaff_ions[] = { "ions_gaff.itp" };
static char norm_ions[] = { "ions.itp" };
/* end water and ion include lines for topology */

/* defaults lines for topology */
static char gaff_default[]    = { "  1               2               yes             0.5     0.8333" };
static char gmx_default[]     = { "  1               1               no              1.0     1.0" };
static char oplsaa_default[]  = { "  1               3               yes             0.5     0.5" };
/* end defaults lines for topology */

/* atom type common initial character string lines for topology */
static char opls_head[] = { "opls_" };
/* end atom type common initial character string lines for topology */

/* This lets me add forcefields later */
typedef struct {
     char forcefield[8];          /* force field name */
     int ff_len;                  /* length of force field name */
     int ff_use;                  /* force field numeric designation */
     int ff_correl;               /* option for later use */
     int ff_adjust;               /* do a units conversion if 1 */
     int ua_ff;                   /* united atoms force field if 1 */
     int subtyped;                /* multiple subtypes of force field if 1 */
} FORCEFIELD;

/* force field definitions */
static FORCEFIELD myforcefields[] = {
     { "gaff",   4, 0, 0, 1, 0, 0 },
     { "amber",  5, 1, 0, 1, 0, 1 },
     { "tripos", 6, 2, 0, 1, 0, 0 },
     { "glycam", 6, 3, 0, 1, 0, 1 },
     { "gmx",    3, 4, 0, 0, 1, 1 },
     { "oplsaa", 6, 5, 0, 0, 0, 0 }
};
/* end force field definitions */

int main(int argc, char *argv[])
{
     char filename[MAXCHAR], inputname[MAXCHAR], itp_name[MAXCHAR];
     char logname[MAXCHAR], outmol[MAXCHAR], leap_name[MAXCHAR];
     char residue[6], dummy[6], amberff[12], ante_name[MAXCHAR];
     char gro_name[MAXCHAR], tp_amber[6], cap_amber[6];
     char *cslash, *dnext, *d, *whichitp;
     int num_atoms, num_bonds, i, j, max_rng, use_mine;
     int angle_max, ang_count, tors_max, tors_count;
     int num, num2, count1, count2, nmbr, do_rename, all_renum;
     int found, gro_only, by_correl, parm_adj, do_ua, k;
     int purge_level, iprtit;
     int FFmol, FFbond, FFangle, FFtors, FFimpro, FFvdw, FFequiv;
     int FFeqlen;
     char *my_options[NUMOPTS] = { "-h", "-H", "-v", "-V", "-n",
                                   "-dir", "-ff", "-r", "-meas",
                                   "-rename", "-gro", "-last",
                                   "-resnum", "-renumall",
                                   "-purge", "-move", "-charge" };

     if(argc < 2){
          fprintf(stderr, "Need parameters to run.\n");
          usage(argv[0]);
          exit(1);
     }

     filename[0] = NULLCHAR;
     ante_name[0] = NULLCHAR;
     amberff[0] = NULLCHAR;
     tp_amber[0] = NULLCHAR;
     residue[0] = NULLCHAR;
     use_mine = 0;
     use_amber = -1;
     do_rename = 0;
     lstart_at = 1;
     res_start_at = 0;
     rstart_at = 0;
     lstart2_at = 1;
     all_renum = 0;
     gro_only = 0;
     purge_level = 1;         /* default to purge of 1 */
     whichitp = NULL;
     do_ua = 0;
     center_it = 0;           /* do not center molecule */
     charge_set = 0;          /* do not use force field based charges */

     for (i = 1; i < argc; i++) {
          found = 0;
          for(j = 0; j < NUMOPTS; j++) 
               if(!strcmp(argv[i], my_options[j])) {
                    found = j + 1;
                    break;
               }
          switch(found) {
               case  1:                           /* -h */
               case  2:  usage(argv[0]);          /* -H */
                         break;

               case  3:                           /* -v */
               case  4:  version(argv[0]);        /* -V */
                         break;

               case  5:                           /* -n */
                         if(strlen(argv[(i + 1)]) > (MAXCHAR - 10)) {
                               sprintf(mess,
                                       "File name, %s, too long.\n",
                                       argv[(i + 1)]);
                               my_fatal(FARGS, mess);
                               exit(1);
                         }
                         strncpy(filename, argv[i+1], (MAXCHAR - 10));
                         filename[(MAXCHAR - 10)] = NULLCHAR;
                         i++;
                         break;

               case  6:                           /* -dir */
                         if(strlen(argv[(i + 1)]) > (MAXCHAR - 40)) {
                              sprintf(mess,
                                "Forcefield home directory name, %s, too long.\n",
                                      argv[(i + 1)]);
                              my_fatal(FARGS, mess);
                              exit(1);
                         }
                         strncpy(ante_name, argv[i+1], (MAXCHAR - 40));
                         ante_name[(MAXCHAR - 40)] = NULLCHAR;
                         i++;
                         break;

               case  7:                           /* -ff */
                         strncpy(amberff, argv[(i + 1)], 10);
                         amberff[10] = NULLCHAR;
                         i++;
                         break;

               case  8:                           /* -r */
                         strncpy(residue, argv[(i + 1)], 5);
                         residue[5] = NULLCHAR;
                         i++;
                         break;

               case  9:                           /* -meas */
                         use_mine = 1;
                         break;

               case 10:                           /* -rename */
                         do_rename = 1;
                         break;

               case 11:                           /* -gro */
                         gro_only = 1;
                         break;

               case 12:                           /* -last */
                         if(strspn(argv[(i + 1)], "0123456789")
                            != strlen(argv[(i + 1)])) {
                                 my_fatal(FARGS,
                                    "Non digit characters in -last specification\n");
                                 exit(1);
                         }
                         sscanf(argv[(i + 1)], "%d", &lstart_at);
                         lstart_at++;
                         i++;
                         break;

               case 13:                           /* -resnum */
                         if(strspn(argv[(i + 1)], "0123456789")
                            != strlen(argv[(i + 1)])) {
                                 my_fatal(FARGS,
                                    "Non digit characters in -resnum specification\n");
                                 exit(1);
                         }
                         sscanf(argv[(i + 1)], "%d", &res_start_at);
                         i++;
                         break;

               case 14:                           /* -renumall */
                         all_renum = 1;
                         break;

               case 15:                           /* -purge [level] */
                         purge_level = 1;
                         if(i == (argc - 1)) break;
                         if(strlen(argv[(i + 1)]) > 1) break;
                         if(strspn(argv[(i + 1)], "012345")
                            != strlen(argv[(i + 1)])) break;
                         purge_level = atoi(argv[(i + 1)]);
                         i++;
                         break;

               case 16:                           /* -move */
                         center_it = 1;
                         break;

               case 17:                           /* -charge */
                         charge_set = 1;
                         break;

               default:  sprintf(mess, "Unrecognized option %s\n",
                                 argv[i]);
                         my_fatal(FARGS, mess);
                         exit(1);
          }
     }

     if(!filename[0] || !ante_name[0] || !amberff[0]) {
          fprintf(stderr, "Need correct parameters\n");
          usage(argv[0]);
          exit(1);
     }

     dnext = &filename[0];
     while(1) {
          if((cslash = strchr(dnext, '/')) == NULL) break;
          dnext = cslash;
          dnext++;
     }
     strcpy(inputname, filename);
     strcat(inputname, ".mol2");
     strcpy(logname, dnext);
     strcat(logname, ".log");
     strcpy(outmol, dnext);
     strcat(outmol, "MOL.mol2");
     strcpy(itp_name, dnext);
     strcat(itp_name, ".top");
     strcpy(gro_name, dnext);
     strcat(gro_name, ".gro");

     if((stdlog = fopen(logname, "w")) == NULL) {
          sprintf(mess, "Cannot open log file %s\n", logname);
          my_fatal(FARGS, mess);
          exit(1);
     }

/* log the run command */
     fprintf(stdlog, "topolbuild %s\nCommand line:\n",
             vers_no);
     j = 0;
     for(i = 0; i < argc; i++)
          j += strlen(argv[i]) + 1;
     j++;
     
     in_command = (char *)allocator(j, sizeof(char), FARGS);
     strcpy(in_command, argv[0]);
     for(i = 1; i < argc; i++) {
          strcat(in_command, " ");
          strcat(in_command, argv[i]);
     }

     fprintf(stdlog, "     %s\n\n", in_command);

/* lets me add forcefields later */
     for(i = 0; i < NUMFORCEFIELDS; i++)
          if(!strncmp(amberff, myforcefields[i].forcefield,
                      myforcefields[i].ff_len)) {
               use_amber = myforcefields[i].ff_use;
               by_correl = myforcefields[i].ff_correl;
               parm_adj  = myforcefields[i].ff_adjust;
               do_ua     = myforcefields[i].ua_ff;
/* check force field subtypes */
               if(myforcefields[i].subtyped) {
                    strncpy(tp_amber,
                            &amberff[(myforcefields[i].ff_len)], 5);
                    tp_amber[5] = NULLCHAR;
               }
               break;
          }

     if(use_amber < 0) {
          sprintf(mess,
                  "Forcefield %s not recognized\n", amberff);
          my_fatal(FARGS, mess);
          exit(1);
     }

     use_define[4] = NULL;
     no_oppose = 0;

     switch(use_amber) {
/* GAFF AA force field */
          case 0:  use_define[0] = gaff_def[0];
                   use_define[1] = gaff_def[1];
                   use_define[2] = &gaff_spce[0];
                   use_define[3] = &gaff_ions[0];
                   use_define[5] = &gaff_default[0];
                   break;

/* Amber AA force field */
          case 1:  use_define[0] = amber_def[0];
                   i = strlen(amber_def[1]) + strlen(tp_amber) + 2;
                   use_define[1] = (char *)allocator(i, sizeof(char),
                                                     FARGS);
                   strcpy(use_define[1], amber_def[1]);
                   strcpy(cap_amber, tp_amber);
                   uppercase(cap_amber);
                   strcat(use_define[1], cap_amber);
                   use_define[2] = &gaff_spce[0];
                   use_define[3] = &gaff_ions[0];
                   use_define[5] = &gaff_default[0];
                   break;

/* Tripos AA force field */
          case 2:  use_define[0] = tripos_def[0];
                   use_define[1] = tripos_def[1];
                   use_define[2] = &gaff_spce[0];
                   use_define[3] = &gaff_ions[0];
                   use_define[5] = &gaff_default[0];
                   break;

/* Glycam AA force field */
          case 3:  use_define[0] = glycam_def[0];
                   use_define[1] = glycam_def[1];
                   use_define[2] = &gaff_spce[0];
                   use_define[3] = &gaff_ions[0];
                   use_define[5] = &gaff_default[0];
                   break;

/* Gromos / Gromacs UA force fields */
          case 4:  use_define[0] = gromos_def[0];
                   i = strlen(gromos_def[1]) + strlen(tp_amber) + 2;
                   use_define[1] = (char *)allocator(i, sizeof(char),
                                                     FARGS);
                   strcpy(use_define[1], gromos_def[1]);
                   strcpy(cap_amber, tp_amber);
                   uppercase(cap_amber);
                   strcat(use_define[1], cap_amber);
                   use_define[2] = &norm_spce[0];
                   use_define[3] = &norm_ions[0];
                   use_define[5] = &gmx_default[0];
                   no_oppose = 1;
                   break;

/* OPLS-AA force field. */
          case 5:  use_define[0] = oplsaa_def[0];
                   use_define[1] = oplsaa_def[1];
                   use_define[2] = &norm_spce[0];
                   use_define[3] = &norm_ions[0];
                   use_define[4] = &opls_head[0];
                   use_define[5] = &oplsaa_default[0];
                   break;

          default: my_fatal(FARGS,
                            "Specified force field not implemented.\n");
                   exit(1);
     }                                /* end switch(use_amber) */
     switch(use_amber) {              /* gaff and amber use antechamber files */
          case 0:
          case 1:                     /* glycam is an amber type */
          case 3:  strcpy(leap_name, ante_name);
                   strcat(ante_name, "/dat/antechamber/");
                   strcat(leap_name,
                          "/dat/leap/parm/");     /* just fall through to break */
          default: break;
     }                                /* end switch(use_amber) */
 
     if((use_amber != 5) && (charge_set == 1)) {
          sprintf(mess,
                  "Option -charge ignored for forcefield %s.",
                  amberff);
          my_warning(FARGS, mess);
          charge_set = 0;
     }
     if((!res_start_at) && (lstart_at != 1)) {
          my_fatal(FARGS,
            "Must give last residue number used when giving last atom number used.\n");
          exit(1);
     }

     if(res_start_at && (lstart_at == 1)) {
          my_fatal(FARGS,
            "Must give last atom number used > 1 when giving last residue number used.\n");
          exit(1);
     }

     if(all_renum && !res_start_at) {
          my_fatal(FARGS,
            "-renumall option requires -last and -resnum options be given.\n");
          exit(1);
     }

     if(all_renum && center_it) {
          my_fatal(FARGS,
            "Options -renumall and -move conflict.\n");
          exit(1);
     }

     if(all_renum) {
          rstart_at = res_start_at;
          lstart2_at = lstart_at;
     }

     readmol2(inputname, &num_atoms, &num_bonds);

     atom_residue_order(num_atoms, residue);

     new_names = atom_name;
     fprintf(stdlog, "Molecule Name %s\n", mol_name);
     fprintf(stdlog, "Model Type %s\n", model_type);
     fprintf(stdlog, "Charge Method %s\n", charge_type);
     fprintf(stdlog, "Has %d atoms and %d bonds\n", num_atoms,
             num_bonds);

     angle_max = num_bonds;
     find_angles(0, angle_max, num_atoms, num_bonds, &ang_count);
     angle_max = ang_count;
     find_angles(1, angle_max, num_atoms, num_bonds, &ang_count);

     tors_max = num_bonds;
     find_torsions(0, tors_max, num_atoms, num_bonds, &tors_count);
     tors_max = tors_count;
     find_torsions(1, tors_max, num_atoms, num_bonds, &tors_count);

     measure_bondlen(num_bonds);
     measure_angles(angle_max);
     measure_torsions(tors_max);

     max_rng = 3 * num_bonds/4;
     detect_ring(num_atoms, num_bonds, max_rng, &num_rngs);
     loc_improper(num_atoms, do_ua);

     measure_improp(improper_num);

     planar_angles(ang_count);

     judge_bond(num_bonds, num_atoms, num_rngs);
     molec_alloc(num_atoms, num_bonds, ang_count, tors_count,
                 improper_num);     /* allocate rest of molecule storage */

     set_mainchain(num_atoms);
     if(do_rename) renameit(num_atoms, num_rngs);

     writemol2(outmol, num_atoms, num_bonds, new_names);

     fprintf(stdlog, "\n Atom Record\n");
     fprintf(stdlog, "     no.    name");
     if(new_names != atom_name)
          fprintf(stdlog, " new name");
     fprintf(stdlog, "   type residue  Sym seg atno            coordinates");
     fprintf(stdlog, "             charge     mass  bonds\n");

     for(i = 0; i < num_atoms; i++) {
          strncpy(dummy, atom_residue[i], 5);
          dummy[5] = NULLCHAR;
          if((cslash = strpbrk(dummy, "0123456789")) != NULL) {
               k = strlen(dummy);
               if((cslash != &dummy[0]) &&
                  ((d = strpbrk(&dummy[k-1], "0123456789")) != NULL))
                         *cslash = NULLCHAR;
          }

          fprintf(stdlog, "%7d %8s", (i + lstart_at), atom_name[i]);
          if(new_names != atom_name)
               fprintf(stdlog, " %8s", new_names[i]);
          fprintf(stdlog,
       " %6s  %5s   %2s  %2d  %3d %10.4f %10.4f %10.4f %9.4f %9.4f  %2d",
                  atom_type[i], dummy, atom_symbl[i], atom_segno[i],
                  atom_atno[i], atom_x[i], atom_y[i], atom_z[i],
                  atom_charge[i], atom_mass[i], bond_count[i]);
          for(j = 0; j < main_num; j++)
               if(main_chain[j] == i) {
                    fprintf(stdlog, "    main");
                    break;
               }
          fprintf(stdlog, "\n");
     }

     if(num_rngs) {
/* write out a rings found report */
          fprintf(stdlog,
        "\n===============================================================\n");
          fprintf(stdlog,
        "----------------------ring property (I)------------------------\n");
	 
          nmbr = 1;
          for(i = 0; i < num_rngs; i++) {
               if(ring_num[i] == 0) continue;
               fprintf(stdlog, "\n        Ring %d\n\n", nmbr++);
               for(j = 0; j < ring_num[i]; j++)
                    fprintf(stdlog, "%5d %5d %5d %5s\n", (i + 1), (j + 1),
                            (ring_atomno[i][j] + lstart_at),
                            new_names[(ring_atomno[i][j])] );
          }              /* end for(i = 0; i < num_rngs; i++) */

          fprintf(stdlog,
       "\n===============================================================\n");
          fprintf(stdlog,
       "---------------------ring property (II)------------------------\n");

          for(i = 0; i < num_rngs; i++) {
               if(ring_num[i] == 0) continue;
               for(j = 0; j < ring_num[i]; j++)
                    fprintf(stdlog,
              "atom[%5d] (%-4s) belongs to one member of %d-membered ring (No %2d)\n",
                            (ring_atomno[i][j] + lstart_at),
                            new_names[(ring_atomno[i][j])],
                            ring_num[i], (i + 1));
          }              /* end for(i = 0; i < num_rngs; i++) */
     }                   /* end if(num_rngs) */

     fprintf(stdlog,
       "\n===============================================================\n");
     fprintf(stdlog,
       "---------------------ring property (III)-----------------------\n");

     for(i = 0; i < num_atoms; i++) {
          for(j = 1; j <= 10; j++) {
               if(arom_atomno[i][j] > 0)
                    fprintf(stdlog,
             "atom[%5d] (%-4s) involved in %d %d-member ring(s)\n",
                            (i + lstart_at), new_names[i],
                            arom_atomno[i][j], j);
          }              /* end for(j = 1; j <= 10; j++) */
          if(arom_num[i] > 0)
               fprintf(stdlog,
                  "atom[%5d] (%-4s) is not in any ring (number[%d]=%d)\n",
                       (i + lstart_at), new_names[i], (i + lstart_at),
                       arom_num[i]);
     }                   /* end for(i = 0; i < num_atoms; i++) */

     if(num_rngs) {
          fprintf(stdlog,
    "\n===============================================================\n");
          fprintf(stdlog,
    "----------------------aromatic property------------------------\n");
          for(i = 0; i < num_atoms; i++) {
               if(ar_set[i][1] >= 1)
                    fprintf(stdlog,
                "atom[%5d] (%-4s) is in %d pure aromatic ring(s) (AR1)\n",
                            (i + lstart_at), new_names[i], ar_set[i][1]);
               if(ar_set[i][2] >= 1)
                    fprintf(stdlog,
                "atom[%5d] (%-4s) is in %d planar ring(s) (AR2)\n",
                            (i + lstart_at), new_names[i], ar_set[i][2]);
               if(ar_set[i][3] >= 1)
                    fprintf(stdlog,
            "atom[%5d] (%-4s) is in %d planar ring(s) that have \"outside\" bonds (AR3)\n",
                            (i + lstart_at), new_names[i],
                            ar_set[i][3]);
               if(ar_set[i][4] >= 1)
                    fprintf(stdlog,
            "atom[%5d] (%-4s) is in %d non-planar ring(s) (AR4)\n",
                            (i + lstart_at), new_names[i],
                            ar_set[i][4]);
               if(ar_set[i][5] >= 1)
                    fprintf(stdlog,
            "atom[%5d] (%-4s) is in %d pure aliphatic ring(s) made of sp3 carbons (AR5)\n",
                            (i + lstart_at), new_names[i],
                            ar_set[i][5]);
          }              /* end for(i = 0; i < num_atoms; i++) */
     }

     fprintf(stdlog,
        "\n===============================================================\n");
     fprintf(stdlog,
        "---------------------electronic property-----------------------\n");
     for(i = 0; i < num_atoms; i++) {
          if(atom_ewd[i])
               fprintf(stdlog,
                       "atom[%5d] (%-4s) is an electron withdrawing atom\n",
                       (i + lstart_at), new_names[i]);
          else
               fprintf(stdlog,
                       "atom[%5d] (%-4s) is not an electron withdrawing atom\n",
                       (i + lstart_at), new_names[i]);
     }                   /* end for(i = 0; i < num_atoms; i++) */

     fprintf(stdlog,
       "\n===============================================================\n");
     fprintf(stdlog,
        "--------------------connectivity property----------------------\n");
     for(i = 0; i < num_atoms; i++) {
          for(j = 0; j <= 5; j++) {
               if(connect[i][j] != -1)
                    fprintf(stdlog,
                            "atom[%5d] (%-4s) %5d %5s\n", (i + lstart_at),
                            new_names[i], (connect[i][j] + lstart_at),
                            new_names[(connect[i][j])]);
          }              /* end for(j = 0; j <= 5; j++) */
     }                   /* end for(i = 0; i < num_atoms; i++) */

     fprintf(stdlog,
        "------------------------------end------------------------------\n\n");

     fprintf(stdlog, "\nFound %d angles, %d torsions", angle_max, tors_max);
     if(improper_num) fprintf(stdlog, ", and %d impropers\n", improper_num);
     else fprintf(stdlog, "\n");

/* If not united atom force field and only want the coordinates file,
   write and exit now.
   For united atoms force field, need to assign atom types and revise
   atom numbering first.
*/
     if(gro_only && !do_ua) {
          write_gro(gro_name, num_atoms, new_names, residue);
          return(0);
     }

     no_constraints = 0;
     switch(use_amber) {
          case 3:              /* amber force fields */
          case 0:
          case 1:  do_an_amber(ante_name, use_amber, leap_name, num_atoms,
                               num_bonds, tp_amber, angle_max, tors_max);
                   amber_flat_rings(tors_max);
                   break;
                        /* end amber force fields */

                        /* Tripos FF.  Already have correct atom and bond types */
          case 2:  for(i = 0; i < num_atoms; i++)
                        strcpy(a_ff_type[i], atom_type[i]);
                   FFbond = 10;
                   FFangle = 10;
                   FFtors = 10;
                   FFimpro = 10;
                   FFvdw = 10;

                   read_tripos(ante_name, 0, &FFbond, &FFangle, &FFtors, &FFimpro,
                               &FFvdw);

                   fprintf(stdlog,
                           "\n\tFound\n\t\t%d bond records\n\t\t%d angle records\n",
                           FFbond, FFangle);
                   fprintf(stdlog,
                           "\t\t%d torsion records\n\t\t%d impropers records\n",
                           FFtors, FFimpro);
                   fprintf(stdlog, "\t\t%d van der Waals records\n", FFvdw);
                   read_tripos(ante_name, 1, &FFbond, &FFangle, &FFtors, &FFimpro,
                               &FFvdw);

                   set_amber_vdwFF(num_atoms, FFvdw);     /* works the same for Tripos FF */
                   set_tripos_bondFF(num_bonds, FFbond);
                   set_tripos_angleFF(angle_max, FFangle);
                   set_tripos_torsFF(tors_max, FFtors);
                   set_tripos_improFF(improper_num, FFimpro);
                   free_triposFF(FFbond, FFangle, FFtors, FFimpro, FFvdw);
                   if(use_charge)
                        adjust_ion_types(num_atoms, a_ff_type);
                   amber_flat_rings(tors_max);
                   break;
                        /* end Tripos FF */

/* Gromacs FF.  Needs to use two levels of atom type definition tables
   because some pairs of atom types use different force constants and
   base values when in different structures.  These really should have
   been made distinct atom types, but the definers of these force fields
   did not do that.
   Measured lengths need external correction before writing the topology.
   It is imperative that the standard gromacs non-bonded parameters file
   be the included non-bonded parameters file in the topology.
*/
          case 4:  do_gmx_FFs(ante_name, "/ATOMTYPE_GMX1.DEF",
                              num_atoms, num_bonds, amberff,
                              angle_max, tors_max, do_ua,
                              residue, by_correl);

                   if(!parm_adj) {             /* adjust measured bond lengths */
                        for(i = 0; i < num_bonds; i++)
                             meas_len[i] = meas_len[i] * 0.1;
                   }
                   correct_dihed_types(tors_max);
                   strcpy(mycorr_name, "ffG");
                   strcat(mycorr_name, &amberff[3]);
                   strcat(mycorr_name, "nb.itp");
                   whichitp = &mycorr_name[0];
                   break;
                        /* end Gromacs FF */

/* OPLS all atoms FF.  Needs to use two levels of atom type definition tables
   because force field parameters are defined according to a modification of
   the Kollman atom types, but the nonbonded parameters are defined based on
   special numeric atom types just for OPLS.
   Measured lengths need external correction before writing the topology.
   It is imperative that the standard OPLS non-bonded parameters file be
   the included non-bonded parameters file in the topology.
*/
          case 5:  do_gmx_FFs(ante_name, "/ATOMTYPE_OPLSAA1.DEF",
                              num_atoms, num_bonds, amberff,
                              angle_max, tors_max, do_ua,
                              residue, by_correl);

                   if(!parm_adj) {             /* adjust measured bond lengths */
                        for(i = 0; i < num_bonds; i++)
                             meas_len[i] = meas_len[i] * 0.1;
                   }
                   correct_dihed_types(tors_max);
                   strcpy(mycorr_name, "ffoplsaanb.itp");
                   whichitp = &mycorr_name[0];
                   break;
                        /* end OPLS all atoms FF */

          default: sprintf(mess,
                           "Force field %s not yet implemented\n",
                           amberff);
                   my_fatal(FARGS, mess);
                   exit(1);
     }                                 /* end switch(use_amber) */

/* If only want the coordinates file for a united atoms force field,
   can write and exit now, finally */
     if(gro_only) {
          write_gro(gro_name, num_atoms, new_names, residue);
          return(0);
     }

     set_wanted_tors(tors_max, purge_level);

/* fix measured torsions for multiplicities as needed and able */
     fix_tor(tors_max, improper_num);
     count_torsions(tors_max);

/* write log of results */
     fprintf(stdlog, "\n\n\tAtom Force Field Results\n");
     fprintf(stdlog,
             "     no.    name   type   pol.   radius     pot.\n");
     for(i = 0; i < num_atoms; i++) {
          fprintf(stdlog, "%7d %8s", (i + lstart_at),
                  new_names[i]);
          if(a_ff_type[i] == NULLCHAR) fprintf(stdlog, "   ***");
          else fprintf(stdlog, " %6s", a_ff_type[i]);
          if(a_pol[i] > 9999.0) fprintf(stdlog, " *******");
          else fprintf(stdlog, " %7.3f", a_pol[i]);
          if(a_vdw_atom_radius[i] > 9999.0) fprintf(stdlog,
                                                    " ******** ********\n");
          else fprintf(stdlog, " %8.4f %8.4f\n", a_vdw_atom_radius[i],
                       a_vdw_atom_pot[i]);
     }                                 /* end for(i = 0; i < num_atoms; i++) */

     fprintf(stdlog, "\n\tBond Force Field Results\n");
     fprintf(stdlog,
             "  no.     bond    type  interp.  force  length   measured\n");
     for(i = 0; i < num_bonds; i++) {
          fprintf(stdlog, "%4d %5s-%5s   %2s   %3d",
                  (i + 1), new_names[(bond_i[i])], new_names[(bond_j[i])],
                  bond_type[i], how_bonded[i]);
          switch(a_bond_type[i]) {
               case 0:  fprintf(stdlog, "     ****** *******");
                        break;

               case 5:  fprintf(stdlog, "                   ");
                        break;

               case 7:
               case 6:
               case 1:
               case 2:  fprintf(stdlog, "     %6.1f %7.3f",
                                a_bond_force[i], a_bond_length[i]);
                        break;

               case 3:
               case 4:  fprintf(stdlog, "     %6.1f %6.1f %7.3f",
                                a_bond_force[i], a_bond_beta[i],
                                a_bond_length[i]);
                        break;

               case 8:
               case 9:  fprintf(stdlog, "     %6.1f %7d",
                                a_bond_force[i], (int)a_bond_length[i]);

               default: fprintf(stdlog, "Bond type error, bond number %d.  Length", i);
                        break;
          }
          fprintf(stdlog, "    %7.3f\n", meas_len[i]);
     }                                 /* end for(i = 0; i < num_bonds; i++) */

     if(no_constraints > 0) {
          iprtit = 0;
          for(i = 0; i < angle_max; i++) {
               if((constrain_i[i] < 0) ||
                  (constrain_j[i] < 0) ||
                  (constr_type[i] < 1) ||
                  (constr_type[i] > 2))
                         continue;
               if(!iprtit) {
                    fprintf(stdlog, "\n\tConstraints Force Field Results\n");
                    fprintf(stdlog,
                            "  no.     Atoms    type  distance\n");
                    iprtit = 1;
               }
               fprintf(stdlog,"%4d %6s %6s   %1d %10.5f\n", (i + 1),
                       new_names[(constrain_i[i])], new_names[(constrain_j[i])],
                       constr_type[i], dist_constr[i]);
          }
     }
     fprintf(stdlog, "\n\tAngles Force Field Results\n");
     fprintf(stdlog, "Angle        Atoms        force    angle     method measured\n");
     for(i = 0; i < angle_max; i++) {
          fprintf(stdlog,"%4d %5s-%5s-%5s", (i + 1), new_names[(angle_tab[i][0])],
                          new_names[(angle_tab[i][1])],
                          new_names[(angle_tab[i][2])]);
          switch(a_angle_type[i]) {
               case -1: fprintf(stdlog,
                                " Deleted planar opposition angle %d.  Theta",
                                (i + 1));
                        break;

               case 0:  fprintf(stdlog, "   ******   ******   ******");
                        break;

               case 1:
               case 2:  fprintf(stdlog, "   %6.1f   %6.2f   %6d",
                                a_angle_force[i], a_angle_angle[i], a_meth[i]);
                        break;

               case 3:  for(j = 0; j < 3; j++)
                             fprintf(stdlog, "   %10.5f", a_angle_quartic[j][i]);
                        fprintf(stdlog, "  %6d", a_meth[i]);
                        break;

               case 4:  for(j = 0; j < 4; j++)
                             fprintf(stdlog, "   %10.5f", a_angle_quartic[j][i]);
                        fprintf(stdlog, "  %6d", a_meth[i]);
                        break;

               case 5:  for(j = 0; j < 3; j++)
                             fprintf(stdlog, "   %10.5f", a_angle_quartic[j][i]);
                        fprintf(stdlog, "   %6.2f   %6d",
                                a_angle_angle[i], a_meth[i]);
                        break;

               case 6:  for(j = 0; j < 5; j++)
                             fprintf(stdlog, "   %10.5f", a_angle_quartic[j][i]);
                        fprintf(stdlog, "   %6.2f   %6d",
                                a_angle_angle[i], a_meth[i]);
                        break;

               case 8:  fprintf(stdlog, "%6.1f  %6d   %6d",
                                a_angle_quartic[1][i], (int)a_angle_quartic[0][i],
                                 a_meth[i]);
                        break;

               default: fprintf(stdlog,
                                "Angle type error, angle number %d.  Theta",
                                (i + 1));
                        break;
          }
          fprintf(stdlog, "    %7.3f\n", meas_ang[i]);
     }                                 /* end for(i = 0; i < angle_max; i++) */

     if((tors_max > 0) && (tors_max != tors_imp)) {
          fprintf(stdlog, "\n\tTorsions Force Field Results\n");
          if(regs_all)
               regs_log(tors_max);
          if(rbs_all)
               rbs_log(tors_max);
          if(fourier_all)
               fourier_log(tors_max);
          if(table_all)
               table_log(tors_max);
     }                                 /* end if((tors_max > 0) && (tors_max != tors_imp)) */

     if(tors_imp) {
          fprintf(stdlog, "\n\tPlanar Ring Improper Force Field Results\n");
          fprintf(stdlog, "  no.       Improper Locus       force  phase   measured\n");
          for(i = 0; i < tors_max; i++) {
               if(a_tors_type[i] != 2)
                    continue;
               fprintf(stdlog, "%4d %5s-%5s-%5s-%5s", (i + 1), new_names[(torsion_tab[i][0])],
                       new_names[(torsion_tab[i][1])], new_names[(torsion_tab[i][2])],
                       new_names[(torsion_tab[i][3])]);
               fprintf(stdlog, "   %6.1f %6.1f   %7.3f\n", a_tors_force[i],
                       a_tors_phase[i], meas_tors[i]);
         }
     }

     if(improper_num) {
          fprintf(stdlog, "\n\tImproper Force Field Results\n");
          fprintf(stdlog, "  no.       Improper Locus       force  phase   measured\n");
          for(i = 0; i < improper_num; i++) {
               fprintf(stdlog, "%4d %5s-%5s-%5s-%5s", (i + 1),
                       new_names[(improperid[i][0])], new_names[(improperid[i][1])],
                       new_names[(improperid[i][2])], new_names[(improperid[i][3])]);
               if(a_impr_mult[i] > 999) {
                    fprintf(stdlog, "   ****** ******");
               }
               else
                    fprintf(stdlog, "   %6.1f %6.1f",
                            a_impr_force[i], a_impr_phase[i]);
               fprintf(stdlog, "   %7.3f", meas_impro[i]);
               if(a_impr_mult[i] == 101)
                    fprintf(stdlog, "   default values\n");
               else  
                    fprintf(stdlog, "\n");
          }                            /* end for(i = 0; i < improper_num; i++) */
     }                                 /* end if(improper_num) */

/* After final arrays are fixed, set up and write the .gro file */

     write_gro(gro_name, num_atoms, new_names, residue);

/* Then set up and write the .itp file */

     write_top(itp_name, residue, num_atoms, num_bonds, tors_max, angle_max,
               tors_max, improper_num, a_ff_type, new_names, use_mine,
               parm_adj);

     write_ff(dnext, a_ff_type, num_atoms, a_vdw_atom_radius,
              a_vdw_atom_pot, parm_adj, whichitp);

     return(0);
}

/* print the usage message
   parameter is the name of the program
*/
void usage(char *name)
{
     version(name);
     printf("Read a Tripos .mol2 file and generate something approximating a gromacs\n");
     printf(".gro and .itp from it.\n");
     printf("\tUses the following parameters and flags:\n");
     printf("\t\t-n                Common name for the .mol2, .gro, .itp, .top, and\n");
     printf("\t\t                  .log files (the path to and name of the input\n");
     printf("\t\t                  file stripped of the .mol2).  A revised .mol2\n");
     printf("\t\t                  file will be output to common nameMOL.mol2 as well.\n");
     printf("\t\t-dir              Absolute path to a directory that contains\n");
     printf("\t\t                  appropriate force field data.\n");
     printf("\t\t                  For amber type force fields, including glycam,\n");
     printf("\t\t                  this would be a path to a directory that contains\n");
     printf("\t\t                  directories ./dat/antechamber and ./dat/leap/parm\n");
     printf("\t\t                  (these in turn must contain appropriate antechamber\n");
     printf("\t\t                  and leap data files)\n");
     printf("\t\t                  For tripos force field this would be the value of\n");
     printf("\t\t                  Tripos variable $TA_ASCTABLES\n");
     printf("\t\t                  For gromacs force fields including oplsaa, this\n");
     printf("\t\t                  would be a path to a directory that contains the\n");
     printf("\t\t                  files ATOMTYPE_GMX1.DEF, ATOMTYPE_gmx####.DEF,\n");
     printf("\t\t                  ffgmx####.dat (where #### designates the particular\n");
     printf("\t\t                  gromacs force field), ATOMTYPE_OPLSAA1.DEF,\n");
     printf("\t\t                  ATOMTYPE_oplsaa.DEF, and ffoplsaa.dat\n");
     printf("\t\t-ff               Force field to use, either tripos, gaff, glycam## or\n");
     printf("\t\t                  amber## (where ## designates the year of the glycam\n");
     printf("\t\t                  or amber file (the number between parm or the glycam\n");
     printf("\t\t                  and the .dat in /dat/leap/parm)), gmx#### (where\n");
     printf("\t\t                  #### designates the particular gromacs force field,\n");
     printf("\t\t                  or oplsaa.\n");
     printf("\t\t-r                Optional residue name to be used for all atoms in the\n");
     printf("\t\t                  output.\n");
     printf("\t\t-last             Optional last number of a topology to which this\n");
     printf("\t\t                  topology is to be appended.  The first number for an\n");
     printf("\t\t                  atom in this topology will be one more than the\n");
     printf("\t\t                  number given.  A -resnum is required along with this\n");
     printf("\t\t                  option.\n");
     printf("\t\t-resnum           Optional last residue number of a topology to which\n");
     printf("\t\t                  this topology is to be appended.  The first number\n");
     printf("\t\t                  for a residue in this topology will be one more than\n");
     printf("\t\t                  the number given.  A -last is required along with\n");
     printf("\t\t                  this option.\n");
     printf("\t\t-renumall         Optional flag to renumber topology as well as\n");
     printf("\t\t                  coordinates to permit appending topology elements to\n");
     printf("\t\t                  a main topology to permit addition of inter-molecular\n");
     printf("\t\t                  restraints\n");
     printf("\t\t-meas             Optional flag to use measured bond lengths, angles\n");
     printf("\t\t                  and torsions.  Otherwise, use the values found in the\n");
     printf("\t\t                  selected forcefield.  Note that neither for oplsaa,\n");
     printf("\t\t                  nor for standard Ryckaert-Bellemans dihedrals is\n");
     printf("\t\t                  recalculation of the torsion parameters to the measured\n");
     printf("\t\t                  values done.\n");
     printf("\t\t-purge [level]    Optional flag to decrease the number of dihedral\n");
     printf("\t\t                  angles included in the topology.  A single digit\n");
     printf("\t\t                  between 0 and 5 can be specified to indicate the\n");
     printf("\t\t                  extent of the decrease with 0 as unpurged and 5 as\n");
     printf("\t\t                  most stringent purge.  Default value is 1.\n");
     printf("\t\t-gro              Only write out gromacs coordinates (.gro) file.  The\n");
     printf("\t\t                  -r, -last, -resnum, -renumall, and -rename options\n");
     printf("\t\t                  will act upon the coordinate file appropriately.  The\n");
     printf("\t\t                  -meas option is inappliacble to gromacs coordinate\n");
     printf("\t\t                  files and will be ignored if specified.\n");
     printf("\t\t-move             Optional flag to translate the molecule so that it is\n");
     printf("\t\t                  centered on its center of mass.  This option conflicts\n");
     printf("\t\t                  with -renumall because it is assumed that -renumall\n");
     printf("\t\t                  is done for a molecule aligned with another to which\n");
     printf("\t\t                  this topology is to be added.\n");
     printf("\t\t-rename           Optional flag to rename atoms for the output.\n");
     printf("\t\t-charge           Optional flag to assign charges from the force field\n");
     printf("\t\t                  atom types to replace the .mol2 charges.  Currently\n");
     printf("\t\t                  only works with oplsaa.\n\n");
     printf("\t\t-h                Prints program version information and this usage\n");
     printf("\t\t                  message.\n");
     printf("\t\t-H                Same as -h.\n");
     printf("\t\t-v                Prints program version information\n");
     printf("\t\t-V                Same as -v.\n");
     return;
}

/* print the version message
   parameter is the name of the program
*/
void version(char *name)
{
     printf("\n\t\t\t%s\n", name);
     printf("\t\t\t%s\n", vers_no);
     printf("\t\t\tWritten by Bruce D. Ray, Ph.D.\n");
     printf("\t\t\tIUPUI Physics Department\n");
     printf("\t\t\t402 N. Blackford St.\n");
     printf("\t\t\tIndianapolis, Indiana  46202\n\t\t\tUSA\n\n");
     return;
}

/* Handle amber type force field conversions.
   parameters are:
        ante_name                  directory of antechamber atom type data files
        use_amber                  flag for amber versus gaff force field
        leap_name                  directory of leap amber and gaff parm files
        num_atoms                  number of atoms in molecule
        num_bonds                  number of bonds in molecule
        tp_amber                   year number of amber force field to use
        angle_max                  number of angles in molecule
        tors_max                   number of torsions in molecule 
*/
void do_an_amber(char *ante_name, int use_amber, char *leap_name, int num_atoms,
                 int num_bonds, char *tp_amber, int angle_max, int tors_max)
{
     int num, num2, count1, count2, i, j;
     int FFmol, FFbond, FFangle, FFtors, FFimpro, FFvdw, FFequiv, FFeqlen;

     strcpy(ambername, ante_name);
     switch(use_amber) {                       /* lets me add forcefields later */
          case 0:  strcat(ambername, "ATOMTYPE_GFF.DEF");
                   break;

          case 1:  strcat(ambername, "ATOMTYPE_AMBER.DEF");
                   break;

          case 3:  strcat(ambername, "ATOMTYPE_GLYCAM.DEF");
                   break;

          default: my_fatal(FARGS, "System error in switch statements\n");
                   exit(100);
                   break;
     }

     read_at_types(&num, &num2, (4 * MAXCHAR), MAXCHAR, 0,
                   ambername);                         /* first count lines */

     fprintf(stdlog, "\nIn %s there are %d WILDATOM lines and %d ATD lines\n",
                     ambername, num2, num);
     count1 = num;
     count2 = num2;
     read_at_types(&num, &num2, count1, count2, 1,
                   ambername);                           /* now read them */

     init_atom_types(num_atoms, num_bonds);
     bond_info(num_atoms, num_bonds);

     if(!use_amber) {
          judge_atom_types(count1, count2, num_atoms, num_bonds, a_ff_type);
          adjustments(num_atoms, num_bonds, a_ff_type);
          adjust_type_cp(num_atoms, num_bonds, a_ff_type);
          check_type_errors(num_atoms, num_bonds, a_ff_type);
     }
     else
          judge_atom_types(count1, count2, num_atoms, num_bonds, a_ff_type);

     free_atom_types(num_atoms);
     release_types(count1, count2);

     if(!use_amber) {
          similars_alloc(num_atoms);          /* set up the similar names arrays */

          strcpy(mycorr_name, ante_name);
          strcat(mycorr_name, "ATCOR.DAT");
          read_at_corr(&num, MAXCHAR, 0, mycorr_name);     /* count them */

          fprintf(stdlog,
                  "Number of correlation entries in file = %d\n", num);

          num2 = num;
          read_at_corr(&num, num2, 1, mycorr_name);        /* now read them */

          set_similars(num_atoms, num2);                 /* use them */

          release_corr(num2);                            /* dump them */
     }

     strcpy(ambername, leap_name);                  /* select parameters file */
     switch(use_amber) {
          case 0:  strcat(ambername, "gaff.dat");
                   break;

          case 1:  strcat(ambername, "parm");
                   strcat(ambername, tp_amber);
                   strcat(ambername, ".dat");
                   break;

          case 3:  strcat(ambername, "glycam");
                   strcat(ambername, tp_amber);
                   strcat(ambername, ".dat");
                   break;

          default: my_fatal(FARGS,
                            "System error in switch statements\n");
                   exit(100);
                   break;
     }

     FFmol = 10;
     FFbond = 10;
     FFangle = 10;
     FFtors = 10;
     FFimpro = 10;
     FFvdw = 10;
     FFequiv = 10;
     FFeqlen = 10;

     read_amber(ambername, 0, &FFmol, &FFbond, &FFangle, &FFtors, &FFimpro,
                &FFvdw, &FFequiv, &FFeqlen);

     fprintf(stdlog,
   "\n\tFound\n\t\t%d atom records\n\t\t%d bond records\n\t\t%d angle records\n",
            FFmol, FFbond, FFangle);
     fprintf(stdlog,
   "\t\t%d torsion records\n\t\t%d impropers records\n\t\t%d van der Waals records\n",
            FFtors, FFimpro, FFvdw);
     fprintf(stdlog,
   "\t\t%d van der Waals equivalence records with %d maximum record length\n",
            FFequiv, FFeqlen);

     FFvdw = FFvdw + (FFequiv * (FFeqlen - 1));
     read_amber(ambername, 1, &FFmol, &FFbond, &FFangle, &FFtors, &FFimpro,
                &FFvdw, &FFequiv, &FFeqlen);

     free_equivs(FFequiv, FFeqlen);

     if(!use_amber) {
          set_gaff_atomFF(num_atoms, FFmol);
          set_gaff_vdwFF(num_atoms, FFvdw);
          set_gaff_bondFF(num_bonds, FFbond);
          set_gaff_angleFF(angle_max, FFangle);
          set_gaff_torsFF(tors_max, FFtors);
          if(improper_num)
               set_gaff_improFF(improper_num, FFimpro);

          for(i = 0; i < num_atoms; i++)             /* make gaff atom types upper
                                                        case for topology output */
               for(j = 0; j < strlen(a_ff_type[i]); j++)
                    a_ff_type[i][j] = toupper(a_ff_type[i][j]);

          release_similars(num_atoms);
     }
     else {
          set_amber_atomFF(num_atoms, FFmol);
          set_amber_vdwFF(num_atoms, FFvdw);
          set_amber_bondFF(num_bonds, FFbond);
          set_amber_angleFF(angle_max, FFangle);
          set_amber_torsFF(tors_max, FFtors);
          if(improper_num)
               set_amber_improFF(improper_num, FFimpro);
     }

     free_amberFF(FFmol, FFbond, FFangle, FFtors, FFimpro, FFvdw);
     if(use_charge)
          adjust_ion_types(num_atoms, a_ff_type);

     return;
}

/* Send the regular dihedral angles to the log.
   parameter is
        i_count                        total number of torsions
*/
void regs_log(int i_count)
{
     int i, j, k;

     fprintf(stdlog, "Simple Torsions\n");
     fprintf(stdlog,
   "Tors           Atoms           mul force   phase   term   measured   purge\n");

     k = regs_all;
     for(i = 0; i < i_count; i++) {
          if((a_tors_type[i] > 1) &&
             (a_tors_type[i] != 9))
                    continue;
          fprintf(stdlog,"%4d %5s-%5s-%5s-%5s", (i + 1),
                          new_names[(torsion_tab[i][0])],
                          new_names[(torsion_tab[i][1])],
                          new_names[(torsion_tab[i][2])],
                          new_names[(torsion_tab[i][3])]);
          switch(a_tors_type[i]) {
               case 0:  fprintf(stdlog,"  ***  ****** ****** ******");
                        break;

               case 9:
               case 1:  fprintf(stdlog, "  %3d %6.2f  %6.1f %6.1f",
                                a_tors_mult[i], a_tors_force[i],
                                a_tors_phase[i], a_tors_term[i]);
               default: break;
          }

          fprintf(stdlog, "   %7.3f", meas_tors[i]);
          if(wanted_tors != NULL) {
               if(wanted_tors[i])
                    fprintf(stdlog, "   %2d", (6 - wanted_tors[i]));
               else
                    fprintf(stdlog, "    X");
          }
          fprintf(stdlog, "\n");

          if(multors_cnt[i])
               for(j = 0; j < multors_cnt[i]; j++) {
                    fprintf(stdlog,"%4d %5s-%5s-%5s-%5s", (i + 1),
                                    new_names[(torsion_tab[i][0])],
                                    new_names[(torsion_tab[i][1])],
                                    new_names[(torsion_tab[i][2])],
                                    new_names[(torsion_tab[i][3])]);
                    switch(a_tors_type[i]) {
                         case 0:  fprintf(stdlog,
                                    "  ***  ****** ****** ******\n");
                                  break;

                         case 9:
                         case 1:  fprintf(stdlog,
                                    "  %3d %6.2f  %6.1f %6.1f\n",
                                          (int) multiples[i][j][0],
                                          multiples[i][j][1],
                                          multiples[i][j][2],
                                          multiples[i][j][3]);
                                  break;
                    }          /* end switch(a_tors_type[i]) */
               }               /* end for(j = 0; j < multors_cnt[i]; j++) */
          k--;
          if(!k) break;
     }                         /* end for(i = 0; i < i_count; i++) */

     return;
}

/* Send RB dihedrals to the log
   parameter is
        i_count                        total number of torsions
*/
void rbs_log(int i_count)
{
     int i, j, k;

     fprintf(stdlog, "RB Dihedrals\n");
     fprintf(stdlog,
   "Tors           Atoms           c0          c1          c2          c3");
     fprintf(stdlog,
        "          c4          c5     measured   purge\n");

     k = rbs_all;
     for(i = 0; i < i_count; i++) {
          if(a_tors_type[i] != 3)
               continue;
          fprintf(stdlog,"%4d %5s-%5s-%5s-%5s", (i + 1),
                  new_names[(torsion_tab[i][0])],
                  new_names[(torsion_tab[i][1])],
                  new_names[(torsion_tab[i][2])],
                  new_names[(torsion_tab[i][3])]);
          for(j = 0; j < 6; j++)           /* print RB coefficients */
               fprintf(stdlog, "  %10.5f", a_tors_set[j][i]);
          fprintf(stdlog, "   %7.3f", meas_tors[i]);
          if(wanted_tors != NULL) {
               if(wanted_tors[i])
                    fprintf(stdlog, "   %2d", (6 - wanted_tors[i]));
               else
                    fprintf(stdlog, "    X");
          }
          fprintf(stdlog, "\n");
          k--;
          if(!k) break;
     }                            /* end for(i = 0; i < i_count; i++) */

     return;
}

/* Send fourier dihedrals to the log
   parameter is
        i_count                        total number of torsions
*/
void fourier_log(int i_count)
{
     int i, j, k;

     fprintf(stdlog, "Fourier Dihedrals\n");
     fprintf(stdlog,
    "Tors           Atoms           c0          c1          c2          c3");
     fprintf(stdlog, "     measured   purge\n");

     k = fourier_all;
     for(i = 0; i < i_count; i++) {
          if(a_tors_type[i] != 5)
               continue;
          fprintf(stdlog,"%4d %5s-%5s-%5s-%5s", (i + 1),
                  new_names[(torsion_tab[i][0])],
                  new_names[(torsion_tab[i][1])],
                  new_names[(torsion_tab[i][2])],
                  new_names[(torsion_tab[i][3])]);
          for(j = 0; j < 4; j++)           /* print fourier coefficients */
               fprintf(stdlog, "  %10.5f", a_tors_set[j][i]);
          fprintf(stdlog, "   %7.3f", meas_tors[i]);
          if(wanted_tors != NULL) {
               if(wanted_tors[i])
                    fprintf(stdlog, "   %2d", (6 - wanted_tors[i]));
               else
                    fprintf(stdlog, "    X");
          }
          fprintf(stdlog, "\n");
          k--;
          if(!k) break;
     }                           /* end for(i = 0; i < i_count; i++) */

     return;
}

/* Send tabular dihedrals to the log
   parameter is
        i_count                        total number of torsions
*/
void table_log(int i_count)
{
     int i, j, k;

     fprintf(stdlog, "Tabular Dihedrals\n");
     fprintf(stdlog,
             "Tors           Atoms           num          num");
     fprintf(stdlog, "     measured   purge\n");

     k = table_all;
     for(i = 0; i < i_count; i++) {
          if(a_tors_type[i] != 8)
               continue;
          fprintf(stdlog,"%4d %5s-%5s-%5s-%5s", (i + 1),
                  new_names[(torsion_tab[i][0])],
                  new_names[(torsion_tab[i][1])],
                  new_names[(torsion_tab[i][2])],
                  new_names[(torsion_tab[i][3])]);
          fprintf(stdlog, "  %3d %6.2f",
                  (int)a_tors_set[0][i], a_tors_set[1][i]);
          fprintf(stdlog, "   %7.3f", meas_tors[i]);
          if(wanted_tors != NULL) {
               if(wanted_tors[i])
                    fprintf(stdlog, "   %2d", (6 - wanted_tors[i]));
               else
                    fprintf(stdlog, "    X");
          }
          fprintf(stdlog, "\n");
          k--;
          if(!k) break;
     }                            /* end for(i = 0; i < i_count; i++) */

     return;
}
