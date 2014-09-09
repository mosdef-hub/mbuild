/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to read and write *.mol2 (Tripos Sybyl) files
        Function to assign atomic number and mass based on symbol given
        caseless strncmp function

   NOTICE: Portions of this are derivative work based on study of:
      1. The routines found GROMACS version 3.3.1 written by David
         van der Spoel, Erik Lindahl, Berk Hess, and others and
         copyright (c) 1991-2000, University of Groningen, The Netherlands.
         Copyright (c) 2001-2004, The GROMACS development team

      2. The routines found in antechamber 1.27 file mol2.c

      3. The routines found in Dock 6.1 file dockmol.cpp
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "utilities.h"
#include "block_memory.h"
#include "readmol2.h"
#include "mol2.h"
#include "tripos.h"

#define MAXNAME 128
#define MAXNAMEM1 127

/* table for atomic masses */

typedef struct {
     float mass;
     int len;
     int at_no;
     char type[3];
} ATOMIC;

static ATOMIC atomicno[] = {
     {   1.00800,  1,   1, "H"  },
     {   4.00260,  2,   2, "He" },
     {   6.94100,  2,   3, "Li" },
     {   9.01218,  2,   4, "Be" },
     {  10.81000,  1,   5, "B"  },
     {  12.01100,  1,   6, "C"  },
     {  14.00670,  1,   7, "N"  },
     {  15.99940,  1,   8, "O"  },
     {  18.99840,  1,   9, "F"  },
     {  20.17970,  2,  10, "Ne" },
     {  22.98977,  2,  11, "Na" },
     {  24.30500,  2,  12, "Mg" },
     {  26.98150,  2,  13, "Al" },
     {  28.08000,  2,  14, "Si" },
     {  30.97376,  1,  15, "P"  },
     {  32.06000,  1,  16, "S"  },
     {  35.45300,  2,  17, "Cl" },
     {  39.94800,  2,  18, "Ar" },
     {  39.09830,  1,  19, "K"  },
     {  40.08000,  2,  20, "Ca" },
     {  44.95590,  2,  21, "Sc" },
     {  47.90000,  2,  22, "Ti" },
     {  50.94140,  1,  23, "V"  },
     {  51.99600,  2,  24, "Cr" },
     {  54.93800,  2,  25, "Mn" },
     {  55.84700,  2,  26, "Fe" },
     {  58.93320,  2,  27, "Co" },
     {  58.71000,  2,  28, "Ni" },
     {  63.54600,  2,  29, "Cu" },
     {  65.37000,  2,  30, "Zn" },
     {  69.72000,  2,  31, "Ga" },
     {  72.59000,  2,  32, "Ge" },
     {  74.92160,  2,  33, "As" },
     {  78.96000,  2,  34, "Se" },
     {  79.90400,  2,  35, "Br" },
     {  83.79800,  2,  36, "Kr" },
     {  85.46780,  2,  37, "Rb" },
     {  87.62000,  2,  38, "Sr" },
     {  88.90590,  1,  39, "Y"  },
     {  91.22000,  2,  40, "Zr" },
     {  92.90640,  2,  41, "Nb" },
     {  95.94000,  2,  42, "Mo" },
     {  98.90620,  2,  43, "Tc" },
     { 101.07000,  2,  44, "Ru" },
     { 102.90550,  2,  45, "Rh" },
     { 106.40000,  2,  46, "Pd" },
     { 107.86800,  2,  47, "Ag" },
     { 112.40000,  2,  48, "Cd" },
     { 114.82000,  2,  49, "In" },
     { 118.69000,  2,  50, "Sn" },
     { 121.75000,  2,  51, "Sb" },
     { 127.60000,  2,  52, "Te" },
     { 126.90450,  1,  53, "I"  },
     { 131.29300,  2,  54, "Xe" },
     { 132.90540,  2,  55, "Cs" },
     { 137.33000,  2,  56, "Ba" },
     { 138.90550,  2,  57, "La" },
     { 140.12000,  2,  58, "Ce" },
     { 140.90770,  2,  59, "Pr" },
     { 144.24000,  2,  60, "Nd" },
     { 145.00000,  2,  61, "Pm" },
     { 150.40000,  2,  62, "Sm" },
     { 151.96000,  2,  63, "Eu" },
     { 157.25000,  2,  64, "Gd" },
     { 158.92540,  2,  65, "Tb" },
     { 162.50000,  2,  66, "Dy" },
     { 164.93030,  2,  67, "Ho" },
     { 167.26000,  2,  68, "Er" },
     { 168.93420,  2,  69, "Tm" },
     { 173.04000,  2,  70, "Yb" },
     { 174.97000,  2,  71, "Lu" },
     { 178.49000,  2,  72, "Hf" },
     { 180.94790,  2,  73, "Ta" },
     { 183.85000,  1,  74, "W"  },
     { 186.20000,  2,  75, "Re" },
     { 190.20000,  2,  76, "Os" },
     { 192.22000,  2,  77, "Ir" },
     { 195.09000,  2,  78, "Pt" },
     { 196.96650,  2,  79, "Au" },
     { 200.59000,  2,  80, "Hg" },
     {   0.00000,  2,   0, "Lp" },
     {   0.00000,  2,   0, "Du" }
};

#define ntypes   82            /* number of atomic mass table entries */

/* end atomic mass information */

/* set atomic number
   parameters are:
        len                    length of atom type string
        type                   atom type to find
        atmass                 returned atomic mass
        atno                   returned atomic number
        atsymb                 returned atom symbol
*/
void set_atomic_no(int len, char *type, float *atmass, int *atno, char *atsymb)
{
     int k, indx;

     indx = 0;
     for(k = 0; k < ntypes; k++ ) {
          if(len != atomicno[k].len) continue;
          if(!(caseless_strncmp(type, atomicno[k].type, len))) {
               *atmass = atomicno[k].mass;
               *atno = atomicno[k].at_no;
               strncpy(atsymb, atomicno[k].type, len);
               atsymb[len] = NULLCHAR;
               indx++;
               break;
          }
     }          /* end for( k = 0; k < ntypes; k++ ) */

     if(!indx ) {
          sprintf(mess,
                  "Atomic symbol %s in atom type not found.\n",
                  type);
          my_fatal(FARGS, mess);
          exit(1);
     }

     return;
}

/* Allocate storage for the mol2 file to be read
   parameters are:
        read_atomnum           number of atoms
        read_bondnum           number of bonds
        feat                   number of features (used for number of restraints)
*/
void alloc_mol2(int read_atomnum, int read_bondnum, int feat)
{
     int i, j;

     atom_segno = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     atom_atno = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     bond_i = (int *)allocator(read_bondnum, sizeof(int), FARGS);
     bond_j = (int *)allocator(read_bondnum, sizeof(int), FARGS);
     how_bonded = (int *)allocator(read_bondnum, sizeof(int), FARGS); 
     bond_count = (int *)allocator(read_atomnum, sizeof(int), FARGS);

     for(i = 0; i < read_atomnum; i++)
          bond_count[i] = 0;

     atom_sat = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     atom_ewd = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     connect = (int **)allocator(read_atomnum, sizeof(int*), FARGS);

     for(i = 0; i < read_atomnum; i++) {
          connect[i] = (int *)allocator(8, sizeof(int), FARGS);
          for(j = 0; j < 8; j++)
               connect[i][j] = -1;
     }

     bonds_list = (int **)allocator(read_atomnum, sizeof(int*), FARGS);

     for(i = 0; i < read_atomnum; i++) {
          bonds_list[i] = (int *)allocator(8, sizeof(int), FARGS);
          for(j = 0; j < 8; j++)
               bonds_list[i][j] = -1;
     }

     atom_x = (float *)allocator(read_atomnum, sizeof(float), FARGS);
     atom_y = (float *)allocator(read_atomnum, sizeof(float), FARGS);
     atom_z = (float *)allocator(read_atomnum, sizeof(float), FARGS);
     atom_charge = (float *)allocator(read_atomnum, sizeof(float), FARGS);
     atom_mass = (float *)allocator(read_atomnum, sizeof(float), FARGS);
     atom_name = (char **)mat_alloc(read_atomnum, 12, FARGS);
     atom_symbl = (char **)mat_alloc(read_atomnum, 5, FARGS);
     atom_type = (char **)mat_alloc(read_atomnum, 12, FARGS);
     atom_residue = (char **)mat_alloc(read_atomnum, 12, FARGS);
     bond_type = (char **)mat_alloc(read_bondnum, 6, FARGS);
     a_ff_type = (char **)mat_alloc(read_atomnum, 9, FARGS);
     improper = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     atom_order = (int *)allocator(read_atomnum, sizeof(int), FARGS);
     back_order = (int *)allocator(read_atomnum, sizeof(int), FARGS);

     if(feat > 0) {
          restr = (int **)allocator(feat, sizeof(int *), FARGS);
          frest = (float **)allocator(feat, sizeof(float *), FARGS);
          for(i = 0; i < feat; i++) {
               restr[i] = (int *)allocator(5, sizeof(int), FARGS);
               frest[i] = (float *)allocator(3, sizeof(float), FARGS);
          }
     }

     rstr_dist = 0;
     rstr_range = 0;
     rstr_angle = 0;
     rstr_torsn = 0;

     return;
}

/* charge adjustment routine to give integer total charge
   parameter is:
        read_atomnum           number of atoms
*/
void adjust_chrg(int read_atomnum)
{
     double tot_charge, delta_chg, chg_intgr, max_chg;
     int i, where_max;

     tot_charge = 0.0;
     where_max = -1;
     max_chg = 0.0;

     for(i = 0; i <  read_atomnum; i++ ) {
          if(atom_charge[i] == floor(atom_charge[i]))
               continue;                /* skip integral charges */
          if(bond_count[i]) {           /* only adjust bonded atoms */
               tot_charge += (double)atom_charge[i];
               if(fabs(max_chg) < fabs((double)atom_charge[i])) {
                    max_chg = atom_charge[i];
                    where_max = i;
               }
          }
     }

     chg_intgr = 1.0;
     if(tot_charge < 0.0) chg_intgr = -1.0;
     delta_chg = tot_charge - (double)((int)(tot_charge + (0.5 * chg_intgr)));

     if(fabs(delta_chg) > 1e-6) {
          if((fabs(delta_chg) < 1e-2) && !(where_max < 0)) {
               atom_charge[where_max] += (float) (chg_intgr * delta_chg);
               return;                  /* best to only adjust maximum charge */
          }
          for(i = 0; i <  read_atomnum; i++ ) {   /* proportional adjustment of charges */
               if(atom_charge[i] == floor(atom_charge[i]))
                    continue;                /* skip integral charges */
               if(bond_count[i])             /* only adjust bonded atoms */
                    atom_charge[i] += (float)((chg_intgr * delta_chg *
                                      (double)atom_charge[i]) / tot_charge);
          }
     }

     return;
}

/* Read MOL2  format
   parameters are:
        filename               name of file to read
        atomnum                returned number of atoms read
        bondnum                returned number of bonds read
*/

void readmol2(const char *filename, int *atomnum, int *bondnum )
{
     int i, j, which, tmpint1, tmpint3, len, feat;
     int read_atomnum, read_bondnum, index, k;
     int allowed, charges, mandatory, has_charge;
     int max, bndcheck;

     char type[MAXCHAR], *here;
     char tmpchar4[MAXCHAR], tmpchar5[MAXCHAR];
     char tmpc1[MAXCHAR], tempt;

     FILE *fpin;

     if ((fpin = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     adj_chg = 1.0;

     which = 0;     /* switch variable for header records handling */

     while(fgets(line, MAXCHAR, fpin) != NULL) {    /* Handle header records */
          if((here = strchr(line, '\n')) != NULL)
               *here = NULLCHAR;
          switch(which) {
               case 0: if(!strncmp(line, "@<TRIPOS>MOLECULE", 17))
                            which = 1;
                       break;
/* this loops until we have found a molecule record which must be the first
   official mol2 format record.  Then go to case 1. */

               case 1: if(!strncmp(line, "@<TRIPOS>", 9)) {
                            my_fatal(FARGS,
              "Missing header information before atom records reading mol2 file.\n");
                            exit(1);
					   }
/* must have name information before the atoms records. */

                       strncpy(mol_name, line, MAXNAMEM1);  /* molecule name line */
                       mol_name[MAXNAMEM1] = NULLCHAR;      /*  with strncpy, must always ensure
                                                                the string is terminated */
                       which = 2;
                       break;

               case 2: if(!strncmp(line, "@<TRIPOS>", 9)) {
                            my_fatal(FARGS,
               "Missing header information before atom records reading mol2 file.\n");
                            exit(1);
                       }
/* must have size information before the atoms records. */

                       index = strspn(line, "0123456789 \t\n ");
                       if( index < strlen(line) ) {
                            my_fatal(FARGS,
               "Illegal characters for size information line in mol2 file.\n");
                            exit(1);
                       }

                       tmpint1 = 0;
                       feat = 0;
                       tmpint3 = 0;
                       j = sscanf(line, "%d%d%d%d%d", &read_atomnum,
                                  &read_bondnum, &tmpint1, &feat, &tmpint3);
                       if(j < 4) feat = 0;
/* Allocate arrays for atoms and bonds */
                       alloc_mol2(read_atomnum, read_bondnum, feat);

                       which = 3;
                       break;

               case 3: if(!strncmp(line, "@<TRIPOS>", 9)) {
                            my_fatal(FARGS,
	       "Missing header information before atom records reading mol2 file.\n");
                            exit(1);
                       }
/* must have model information before the atoms records. */

                       strncpy(model_type, line, MAXNAMEM1);  /* model type information */
                       model_type[MAXNAMEM1] = NULLCHAR;      /*  with strncpy, must always ensure
                                                                  the string is terminated */
                       which = 4;
                       break;


               case 4: if(!strncmp(line, "@<TRIPOS>", 9)) {
                            my_fatal(FARGS,
               "Missing header information before atom records reading mol2 file.\n");
                            exit(1);
                       }
/* must have charge information before the atoms records. */

                       strncpy(charge_type, line, MAXNAMEM1);  /* charge type information */
                       charge_type[MAXNAMEM1] = NULLCHAR;      /*  with strncpy, must always ensure
                                                                   the string is terminated */

                       if(!strncmp(charge_type, "INVALID", 7) ||
                          !strncmp(charge_type, "NO_CHARGE", 9))
                                 use_charge = 0;
                       else use_charge = 1;
                       if(!strncmp(charge_type, "AMBER", 5))
                            adj_chg = 2.0;
                       which = 5;
                       break;

               case 5: if(!strncmp(line, "@<TRIPOS>ATOM", 13))
                            which = 6;
                       break;
/* this loops until we have found an atom record which follows the molecule record,
   title, size, model type, charge type, and energy, and comments lines.  We ignore
   the energy and comments, but save title, size, model, and charges information.
   case 6 will break out of this loop to a new loop where we process atom records
   until we reach a bond record. */

               default: my_fatal(FARGS,
               "Major system error.  Reached default in switch in reading mol2 file.\n");
                        exit(1);       /* not strictly necessary, but we have major
                                          problems if reach here */

          }        /* end switch(which) */

          if (which == 6) break;     /* escape the while loop */

     }        /* end while(fgets(line, MAXCHAR, fpin) != NULL) */

     if (which != 6) {
          my_fatal(FARGS,
            "File read to end without finding molecule record reading mol2 file.\n");
          exit(1);
     }

     for(i = 0; i < read_atomnum; i++ ) {     /* handle atom records */
          if(fgets(line, MAXCHAR, fpin) == NULL) {
               my_fatal(FARGS,
                 "Premature termination of mol2 file while reading atoms.\n");
               exit(1);          /* premature termination of file */
          }

          if(!strncmp(line, "@<TRIPOS>", 9)) {
               my_fatal(FARGS,
                 "Premature end of atoms record entry in mol2 file.\n");
               exit(1);          /* not enough atoms */
          }

          sscanf(line, "%d%s%f%f%f%s%d%s%f", &tmpint1, tmpc1,
                        &atom_x[i], &atom_y[i], &atom_z[i], tmpchar4,
                        &atom_segno[i], tmpchar5, &atom_charge[i]);

          if(!use_charge) atom_charge[i] = 0.0;

          len = strcspn(tmpchar4, ".");

          if(len > 2) {
               sprintf(mess,
                  "Atomic symbol %s in atom type in mol2 file too long.\n",
                  type);
               my_fatal(FARGS, mess);
               exit(1);
          }

          strncpy(type, tmpchar4, len);      /* set type for mass determination */
          type[len] = NULLCHAR;              /* with strncpy, must always ensure
                                                the string is terminated */
/* set atom mass based on atom type */

          set_atomic_no(len, type, &atom_mass[i], &atom_atno[i], atom_symbl[i]);

/* save full type and residue name for atom */

          strncpy(&atom_type[i][0], tmpchar4, 10);
          atom_type[i][10] = NULLCHAR;      /* with strncpy, must always ensure
                                               the string is terminated */

          strncpy(&atom_residue[i][0], tmpchar5, 10);
          atom_residue[i][10] = NULLCHAR;     /* with strncpy, must always ensure
                                                 the string is terminated */

/* Need alternative name if the atom name column lacks good names.  Take
   alternative name from the atom type */

         if(tmpc1[0] == '*')  strncpy(&atom_name[i][0], tmpchar4, len);
         else  strncpy(&atom_name[i][0], tmpc1, 4);

         atom_name[i][4] = NULLCHAR;          /* with strncpy, must always ensure
                                                 the string is terminated */

/* Make sure the atom name does not begin with a digit */

         index = 0;
         k = strlen(&atom_name[i][0]);

         while(isdigit(atom_name[i][0])) {
              tempt = atom_name[i][0];
              index++;

              if(index > k ) {
                   my_fatal(FARGS,
                     "Atom name in mol2 file contains only digits.\n");
                   exit(1);
              }

              for(j = 1; j < k; j++ ) {
                   atom_name[i][j-1] = atom_name[i][j];
              }      /* end for( j = 1; j < strlen(&atom_name[i][0]) */

              atom_name[i][k-1] = tempt;
          }          /* end while(isdigit(atom_name[i][0])) */

/* now we should have an acceptable atom name */

     }              /* end  for( i = 0; i < read_atomnum, i++ ) */

     if(fgets(line, MAXCHAR, fpin) == NULL) {
          my_fatal(FARGS,
            "Premature termination of mol2 file before bonds record.\n");
          exit(1);          /* premature termination of file */
     }

     if(strncmp(line, "@<TRIPOS>BOND", 13)) {
          my_fatal(FARGS,
            "Line after last atom record not @<TRIPOS>BOND\n");
          exit(1);          /* premature termination of file */
     }

     for(i = 0; i < read_bondnum; i++ ) {     /* handle bond records */
          if(fgets(line, MAXCHAR, fpin) == NULL) {
               my_fatal(FARGS,
                 "Premature termination of mol2 file while reading bonds.\n");
               exit(1);          /* premature termination of file */
          }

          if(!strncmp(line, "@<TRIPOS>", 9)) {
               my_fatal(FARGS,
                 "Premature end of bond records in mol2 file.\n");
               exit(1);          /* not enough bonds */
          }

          sscanf(line, "%d%d%d%s", &tmpint3, &bond_i[i],
                 &bond_j[i], type);

          if((bond_i[i] < 1) || (bond_i[i] > read_atomnum)) {
               sprintf(mess,
                 "For bond number %d, first bonded atom number, %d, out of range\n.",
                       (i+1), bond_i[i]);
               my_fatal(FARGS, mess);
               exit(1);
          }

          if((bond_j[i] < 1) || (bond_j[i] > read_atomnum)) {
               sprintf(mess,
                 "For bond number %d, second bonded atom number, %d, out of range\n.",
                       (i+1), bond_j[i]);
               my_fatal(FARGS, mess);
               exit(1);
          }

          bond_i[i]--;        /* correct to start with 0 instead of 1 */
          bond_j[i]--;

          if(bond_count[(bond_i[i])] > 6) {
               sprintf(mess,
                 "Too many bonds to atom number %d with bond number %d.\n",
                       (bond_i[i] + 1), (i+1));
               my_fatal(FARGS, mess);
               exit(1);
          }

          if(bond_count[(bond_j[i])] > 6) {
               sprintf(mess,
                 "Too many bonds to atom number %d with bond number %d.\n",
                       (bond_j[i] + 1), (i+1));
               my_fatal(FARGS, mess);
               exit(1);
          }

/* make connection tables */
/* atoms bound to each other side 1 */
          connect[(bond_i[i])][(bond_count[(bond_i[i])])] = bond_j[i];
/* bonds atoms participate in side 1 */
          bonds_list[(bond_i[i])][(bond_count[(bond_i[i])])] = i;
          bond_count[(bond_i[i])]++;

/* atoms bound to each other side 2 */
          connect[(bond_j[i])][(bond_count[(bond_j[i])])] = bond_i[i];
/* bonds atoms participate in side 2 */
          bonds_list[(bond_j[i])][(bond_count[(bond_j[i])])] = i;
          bond_count[(bond_j[i])]++;

/* get the bond type */
          if(len = strlen(type) > 2) {
               if(!strncmp(type, "SINGLE", 6))
                    strcpy(type, "1");
               if(!strncmp(type, "DOUBLE", 6))
                    strcpy(type, "2");
               if(!strncmp(type, "TRIPLE", 6))
                    strcpy(type, "3");
          }

          strncpy(&bond_type[i][0], type, 5);
          bond_type[i][5] = NULLCHAR;   /*  with strncpy, must always ensure
                                            the string is terminated */

          len = strlen(&bond_type[i][0]);
          if(len > 2) {
               sprintf(mess,
                       "Untranslatable bond type %s\n",
                       bond_type[i]);
               my_fatal(FARGS, mess);
               exit(1);
          }
          if(len == 1) {
               switch(bond_type[i][0]) {
                    case '1':  how_bonded[i] = 1;
                               break;
                    case '2':  how_bonded[i] = 2;
                               break;
                    case '3':  how_bonded[i] = 3;
                               break;
                    default:  sprintf(mess,
                                      "Untranslatable bond type %s\n",
                                       bond_type[i]);
                              my_fatal(FARGS, mess);
                              exit(1);
               }
               continue;
          }

          if(!strncmp(&bond_type[i][0], "ar", 2)) {
               how_bonded[i] = 10;
               continue;
          }
          if(!strncmp(&bond_type[i][0], "am", 2)) {
               how_bonded[i] = 1;
               continue;
          }

/* failed to match anything */
          sprintf(mess,
                  "Untranslatable bond type %s\n",
                  bond_type[i]);
          my_fatal(FARGS, mess);
          exit(1);
 
     }        /* end for( i = 0; i < read_bondnum, i++ ) */

     if(fgets(line, MAXCHAR, fpin) == NULL) {
          my_fatal(FARGS,
                   "Incomplete mol2 file after bonds record.\n");
          exit(1);          /* premature termination of file */
     }

     if(strncmp(line, "@<TRIPOS>", 9)) {
          my_warning(FARGS,
            "Bond record does not end with @<TRIPOS>\n");
          feat = 0;                     /* no restraints if bad ending */
     }

     if(use_charge) {                   /* adjust charges to be integral */
          adjust_chrg(read_atomnum);
          adjust_chrg(read_atomnum);    /* second time in case to round off errors */
     }                                  /* end if(use_charge) */

/* now need to analyze the bonding for consistency */
     for(i = 0; i < read_atomnum; i++) {
          allowed = 0;
          charges = 0;
          mandatory = 0;
          bndcheck = 0;

          switch(atom_atno[i]) {
               case 0:   allowed = 1;       /* LP */
                         mandatory = 0;     /* sum of bonds and charges must be at least */
                         bndcheck = 1;
                         break;

               case 1:   allowed = 1;       /* H */
                         mandatory = 1;     /* sum of bonds and charges must be at least */
                         bndcheck = 1;
                         break;

               case 9:                      /* F */
               case 17:                     /* Cl */
               case 35:                     /* Br */
               case 53:  allowed = 1;       /* I */
                         charges = -1;
                         mandatory = 1;     /* sum of bonds and charges must be at least */
                         bndcheck = 1;
                         break;

               case 6:   for(j = CSTART; j < NSTART; j++) {       /* C */
                              if(!strcmp(atom_type[i], atom_desc[j].name)) {
                                   mandatory = atom_desc[j].minbonds;
                                   allowed = atom_desc[j].maxbonds;
                                   break;
                              }
                         }
                         if((bond_count[i] == 3) &&
                             !strcmp(atom_type[i], "C.2"))
                                   if(!strncmp(atom_type[(connect[i][0])], "N", 1) ||
                                      !strncmp(atom_type[(connect[i][1])], "N", 1) ||
                                      !strncmp(atom_type[(connect[i][2])], "N", 1))
                                           allowed = 6;
                         bndcheck = 1;
                         break;

               case 7:   for(j = NSTART; j < OSTART; j++) {       /* N */
                              if(!strcmp(atom_type[i], atom_desc[j].name)) {
                                   mandatory = atom_desc[j].minbonds;
                                   allowed = atom_desc[j].maxbonds;
                                   break;
                              }
                         }
                         if((!strcmp(atom_type[i], "N.4") ||
                             !strcmp(atom_type[i], "N.pl3") ||
                             !strcmp(atom_type[i], "N.am") ||
                             !strcmp(atom_type[i], "N.2")) &&
                            (atom_charge[i] < 0.5)) {
                                   allowed = 4;
                                   if(atom_charge[i] < 0.0)
                                        allowed = 5;
                         }
                         charges = 1;
                         bndcheck = 1;
                         break;

               case 8:   for(j = OSTART; j < PSTART; j++) {       /* O */
                              if(!strcmp(atom_type[i], atom_desc[j].name)) {
                                   mandatory = atom_desc[j].minbonds;
                                   allowed = atom_desc[j].maxbonds;
                                   break;
                              }
                         }
                         if((!strcmp(atom_type[i], "O.co2") ||
                             !strcmp(atom_type[i], "O.2") ||
                             !strcmp(atom_type[i], "O.3")) &&
                            (atom_charge[i] < -0.5))
                                   allowed = 3;
                         charges = -1;
                         bndcheck = 1;
                         break;

               case 15:  allowed = 5;       /* P */
                         mandatory = 3;     /* sum of bonds and charges must be at least */
                         bndcheck = 1;
                         break;

               case 16:  for(j = SSTART; j < (SSTART + SLEN); j++) {       /* S */
                              if(!strcmp(atom_type[i], atom_desc[j].name)) {
                                   mandatory = atom_desc[j].minbonds;
                                   allowed = atom_desc[j].maxbonds;
                                   break;
                              }
                         }
                         bndcheck = 1;
                         break;

               default:  bndcheck = 0;
                         break;             /* all others default to un-bonded */
          }

          if(bndcheck) {
               if(bond_count[i] > allowed) {
                    sprintf(mess,
                       "Atom %d (%s) has %d connections when allowed %d\n",
                            (i+1), atom_name[i], bond_count[i], allowed);
                    my_fatal(FARGS, mess);
                    exit(1);
               }
               if(bond_count[i] < mandatory) {
                    sprintf(mess,
                       "Atom %d (%s) has %d connections when required to have %d\n",
                            (i+1), atom_name[i], bond_count[i], mandatory);
                    my_fatal(FARGS, mess);
                    exit(1);
               }
          }
     }          /* end for(i = 0; i < read_atomnum; i++) */

/* last need to read the different kinds of restraints */
     if(feat > 0) {                        /* possible restraints? */
          which = 0;
          max = feat;
          tmpint1 = 0;
          while(fgets(line, MAXCHAR, fpin) != NULL) {

/* simple distance restraints: get   atom1, atom2, dist, force_const */
               if(!strncmp(line, "@<TRIPOS>FFCON_DIST", 19)) {
                    if(which & 1) {
                         my_fatal(FARGS,
                           "Second FFCON_DIST tag in file\n");
                         exit(1);
                    }
                    which = (which || 1);
                    i = read_restr(fpin, &rstr_dist, 1, read_atomnum,
                                   max, &tmpint1);
                    if(!i) break;
               }

/* distance range restraints: get  atom1, atom2, low_dist, high_dist, force_const */
               if(!strncmp(line, "@<TRIPOS>FFCON_RANGE", 20)) {
                    if(which & 2) {
                         my_fatal(FARGS,
                           "Second FFCON_RANGE tag in file\n");
                         exit(1);
                    }
                    which = (which || 2);
                    i = read_restr(fpin, &rstr_range, 2, read_atomnum,
                                   max, &tmpint1);
                    if(!i) break;
               }

/* angle restraints: get  atom1, atom2, atom3, theta, force_const */
               if(!strncmp(line, "@<TRIPOS>FFCON_ANGLE", 20)) {
                    if(which & 4) {
                         my_fatal(FARGS,
                           "Second FFCON_ANGLE tag in file\n");
                         exit(1);
                    }
                    which = (which || 4);
                    i = read_restr(fpin, &rstr_angle, 3, read_atomnum,
                                   max, &tmpint1);
                    if(!i) break;
               }

/* torsion restraints: get  atom1, atom2, atom3, atom4, phase, force_const */
               if(!strncmp(line, "@<TRIPOS>FFCON_TORSION", 22)) {
                    if(which & 8) {
                         my_fatal(FARGS,
                           "Second FFCON_TORSION tag in file\n");
                         exit(1);
                    }
                    which = (which || 8);
                    i = read_restr(fpin, &rstr_torsn, 4, read_atomnum,
                                   max, &tmpint1);
                    if(!i) break;
               }
          }                        /* end while(fgets(line, MAXCHAR, fpin) != NULL) */
     }                             /* end if(feat > 0) */

     *atomnum = read_atomnum;
     *bondnum = read_bondnum;

     fclose(fpin);          /* our work is done here */

/* exits with atom names, types, masses, charges, coordinates,
   residue names, atomic symbol, atomic number, and residue numbers,
   in atom_name, atom_type, atom_mass, atom_charge, (atom_x, atom_y,
   atom_z), atom_residue, atom_atno, atom_symbl, and atom_segno
   respectively; with bond information in bond_i, bond_j, bond_type,
   and how_bonded; and with a count of atoms bound to each atom in
   bond_count, a table of atom connections to each other in connect,
   and a table of bonds in which each atom participates in bonds_list */

     return;                  /* I like explicit returns */

/* WARNING     WARNING     WARNING     WARNING     WARNING     WARNING     WARNING */
/* The following memory allocations were made and not freed in this code:

        atom_name, atom_type, atom_residue, bond_type, atom_name[i],
        atom_type[i], atom_residue[i], bond_type[i], bond_i, bond_j,
        bond_count, connect, connect[i], bond_list, bonds_list[i],
        atom_segno, atom_x, atom_y, atom_z, atom_charge, atom_mass,
		how_bonded, atom_atno, atom_symbl, and atom_symbl[i]

   These will all need to be freed at some point. */
/* WARNING     WARNING     WARNING     WARNING     WARNING     WARNING     WARNING */
}


/* write MOL2 format
   parameters are:
        filename             name of .mol2 file to write
        atomnum              number of atoms
        bondnum              number of bonds
        names                list of names to write
*/
void writemol2(const char *filename, int atomnum, int bondnum, char **names)
{
     int i;
     FILE *fpout;

     if((fpout = fopen(filename, "w")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     fprintf(fpout, "@<TRIPOS>MOLECULE\n");

     fprintf(fpout, "%s\n", mol_name);

     fprintf(fpout, "%5d%6d%6d%6d%6d\n", atomnum, bondnum, 1, 0, 0);

     fprintf(fpout, "%s\n", model_type);
 
     fprintf(fpout, "%s\n\n\n", charge_type);

     fprintf(fpout, "@<TRIPOS>ATOM\n");

     for(i = 0; i < atomnum; i++) {
          fprintf(fpout,
                  "%7d %-8s%10.4f%10.4f%10.4f %-6s%5d %-8s%10.4f\n",
                  (i + 1), names[i], atom_x[i], atom_y[i], atom_z[i],
                  atom_type[i], atom_segno[0], atom_residue[0],
                  atom_charge[i]);
     }        /* end for (i = 0; i < atomnum; i++) */

     fprintf(fpout, "@<TRIPOS>BOND\n");

     for(i = 0; i < bondnum; i++) {
          fprintf(fpout, "%6d%5d%5d %s\n", (i + 1),
                  (bond_i[i] + 1), (bond_j[i] + 1), bond_type[i]);
     }                            /* end for (i = 0; i < bondnum; i++) */

     fprintf(fpout, "@<TRIPOS>SUBSTRUCTURE\n");
     fprintf(fpout, "%s%-12s%s", "     1 ", atom_residue[0],
	     "1 TEMP              0 ****  ****    0 ROOT\n");

     fclose(fpout);           /* its all out now */

     return;                  /* I like explicit returns */
}

/* reads simple distance restraints, range restraints, angle restraints
   and torsion restraints.  Returns 0 if reached end of file & 1 otherwise.
   parameters are:
        fpin                         file pointer to input file
        num_restr                    return of number of restraints
        restr_type                   type of restraint being read
        num_atoms                    number of atoms
*/
int read_restr(FILE *fpin, int *num_restr, int restr_type, int num_atoms,
               int max, int *incount)
{
     int count, start;

     count = *incount;
     start = *incount;
     while(fgets(line, MAXCHAR, fpin) != NULL) {
          if(!strncmp(line, "@<TRIPOS>", 9)) {
               *num_restr = count - start;
               *incount = count;
               return(1);
          }
          switch(restr_type) {
/* simple distance restraint: only need atom1, atom2, dist, force_const */
               case 1:  if(count >= max) {
                             my_fatal(FARGS,
                               "Maximum allowance for FFCON_DIST exceeded\n");
                             exit(1);
                        }
                        sscanf(line, "%d%d%f%f", &restr[count][1], &restr[count][2],
                               &frest[count][0], &frest[count][1]);
                        if((restr[count][1] > num_atoms) ||
                           (restr[count][2] > num_atoms)) {
                                  my_fatal(FARGS,
                                    "Atom number out of range in distance restraint\n");
                                  exit(1);
                        }
                        restr[count][0] = 1;
                        count++;
                        break;

/* range of distances restraint: only need  atom1, atom2, low_dist,
   high_dist, force_const */
               case 2:  if(count >= max) {
                             my_fatal(FARGS,
                               "Maximum allowance for FFCON_RANGE exceeded\n");
                             exit(1);
                        }
                        sscanf(line, "%d%d%f%f%f", &restr[count][1],
                               &restr[count][2], &frest[count][0],
                               &frest[count][1], &frest[count][2]);
                        if((restr[count][1] > num_atoms) ||
                           (restr[count][2] > num_atoms)) {
                                  my_fatal(FARGS,
                                    "Atom number out of range in range restraint\n");
                                  exit(1);
                        }
                        restr[count][0] = 2;
                        count++;
                        break;

/* angle restraint:  atom1, atom2, atom3, theta, force_const */
               case 3:  if(count >= max) {
                             my_fatal(FARGS,
                               "Maximum allowance for FFCON_ANGLE exceeded\n");
                             exit(1);
                        }
                        sscanf(line, "%d%d%d%f%f", &restr[count][1],
                               &restr[count][2], &restr[count][3],
                               &frest[count][0], &frest[count][1]);
                        if((restr[count][1] > num_atoms) ||
                           (restr[count][2] > num_atoms) ||
                           (restr[count][3] > num_atoms)) {
                                  my_fatal(FARGS,
                                    "Atom number out of range in angle restraint\n");
                                  exit(1);
                        }
                        restr[count][0] = 3;
                        count++;
                        break;

/* torsion restraint:   atom1, atom2, atom3, atom4, phase, force_const */
               case 4:  if(count >= max) {
                             my_fatal(FARGS,
                               "Maximum allowance for FFCON_TORSION exceeded\n");
                             exit(1);
                        }
                        sscanf(line, "%d%d%d%d%f%f", &restr[count][1],
                               &restr[count][2], &restr[count][3],
                               &restr[count][4], &frest[count][0],
                               &frest[count][1]);
                        if((restr[count][1] > num_atoms) ||
                           (restr[count][2] > num_atoms) ||
                           (restr[count][3] > num_atoms) ||
                           (restr[count][4] > num_atoms)) {
                                  my_fatal(FARGS,
                                    "Atom number out of range in torsion restraint\n");
                                  exit(1);
                        }
                        restr[count][0] = 4;
                        count++;
                        break;

               default: my_fatal(FARGS,
                          "Unrecognized restraint type in input\n");
                        exit(1);
                        break;
          }
     }

     *num_restr = count - start;
     *incount = count;

     return(0);
}
