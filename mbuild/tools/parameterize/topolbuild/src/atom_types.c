/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Functions to read antechamber atom types and atom correlation
        data files and to free the storage allocated in reading them.

   NOTICE: This is a derivative work based on study of the routines
   found in the antechamber 1.27 program atomtype, version 1.0,
   dated October, 2001 by Junmei Wang, Department of Pharmaceutical
   Chemistry, School of Pharmacy, University of California,
   San Francisco, CA  94143

   Portions of this work simplify storage allocation, and clarify
   decision trees compared to the work cited above.
*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "block_memory.h"
#include "atom_types.h"
#include "readmol2.h"
#include "utilities.h"

/* read an antechamber style atom type definition file
   parameters are:
        at_type_cnt                   returned count of atom types
        wild_cnt                      returned count of wild atoms
        max_type                      maximum number of atom types allowed
        max_wild                      maximum number of wild atoms allowed
        just_count                    count records and return the counts
                                      but do not allocate storage or store data
        filename                      name of the atom type definition file
                                      to read
*/
void read_at_types(int *at_type_cnt, int *wild_cnt, int max_type, int max_wild,
                   int just_count, char *filename)
{
     int i, j, k, num_wild, num_type, wild_allow;
     int lev0, lvls;
     float dummy;
     FILE *def;
     char *dd;
     char atype[10];

     if(just_count) {
          wild_name = (char **)mat_alloc(max_wild, 8, FARGS);
          wild_count = (int *)allocator(max_wild, sizeof(int), FARGS);
          at_type_name = (char **)mat_alloc(max_type, 12, FARGS);
          at_type_residue = (int *)allocator(max_type, sizeof(int), FARGS);
          residue_name = (char **)allocator(max_type, sizeof(char *), FARGS);
          at_type_no = (int *)allocator(max_type, sizeof(int), FARGS);
          at_type_attached = (int *)allocator(max_type, sizeof(int), FARGS);
          at_type_attached_H = (int *)allocator(max_type, sizeof(int), FARGS);
          at_type_ewd = (int *)allocator(max_type, sizeof(int), FARGS);
          at_type_prop = (char **)allocator(max_type, sizeof(char *), FARGS);
          at_type_env = (char **)allocator(max_type, sizeof(char *), FARGS);
          at_env_bonds = (char **)allocator(max_type, sizeof(char *), FARGS);
          wild_elements = (char ***)allocator(max_wild, sizeof(char*), FARGS);
          wild_atno = (int **)allocator(max_wild, sizeof(int*), FARGS);
          for(i = 0; i < max_type; i++) {
               at_type_residue[i] = 0;
               at_type_no[i] = 0;
               at_type_attached[i] = 999;
               at_type_attached_H[i] = 999;
               at_type_ewd[i] = 999;
               at_type_prop[i] = NULL;
               at_type_env[i] = NULL;
               at_env_bonds[i] = NULL;
               residue_name[i] = NULL;
          }
          for(i = 0; i < max_wild; i++) {
               wild_count[i] = 0;
               wild_elements[i] = 0;
          }
     }

     num_wild = 0;
     num_type = 0;

     if((def = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     wild_allow = 1;

     while(fgets(line, ((2 * MAXCHAR) - 2), def) !=NULL) {
          line[((2 * MAXCHAR) - 1)] = NULLCHAR;
          if((!strncmp(line, "WILDATOM", 8)) && wild_allow) {
               if(num_wild >= max_wild) {
                    sprintf(mess,
                            "More WILDATOM records (%d) than allowed (%d)\n",
                            (num_wild+1), max_wild);
                    my_fatal(FARGS, mess);
                    exit(1);
               }
               if((i = make_tokens(name, &line[9], 11)) < 2) {
                    my_fatal(FARGS, "Faulty WILDATOM definition\n");
                    exit(1);
               }

/* only allocate the memory I need and store a count of elements for each WILDATOM */

               if(just_count) {
                    strncpy(&wild_name[num_wild][0], name[0], 4);
                    wild_name[num_wild][4] = NULLCHAR;
                    wild_elements[num_wild] = (char **)mat_alloc((i-1), 6, FARGS);
                    wild_atno[num_wild] = (int *)allocator((i-1), sizeof(int), FARGS);
                    wild_count[num_wild] = i-1;
                    for(j = 0; j < (i-1); j++) {
                         strncpy(wild_elements[num_wild][j], name[(j+1)], 4);
                         wild_elements[num_wild][j][4] = NULLCHAR;
                         k = strcspn(name[(j+1)], "0123456789");
                         strncpy(atype, name[(j+1)], k);
                         atype[k] = NULLCHAR;
                         set_atomic_no(k, &atype[0], &dummy, &wild_atno[num_wild][j],
                                       mess);
                    }
               }
               num_wild++;
               continue;
          }                    /* end if((!strncmp(line, "WILDATOM", 8)) && wild_allow) */

          if(!strncmp(line, "ATD", 3)) {           /* atom type definitions */
               wild_allow = 0;                     /* no more WILDATOM entries */
               if(num_type >= max_type) {
                    sprintf(mess,
                            "More atom type definitions (%d) than allowed (%d)\n",
                            (num_type+1), max_type);
                    my_fatal(FARGS, mess);
                    exit(1);
               }

               if((dd = strchr(line, '&')) == NULL) {
                    my_fatal(FARGS, "Faulty ATD line\n");
                    exit(1);
               }
               dd++;
               *dd = NULLCHAR;
               i = make_tokens(&name[2], &line[4], 9);
               if((!i) || (!strcmp(name[2], "&"))) {
                    my_fatal(FARGS, "Faulty ATD line\n");
                    exit(1);
               }
               num_type++;
               if(just_count) {
                    strncpy(&at_type_name[(num_type - 1)][0], name[2], 10);
                    at_type_name[(num_type - 1)][10] = NULLCHAR;
                    if(!strcmp(name[3], "&")) continue;
                    at_type_residue[(num_type - 1)] = -1;
                    if(strcmp(name[3], "*")) {
                         if(!strcmp(name[3], "AA"))
                              at_type_residue[(num_type - 1)] = 1;
                         if(!strcmp(name[3], "NA"))
                              at_type_residue[(num_type - 1)] = 2;
                         if(!strcmp(name[3], "BIO"))
                              at_type_residue[(num_type - 1)] = 3;
                         if(strcmp(name[3], "AA") && strcmp(name[3], "NA") &&
                            strcmp(name[3], "BIO")) {
                                at_type_residue[(num_type - 1)] = 4;
                                residue_name[(num_type - 1)] =
                                     (char *)allocator((strlen(name[3]) + 2),
                                                       sizeof(char), FARGS);
                                strcpy(&residue_name[(num_type - 1)][0], name[3]);
                         }
                    }                               /* end if(strcmp(name[3], "*")) */
                    if(!strcmp(name[4], "&")) continue;
                    if(strcmp(name[4], "*")) {
                         at_type_no[(num_type - 1)] = make_int(name[4],
                                                               "Faulty atomic number %s\n",
                                                               FARGS);
                    }
                    else at_type_no[(num_type - 1)] = -1;
                    if(!strcmp(name[5], "&")) continue;
                    if(strcmp(name[5], "*")) {
                         at_type_attached[(num_type - 1)] = make_int(name[5],
                                                    "Faulty number of attached atoms %s\n",
                                                                    FARGS);
                    }
                    else at_type_attached[(num_type - 1)] = -1;
                    if(!strcmp(name[6], "&")) continue;
                    if(strcmp(name[6], "*")) {
                         at_type_attached_H[(num_type - 1)] = make_int(name[6],
                                                        "Faulty H attachments number %s\n",
                                                                       FARGS);
                    }
                    else at_type_attached_H[(num_type - 1)] = -1;
                    if(!strcmp(name[7], "&")) continue;
                    if(strcmp(name[7], "*")) {
                         at_type_ewd[(num_type - 1)] = make_int(name[7],
                                                  "Faulty electron withdrawal number %s\n",
                                                                FARGS);
                    }
                    else at_type_ewd[(num_type - 1)] = -1;
                    if(!strcmp(name[8], "&")) {
                         at_type_prop[(num_type - 1)] =
                                                (char *)allocator(2, sizeof(char),
                                                                  FARGS);
                         at_type_prop[(num_type - 1)][0] = NULLCHAR;
                         continue;
                    }
                    if(strcmp(name[8], "*")) {
                         if(str_balanced('[', ']', name[8])) {
                              sprintf(mess,
                                      "ATD %s [ and ] do not match in %s.\n",
                                      name[2], name[8]);
                              my_fatal(FARGS, mess);
                              exit(1);
                         }
                         at_type_prop[(num_type - 1)] =
                                                (char *)allocator((strlen(name[8]) + 2),
                                                                  sizeof(char), FARGS);
                         strcpy(&at_type_prop[(num_type - 1)][0], name[8]);
                    }
                    else {
                         at_type_prop[(num_type - 1)] =
                                                (char *)allocator(2, sizeof(char),
                                                                  FARGS);
                         strcpy(&at_type_prop[(num_type - 1)][0], "*");
                    }
                    if(!strcmp(name[9], "&")) {
                         at_type_env[(num_type - 1)] =
                                                (char *)allocator(2, sizeof(char),
                                                                  FARGS);
                         at_type_env[(num_type - 1)][0] = NULLCHAR;
                         continue;
                    }
                    if(strcmp(name[9], "*")) {
                         if(str_balanced('(', ')', name[9])) {
                              sprintf(mess,
                                      "ATD %s ( and ) do not match in %s.\n",
                                      name[2], name[9]);
                              my_fatal(FARGS, mess);
                              exit(1);
                         }
                         if(str_balanced('[', ']', name[9])) {
                              sprintf(mess,
                                      "ATD %s [ and ] do not match in %s.\n",
                                      name[2], name[9]);
                              my_fatal(FARGS, mess);
                              exit(1);
                         }
                         if(str_balanced('<', '>', name[9])) {
                              sprintf(mess,
                                      "ATD %s < and > do not match in %s.\n",
                                      name[2], name[9]);
                              my_fatal(FARGS, mess);
                              exit(1);
                         }
                         lev0 = 0;
                         for(lvls = 0; lvls < strlen(name[9]); lvls++) {
                              if(lev0) {
                                   if((name[9][lvls] == '(') ||
                                      (name[9][lvls] == ')') ||
                                      (name[9][lvls] == '[') ||
                                      (name[9][lvls] == '<') ||
                                      (name[9][lvls] == '>')) {
                                             sprintf(mess,
                                              "ATD %s illegal placement of %c in %s.\n",
                                                     name[2], name[9][lvls], name[9]);
                                             my_fatal(FARGS, mess);
                                             exit(1);
                                   }
                                   if(name[9][lvls] == ']')
                                        lev0 = 0;
                              }
                              else
                                   if(name[9][lvls] == '[')
                                        lev0 = lvls + 1;
                         }         /* end for(lvle = 0; lvls < strlen(name[9]; lvls++) */

                         at_type_env[(num_type - 1)] =
                                                (char *)allocator((strlen(name[9]) + 2),
                                                                  sizeof(char), FARGS);
                         strcpy(&at_type_env[(num_type - 1)][0], name[9]);
                    }
                    else {
                         at_type_env[(num_type - 1)] =
                                                (char *)allocator(2, sizeof(char),
                                                                  FARGS);
                         strcpy(&at_type_env[(num_type - 1)][0], "*");
                    }
                    if(!strcmp(name[10], "&")) {
                         at_env_bonds[(num_type - 1)] =
                                                (char *)allocator(2, sizeof(char),
                                                                  FARGS);
                         at_env_bonds[(num_type - 1)][0] = NULLCHAR;
                         continue;
                    }
                    if(str_balanced('[', ']', name[10])) {
                         sprintf(mess,
                                 "ATD %s [ and ] do not match in %s.\n",
                                 name[2], name[10]);
                         my_fatal(FARGS, mess);
                         exit(1);
                    }
                    at_env_bonds[(num_type - 1)] =
                                           (char *)allocator((strlen(name[10]) + 2),
                                                             sizeof(char), FARGS);
                    strcpy(&at_env_bonds[(num_type - 1)][0], name[10]);
               }                   /* end if(just_count) */
          }                        /* end if(!strncmp(line, "ATD", 3)) */
     }                             /* end while(fgets(line, 2*MAXCHAR, def) !=NULL) */

     *wild_cnt = num_wild;
     num_wilds = num_wild;
     *at_type_cnt = num_type;
     fclose(def);

     return;
}

/* read an antechamber style atom types correlation file
   parameters are:
        corr_cnt                  returned count of number of atom types
                                  correlation records
        max_corr                  maximum number of atom types correlations
                                  allowed
        just_count                count records and return the counts
                                  but do not allocate storage or store data
        filename                  name of the atom types correlations file
                                  to read
*/
void read_at_corr(int *corr_cnt, int max_corr, int just_count,
                  char *filename)
{
     int i, j, k, count;
     FILE *corr;

     if(just_count) {
          corr_name = (char **)mat_alloc(max_corr, 4, FARGS);
          index_improper = (int *)allocator(max_corr, sizeof(int), FARGS);
          corr_to = (char ***)tabl_alloc(max_corr, 3, 4, FARGS);
          corr_num = (int *)allocator(max_corr, sizeof(int), FARGS);
          for(i = 0; i < max_corr; i++)
               corr_num[i] = 0;
     }

     if((corr = fopen(filename, "r")) == NULL) {
          sprintf(mess, "Cannot open file %s\n", filename);
          my_fatal(FARGS, mess);
          exit(1);
     }

     count = 0;
     while(fgets(line, MAXCHAR, corr) !=NULL) {
          if(!strncmp(line, "CORR", 4)) {
               if((i = make_tokens(name, &line[4], 5)) < 2) {
                    my_fatal(FARGS,
                             "Faulty atom correlation definition\n");
                    exit(1);
               }
               if(i == 3)
                    if(!strcmp(name[0], name[2])) {
                         name[2][0] = NULLCHAR;
                         i = 2;
                    }

               if(just_count) {
                    strncpy(&corr_name[count][0], name[0], 3);
                    corr_name[count][3] = NULLCHAR;
                    if(strlen(corr_name[count]) == 1) {
                         corr_name[count][1] = ' ';
                         corr_name[count][2] = NULLCHAR;
                    }
                    corr_num[count] = i - 2;
                    index_improper[count] = make_int(name[1],
                                                  "Faulty atom impropers index %s\n",
                                                     FARGS);
                    for(j = 0; j < (i-2); j++) {
                         strncpy(&corr_to[count][j][0], name[(j+2)], 3);
                         corr_to[count][j][3] = NULLCHAR;
                         if(strlen(corr_to[count][j]) == 1) {
                              corr_to[count][j][1] = ' ';
                              corr_to[count][j][2] = NULLCHAR;
                         }
                    }              /* end for(j = 0; j < (i-2); j++) */
               }                   /* end if(just_count) */
          count++;
          }                        /* end if(!strncmp(line, "CORR", 4)) */
     }                             /* end while(fgets(line, MAXCHAR, corr) !=NULL) */

     *corr_cnt = count;
     fclose(corr);
     return;
}

/* free storage allocated for atom types data
   parameters are:
        at_type_cnt               number of atom types
        wild_cnt                  number of wild atom records
*/
void release_types(int at_type_cnt, int wild_cnt)
{
     int i, j;

/* free wilds arrays */
     for(i = 0; i < wild_cnt; i++){
          for(j = 0; j < wild_count[i]; j++)
               free_me(wild_elements[i][j], FARGS);
          free_me(wild_name[i], FARGS);
          free_me(wild_elements[i], FARGS);
          free_me(wild_atno[i], FARGS);
     }
     free_me(wild_name, FARGS);
     free_me(wild_elements, FARGS);
     free_me(wild_atno, FARGS);
     free_me(wild_count, FARGS);

/* free atom type arrays */
     free_me(at_type_residue, FARGS);
     free_me(at_type_no, FARGS);
     free_me(at_type_attached, FARGS);
     free_me(at_type_attached_H, FARGS);
     free_me(at_type_ewd, FARGS);

     for(i = 0; i < at_type_cnt; i++) {
          free_me(at_type_name[i], FARGS);
          if(at_type_prop[i] != NULL)
               free_me(at_type_prop[i], FARGS);
          if(at_type_env[i] != NULL)
               free_me(at_type_env[i], FARGS);
          if(at_env_bonds[i] != NULL)
               free_me(at_env_bonds[i], FARGS);
          if(residue_name[i] != NULL)
               free_me(residue_name[i], FARGS);
     }
     free_me(at_type_name, FARGS);
     free_me(at_type_prop, FARGS);
     free_me(at_type_env, FARGS);
     free_me(at_env_bonds, FARGS);
     free_me(residue_name, FARGS);

     return;
}

/* free storage allocated for atom types correlations data
   parameter is:
        corr_count                number of atom types correlations stored
*/
void release_corr(int corr_count)
{
     int i, j;

     for(i = 0; i < corr_count; i++) {
          free_me(corr_name[i], FARGS);
          free_me(corr_to[i][0], FARGS);
          free_me(corr_to[i][1], FARGS);
          free_me(corr_to[i][2], FARGS);
          free_me(corr_to[i], FARGS);
     }

     free_me(corr_to, FARGS);
     free_me(corr_name, FARGS);
     free_me(index_improper, FARGS);
     free_me(corr_num, FARGS);

     return;
}
