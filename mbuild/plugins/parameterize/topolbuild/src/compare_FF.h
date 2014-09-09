/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        force field comparison header  compare_FF.h
*/

#include "similars.h"
#include "mol2.h"
#include "amber_FF.h"
#include "GMX_FF.h"
#include "block_memory.h"

#define MAXCHAR 256

extern char mess[MAXCHAR];                 /* for error messages */

extern int compare_bonds(char *bond_atom1, char *bond_atom2,
                         int bind_cntFF, int where);
extern int angles_comp(char *angle_atom1, char *angle_atom2,
                       char *angle_atom3, int angle_cntFF,
                       int where, int method);
extern int compare_angles(char *angle_atom1, char *angle_atom2,
                          char *angle_atom3, int angle_cntFF,
                          float *force, float *angle, int *method);
extern void calc_angle(int ang1, int ang2, int ang3, float *force,
                       float *angle, float *bnd_length, int *method,
                       int *sim_atom_indx, char **corr_name1,
                       char **corr_name2, char **corr_name3,
                       int angle_cntFF, char **assigned_type);
extern int compare_tors(char *tors_atom1, char *tors_atom2,
                        char *tors_atom3, char *tors_atom4,
                        int tors_cntFF, int *mult, float *force,
                        float *phase, float *term);
extern double force_param(int atomic);
extern int compare_improp(char *impro_atom1, char *impro_atom2,
                          char *impro_atom3, char *impro_atom4,
                          int impro_cntFF, int *mult, float *force,
                          float *phase, float *term);
extern int atom_cmp(char *my_atom, char *at_type);
