/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define storage and functions for judging atom types on the basis
        of bonding and of chemical environment      judge_types.h
*/

#include "block_memory.h"
#include "atom_types.h"
#include "mol2.h"
#include "rings.h"

#define MAXCHAR 256

#define MAXLEVELS 16
#define MAXCHAIN 16
#define MAXENVBOND 16
#define MAXENVSTR 48
#define MAXSCAN 64
#define NULLCHAR 0
#define NUM_FEAT 25

extern void init_atom_types(int atom_num, int num_bonds);
extern void free_atom_types(int atom_num);
extern void judge_atom_types(int num_types, int num_wild, int atom_num,
                             int num_bonds, char **assigned_type);
extern void bond_info(int atom_num, int num_bonds);
extern int check_AA(int number);
extern int check_NA(int number);
extern int check_BIO(int number);
extern int prop_check(int num, int mark, int num_bonds, char *prop_str);
extern int chem_environ(int num, char *env_str, char *bonds_str, int bndindx,
                        int atom_num, int num_bonds);
extern void check_feat(int featindx, int status, int *setindx);
extern void bond_check(char *thestr, int featindx, int status,
                       int *setindx, int mark, int num_bonds, int num,
                       int bndtyp);
extern int junctbond(int num, int mark, int num_bonds, int bndtyp);
extern void balanced_str(char cin, char cout, char *str);
extern int str_balanced(char cin, char cout, char *str);
extern void enviromatch(int selection, int startnum, int chain, int maximum,
                        int num_bonds);
extern int wildmatch(char *my_atom, int my_id);
extern void dccheck(int startpoint, int *success, int chain,  int num_atoms,
                    int num_bonds);
extern int check_chain(int chain, int num_atoms, int num_bonds);
extern int check_bonds(int ida, int idb, char *bondstr, int num_bonds);

extern char mess[MAXCHAR];

extern int scan_num;
extern int initial;

/* bond information arrays */

extern int *sb;
extern int *SB;
extern int *db;
extern int *DB;
extern int *tb;
extern int *TB;
extern int *AB;
extern int *DL;

/* end bond information arrays */

/* chemical environment arrays */

extern char **apch;
extern char **env_bonds;
extern char **env_atom1;
extern char **env_atom2;
extern char **atom_chem;
extern char **envblockstr;
extern char ***env_atom_name;
extern char ***env_ap;
extern char ***env_strname;

extern int *env_len;
extern int *chem_indx;
extern int **env_con;
extern int **env_index;

extern int *selchain;
extern int *selindx;
extern int *sel_chain_indx;

extern int *schain_id;
extern int *schain_num;
extern int **schain_atom;

extern int env_bond_num;

/* end chemical environment arrays */
