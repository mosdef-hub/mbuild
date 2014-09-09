/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        write forcefield topology and include files for a topology header
        gromacs_FF.h
*/

extern void write_ff(char *commonname, char **assigned_types, int numatom,
                     float *vdw_rad, float *vdw_pot, int amber, char *altnb);
extern void write_gro (char *filename, int atomnum, char **names, char *residue);
extern void write_top(char *filename, char *residue, int numatoms, int numbonds,
                      int numpairs, int numangles, int numdihed, int numimpro,
                      char **assigned_types, char **names, int usage, int amber);
extern void cgnr_assgn(int num_atoms, char *residue);
extern void write_atoms(FILE *fptop, char *residue, int numatoms,
                        char **assigned_types, char **names);
extern void write_bonds(FILE *fptop, int numbonds, char **names, int amber,
                        double conv, int usage, int tot_rstr);
extern void write_constr(FILE *fptop, int numangles, char **names,
                         int amber, int usage);
extern void write_pairs(FILE *fptop, int numpairs, char **names);
extern void write_angles(FILE *fptop, int numangles, char **names,
                         int usage, int amber, double conv);
extern void write_dihed(FILE *fptop, int numdihed, char **names,
                        int usage, int amber, double conv);
extern void write_impropers(FILE *fptop, int numimpro, int numdihed,
                            char **names, int usage, int amber, double conv);
extern void write_restraints(FILE *fptop, int tot_rstr, char **names);
extern float get_dist(int point1, int point2);
