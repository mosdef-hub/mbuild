/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        read amber forcefield parm files header   read_amber.h
*/

extern void read_amber(char *parmfile, int just_count, int *num_mol,
                       int *num_bind, int *num_angle, int *num_tors,
                       int *num_impro, int *num_vdw, int *num_equiv,
                       int *max_equiv);
extern void read_molparms(int just_count, int mol_max, int *num_mol);
extern void read_bindparms(int just_count, int bind_max, int *num_bind);
extern void read_angleparms(int just_count, int tors_angle, int *num_angle);
extern void read_torsparms(int just_count, int tors_max, int *num_tors);
extern void read_improperparms(int just_count, int impro_max, int *num_impro);
extern void read_vdwparms(int just_count, int vdw_max, int *num_vdw);
void read_vdwequivs(int just_count, int equiv_max, int eqlen_max, int *num_equiv,
                    int *max_equiv);
extern void free_amberFF(int num_mol, int num_bind, int num_angle,
                         int num_tors, int num_impro, int num_vdw);
extern void free_equivs(int num_equiv, int max_equiv);
extern void flag_more_tors(int cntr);
