/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        angle and torsion tables set up header    param_tables.h
*/

extern void find_angles(int just_count, int max_ang, int num_atoms,
                        int num_bonds, int *num_angs);
extern void find_torsions(int just_count, int max_tors, int num_atoms,
                          int num_bonds, int *num_tors);
extern void measure_angles(int num_angs);
extern double angle_comp(int i, int j, int k);
extern void measure_torsions(int num_tors);
extern void measure_bondlen(int num_bonds);
extern void measure_improp(int num_improp);
extern void fix_tor(int count, int num_impr);
extern void correct_dihed_types(int num_tors);
extern double dihed(int i, int j, int k, int l);
extern void set_wanted_tors(int tors_max, int purge_level);
extern void count_torsions(int num_tors);
extern int single_flat_ring(int i);
extern int cross_flat_ring(int i);
extern void amber_flat_rings(int num_tors);
extern void planar_angles(int num_angles);
