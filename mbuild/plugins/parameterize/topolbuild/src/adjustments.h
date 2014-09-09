/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        adjust general amber forcefield atom type assignments routines header
        adjustments.h
*/

extern void adjustments(int num_atoms, int num_bonds, char **assigned_type);
extern void adjust_type_cp(int num_atoms, int num_bonds, char **assigned_type);
extern void check_type_errors(int num_atoms, int num_bonds, char **assigned_type);
extern void warn_me(int id1, int id2, char **assigned_type);
extern void adjust_ion_types(int num_atoms, char **assigned_type);
