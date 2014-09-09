/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        use the amber force field data to set parameters header   use_amber.h
*/

extern void set_amber_atomFF(int num_atoms, int num_mol);
extern void set_amber_vdwFF(int num_atoms, int num_vdw);
extern void set_amber_angleFF(int num_angles, int angle_cntFF);
extern void set_amber_bondFF(int num_bonds, int bind_cntFF);
extern void set_amber_torsFF(int num_tors, int tors_cntFF);
extern void set_amber_improFF(int num_impr, int impr_cntFF);
