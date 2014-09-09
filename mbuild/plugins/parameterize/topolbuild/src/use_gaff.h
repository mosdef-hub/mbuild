/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        use the gaff force field data to set parameters header    use_gaff.h
*/

extern void set_gaff_atomFF(int num_atoms, int num_mol);
extern void set_gaff_bondFF(int num_bonds, int bind_cntFF);
extern void set_gaff_vdwFF(int num_atoms, int num_vdw);
extern void set_gaff_angleFF(int num_angles, int angle_cntFF);
extern void set_gaff_torsFF(int num_tors, int tors_cntFF);
extern void set_gaff_improFF(int num_impr, int impr_cntFF);
extern void check_more_tors(int where, int what);
