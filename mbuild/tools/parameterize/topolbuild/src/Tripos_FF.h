/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Define additional storage elements for Tripos type forcefields
        and read tripos forcefield TAFF files
        and use Tripos force field data to set parameters header
        Tripos_FF.h
*/

extern void read_tripos(char *dir, int just_count, int *num_bind, int *num_angle,
                        int *num_tors, int *num_impro, int *num_vdw);
extern void free_triposFF(int num_bind, int num_angle, int num_tors, int num_impro,
                          int num_vdw);

extern void set_tripos_bondFF(int num_bonds, int bind_cntFF);
extern void set_tripos_angleFF(int num_angles, int angle_cntFF);
extern void set_tripos_torsFF(int num_tors, int tors_cntFF);
extern void set_tripos_improFF(int num_impr, int impr_cntFF);
extern int compare_types(int id_atom, char *name_of_type);

/* can use set_amber_vdwFF for the van der Waals */

/* Tripos type force field additional parameters */

extern char **trp_bnd_type;
extern char **tors_bond;

/* end Tripos type force field additional parameters */
