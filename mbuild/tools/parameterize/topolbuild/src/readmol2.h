/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        read and write *.mol2 (Tripos Sybyl) files
        assign atomic number and mass based on symbol given
        and caseless strncmp function header     readmol2.h

*/

#define MAXCHAR 256

extern void readmol2(const char *filename, int *atomnum, int *bondnum );
extern void alloc_mol2(int read_atomnum, int read_bondnum, int feat);
extern void writemol2(const char *filename, int atomnum, int bondnum,
                      char **names);
extern int read_restr(FILE *fpin, int *num_restr, int restr_type, int num_atoms,
                      int max, int *incount);
extern void set_atomic_no(int len, char *type, float *atmass, int *atno,
                          char *atsymb);
extern void molec_alloc(int num_atoms, int num_bonds, int num_angles,
                        int num_tors, int num_impr);
extern void atom_residue_order(int num_atoms, char *resid);

extern char line[2*MAXCHAR];               /* input line */
extern char mess[MAXCHAR];                 /* for error messages */
