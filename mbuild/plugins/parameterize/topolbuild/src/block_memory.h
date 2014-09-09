/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        memory allocation and error handling routines header  block_memory.h
*/

extern char *allocator(int quantity, int size, const char *file, int line );
extern char **mat_alloc(int quant1, int quant2, const char *file, int line );
extern char ***tabl_alloc(int quant1, int quant2, int quant3, const char *file,
                          int line );
extern char ****align_alloc(int quant1, int quant2, int quant3, int quant4,
                            const char *file, int line );
extern void free_me(void *my_alloc, const char *file, int line );
extern void my_fatal(const char *file, int line, const char *msg);
extern void my_warning(const char *file, int line, const char *msg);

#define FARGS __FILE__,__LINE__

extern FILE *stdlog;
