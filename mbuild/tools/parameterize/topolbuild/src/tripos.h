/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

        tripos.h     --   tripos atom descriptors
*/

typedef struct {
     int atno;
     char name[12];
     int maxbonds;
     int minbonds;
}    TRIPOS;

static TRIPOS atom_desc[] = {
     {  6, "C.3",   4, 4 },
     {  6, "C.2",   3, 3 },
     {  6, "C.1",   2, 2 },
     {  6, "C.ar",  3, 3 },
     {  6, "C.cat", 3, 2 },
     {  7, "N.3",   4, 3 },
     {  7, "N.2",   3, 2 },
     {  7, "N.1",   2, 1 },
     {  7, "N.am",  3, 2 },
     {  7, "N.pl3", 3, 3 },
     {  7, "N.4",   4, 2 },
     {  7, "N.ar",  3, 1 },
     {  8, "O.3",   4, 1 },
     {  8, "O.2",   3, 1 },
     {  8, "O.co2", 3, 1 },
     {  8, "O.spc", 4, 2 },
     {  8, "O.t3p", 4, 2 },
     { 15, "P.3",   5, 3 },
     { 16, "S.3",   4, 2 },
     { 16, "S.2",   3, 1 },
     { 16, "S.o",   4, 2 },
     { 16, "S.o2",  4, 2 }
};
 
#define CSTART  0
#define CLEN    5
#define NSTART  5
#define NLEN    7
#define OSTART 12
#define OLEN    5
#define PSTART 17
#define PLEN    1
#define SSTART 18
#define SLEN    4
