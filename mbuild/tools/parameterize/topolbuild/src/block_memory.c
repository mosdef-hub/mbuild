/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA
        Memory allocation and error handling functions    block_memory.c

   NOTICE: The error handler is derivative work based on study of
   the routines found GROMACS version 3.3.1 written by David
   van der Spoel, Erik Lindahl, Berk Hess, and others and
   copyright (c) 1991-2000, University of Groningen, The Netherlands.
   Copyright (c) 2001-2004, The GROMACS development team
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "block_memory.h"

#define NULLCHAR 0

/* do vector allocations and handle allocation errors
   parameters are:
        quantity             how much to allocate
        size                 size of each memory unit
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
*/
char *allocator(int quantity, int size, const char *file, int line)
{
     char *temp;

     if(!(temp = (char *)malloc((unsigned) quantity * size))) {
          my_fatal(file, line ,
                   "Memory allocation failure.\n");
          exit(1);
     }

     return(temp);
}

/* allocate character matrices
   parameters are:
        quant1               number of strings to allocate
        quant2               length of each character string
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
*/
char **mat_alloc(int quant1, int quant2, const char *file, int line)
{
     char **temp;
     int i;

     temp = (char **)allocator(quant1, sizeof(char*), file, line);

     for(i = 0; i < quant1; i++) {
          temp[i] = (char *)allocator(quant2, sizeof(char), file, line);
          temp[i][0] = NULLCHAR;
     }        /* end for( i = 0; i < quant1; i++) */

     return(temp);
}

/* allocate a table of pointers to pointers to char
   parameters are:
        quant1               number of blocks of strings needed
        quant2               number of strings to allocate
        quant3               length of each character string
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
*/
char ***tabl_alloc(int quant1, int quant2, int quant3, const char *file,
                   int line )
{
     char ***temp;
     int i;

     temp = (char ***)allocator( quant1, sizeof(char*), file, line);

     for(i = 0; i < quant1; i++) {
          temp[i] = (char **)mat_alloc( quant2, quant3, file, line);
          temp[i][0][0] = NULLCHAR;
     }        /* end for( i = 0; i < quant1; i++) */

     return(temp);
}

/* allocate an alignment of pointers to pointers to pointers to char
   parameters are:
        quant1               number of pointers to blocks of strings needed
        quant2               number of blocks of strings needed
        quant3               number of strings to allocate
        quant4               length of each character string
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
*/
char ****align_alloc(int quant1, int quant2, int quant3, int quant4,
                     const char *file, int line )
{
     char ****temp;
     int i;

     temp = (char ****)allocator( quant1, sizeof(char*), file, line);

     for(i = 0; i < quant1; i++) {
          temp[i] = (char ***)tabl_alloc( quant2, quant3, quant4, file, line);
          temp[i][0][0][0] = NULLCHAR;
     }        /* end for( i = 0; i < quant1; i++) */

     return(temp);
}

/* free allocated memory with check that it was allocated
   parameters are:
        my_alloc               pointer to allocated memory
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
*/
void free_me(void *my_alloc, const char *file, int line )
{
     if(!my_alloc) {
          my_fatal(file, line ,
                   "Attempt to free unallocated memory.\n");
          exit(1);
     }
     free(my_alloc);
     return;
}

/* error handler
   parameters are:
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
        msg                  error message to write
*/
void my_fatal(const char *file, int line, const char *msg)
{
     if(stdlog) {
          fprintf(stdlog, "Fatal error.\n");
          fprintf(stdlog, "Source code file: %s, line: %d\n", file, line);
          fprintf(stdlog,"%s\n",msg);
          fflush(stdlog);
     }
     fprintf(stderr, "Fatal error.\n");
     fprintf(stderr, "Source code file: %s, line: %d\n", file, line);
     fprintf(stderr,"%s\n",msg);
     exit(-1);
}

/* warning handler
   parameters are:
        file                 the name of the calling code file
                             for errors
        line                 line number in the calling code
                             for errors
        msg                  error message to write
*/
void my_warning(const char *file, int line, const char *msg)
{
     if(stdlog) {
          fprintf(stdlog, "Warning.\n");
          fprintf(stdlog, "Source code file: %s, line: %d\n", file, line);
          fprintf(stdlog,"%s\n",msg);
          fflush(stdlog);
     }
     fprintf(stderr, "Warning.\n");
     fprintf(stderr, "Source code file: %s, line: %d\n", file, line);
     fprintf(stderr,"%s\n",msg);
     return;
}
