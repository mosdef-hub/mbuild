/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

   utility functions        utilities.c
   Functions are:
        get_name               Get a name and return progress in
                               the string from which it came.
        make_tokens            Separate a string into tokens.
        caseless_strncmp       Case insensitive string comparison.
        spaceline              Check line to be sure it contains
                               graphical characters.
        uppercase              Convert string to upper case letters.
        CubeRoot               Fast cube root implementation, positive
                               numbers only.
        fast_InvSqrt           Double precision version of fast
                               inverse square root.
        get_a_line             Read lines, strip comments and end of
                               line, and return non-blank lines.
        make_int               check string for legal characters and
                               make an integer from it.
*/

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "utilities.h"

char not_token[7] = " \t\v\b\r\f\n";

/* Get a name and return progress in the string from which it came.
   Skip over spaces and tabs before the name mark first space or
   tab after the name.
   parameters are:
        this               string from which the name is to come
        place              where to store the name
        length             maximum for token length
                           (default value is 5)
*/
char *get_name(char *this, char *place, int length)
{
     int k, l, m;
     char *d;

     if(length < 1) m = 5;
     else m = length;
     place[0] = NULLCHAR;
     if((l = strspn(this, not_token)) == strlen(this))
          return(NULL);
     d = strpbrk(&this[l], not_token);
     k = strcspn(&this[l], not_token);
     if(k > m) k = m;
     strncpy(place, &this[l], k);
     place[k] = NULLCHAR;
     if(strlen(place) == 1) {
          place[1] = ' ';
          place[2] = NULLCHAR;
     }

     return(d);
}

/* Separate a line into tokens and return count of tokens.
   parameters are:
        list                 array of pointers to tokens.
        line                 line to be tokenized
        count                maximum number of tokens allowed.
*/
int make_tokens(char **list, char *line, int count)
{
     int i, k;
     char *d;

     i = 0;
     d = line;
     while((k = strspn(d, not_token)) != strlen(d)) {
          list[i] = &d[k];
          i++;
          if((d = strpbrk(&d[k], not_token)) == NULL)
               break;
          *d = NULLCHAR;
          d++;
          if(i == count)
               break;
     }
     
     return(i);
}

/* Case insensitive string comparison
   parameters are:
        string1                 string to compare
        string2                 string to compare
        length                  number of characters to compare
*/
int caseless_strncmp(const char *string1, const char *string2, int length)
{
     int i;
     char lowcase1[MAXCHAR], lowcase2[MAXCHAR];

/* if either string is less than length and the string lengths are
   not equal, then they do not compare */

     if(((strlen(string1) < length) || (strlen(string2) < length))
            && (strlen(string1) != strlen(string2))) return(1);

/* all NULLCHAR strings are equal to each other */

     if(string1[0] == NULLCHAR ) return(0);

     for(i = 0; i < length; i++) {
          if(string1[i] == NULLCHAR) break;     /* quit when hit end of strings */

          lowcase1[i] = tolower(string1[i]);
          lowcase2[i] = tolower(string2[i]);
     }          /* end for( i = 0; i < length; i++) */

     lowcase1[i] = NULLCHAR;     /* terminate strings before comparison */
     lowcase2[i] = NULLCHAR;

     i = strncmp(lowcase1, lowcase2, length);

     return(i);
}

/* check line to be sure it contains graphical characters
   only parameter is line to check
*/
int spaceline(char *line)
{
     int i;

     for(i = 0; i < strlen(line); i++)
          if(isgraph(line[i])) return(0);     /* found a graphic character in the line */

     return(1);                    /* This line looks blank when printed */
}

/* convert string to upper case letters
   parameter is the string to convert
*/
void uppercase(char *inname)
{
     int i;

     for(i = 0; i < strlen(inname); i++)
          inname[i] = toupper(inname[i]);

     return;
}

/* Fast cube root implementation for positive numbers only
   parameter is number for which cubed root is desired
   Of course, for a negative number's cubed root one can
   always use
        -CubeRoot(fabs(x));

   Derived from Turkowski, K. (1998) Computing the Cube Root,
   Apple Computer Technical Report #KT-32
*/
double CubeRoot(double x)
{
     double fr, r;
     int ex, shx;

     fr = frexp(x, &ex);
     shx = ex % 3;
     if(shx > 0)
          shx -= 3;
     ex = (ex - shx)/3;
     fr = ldexp(fr, shx);
/* Use quartic rational polynomial with error < 2^(-24) */
     fr = (((( 45.2548339756803022511987494   * fr +
              192.2798368355061050458134625)  * fr +
              119.1654824285581628956914143)  * fr +
               13.43250139086239872172837314) * fr +
                0.1636161226585754240958355063)
        /
          (((( 14.80884093219134573786480845  * fr +
              151.9714051044435648658557668)  * fr +
              168.5254414101568283957668343)  * fr +
               33.9905941350215598754191872)  * fr +
                1.0);
     r = ldexp(fr, ex);                  /* 24 bits of precision */
     r = (2.0 * r + x/(r * r))/3.0;      /* 48 bits of precision */
     r = (2.0 * r + x/(r * r))/3.0;      /* 96 bits of precision */
     return(r);
}

/* double precision version of fast inverse sqrt,
   parameter is number for which inverse sqrt is desired.
*/
double fast_InvSqrt(double x)
{
#ifdef NOIEEE
/* This is inserted for the rare case that a processor does not
   use IEEE arithmetic.
*/
     double r;
     r = 1.0/sqrt(x);
     return(r);
#else
/* Taken from Charles McEniry (2007)
   The Mathematics Benind the Fast Inverse Square Root Function
   Code

   I use two Newton Raphson steps to make sure I've got good
   double precision conversion.
*/
     const unsigned long long int R = 0xBFCDD6A18F6A6F55ULL;
     union { double f; unsigned long long int ul; } y;
     y.f = x;
/* hidden initial guess, fast */
     y.ul = (R - y.ul ) >> 1;
     y.f = 0.5 * y.f * (3.0 - x * y.f * y.f);
     y.f = 0.5 * y.f * (3.0 - x * y.f * y.f);
     return y.f;
#endif
}

/* Read lines, strip comments and end of line, and return
   non-blank lines.
   parameters are:
        where          pointer to input buffer
        ffptr          pointer of file to read
        len            maximum # characters to read
*/
char *get_a_line(char *where, int len, FILE *ffptr)
{
     char *retrn, *c;
     int j;

     while((retrn = fgets(where, (len - 1), ffptr)) != NULL) {
          where[(len - 1)] = NULLCHAR;
          if((c = strchr(where, '\n')) != NULL)
               *c = NULLCHAR;                     /* strip returns */
          if((c = strchr(where, ';')) != NULL)
               *c = NULLCHAR;                    /* strip comments */
          if((j = strspn(where, " \t")) == strlen(where))
               continue;                       /* drop blank lines */
          if(strlen(where))
               return(retrn);                /* return good string */
     }

     return(retrn);
}

/* check string for legal characters and make an
   integer from it.
   parameters are:
        stuff          string to be converted.
        err_mess       error message in case of illegal
                       characters in the string
        from           calling file name for errors.
        this           calling line number for errors.
*/
int make_int(char *stuff, char *err_mess, const char *from,
             int this)
{
     int k;

     if(strspn(stuff, "0123456789") != strlen(stuff)) {
          sprintf(mess, err_mess, stuff);
          my_fatal(from, this, mess);
          exit(1);
     }
     k = atoi(stuff);
     return(k);
}
