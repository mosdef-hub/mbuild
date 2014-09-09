/* Topology Builder        Bruce D. Ray
   IUPUI Physics Dept. NMR Center
   402 N. Blackford St.
   Indianapolis, IN   46202
   USA

   utility functions header        utilities.h
*/

#define NULLCHAR 0
#define MAXTOKEN 128
#define MAXCHAR 256

extern char mess[MAXCHAR];                         /* for error messages */
extern char *get_name(char *this, char *place,
                      int length);                 /* Get a name and return  */
                                                   /* progress in the string */
                                                   /* from which it came.    */

extern int make_tokens(char **list, char *line,
                       int count);                 /* Separate a string into */
                                                   /* tokens.                */
extern int caseless_strncmp(const char *string1,
                            const char *string2,
                            int length);           /* Case insensitive       */
                                                   /* string comparison.     */

extern int spaceline(char *line);                  /* Check line to be sure  */
                                                   /* it contains graphical  */
                                                   /* characters             */

extern void uppercase(char *inname);               /* Convert string to      */
                                                   /* upper case letters.    */
extern double CubeRoot(double x);                  /* Fast cube root         */
                                                   /* implementation         */
extern double fast_InvSqrt(double x);              /* Double precision fast  */
                                                   /* inverse square root.   */
extern char *get_a_line(char *where, int len,
                        FILE *ffptr);              /* Read lines, strip      */
                                                   /* comments and end of    */
                                                   /* line, and return non-  */
                                                   /* blank lines.           */
extern int make_int(char *stuff, char *err_mess,
                    const char *from, int this);   /* check string for legal */
                                                   /* characters and make an */
                                                   /* integer from it.       */
