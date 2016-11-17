#ifndef DIM_H
#define DIM_H 1 
#include <stdlib.h>
/* define base for big integers */
# define base 100 
/* maximal length for partitions +1 */
# define maxdim 201
/* */
# define maxdimqframe 102

/* initial max size of bignums */
# define maxl 200

# define maxsize 500
# define nmax 150
/*# define max 100*/ 
/* logging option - please do not touch */
# define pcol 100
/* size of string0 - please use carrefuly */
# define bcol 200
/* number of user defined functions */
# define maxfn 50
/* max number of groups (DPmode) */
# define maxprod 10
/* don't modify please... Related to dat files */
# define brtabmax 23000
/* following values count the number of variables in SFN, REP and DPM modes */
# define svarlimit 50
# define rvarlimit 50
# define pvarlimit 50
/** max number of structures for the garbage collector (unsigned short)*/
# define srjctlimit 500
# define rjctlimit 500
# define prjctlimit 500
/** used for outer sfn product and size of tableau, was 50, extended to maxdim/2 for O_sfnProduct to work correctly */
# define worksize (maxdim/2)

# define skipbrackets true
/** Continuation character for input that consists in multiple lines */
# define cont '&'
# define delim ','
# define qt '\''
# define spc ' '

/** for comparison functions */
# define LESS -1
# define EQUAL 0
# define GREATER 1

# define zerochar '0'
/*# define vducol 80*/
# define vdurow 40
# define lastlines 5
/*# define maxelements maxdim*/

/** it has to be less than maxdim, at the begining of SCHUR it was 30 */
# define POLYMAXDEGREE 90

/** needed by all strings */
# define MAXSTRING	bcol

# define tlines     (getenv("LINES")==NULL ? 37 : strtol(getenv("LINES"), (char **)NULL, 10))
/** \def display limit - please do not touch */
//# define tcol 100
# define tcol     (getenv("COLUMNS")==NULL ? 79 : strtol(getenv("COLUMNS"), (char **)NULL, 10))

#endif /* ! DIM_H */
