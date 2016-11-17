#include <stdbool.h>
#include "standard.h"
#include "dim.h"
#include "type.h"
#include "var.h"
#include "utils.h"
#include "s1.h"
#include "s2.h"
#include "outerskew.h"

/** \typedef Tab_Segptr */
typedef struct {
        termptr A[SEGPTRMAX - 1 + 1];
} Tab_Segptr;

/* Beware : in the following, there is dirty tricks to be able to cope
 * with compilers that do not allow nested functions */

#define STRUCT_CONTENT frame rho;bool firsttime;int j;int i;Tableau tableau;Tab_Segptr segptrs

#ifndef HAVE_NESTED_FUNCTIONS
  /* this  must be exactly the same content than below */ 
  typedef struct {
	  STRUCT_CONTENT;
	  int coef;
	  frame frme;
	  frame offset;
	  int corelen;
	  int framelen;
	  int limitx;
  } MyEnv;

  /* static : functions local to this file */
  static int qcf (termptr);
  static void segsort (void);
  static void h_count (int );

  static MyEnv myEnv;
  # define MYENV(m)  myEnv.m
  # define MYSTATIC  static
  # include "outerskew_inc.c"
#endif // ! HAVE_NESTED_FUNCTIONS

termptr
outerskew (int coef, frame core, frame frme, frame offset, int corelen,
           int framelen, int limitx)
{
    register int i,j;

#ifdef HAVE_NESTED_FUNCTIONS
    // local functions' prototypes
    auto int qcf (termptr);
    auto void segsort (void);
    auto void h_count (int);

/* This is dirty, but when we have nested functions we do not want to go through a struct.
 * So this has to be the same as the above struct */
    STRUCT_CONTENT;

    # define MYENV(m) m
    # define MYSTATIC
    # include "outerskew_inc.c"
#else
    MYENV(coef)=coef;
    MYENV(frme)=frme;
    MYENV(offset)=offset;
    MYENV(corelen)=corelen;
    MYENV(framelen)=framelen;
    MYENV(limitx)=limitx;
#endif

    /* outerskew main part */
    if (MYENV(limitx) >= maxdim)
        MYENV(limitx) = maxdim - 1;
    if (MYENV(frme).A[1] <= maxdim) {
        for (i = 0; i < SEGPTRMAX; i++)
            MYENV(segptrs).A[i] = NULL;
        for (i = 1; i <= worksize; i++)
            MYENV(tableau).A[0].A[i - 1] = 0;
        for (i = 1; i <= MYENV(framelen); i++)
            for (j = 1; j <= (MYENV(offset).A[i] + 1); j++)
                MYENV(tableau).A[i].A[j - 1] = 0;
        MYENV(rho) = nolls;
        MYENV(firsttime) = true;
        MYENV(rho) = core;
        h_count (1);
        return MYENV(segptrs).A[SEGPTRMAX - 1];
    } else {
        error (BAD_OUTER_SKEW, 0);
        return NULL;
    }
}                               //outerskew
