#ifndef TYPE_H
#define TYPE_H 1

/** \typedef frame
 * A partition, usable at the limit of maxdim -1, since A[maxdim] is sometimes used specifically
 */
typedef struct
{
  int A[maxdim + 1];
// added by FB&FT 20060301 to take care of zeros...
  /** partition length */
  unsigned short int length; 
} frame;

/** \typedef qframe */
typedef struct
{
  frame A[maxdimqframe - 1 + 1];
} qframe;

/** \typedef bframe
 * This structure is used for bignumbers in Binary Coded Decimal (2 decimal digits by table entry) 
 * or as prime factors (table entry i is the power of the ith prime).
 */
typedef struct
{
  int A[maxl + 1];
} bframe;			


/** \typedef lbframe
 * struct used for polynoms
 */
typedef struct
{
  int A[maxdim + 1];
} lbframe;

/** \typedef barr*/
typedef struct
{
  int A[maxsize + 1];
} barr;

typedef frame *arrayptr;

/** \typedef smalltab*/
typedef struct
{
  char A[brtabmax - 1 + 1];
} smalltab;

/** \typedef string0
 * char strings used for example for file names and commands input
 * These strings are special : right filled with white spaces.
 */
typedef struct
{
  char A[bcol - 1 + 1];
} string0;

/** \typedef hktype */
typedef enum
{ rrow, ccol, rc } hktype;

/** \typedef caystype */
typedef enum
{ norm, pe6, e6bra, prod, triple } caystype;

/** \typedef grptype */
typedef enum
{ sung, un, son, on, spn, sn, an, g2, f4, en, e6, e7, e8, ospnm, unm, nill,
    l168, sunm, spnc, mp, sonc, unc } grptype;

/** \typedef modes 
 * the different modes available : sfn, branch mode, rep mode, DPrep mode*/
typedef enum
{ sfnm, brm, repm, dpm } modes;

/** \typedef term
 * used for sfn mode, a chained list of partitions with multiplicities
 */
typedef struct term
{				
  struct term *next;		// must be in first place
  /** partition multiplicity */
  int mult;
  /** partition by itself */
  frame val;
  /** must be ' ' (for scalars) or '#'*/
  char slab;	
} term;

/** \typedef termptr */
typedef term *termptr;

/** \typedef termarray  an array of partition lists*/
typedef struct
{
  termptr A[maxdim + 1];
} termarray;

/** \typedef ocharac
 * struct used for repmode (and encapsulated by dpmode)
 */
typedef struct ocharac
{				// was S63
  struct ocharac *next;		// must be in first place
  int mult;
  frame val;
  /** can be ' ', '+', '-' or '#' */
  char lab;		
  bool spin;
  /** if C6_double is true conval and conlab are meaningful */
  bool C6_double;		
  frame conval;
  char conlab;
} ocharac;
typedef ocharac *ocharptr;

/** \typedef groop */
typedef struct groop
{		
  grptype name;
  int rank;
  int rank2;
} groop;

/** \typedef groupArray*/
typedef struct
{
  groop A[maxprod - 1 + 1];
} groopArray;

/** \typedef prodrec
 * used for dpmode, a chained list of ocharptr, \see ocharac
 */
typedef struct prodrec
{				// was S64
  struct prodrec *next;		// must be in first place
  struct
  {
    ocharptr A[maxprod - 1 + 1];
  } prods;
  int mult;
} prodrec;
/** \typedef prodtype */
typedef prodrec *prodtype;

/*** \typedef fns a function is a chained list of strings (its lines) */
typedef struct fns
{
  string0 gbuff;
  struct fns *next;
} fns;

/** \typedef fnptrs list of functions */
typedef fns *fnptrs;

/** \typedef branchtab */
typedef struct branchtab
{
  smalltab tab;
} branchtab;

/** \typedef tabptr list of branchtab */
typedef branchtab *tabptr;

/** \typedef byte a byte on 8 bits */
typedef unsigned char byte;

/** \typedef Tableau */
typedef struct
{
  struct
  {
    int A[worksize - 1 + 1];
  } A[maxdim + 1];
} Tableau;

/** \struct fn table of user-defined functions */
struct
{
  fnptrs A[maxfn - 1 + 1];
} fn;

/** \struct svar Variables in sfn mode */
struct
{
  termptr A[svarlimit - 1 + 1];
} svar;

/** \struct vari Variables in rep mode */
struct
{
  ocharptr A[rvarlimit - 1 + 1];
} vari;

/** \struct pvar Variables in DP mode */
struct
{
  prodtype A[pvarlimit - 1 + 1];
} pvar;

/** \struct preject memory management: keep an eye on what to clean after computation 
* of a command line \see stackhandle, collectgarbage */
struct
{
  prodtype A[prjctlimit - 1 + 1];
} preject;

/** \struct reject memory management rep mode */
struct
{
  ocharptr A[rjctlimit - 1 + 1];
} reject;

/** \struct reject memory management rep mode */
struct
{
  ocharptr A[rjctlimit - 1 + 1];
} reject2;

/** \struct sreject memory management sfn mode */
struct
{
  termptr A[srjctlimit - 1 + 1];
} sreject;

#include "sets_mgmt.h"
/** these sets are used to accelerate the checks on special characters*/
struct
{
  setword S[MAXSETSIZE+1];
} brackets, plus_minus_etc, numbers, numbersEtc, letterset;

#endif /* ! TYPE_H */
