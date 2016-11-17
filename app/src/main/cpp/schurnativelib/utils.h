# ifndef UTILS_H
# define UTILS_H 1

# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
//# include "var.h"
# include "mymalloc.h"
# include "sets_mgmt.h"

/* Error codes */
# define MISTAKE      		1
# define COMMAND_NOT_FOUND 	2
# define MISSING_COMMA 		3
# define MISSING_PARAMETER 	4
# define MISSING_INTEGER	5
# define TOO_MANY_PARTS		6
# define WRONG_SERIES		7
# define WRONG_DIGIT 		8
# define NOT_IMPLEMENTED        9
# define CANNOT_BE_ZERO		10
# define INVALID_SERIES		11
# define FILENAME_ERR		12
# define NO_LOG_FILE            13
# define BAD_CONJUGATE		14
# define BAD_OUTER_SKEW		15
# define NO_SUCH_FILE		16
# define BAD_FUNCTION		17
# define NOT_IN_FUNCTION	18
# define MISSING_QUOTE		19
# define GROUP_NOT_SET		20
# define ZERO_DIVISOR		21
# define BAD_IRREP		22
# define PART_TOO_BIG		23
# define LINE_TOO_LONG 		24

# define WARNING_LIMIT		30

# define PARENTHESIS_NOT_NEEDED 31
# define ARE_YOU_SURE		32

# define boolstr(b) ( b ? "true" : "false" )

/* string functions */
#ifndef __GNUC__
char *strdup(char *); 
#endif
char locase(char );
void strcpylowercase(char *, char *);
void strcpyuppercase(char *, char*);
int charval(char);
char boolchar(bool);
bool interp(char *, char *, int);
char *stripwhite(char *);
char *findLastSpace (char *);

void print(char *, ...);

int sgn(int);

void progress(void);

int minusoneto(int);

void sound(int);
  
void inform(char *, char);
void warn(char *, char);
void error(int, int);
void Caseerror (int);
void aaargh(char *, bool);

unsigned short len (frame *);
unsigned sumOfPartition(frame *, unsigned);
unsigned maxOfPartition (frame *);
void free_all(void);
#endif /* ! UTILS_H */

