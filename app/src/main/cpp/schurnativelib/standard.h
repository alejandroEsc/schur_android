# ifndef STANDARD_H
# define STANDARD_H 1

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdbool.h>

typedef struct
{
  FILE *fp;
  unsigned eoln:1, eof:1, out:1, init:1,:12;
  char buf;
} text;

/*
**	Definitions for standard types
*/

/*
**  	Definitions for strings
*/

extern void Getl(text *);

# define Trunc(m)	((int)(m))
# define UpperCase(m)   ((m>='a' && m <='z') ? (m-'a'+'A') : m)

/* This is used when gcc complains too much about unused parameter */
#define UNUSED(x) (void)(x)

#endif /* !STANDARD_H */
