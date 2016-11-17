# ifndef SETS_MGMT_H
# define SETS_MGMT_H 1

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
  
# ifndef SETSPACE
// Max. memory allocated to sets management
# define  SETSPACE 256
# endif

// Max. number of words used in sets. For now on, 9 is a minimum.
#define MAXSETSIZE 9

/*
**	Definitions for set-operations
*/

/* An integer set is encoded as a chained list of unsigned short integer
 * first value is the length of the list in number of unsigned short integers
 * next values are "setword" a binary representation of this part of the set.
 * For example, a set consisting of the 2 elements 32 and 33 is encoded as 3, 0x0000, 0x0000, 0x0003
 */
# define  CLAIMSET     0
# define  NEWTMPSET       1
# define  SAVESET      2
# define  Claimset() 	(void)Currset(CLAIMSET, (setptr)0)
# define  NewTmpSet() 	Currset(NEWTMPSET, (setptr)0)
# define  Saveset(s) 	Currset(SAVESET, s)

typedef unsigned short setword;
# define  NBSETBITS	(8*sizeof(setword))

typedef setword *setptr;

setptr Currset (int, setptr);
setptr Union (setptr , setptr);
bool Member (unsigned int, setptr);
//void Setncpy (setptr, setptr, unsigned int);
void displayCharsOfSet (setptr);
void addtoSet (unsigned int, setptr);
void strtoSet (char *, setptr);
# endif /* ! SETS_MGMT_H */
