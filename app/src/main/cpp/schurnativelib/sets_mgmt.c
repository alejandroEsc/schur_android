/*
 * This file is part of SCHUR.
 *
 * SCHUR - an interactive program for calculating properties of Lie
 * groups and symmetric functions.
 * Copyright (C) 2006  Franck BUTELLE, Frédéric Toumazet
 * 
 * SCHUR is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * SCHUR is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with SCHUR; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "sets_mgmt.h"

setptr
Currset (int n, setptr sp)
{
  static setword Space[SETSPACE];
  static setptr Top = Space;

  switch (n)
    {
    case CLAIMSET:
      Top = Space;
      return (0);
    case NEWTMPSET:
      if (Space + SETSPACE - Top <= MAXSETSIZE)
	{
	  fprintf (stderr, "Set-space exhausted\n");
	  exit (1);
	}
      Top[0] = 0; // initialize an empty set
      return (Top);
    case SAVESET: // the set is already in Space.
	Top = sp + sp[0] + 1;
      return (sp);
    }
  return (0);
}

setptr
Union (register setptr p1, register setptr p2)
{
  register int ig, jg, k;
  register setptr sp = NewTmpSet (), p3 = sp;

  jg = *p1;
  *p3 = jg;
  if (jg > *p2)
    jg = *p2;
  else
    *p3 = *p2;
  k = *p1 - *p2;
  p1++, p2++, p3++;
  for (ig = 0; ig < jg; ig++)
    *p3++ = (*p1++ | *p2++);
  while (k > 0)
    {
      *p3++ = *p1++;
      k--;
    }
  while (k < 0)
    {
      *p3++ = *p2++;
      k++;
    }
  return (Saveset (sp));
}

bool
Member (register unsigned int m, register setptr sp)
{
  register unsigned int ig = m / NBSETBITS + 1;

  return ( (ig <= *sp) && (sp[ig] & (1 << (m % NBSETBITS))) );
}

void
addtoSet (unsigned int m, register setptr sp)
{
  register unsigned int ig = m / NBSETBITS + 1;

  sp[ig] |=  1 << (m % NBSETBITS);
}

/* copy the N higher bytes of S2 in S1 */
/*
void
Setncpy (register setptr S1, register setptr S2, register unsigned int N)
{
  register unsigned int m;

  N /= NBSETBITS;
  *S1++ = --N;
  m = *S2++;
  while (m != 0 && N != 0)
    {
      *S1++ = *S2++;
      --N;
      --m;
    }
  while (N-- != 0)
    *S1++ = 0;
} */

//This function is defined mainly for debugging existing pre-defined sets.
void
displayCharsOfSet (register setptr S)
{
  register unsigned i;
  setword n = *S;

  for (i = 0; i <= n*NBSETBITS-1; i++)
     if (Member(i,S)) {
	if (i >= ' ' && i < 127)	// a character that can be displayed
	  printf ("%c", i);
	else
	  printf ("[%d]", i);
     }
  printf(".\n");
}

// S must be already initialized.
void
strtoSet (char *src, register setptr S)
{
	register unsigned i;

	S[0]=MAXSETSIZE;
	for (i=1; i< MAXSETSIZE;i++)
		S[i]=0;

	for (i=0 ; i <strlen(src) ; i++)
		addtoSet(src[i], S);
}
