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

/* \file s.c
 */

/*
**	Definitions for i/o
*/
# include <stdio.h>
# include <math.h>
# include <limits.h>

# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "s2.h"
# include "g.h"
# include "skew.h"
# include "s.h"

# define Line __LINE__
prodtype discardedprds, pcurrent;

termptr
spseries (int n, int k)
{
  register termptr R165;
  termptr list, temp;

  temp = seriesx ('f', setlimit, full);
  schur_restrict (&temp, k, 'l');
  list = lconjgte (temp);
  schur_restrict (&list, n, 'f');
  schur_restrict (&list, k, 'z');
  temp = lconjgte (list);
  R165 = temp;
  return R165;
}

termptr
rseries (int n, char tag)
{
  register termptr R166;
  termptr temp;

  temp = seriesx (tag, setlimit, full);
  schur_restrict (&temp, n, 'l');
  R166 = temp;
  return R166;
}

termptr
makelist (int n)
{
  register termptr R167;
  termptr result, ptr;
  frame rho, nrho;
  int a, m, p, x;
  register int i;

  result = NULL;
  if ((n >= 0))
    {
      a = 1;
      if ((n == 0))
	{
	  snu (&ptr);
	  ptr->mult = 1;
	  ptr->val = nolls;
	  add (&result, &ptr);
	  dispsfn (&ptr);
	}
      else if ((n == 2))
	{
	  snu (&ptr);
	  ptr->mult = -1;
	  ptr->val = nolls;
	  add (&result, &ptr);
	  dispsfn (&ptr);
	  snu (&ptr);
	  ptr->mult = -1;
	  ptr->val = nolls;
	  ptr->val.A[1] = 1;
	  add (&result, &ptr);
	  dispsfn (&ptr);
	  snu (&ptr);
	  ptr->mult = 1;
	  ptr->val = nolls;
	  ptr->val.A[1] = 2;
	  add (&result, &ptr);
	  dispsfn (&ptr);
	}
      else if ((n == 1))
	{
	  snu (&ptr);
	  ptr->mult = 1;
	  ptr->val = nolls;
	  ptr->val.A[1] = 1;
	  add (&result, &ptr);
	  dispsfn (&ptr);
	}
      else
	while ((a <= n))
	  {
	    rho = nolls;
	    p = (n - a) / 2;
	    if ((bool) ((n) & 1))
	      x = (n + a + 2) / 2;
	    else
	      x = (n + a) / 2;
	    if ((bool) ((x) & 1))
	      m = -1;
	    else
	      m = 1;
	    rho.A[1] = a;
	    if ((p >= 1))
	      for (i = 2; i <= p + 1; i++)
		{
		  rho.A[i] = 1;
		}
	    snu (&ptr);
	    ptr->mult = m;
	    ptr->val = rho;
	    add (&result, &ptr);
	    ldisp (&ptr);
	    if ((p >= 1))
	      {
		nrho = rho;
		nrho.A[p + 1] = 0;
		snu (&ptr);
		ptr->mult = m;
		ptr->val = nrho;
		add (&result, &ptr);
		ldisp (&ptr);
	      }
	    a = a + 1;
	  }
      sort (&result, true);
    }
  else
    (void) fprintf (output.fp, "ERROR: n must be a non-negative int\n"),
      Putl (output, 1);
  R167 = result;
  return R167;
}

termptr
allseriesbut_t (char series, int size, frame shape, unsigned maxlng)
{
  register termptr R168;
  bool first, sb_conj, strip, self, notodd, shaped, oddparts, evenparts, gap;
  int r, edge, last, left=0, take, wt, coe;
  register int i;
  frame rho, part, help;
  termptr listhead, listtail;

  listhead = NULL;
  listtail = NULL;
  //for (i = len (shape) + 1; i <= maxdim; i++)
  for (i = MAX(shape.length,len(&shape)) + 1; i <= maxdim; i++)
    {
      shape.A[i] = 0;
    }
  help=nolls;
  part=nolls;
  rho=nolls;
  if (((series == 'b') || (series == 'c') || (series == 'l')
       || (series == 'q') || (series == 'v') || (series == 'x')))
    sb_conj = true;
  else
    sb_conj = false;

  if (((series == 'a') || (series == 'c')))
    strip = true;
  else
    strip = false;
  if (((series == 'a') || (series == 'c') || (series == 'e')
       || (series == 'g')))
    self = true;
  else
    self = false;
  if (((series == 'v') || (series == 'w') || (series == 'x')
       || (series == 'y')))
    notodd = true;
  else
    notodd = false;
  gap = self;
  shaped = (bool) (size < 0);
  oddparts = (bool) (self && !strip);
  evenparts = strip;
  if ((series == 'b') || (series == 'd'))
    evenparts = true;
  if ((size == 2) && oddparts)
    size = 1;
  r = frank (shape);
  if (self && shaped && sb_conj && strip)
    for (i = 1; i <= r; i++)
      {
	shape.A[i] = shape.A[i] - 1;
      }
  if (sb_conj && shaped)
    conjgte (&shape);
  if (self && shaped && (!sb_conj) && strip)
    {
      conjgte (&shape);
      for (i = 1; i <= r; i++)
	{
	  shape.A[i] = shape.A[i] - 1;
	}
      conjgte (&shape);
    }
  if (((series == 'b') || (series == 'd')))
    if (shaped)
      {
	int lastStep = len (&shape);	//protect termination
	for (i = 1; i <= lastStep; i++)
	  {
	    shape.A[i] = shape.A[i] - (unsigned) ((bool) ((shape.A[i]) & 1));
	  }
      }
    else
      size = size - (unsigned) ((bool) ((size) & 1));
  if (((series == 'r') || (series == 's'))) {
    int lastStep = len (&shape);     //protect termination
    for (i = 2; i <= lastStep; i++)
      {
	shape.A[i] = 1;
      }
  }
  if (notodd)
    for (i = 3; i <= maxdim; i++)
      {
	shape.A[i] = 0;
      }

  if (((series == 'l') || (series == 'm') || (series == 'p')
       || (series == 'q')))
    for (i = 2; i <= maxdim; i++)
      {
	shape.A[i] = 0;
      }
  r = frank (shape);
  if (self && shaped)
    {
      part = shape;
      conjgte (&part);
      for (i = 1; i <= r; i++)
	{
	  shape.A[i] =
	    ( MIN (shape.A[i], part.A[i]) - i ) * 2 + 1 + (unsigned) (strip);
	}
      for (i = r + 1; i <= maxdim; i++)
	{
	  shape.A[i] = 0;
	}
      shape.length=r;
    }
  first = true;
  edge = 1;
  while (edge > 0)
    {
      if (first)
	{
	  if (shaped)
	    part = shape;
	  else if (size == 0)
	    part.A[1] = 0;
	  else
	    {
	      part.A[3] = 0;
	      part.A[1] = size - (unsigned) ((bool) ((size) & 1)
					     && evenparts) -
		(unsigned) (!(bool) ((size) & 1) && oddparts);
	      if (!(bool) ((size) & 1) && oddparts)
		part.A[2] = 1;
	      else
		part.A[2] = 0;
	    }
	  first = false;
	  edge = len (&part);
	  left = 0;
	}
      else
	{
	  take = 1 + (int) ((part.A[edge] > 1) && (oddparts || evenparts));
	  last = part.A[edge] - take;
	  part.A[edge] = last;
	  left = left + take;
	  while (last != 0)
	    {
	      edge = edge + 1;
	      take = (int) (!(bool) ((left) & 1) && oddparts) +
			(int) ((bool) ((left) & 1) && evenparts);
	      if (left == 0)
		take = 0;
	      last =
		MIN (shape.A[edge],
		     MIN (last -
			  (int) ((oddparts || evenparts) && gap) * ((int) (last > 1) + 1),
			  left - take));
	      part.A[edge] = last;
	      left = left - last;
	    }
	  edge = edge - 1;
	}
      if (self)
	{
	  for (i = 1; i <= edge; i++)
	    {
	      help.A[i] = part.A[i] / 2 + i;
	    }
	  help.A[edge + 1] = 0;
	  help.length=edge;
	  rho = help;
	  conjgte (&rho); 
	  for (i = 1; i <= edge; i++)
	    {
	      rho.A[i] = help.A[i];
	    }
	}
      else
	rho = part;
      if (strip)
	for (i = 1; i <= edge; i++)
	  {
	    rho.A[i] = rho.A[i] - 1;
	  }
      if (sb_conj)
	conjgte (&rho);
      wt = wtfrm (&rho);
      if (!((bool) ((wt) & 1) && notodd))
	{
	  coe = 1;
	  if (strip)
	    coe = minusoneto (wt / 2);
	  if (((series == 'h') || (series == 'l') || (series == 'p')))
	    coe = minusoneto (wt);
	  if (series == 'e')
	    coe = minusoneto ((wt + edge) / 2);
	  if (series == 'g')
	    coe = minusoneto ((wt - edge) / 2);
	  if (notodd && ((series != 'x') && (series != 'y')))
	    coe = minusoneto (part.A[2]);
	  if ((series == 'r') && (wt != 0))
	    coe = 2 * minusoneto (wt);
	  if ((series == 's') && (wt != 0))
	    coe = 2;
	  if (maxlng==0  || rho.length<=maxlng)
	    insort (coe, rho, &listhead, &listtail);
	  progress ();
	}
    }
  R168 = listhead;
  return R168;
}

termptr
onlyseries_t (int size, frame shape)
{
  register termptr R169;
  int top;
  register int i;
  frame rho;
  termptr listhead, listtail;

  listhead = NULL;
  listtail = NULL;
  for (i = 0; i <= maxdim; i++)
    {
      rho.A[i] = 0;
    }
  top = INT_MAX;
  if (size >= 0)
    top = Trunc (sqrt (2 * size + 0.25) - 0.5);
  else
    for (i = 0; i <= len (&shape) + 1; i++)
      {
	top = MIN (top, shape.A[i + 1] + i);
      }
  for (i = 0; i <= top; i++)
    {
      rho.A[i + 1] = top - i;
    }
  insort (1, rho, &listhead, &listtail);
  while (rho.A[1] != 0)
    {
      for (i = 2; i <= rho.A[1]; i++)
	{
	  rho.A[i] = rho.A[i] - 1;
	}
      rho.A[1] -= 1;
      insort (1, rho, &listhead, &listtail);
      progress ();
    }
  R169 = listhead;
  return R169;
}

termptr
seriesx (char ser, int size, frame shape)
{
  register termptr result;

  if (ser == 't')
    result = onlyseries_t (size, shape);
  else
    result = allseriesbut_t (ser, size, shape, 0);
  return result;
}

bool
validser (char ser)
{
  register bool R171 = false;	//modified by FB, was uninitialized

  switch ((int) (locase (ser)))
    {
    case 'a':
    case 'b':
    case 'c':
    case 'd':
    case 'e':
    case 'f':
    case 'g':
    case 'h':
    case 'l':
    case 'm':
    case 'p':
    case 'q':
    case 'r':
    case 's':
    case 't':
    case 'v':
    case 'w':
    case 'x':
    case 'y':
      R171 = true;
      break;
    case 'i':
    case 'j':
    case 'k':
    case 'n':
    case 'o':
    case 'u':
    case 'z':
      R171 = false;
      break;
    default:
      print ("ERROR: Series unknown or syntax error\n");
    }
  return R171;
}

termptr
useseries (char ser, termptr list, bool byweight, bool doskew, int p)
{
  register termptr R172;
  termptr result, sublist1, sublist2, locall;

  result = NULL;
  if (validser (ser))
    {
      snu (&locall);
      locall->next = NULL;
      while (list != NULL)
	{
	  if (byweight)
	    sublist1 = seriesx (ser, wtfrm (&list->val), full);
	  else
	    sublist1 = seriesx (ser, -1, list->val);
	  if (doskew)
	    {
	      locall->val = list->val;
	      locall->mult = list->mult;
	      sublist2 = lskew (locall, sublist1);
	      ldisp (&sublist1);
	      sublist1 = sublist2;
	    }
	  merge (&result, &sublist1, true, doskew);
	  list = list->next;
	}
      dispsfn (&locall);
    }
  else
    error (INVALID_SERIES, p);
  R172 = result;
  return R172;
}

void
cleave (termptr * list)
{
  termptr ptr1, ptr2;
  int s, n;
  register int k;

  ptr1 = (*list);
  while (ptr1 != NULL)
    {
      s = ptr1->mult;
      n = abs (s);
      s = sgn (s);
      if (n != 1)
	{
	  ptr1->mult = s;
	  for (k = 1; k <= n - 1; k++)
	    {
	      snu (&ptr2);
	      (*ptr2) = (*ptr1);
	      ptr1->next = ptr2;
	      ptr1 = ptr2;
	    }
	}
      ptr1 = ptr1->next;
    }
}

termptr
seriesm (int wt, bool exclude)
{
  register termptr R173;
  termptr result, ptr;
  frame rho;

  rho = nolls;
  rho.A[1] = wt;
  result = NULL;
  if (!((wt == 0) && exclude))
    {
      snu (&result);
      result->val = rho;
      result->mult = 1;
      ptr = result;
      while (rho.A[1] != (int) (exclude))
	{
	  rho.A[1] = rho.A[1] - 1;
	  snu (&ptr->next);
	  ptr = ptr->next;
	  ptr->val = rho;
	  ptr->mult = 1;
	}
      dispsfn (&ptr->next);
      ptr->next = NULL;
    }
  R173 = result;
  return R173;
}


termptr
rhook (int n)
{
  register termptr R174;
  termptr newlist, temp, temp1;
  int p, x, d;
  register int i;
  register int a;

  newlist = NULL;
  if ((bool) ((n) & 1))
    d = 2;
  else
    d = 0;
  for (a = 1; a <= n; a++)
    {
      p = (n + a + d) / 2;
      snu (&temp);
      if ((bool) ((p) & 1))
	temp->mult = -1;
      else
	temp->mult = 1;
      temp->val = nolls;
      temp->val.A[1] = a;
      x = (n - a) / 2;
      if (x == 0)
	add (&newlist, &temp);
      if (x >= 1)
	{
	  for (i = 1; i <= x; i++)
	    {
	      temp->val.A[i + 1] = 1;
	    }
	  snu (&temp1);
	  temp1->val = temp->val;
	  temp1->mult = temp->mult;
	  temp1->val.A[x + 1] = 0;
	  temp1->next = NULL;
	  add (&newlist, &temp);
	  if (temp1 != NULL)
	    {
	      add (&newlist, &temp1);
	      dispsfn (&temp1);
	    }
	}
      dispsfn (&temp);
    }
  sort (&newlist, true);
  R174 = newlist;
  return R174;
}
