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

#include "standard.h"
#include "define.h"
#include "dim.h"
#include "type.h"
#include "utils.h"

// for sslab
#include "var.h"

// for ldisp
#include "s1.h"

// for cslabel
#include "label.h"

// for outerskew
#include "s2.h"
#include "outerskew.h"
#include "skew.h"

bool
skewcompatable (frame a, frame b)
{
  int i;

  i = 1;
  while ((a.A[i] >= b.A[i]) && (b.A[i] != 0))
    i = i + 1;

  return (bool) (b.A[i] == 0);
}

termptr
skew (frame sfn1, frame sfn2)
{
  return skewx (1, sfn1, sfn2, maxdim);
}

termptr
lskew (termptr list1, termptr list2)
{
  termptr templist, sublist, skewlist, x;
  int multy;
  char ch;

  skewlist = NULL;
  ch = ' ';
  while (list1 != NULL)
    {
      templist = list2;
      while (templist != NULL)
	{
	  sublist = skew (list1->val, templist->val);
	  multy = list1->mult * templist->mult;
	  if (sslab && (group == on))
	    {
	      if (((templist->slab == '#') && (list1->slab == ' '))
		  || ((templist->slab == ' ') && (list1->slab == '#')))
		ch = '#';
	      else
		ch = ' ';
	      cslabel (sublist, ch);
	    }
	  x = sublist;
	  if (multy != 1)
	    while (x != NULL)
	      {
		x->mult = x->mult * multy;
		x = x->next;
	      }
	  add (&skewlist, &sublist);
	  templist = templist->next;
	}
      sort (&skewlist, true);
      list1 = list1->next;
    }
  return skewlist;
}

termptr
pskew (termptr point1, termptr point2)
{
  register termptr top;
  termptr skewlist;
  int multy;
  register termptr temp;

  if ((point1 == NULL) || (point2 == NULL))
    return NULL;
  else
    {
      skewlist = skew (point1->val, point2->val);
      top = skewlist;
      multy = point1->mult * point2->mult;
      if (multy != 1)
	while (skewlist != NULL)
	  {
	    temp = skewlist;
	    temp->mult  *= multy;
	    skewlist = temp->next;
	  }
    }
  return top;
}

termptr
allskews (frame shape, bool exclude)
{
  int edge, previous;
  termptr result, ptr;
  frame rho;

  rho = nolls;
  if (exclude && (shape.A[1] == 0))
    return NULL;
  else
    {
      rho = shape;
      edge = len (&rho);
      snu (&result);
      result->val = rho;
      result->mult = 1;
      ptr = result;
      while (edge > 0)
	{
	  previous = rho.A[edge] - 1;
	  rho.A[edge] = previous;
	  while (previous != 0)
	    {
	      edge = edge + 1;
	      previous = MIN (previous, shape.A[edge]);
	      rho.A[edge] = previous;
	    }
	  edge = edge - 1;
	  if (!(exclude && (edge == 0)))
	    {
	      snu (&ptr->next);
	      ptr = ptr->next;
	      ptr->val = rho;
	      ptr->mult = 1;
	    }
	}
      ptr->next = NULL;
      return result;
    }
}

termptr
skewx (int c, frame alpha, frame beta, int maxx)
{
  int la, lb;
  register int i;
  termptr result;
  if ((cffrmfrm (alpha, beta) == EQUAL) || (beta.A[1] == 0))
    {
      snu (&result);
      result->mult = c;
      result->next = NULL;
      if (sign == EQUAL)
	result->val = nolls;
      else
	result->val = alpha;
    }
  else
    {
      la = len (&alpha);
      lb = len (&beta);
      if (skewcompatable (alpha, beta))
	{
	  if (la == lb)
	    while (alpha.A[la] == beta.A[la])
	      la = la - 1;
	  else
	    for (i = lb; i <= la; i++)
	      {
		beta.A[i + 1] = 0;
	      }
	  result = outerskew (c, nolls, alpha, beta, 0, la, maxx);
	}
      else
	result = NULL;
    }
  return result;
}

termptr
nskew (termptr list, termptr sklist, int n)
{
  termptr newlist, temp;
  register int i;

  if (n <= 0)
    newlist = NULL;
  else
    newlist = list;
  for (i = 1; i <= n; i++)
    {
      temp = lskew (newlist, sklist);
      if (i > 1)
	ldisp (&newlist);
      newlist = temp;
    }
  return newlist;
}
