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

/*
**	Definitions for i/o
*/
# include <stdio.h>
# include "standard.h"
# include "define.h"

/*
**	Start of program definitions
*/
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "skew.h"
# include "s1.h"

void
qup (termptr * list)
{
  termptr newlist, temp;
  bool startx;
  int dummy;
  register int i;

  newlist = NULL;
  while ((*list) != NULL)
    {
      register termptr W1 = *list;

      do
	{
	  startx = true;
	  for (i = 1; i <= maxdim - 1; i++)
	    {
	      if (abs (W1->val.A[i]) < abs (W1->val.A[i + 1]))
		{
		  startx = false;
		  dummy = W1->val.A[i];
		  W1->val.A[i] = W1->val.A[i + 1];
		  W1->val.A[i + 1] = dummy;
		  W1->mult = -W1->mult;
		}
	    }
	}
      while (startx == false);	/*8/6/99 */
      temp = W1->next;
      W1->next = newlist;
      newlist = (*list);
      (*list) = temp;
    }
  (*list) = newlist;
}

int
qlen (frame f)
{
  register int i;

  i = maxdim - 1;
  do
    {
      if (f.A[i] == 0)
	i = i - 1;
    }
  while (!((f.A[i] != 0) || (i == 1)));
  if (f.A[i] == 0)
    return (0);
  return i;
}

void
Qsort (termptr * head, bool reverse)
{
  termptr lastptr, sortedlist, test, temp, tail;
  int ord;
  bool inserted;

  progress ();
  if (reverse)
    {
      lastptr = NULL;
      while ((*head) != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = (*head);
	  (*head) = tail;
	}
      (*head) = lastptr;
    }
  snu (&tail);
  tail->next = NULL;
  snu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while ((*head) != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {
	      ord = qtestord (&test->val, &temp->val);	/*q */
	      switch ((int) (ord))
		{
		case -1:
		  lastptr = temp;
		  temp = temp->next;
		  break;
		case 0:
		  temp->mult = temp->mult + test->mult;
		  dispsfn (&test);
		  inserted = true;
		  break;
		case 1:
		  lastptr->next = test;
		  test->next = temp;
		  inserted = true;
		  break;
		default:
		  Caseerror (Line);
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispsfn (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispsfn (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispsfn (&sortedlist);
}

void
sort (termptr * head, bool reverse)
{
  termptr lastptr, sortedlist, test, temp, tail;
  int w1, w2, ord;
  bool inserted;

  progress ();
  if (reverse)
    {
      lastptr = NULL;
      while ((*head) != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = (*head);
	  (*head) = tail;
	}
      (*head) = lastptr;
    }
  snu (&tail);
  tail->next = NULL;
  snu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while ((*head) != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {
	      ord = qtestord (&test->val, &temp->val);
	      if (weight)
		{
		  w1 = wtfrm (&test->val);
		  w2 = wtfrm (&temp->val);
		  if ((w1 == w2))
		    {
		      switch ((int) (ord))
			{
			case -1:
			  lastptr = temp;
			  temp = temp->next;
			  break;
			case 0:
			  temp->mult = temp->mult + test->mult;
			  dispsfn (&test);
			  inserted = true;
			  break;
			case 1:
			  lastptr->next = test;
			  test->next = temp;
			  inserted = true;
			  break;
			default:
			  Caseerror (Line);
			}
		    }
		  else if ((w1 > w2))
		    {
		      lastptr->next = test;
		      test->next = temp;
		      inserted = true;
		    }
		  else if ((w1 < w2))
		    {
		      lastptr = temp;
		      temp = temp->next;
		    }
		}
	      else
		{
		  switch ((int) (ord))
		    {
		    case -1:
		      lastptr = temp;
		      temp = temp->next;
		      break;
		    case 0:
		      temp->mult = temp->mult + test->mult;
		      dispsfn (&test);
		      inserted = true;
		      break;
		    case 1:
		      lastptr->next = test;
		      test->next = temp;
		      inserted = true;
		      break;
		    default:
		      Caseerror (Line);
		    }
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispsfn (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispsfn (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispsfn (&sortedlist);
}

void
osort (ocharptr * head, bool reverse)
{
  ocharptr lastptr, sortedlist, test, temp, tail;
  bool inserted;
  int ord;

  progress ();

  lastptr = NULL;
  if (reverse)			// if reverse is true, reverse the list to start from the end
    {
      while ((*head) != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = (*head);
	  (*head) = tail;
	}
      (*head) = lastptr;
    }

  cnu (&tail);
  tail->next = NULL;
  cnu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while ((*head) != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {

	      ord = otestord (test, temp);
	      switch ((int) (ord))
		{
		case -1:
		  lastptr = temp;
		  temp = temp->next;
		  break;
		case 0:
		  temp->mult = temp->mult + test->mult;
		  dispchr (&test);
		  inserted = true;
		  break;
		case 1:
		  lastptr->next = test;
		  test->next = temp;
		  inserted = true;
		  break;
		default:
		  Caseerror (Line);
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispchr (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispchr (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispchr (&sortedlist);
}


void
schur_psort (prodtype * head, bool reverse)
{
  prodtype lastptr, sortedlist, test, temp, tail;
  int ord;
  bool inserted;
  progress ();
  (*head) = pcheck ((*head));
  temp = (*head);
  (*head) = prodexpand ((*head));
  pdisp (&temp);
  if (reverse)
    {
      lastptr = NULL;
      while ((*head) != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = (*head);
	  (*head) = tail;
	}
      (*head) = lastptr;
    }
  pnu (&tail);
  tail->next = NULL;
  pnu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while ((*head) != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {
	      ord = ptestord (test, temp);
	      switch ((int) (ord))
		{
		case -1:
		  lastptr = temp;
		  temp = temp->next;
		  break;
		case 0:
		  temp->mult = temp->mult + test->mult;
		  dispprod (&test);
		  inserted = true;
		  break;
		case 1:
		  lastptr->next = test;
		  test->next = temp;
		  inserted = true;
		  break;
		default:
		  Caseerror (Line);
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispprod (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispprod (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispprod (&sortedlist);

}

void
gsort (ocharptr * head, bool reverse, caystype cays)
{
  ocharptr lastptr, sortedlist, test, temp, tail;
  bool inserted;
  int odr;

  lastptr = NULL;
  if (reverse)
    {
      while ((*head) != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = (*head);
	  (*head) = tail;
	}
      (*head) = lastptr;
    }
  cnu (&tail);
  tail->next = NULL;
  cnu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while ((*head) != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {
	      odr = testord (&test->val, &temp->val);
	      if (odr == EQUAL)
		{
		  if (test->C6_double)
		    if (temp->C6_double)
		      {
			odr = testord (&test->conval, &temp->conval);
			if ((odr == EQUAL) && (cays == prod))
			  {
			    if ((unsigned) (test->conlab) >
				(unsigned) (temp->conlab))
			      odr = GREATER;
			    else
			      if ((unsigned) (test->conlab) <
				  (unsigned) (temp->conlab))
			      odr = LESS;
			  }
		      }
		    else
		      odr = GREATER;
		  else if (temp->C6_double)
		    odr = LESS;
		}
	      if (odr == EQUAL)
		{
		  if (test->spin != temp->spin)
		    if (temp->spin)
		      odr = GREATER;
		    else
		      odr = LESS;
		  else if (test->lab != temp->lab)
		    {
		      if (test->lab > temp->lab)
			odr = GREATER;
		      else
			odr = LESS;
		    }
		}
	      switch ((int) (odr))
		{
		case -1:
		  lastptr = temp;
		  temp = temp->next;
		  break;
		case 0:
		  temp->mult = temp->mult + test->mult;
		  dispchr (&test);
		  inserted = true;
		  break;
		case 1:
		  lastptr->next = test;
		  test->next = temp;
		  inserted = true;
		  break;
		default:
		  Caseerror (Line);
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispchr (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispchr (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispchr (&sortedlist);
}

void
merge (termptr * list1, termptr * list2, bool adding, bool inc)
{
  termptr ptr1, ptr2, ptr3;
  int coe;

  if (!adding && !inc)
    aaargh ("bad merge?", true);
  if (*list2 != NULL)
    if (*list1 == NULL)
      {
	*list1 = *list2;
	*list2 = NULL;
	if (!adding)
	  coeffset (list1, 0, '-');	// opposite 
      }
    else
      {
	ptr1 = *list1;
	ptr2 = *list1;
	while (*list2 != NULL)
	  {
	    while (cfptrptr (ptr1, *list2) == 1)
	      {
		ptr2 = ptr1;
		ptr1 = ptr1->next;
	      }
	    if (sign == 0)
	      {
		if (inc)
		  {
		    coe = ptr1->mult;
		    if (adding)
		      coe = coe + (*list2)->mult;
		    else
		      coe = coe - (*list2)->mult;
		    ptr1->mult = coe;
		    if (coe == 0)
		      if (ptr1 == *list1)
			{
			  *list1 = (*list1)->next;
			  dispsfn (&ptr1);
			  ptr1 = *list1;
			  ptr2 = *list1;
			}
		      else
			{
			  ptr2->next = ptr1->next;
			  dispsfn (&ptr1);
			  ptr1 = ptr2;
			}
		  }
		ptr3 = *list2;
		*list2 = (*list2)->next;
		dispsfn (&ptr3);
	      }
	    else
	      {
		if (!adding)
		  (*list2)->mult = -(*list2)->mult;
		ptr3 = *list2;
		*list2 = (*list2)->next;
		if (ptr1 == *list1)
		  {
		    *list1 = ptr3;
		    ptr3->next = ptr1;
		    ptr1 = *list1;
		    ptr2 = *list1;
		  }
		else
		  {
		    ptr2->next = ptr3;
		    ptr3->next = ptr1;
		    ptr1 = ptr3;
		    ptr2 = ptr3;
		  }
	      }
	  }
      }
}

void
cmerge (ocharptr * list1, ocharptr * list2, bool adding, bool inc)
{
  ocharptr ptr1, ptr2, ptr3;
  int coe;

  if ((!adding) && (!inc))
    aaargh ("bad merge?", true);
  if ((*list2) != NULL)
    if ((*list1) == NULL)
      {
	(*list1) = (*list2);
	(*list2) = NULL;
      }
    else
      {
	ptr1 = (*list1);
	ptr2 = (*list1);
	while ((*list2) != NULL)
	  {
	    while (ccfptrptr (ptr1, (*list2)) == 1)
	      {
		ptr2 = ptr1;
		ptr1 = ptr1->next;
	      }
	    if (sign == 0)
	      {
		if (inc)
		  {
		    coe = ptr1->mult;
		    if (adding)
		      coe = coe + (*list2)->mult;
		    else
		      coe = coe - (*list2)->mult;
		    ptr1->mult = coe;
		    if (coe == 0)
		      if (ptr1 == (*list1))
			{
			  (*list1) = (*list1)->next;
			  dispchr (&ptr1);
			  ptr1 = (*list1);
			  ptr2 = (*list1);
			}
		      else
			{
			  ptr2->next = ptr1->next;
			  dispchr (&ptr1);
			  ptr1 = ptr2;
			}
		  }
		ptr3 = (*list2);
		(*list2) = (*list2)->next;
		dispchr (&ptr3);
	      }
	    else
	      {
		if (!adding)
		  (*list2)->mult = -(*list2)->mult;
		ptr3 = (*list2);
		(*list2) = (*list2)->next;
		if (ptr1 == (*list1))
		  {
		    (*list1) = ptr3;
		    ptr3->next = ptr1;
		    ptr1 = (*list1);
		    ptr2 = (*list1);
		  }
		else
		  {
		    ptr2->next = ptr3;
		    ptr3->next = ptr1;
		    ptr1 = ptr3;
		    ptr2 = ptr3;
		  }
	      }
	  }
      }
}

void
insort (int c, frame f, termptr * head, termptr * tail)
{
  termptr ptr1, ptr2, ptr3;

  progress ();
  if (c != 0)
    if ((*head) == NULL)
      {
	snu (&(*head));
	(*head)->val = f;
	(*head)->mult = c;
	(*head)->next = NULL;
	(*tail) = (*head);
      }
    else if (cfptrfrm ((*tail), f) != -1)
      if (sign == 1)
	{
	  snu (&(*tail)->next);
	  (*tail) = (*tail)->next;
	  (*tail)->val = f;
	  (*tail)->mult = c;
	  (*tail)->next = NULL;
	}
      else
	{
	  (*tail)->mult = (*tail)->mult + c;
	  if ((*tail)->mult == 0)
	    {
	      ptr1 = (*head);
	      ptr2 = (*head);
	      dispsfn (&(*tail));
	      while (ptr1 != (*tail))
		{
		  ptr2 = ptr1;
		  ptr1 = ptr1->next;
		}
	      if (ptr1 != ptr2)
		{
		  ptr2->next = NULL;
		  (*tail) = ptr2;
		}
	      else
		{
		  (*head) = NULL;
		  (*tail) = NULL;
		}
	    }
	}
    else
      {
	ptr1 = (*head);
	ptr2 = (*head);
	while (cfptrfrm (ptr1, f) == 1)
	  {
	    ptr2 = ptr1;
	    ptr1 = ptr1->next;
	  }
	if (sign == 0)
	  {
	    ptr1->mult = ptr1->mult + c;
	    if (ptr1->mult == 0)
	      if (ptr1 != (*head))
		{
		  ptr2->next = ptr1->next;
		  dispsfn (&ptr1);
		  if (ptr2->next == NULL)
		    (*tail) = ptr2;
		}
	      else
		{
		  (*head) = (*head)->next;
		  dispsfn (&ptr1);
		  if ((*head) == NULL)
		    (*tail) = NULL;
		}
	  }
	else
	  {
	    snu (&ptr3);
	    ptr3->val = f;
	    ptr3->mult = c;
	    if (ptr1 == (*head))
	      {
		ptr3->next = ptr2;
		(*head) = ptr3;
	      }
	    else
	      {
		ptr3->next = ptr2->next;
		ptr2->next = ptr3;
	      }
	  }
      }
}

/** divide all multiplicities by c*/
void
pdiv (prodtype * list, int c)
{
  prodtype ptr;

  ptr = (*list);
  if (c == 0)
    print ("division by zero not allowed\n");
  else
    {
      while (ptr != NULL)
	{
	  if (ptr->mult % c == 0)
	    ptr->mult = ptr->mult / c;
	  else
	    print ("error: inappropriate coefficients\n");
	  ptr = ptr->next;
	}
    }
}

/** divide all multiplicities by c*/
void
rdiv (ocharptr * list, int c)
{
  ocharptr ptr;

  ptr = *list;
  if (c == 0)
    print ("division by zero not allowed\n");
  else
    {
      while (ptr != NULL)
	{
	  if (ptr->mult % c == 0)
	    ptr->mult = ptr->mult / c;
	  else
	    print ("error: inappropriate coefficients\n");
	  ptr = ptr->next;
	}
    }
}

/** memory allocation of an ocharac (used in rep mode) and initialize it */
void
cnu (ocharptr * chrc)
{
  *chrc = (ocharptr) malloc (sizeof (**chrc));
  (*chrc)->next = NULL;
  (*chrc)->conval = nolls;	// nolls is an array full of zeros
  (*chrc)->val = nolls;
  (*chrc)->lab = ' ';
  rcount = rcount + 1;
}

/** memory allocation of a sfn (used in sfn mode) and initialize it */
void
snu (termptr * sfn)
{
  *sfn = (termptr) malloc (sizeof (**sfn));
  (*sfn)->val = nolls;
  (*sfn)->slab = ' ';
  (*sfn)->next = NULL;
  count = count + 1;
}

/** memory allocation of a prodtype (used in DP mode) and initialize it */
void
pnu (prodtype * prd)
{
  register int j;

  pcount = pcount + 1;
  *prd = (prodtype) malloc ( sizeof (**prd) );
  for (j = 1; j <= maxprod; j++)
    {
      (*prd)->prods.A[j - 1] = NULL;
    }
  (*prd)->next = NULL;
}

/** free memory pointed to by sfn (a partition list)*/
void
ldisp (termptr * sfn)
{
  termptr temp;

  while (*sfn != NULL)
    {
      temp = *sfn;
      *sfn = temp->next;
      count = count - 1;
      free (temp);
    }
}

/** free memory pointed to by chrc (a list of ocharac) */
void
odisp (ocharptr * chrc)
{
  ocharptr temp;

  while (*chrc != NULL)
    {
      rcount = rcount - 1;
      temp = *chrc;
      *chrc = temp->next;
      free (temp);
    }
}

/** free memory pointed to by prd (a list of prodtype) */
void
pdisp (prodtype * prd)
{
  prodtype temp;
  register int j;

  while ((*prd) != NULL)
    {
      pcount = pcount - 1;
      temp = *prd;
      *prd = (*prd)->next;
      for (j = 1; j <= maxprod; j++)
	{
	  odisp (&temp->prods.A[j - 1]);
	}
      free (temp);
    }
}

/** free memory pointed to by first element of chrc, dangerous !
 * \see odisp */
void
dispchr (ocharptr * chrc)
{
  if (*chrc != NULL)
    {
      free (*chrc);
      *chrc = NULL;
      rcount = rcount - 1;
    }
}

void
dispsfn (termptr * sfn)
{
  if (*sfn != NULL)
    {
      free (*sfn);
      *sfn = NULL;
      count = count - 1;
    }
}

void
dispprod (prodtype * prd)
{
  register int j;

  if (*prd != NULL)
    {
      pcount = pcount - 1;
      for (j = 1; j <= maxprod; j++)
	{
	  odisp (&(*prd)->prods.A[j - 1]);
	}
      free (*prd);
      *prd = NULL;
    }
}

void
fnu (fnptrs * fnx)
{
  (*fnx) = (fnptrs) malloc ((unsigned) (sizeof (*(*fnx))));
  (*fnx)->next = NULL;
}

void
fdisp (fnptrs * fnx)
{
  fnptrs temp, temp1;

  if ((*fnx) != NULL)
    {
      temp = (*fnx);
      while (temp->next != NULL)
	{
	  temp1 = temp->next;
	  free (temp);
	  temp = temp1;
	}
      (*fnx) = NULL;
    }
}

/** return length of a string0 */
int
blen (string0 * f)
{
  register int i;

  i = bcol - 1;
  do
    {
      if ((f->A[i - 1] == ' '))
	i = i - 1;
    }
  while (!((f->A[i - 1] != ' ') || (i == 1)));
  return i;
}

/** returns weight of partition */
int
wtfrm (frame * partition)
{
  register int weightx;
  register int i;

  weightx = 0;
  for (i = 1; i <= maxdim - 1; i++)
    {
      weightx = weightx + partition->A[i];
    }
  return weightx;
}

/** return comparison between partition a and b, LESS, EQUAL or GREATER 
 * stops on first 0 encountered 
 * \see qcffrmfrm */
int
cffrmfrm (frame a, frame b)
{
  register int i;

  if (a.A[1] == b.A[1])
    i = 2;
  else
    i = 1;
  while ((a.A[i] == b.A[i]) && (a.A[i] != 0) && (i < maxdim - 1))
    i = i + 1;
  return sgn (a.A[i] - b.A[i]);
}

/** return comparison between first partition of p and b, LESS, EQUAL or GREATER 
 * if p is NULL it returns LESS*/
int
cfptrfrm (termptr p, frame b)
{
  register int result;

  if (p != NULL)
    result = cffrmfrm (p->val, b);
  else
    {
      result = LESS;
      sign = LESS;
    }
  return result;
}

/** return comparison between first elements of list of partitions p and q, LESS, EQUAL or GREATER */
int
cfptrptr (termptr p, termptr q)
{
  if ((p == NULL) && (q != NULL))
    {
      sign = LESS;
      return LESS;
    }
  else if ((p != NULL) && (q == NULL))
    {
      sign = GREATER;
      return GREATER;
    }
  else
    return cffrmfrm (p->val, q->val);
}

/** return comparison between first elements ocharptr p and q, LESS, EQUAL or GREATER */
int
ccfptrptr (ocharptr p, ocharptr q)
{
  if ((p == NULL) && (q != NULL))
    {
      sign = LESS;
      return LESS;
    }
  else if ((p != NULL) && (q == NULL))
    {
      sign = GREATER;
      return GREATER;
    }
  else
    return cffrmfrm (p->val, q->val);
}

/** modify all multiplicities of a list of partitions
 * \param list list of partitions
 * \param c factor to apply
 * \param action : "-" (opposite), "*" (multiply by c), "=" (force to be c), "^" , "|" (absolute value),
 * 	"\" divide by c */
void
coeffset (termptr * list, int c, char action)
{
  termptr ptr;

  ptr = (*list);
  while (ptr != NULL)
    {
      switch ((int) (action))
	{
	case '-':		//added by FB to accelerate some computation
	  ptr->mult = -ptr->mult;
	case '*':
	  ptr->mult = ptr->mult * c;
	  break;
	case '=':
	  ptr->mult = c;
	  break;
	case '^':
	  ptr->mult = ptr->mult * minusoneto (wtfrm (&ptr->val));
	  break;
	case '|':
	  ptr->mult = abs (ptr->mult);
	  break;
	case '\\':
	  if ((ptr->mult % c == 0))
	    ptr->mult = ptr->mult / c;
	  else
	    print ("error: inappropriate coefficients\n");
	  break;
	default:
	  Caseerror (Line);
	}
      ptr = ptr->next;
    }
}

/** returns comparison of partition a and b
 * it differs from cffrmfrm : do not stop on zeros.
 * \see cffrmfrm */
int
qcffrmfrm (frame a, frame b)
{
  register int i;

  if (a.A[1] == b.A[1])
    i = 2;
  else
    i = 1;
  while ((a.A[i] == b.A[i]) && (i < maxdim))
    i = i + 1;
  return sgn (a.A[i] - b.A[i]);
}

/** returns comparison of first element of partitions' lists  first and last */
int
qtestord (frame * first, frame * last)
{
  if (first == NULL || last == NULL)	// added validity check
    aaargh ("NULL partitions list in qtestord", true);
  return qcffrmfrm ((*first), (*last));
}

/** returns comparison of first element of partitions' lists  first and last 
 * stops at first zero encountered*/
int
testord (frame * first, frame * last)
{
  return cffrmfrm ((*first), (*last));
}


/** test if first is < = or > to last, return LESS,EQUAL or GREATER  respectively*/
int
otestord (ocharptr first, ocharptr last)
{
  register int ord = EQUAL;	//modified by FB, was uninitialized

  if ((first->C6_double && last->C6_double))
    {
      ord = testord (&first->conval, &last->conval);
      if (ord == EQUAL)
	if (first->spin != last->spin)
	  if (first->spin)
	    ord = GREATER;
	  else
	    ord = LESS;
      if (ord == EQUAL)
	ord = testord (&first->val, &last->val);
    }
  else if ((first->C6_double && (!last->C6_double)))
    ord = GREATER;
  else if (((!first->C6_double) && last->C6_double))
    ord = LESS;
  if (((!first->C6_double) && (!last->C6_double)))
    {
      ord = testord (&first->val, &last->val);
      if (first->val.A[1] == last->val.A[1])
	if (first->spin != last->spin)
	  if (first->spin)
	    ord = GREATER;
	  else
	    ord = LESS;
      if ((ord == EQUAL))
	if (first->lab != last->lab)
	  if ((((first->lab == '+') || (first->lab == '-'))
	       && ((last->lab == '+') || (last->lab == '-'))))
	    {
	      if (first->lab < last->lab)
		ord = GREATER;
	      else
		ord = LESS;
	    }
	  else if (first->lab > last->lab)
	    ord = GREATER;
	  else
	    ord = LESS;
    }	 /**2/1/99**/
  if (ord == EQUAL)
    {
      if (((first->lab == '#') && (last->lab == ' ')))
	ord = GREATER;
      if (((first->lab == ' ') && (last->lab == '#')))
	ord = LESS;
    }
  return ord;
}				// otestord

/* returns the rest of the integer division of a by b*  unused
int
stmod (int a, int b)
{
  return (a - (a / b) * b);
}*/

prodtype
pcheck (prodtype pr)
{
  register prodtype sublist;
  prodtype help;
  register int j;
  bool test;

  sublist = NULL;
  while (pr != NULL)
    {
      test = false;
      for (j = 1; j <= nprod; j++)
	{
	  if (pr->prods.A[j - 1] == NULL)
	    test = true;
	}
      if (test)
	{
	  help = pr;
	  pr = pr->next;
	  dispprod (&help);
	}
      else
	{
	  help = pr->next;
	  pr->next = sublist;
	  sublist = pr;
	  pr = help;
	}
    }
  return sublist;
}

int
ptestord (prodtype first, prodtype last)
{
  register int j;
  int ord;
  j = 0;
  ord = EQUAL;
  do
    {
      j = j + 1;
      if ((currgrp.A[j - 1].name == spnc) || (currgrp.A[j - 1].name == sonc))
	ord = sotestord (first->prods.A[j - 1], last->prods.A[j - 1]);
      else			/*25/4/98 */
	ord = otestord (first->prods.A[j - 1], last->prods.A[j - 1]);
    }
  while (!((ord != EQUAL) || (j == nprod)));
  return ord;
}

void
pexpand (prodtype * sublist, prodtype * pr, int j)
{
  prodtype entry;
  ocharptr chrc, top;
  int m;
  register int k;

  top = (*pr)->prods.A[j - 1];
  while ((*pr)->prods.A[j - 1] != NULL)
    {
      if (j != nprod)
	pexpand (&(*sublist), &(*pr), j + 1);
      else
	{
	  pnu (&entry);
	  m = (*pr)->mult;
	  for (k = 1; k <= nprod; k++)
	    {
	      cnu (&chrc);
	      (*chrc) = (*(*pr)->prods.A[k - 1]);
	      m = m * chrc->mult;
	      chrc->next = NULL;
	      chrc->mult = 1;
	      entry->prods.A[k - 1] = chrc;
	      chrc = NULL;
	    }
	  entry->mult = m;
	  entry->next = (*sublist);
	  (*sublist) = entry;
	}
      (*pr)->prods.A[j - 1] = (*pr)->prods.A[j - 1]->next;
    }
  (*pr)->prods.A[j - 1] = top;
}

prodtype
prodexpand (prodtype pr)
{
  prodtype sublist;

  sublist = NULL;
  while (pr != NULL)
    {
      pexpand (&sublist, &pr, 1);
      pr = pr->next;
    }

  return sublist;
}

termptr
swt (termptr list)
{
  termptr newlist, temp = NULL;
  int w;
  register int i;

  newlist = NULL;
  i = 1;
  w = 0;
  while (list != NULL)
    {
      if (newlist == NULL)
	{
	  snu (&newlist);
	  temp = newlist;
	}
      else
	{
	  snu (&temp->next);
	  temp = temp->next;
	}
      (*temp) = (*list);
      while (temp->val.A[i] != 0)
	{
	  w = w + temp->val.A[i];
	  i = i + 1;
	}
      for (i = 1; i <= maxdim; i++)
	{
	  temp->val.A[i] = 0;
	}
      temp->val.A[1] = w;
      w = 0;
      i = 1;
      list = list->next;
    }
  sort (&newlist, false);
  return newlist;
}

ocharptr
rwt (ocharptr list)
{
  ocharptr newlist, temp = NULL;	//modified by FB, was uninitialized
  int w;
  register int i;

  newlist = NULL;
  i = 1;
  w = 0;
  while (list != NULL)
    {
      if (newlist == NULL)
	{
	  cnu (&newlist);
	  temp = newlist;
	}
      else
	{
	  cnu (&temp->next);
	  temp = temp->next;
	}
      (*temp) = (*list);
      while (temp->val.A[i] != 0)
	{
	  w = w + temp->val.A[i];
	  i = i + 1;
	}
      for (i = 1; i <= maxdim; i++)
	{
	  temp->val.A[i] = 0;
	}
      temp->val.A[1] = w;
      w = 0;
      i = 1;
      list = list->next;
    }
  osort (&newlist, false);
  return newlist;
}

//Compute the conjugate of a partition f
void
conjgte (frame * f)
{
  int partc, edge = 0, saveheight;
  register int i;
  frame rho;

  partc = 1;
  rho = nolls;
  saveheight = f->A[1];
  if (saveheight < maxdim)
    while (f->A[partc] != 0)
      {
	edge = f->A[partc];
	do
	  {
	    partc = partc + 1;
	  }
	while (!(f->A[partc] != edge));
	for (i = f->A[partc] + 1; i <= edge; i++)
	  {
	    rho.A[i] = partc - 1;
	  }
      }
  else
    error (BAD_CONJUGATE, 0);
  rho.length = saveheight;
  (*f) = rho;
}

termptr
lconjgte (termptr x)
{
  termptr conjx;

  conjx = x;
  while (x != NULL)
    {
      conjgte (&x->val);
      x = x->next;
    }
  sort (&conjx, false);
  return conjx;
}

bool
sconjgte (frame * x)
{
  int lx, ly;
  register int i;
  bool test;
  frame y;

  y = nolls;
  lx = len (x);
  y = *x;
  conjgte (&y);
  ly = len (&y);
  if (lx != ly)
    test = false;
  else
    {
      test = true;
      for (i = 1; i <= lx; i++)
	{
	  if (x->A[i] != y.A[i])
	    test = false;
	}
    }
  return test;
}

void
add (termptr * top, termptr * extra)
{
  termptr last;

  if ((*extra) != NULL)
    {
      last = (*extra);
      while (last->next != NULL)
	last = last->next;
      last->next = (*top);
      (*top) = (*extra);
      (*extra) = NULL;
    }
}

void
oadd (ocharptr * top, ocharptr * extra)
{
  ocharptr last;

  if (*extra != NULL)
    {
      last = *extra;
      while (last->next != NULL)
	last = last->next;
      last->next = *top;
      *top = *extra;
      *extra = NULL;
    }
}

termptr
ladd (termptr a, termptr b)
{
  termptr head, lastptr, point;

  snu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (a != NULL)
    {
      snu (&point);
      *point = *a;
      a = a->next;
      lastptr->next = point;
      lastptr = point;
    }
  while (b != NULL)
    {
      snu (&point);
      *point = *b;
      b = b->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  sort (&head, true);
  dispsfn (&lastptr);
  return head;
}

termptr
sfncopy (termptr sfn)
{
  termptr head, lastptr, point;

  snu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (sfn != NULL)
    {
      snu (&point);
      *point = *sfn;
      sfn = sfn->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  sort (&head, true);
  dispsfn (&lastptr);
  return head;
}

termptr
absval (termptr sf)
{
  termptr t;

  t = sfncopy (sf);
  if (t != NULL)		// sfncopy can return an empty list=NULL
    coeffset (&t, 0, '|');
  return t;
}

ocharptr
chrccopy (ocharptr chrc)
{
  ocharptr head, lastptr, point;

  cnu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (chrc != NULL)
    {
      cnu (&point);
      *point = *chrc;
      chrc = chrc->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  if (currgrp.A[1 - 1].name == spnc)
    sposort (&head, true);
  else
    osort (&head, true);

  /*osort(&head, true); */
  dispchr (&lastptr);
  return head;
}

prodtype
prodcopy (prodtype prd)
{
  prodtype head, lastptr, point;
  register int j;
  pnu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (prd != NULL)
    {
      pnu (&point);
      *point = *prd;
      for (j = 0; j < maxprod; j++)
	{
	  point->prods.A[j] = gchrccopy (prd->prods.A[j]);
	}
      prd = prd->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  schur_psort (&head, true);
  dispprod (&lastptr);
  return head;
}

void
mergemin (frame one, frame two, frame * three)
{
  register int i;
  for (i = 0; i <= maxdim; i++)
      three->A[i] = 0;
  i = 1;
  while ((one.A[i] != 0) || (two.A[i] != 0))
    {
      three->A[i] = MIN (one.A[i], two.A[i]);
      i = i + 1;
    }
}

void
mergemax (frame one, frame two, frame * three)
{
  register int i;
  for (i = 0; i <= maxdim; i++)
      three->A[i] = 0;
  i = 1;
  while ((one.A[i] != 0) || (two.A[i] != 0))
    {
      three->A[i] = MAX (one.A[i], two.A[i]);
      i = i + 1;
    }
}

termptr
sameweight (int wtfrmx, int step)
{
  int left, edge1, edge2, sum, prev;
  register int i;
  frame rho;
  termptr result, ptr1;

  edge2 = 1;
  sum = 0;
  rho = nolls;
  rho.A[1] = wtfrmx;
  snu (&result);
  result->val = rho;
  result->mult = 1;
  ptr1 = result;
  while (rho.A[1] > step)
    {
      edge1 = edge2;
      edge2 = -1;
      left = sum + rho.A[edge1];
      prev = rho.A[edge1] - step;
      do
	{
	  prev = MIN (left, prev);
	  if ((edge2 < 0) && (prev <= step))
	    {
	      edge2 = edge1 - 1;
	      sum = left;
	    }
	  rho.A[edge1] = prev;
	  left = left - prev;
	  edge1 = edge1 + 1;
	}
      while (!(prev == 0));
      snu (&ptr1->next);
      ptr1 = ptr1->next;
      ptr1->val = rho;
      ptr1->mult = 1;
      for (i = edge1 - 1; i <= maxdim; i++)
	{
	  ptr1->val.A[i] = 0;
	}
    }
  dispsfn (&ptr1->next);
  ptr1->next = NULL;
  return result;
}

void
limit (termptr * list, int n)
{
  termptr newlist, temp;
  int i;
  bool test;

  newlist = NULL;
  while (*list != NULL)
    {
      register termptr W32 = *list;

      i = 1;
      test = true;
      if ((n < 0))
	do
	  {
	    if ((((n == -1) && (bool) ((W32->val.A[i]) & 1))
		 || ((n == -2) && (!(bool) ((W32->val.A[i]) & 1)))))
	      test = false;
	    else
	      i = i + 1;
	  }
	while (!((W32->val.A[i] == 0) || (i == maxdim) || (test == false)));
      else if ((n > 0))
	do
	  {
	    if ((W32->val.A[i] == n))
	      test = false;
	    else
	      i = i + 1;
	  }
	while (!((W32->val.A[i] == 0) || (i == maxdim) || (test == false)));
      temp = W32->next;
      if ((test == false))
	dispsfn (&(*list));
      else
	{
	  W32->next = newlist;
	  newlist = *list;
	}
      *list = temp;
    }
  sort (&newlist, true);
  *list = newlist;
}

termptr
skcompat (frame x)
{
  return allskews (x, false);
}

termptr
lskcmpt (termptr x)
{
  register termptr temp;
  frame big;

  big = nolls;
  temp = NULL;
  if (x != NULL)
    {
      big = x->val;
      x = x->next;
      while (x != NULL)
	{
	  mergemax (big, x->val, &big);
	  x = x->next;
	}
      temp = skcompat (big);
    }
  return temp;
}

termptr
eqwt (int n)
{
  return sameweight (n, 1);
}

termptr
sfnmult (int factor, termptr sfn)
{
  termptr t;

  if (factor == 0)
    return NULL;
  else
    {
      t = sfncopy (sfn);
      if (t != NULL)		// sfncopy can return an empty list=NULL
	coeffset (&t, factor, '*');
    }
  return t;
}

ocharptr
chrcadd (ocharptr a, ocharptr b)
{
  ocharptr head, lastptr, point;

  cnu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (a != NULL)
    {
      cnu (&point);
      (*point) = (*a);
      a = a->next;
      lastptr->next = point;
      lastptr = point;
    }
  while (b != NULL)
    {
      cnu (&point);
      (*point) = (*b);
      b = b->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == sonc))
    sposort (&head, true);
  else
    osort (&head, true);
  /*osort(&head, true); */
  dispchr (&lastptr);
  return head;
}

ocharptr
chrcmult (int factor, ocharptr chrc)
{
  register ocharptr point;

  if (factor != 0)
    {
      if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == sonc))
	chrc = gchrccopy (chrc);
      else
	chrc = chrccopy (chrc);
      point = chrc;
      while (point != NULL)
	{
	  point->mult *= factor;
	  point = point->next;
	}
      return chrc;
    }
  else
    return NULL;
}

prodtype
prodadd (prodtype a, prodtype b)
{
  prodtype head, lastptr, point;
  register int i;

  pnu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (a != NULL)
    {
      pnu (&point);
      point->mult = a->mult;
      for (i = 1; i <= nprod; i++)
	{
	  point->prods.A[i - 1] = gchrccopy (a->prods.A[i - 1]);
	}
      a = a->next;
      lastptr->next = point;
      lastptr = point;
    }
  while (b != NULL)
    {
      pnu (&point);
      point->mult = b->mult;
      for (i = 1; i <= nprod; i++)
	{
	  point->prods.A[i - 1] = chrccopy (b->prods.A[i - 1]);
	}
      b = b->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  schur_psort (&head, true);
  dispprod (&lastptr);
  return head;
}

prodtype
prodmult (int factor, prodtype prd)
{
  register prodtype point;

  if (factor != 0)
    {
      prd = prodcopy (prd);
      point = prd;
      while (point != NULL)
	{
	  point->mult *= factor;
	  point = point->next;
	}
      return prd;
    }
  else
    return NULL;
}

ocharptr
sfntochrc (termptr sfn, bool spinor, char tag)
{
  ocharptr lastchrc, headchrc, listchrc;

  cnu (&lastchrc);
  headchrc = lastchrc;
  headchrc->next = NULL;
  while (sfn != NULL)
    {
      cnu (&listchrc);
      {
	listchrc->mult = sfn->mult;
	listchrc->val = sfn->val;
	if (sfn->slab == '#')
	  listchrc->lab = '#';
	else
	  listchrc->lab = tag;
	listchrc->spin = spinor;
	listchrc->next = NULL;
	listchrc->C6_double = false;
	lastchrc->next = listchrc;
	lastchrc = listchrc;
	sfn = sfn->next;
      }
    }
  lastchrc = headchrc;
  headchrc = headchrc->next;
  dispchr (&lastchrc);
  return headchrc;
}

void
sstkhandle (termptr sfn)
{
  if (srjctindex < srjctlimit)
    srjctindex = srjctindex + 1;
  else
    {
      warn ("stack overloaded", cont);
      inform (":non fatal;", cr);
    }
  sreject.A[srjctindex - 1] = sfn;
}

void
stackhandle (ocharptr chrc, bool stack)
{
  if (stack)
    {
      if (rjctindex < rjctlimit)
	rjctindex++;
      else
	{
	  warn ("stack overloaded", cont);
	  inform (":non fatal;", cr);
	}
      reject.A[rjctindex - 1] = chrc;
    }
  else
    {
      if (rjctindex2 < rjctlimit)
	rjctindex++;
      else
	{
	  warn ("stack overloaded", cont);
	  inform (":non fatal;", cr);
	}
      reject2.A[rjctindex - 1] = chrc;
    }
}

void
pstackhandle (prodtype prd)
{
  if (prjctindex < prjctlimit)
    prjctindex = prjctindex + 1;
  else
    {
      warn ("stack overloaded", cont);
      inform (":non fatal;", cr);
    }
  preject.A[prjctindex - 1] = prd;
}

int
qqlen (frame f)
{
  register int i;

  i = maxdim;
  do
    {
      if (f.A[i] == 0)
	i = i - 1;
    }
  while (!((f.A[i] != 0) || (i == 1)));
  return i;
}


int
sotestord (ocharptr first, ocharptr last)
{
  int ord;

  {
    if ((first->val.A[maxdim] > last->val.A[maxdim]))
      ord = GREATER;
    else if ((first->val.A[maxdim] < last->val.A[maxdim]))
      ord = LESS;
    else
      ord = EQUAL;
  }
  if (ord == EQUAL)
    {
      if (first->spin != last->spin)
	if (first->spin)
	  ord = GREATER;
	else
	  ord = LESS;
      if (first->spin == last->spin)
	ord = EQUAL;
    }
  if (ord == EQUAL)
    ord = testord (&first->val, &last->val);	 /**2/1/99**/
  if (ord == EQUAL)
    {
      if (((first->lab == '#') && (last->lab == ' ')))
	ord = GREATER;
      if (((first->lab == ' ') && (last->lab == '#')))
	ord = LESS;
    }
  return ord;
}

void
sposort (ocharptr * head, bool reverse)
{
  ocharptr lastptr, sortedlist, test, temp, tail;
  bool inserted;
  int ord;

  lastptr = NULL;
  progress ();
  if (reverse)
    {
      while (*head != NULL)
	{
	  tail = (*head)->next;
	  (*head)->next = lastptr;
	  lastptr = *head;
	  *head = tail;
	}
      *head = lastptr;
    }
  cnu (&tail);
  tail->next = NULL;
  cnu (&lastptr);
  lastptr->next = tail;
  sortedlist = lastptr;
  tail = lastptr;
  while (*head != NULL)
    {
      test = (*head);
      (*head) = (*head)->next;
      lastptr = sortedlist;
      temp = sortedlist->next;
      inserted = false;
      do
	{
	  if (temp->next != NULL)
	    {
	      ord = sotestord (test, temp);
	      switch ((int) (ord))
		{
		case -1:
		  lastptr = temp;
		  temp = temp->next;
		  break;
		case 0:
		  temp->mult = temp->mult + test->mult;
		  dispchr (&test);
		  inserted = true;
		  break;
		case 1:
		  lastptr->next = test;
		  test->next = temp;
		  inserted = true;
		  break;
		default:
		  Caseerror (Line);
		}
	    }
	  else
	    {
	      lastptr->next = test;
	      test->next = temp;
	      inserted = true;
	      tail = test;
	    }
	}
      while (!(inserted));
    }
  dispchr (&tail->next);
  tail->next = NULL;
  (*head) = sortedlist;
  tail = (*head)->next;
  while (tail != NULL)
    {
      if (tail->mult == 0)
	{
	  (*head)->next = tail->next;
	  dispchr (&tail);
	  tail = (*head)->next;
	}
      else
	{
	  (*head) = tail;
	  tail = tail->next;
	}
    }
  (*head) = sortedlist->next;
  dispchr (&sortedlist);
}


ocharptr
gchrccopy (ocharptr chrc)
{
  ocharptr head, lastptr, point;

  cnu (&lastptr);
  head = lastptr;
  head->next = NULL;
  while (chrc != NULL)
    {
      cnu (&point);
      (*point) = (*chrc);
      chrc = chrc->next;
      lastptr->next = point;
      lastptr = point;
    }
  lastptr = head;
  head = head->next;
  if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == sonc))
    sposort (&head, true);
  else
    osort (&head, true);

  dispchr (&lastptr);
  return head;
}

void
coeffsetch (ocharptr * list, int c, char action)
{
  ocharptr ptr;

  ptr = (*list);
  while (ptr != NULL)
    {
      switch ((int) (action))
	{
	case '*':
	  ptr->mult = ptr->mult * c;
	  break;
	case '=':
	  ptr->mult = c;
	  break;
	case '^':
	  ptr->mult = ptr->mult * minusoneto (wtfrm (&ptr->val));
	  break;
	case '|':
	  ptr->mult = abs (ptr->mult);
	  break;
	case '\\':
	  if ((ptr->mult % c == 0))
	    ptr->mult = ptr->mult / c;
	  else
	    print ("error: inappropriate coefficients\n");
	  break;
	default:
	  Caseerror (Line);
	}
      ptr = ptr->next;
    }
}

ocharptr
absvalch (ocharptr sf)
{
  ocharptr t;

  t = chrccopy (sf);
  coeffsetch (&t, 0, '|');
  return t;
}

bool
conjtest (frame * x)
{
  int lx, ly;
  register int i;
  bool test;
  frame y;

  lx = len (x);
  y = (*x);
  conjgte (&y);
  ly = len (&y);
  test = true;
  if (lx > ly)
    test = false;
  else if (lx == ly)
    for (i = 1; i <= lx; i++)
      {
	if (x->A[i] < y.A[i])
	  test = false;
      }
  return test;
}


termptr
conjadd (termptr list)
{
  termptr newlist, temp, tlist;
  bool test;

  newlist = NULL;
  tlist = list;
  while (list != NULL)
    {

      snu (&temp);
      temp->mult = list->mult;
      temp->val = list->val;
      test = conjtest (&temp->val);
      if (test == false)
	conjgte (&temp->val);
      add (&newlist, &temp);
      dispsfn (&temp);
      list = list->next;
    }
  ldisp (&tlist);
  sort (&newlist, true);
  return newlist;
}
