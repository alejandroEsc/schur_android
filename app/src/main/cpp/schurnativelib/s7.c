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

# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "gr.h"
# include "s2.h"
# include "s5.h"
# include "bignums.h"
# include "g.h"
# include "s7.h"


termptr
rsameweight (frame rho)
{
  register termptr R169;
  int left, edge1, edge2, sum, prev, k;
  register int i;
  termptr result, ptr1;
  frame temp;
  temp = rho;
  sum = 0;
  k = len (&rho);
  edge2 = 0;
  for (i = 1; i <= k; i++)
    {
      if (rho.A[i] > 1)
	edge2 = edge2 + 1;
    }
  for (i = 1; i <= k; i++)
    {
      if (rho.A[i] == 1)
	sum = sum + 1;
    }
  snu (&result);
  result->val = temp;
  result->mult = 1;
  ptr1 = result;
  while (temp.A[1] > 1)
    {
      edge1 = edge2;
      edge2 = -1;
      left = sum + temp.A[edge1];
      prev = temp.A[edge1] - 1;

      do
	{
	  prev = MIN (left, prev);
	  if ((edge2 < 0) && (prev <= 1))
	    {
	      edge2 = edge1 - 1;
	      sum = left;
	    }
	  temp.A[edge1] = prev;
	  left = left - prev;
	  edge1 = edge1 + 1;
	}
      while (!(prev == 0));
      snu (&ptr1->next);
      ptr1 = ptr1->next;
      ptr1->val = temp;
      ptr1->mult = 1;
      for (i = edge1 - 1; i <= maxdim; i++)
	{
	  ptr1->val.A[i] = 0;
	}
    }
  dispsfn (&ptr1->next);
  ptr1->next = NULL;
  R169 = result;
  return R169;
}

termptr
lsfntomono (termptr list)
{
  register termptr R170;
  termptr newlist, tlist, temp;

  /*int       m; */
  bframe bigg;
  newlist = NULL;
  while (list != NULL)
    {

      /*m = list->mult; */
      tlist = rsameweight (list->val);
      while (tlist != NULL)

	{
	  register termptr W7 = &(*tlist);
	  snu (&temp);
	  kostka (list, tlist, &bigg);
	  temp->mult = frbig (bigg);
	  temp->val = W7->val;
	  add (&newlist, &temp);
	  dispsfn (&temp);
	  tlist = W7->next;
	} sort (&newlist, true);
      list = list->next;
    } R170 = newlist;
  return R170;
}

void
msort (arrayptr * l, int n)
{
  int q, i;
  register int j;
  bool done;
  i = 1;

  do
    {
      done = true;
      for (j = n; j >= i + 1; j--)
	{
	  if ((*l)->A[j] > (*l)->A[j - 1])
	    {
	      q = (*l)->A[j];
	      (*l)->A[j] = (*l)->A[j - 1];
	      (*l)->A[j - 1] = q;
	      done = false;
	    }
	}
      i = i + 1;
    }
  while (!(done || (i == n)));
}

void
ifzero (arrayptr * l, int n, int *rs)
{
  int i;
  (*rs) = 0;
  i = 1;
  while (((*rs) == 0) && (i < n))
    {
      if ((*l)->A[i] == (*l)->A[i + 1] - 1)
	(*rs) = 1;
      i = i + 1;
    }
}
void
cutzero (arrayptr * l, int n, int *m)
{
  int i;
  bool done;
  i = n;
  done = true;
  while (done && (i > 0))
    if ((*l)->A[i] == 0)
      i = i - 1;

    else
      done = false;
  (*m) = i;
  if ((*m) > 0)
    {
      i = 1;
      done = true;
      while (done)
	if ((*l)->A[i] == 0)
	  i = i + 1;

	else
	  done = false;
      if (i - 1 >= (*l)->A[i])
	(*m) = 0;
    }
}
void
xstndise (arrayptr * l, int n, int *m, int *znak)
{
  arrayptr list, temp;
  int q, t, s, rs;
  register int i, j;
  bool done;
  list = (frame *) malloc ((unsigned) (sizeof (*list)));
  temp = (frame *) malloc ((unsigned) (sizeof (*temp)));
  (*znak) = 1;
  cutzero (&(*l), n, &(*m));
  if ((*m) > 0)
    {
      ifzero (&(*l), (*m), &rs);
      if (rs != 1)
	{
	  j = 0;
	  while ((rs != 1) && ((*m) - j > 1))
	    {
	      j = j + 1;
	      t = (*m) - j + 1;
	      for (i = 1; i <= t; i++)
		{
		  q = j + i - 1;
		  list->A[i] = (*l)->A[q] - i + 1;
		  temp->A[i] = list->A[i];
		}
	      msort (&temp, t);
	      s = 1;

	      do
		{
		  done = false;
		  if (temp->A[1] == list->A[s])
		    done = true;
		  s = s + 1;
		}
	      while (!(done && (s <= t + 1)));
	      s = s - 1;
	      if (s > 1)
		for (i = 1; i <= s - 1; i++)
		  {
		    (*znak) = (*znak) * (-1);
		  }
	      for (i = j + s - 1; i >= j + 1; i--)
		{
		  (*l)->A[i] = (*l)->A[i - 1] + 1;
		}
	      (*l)->A[j] = list->A[s];
	      s = j;
	      if (t < (*m))
		{
		  t = t + 1;
		  s = j - 1;
		}
	      for (i = 1; i <= t; i++)
		{
		  q = s + i - 1;
		  list->A[i] = (*l)->A[q];
		}
	      ifzero (&list, t, &rs);
	      if (rs == 1)
		(*m) = 0;
	    }
	}
      else
	(*m) = 0;
    }
  free (list);
  free (temp);
}

void
smallcs (arrayptr * l, arrayptr * s, int n, int losp, termptr * wtp,
	 termptr * shoot, int multy)
{
  int znak, m;
  register int q;
  if (n < 0)
    print ( "Gordan not defined on %10d elements.\n");

  else if (n == 2)
    {
      if ((*l)->A[1] != (*l)->A[2])
	{
	  q = (*l)->A[2];
	  (*l)->A[2] = (*l)->A[1];
	  (*l)->A[1] = q;
	  for (q = 1; q <= losp; q++)
	    {
	      (*l)->A[q] = (*l)->A[q] + (*s)->A[q];
	    }
	  xstndise (&(*l), n, &m, &znak);
	  if (m > 0)
	    {
	      snu (&(*wtp));
	      (*wtp)->val = (*(*l));
	      (*wtp)->mult = znak * multy;
	      (*wtp)->next = (*shoot);
	      (*shoot) = (*wtp);
	    }
	}
    }
}
void
permute (arrayptr * l, arrayptr s, int n, int losp, bool * stop,
	 termptr * wtp, termptr * shoot, int multy)
{
  arrayptr r;
  int q, c, k, t, znak;
  register int j;
  bool done;
  r = (frame *) malloc ((unsigned) (sizeof (*r)));
  k = n;
  q = k - 1;
  done = true;
  while (((*l)->A[k] >= (*l)->A[q]) && done)
    {
      k = k - 1;
      q = q - 1;
      if (k == 1)
	{
	  q = 1;
	  done = false;
	}
    }
  (*stop) = false;
  if (done)
    {
      c = 1;
      t = n - k + 1;
      for (j = k; j <= n; j++)
	{
	  r->A[c] = (*l)->A[j];
	  c = c + 1;
	}
      msort (&r, t);
      c = 1;
      while ((r->A[c] >= (*l)->A[k - 1]) && (c <= t))
	{
	  c = c + 1;
	}
      q = (*l)->A[k - 1];
      (*l)->A[k - 1] = r->A[c];
      r->A[c] = q;
      c = 1;
      for (j = k; j <= n; j++)
	{
	  (*l)->A[j] = r->A[c];
	  c = c + 1;
	}
      for (j = 0; j <= maxdim; j++)
	{
	  r->A[j] = 0;
	}
      for (j = 1; j <= n; j++)
	{
	  r->A[j] = (*l)->A[j];
	  if (j <= losp)
	    r->A[j] = r->A[j] + s->A[j];
	}
      xstndise (&r, n, &q, &znak);
      if (q > 0)
	{
	  snu (&(*wtp));
	  (*wtp)->val = (*r);
	  (*wtp)->mult = znak * multy;
	  (*wtp)->next = (*shoot);
	  (*shoot) = (*wtp);
	}
      free (r);
      (*stop) = true;
    }
}
void
gordan3 (arrayptr * l, arrayptr * s, int n, int q, termptr * wtp, int multy)
{
  arrayptr r;
  termptr shoot;
  int c, i;
  bool stop;
  r = (frame *) malloc ((unsigned) (sizeof (*r)));
  for (i = 1; i <= maxdim; i++)
    {
      r->A[i] = 0;
    }
  c = 1;
  r->A[c] = (*l)->A[c];
  if (c <= q)
    r->A[c] = r->A[c] + (*s)->A[c];
  while ((c <= n) && (r->A[c] > 0))
    {
      c = c + 1;
      r->A[c] = (*l)->A[c];
      if (c <= q)
	r->A[c] = r->A[c] + (*s)->A[c];
    }
  c = c - 1;
  shoot = NULL;
  (*wtp)->val = (*r);
  (*wtp)->mult = multy;
  (*wtp)->next = shoot;
  shoot = (*wtp);
  free (r);
  stop = true;
  if (n < 3)
    smallcs (&(*l), &(*s), n, q, &(*wtp), &shoot, multy);

  else
    while (stop)
      permute (&(*l), (*s), n, q, &stop, &(*wtp), &shoot, multy);
}

termptr
gordan2 (termptr sfn1, termptr sfn2)
{
  arrayptr l, s;
  termptr wtp;
  int l2, n, multy;
  register int i;
  //l1 = len (&sfn1->val);
  l2 = len (&sfn2->val);
  multy = sfn1->mult * sfn2->mult;
  n = l2 + wtfrm (&sfn1->val);
  l = (frame *) malloc ((unsigned) (sizeof (*l)));
  s = (frame *) malloc ((unsigned) (sizeof (*s)));
  for (i = 1; i <= maxdim; i++)
    {
      l->A[i] = 0;
      s->A[i] = 0;
    }
  (*l) = sfn1->val;
  (*s) = sfn2->val;
  msort (&l, n);
  msort (&s, l2);
  snu (&wtp);
  gordan3 (&l, &s, n, l2, &wtp, multy);
  free (l);
  free (s);
  sort (&wtp, true);
  return wtp;
}

termptr
gordan (termptr list1, termptr list2)
{
  register termptr R172;
  termptr templist, sublist, prodlist;
  prodlist = NULL;
  while (list1 != NULL)
    {
      templist = list2;
      while (templist != NULL)
	{
	  sublist = gordan2 (list1, templist);
	  add (&prodlist, &sublist);
	  templist = templist->next;
	}
      sort (&prodlist, true);
      list1 = list1->next;
    }
  R172 = prodlist;
  return R172;
}

void
multpart (termptr list, int pn)
{
  int i, p, q, l;
  p = 1;
  while (list != NULL)

    {
      register termptr tempo = list;
      i = 1;
      l = len (&list->val);
      while (i <= l)
	{
	  q = tempo->val.A[i] * pn;
	  if (q >=maxdim)
	    error (PART_TOO_BIG, p);
	  tempo->val.A[i] = q;
	  i = i + 1;
	}
      list = list->next;
    }
}

termptr
monotosfn (termptr mono)
{
  register termptr R173;
  termptr temp;
  register int i;
  snu (&temp);
  temp->mult = 1;
  for (i = 1; i <= maxdim; i++)
    {
      temp->val.A[i] = 0;
    }
  R173 = gordan (mono, temp);
  dispsfn (&temp);
  return R173;
}

termptr
elemtomono (termptr emlist)
{
  register termptr R174;
  termptr temp;
  temp = elemtosfn (emlist);
  R174 = lsfntomono (temp);
  return R174;
}

termptr
monotoelem (termptr melist)
{
  register termptr R175;
  termptr temp;
  temp = monotosfn (melist);
  R175 = sfntoelem (temp);
  return R175;
}

termptr
homotomono (termptr hmlist)
{
  register termptr R176;
  termptr temp;
  temp = homotosfn (hmlist);
  R176 = lsfntomono (temp);
  return R176;
}

termptr
monotohomo (termptr mhlist)
{
  register termptr R177;
  termptr temp;
  temp = monotosfn (mhlist);
  R177 = sfntohomo (temp);
  return R177;
}

termptr
pntosfn (int n)
{
  register termptr R178;
  termptr temp, result;
  register int k;
  register int i;
  result = NULL;
  for (i = 1; i <= n; i++)
    {
      snu (&temp);
      for (k = 1; k <= maxdim; k++)
	{
	  temp->val.A[k] = 0;
	}
      if ((bool) ((i + 1) & 1))
	temp->mult = -1;

      else
	temp->mult = 1;
      temp->val.A[1] = n - i + 1;
      if (i > 1)
	for (k = 2; k <= i; k++)
	  {
	    temp->val.A[k] = 1;
	  }
      merge (&result, &temp, true, true);
    }
  R178 = result;
  return R178;
}

termptr
powersumtosfn (termptr sfn)
{
  register termptr R179;
  termptr result, onesfn, temp, list;
  int k, m;
  register int i;
  result = NULL;
  while (sfn != NULL)
    {
      /*register termptr W41 = &(*sfn); *//*12/12/95 */
      k = len (&sfn->val);
      snu (&onesfn);
      onesfn->mult = 1;
      m = sfn->mult;
      for (i = 1; i <= maxdim; i++)
	{
	  onesfn->val.A[i] = 0;
	}
      for (i = 1; i <= k; i++)
	{
	  list = pntosfn (sfn->val.A[i]);
	  temp = louter (onesfn, list);
	  ldisp (&onesfn);
	  onesfn = temp;
	  ldisp (&list);
	}
      temp = sfnmult (m, onesfn);
      merge (&result, &temp, true, true);
      sort (&result, true);
      ldisp (&onesfn);
      sfn = sfn->next;
    }
  R179 = result;
  return R179;
}
