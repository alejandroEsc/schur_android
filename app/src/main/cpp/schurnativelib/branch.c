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

/** \file branch.c
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
# include "s1.h"
# include "s2.h"
# include "s4.h"
# include "s6.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "q1.h"
# include "q2.h"
# include "g.h"
# include "skew.h"
# include "label.h"
# include "branch.h"

byte *G275_n;
byte *G273_l;
ocharptr *G271_lambda1;
rep_ptr *G269_temp_ptr;
rep_ptr *G267_g2_head;
rep_ptr *G265_g2_currentb;
rep_ptr *G263_so7_head;
rep_ptr *G261_so7_currentb;
ocharptr *G259_lambda;


ocharptr
uadd (ocharptr list, int extra, int index)
{
  register ocharptr R241;
  ocharptr newlist, temp;

  newlist = NULL;
  temp = NULL;
  while (list != NULL)
    {
      /*register ocharptr W2 = &(*list); *//*12/12/95 */

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
      if (index == 1)
	temp->val.A[1] = temp->val.A[1] + extra;
      else if (index == 2)
	{
	  if (extra == 0)
	    (void) fprintf (output.fp, "division by zero not allowed\n"),
	      Putl (output, 1);
	  else
	    {
	      if ((temp->val.A[1] % extra == 0))
		temp->val.A[1] = temp->val.A[1] / extra;
	      else
		(void) fprintf (output.fp, "error: inappropriate divisor\n"),
		  Putl (output, 1);
	    }
	}
      list = list->next;
    }
  R241 = newlist;
  return R241;
}

ocharptr
snanbrnch (ocharptr snrep, int n)
{
  register ocharptr R242;
  ocharptr temp, newlist;
  bool sj;
  frame dummy;
  int pts;

  dummy = nolls;
  newlist = NULL;
  while (snrep != NULL)
    {
      register ocharptr W3 = &(*snrep);

      cnu (&temp);
      temp->mult = W3->mult;
      temp->val = W3->val;
      temp->spin = W3->spin;
      temp->lab = W3->lab;
      temp->C6_double = false;
      if ((!W3->spin))
	{
	  sj = sconjgte (&W3->val);
	  if (sj)
	    anexp (&temp);
	  else
	    {
	      dummy = W3->val;
	      conjgte (&dummy);
	      if (dummy.A[1] > W3->val.A[1])
		temp->val = dummy;
	    }
	}
      else
	{
	  pts = n - len (&W3->val);
	  if ((!(bool) ((pts) & 1)))
	    dlabel (temp);
	  if ((((bool) ((pts) & 1)) && (W3->lab == ' ')))
	    anexp (&temp);
	}
      oadd (&newlist, &temp);
      snrep = W3->next;
    }
  osort (&newlist, true);
  R242 = newlist;
  return R242;
}

ocharptr
onsnbr1 (ocharptr onrep, int n, bool t)
{
  register ocharptr R243;
  termptr dummy, ss, sub1, temp, zero, undum;
  ocharptr brnch;
  bool test;
  char ch;
  register int ib;

  progress ();
  test = false;
  if (t)
    ch = 'g';
  else
    ch = 'c';
  snu (&zero);
  {
    register termptr W4 = &(*zero);
    for (ib = 0; ib <= maxdim; ib++)
      {
	W4->val.A[ib] = 0;
      }
    W4->next = NULL;
  }
  snu (&undum);
  undum->mult = 1;
  undum->next = NULL;
  temp = NULL;
  zero->mult = onrep->mult;
  undum->val = onrep->val;
  if (onrep->spin)
    {
      if (t)
	ss = seriesx ('a', -1, onrep->val);
      else
	ss = seriesx ('e', -1, onrep->val);
    }
  else
    ss = seriesx (ch, -1, onrep->val);
  sub1 = lskew (undum, ss);
  ldisp (&ss);
  dispsfn (&undum);
  ss = plethonerinner (sub1);
  ldisp (&sub1);
  dummy = ladd (temp, ss);
  temp = dummy;
  ldisp (&ss);
  ss = makeweight (n, temp);
  ldisp (&temp);
  if (onrep->lab == '#')
    ss = lconjgte (ss);
  if (onrep->spin)
    {
      if ((bool) ((n) & 1))
	qsn = 1;
      else
	qsn = 0;
      if (((!(bool) ((n) & 1)) && (t == false)))
	qsn = 2;
      test = true;
      temp = craise (ss, 1);
      ldisp (&ss);
      ss = reducedq (temp);
      ldisp (&temp);
      redu = true;
      temp = qmult (ss);
      redu = false;
      ldisp (&ss);
      ss = temp;
    }
  brnch = formbb (zero, ss, true, test, ' ');
  ldisp (&ss);
  snselect (&brnch, n, qspecial);
  osort (&brnch, true);
  dispsfn (&zero);
  R243 = brnch;
  return R243;
}

ocharptr
onsnbr (ocharptr lambda, int n, bool t)
{
  register ocharptr R244;
  ocharptr brnch, subbrnch;
  int w;

  brnch = NULL;
  if (t)
    w = n;
  else
    w = n + 1;
  while (lambda != NULL)
    {
      register ocharptr W7 = &(*lambda);

      subbrnch = onsnbr1 (lambda, w, t);
      oadd (&brnch, &subbrnch);
      lambda = W7->next;
    }
  osort (&brnch, true);
  R244 = brnch;
  return R244;
}



ocharptr
so2k1so3brnch (char br, ocharptr unirr, int n, bool test)
{
  register ocharptr R245;
  termptr ss, sub1, sub, zero, undum;
  ocharptr brnch, osub;
  register int ib;

  snu (&zero);
  {
    register termptr W8 = &(*zero);
    for (ib = 0; ib <= maxdim; ib++)
      {
	W8->val.A[ib] = 0;
      }
    W8->next = NULL;
  }
  snu (&undum);
  undum->mult = 1;
  undum->next = NULL;
  brnch = NULL;
  while (unirr != NULL)
    {
      register ocharptr W11 = &(*unirr);

      if (((br == 'c') && W11->spin))
	(void) fprintf (output.fp, "not implemented for spin irreps\n"),
	  Putl (output, 1);
      else
	{
	  zero->mult = W11->mult;
	  undum->val = W11->val;
	  ss = seriesx (br, -1, W11->val);
	  sub1 = lskew (undum, ss);
	  ldisp (&ss);
	  sub = sub1;
	  osub = formbb (zero, sub, true, false, ' ');
	  oadd (&brnch, &osub);
	  ldisp (&sub);
	}
      unirr = W11->next;
    }
  ldisp (&zero);
  ldisp (&undum);
  if (test)
    R245 = unso3brnch (brnch, n);
  else
    R245 = gmodify (brnch, currgrp.A[1 - 1]);
  odisp (&brnch);
  return R245;
}

bool
get_ln (void)
{
  register bool R247;

  {
    register ocharptr W12 = &(*(*G259_lambda));

    (*G273_l) = W12->val.A[1] - W12->val.A[2];
    (*G275_n) = W12->val.A[2];
    if ((*G273_l) < (*G275_n))
      R247 = false;
    else
      R247 = true;
  }
  return R247;
}

void
so7_to_lambda (void)
{
  ocharptr lambda_currentb, new_lambda;
  register int ib;

  (*G261_so7_currentb) = (*G263_so7_head);
  (*G271_lambda1) =
    (ocharptr) malloc ((unsigned) (sizeof (*(*G271_lambda1))));
  rcount = rcount + 1;
  lambda_currentb = (*G271_lambda1);
  while ((*G261_so7_currentb) != NULL)
    {
      {
	register ocharptr W13 = &(*lambda_currentb);

	W13->mult = (*G261_so7_currentb)->mult;
	W13->C6_double = false;
	for (ib = 0; ib <= maxdim; ib++)
	  {
	    W13->val.A[ib] = 0;
	  }
	W13->val.A[1] = (*G261_so7_currentb)->l;
	W13->val.A[2] = (*G261_so7_currentb)->n;
	W13->lab = ' ';
	W13->conlab = ' ';
	W13->spin = false;
	new_lambda = (ocharptr) malloc ((unsigned) (sizeof (*new_lambda)));
	rcount = rcount + 1;
	W13->next = new_lambda;
      }
      (*G261_so7_currentb) = (*G261_so7_currentb)->next;
      if ((*G261_so7_currentb) != NULL)
	lambda_currentb = lambda_currentb->next;
      else
	{
	  free (new_lambda);
	  rcount = rcount - 1;
	  lambda_currentb->next = NULL;
	}
    }
}

void
init (void)
{
  (*G261_so7_currentb) =
    (representation *) malloc ((unsigned) (sizeof (*(*G261_so7_currentb))));
  rcount = rcount + 1;
  {
    register representation *W16 = &(*(*G261_so7_currentb));

    W16->prev = NULL;
    W16->next = NULL;
    W16->l = 0;
    W16->n = 0;
    W16->mult = 0;
  }
  (*G263_so7_head) = (*G261_so7_currentb);
  (*G265_g2_currentb) =
    (representation *) malloc ((unsigned) (sizeof (*(*G265_g2_currentb))));
  rcount = rcount + 1;
  {
    register representation *W17 = &(*(*G265_g2_currentb));

    W17->prev = NULL;
    W17->next = NULL;
    W17->l = 0;
    W17->n = 0;
    W17->mult = (*G259_lambda)->mult;
  }
  (*G267_g2_head) = (*G265_g2_currentb);
}

void
dispose_heap (void)
{
  (*G261_so7_currentb) = (*G263_so7_head);
  while ((*G261_so7_currentb) != NULL)
    {
      (*G269_temp_ptr) = (*G261_so7_currentb);
      (*G261_so7_currentb) = (*G261_so7_currentb)->next;
      free ((*G269_temp_ptr));
      rcount = rcount - 1;
    }
  (*G265_g2_currentb) = (*G267_g2_head);
  while ((*G265_g2_currentb) != NULL)
    {
      (*G269_temp_ptr) = (*G265_g2_currentb);
      (*G265_g2_currentb) = (*G265_g2_currentb)->next;
      free ((*G269_temp_ptr));
      rcount = rcount - 1;
    }
}

void
insert_node (rep_ptr * head, rep_ptr * currentb, byte l, byte n, int mult)
{
  rep_ptr new_node;

  new_node = (representation *) malloc ((unsigned) (sizeof (*new_node)));
  rcount = rcount + 1;
  if ((*head) != NULL)
    {
      if ((*currentb)->prev != NULL)
	{
	  new_node->next = (*currentb);
	  new_node->prev = (*currentb)->prev;
	  (*currentb)->prev->next = new_node;
	  (*currentb)->prev = new_node;
	}
      else
	{
	  new_node->prev = NULL;
	  new_node->next = (*currentb);
	  (*currentb)->prev = new_node;
	  (*head) = new_node;
	}
    }
  else
    {
      new_node->prev = NULL;
      new_node->next = NULL;
      (*head) = new_node;
    }
  new_node->l = l;
  new_node->n = n;
  new_node->mult = mult;
  (*currentb) = new_node;
}

void
delete_node (rep_ptr * head, rep_ptr * currentb)
{
  rep_ptr temp;

  temp = (*currentb);
  if ((*currentb) != NULL)
    {
      if (((*currentb)->prev == NULL) && ((*currentb)->next != NULL))
	{
	  (*head) = (*currentb)->next;
	  (*currentb) = (*head);
	  (*currentb)->prev = NULL;
	  free (temp);
	  rcount = rcount - 1;
	}
      else if (((*currentb)->next != NULL) && ((*currentb)->prev != NULL))
	{
	  (*currentb)->prev->next = (*currentb)->next;
	  (*currentb)->next->prev = (*currentb)->prev;
	  (*currentb) = (*currentb)->prev;
	  free (temp);
	  rcount = rcount - 1;
	}
      else if ((*currentb)->prev != NULL)
	{
	  (*currentb)->prev->next = NULL;
	  (*currentb) = (*currentb)->prev;
	  free (temp);
	  rcount = rcount - 1;
	}
      else if (((*currentb)->prev == NULL) && ((*currentb)->next == NULL))
	{
	  free ((*currentb));
	  rcount = rcount - 1;
	  (*currentb) = NULL;
	  (*head) = NULL;
	}
    }
}

bool
found (rep_ptr * head, rep_ptr * currentb, int l, int n)
{
  register bool R248=false;
  bool finish;

  (*currentb) = (*head);
  finish = false;
  if ((*head) != NULL)
    {
      do
	{
	  if ((*currentb)->l > l)
	    {
	      if ((*currentb)->next != NULL)
		(*currentb) = (*currentb)->next;
	      else
		finish = true;
	    }
	  else
	    finish = true;
	}
      while (!(finish));
      if ((*currentb)->l == l)
	{
	  finish = false;
	  do
	    {
	      if (((*currentb)->n > n) && ((*currentb)->l == l))
		{
		  if ((*currentb)->next != NULL)
		    (*currentb) = (*currentb)->next;
		  else
		    {
		      finish = true;
		      R248 = false;
		    }
		}
	      else if (((*currentb)->n == n) && ((*currentb)->l == l))
		{
		  finish = true;
		  R248 = true;
		}
	      else
		{
		  finish = true;
		  R248 = false;
		}
	    }
	  while (!(finish));
	}
      else
	R248 = false;
    }
  else
    R248 = false;
  return R248;
}

void
first_term (int l, int n, int mult)
{
  register int u;

  if ((l - n) / ((double) (n + 1)) >= 0)
    {
      (*G269_temp_ptr) = (*G263_so7_head);
      {
	int lastStep = Trunc ((l - n) / ((double) (n + 1)));
	for (u = 0; u <= lastStep; u++)
	  {
	    if (found
		(&(*G269_temp_ptr), &(*G261_so7_currentb), l - (n + 1) * u,
		 n))
	      (*G261_so7_currentb)->mult = (*G261_so7_currentb)->mult + mult;
	    else
	      insert_node (&(*G263_so7_head), &(*G261_so7_currentb),
			   l - (n + 1) * u, n, mult);
	    if ((*G261_so7_currentb)->mult == 0)
	      delete_node (&(*G263_so7_head), &(*G261_so7_currentb));
	    (*G269_temp_ptr) = (*G261_so7_currentb);
	    if (found
		(&(*G269_temp_ptr), &(*G261_so7_currentb), l - (n + 1) * u,
		 0))
	      (*G261_so7_currentb)->mult = (*G261_so7_currentb)->mult - mult;
	    else
	      insert_node (&(*G263_so7_head), &(*G261_so7_currentb),
			   l - (n + 1) * u, 0, -mult);
	    if ((*G261_so7_currentb)->mult == 0)
	      delete_node (&(*G263_so7_head), &(*G261_so7_currentb));
	    (*G269_temp_ptr) = (*G261_so7_currentb);
	  }
      }
    }
}

void
second_term (int l, int n, int mult)
{
  register int u;

  if ((l - n - 1) / ((double) (n + 1)) >= 0)
    {
      (*G269_temp_ptr) = (*G263_so7_head);
      {
	int lastStep = Trunc ((l - n - 1) / ((double) (n + 1)));
	for (u = 0; u <= lastStep; u++)
	  {
	    if (found
		(&(*G269_temp_ptr), &(*G261_so7_currentb),
		 l - 1 - (n + 1) * u, n))
	      (*G261_so7_currentb)->mult = (*G261_so7_currentb)->mult - mult;
	    else
	      insert_node (&(*G263_so7_head), &(*G261_so7_currentb),
			   l - 1 - (n + 1) * u, n, -mult);
	    if ((*G261_so7_currentb)->mult == 0)
	      delete_node (&(*G263_so7_head), &(*G261_so7_currentb));
	    (*G269_temp_ptr) = (*G261_so7_currentb);
	    if (found
		(&(*G269_temp_ptr), &(*G261_so7_currentb),
		 l - 1 - (n + 1) * u, 0))
	      (*G261_so7_currentb)->mult = (*G261_so7_currentb)->mult + mult;
	    else
	      insert_node (&(*G263_so7_head), &(*G261_so7_currentb),
			   l - 1 - (n + 1) * u, 0, mult);
	    if ((*G261_so7_currentb)->mult == 0)
	      delete_node (&(*G263_so7_head), &(*G261_so7_currentb));
	    (*G269_temp_ptr) = (*G261_so7_currentb);
	  }
      }
    }
}

void
third_term (int l, int n, int mult)
{
  register int u;
  register int a;
  for (a = 1; a <= (n - 1); a++)
    {
      if ((l - n) / ((double) (n + 1)) >= 0)
	{
	  (*G269_temp_ptr) = (*G267_g2_head);
	  {
	    int lastStep = Trunc ((l - n) / ((double) (n + 1)));
	    for (u = 0; u <= lastStep; u++)
	      {
		if (found
		    (&(*G269_temp_ptr), &(*G265_g2_currentb),
		     l - (n + 1) * u, a))
		  (*G265_g2_currentb)->mult =
		    (*G265_g2_currentb)->mult - mult;
		else
		  insert_node (&(*G267_g2_head), &(*G265_g2_currentb),
			       l - (n + 1) * u, a, -mult);
		if ((*G265_g2_currentb)->mult == 0)
		  delete_node (&(*G267_g2_head), &(*G265_g2_currentb));
		(*G269_temp_ptr) = (*G265_g2_currentb);
	      }
	  }
	}
    }
}

void
fourth_term (int l, int n, int mult)
{
  register int u;
  register int a;
  for (a = 1; a <= (n - 1); a++)
    {
      if (((l - a) / ((double) (n + 1)) - 1) >= 0)
	{
	  (*G269_temp_ptr) = (*G267_g2_head);
	  {
	    int lastStep = Trunc ((l - a) / ((double) (n + 1)) - 1);
	    for (u = 0; u <= lastStep; u++)
	      {
		if (found
		    (&(*G269_temp_ptr), &(*G265_g2_currentb),
		     l - (a + 1) - (n + 1) * u, a))
		  (*G265_g2_currentb)->mult =
		    (*G265_g2_currentb)->mult + mult;
		else
		  insert_node (&(*G267_g2_head), &(*G265_g2_currentb),
			       l - (a + 1) - (n + 1) * u, a, mult);
		if ((*G265_g2_currentb)->mult == 0)
		  delete_node (&(*G267_g2_head), &(*G265_g2_currentb));
		(*G269_temp_ptr) = (*G265_g2_currentb);
	      }
	  }
	}
    }
}

void
fifth_term (int l, int mult)
{
  insert_node (&(*G263_so7_head), &(*G261_so7_currentb), l, 0, mult);
}

void
P62_branch (int l, int n, int mult)
{
  first_term (l, n, mult);
  second_term (l, n, mult);
  third_term (l, n, mult);
  fourth_term (l, n, mult);
}

ocharptr
g2_so7a (ocharptr lambda)
{
  register ocharptr R246;
  rep_ptr so7_currentb, so7_head;
  rep_ptr g2_currentb, g2_head;
  rep_ptr temp_ptr;
  ocharptr lambda1g;
  byte l, n;
  ocharptr *F260;
  rep_ptr *F262;
  rep_ptr *F264;
  rep_ptr *F266;
  rep_ptr *F268;
  rep_ptr *F270;
  ocharptr *F272;
  byte *F274;
  byte *F276;

  F276 = G275_n;
  G275_n = &n;
  F274 = G273_l;
  G273_l = &l;
  F272 = G271_lambda1;
  G271_lambda1 = &lambda1g;
  F270 = G269_temp_ptr;
  G269_temp_ptr = &temp_ptr;
  F268 = G267_g2_head;
  G267_g2_head = &g2_head;
  F266 = G265_g2_currentb;
  G265_g2_currentb = &g2_currentb;
  F264 = G263_so7_head;
  G263_so7_head = &so7_head;
  F262 = G261_so7_currentb;
  G261_so7_currentb = &so7_currentb;
  F260 = G259_lambda;
  G259_lambda = &lambda;
  {
    /*register ocharptr W30 = &(*(*G259_lambda)); *//*12/12/95 */

    if (get_ln () == false)
      {
	R246 = NULL;
      }
  }
  init ();
  if ((*G275_n) != 0)
    {
      P62_branch ((*G273_l), (*G275_n), (*G259_lambda)->mult);
      while ((*G267_g2_head) != NULL)
	{
	  (*G265_g2_currentb) = (*G267_g2_head);
	  {
	    register representation *W31 = &(*(*G265_g2_currentb));

	    P62_branch (W31->l, W31->n, W31->mult);
	  }
	  delete_node (&(*G267_g2_head), &(*G267_g2_head));
	}
    }
  else
    fifth_term ((*G273_l), (*G259_lambda)->mult);
  so7_to_lambda ();
  dispose_heap ();
  R246 = (*G271_lambda1);
  G259_lambda = F260;
  G261_so7_currentb = F262;
  G263_so7_head = F264;
  G265_g2_currentb = F266;
  G267_g2_head = F268;
  G269_temp_ptr = F270;
  G271_lambda1 = F272;
  G273_l = F274;
  G275_n = F276;
  return R246;
}

ocharptr
g2_so7 (ocharptr lambda)
{
  register ocharptr R249;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (lambda != NULL)
    {
      register ocharptr W32 = &(*lambda);

      subbrnch = g2_so7a (lambda);
      oadd (&brnch, &subbrnch);
      lambda = W32->next;
    }
  osort (&brnch, true);
  R249 = brnch;
  return R249;
}

ocharptr
g2so3brnch (ocharptr lambda)
{
  register ocharptr R250;
  ocharptr brnch, subbrnch;

  subbrnch = g2_so7 (lambda);
  brnch = so2k1so3brnch ('c', subbrnch, 7, true);
  R250 = brnch;
  odisp (&subbrnch);
  return R250;
}

ocharptr
so7g2brnch (ocharptr lambda)
{
  register ocharptr R251;
  ocharptr list, sublist, indiv;
  int imax, jmax, kmax, lmax;
  register int p;
  register int jb;
  register int k;
  register int ib;

  list = NULL;
  while (lambda != NULL)
    {
      register ocharptr W33 = &(*lambda);

      sublist = NULL;
      imax = 2 * W33->val.A[3];
      if (lambda->spin)
	imax = imax + 1;
      kmax = W33->val.A[2] - W33->val.A[3];
      jmax = W33->val.A[1] - W33->val.A[2];
      lmax = W33->val.A[1] - W33->val.A[3];
      for (ib = 0; ib <= imax; ib++)
	{
	  for (k = 0; k <= kmax; k++)
	    {
	      {
		int lastStep = MIN (ib + k, jmax);
		for (jb = 0; jb <= lastStep; jb++)
		  {
		    cnu (&indiv);
		    for (p = 0; p <= maxdim; p++)
		      {
			indiv->val.A[p] = 0;
		      }
		    indiv->val.A[1] = kmax + lmax + ib - k;
		    indiv->val.A[2] = kmax + jb - k;
		    indiv->mult = W33->mult;
		    indiv->spin = false;
		    indiv->C6_double = false;
		    indiv->lab = ' ';
		    indiv->next = sublist;
		    sublist = indiv;
		  }
	      }
	    }
	}
      lambda = W33->next;
      oadd (&list, &sublist);
    }
  g2modify (&list, 0);
  osort (&list, false);
  R251 = list;
  return R251;
}

void
u1dtrmne (ocharptr * list, int n, int m, int n2, int m2, int ss)
{
  int u1value, valm;
  register int ib;
  ocharptr toplist;

  toplist = (*list);
  while ((*list) != NULL)
    {
      register ocharptr W42 = &(*(*list));

      if (ss == 1)
	{
	  u1value = m * wtfrm (&W42->conval) - n * wtfrm (&W42->val);
	  if (((n == 1) && (m > 1)) || ((n > 1) && (m == 1)))
	    u1value = -u1value;
	  valm = W42->val.A[m];

	  for (ib = m; ib >= 2; ib--)
	    {
	      W42->val.A[ib] = W42->val.A[ib - 1] - valm;
	    }
	  W42->val.A[1] = u1value;
	  for (ib = 1; ib <= n; ib++)
	    {
	      W42->conval.A[ib] = W42->conval.A[ib] - W42->conval.A[n];
	    }
	}
      else if (ss == 2)
	{
	  u1value = m * wtfrm (&W42->conval) + n * wtfrm (&W42->val);
	  valm = W42->val.A[m];

	  for (ib = m; ib >= 2; ib--)
	    {
	      W42->val.A[ib] = W42->val.A[ib - 1] - valm;
	    }
	  W42->val.A[1] = u1value;
	  for (ib = 1; ib <= n; ib++)
	    {
	      W42->conval.A[ib] = W42->conval.A[ib] - W42->conval.A[n];
	    }
	}
      else
	{
	  u1value =
	    (n2 - m2) * wtfrm (&W42->conval) - (n - m) * wtfrm (&W42->val);

	  for (ib = maxdim; ib >= 2; ib--)
	    {
	      W42->val.A[ib] = W42->val.A[ib - 1];
	    }
	  W42->val.A[1] = u1value;
	}
      (*list) = W42->next;
    }
  osort (&toplist, true);
  (*list) = toplist;
}

ocharptr
unospnbrnch (char br, ocharptr unirr, int n)
{
  register ocharptr R252;
  termptr ss, sub1, sub2, sub, zero, undum;
  ocharptr brnch, osub;
  register int ib;

  snu (&zero);
  {
    register termptr W53 = &(*zero);
    for (ib = 0; ib <= maxdim; ib++)
      {
	W53->val.A[ib] = 0;
      }
    W53->next = NULL;
  }
  snu (&undum);
  undum->mult = 1;
  undum->next = NULL;
  brnch = NULL;
  while (unirr != NULL)
    {
      register ocharptr W56 = &(*unirr);

      zero->mult = W56->mult;
      undum->val = W56->val;
      ss = seriesx (br, -1, W56->val);
      sub1 = lskew (undum, ss);
      ldisp (&ss);
      if (W56->C6_double)
	{
	  undum->val = W56->conval;
	  ss = seriesx (br, -1, W56->conval);
	  sub2 = lskew (undum, ss);
	  if ((br != 'm'))
	    sub = louter (sub1, sub2);
	  else
	    osub = formbb (sub1, sub2, false, false, ' ');
	  ldisp (&sub1);
	  ldisp (&sub2);
	}
      else
	sub = sub1;
      if (((br != 'm') || ((br == 'm') && (!W56->C6_double))))
	{
	  osub = formbb (zero, sub, true, false, ' ');
	  ldisp (&sub);
	}
      oadd (&brnch, &osub);
      osort (&brnch, true);
      unirr = W56->next;
    }

  dispsfn (&zero);
  dispsfn (&undum);
  if (br == 'd')
    {
       /**/ gmodify (brnch, currgrp.A[1 - 1]);
      omodify (&brnch, n);
       /**/
	/* gumodify(&brnch, n); *//*21 3 1999 */
    }

  if (br == 'b')
    spmodify (&brnch, n);
  if (br == 'm')
    gumodify (&brnch, n);
  R252 = brnch;
  return R252;
}

ocharptr
supqu1br (ocharptr irrep, int n, int m, int n2, int m2, int ss)
{
  ocharptr odum, co, cu, list, sublist;
  termptr zeta, tzeta, dum, ort, unnit, ks;
  frame zero;
  register int ib;
  bool spinn;
  /*char  ch; *//*12/12/95 */

  list = NULL;
  for (ib = 0; ib <= maxdim; ib++)
    {
      zero.A[ib] = 0;
    }
  snu (&dum);
  {
    register termptr W59 = &(*dum);

    W59->val = zero;
    W59->mult = 1;
    W59->next = NULL;
  }
  cnu (&odum);
  {
    register ocharptr W60 = &(*odum);

    W60->mult = 1;
    W60->val = zero;
    W60->spin = false;
    W60->C6_double = false;
    W60->lab = ' ';
    W60->next = NULL;
  }
  while (irrep != NULL)
    {
      spinn = irrep->spin;
      zeta = skcompat (irrep->val);
      if ((ss == 0))
	schur_restrict (&zeta, -m, 'w');
      while (zeta != NULL)
	{
	  tzeta = zeta;
	  zeta = zeta->next;
	  tzeta->next = NULL;
	  unnit = skew (irrep->val, tzeta->val);
	  if (ss == 2)
	    ort = lconjgte (tzeta);
	  else
	    ort = tzeta;
	  if (ss == 5)
	    {
	      ks = seriesx ('b', -1, ort->val);
	      ort = lskew (ort, ks);
	      ldisp (&ks);
	    }
	  if (ss == 6)
	    {
	      ks = seriesx ('d', -1, ort->val);
	      ort = lskew (ort, ks);
	      ldisp (&ks);
	    }
	  if (((ss == 5) || (ss == 6)))
	    dispsfn (&tzeta);
	  cu = formbb (dum, unnit, true, spinn, ' ');
	  co = formbb (dum, ort, true, spinn, ' ');
	  ldisp (&unnit);
	  ldisp (&ort);
	  if (ss != 0)
	    {
	      if (((ss != 3) && (ss != 5) && (ss != 6)))
		{
		  gumodify (&cu, n);
		  gumodify (&co, m);
		}
	      if ((ss == 3))
		{
		  unmmodify (&cu, n, m);
		  unmmodify (&co, n2, m2);
		}
	      if ((ss == 5))
		{
		  spmodify (&cu, n);
		  spmodify (&co, m);
		}
	      if ((ss == 6))
		{
		  omodify (&cu, n);
		  omodify (&co, m);
		}
	    }
	  else
	    {
	      snselect (&cu, n, qspecial);
	      snselect (&co, m, qspecial);
	    }
	  sublist = oformbb (cu, co, false, spinn);
	  odisp (&cu);
	  odisp (&co);
	  oadd (&list, &sublist);
	}
      irrep = irrep->next;
    }
  dispsfn (&dum);
  dispchr (&odum);
  osort (&list, false);
  if (((ss != 0) && (ss != 4) && (ss != 5) && (ss != 6)))
    u1dtrmne (&list, n, m, n2, m2, ss);
  
  return list;
}



ocharptr
unsubgrbr (ocharptr irrep, int n, int ss, char br)
{
  register ocharptr R254;
  ocharptr list, sublist;
  termptr temp, tlist, nlist, xlist, sx, ttemp, dummy, dumm, dums,
    zeta, tzeta, unnit;
  int jb;
  jb = n / 2;
  temp = NULL;

  while (irrep != NULL)
    {
      register ocharptr W61 = &(*irrep);

      snu (&tlist);
      tlist->mult = W61->mult;
      tlist->val = W61->val;
      sx = seriesx (br, -1, W61->val);
      nlist = lskew (tlist, sx);
      ldisp (&sx);
      dispsfn (&tlist);
      if ((ss == 0))
	dums = nlist;

      if ((ss == 1))
	{
	  {
	    dumm = nlist;
	  }			/*31/12/97 */
	  ttemp = NULL;
	  while (nlist != NULL)
	    {
	      register termptr W62 = &(*nlist);

	      snu (&tlist);
	      tlist->mult = W62->mult;
	      tlist->val = W62->val;
	      sx = seriesx ('d', -1, W62->val);
	      xlist = lskew (tlist, sx);

	      ldisp (&sx);
	      dispsfn (&tlist);

	      add (&ttemp, &xlist);

	      nlist = W62->next;
	    }
	  {
	    ldisp (&dumm);
	  }			/*31/12/97 */
	  nlist = ttemp;
	}

      add (&temp, &nlist);
      irrep = W61->next;
    }

  list = NULL;
  while (temp != NULL)
    {
      register termptr W63 = &(*temp);

      zeta = skcompat (temp->val);
      dummy = zeta;

      while (zeta != NULL)
	{
	  register termptr W64 = &(*zeta);

	  snu (&tzeta);
	  tzeta->mult = W63->mult;	/*20/4/98 */
	  tzeta->val = W64->val;
	  tzeta->next = NULL;
	  unnit = skew (temp->val, tzeta->val);

	  sublist = formbb (tzeta, unnit, true, false, ' ');

	  ldisp (&unnit);	/*31/12/97 */
	  dispsfn (&tzeta);
	  oadd (&list, &sublist);
	  zeta = zeta->next;

	}

      ldisp (&dummy);		/*31/12/97 */
      temp = W63->next;
    }

  if ((ss == 1))
    {
      ldisp (&ttemp);
    }

  if ((ss == 0))
    ldisp (&dums);
  gumodify (&list, jb);

  R254 = list;
  return R254;
}

ocharptr
onmonombrnch (ocharptr irrep, int n, int m, int ss)
{
  register ocharptr R255;
  termptr zeta, tzeta, dum, left, right, dleft, dright, series;
  ocharptr branch, subbrnch, odum, cleft, cright, dummy,
    extra, diffbrnch, dsub, ddiff;
  frame zero;
  int lamwt, mult1, mult2;
  register int ib;
  bool same, neg;
  char labl;

  branch = NULL;
  zero = nolls;
  while (irrep != NULL)
    {
      register ocharptr W65 = &(*irrep);

      snu (&dum);
      {
	register termptr W66 = &(*dum);

	W66->val = zero;
	W66->mult = 1;
	W66->next = NULL;
      }
      cnu (&odum);
      {
	register ocharptr W67 = &(*odum);

	W67->mult = 1;
	W67->val = zero;
	W67->spin = false;
	W67->C6_double = false;
	W67->lab = ' ';
	W67->next = NULL;
      }
      labl = W65->lab;
      if ((bool) ((n + m) & 1) || ((bool) ((n) & 1) && (bool) ((m) & 1)))
	W65->lab = ' ';
      if (((W65->lab == '+') || (W65->lab == '-'))
	  && !((bool) ((n) & 1) || (bool) ((m) & 1)))
	{
	  cnu (&dummy);
	  (*dummy) = (*irrep);
	  dummy->lab = ' ';
	  dummy->next = NULL;
	  subbrnch = onmonombrnch (dummy, n, m, ss);
	  if (W65->spin)
	    {
	      lamwt = wtfrm (&W65->val);
	      dsub = subbrnch;
	      while (dsub != NULL)
		{
		  same = (bool) (dsub->lab == dsub->conlab);
		  neg =
		    (bool) ((wtfrm (&dsub->val) - wtfrm (&dsub->conval) -
			     lamwt) & 1);
		  if ((neg && same) || (!neg && !same))
		    {
		      if (W65->lab == '+')
			dsub->mult = 0;
		    }
		  else if (W65->lab == '-')
		    dsub->mult = 0;
		  dsub = dsub->next;
		}
	    }
	  else
	    {
	      dummy->lab = '"';
	      diffbrnch = onmonombrnch (dummy, n, m, ss);
	      dsub = subbrnch;
	      ddiff = diffbrnch;
	      if (diffbrnch != NULL)
		{
		  while (ddiff->next != NULL)
		    ddiff = ddiff->next;
		  cnu (&extra);
		  (*extra) = (*ddiff);
		  extra->val.A[1] = -1;
		  ddiff->next = extra;
		}
	      else
		{
		  cnu (&diffbrnch);
		  diffbrnch->val.A[1] = -1;
		  diffbrnch->next = NULL;
		}
	      ddiff = diffbrnch;
	      while (dsub != NULL)
		if (dsub->val.A[m / 2] != 0)
		  if (dsub->conval.A[n / 2] != 0)
		    if ((testord (&ddiff->val, &dsub->val) == EQUAL)
			&& (testord (&ddiff->conval, &dsub->conval) == EQUAL))
		      if (dsub->mult == ddiff->mult)
			{
			  for (ib = 1; ib <= 4; ib++)
			    {
			      same = (bool) (dsub->lab == dsub->conlab);
			      if (((W65->lab == '+') && !same)
				  || ((W65->lab == '-') && same))
				dsub->mult = 0;
			      dsub = dsub->next;
			    }
			  ddiff = ddiff->next;
			}
		      else
			{
			  mult1 = (dsub->mult + ddiff->mult) / 2;
			  mult2 = (dsub->mult - ddiff->mult) / 2;
			  for (ib = 1; ib <= 4; ib++)
			    {
			      same = (bool) (dsub->lab == dsub->conlab);
			      if (((W65->lab == '+') && !same)
				  || ((W65->lab == '-') && same))
				dsub->mult = mult2;
			      else
				dsub->mult = mult1;
			      dsub = dsub->next;
			    }
			  ddiff = ddiff->next;
			}
		    else
		      for (ib = 1; ib <= 4; ib++)
			{
			  same = (bool) (dsub->lab == dsub->conlab);
			  if (((W65->lab == '+') && !same)
			      || ((W65->lab == '-') && same))
			    dsub->mult = 0;
			  else
			    dsub->mult = dsub->mult / 2;
			  dsub = dsub->next;
			}
		  else
		    {
		      dsub->mult = dsub->mult / 2;
		      dsub = dsub->next;
		      dsub->mult = dsub->mult / 2;
		      dsub = dsub->next;
		    }
		else
		  {
		    if (dsub->conval.A[n / 2] != 0)
		      {
			dsub->mult = dsub->mult / 2;
			dsub = dsub->next;
			dsub->mult = dsub->mult / 2;
			dsub = dsub->next;
		      }
		    else
		      {
			dsub->mult = dsub->mult / 2;
			dsub = dsub->next;
		      }
		  }
	      odisp (&diffbrnch);
	    }
	  dispchr (&dummy);
	}
      else
	{
	  if (W65->lab == '"')
	    {
	      int lastStep = (n + m) / 2;
	      for (ib = 1; ib <= lastStep; ib++)
		{
		  W65->val.A[ib] = W65->val.A[ib] - 1;
		}
	    }
	  subbrnch = NULL;
	  zeta = skcompat (W65->val);
	  while (zeta != NULL)
	    {
	      tzeta = zeta;
	      zeta = zeta->next;
	      tzeta->next = NULL;
	      if (W65->spin)
		series = skcompat (tzeta->val);
	      else if (W65->lab == '"')
		series = seriesx ('b', -1, tzeta->val);
	      else
		series = seriesx ('d', -1, tzeta->val);
	      tzeta->mult = tzeta->mult * W65->mult;
	      left = lskew (tzeta, series);
	      ldisp (&series);
	      right = skew (W65->val, tzeta->val);
	      dispsfn (&tzeta);
	      if (W65->lab == '"')
		{
		  left = lconjgte (left);
		  right = lconjgte (right);
		  dleft = left;
		  while (dleft != NULL)
		    {
		      register termptr W76 = &(*dleft);

		      for (ib = maxdim - 1; ib >= 1; ib--)
			{
			  W76->val.A[ib + 1] = W76->val.A[ib];
			}
		      W76->val.A[1] = n / 2;
		      dleft = W76->next;
		    }
		  dright = right;
		  while (dright != NULL)
		    {
		      register termptr W79 = &(*dright);

		      for (ib = maxdim - 1; ib >= 1; ib--)
			{
			  W79->val.A[ib + 1] = W79->val.A[ib];
			}
		      W79->val.A[1] = m / 2;
		      dright = W79->next;
		    }
		  stndise (&left);
		  stndise (&right);
		  left = lconjgte (left);
		  right = lconjgte (right);
		}
	      if (ss == 2)
		left = lconjgte (left);
	      cleft = formbb (dum, left, true, W65->spin, W65->lab);
	      if (ss == 1)
		omodify (&cleft, n);
	      else
		spmodify (&cleft, n);
	      if (!(bool) ((n) & 1) && !(W65->lab == '"') && (ss == 1))
		so2nexp (&cleft, n);
	      cright = formbb (dum, right, true, W65->spin, W65->lab);
	      omodify (&cright, m);
	      if (!(bool) ((m) & 1) && !(W65->lab == '"') && (ss == 1))
		so2nexp (&cright, m);
	      dsub = oformbb (cleft, cright, false, W65->spin);
	      ldisp (&left);
	      ldisp (&right);
	      odisp (&cleft);
	      odisp (&cright);
	      oadd (&subbrnch, &dsub);
	    }
	}
      W65->lab = labl;
      if ((bool) ((n) & 1) && (bool) ((m) & 1))
	{
	  dsub = subbrnch;
	  if (W65->spin && (W65->lab == ' '))
	    while (dsub != NULL)
	      {
		dsub->mult = dsub->mult * 2;
		dsub = dsub->next;
	      }
	  else if (!W65->spin && ((W65->lab == '+') || (W65->lab == '-')))
	    {
	      osort (&subbrnch, true);
	      dsub = subbrnch;
	      while (dsub != NULL)
		{
		  dsub->mult = dsub->mult / 2;
		  dsub = dsub->next;
		}
	    }
	}
      oadd (&branch, &subbrnch);
      dispsfn (&dum);
      dispchr (&odum);
      if (W65->lab == '"')
	{
	  int lastStep = (n + m) / 2;
	  for (ib = 1; ib <= lastStep; ib++)
	    {
	      W65->val.A[ib] = W65->val.A[ib] + 1;
	    }
	}
      if (ss == 2)
	gsort (&branch, true, prod);
      irrep = irrep->next;
    }
  gsort (&branch, true, prod);
  R255 = branch;
  return R255;
}


ocharptr
sonunbr1 (ocharptr lambda, int n)
{
  register ocharptr R257;
  ocharptr list, sublist;
  termptr slam, xi, xit, temp1, temp2, temp, bb;
  int u1val, k;
  register int ib;

  list = NULL;
  k = n / 2;
  snu (&slam);
  {
    register termptr W84 = &(*slam);

    W84->val = lambda->val;
    W84->mult = lambda->mult;
    W84->next = NULL;
  }
  xi = skcompat (slam->val);
  while (xi != NULL)
    {
      temp1 = pskew (slam, xi);
      xit = xi;
      xi = xi->next;
      xit->next = NULL;
      while (temp1 != NULL)
	{
	  if ((bool) ((n) & 1))
	    bb = seriesx ('f', -1, temp1->val);
	  else
	    bb = seriesx ('b', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  ldisp (&bb);
	  dispsfn (&temp);
	  sublist = formbb (xit, temp2, false, false, ' ');
	  gumodify (&sublist, k);
	  oadd (&list, &sublist);
	  ldisp (&temp2);
	}
      dispsfn (&xit);
    }
  sublist = list;
  while (sublist != NULL)
    {
      register ocharptr W85 = &(*sublist);

      if (W85->C6_double)
	{
	  u1val = 2 * (wtfrm (&W85->val) - wtfrm (&W85->conval));
	  for (ib = 1; ib <= k; ib++)
	    {
	      W85->val.A[ib] = W85->val.A[ib] + W85->conval.A[1] - W85->
		conval.A[k - ib + 1];
	    }
	  W85->C6_double = false;
	}
      else
	{
	  u1val = 2 * wtfrm (&W85->val);
	  for (ib = 1; ib <= k; ib++)
	    {
	      W85->val.A[ib] = W85->val.A[ib] - W85->val.A[k];
	    }
	}

      for (ib = maxdim - 1; ib >= 1; ib--)
	{
	  W85->val.A[ib + 1] = W85->val.A[ib];
	}
      W85->val.A[1] = u1val;
      sublist = W85->next;
    }
  dispsfn (&slam);
  R257 = list;
  return R257;
}

ocharptr
sonunbr2 (ocharptr lambda, int n)
{
  register ocharptr R258;
  ocharptr list, sub1;
  termptr slam, xi, xit, temp1, temp2, temp, bb, temp3, ser;
  frame lfrm;
  int u1val, limitb, k;
  register int ib;
  char labyl = '+';

  lfrm = nolls;
  list = NULL;
  k = n / 2;
  snu (&slam);
  {
    register termptr W92 = &(*slam);

    if (lambda->spin || (lambda->lab == ' '))
      W92->val = lambda->val;
    else
      for (ib = 0; ib <= maxdim; ib++)
	{
	  W92->val.A[ib] = MAX (0, lambda->val.A[ib] - 1);
	}
    W92->mult = lambda->mult;
    W92->next = NULL;
  }
  if (lambda->lab != ' ')
    {
      if ((bool) ((n) & 1))
	if (lambda->lab == '+')
	  labyl = '-';
	else
	  labyl = '+';
      else if (lambda->lab == '+')
	labyl = '+';
      else
	labyl = '-';
    }
  xi = skcompat (slam->val);
  while (xi != NULL)
    {
      temp1 = pskew (slam, xi);
      xit = xi;
      xi = xi->next;
      xit->next = NULL;
      limitb = xit->val.A[1] + n;
      if (lambda->lab == ' ')
	{
	  for (ib = 1; ib <= limitb; ib++)
	    {
	      lfrm.A[ib] = 1;
	    }
	  for (ib = limitb + 1; ib <= maxdim; ib++)
	    {
	      lfrm.A[ib] = 0;
	    }
	}
      while (temp1 != NULL)
	{
	  if ((bool) ((n) & 1))
	    bb = seriesx ('f', -1, temp1->val);
	  else
	    bb = seriesx ('b', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  if (lambda->spin)
	    if (lambda->lab == ' ')
	      ser = seriesx ('q', -1, lfrm);
	    else
	      ser = ql (labyl, limitb);
	  else
	    ser = xv (labyl, limitb);
	  temp3 = louter2 (temp2, ser, limitb);
	  ldisp (&bb);
	  ldisp (&ser);
	  dispsfn (&temp);
	  sub1 = formbb (xit, temp3, true, false, ' ');
	  gumodify (&sub1, k);
	  oadd (&list, &sub1);
	  ldisp (&temp2);
	  ldisp (&temp3);
	}
      dispsfn (&xit);
    }
  sub1 = list;
  while (sub1 != NULL)
    {
      register ocharptr W99 = &(*sub1);

      if (W99->C6_double)
	{
	  u1val = 2 * (wtfrm (&W99->val) - wtfrm (&W99->conval));
	  for (ib = 1; ib <= k; ib++)
	    {
	      W99->val.A[ib] = W99->val.A[ib] + W99->conval.A[1] - W99->
		conval.A[k - ib + 1];
	    }
	  W99->C6_double = false;
	}
      else
	{
	  u1val = 2 * wtfrm (&W99->val);
	  for (ib = 1; ib <= k; ib++)
	    {
	      W99->val.A[ib] = W99->val.A[ib] - W99->val.A[k];
	    }
	}
      if (lambda->spin)
	u1val = u1val - k;
      else
	u1val = u1val - n;

      for (ib = maxdim - 1; ib >= 1; ib--)
	{
	  W99->val.A[ib + 1] = W99->val.A[ib];
	}
      W99->val.A[1] = u1val;
      sub1 = W99->next;
    }
  dispsfn (&slam);
  R258 = list;
  return R258;
}

ocharptr
sonunbrnch (ocharptr lambda, int n)
{
  register ocharptr R256;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (lambda != NULL)
    {
      register ocharptr W106 = &(*lambda);

      if (W106->spin || (W106->lab != ' '))
	subbrnch = sonunbr2 (lambda, n);
      else
	subbrnch = sonunbr1 (lambda, n);
      oadd (&brnch, &subbrnch);
      lambda = W106->next;
    }
  osort (&brnch, true);
  R256 = brnch;
  return R256;
}

ocharptr
unsnbr1 (ocharptr onrep, int n)
{
  register ocharptr R257;
  termptr dummy, ss, sub1, temp, zero, undum;
  ocharptr brnch;
  register int ib;

  progress ();
  /* putchrc(&output,onrep,true); */
  snu (&zero);
  {
    register termptr W107 = &(*zero);
    for (ib = 0; ib <= maxdim; ib++)
      {
	W107->val.A[ib] = 0;
      }
    W107->next = NULL;
  }
  snu (&undum);
  undum->mult = 1;
  undum->next = NULL;
  temp = NULL;
  zero->mult = onrep->mult;
  undum->val = onrep->val;
  ss = seriesx ('m', -1, onrep->val);
  /*putsfn(&output,ss,true);
     putsfn(&output,undum,true); */
  sub1 = lskew (undum, ss);
  ldisp (&ss);
  dispsfn (&undum);
  ss = plethonerinner (sub1);
  ldisp (&sub1);
  dummy = ladd (temp, ss);
  temp = dummy;
  ldisp (&ss);
  ss = makeweight (n, temp);
  ldisp (&temp);
  if (onrep->lab == '#')
    ss = lconjgte (ss);
  brnch = formbb (zero, ss, true, false, ' ');
  ldisp (&ss);
  snselect (&brnch, n, qspecial);
  osort (&brnch, true);
  dispsfn (&zero);
  R257 = brnch;
  return R257;
}

ocharptr
unsnbr (ocharptr lambda, int n)
{
  register ocharptr R258;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (lambda != NULL)
    {
      register ocharptr W108 = &(*lambda);

      subbrnch = unsnbr1 (lambda, n);
      oadd (&brnch, &subbrnch);
      lambda = W108->next;
    }
  osort (&brnch, true);
  R258 = brnch;
  return R258;
}
