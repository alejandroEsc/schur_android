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
# include <stdlib.h>
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
# include "s5.h"
# include "s6.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "q1.h"
# include "q2.h"
# include "g.h"
# include "gr.h"
# include "branch.h"
# include "skew.h"
# include "label.h"
# include "s3.h"

#define INAPPROPRIATE_GROUP()  {print("Inappropriate group\n");return(NULL);}


ocharptr
tspncprod (ocharptr left, ocharptr right, int n)
{
  ocharptr temp, newlist, tempc;
  termptr munu, /* mu, */ nu, nuu, nudd, ds,
    templ, templist, nud, dummy, dummyy;
  int pr, pl, k, l,  kl, kk, ll, nnn, xk, restx;
  bool spinkl, ppp;
  char tag;
  newlist = NULL;
  nnn = n / 2;
  k = left->val.A[maxdim];
  l = right->val.A[maxdim];
  if (left->spin)
    kk = 2 * k + 1;
  else
    kk = 2 * k;
  if (right->spin)
    ll = 2 * l + 1;
  else
    ll = 2 * l;
  if ((ll > kk) || ((ll == kk) && (left->val.A[1] < right->val.A[1])))
    {
      tempc = left;
      left = right;
      right = tempc;
      xk = k;
      l = k;
      k = xk;
      xk = ll;
      ll = kk;
      kk = xk;
    }
  kl = (kk + ll) / 2;
  if ((bool) ((kk + ll) & 1))
    spinkl = true;
  else
    spinkl = false;
  pl = len (&left->val);
  pr = len (&right->val);
  //nn = MIN (nnn, ll);
  /*snu(&mu);
     mu->mult = left->mult;
     mu->val = left->val;
     mu->val.A[maxdim] = 0; */
  ds = rseries (nnn, 'd');

  modspnsfn (&ds, kk + ll);
  restx = setlimit - pr - pl;
  schur_restrict (&ds, restx, 'w');

  if ((bool) ((ll) & 1))
    tag = 'g';
  else
    tag = 'c';
  snu (&templist);
  templist->mult = right->mult;
  templist->val = right->val;
  templist->val.A[maxdim] = 0;
  templist->val.length=len(&templist->val);
  dummy = templist;
  if (ll == nnn)
    ppp = true;
  else
    ppp = false;
  nu = signseq (&dummy, ll, nnn, tag, ppp, false);

  restx = setlimit - pr;
  modspnsfn (&nu, kk + ll);
  schur_restrict (&nu, restx, 'w');
  if ((bool) ((kk) & 1))
    tag = 'g';
  else
    tag = 'c';
  snu (&templ);
  templ->mult = left->mult;
  templ->val = left->val;
  templ->val.A[maxdim] = 0;
  templ->val.length=len(&templ->val);
  dummyy = templ;
  if (kk == nnn)
    ppp = true;
  else
    ppp = false;
  nuu = signseq (&dummyy, kk, nnn, tag, ppp, false);
  modspnsfn (&nuu, kk + ll);
  restx = setlimit - pl;
  schur_restrict (&nuu, restx, 'w');
  nud = louter2 (nu, nuu, nnn);
  modspnsfn (&nud, kk + ll);
  ldisp (&nu);
  ldisp (&nuu);
  dispsfn (&dummy);
  dispsfn (&dummyy);
  munu = louter2 (ds, nud, nnn);
  schur_restrict (&munu, setlimit, 'w');
  ldisp (&ds);
  ldisp (&nud);
  modspnsfn (&munu, kk + ll);
  nudd = munu;
  while (nudd != NULL)
    {
      register termptr W4 = &(*nudd);

      cnu (&temp);
      temp->mult = W4->mult;
      temp->val = W4->val;
      temp->C6_double = true;
      temp->conval = nolls;
      temp->conval.A[1] = kl;
      temp->val.A[maxdim] = kl;
      temp->spin = spinkl;
      temp->lab = ' ';
      temp->conlab = ' ';
      oadd (&newlist, &temp);
      nudd = W4->next;
    }
  ldisp (&munu);
  sposort (&newlist, true);
  return newlist;
}

/*ocharptr onmbrnch();*/

ocharptr
onmbrnch1 (ocharptr irrep, int n, int m, int s)
{
  register ocharptr R270;
  ocharptr co, cu, list, sublist;
  termptr temp3, temp2, sub1, ss, zeta, tzeta,
    undum, dum, ort, unnit, odum, temp, temp1;
  frame zero;
  int z, tz;
  register int i;
  char ch;

  list = NULL;
  for (i = 0; i <= maxdim; i++)
    {
      zero.A[i] = 0;
    }
  while (irrep != NULL)
    {
      if (((s == 1) || (s == 2)))
	ch = 'c';
      else
	ch = 'a';
      ss = seriesx (ch, -1, irrep->val);
      snu (&undum);
      undum->mult = 1;
      undum->next = NULL;
      undum->val = irrep->val;
      z = irrep->mult;
      sub1 = lskew (undum, ss);
      ldisp (&ss);
      dispsfn (&undum);
      temp2 = sub1;
      while (sub1 != NULL)
	{
	  zeta = eqwt (wtfrm (&sub1->val));
	  while (zeta != NULL)
	    {
	      tzeta = zeta;
	      zeta = zeta->next;
	      tzeta->next = NULL;
	      ort = tzeta;
	      if (((s == 1) || (s == 3)))
		ch = 'd';
	      else
		ch = 'b';
	      ss = seriesx (ch, -1, ort->val);
	      snu (&odum);
	      odum->mult = ort->mult;
	      odum->next = NULL;
	      odum->val = ort->val;
	      temp = lskew (odum, ss);
	      ldisp (&ss);
	      dispsfn (&odum);
	      snu (&dum);
	      {
		register termptr W7 = &(*dum);

		W7->val = zero;
		W7->next = NULL;
		W7->mult = 1;
	      }
	      co = formbb (dum, temp, true, false, ' ');
	      ldisp (&temp);
	      dispsfn (&dum);
	      if (((s == 1) || (s == 3)))
		omodify (&co, m);
	      else
		spmodify (&co, m);
	      unnit = inner (sub1->val, tzeta->val);
	      ldisp (&ort);
	      temp3 = unnit;
	      while (unnit != NULL)
		{
		  snu (&dum);
		  {
		    register termptr W8 = &(*dum);

		    W8->val = zero;
		    W8->next = NULL;
		    W8->mult = 1;
		  }
		  snu (&odum);
		  odum->mult = unnit->mult;
		  odum->next = NULL;
		  odum->val = unnit->val;
		  if (((s == 2) || (s == 3)))
		    ch = 'b';
		  else
		    ch = 'd';
		  ss = seriesx (ch, -1, unnit->val);
		  temp1 = lskew (odum, ss);
		  ldisp (&ss);
		  dispsfn (&odum);
		  cu = formbb (dum, temp1, true, false, ' ');
		  ldisp (&temp1);
		  dispsfn (&dum);
		  if ((s == 1))
		    omodify (&cu, n);
		  else
		    spmodify (&cu, n);
		  sublist = oformbb (cu, co, false, false);
		  odisp (&cu);
		  tz = z * sub1->mult;
		  cu = chrcmult (tz, sublist);
		  odisp (&sublist);
		  oadd (&list, &cu);
		  unnit = unnit->next;
		}
	      odisp (&co);
	      ldisp (&temp3);
	    }
	  sub1 = sub1->next;
	}
      ldisp (&temp2);
      irrep = irrep->next;
    }
  osort (&list, false);
  R270 = list;
  return R270;
}

ocharptr
onmbrnch (ocharptr irrep, int n, int m, int s)
{
  register ocharptr R269;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (irrep != NULL)
    {
      register ocharptr W9 = &(*irrep);

      subbrnch = onmbrnch1 (irrep, n, m, s);
      oadd (&brnch, &subbrnch);
      irrep = W9->next;
    }
  osort (&brnch, true);
  R269 = brnch;
  return R269;
}

ocharptr
rcompare (ocharptr list1, ocharptr list2)
{
  register ocharptr R271;
  ocharptr temp, temp2, templist;
  int i;
  bool test;

  temp = NULL;
  while (list1 != NULL)
    {
      temp2 = list2;
      while (temp2 != NULL)
	{
	  test = true;
	  i = 1;
	  if ((list1->lab != temp2->lab) || (list1->spin != temp2->spin))
	    test = false;
	  else
	    do
	      {
		if ((list1->val.A[i] != temp2->val.A[i]))
		  test = false;
		i = i + 1;
	      }
	    while (!(((i > maxdim - 1) || (test == false))));
	  if (test == true)
	    {
	      cnu (&templist);
	      templist->mult = list1->mult * temp2->mult;
	      templist->spin = list1->spin;
	      templist->lab = list1->lab;
	      templist->val = list1->val;
	      oadd (&temp, &templist);
	    }
	  temp2 = temp2->next;
	}
      list1 = list1->next;
    }
  osort (&temp, true);
  R271 = temp;
  return R271;
}

termptr
compare (termptr list1, termptr list2)
{
  register termptr R272;
  termptr temp, temp2, templist;
  int i;
  bool test;

  temp = NULL;
  while (list1 != NULL)
    {
      temp2 = list2;
      while (temp2 != NULL)
	{
	  test = true;
	  i = 1;
	  do
	    {
	      if ((list1->val.A[i] != temp2->val.A[i]))
		test = false;
	      i = i + 1;
	    }
	  while (!(((i > maxdim - 1) || (test == false))));
	  if (test == true)
	    {
	      snu (&templist);
	      templist->mult = list1->mult * temp2->mult;
	      templist->val = list1->val;
	      add (&temp, &templist);
	    }
	  temp2 = temp2->next;
	}
      list1 = list1->next;
    }
  sort (&temp, true);
  R272 = temp;
  return R272;
}

ocharptr
kronk (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr temp = NULL;
  ocharptr kronk1;

  echo = false;
  group = grp.name;
  qspecial = false;

  chrc1 = gmodify (chrc1, grp);
  chrc2 = gmodify (chrc2, grp);
  if ((chrc1 != NULL) && (chrc2 != NULL))
    {
      qspecial = true;
	switch (grp.name)
	  {
	  case sung:
	    temp = unitprod (chrc1, chrc2, true, grp.rank);
	    break;
	  case un:
	    temp = unitprod (chrc1, chrc2, false, grp.rank);
	    break;
	  case son:
	    temp = orthprod (chrc1, chrc2, true, grp.rank);
	    break;
	  case on:
	    temp = orthprod (chrc1, chrc2, false, grp.rank);
	    break;
	  case spn:
	    temp = sympprod (chrc1, chrc2, grp.rank);
	    break;
	  case spnc:
	    if (sb_conj)
	      kronk1 = ltspncprod (chrc1, chrc2, grp.rank);
	    else
	      kronk1 = lspncprod (chrc1, chrc2, grp.rank);
	    rrestrict (&kronk1, setlimit, 'w');
	    temp = kronk1;
	    break;
	  case sonc:
	    temp = lsoncprod (chrc1, chrc2, grp.rank);
	    break;
	  case sn:
	    temp = symmprod (chrc1, chrc2, grp.rank);
	    break;
	  case g2:
	    temp = g2prod (chrc1, chrc2);
	    break;
	  case f4:
	    temp = f4prod (chrc1, chrc2);
	    break;
	  case e6:
	    temp = e6prod (chrc1, chrc2);
	    break;
	  case e7:
	    temp = e7prod (chrc1, chrc2);
	    break;
	  case e8:
	    temp = e8prod (chrc1, chrc2);
	    break;
	  case unm:
	  case sunm:
	  case ospnm:
	  case an:
	  case mp:
	  case nill:
	    inform ("not implemented", cr);
	    break;
	    /* case nill:
	       warn("group not set:no", cont);
	       inform(" action taken", cr);
	       temp = NULL;
	       break ; */
	  default:
	    Caseerror (Line);
	  }
      echo = true;
      odisp (&chrc1);
      odisp (&chrc2);
    }
  return temp;
}


prodtype
pkronk (prodtype p1, prodtype p2)
{
  register prodtype R274;
  prodtype sublist, p11, p22, entry;
  register int j;
  bool okay = false, ktest;
  for (j = 1; j <= nprod; j++)
    {
      switch ((int) (currgrp.A[j - 1].name))
	{
	case an:
	case unm:
	case sunm:
	case ospnm:
	case nill:
	case unc:
	case mp:
	  okay = false;
	  break;
	case sn:
	case un:
	case sung:
	case son:
	case on:
	case spn:
	case en:
	case g2:
	case f4:
	case e6:
	case e7:
	case e8:
	case spnc:
	case sonc:
	  okay = true;
	  break;
	default:
	  Caseerror (Line);
	}
    }
  if (okay)
    {
      p11 = p1;
      p22 = p2;
      sublist = NULL;
      while (p11 != NULL)
	{
	  while (p22 != NULL)
	    {
	      ktest = true;
	      j = 1;
	      pnu (&entry);
	      do
		{
		  entry->prods.A[j - 1] =
		    kronk (p11->prods.A[j - 1], p22->prods.A[j - 1],
			   currgrp.A[j - 1]);
		  if (entry->prods.A[j - 1] == NULL)
		    ktest = false;
		  j = j + 1;
		}
	      while (!((j > nprod) || (ktest == false)));
	      if (ktest == true)
		{
		  entry->mult = p11->mult * p22->mult;
		  entry->next = sublist;
		  sublist = entry;
		}
	      p22 = p22->next;
	    }
	  p22 = p2;
	  p11 = p11->next;
	}
      schur_psort (&sublist, true);
      R274 = sublist;
      sublist = NULL;
    }
  else
    {
      inform ("invalid group,no", cont);
      inform ("t implemented.", cr);
      R274 = NULL;
    }
  return R274;
}

prodtype
contract (prodtype pr, int n, int m, char *wrd)
{
  register prodtype R275;
  int k;
  register int j;
  prodtype top;
  ocharptr h1, h2;
  termptr l;
  top = pr;
  if (interp ("i_sfnproduct", wrd, 1))
    k = 1;
  else if (interp ("o_sfnproduct", wrd, 1))
    k = 2;
  else if (interp ("sk_sfn", wrd, 2))
    k = 3;
  else if (interp ("mixedtensorreps", wrd, 2))
    k = 5;
  else if (interp ("plethysmoutersfn", wrd, 2))
    k = 6;
  else if (interp ("rd_i_sfnproduct", wrd, 4))
    k = 7;
  else
    k = 4;
  while (top != NULL) {
      register prodtype W12 = &(*top);

      h1 = W12->prods.A[n - 1];
      h2 = W12->prods.A[m - 1];
      if (h1==NULL||h2==NULL) {
	      error(GROUP_NOT_SET,MAXSTRING);
	      return(pr);
      }
      
      switch ((int) (k)) {
	case 1:
	  l = inner (W12->prods.A[n - 1]->val, W12->prods.A[m - 1]->val);
	  break;
	case 2:
	  l = outer (W12->prods.A[n - 1]->val, W12->prods.A[m - 1]->val);
	  break;
	case 3:
	  l = skew (W12->prods.A[n - 1]->val, W12->prods.A[m - 1]->val);
	  break;
	case 4:
	  W12->prods.A[n - 1] =
	    kronk (W12->prods.A[n - 1], W12->prods.A[m - 1],
		   currgrp.A[n - 1]);
	  break;
	case 5:
	  W12->prods.A[n - 1] = oformbb (h1, h2, false, h1->spin);
	  break;
	case 6:
	  W12->prods.A[n - 1] =
	    gpleth (W12->prods.A[n - 1], W12->prods.A[m - 1],
		    currgrp.A[n - 1]);
	  break;
	case 7:
	  l = reduinnprd (W12->prods.A[n - 1]->val, W12->prods.A[m - 1]->val);
	  break;
	default:
	  Caseerror (Line);
	}
      switch ((int) (k)) {
	case 1:
	case 2:
	case 3:
	case 7:
	  W12->prods.A[n - 1] = sfntochrc (l, h1->spin, h1->lab);
	  ldisp (&l);
	  break;
	case 4:
	case 5:
	case 6:
	  break;
	default:
	  Caseerror (Line);
	}
      W12->prods.A[m - 1] = NULL;
      odisp (&h1);
      odisp (&h2);
      for (j = m; j <= nprod - 1; j++)
	{
	  W12->prods.A[j - 1] = W12->prods.A[j + 1 - 1];
	}
      W12->prods.A[nprod - 1] = NULL;
      top = top->next;
    }
  for (j = m; j <= nprod - 1; j++) {
      currgrp.A[j - 1] = currgrp.A[j + 1 - 1];
    }
  nprod = nprod - 1;
  putgroup (currgrp);
  top = pr;
  schur_psort (&top, true);
  R275 = top;
  return R275;
}


ocharptr
u1prod (ocharptr left, ocharptr right)
{
  ocharptr temp, rightemp, result;
  register int i;

  result = NULL;
  while (left != NULL)
    {
      rightemp = right;
      while (rightemp != NULL)
	{
	  if ((left->val.A[1] + rightemp->val.A[1] <= plwt))
	    {
	      cnu (&temp);
		temp->mult = 1;
		temp->C6_double = false;
		for (i = 1; i <= maxdim; i++)
		    temp->val.A[i] = 0;
		temp->val.A[1] = left->val.A[1] + rightemp->val.A[1];
		temp->lab = ' ';
		temp->spin = false;
		temp->next = result;
		result = temp;
	    }
	  rightemp = rightemp->next;
	}
      left = left->next;
    }
  osort (&result, false);
  return result;
}

ocharptr
unddprod (ocharptr firstx, ocharptr last, int n)
{
  register ocharptr R277;
  termptr int1, int2, int3, int4, alpha, beta,
    betop, fir, las, lhs, rhs, temp;
  ocharptr prodlist, sublist;
  frame core;
  bool spinor;

  core = nolls;
  prodlist = NULL;
  snu (&fir);
  snu (&las);
  fir->val = firstx->val;
  las->val = last->val;
  fir->mult = firstx->mult;
  las->mult = last->mult;
  fir->next = NULL;
  las->next = NULL;
  mergemin (firstx->conval, last->val, &core);
  alpha = skcompat (core);
  mergemin (firstx->val, last->conval, &core);
  betop = skcompat (core);
  if (firstx->spin == last->spin)
    spinor = false;
  else
    spinor = true;
  while (alpha != NULL)
    {
      int1 = skew (firstx->conval, alpha->val);
      int4 = pskew (las, alpha);
      beta = betop;
      while (beta != NULL)
	{
	  int2 = skew (last->conval, beta->val);
	  int3 = pskew (fir, beta);
	  lhs = louter (int1, int2);
	  rhs = louter (int3, int4);
	  sublist = formbb (lhs, rhs, true, spinor, ' ');
	  gumodify (&sublist, n);
	  ldisp (&lhs);
	  ldisp (&rhs);
	  oadd (&prodlist, &sublist);
	  ldisp (&int2);
	  ldisp (&int3);
	  beta = beta->next;
	}
      ldisp (&int1);
      ldisp (&int4);
      temp = alpha->next;
      dispsfn (&alpha);
      alpha = temp;
    }
  ldisp (&fir);
  ldisp (&las);
  osort (&prodlist, true);
  ldisp (&betop);
  R277 = prodlist;
  return R277;
}

/*ocharptr unitprod();*/

ocharptr
unssprod (ocharptr firstx, ocharptr last, int n)
{
  termptr prodx, temp;
  char ch;
  ocharptr newprod, lastptr;
  int multy;

  if (sslab)
    ch = last->lab;
  else
    ch = ' ';
  prodx = outer2 (firstx->val, last->val, n);
  multy = firstx->mult * last->mult;
  temp = prodx;
  lastptr = NULL;
  while (temp != NULL)
    {
      cnu (&newprod);
      {
	newprod->C6_double = false;
	newprod->lab = ch;
	if (firstx->spin == last->spin)
	  newprod->spin = false;
	else
	  newprod->spin = true;
	newprod->val = temp->val;
	newprod->mult = temp->mult * multy;
	newprod->next = lastptr;
	lastptr = newprod;
	temp = temp->next;
      }
    }
  ldisp (&prodx);
  return lastptr;
}

ocharptr
undsprod (ocharptr firstx, ocharptr last, int n)
{
  register ocharptr R280;
  ocharptr sublist, prodlist;
  termptr fir, las, xi, lhs, rhs1, rhs, temp;
  frame core;
  bool spinor;

  core = nolls;
  prodlist = NULL;
  snu (&fir);
  snu (&las);
  fir->val = firstx->val;
  las->val = last->val;
  fir->mult = firstx->mult;
  las->mult = last->mult;
  fir->next = NULL;
  las->next = NULL;
  mergemin (firstx->conval, last->val, &core);
  xi = skcompat (core);
  if (firstx->spin == last->spin)
    spinor = false;
  else
    spinor = true;
  while (xi != NULL)
    {
      lhs = skew (firstx->conval, xi->val);
      rhs1 = pskew (las, xi);
      rhs = louter2 (rhs1, fir, n + firstx->conval.A[1]);
      sublist = formbb (lhs, rhs, true, spinor, ' ');
      gumodify (&sublist, n);
      ldisp (&lhs);
      ldisp (&rhs);
      ldisp (&rhs1);
      oadd (&prodlist, &sublist);
      temp = xi;
      xi = xi->next;
      dispsfn (&temp);
    }
  osort (&prodlist, true);
  dispsfn (&fir);
  dispsfn (&las);
  R280 = prodlist;
  return R280;
}

ocharptr
unitprod (ocharptr left, ocharptr right, bool special, int n)
{
  ocharptr sublist, prodlist, rr;

  if (n == 1)
    return (u1prod (left, right));
  else
    {
      prodlist = NULL;
      while (left != NULL)
	{
	  rr = right;
	  while (rr != NULL)
	    {
	      if ((wprod && ((wtfrm (&rr->val) + wtfrm (&left->val)) <= plwt)
		   && !rr->C6_double && !left->C6_double) || (!wprod))
		{
		  if (left->C6_double)
		    if (rr->C6_double)
		      sublist = unddprod (left, rr, n);
		    else
		      sublist = undsprod (left, rr, n);
		  else if (rr->C6_double)
		    sublist = undsprod (rr, left, n);
		  else
		    sublist = unssprod (left, rr, n);
		  if (special)
		    spclise (&sublist, n);
		  oadd (&prodlist, &sublist);
		}
	      rr = rr->next;
	    }
	  osort (&prodlist, true);
	  left = left->next;
	}
      return prodlist;
    }
}

/*ocharptr sunmbrnch();*/

ocharptr
sunmbrnch1 (ocharptr irrep, int n, int m, int n2, int m2, int ss)
{
  register ocharptr R282;
  ocharptr co, cu, list, sublist;
  termptr zeta, tzeta, dum, ort, unnit;
  frame zero;


  list = NULL;
  zero = nolls;
  snu (&dum);
  {
    register termptr W24 = &(*dum);

    W24->val = zero;
    W24->next = NULL;
    W24->mult = 1;
  }
  while (irrep != NULL)
    {
      zeta = eqwt (wtfrm (&irrep->val));
      while (zeta != NULL)
	{
	  tzeta = zeta;
	  zeta = zeta->next;
	  tzeta->next = NULL;
	  ort = tzeta;
	  unnit = inner (irrep->val, tzeta->val);
	  cu = formbb (dum, unnit, true, false, ' ');
	  co = formbb (dum, ort, true, false, ' ');
	  if (ss == 1)
	    {
	      gumodify (&cu, n);
	      gumodify (&co, m);
	    }
	  else
	    {
	      unmmodify (&cu, n, m);
	      unmmodify (&co, n2, m2);
	    }
	  sublist = oformbb (cu, co, false, false);
	  odisp (&cu);
	  odisp (&co);
	  cu = chrcmult (irrep->mult, sublist);
	  odisp (&sublist);
	  oadd (&list, &cu);
	  ldisp (&unnit);
	  ldisp (&ort);
	}
      dispsfn (&dum);
      irrep = irrep->next;
    }
  osort (&list, false);
  R282 = list;
  return R282;
}

ocharptr
sunmbrnch (ocharptr irrep, int n, int m, int n2, int m2, int ss)
{
  register ocharptr R281;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (irrep != NULL)
    {
      register ocharptr W25 = &(*irrep);

      subbrnch = sunmbrnch1 (irrep, n, m, n2, m2, ss);
      oadd (&brnch, &subbrnch);
      irrep = W25->next;
    }
  osort (&brnch, true);
  R281 = brnch;
  return R281;
}

/*ocharptr orthprod();*/

ocharptr
onttprd (int n, ocharptr firstx, ocharptr last)
{
  register ocharptr R284;
  termptr tens1, tens2, zeta, subprod, subprod1, subprod2, prodx, temp;
  frame core;
  ocharptr lastptr, oprod;
  char ch;
  core = nolls;
  snu (&tens1);
  snu (&tens2);
  ch = ' ';
  if (((firstx->lab == ' ') && (last->lab == '#'))
      || ((firstx->lab == '#') && (last->lab == ' ')))
    ch = '#';
  tens1->val = firstx->val;
  tens2->val = last->val;
  tens1->mult = firstx->mult;
  tens2->mult = last->mult;
  tens1->slab = firstx->lab;
  tens2->slab = last->lab;
  tens1->next = NULL;
  tens2->next = NULL;
  mergemin (tens1->val, tens2->val, &core);
  prodx = NULL;
  zeta = skcompat (core);
  while (zeta != NULL)
    {
      subprod1 = pskew (tens1, zeta);
      subprod2 = pskew (tens2, zeta);
      subprod = louter (subprod1, subprod2);
      add (&prodx, &subprod);
      ldisp (&subprod1);
      ldisp (&subprod2);
      temp = zeta;
      zeta = zeta->next;
      dispsfn (&temp);
    }
  sort (&prodx, true);
  lastptr = NULL;
  temp = prodx;
  cslabel (temp, ch);
  while (temp != NULL)
    {
      cnu (&oprod);
      {
	register ocharptr W26 = &(*oprod);

	W26->mult = temp->mult;
	W26->val = temp->val;
	W26->lab = temp->slab;
	W26->C6_double = false;
	W26->spin = false;
	W26->next = lastptr;
	lastptr = oprod;
	temp = temp->next;
      }
    }
  ldisp (&tens1);
  ldisp (&tens2);
  ldisp (&prodx);
  omodify (&lastptr, n);
  R284 = lastptr;
  return R284;
}

ocharptr
onstprd (int n, ocharptr firstx, ocharptr last)
{
  register ocharptr R285;
  termptr div1, temp1, temp2, temp, zeta, spinor,
    tensor, prodx, subprod1, subprod2, subprod;
  frame one, two, three, core;
  register int i;
  ocharptr lastptr, oprod;
  char ch;
  one = nolls;
  two = nolls;
  three = nolls;
  core = nolls;
  sslab = false;
  snu (&spinor);
  snu (&tensor);
  spinor->mult = firstx->mult;
  tensor->mult = last->mult;
  spinor->val = firstx->val;
  spinor->slab = firstx->lab;
  tensor->val = last->val;
  tensor->slab = last->lab;
  spinor->next = NULL;
  tensor->next = NULL;
  div1 = seriesx ('q', -1, tensor->val);
  if ((bool) ((n) & 1))
    {
      sslab = true;
      slabel (div1);
      temp1 = lskew (tensor, div1);
      sslab = false;
    }
  else
    temp1 = lskew (tensor, div1);
  ldisp (&div1);
  for (i = 1; i <= maxdim; i++)
    {
      three.A[i] = 0;
    }
  temp2 = temp1;
  if (temp1 != NULL)
    {
      three = temp1->val;
      while (temp2->next != NULL)
	{
	  temp2 = temp2->next;
	  one = three;
	  two = temp2->val;
	  mergemax (one, two, &three);
	}
    }
  mergemin (spinor->val, three, &core);
  prodx = NULL;
  zeta = skcompat (core);
  while (zeta != NULL)
    {
      if (((bool) ((n) & 1)) && ((group == on)))
/*(Member((unsigned)(group), Conset[2]))))*/
	sslab = true;
      else
	sslab = false;
      zeta->slab = ' ';
      ch = ' ';
      subprod1 = pskew (spinor, zeta);
      temp = zeta;
      zeta = zeta->next;
      temp->next = NULL;
      cslabel (subprod1, ch);
      subprod2 = lskew (temp1, temp);
      subprod = louter (subprod1, subprod2);
      ldisp (&subprod1);
      ldisp (&subprod2);
      dispsfn (&temp);
      add (&prodx, &subprod);
      progress ();
    }
  ldisp (&temp1);
  sort (&prodx, true);
  lastptr = NULL;
  temp = prodx;
  while (temp != NULL)
    {
      cnu (&oprod);
      {
	register ocharptr W29 = &(*oprod);

	W29->C6_double = false;
	W29->mult = temp->mult;
	W29->val = temp->val;
	if (((bool) ((n) & 1)) && ((group == on)))
	  /*if (((bool)((n) & 1)) && (Member((unsigned)(group), Conset[3]))) */
	  W29->lab = temp->slab;
	else
	  W29->lab = ' ';
	W29->next = lastptr;
	W29->spin = true;
	lastptr = oprod;
	temp = temp->next;
      }
    }
  ldisp (&spinor);
  ldisp (&tensor);
  ldisp (&prodx);
  omodify (&lastptr, n);
  if (((bool) ((n) & 1)) && ((group == on)))
    /*if (((bool)((n) & 1)) && (Member((unsigned)(group), Conset[4]))) */
    olabel (lastptr);
  R285 = lastptr;
  sslab = false;
  return R285;
}

ocharptr
onssprd (int n, ocharptr firstx, ocharptr last)
{
  register ocharptr R286;
  termptr spin1, spin2, zeta, subprod, subprod1, subprod2,
    prodx, newprod, nusbprd, temp, qq;
  frame core, limitx;
  ocharptr tchar, lastptr, oprod, tempf, templ;
  register int i;
  register int m;
  char ch;
  core = nolls;
  limitx = nolls;
  cnu (&tempf);
  (*tempf) = (*firstx);
  tempf->next = NULL;
  cnu (&templ);
  (*templ) = (*last);
  templ->next = NULL;
  sslab = false;
  ch = ' ';
  /*if (((bool)((n) & 1) && (Member((unsigned)(group), Conset[5]))) || (!(bool)((n) & 1)) && (Member((unsigned)(group), Conset[6]))) { */
  if ((((bool) ((n) & 1) && ((group == on) || (group == son)))
       || ((!(bool) ((n) & 1)) && (group == on))))
    {
      i = n / 2;
      tempf->spin = false;
      templ->spin = false;
      for (m = 1; m <= i; m++)
	{
	  limitx.A[m] = 1;
	}
      m = 2 * i;
      fixg (&currgrp.A[1 - 1], spn, m, 0);
      oprod = sympprod (tempf, templ, m);
      tempf->spin = true;
      templ->spin = true;
      qq = seriesx ('q', i, limitx);
      tchar = sfntochrc (qq, false, ' ');
      ldisp (&qq);
      /*if (((bool)((n) & 1) && (Member((unsigned)(group), Conset[7])))) */
      if (((bool) ((n) & 1)) && ((group == on)))

	{
	  sslab = true;
	  oclabel (tchar);
	}
      fixg (&currgrp.A[1 - 1], un, i, 0);
      lastptr = unitprod (oprod, tchar, false, i);
      odisp (&oprod);
      odisp (&tchar);
      /*if ((Member((unsigned)(group), Conset[8]))) */
      if ((group == son))
	fixg (&currgrp.A[1 - 1], son, n, 0);
      else
	fixg (&currgrp.A[1 - 1], on, n, 0);
      if ((!(bool) ((n) & 1)))
	onexp (&lastptr, i);
    }
  else
    {
      /*if (((bool)((n) & 1) && (Member((unsigned)(group), Conset[9])))) */
      if (((bool) ((n) & 1)) && ((group == on)))
	if (((tempf->lab == '#') && (templ->lab == ' '))
	    || ((tempf->lab == ' ') && (templ->lab == '#')))
	  ch = '#';
      snu (&spin1);
      snu (&spin2);
      spin1->val = tempf->val;
      spin2->val = templ->val;
      spin1->mult = tempf->mult;
      spin2->mult = templ->mult;
      spin1->slab = ' ';
      spin2->slab = ' ';
      spin1->next = NULL;
      spin2->next = NULL;
      mergemin (spin1->val, spin2->val, &core);
      zeta = skcompat (core);
      prodx = NULL;
      while (zeta != NULL)
	{
	  subprod1 = pskew (spin1, zeta);
	  subprod2 = pskew (spin2, zeta);
	  subprod = louter (subprod1, subprod2);
	  add (&prodx, &subprod);
	  ldisp (&subprod1);
	  ldisp (&subprod2);
	  temp = zeta;
	  zeta = zeta->next;
	  dispsfn (&temp);
	}
      sort (&prodx, true);
      progress ();
      newprod = NULL;
      while (prodx != NULL)
	{
	  for (i = 0; i <= maxdim; i++)
	    {
	      limitx.A[i] = 0;
	    }
	  m = n + prodx->val.A[1];
	  prodx->slab = ch;
	  for (i = 1; i <= m; i++)
	    {
	      limitx.A[i] = 1;
	    }
	  qq = seriesx ('q', -1, limitx);
	  /*if ((bool)((n) & 1) && (Member((unsigned)(group), Conset[10]))) */
	  if (((bool) ((n) & 1)) && ((group == on)))
	    sslab = true;
	  if (sslab == true)
	    slabel (qq);
	  temp = prodx;
	  prodx = prodx->next;
	  temp->next = NULL;
	  nusbprd = louter2 (qq, temp, m + 1);
	  ldisp (&qq);
	  dispsfn (&temp);
	  add (&newprod, &nusbprd);
	}
      lastptr = NULL;
      while (newprod != NULL)
	{
	  cnu (&oprod);
	  {
	    register ocharptr W36 = &(*oprod);

	    W36->mult = newprod->mult;
	    W36->val = newprod->val;
	    progress ();
	    if (((bool) ((n) & 1)) && ((group == on)))
	      /*if ((bool)((n) & 1) && (Member((unsigned)(group), Conset[11]))) */
	      W36->lab = newprod->slab;
	    else
	      W36->lab = ' ';
	    W36->C6_double = false;
	    W36->spin = false;
	    W36->next = lastptr;
	    lastptr = oprod;
	    temp = newprod;
	    newprod = newprod->next;
	    dispsfn (&temp);
	  }
	}
      omodify (&lastptr, n);
      oprod = lastptr;
      if ((bool) ((n) & 1))
	while (oprod != NULL)
	  {
	    register ocharptr W37 = &(*oprod);

	    W37->mult = W37->mult / 2;
	    oprod = W37->next;
	  }
      ldisp (&spin1);
      ldisp (&spin2);
    }
  R286 = lastptr;
  sslab = false;
  dispchr (&templ);
  dispchr (&tempf);
  return R286;
}

ocharptr
so2nprd1 (ocharptr lambda, ocharptr mu, int n)
{
  register ocharptr R287;
  ocharptr list, stornext, sub1, sublist;
  termptr slam, xi, xit, temp1, temp2, temp, bb;
  stornext = mu->next;
  mu->next = NULL;
  list = NULL;
  snu (&slam);
  {
    register termptr W38 = &(*slam);

    W38->val = lambda->val;
    W38->mult = lambda->mult;
    W38->next = NULL;
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
	  bb = seriesx ('b', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  ldisp (&bb);
	  dispsfn (&temp);
	  sub1 = formbb (xit, temp2, true, false, ' ');
	  ldisp (&temp2);
	  gumodify (&sub1, n / 2);
	  sublist = unitprod (sub1, mu, false, n / 2);
	  odisp (&sub1);
	  oadd (&list, &sublist);
	}
      dispsfn (&xit);
    }
  dispsfn (&slam);
  mu->next = stornext;
  so2nmod (&list, mu->lab, n);
  R287 = list;
  return R287;
}

ocharptr
so2nprd2 (ocharptr lambda, ocharptr mu, int n)
{
  register ocharptr R288;
  ocharptr list, stornext, sub1, sublist;
  termptr slam, xi, xit, temp1, temp2, temp, bb, temp3, ser;
  int limitx;
  register int i;
  char labyl;
  frame dum, dumm;
  list = NULL;
  dum = nolls;
  dumm = nolls;
  stornext = mu->next;
  mu->next = NULL;
  snu (&slam);
  {
    register termptr W39 = &(*slam);

    if (mu->spin)
      W39->val = lambda->val;
    else
      {
	for (i = 0; i <= maxdim - 1; i++)
	  {
	    W39->val.A[i] = MAX (0, lambda->val.A[i] - 1);
	    dum.A[i] = MAX (0, mu->val.A[i] - 1);
	  }
	dumm = mu->val;
	mu->val = dum;
      }
    W39->mult = lambda->mult;
    W39->next = NULL;
  }
  if ((bool) ((n / 2) & 1))
    if (mu->lab == lambda->lab)
      labyl = '-';
    else
      labyl = '+';
  else if (mu->lab == lambda->lab)
    labyl = '+';
  else
    labyl = '-';
  xi = skcompat (slam->val);
  while (xi != NULL)
    {
      temp1 = pskew (slam, xi);
      xit = xi;
      xi = xi->next;
      xit->next = NULL;
      limitx = xit->val.A[1] + n / 2;
      while (temp1 != NULL)
	{
	  bb = seriesx ('b', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  if (mu->spin)
	    ser = ql (labyl, limitx);
	  else
	    ser = xv (labyl, limitx);
	  temp3 = louter2 (temp2, ser, limitx);
	  ldisp (&bb);
	  ldisp (&ser);
	  dispsfn (&temp);
	  sub1 = formbb (xit, temp3, true, lambda->spin, ' ');
	  gumodify (&sub1, n / 2);
	  sublist = unitprod (sub1, mu, false, n / 2);
	  odisp (&sub1);
	  oadd (&list, &sublist);
	  ldisp (&temp2);
	  ldisp (&temp3);
	}
      dispsfn (&xit);
    }
  if (!mu->spin)
    mu->val = dumm;
  dispsfn (&slam);
  mu->next = stornext;
  so2nmod (&list, mu->lab, n);
  R288 = list;
  return R288;
}

ocharptr
so2nprd3 (ocharptr lambda, ocharptr mu, int n)
{
  register ocharptr R289;
  ocharptr list, stornext, sub1, sublist;
  termptr slam, xi, temp1, temp2, temp, bb, temp3, ser, xit;
  int limitx;
  char labyl;
  stornext = mu->next;
  mu->next = NULL;
  list = NULL;
  snu (&slam);
  {
    register termptr W42 = &(*slam);

    W42->val = lambda->val;
    W42->mult = lambda->mult;
    W42->next = NULL;
  }
  xi = skcompat (slam->val);
  limitx = slam->val.A[1] + n / 2;
  if (mu->lab == lambda->lab)
    labyl = '+';
  else
    labyl = '-';
  while (xi != NULL)
    {
      xit = xi;
      xi = xi->next;
      xit->next = NULL;
      temp1 = pskew (slam, xit);
      while (temp1 != NULL)
	{
	  bb = seriesx ('b', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  dispsfn (&temp);
	  ser = ql (labyl, limitx);
	  temp3 = louter2 (xit, ser, limitx);
	  ldisp (&bb);
	  ldisp (&ser);
	  sub1 = formbb (temp3, temp2, true, true, ' ');
	  gumodify (&sub1, n / 2);
	  sublist = unitprod (sub1, mu, false, n / 2);
	  odisp (&sub1);
	  oadd (&list, &sublist);
	  ldisp (&temp2);
	  ldisp (&temp3);
	}
      dispsfn (&xit);
    }
  dispsfn (&slam);
  mu->next = stornext;
  so2nmod (&list, mu->lab, n);
  R289 = list;
  return R289;
}

ocharptr
orthprod (ocharptr left, ocharptr right, bool special, int n)
{
  register ocharptr R283;
  ocharptr prodlist, sublist, rr;
  int lwt, rwtx;
  prodlist = NULL;
  if (special && !(bool) ((n) & 1))
    while (left != NULL)
      {
	rr = right;
	lwt = wtfrm (&left->val);
	while (rr != NULL)
	  {
	    rwtx = wtfrm (&rr->val);
	    if (left->spin)
	      if (rr->spin)
		if (lwt < rwtx)
		  sublist = so2nprd2 (left, rr, n);
		else
		  sublist = so2nprd2 (rr, left, n);
	      else if (rr->lab == ' ')
		sublist = so2nprd1 (rr, left, n);
	      else
		sublist = so2nprd3 (left, rr, n);
	    else if (rr->spin)
	      if (left->lab == ' ')
		sublist = so2nprd1 (left, rr, n);
	      else
		sublist = so2nprd3 (rr, left, n);
	    else if (left->lab == ' ')
	      if (rr->lab == ' ')
		{
		  sublist = onttprd (n, left, rr);
		  so2nexp (&sublist, n);
		}
	      else
		sublist = so2nprd1 (left, rr, n);
	    else if (rr->lab == ' ')
	      sublist = so2nprd1 (rr, left, n);
	    else if (lwt < rwtx)
	      sublist = so2nprd2 (left, rr, n);
	    else
	      sublist = so2nprd2 (rr, left, n);
	    oadd (&prodlist, &sublist);
	    rr = rr->next;
	  }
	left = left->next;
      }
  else
    while (left != NULL)
      {
	rr = right;
	while (rr != NULL)
	  {
	    if (left->spin)
	      if (rr->spin)
		sublist = onssprd (n, left, rr);
	      else
		sublist = onstprd (n, left, rr);
	    else if (rr->spin)
	      sublist = onstprd (n, rr, left);
	    else
	      sublist = onttprd (n, left, rr);
	    oadd (&prodlist, &sublist);
	    rr = rr->next;
	  }
	left = left->next;
      }
  osort (&prodlist, true);
  R283 = prodlist;
  return R283;
}

void
onttscalar (ocharptr list1, ocharptr list2, groop grp)
{
  ocharptr temp2;
  int i, total;
  bool test;

  total = 0;
  if (grp.name == on)
    {
      while (list1 != NULL)
	{
	  temp2 = list2;
	  while (temp2 != NULL)
	    {
	      test = true;
	      i = 1;
	      if ((list1->lab != temp2->lab))
		test = false;
	      else
		do
		  {
		    if ((list1->val.A[i] != temp2->val.A[i]))
		      test = false;
		    i = i + 1;
		  }
		while (!(((i > maxdim - 1) || (test == false))));
	      if (test == true)
		total = total + list1->mult * temp2->mult;
	      temp2 = temp2->next;
	    }
	  list1 = list1->next;
	}
      inform ("# of scalars =;", cont);
      print ( "%10d\n", total);
    }
  else
    print ( "Group not set as an orthogonal group O(n)\n");
}

ocharptr
lspncprod (ocharptr list1, ocharptr list2, int n)
{
  register ocharptr R290;
  ocharptr newlist, temp, ll2;

  newlist = NULL;
  while (list1 != NULL)
    {
      register ocharptr W43 = &(*list1);

      ll2 = list2;
      while (ll2 != NULL)
	{
	  register ocharptr W44 = &(*ll2);

	  temp = spncprod (list1, ll2, n);
	  oadd (&newlist, &temp);
	  ll2 = W44->next;
	}
      list1 = W43->next;
    }
  spncmodify (&newlist, n);
  sposort (&newlist, true);
  R290 = newlist;
  return R290;
}

ocharptr
spncprod (ocharptr left, ocharptr right, int n)
{
  register ocharptr R291;
  ocharptr temp, newlist, tempc;
  termptr munu, mu, nu, ds, templist, nud, dummy;
  int k, l, nn, kl, kk, ll, nnn, xk, restx;
  bool spinkl, ppp;
  char tag;

  newlist = NULL;
  if ((wtfrm (&left->val) + wtfrm (&right->val) <= (int) (setlimit)))
    {
      nnn = n / 2;
      restx = setlimit - 1;
      k = left->val.A[maxdim];
      l = right->val.A[maxdim];
      if (left->spin)
	kk = 2 * k + 1;
      else
	kk = 2 * k;
      if (right->spin)
	ll = 2 * l + 1;
      else
	ll = 2 * l;
      if ((ll > kk) || ((ll == kk) && (left->val.A[1] < right->val.A[1])))
	{
	  tempc = left;
	  left = right;
	  right = tempc;
	  xk = k;
	  l = k;
	  k = xk;
	  xk = ll;
	  ll = kk;
	  kk = xk;
	}
      kl = (kk + ll) / 2;
      if ((bool) ((kk + ll) & 1))
	spinkl = true;
      else
	spinkl = false;
      nn = MIN (nnn, ll);
      snu (&mu);
      mu->mult = left->mult;
      mu->val = left->val;
      mu->val.A[maxdim] = 0;
      ds = rseries (nn, 'd');
      if ((bool) ((ll) & 1))
	tag = 'g';
      else
	tag = 'c';
      snu (&templist);
      templist->mult = right->mult;
      templist->val = right->val;
      templist->val.A[maxdim] = 0;
      dummy = templist;
      if (ll == nn)
	ppp = true;
      else
	ppp = false;
      nu = signseq (&dummy, ll, nn, tag, ppp, false);
      dispsfn (&dummy);
      nud = louter2 (ds, nu, nn);
      ldisp (&ds);
      ldisp (&nu);
      munu = louter2 (mu, nud, nnn);
      schur_restrict (&munu, restx, 'w');
      dispsfn (&mu);
      ldisp (&nud);
      nud = munu;
      while (nud != NULL)
	{
	  register termptr W45 = &(*nud);

	  cnu (&temp);
	  temp->mult = W45->mult;
	  temp->val = W45->val;
	  temp->C6_double = true;
	  temp->conval = nolls;
	  temp->conval.A[1] = kl;
	  temp->val.A[maxdim] = kl;
	  temp->spin = spinkl;
	  temp->lab = ' ';
	  temp->conlab = ' ';
	  oadd (&newlist, &temp);
	  nud = W45->next;
	}
      ldisp (&munu);
      sposort (&newlist, true);
    }
  R291 = newlist;
  return R291;
}

/*ocharptr sympprod();*/

ocharptr
spnprd (int n, ocharptr firstx, ocharptr last)
{
  register ocharptr R293;
  termptr tens1, tens2, zeta, subprod, subprod1, subprod2, prodx, temp;
  frame core;
  ocharptr lastptr, spprod;

  core = nolls;
  snu (&tens1);
  snu (&tens2);
  tens1->val = firstx->val;
  tens2->val = last->val;
  tens1->mult = firstx->mult;
  tens2->mult = last->mult;
  tens1->next = NULL;
  tens2->next = NULL;
  mergemin (tens1->val, tens2->val, &core);
  zeta = skcompat (core);
  prodx = NULL;
  while (zeta != NULL)
    {
      subprod1 = pskew (tens1, zeta);
      subprod2 = pskew (tens2, zeta);
      subprod = louter (subprod1, subprod2);
      add (&prodx, &subprod);
      ldisp (&subprod1);
      ldisp (&subprod2);
      temp = zeta;
      zeta = zeta->next;
      dispsfn (&temp);
    }
  sort (&prodx, true);
  lastptr = NULL;
  temp = prodx;
  while (temp != NULL)
    {
      cnu (&spprod);
      {
	register ocharptr W46 = &(*spprod);

	W46->mult = temp->mult;
	W46->val = temp->val;
	W46->lab = ' ';
	W46->C6_double = false;
	W46->spin = false;
	W46->next = lastptr;
	lastptr = spprod;
	temp = temp->next;
      }
    }
  ldisp (&prodx);
  ldisp (&tens1);
  ldisp (&tens2);
  spmodify (&lastptr, n);
  R293 = lastptr;
  return R293;
}

ocharptr
sympprod (ocharptr left, ocharptr right, int n)
{
  register ocharptr R292;
  ocharptr sublist, prodlist, rr;

  prodlist = NULL;
  while (left != NULL)
    {
      rr = right;
      while (rr != NULL)
	{
	  sublist = spnprd (n, left, rr);
	  oadd (&prodlist, &sublist);
	  rr = rr->next;
	}
      left = left->next;
    }
  osort (&prodlist, true);
  R292 = prodlist;
  return R292;
}

ocharptr
su2su6prod (ocharptr left, ocharptr right)
{
  register ocharptr R294;
  ocharptr product, rr, last, subpr;
  termptr su6pr, dsu6pr;
  frame ldum, rdum;
  int s1, s2, j;
  register int i;

  ldum = nolls;
  rdum = nolls;
  product = NULL;
  while (left != NULL)
    {
      rr = right;
      while (rr != NULL)
	{
	  s1 = left->val.A[1];
	  s2 = rr->val.A[1];
	  for (i = 1; i <= maxdim - 1; i++)
	    {
	      ldum.A[i] = left->val.A[i + 1];
	      rdum.A[i] = rr->val.A[i + 1];
	    }
	  ldum.A[0] = 0;
	  ldum.A[maxdim] = 0;
	  rdum.A[0] = 0;
	  rdum.A[maxdim] = 0;
	  su6pr = outer2 (ldum, rdum, 6);
	  dsu6pr = su6pr;
	  while (dsu6pr != NULL)
	    {
	      register termptr W49 = &(*dsu6pr);

	      if (W49->val.A[6] != 0)
		for (i = 1; i <= 6; i++)
		  {
		    W49->val.A[i] = W49->val.A[i] - W49->val.A[6];
		  }
	      dsu6pr = W49->next;
	    }
	  sort (&su6pr, true);
	  last = NULL;
	  j = abs (s1 - s2);
	  do
	    {
	      dsu6pr = su6pr;
	      while (dsu6pr != NULL)
		{
		  cnu (&subpr);
		  {
		    register ocharptr W52 = &(*subpr);

		    W52->val.A[0] = 0;
		    W52->val.A[1] = j;
		    for (i = 2; i <= maxdim - 1; i++)
		      {
			W52->val.A[i] = dsu6pr->val.A[i - 1];
		      }
		    W52->mult = dsu6pr->mult * left->mult * rr->mult;
		    W52->spin = false;
		    W52->lab = ' ';
		    W52->C6_double = false;
		    W52->next = last;
		    last = subpr;
		  }
		  dsu6pr = dsu6pr->next;
		}
	      j = j + 2;
	    }
	  while (!(j > s1 + s2));
	  oadd (&product, &subpr);
	  ldisp (&su6pr);
	  rr = rr->next;
	}
      left = left->next;
    }
  osort (&product, false);
  R294 = product;
  return R294;
}

ocharptr
ssymmprod (ocharptr list1, ocharptr list2, int n)
{
  register ocharptr R295;
  termptr sublist, x, newlist, stemp;
  ocharptr templist, temp, templ, tlist, sulist, temp1,
    temp2, prodlist, templist1;
  bool test1, test2, test;
  frame nu;
  int multy, s, t, m, l, p;
  register int z;

  prodlist = NULL;
  sublist = NULL;
  p = len (&list1->val);
  cnu (&temp1);
  (*temp1) = (*list1);
  temp1->next = NULL;
  if (temp1->spin)
    test1 = true;
  else
    test1 = false;
  templist = list2;
  if (temp1->lab == '-')
    s = -1;
  else if (temp1->lab == '+')
    s = 1;
  else
    s = 0;
  cnu (&temp2);
  (*temp2) = (*templist);
  temp2->next = NULL;
  if (temp2->spin)
    test2 = true;
  else
    test2 = false;
  if (temp2->lab == '-')
    t = -1;
  else if (temp2->lab == '+')
    t = 1;
  else
    t = 0;
  if (((!test1) && test2))
    {
      templ = temp1;
      temp1 = temp2;
      temp2 = templ;
      test = test1;
      test1 = test2;
      test2 = test;
    }
  if ((test1) && (test2))
    sublist = lsinner (&temp1->val, &temp2->val, s, t, n);
  else if ((test1) && (!test2))
    sublist = lsqinner (temp1->val, temp2->val, n);
  else if ((!test1) && (!test2))
    sublist = inner (temp1->val, temp2->val);
  if (((sublist != NULL)
       && ((((s != 0) && (t != 0)) || ((s == 0) && (t == 0)))
	   || ((((s != 0) && (t == 0)) || ((s == 0) && (t != 0)))
	       && (test1 == test2)))))
    {
      multy = temp1->mult * temp2->mult;
      x = sublist;
      if (multy != 1)
	while (x != NULL)
	  {
	    x->mult = x->mult * multy;
	    x = x->next;
	  }
    }
  if ((test1 && (!test2)))
    {
      newlist = reducedq (sublist);
      ldisp (&sublist);
      tlist = NULL;
      sublist = newlist;
      while (sublist != NULL)
	{
	  register termptr W55 = &(*sublist);
	  if (((n % 2 == 0) && (p % 2 == 1))
	      || ((n % 2 == 1) && (p % 2 == 0)))
	    {
	      cnu (&temp);
	      temp->mult = W55->mult;
	      temp->val = W55->val;
	      temp->spin = true;
	      temp->lab = ' ';
	      temp->C6_double = false;
	      oadd (&tlist, &temp);
	      dispchr (&temp);
	    }
	  else
	    {
	      cnu (&temp);
	      temp->mult = W55->mult / 2;
	      temp->val = W55->val;
	      temp->lab = ' ';
	      temp->spin = true;
	      temp->C6_double = false;
	       /**/ oadd (&tlist, &temp);
	      dispchr (&temp);
	      if ((W55->mult % 2 == 1))
		{
		  m = n - wtfrm (&temp1->val);
		  l = len (&temp1->val);
		  nu = nolls;
		  nu.A[1] = m;
		  for (z = 2; z <= l + 1; z++)
		    {
		      nu.A[z] = temp1->val.A[z - 1];
		    }
		  cnu (&temp);
		  temp->mult = 1;
		  temp->val = W55->val;
		  temp->spin = true;
		  temp->lab = ' ';
		  temp->C6_double = false;
		  snu (&stemp);
		  stemp->mult = 1;
		  stemp->val = temp2->val;
		  s = s * dcoeffq (&stemp, &nu);
		  dispsfn (&stemp);
		  if (s == 0)
		    temp->lab = ' ';
		  else if (s == 1)
		    temp->lab = '+';
		  else if (s == -1)
		    temp->lab = '-';
		  oadd (&tlist, &temp);
		  dispchr (&temp);
		}
	    }
	  sublist = W55->next;
	}
      ldisp (&newlist);
      multy = temp1->mult * temp2->mult;
      sulist = tlist;
      if (multy != 1)
	while (sulist != NULL)
	  {
	    register ocharptr W58 = &(*sulist);

	    sulist->mult = W58->mult * multy;
	    sulist = W58->next;
	  }
      templist1 = tlist;
    }
  else
    {
      templist1 = sfntochrc (sublist, false, ' ');
      ldisp (&sublist);
    }
  oadd (&prodlist, &templist1);
  odisp (&templist1);
  odisp (&temp1);
  odisp (&temp2);
  osort (&prodlist, false);
  snselect (&prodlist, n, qspecial);
  R295 = prodlist;
  return R295;
}



ocharptr
symmprod (ocharptr list1, ocharptr list2, int n)
{
  register ocharptr R296;
  ocharptr newlist, temp, rr;

  newlist = NULL;
  snselect (&list1, n, qspecial);
  snselect (&list2, n, qspecial);
   /**/ while (list1 != NULL)
    {
      rr = list2;
      while (rr != NULL)
	{
	  temp = ssymmprod (list1, rr, n);
	  oadd (&newlist, &temp);
	  rr = rr->next;
	}
      list1 = list1->next;
    }
  osort (&newlist, false);
  R296 = newlist;
  return R296;
}


ocharptr
ebrnch (ocharptr ssuper, ocharptr index, tabptr brtab)
{
  ocharptr subbr, tindex, br;
  bool located;

  br = NULL;
  while (ssuper != NULL)
    {
      register ocharptr W59 = &(*ssuper);

      located = false;
      tindex = index;
      while (!located && (tindex != NULL))
	if (otestord (ssuper, tindex) == EQUAL)
	  located = true;
	else
	  tindex = tindex->next;
      if (located)
	{
	  subbr = extract (tindex->mult, W59->mult, brtab);
	  oadd (&br, &subbr);
	  ssuper = W59->next;
	}
      else
	{
	  print ( "branching not found\n");
	  ssuper = NULL;
	}
    }
  return br;
}

ocharptr
g2brnch (ocharptr list)
{
  register ocharptr R298;
  ocharptr bran, brnch, indiv;
  register int i;
  register int t;
  register int s;
  register int r, lastStep, lastStep2;

  bran = NULL;
  while (list != NULL)
    {
      register ocharptr W60 = &(*list);

      brnch = NULL;
      lastStep = W60->val.A[2];
      for (r = 0; r <= lastStep; r++)
	{
	  lastStep2 = W60->val.A[1] - 2 * W60->val.A[2];
	  for (s = 0; s <= lastStep2; s++)
	    {
	      for (t = 0; t <= r + s; t++)
		{
		  cnu (&indiv);
		  for (i = 0; i <= maxdim; i++)
		    {
		      indiv->val.A[i] = 0;
		    }
		  indiv->val.A[1] = W60->val.A[1] - W60->val.A[2] + r - t;
		  indiv->val.A[2] = W60->val.A[2] + s - t;
		  indiv->mult = W60->mult;
		  indiv->spin = false;
		  indiv->C6_double = false;
		  indiv->lab = ' ';
		  indiv->next = brnch;
		  brnch = indiv;
		}
	    }
	}
      list = W60->next;
      oadd (&bran, &brnch);
    }
  osort (&bran, false);
  R298 = bran;
  return R298;
}

ocharptr
g2prod (ocharptr left, ocharptr right)
{
  register ocharptr R299;
  ocharptr brnch, dum, prodx;
  /*int       eps, i; *//*12/12/95 */

  brnch = g2brnch (left);
  dum = right;
  while (dum != NULL)
    {
      dum->val.A[1] = dum->val.A[1] + 1;
      dum = dum->next;
    }
  prodx = unitprod (brnch, right, true, 3);
  dum = right;
  odisp (&brnch);
  while (dum != NULL)
    {
      dum->val.A[1] = dum->val.A[1] - 1;
      dum = dum->next;
    }
  g2modify (&prodx, 1);
  R299 = prodx;
  return R299;
}

ocharptr
e8prod (ocharptr left, ocharptr right)
{
  register ocharptr R300;
  ocharptr temp, list, lastptr, newlist, branch;

  newlist = NULL;
  if (!e8load)
    brload (49);
  while (left != NULL)
    {
      while (right != NULL)
	{
	  if ((left->val.A[1] > right->val.A[1])
	      || ((left->val.A[1] == right->val.A[1])
		  && (len (&left->val) < len (&right->val))))
	    {
	      temp = right;
	      right = left;
	      left = temp;
	    }
	  branch = ebrnch (left, e8index, e8tab);
	  temp = branch;
	  while (temp != NULL)
	    {
	      lastptr = right;
	      lastptr->val.A[1] = lastptr->val.A[1] + 21;
	      list = unddprod (right, temp, 9);
	      spclise (&list, 9);
	      lastptr = right;
	      lastptr->val.A[1] = lastptr->val.A[1] - 21;
	      e8modify (&list, 21);
	      oadd (&newlist, &list);
	      osort (&newlist, true);
	      temp = temp->next;
	    }
	  right = right->next;
	}
      left = left->next;
    }
  odisp (&branch);
  R300 = newlist;
  return R300;
}

ocharptr
e7prod (ocharptr left, ocharptr right)
{
  register ocharptr R301;
  ocharptr temp, list, lastptr, newlist, branch;

  newlist = NULL;
  if (!e7load)
    brload (47);
  while (left != NULL)
    {
      while (right != NULL)
	{
	  if ((left->val.A[1] > right->val.A[1])
	      || ((left->val.A[1] == right->val.A[1])
		  && (len (&left->val) < len (&right->val))))
	    {
	      temp = right;
	      right = left;
	      left = temp;
	    }
	  branch = ebrnch (left, e7index, e7tab);
	  lastptr = right->next;
	  right->next = NULL;
	  right->val.A[1] = right->val.A[1] + 10;
	  list = unitprod (right, branch, true, 8);
	  right->val.A[1] = right->val.A[1] - 10;
	  right->next = lastptr;
	  odisp (&branch);
	  e7modify (&list, 10);
	  oadd (&newlist, &list);
	  right = right->next;
	}
      left = left->next;
      osort (&newlist, true);
    }
  R301 = newlist;
  return R301;
}

ocharptr
e6prod (ocharptr left, ocharptr right)
{
  register ocharptr R302;
  ocharptr newlist, list, lastptr, temp, branch;

  newlist = NULL;
  if (!e6load)
    brload (44);
  while (left != NULL)
    {
      while (right != NULL)
	{
	  if ((left->val.A[1] > right->val.A[1])
	      || ((left->val.A[1] == right->val.A[1])
		  && (len (&left->val) < len (&right->val))))
	    {
	      temp = right;
	      right = left;
	      left = temp;
	    }
	  branch = ebrnch (left, e6index, e6tab);
	  lastptr = right->next;
	  right->next = NULL;
	  right->val.A[1] = right->val.A[1] + 10;
	  list = su2su6prod (branch, right);
	  right->val.A[1] = right->val.A[1] - 10;
	  right->next = lastptr;
	  odisp (&branch);
	  e6modify (&list, 10);
	  oadd (&newlist, &list);
	  right = right->next;
	}
      left = left->next;
      osort (&newlist, true);
    }
  R302 = newlist;
  return R302;
}

ocharptr
f4prod (ocharptr left, ocharptr right)
{
  register ocharptr R303;
  ocharptr temp, list, lastptr, newlist, branch;

  newlist = NULL;
  if (!f4load)
    brload (43);
  while (left != NULL)
    {
      while (right != NULL)
	{
	  if ((left->val.A[1] > right->val.A[1])
	      || ((left->val.A[1] == right->val.A[1])
		  && (len (&left->val) < len (&right->val))))
	    {
	      temp = right;
	      right = left;
	      left = temp;
	    }
	  branch = ebrnch (left, f4index, f4tab);
	  lastptr = right->next;
	  right->next = NULL;
	  right->val.A[1] = right->val.A[1] + 2;
	  list = orthprod (branch, right, true, 9);
	  right->val.A[1] = right->val.A[1] - 2;
	  right->next = lastptr;
	  odisp (&branch);
	  f4modify (&list, 2);
	  oadd (&newlist, &list);
	  right = right->next;
	}
      left = left->next;
      osort (&newlist, true);
    }
  R303 = newlist;
  return R303;
}

ocharptr
ltspncprod (ocharptr list1, ocharptr list2, int n)
{
  register ocharptr R304;
  ocharptr temp, newlist, tlist2;

  newlist = NULL;
  while (list1 != NULL)
    {
      tlist2 = list2;
      while (tlist2 != NULL)
	{
	  temp = tspncprod (list1, tlist2, n);
	  oadd (&newlist, &temp);
	  tlist2 = tlist2->next;
	}
      list1 = list1->next;
    }
  osort (&newlist, true);
  R304 = newlist;
  return R304;
}


ocharptr
lsoncprod (ocharptr list1, ocharptr list2, int n)
{
  register ocharptr R305;
  ocharptr newlist, temp, ll2;

  newlist = NULL;
  while (list1 != NULL)
    {
      register ocharptr W43 = &(*list1);

      ll2 = list2;
      while (ll2 != NULL)
	{
	  register ocharptr W44 = &(*ll2);

	  temp = soncprod (list1, ll2, n);
	  oadd (&newlist, &temp);
	  ll2 = W44->next;
	}
      list1 = W43->next;
    }
  /*spncmodify(&newlist, n); */
  osort (&newlist, true);
  R305 = newlist;
  return R305;
}


ocharptr
soncprod (ocharptr left, ocharptr right, int n)
{
  ocharptr temp, tempc, newlist;
  termptr munu, mu, nu, nuu, nudd, ds, templ, templist, nud, dummy, dummyy;
  int pr, pl, k, l, /*nn,*/ kl, kk, ll, nnn, xk, restx;
  bool ppp;
  //char tag;

  newlist = NULL;
  nnn = n / 2;
  k = left->val.A[maxdim];
  l = right->val.A[maxdim];
  kk = 2 * k;
  ll = 2 * l;
  if ((ll > kk) || ((ll == kk) && (left->val.A[1] < right->val.A[1])))
    {
      tempc = left;
      left = right;
      right = tempc;
      xk = k;
      l = k;
      k = xk;
      xk = ll;
      ll = kk;
      kk = xk;
    }
  kl = (kk + ll) / 2;
  pl = len (&left->val);
  pr = len (&right->val);
  //nn = MIN (nnn, ll);
  nnn = MIN (nnn, kl);
  snu (&mu);
  mu->mult = left->mult;
  mu->val = left->val;
  mu->val.A[maxdim] = 0;
  ds = rseries (nnn, 'b');
  restx = setlimit - pr - pl;
  schur_restrict (&ds, restx, 'w');
  /*if (ll & 1)
    tag = 'a'; */
  snu (&templist);
  templist->mult = right->mult;
  templist->val = right->val;
  templist->val.A[maxdim] = 0;
  dummy = templist;
  if (ll == nnn)
    ppp = true;
  else
    ppp = false;
  nu = signseq (&dummy, ll, nnn, 'a', ppp, false);
  restx = setlimit - pr;
  schur_restrict (&nu, restx, 'w');
  snu (&templ);
  templ->mult = left->mult;
  templ->val = left->val;
  templ->val.A[maxdim] = 0;
  dummyy = templ;
  if (kk == nnn)
    ppp = true;
  else
    ppp = false;
  nuu = signseq (&dummyy, kk, nnn, 'a', ppp, false);
  restx = setlimit - pl;
  schur_restrict (&nuu, restx, 'w');
  nud = louter2 (nu, nuu, nnn);
  ldisp (&nu);
  ldisp (&nuu);
  dispsfn (&dummy);
  dispsfn (&dummyy);
  munu = louter2 (ds, nud, nnn);
  schur_restrict (&munu, setlimit, 'w');
  ldisp (&ds);
  ldisp (&nud);
  nudd = munu;
  while (nudd != NULL)
    {
      register termptr W4 = &(*nudd);

      cnu (&temp);
      temp->mult = W4->mult;
      temp->val = W4->val;
      temp->C6_double = true;
      temp->conval = nolls;
      temp->conval.A[1] = kl;
      temp->val.A[maxdim] = kl;
      temp->spin = false;
      temp->lab = ' ';
      temp->conlab = ' ';
      oadd (&newlist, &temp);
      nudd = W4->next;
    }
  ldisp (&munu);
  osort (&newlist, true);
  return newlist;
}

bool
generic (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  bool test;
  int n, k, s;

  test = false;
  if (((grp.name == spnc) || (grp.name == sonc)))
    {
      n = grp.rank / 2;
      k = chrc2->val.A[maxdim];
      /*if (chrc2->spin)
         k = 2 * k + 1;
         else
         k = 2 * k; */

      s = n + chrc1->val.A[1];

      if ((k >= s))
	test = true;
    }
  else
    INAPPROPRIATE_GROUP();
  return test;
}

ocharptr
covariant (ocharptr chrc, groop grp)
{
  ocharptr newlist, covar;
  int n, /*p,*/ q;
  register int i;
  newlist = NULL;
  if (grp.name == un)
    {
      n = grp.rank;
      while (chrc != NULL)
	{
	  register ocharptr W67 = &(*chrc);
	  cnu (&covar);
	  covar->mult = W67->mult;
	  covar->C6_double = true;
	  covar->val = nolls;
	  covar->conval = nolls;
	  covar->spin = false;
	  if (chrc->C6_double)
	    {
	      q = len (&W67->val);
	      //p = len (&W67->conval);
	      for (i = 1; i <= q; i++)
		{
		  covar->val.A[i] = W67->val.A[i] + W67->conval.A[1];
		}
	      for (i = 1; i <= n - q - 1; i++)
		{
		  covar->val.A[q + i] =
		    W67->conval.A[1] - W67->conval.A[n - q - i + 1];
		}
	      covar->conval.A[1] = -W67->conval.A[1];
	    }
	  else
	    covar->val = W67->val;
	  oadd (&newlist, &covar);
	  dispchr (&covar);
	  chrc = W67->next;
	}
    }
  else
    INAPPROPRIATE_GROUP ();
  return newlist;
}

ocharptr
fprod (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr R309;
  ocharptr newlist, temp, tchrc1, tchrc2, ntemp, otemp;
  int p, k, n;
  bool test;
  if (((grp.name == spnc) || (grp.name == sonc)))
    {
      newlist = NULL;

      n = grp.rank / 2;
      while (chrc1 != NULL)
	{
	  register ocharptr W72 = &(*chrc1);

	  p = 2 * chrc1->conval.A[1];
	  k = chrc2->val.A[maxdim];
	  if (chrc2->spin)
	    k = 2 * k + 1;
	  else
	    k = 2 * k;
	  k = k + p;
	  if ((bool) ((k) & 1))
	    test = true;
	  else
	    test = false;
	  k = k / 2;
	  cnu (&tchrc1);
	  tchrc1->mult = W72->mult;
	  tchrc1->val = W72->val;
	  tchrc1->C6_double = false;
	  cnu (&tchrc2);
	  tchrc2->mult = chrc2->mult;
	  tchrc2->val = chrc2->val;
	  tchrc2->C6_double = false;
	  temp = unitprod (tchrc1, tchrc2, false, n);
	  dispchr (&tchrc1);
	  dispchr (&tchrc2);
	  otemp = temp;
	  while (temp != NULL)
	    {
	      register ocharptr W73 = &(*temp);

	      cnu (&ntemp);
	      ntemp->mult = W73->mult;
	      ntemp->C6_double = true;
	      ntemp->spin = test;
	      ntemp->val = W73->val;
	      ntemp->conval = nolls;
	      ntemp->val.A[maxdim] = k;
	      ntemp->conval.A[1] = k;
	      oadd (&newlist, &ntemp);
	      dispchr (&ntemp);
	      temp = temp->next;
	    }
	  odisp (&otemp);
	  sposort (&newlist, true);
	  chrc1 = W72->next;
	}
    }
  else
    INAPPROPRIATE_GROUP ();

  R309 = newlist;
  return R309;
}

ocharptr
ffprod (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr R310;
  ocharptr brlist, colist, flist, tlist;
  int m, n;
  flist = NULL;
  if (((grp.name == spnc) || (grp.name == sonc)))
    {
      n = grp.rank;
      m = n / 2;
      /*if (generic(chrc1, chrc2, grp)) { */
      if ((grp.name == spnc))
	brlist = unsubgrbr (chrc1, n, 0, 'd');
      else
	brlist = unsubgrbr (chrc1, n, 0, 'b');

      fixcg (&currgrp.A[1 - 1], un, m, 0);
      colist = covariant (brlist, currgrp.A[1 - 1]);

      if ((grp.name == spnc))
	fixcg (&currgrp.A[1 - 1], spnc, n, 0);
      else
	fixcg (&currgrp.A[1 - 1], sonc, n, 0);
      odisp (&brlist);
      /* flist = sprch(fprod(colist, chrc2, grp)); */
      tlist = fprod (colist, chrc2, grp);
      flist = sprch (tlist);
      odisp (&tlist);
      odisp (&colist);
      /*} else {
         (void)fprintf(output.fp, "Non-generic product: product set to zero\n"), Putl(output, 1);
         if (logging)
         (void)fprintf(logfile.fp, "Non-generic product: product set to zero\n"), Putl(logfile, 1);
         } */
    }
  else
    INAPPROPRIATE_GROUP ();

  R310 = flist;
  return R310;
}

ocharptr
sprch (ocharptr list)
{
  register ocharptr R311;
  ocharptr newlist, temp;
  int n, w;
  register int k;

  newlist = NULL;
  n = currgrp.A[1 - 1].rank / 2;

  if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == sonc))
    {
      while (list != NULL)
	{
	  register ocharptr W74 = &(*list);

	  cnu (&temp);
	  temp->mult = W74->mult;
	  temp->spin = W74->spin;
	  temp->C6_double = W74->C6_double;
	  temp->val = W74->val;
	  temp->conval = W74->conval;
	  w = W74->val.A[n];

	  if (w > 0)
	    {
	      for (k = 1; k <= n; k++)
		{
		  temp->val.A[k] = W74->val.A[k] - w;
		}
	      temp->val.A[maxdim] = W74->val.A[maxdim] + w;
	      temp->conval.A[1] = W74->conval.A[1] + w;
	    }
	  oadd (&newlist, &temp);
	  dispchr (&temp);
	  list = W74->next;
	}
      sposort (&newlist, true);
    }
  else
    INAPPROPRIATE_GROUP ();
  R311 = newlist;
  return R311;
}

bool
highlystandard (ocharptr chrc, groop grp)
{
  register bool R312;
  bool test = false;		// modified by FB, was uninitialized
  int k, n;
  termptr temp;

  if (((grp.name == spnc) || (grp.name == sonc)))
    {
      test = false;
      n = grp.rank / 2;
      k = 2 * chrc->val.A[maxdim];
      snu (&temp);
      temp->mult = 1;
      temp->val = chrc->val;
      temp->val.A[maxdim] = 0;
      conjgte (&temp->val);
      if ((grp.name == spnc))
	{
	  if (chrc->spin)
	    k = k + 1;

	  if ((((k < n) && (temp->val.A[2] == 0))
	       || (n <= k - temp->val.A[2])))
	    /*2/5/99 */
	    test = true;
	}
      if ((grp.name == sonc))
	if ((k + 2 - temp->val.A[1] > MIN (n, k)))	/*2/5/99 */
	  test = true;
      dispsfn (&temp);
    }
  else
    INAPPROPRIATE_GROUP ();
  R312 = test;
  return R312;
}

bool
hstd (ocharptr chrc, groop grp)
{
  bool test = false;

  if (((grp.name == spnc) || (grp.name == sonc)))
    {
      test = true;
      while (chrc != NULL)
	{
	  register ocharptr W77 = &(*chrc);

	  if ((highlystandard (chrc, grp) == false))
	    test = false;
	  chrc = W77->next;
	}
    }
  else
    INAPPROPRIATE_GROUP ();
  return test;
}
