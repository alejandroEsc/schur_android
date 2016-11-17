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
# include <string.h>
# include <unistd.h>
# include <limits.h>

# include "standard.h"
# include "define.h"
# include "ReadWrite.h"
# include "sets_mgmt.h"
# include "Scanck.h"
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "s2.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "g.h"
# include "skew.h"
# include "s6.h"

static FILE *Tmpfil;
static long Tmplng;
groop *G219_tog;
int *G217_op;

termptr
osnmod (termptr sfnn, int parts, int k)
{
  enum
  { unmodified, modified, null } partition;

  register termptr R197;
  int depth, hook, hookpart;

  hook = 2 * parts - k;
  depth = 1;
  partition = modified;
  if (hook > 0)
    do
      {
	partition = unmodified;
	hookpart = sfnn->val.A[parts] - depth + 1;
	if (hook < hookpart)
	  {
	    partition = null;
	    sfnn = NULL;
	  }
	else if (hook == hookpart)
	  {
	    if (!(bool) ((sfnn->val.A[parts]) & 1))
	      sfnn->mult = -sfnn->mult;
	    sfnn->val.A[parts] = sfnn->val.A[parts] - hookpart;
	    if (sfnn->slab == ' ')
	      sfnn->slab = '#';
	    else
	      sfnn->slab = ' ';
	    partition = modified;
	  }
	else
	  {
	    depth = sfnn->val.A[parts];
	    hook = hook - hookpart;
	    sfnn->val.A[parts] = sfnn->val.A[parts] - hookpart;
	    parts = parts - 1;
	    if (parts < 1)
	      {
		partition = null;
		sfnn = NULL;
	      }
	  }
      }
    while (!(partition != unmodified));
  R197 = sfnn;
  return R197;
}


termptr
signseq (termptr * list, int k, int n, char tag, bool ppp, bool song)
{
  enum
  { zero, self, star, nostar } assoc;

  register termptr R198;
  termptr tlist, ones, temp, slist, plist;
  int parts, mm, nn, kk;
  bool tboole;
  register int i;

  if (ppp)
    nn = min (n, k);
  else
    nn = n;
  tboole = iosup;		/*24-10-97 */
  kk = k / 2;
  iosup = false;		/*true */
  mm = (*list)->mult;
  assoc = zero;
  if ((((bool) ((k) & 1) && (tag == 'a')) || (tag == 'c') || (tag == 'g')
       || (tag == 'e')))
    {
      if (((*list)->slab == '#'))
	assoc = star;
      if (((*list)->slab == ' '))
	assoc = nostar;
      parts = len (&(*list)->val);
      if (((!(bool) ((k) & 1)) && (parts < kk) && ((*list)->slab == ' ')))
	assoc = nostar;
      if (((!(bool) ((k) & 1)) && (parts == kk)))
	assoc = self;
      if ((parts > kk))
	{
	  (*list)->slab = ' ';
	  (*list) = osnmod ((*list), parts, k);
	  if ((*list) == NULL)
	    print ( "ERROR: Input Sfn modifies to nil - no signed sequence\n");
	  else
	    {
	      parts = len (&(*list)->val);
	      if ((*list)->slab == '#')
		assoc = star;
	      else
		if ((((*list)->slab == ' ') && (!(bool) ((k) & 1))
		     && (parts < kk)))
		assoc = nostar;
	      else if (((!(bool) ((k) & 1)) && (parts == kk)))
		assoc = self;
	      if (((bool) ((k) & 1) && ((*list)->slab == ' ')))
		assoc = nostar;
	    }
	}
    }
  tlist = NULL;
  if ((*list) != NULL)
    {
	if (tag=='t')
          plist = seriesx (tag, n*(n+1), full); // n(n+1) was setlimit FB 2006.
	else
	  plist = allseriesbut_t (tag, n*(n+1), full, n);

      slist = sfnmult (mm, plist);
      ldisp (&plist);
      schur_restrict (&slist, nn, 'l');
      if ((!song))
	{
	  if (((tag == 'c') || ((bool) ((k) & 1) && (tag == 'a'))))
	    {
	      if ((assoc == star))
		schur_restrict (&slist, 1, 'q');
	      else if ((assoc == nostar))
		schur_restrict (&slist, 1, 'p');
	    }
	  if ((tag == 'g'))
	    {
	      if ((assoc == star))
		schur_restrict (&slist, 1, 'o');
	      else
		schur_restrict (&slist, 1, 'e');
	    }
	}
      snu (&ones);
      ones->mult = 1;
      ones->val = nolls;
      for (i = 1; i <= kk; i++)
	{
	  ones->val.A[i] = 1;
	}
      tlist = attach ((*list), ones, true);
      ldisp (&(*list));
      temp = insert (tlist, slist);
      ldisp (&tlist);
      ldisp (&slist);
      tlist = attach (ones, temp, false);
      ldisp (&temp);
      dispsfn (&ones);
      stndise (&tlist);
      schur_restrict (&tlist, nn, 'l');

      sort (&tlist, true);
    }
  R198 = tlist;
  iosup = tboole;		/*false; true; */
  return R198;
}

ocharptr
signseqgr (ocharptr list, int parts, int q)
{
  register ocharptr R199;
  ocharptr newlist, temp, tlist;

  newlist = NULL;
  while (list != NULL)
    {
      register ocharptr W5 = &(*list);

      cnu (&temp);
      temp->mult = W5->mult;
      temp->val = W5->val;
      temp->spin = W5->spin;
      temp->lab = W5->lab;
      temp->conlab = W5->conlab;
      temp->C6_double = W5->C6_double;
      temp->conval = W5->conval;
      if (currgrp.A[1 - 1].name == un)
	tlist = upqseq (temp, currgrp.A[1 - 1].rank, parts, q);
      else if (currgrp.A[1 - 1].name == spn)
	tlist = signseqsp2k (temp, parts, currgrp.A[1 - 1].rank);
      else if (currgrp.A[1 - 1].name == son)
	tlist = signseqson (temp, parts, currgrp.A[1 - 1].rank);
      else if (currgrp.A[1 - 1].name == on)
	tlist = signseqon (temp, parts, currgrp.A[1 - 1].rank);
      else
	{
	  print ("ERROR: signseq not valid for this group\n");	// added FB&FT
	  return NULL;
	}

      oadd (&newlist, &tlist);
      dispchr (&temp);
      list = W5->next;
    }
  osort (&newlist, true);
  R199 = newlist;
  return R199;
}

ocharptr
signseqon (ocharptr list, int parts, int n)
{
  ocharptr newlist, nlist;
  termptr slist, temp, newtemp;
  char tag;

  newlist = NULL;
  if (list->spin)
    if ((bool) ((n) & 1))
      tag = 'a';
    else
      tag = 'e';
  if ((!list->spin))
    if ((bool) ((n) & 1))
      tag = 'g';
    else
      tag = 'c';
  snu (&slist);
  slist->mult = list->mult;
  slist->val = list->val;
  if (list->lab == '#')
    slist->slab = '#';
  else
    slist->slab = ' ';
  temp = signseq (&slist, n, parts, tag, false, false);
  dispsfn (&slist);
  newtemp = temp;
  while (temp != NULL)
    {
      register termptr W6 = &(*temp);

      cnu (&nlist);
      nlist->mult = W6->mult;
      nlist->val = W6->val;
      nlist->spin = list->spin;
      nlist->C6_double = false;
      nlist->lab = ' ';
      oadd (&newlist, &nlist);
      temp = W6->next;
    }
  ldisp (&newtemp);
  osort (&newlist, true);
  return newlist;
}

ocharptr
signseqsp2k (ocharptr list, int parts, int n)
{
  register ocharptr R201;
  termptr slist, temp;
  ocharptr newlist, nlist;

  newlist = NULL;
  snu (&slist);
  slist->mult = list->mult;
  slist->val = list->val;
  if (list->lab == '#')
    slist->slab = '#';
  else
    slist->slab = ' ';
  temp = signseq (&slist, n, parts, 'a', false, false);
  dispsfn (&slist);
  while (temp != NULL)
    {
      register termptr W7 = &(*temp);

      cnu (&nlist);
      nlist->mult = W7->mult;
      nlist->val = W7->val;
      nlist->spin = false;
      nlist->C6_double = false;
      nlist->lab = ' ';
      oadd (&newlist, &nlist);
      temp = W7->next;
    }
  ldisp (&temp);
  osort (&newlist, true);
  R201 = newlist;
  return R201;
}

ocharptr
signseqso2k (ocharptr list, int parts, int n)
{
  register ocharptr R202;
  ocharptr newlist, templist, temp, newtemp, mlist, ttemp;
  char ch;

  newlist = NULL;
  {
    register ocharptr W8 = &(*list);

    cnu (&templist);
    templist->mult = W8->mult;
    templist->lab = ' ';
    templist->val = W8->val;
    templist->spin = false;
    templist->C6_double = false;
    ch = W8->lab;
  }
  temp = signseqon (templist, parts, n);
  dispchr (&templist);
  ttemp = temp;
  while (temp != NULL)
    {
      register ocharptr W9 = &(*temp);

      cnu (&newtemp);
      newtemp->mult = list->mult;
      newtemp->val = W9->val;
      newtemp->lab = '+';
      newtemp->spin = false;
      newtemp->C6_double = false;
      mlist = gmodify (newtemp, currgrp.A[1 - 1]);
      dispchr (&newtemp);
      mlist->val = W9->val;
      mlist->spin = false;
      mlist->C6_double = false;
      if ((mlist->lab == '+'))
	mlist->lab = ch;
      else if ((ch == '+'))
	mlist->lab = '-';
      else
	mlist->lab = '+';
      if ((mlist->mult != list->mult))
	mlist->mult = -list->mult;
      else
	mlist->mult = list->mult;
      oadd (&newlist, &mlist);
      temp = W9->next;
    }
  odisp (&ttemp);
  osort (&newlist, true);
  R202 = newlist;
  return R202;
}

ocharptr
signseqson (ocharptr list, int parts, int n)
{
  register ocharptr R203;
  int signx, wsign, w, wtemp;
  ocharptr newlist, nlist;
  termptr slist, temp, newtemp;
  char tag;

  newlist = NULL;
  if (((len (&list->val) == (n / 2)) && (list->spin == false)
       && (!(bool) ((n) & 1))))
    newlist = signseqso2k (list, parts, n);
  else
    {
      w = wtfrm (&list->val);
      if ((bool) ((w) & 1))
	wsign = 1;
      else
	wsign = -1;
      if (list->lab == '-')
	wsign = -wsign;
      if (list->spin)
	if ((bool) ((n) & 1))
	  tag = 'a';
	else
	  tag = 'e';
      if ((!list->spin))
	if ((bool) ((n) & 1))
	  tag = 'g';
	else
	  tag = 'c';
      snu (&slist);
      slist->mult = list->mult;
      slist->val = list->val;
      if (list->lab == '#')
	slist->slab = '#';
      else
	slist->slab = ' ';
      temp = signseq (&slist, n, parts, tag, false, true);
      dispsfn (&slist);
      newtemp = temp;
      while (temp != NULL)
	{
	  register termptr W10 = &(*temp);

	  cnu (&nlist);
	  nlist->mult = W10->mult;
	  nlist->val = W10->val;
	  nlist->spin = list->spin;
	  nlist->C6_double = false;
	  nlist->lab = ' ';
	  wtemp = wtfrm (&W10->val);
	  if ((!(bool) ((n) & 1)) && (list->spin))
	    {
	      if ((bool) ((wtemp) & 1))
		signx = wsign;
	      else
		signx = -wsign;
	      if (signx == -1)
		nlist->lab = '-';
	      else
		nlist->lab = '+';
	    }
	  oadd (&newlist, &nlist);
	  temp = W10->next;
	}
      ldisp (&newtemp);
      osort (&newlist, true);
    }
  R203 = newlist;
  return R203;
}

void
multsplit (termptr list)
{
  termptr dummy, ptr, oddt, teven, etemp, otemp;

  oddt = NULL;
  teven = NULL;
  etemp = NULL;
  otemp = NULL;
  while (list != NULL)
    {
      dummy = list->next;
      list->next = NULL;
      ptr = sfncopy (list);
      list->next = dummy;
      if (list->mult & 1) 
	{
	  if (oddt == NULL)
	    {
	      oddt = ptr;
	      otemp = oddt;
	    }
	  else
	    {
	      otemp->next = ptr;
	      otemp = otemp->next;
	    }
	}
      else
	{
	  if (teven == NULL)
	    {
	      teven = ptr;
	      etemp = teven;
	    }
	  else
	    {
	      etemp->next = ptr;
	      etemp = etemp->next;
	    }
	}
      ptr = NULL;
      list = list->next;
    }
  sort (&oddt, false);
  sort (&teven, false);
  ldisp (&svar.A[17 - 1]);
  ldisp (&svar.A[18 - 1]);
  svar.A[17 - 1] = oddt;
  svar.A[18 - 1] = teven;
}

void
chkskcmpt (termptr list, termptr * sfn)
{
  termptr top, topl, llist, helpful = NULL, help;
  bool okx, ok1;
  int j;

  top = (*sfn);
  topl = list;
  llist = list;
  while ((*sfn) != NULL)
    {
      okx = false;
      llist = topl;
      while ((llist != NULL) && (!okx))
	{
	  ok1 = true;
	  j = 1;
	  while (ok1 && (j <= maxdim))
	    {
	      ok1 = (bool) ((*sfn)->val.A[j] <= llist->val.A[j]);
	      j = j + 1;
	    }
	  okx = ok1;
	  llist = llist->next;
	}
      if (okx)
	{
	  helpful = (*sfn);
	  (*sfn) = (*sfn)->next;
	}
      else if ((*sfn) == top)
	{
	  help = (*sfn);
	  (*sfn) = (*sfn)->next;
	  top = (*sfn);
	  dispsfn (&help);
	}
      else
	{
	  help = (*sfn);
	  helpful->next = (*sfn)->next;
	  (*sfn) = (*sfn)->next;
	  dispsfn (&help);
	}
    }
  (*sfn) = top;
}

void
spinsplit (ocharptr list)
{
  ocharptr dummy, ptr, tensor, spinor, stemp, ttemp;

  tensor = NULL;
  spinor = NULL;
  stemp = NULL;
  ttemp = NULL;
  while (list != NULL)
    {
      dummy = list->next;
      list->next = NULL;
      ptr = chrccopy (list);
      list->next = dummy;
      if (list->spin)
	{
	  if (spinor == NULL)
	    {
	      spinor = ptr;
	      stemp = spinor;
	    }
	  else
	    {
	      stemp->next = ptr;
	      stemp = stemp->next;
	    }
	}
      else
	{
	  if (tensor == NULL)
	    {
	      tensor = ptr;
	      ttemp = tensor;
	    }
	  else
	    {
	      ttemp->next = ptr;
	      ttemp = ttemp->next;
	    }
	}
      ptr = NULL;
      list = list->next;
    }
  osort (&tensor, false);
  osort (&spinor, false);
  odisp (&vari.A[1 - 1]);
  odisp (&vari.A[2 - 1]);
  vari.A[1 - 1] = tensor;
  vari.A[2 - 1] = spinor;
}

termptr
attach (termptr partat, termptr list, bool addx)
{
  register termptr R204;
  frame f;
  termptr result, ptr = NULL;	// was uninitialized. FB
  int j;
  register int i;

  result = NULL;
  j = len (&partat->val);
  while (list != NULL)
    {
      if (result == NULL)
	{
	  snu (&result);
	  ptr = result;
	}
      else
	{
	  snu (&ptr->next);
	  ptr = ptr->next;
	}
      ptr->next = NULL;
      f=nolls;
      i = 1;
      while (list->val.A[i] != 0)
	{
	  f.A[i] = list->val.A[i];
	  i = i + 1;
	}
      if (addx)
	for (i = 1; i <= j; i++)
	  {
	    f.A[i] = f.A[i] + partat->val.A[i];
	  }
      else
	for (i = 1; i <= j; i++)
	  {
	    f.A[i] = f.A[i] - partat->val.A[i];
	  }
      ptr->val = f;
      ptr->mult = list->mult;
      list = list->next;
    }
  R204 = result;
  return R204;
}

ocharptr
rattach (ocharptr partra, ocharptr list, bool addx)
{
  frame f;
  ocharptr result, ptr = NULL;
  int j;
  register int i;

  result = NULL;
  j = len (&partra->val);
  while (list != NULL)
    {
      if (result == NULL)
	{
	  cnu (&result);
	  ptr = result;
	}
      else
	{
	  cnu (&ptr->next);
	  ptr = ptr->next;
	}
      ptr->next = NULL;
      for (i = 1; i <= maxdim; i++)
	{
	  f.A[i] = 0;
	}
      i = 1;
      while (list->val.A[i] != 0)
	{
	  f.A[i] = list->val.A[i];
	  i = i + 1;
	}
      if (addx)
	for (i = 1; i <= j; i++)
	  {
	    f.A[i] = f.A[i] + partra->val.A[i];
	  }
      else
	for (i = 1; i <= j; i++)
	  {
	    f.A[i] = f.A[i] - partra->val.A[i];
	  }
      ptr->val = f;
      ptr->mult = list->mult;
      ptr->spin = list->spin;
      ptr->lab = list->lab;
      ptr->C6_double = list->C6_double;
      list = list->next;
    }
  return result;
}

termptr
insert (termptr part, termptr list)
{
  int j, l;
  register int i;
  frame f;
  termptr result, ptr = NULL,	// was uninitialized. FB
    newlist, alist;

  newlist = NULL;
  while (part != NULL)
    {
      alist = list;
      //j = len (part->val); //the old way
      j = MAX( part->val.length, qlen(part->val) );
      result = NULL;
      while (alist != NULL)
	{
	  if (result == NULL)
	    {
	      snu (&result);
	      ptr = result;
	    }
	  else
	    {
	      snu (&ptr->next);
	      ptr = ptr->next;
	    }
	  ptr->next = NULL;
	  //l = len (alist->val); // same
	  l = MAX(alist->val.length,qlen(alist->val));
	  f = nolls;
	  if (l + j > maxdim - 1)
	    error (TOO_MANY_PARTS, 1);
	  else
	    {
	      for (i = 1; i <= j; i++)
		{
		  f.A[i] = part->val.A[i];
		}
	      for (i = 1; i <= l; i++)
		{
		  f.A[i + j] = alist->val.A[i];
		}
	    }
	  f.length=l+j; // added 20060301
	  ptr->val = f;
	  ptr->mult = alist->mult;
	  alist = alist->next;
	}
      add (&newlist, &result);
      ldisp (&alist);
      part = part->next;
    }
  return newlist;
}

termptr
deskewm (termptr * intermediate, int plethweight, int loweight)
{
  termptr extract, extr, ptr1, ptr2, tail = NULL, plethtail = NULL,
    plethresult;
  int w, nxt, old, nextwt, wt;

  wt = loweight;
  plethresult = NULL;
  do
    {
      extract = NULL;
      ptr1 = (*intermediate);
      ptr2 = ptr1;
      nextwt = INT_MAX;
      while (ptr1 != NULL)
	{
	  w = wtfrm (&ptr1->val);
	  if (w == wt)
	    {
	      extr = ptr1;
	      if (ptr1 != ptr2)
		{
		  ptr2->next = ptr1->next;
		  ptr1 = ptr2;
		}
	      else
		{
		  (*intermediate) = ptr1->next;
		  ptr1 = (*intermediate);
		  ptr2 = ptr1;
		}
	      if (extract != NULL)
		tail->next = extr;
	      else
		extract = extr;
	      tail = extr;
	      nxt = plethweight - w;
	      w = 0;
	      {
		register termptr W28 = &(*extr);

		do
		  {
		    w = w + 1;
		    old = nxt;
		    nxt = W28->val.A[w];
		    W28->val.A[w] = old;
		  }
		while (!(old == 0));
	      }
	    }
	  else
	    nextwt = min (w, nextwt);
	  ptr2 = ptr1;
	  if (ptr1 != NULL)
	    ptr1 = ptr1->next;
	}
      if (extract != NULL)
	{
	  tail->next = NULL;
	  ptr2 = seriesm (plethweight - wt - 1, true);
	  ptr1 = lskew (extract, ptr2);
	  merge (&(*intermediate), &ptr1, false, true);
	  ldisp (&ptr2);
	  if (plethresult != NULL)
	    plethtail->next = extract;
	  else
	    plethresult = extract;
	  plethtail = tail;
	}
      wt = nextwt;
    }
  while (!(*intermediate == NULL));
  return plethresult;
}

termptr
plethfactor (termptr list, termptr tterm)
{
  int molt, i, ww;
  register int j;
  termptr dummy, factlist, temp0, temp1, temp2, temp;
  frame mu;
  ww = wtfrm (&list->val) * wtfrm (&tterm->val);
  if (plwt < ww)
    return NULL;
  else
    {
      i = len (&tterm->val);
      mu = tterm->val;
      conjgte (&mu);
      j = len (&mu);
      if ((i <= j))
	{
	  factlist = factorsof (tterm->mult, tterm->val);
	  temp2 = NULL;
	  temp = factlist;
	  while (factlist != NULL)
	    {
	      snu (&temp0);
	      temp0->mult = 1;
	      i = 1;
	      molt = factlist->mult;
	      while (factlist->val.A[i] != 0)
		{
		  snu (&temp1);
		  temp1->mult = 1;
		  temp1->val.A[1] = factlist->val.A[i];
		  dummy = plethysm (list, temp1);
		  dispsfn (&temp1);
		  temp1 = temp0;
		  temp0 = louter (temp1, dummy);
		  schur_restrict (&temp0, plwt, 'w');
		  ldisp (&temp1);
		  ldisp (&dummy);
		  i = i + 1;
		}
	      if (temp0 != NULL)
		{
		  dummy = temp0;
		  temp0 = sfnmult (molt, dummy);
		  ldisp (&dummy);
		  dummy = temp2;
		  temp2 = ladd (temp0, dummy);
		  ldisp (&dummy);
		  ldisp (&temp0);
		}
	      factlist = factlist->next;
	    }
	  ldisp (&temp);
	  return temp2;
	}
      else
	  return plethysm (list, tterm);
    }
}

termptr
plethysm (termptr l, termptr m)
{
  bool found, conjug;
  int c, ones, twos, n, k, lambdawt, lambdaln, ww;
  register int termp;
  termptr result, tail = NULL,
    ptr, mser, expansion, lambdap, mup, skewed;
  frame lambda, mu, rho;
  ww = wtfrm (&l->val) * wtfrm (&m->val);
  if ((plwt < ww))
    result = NULL;
  else
    {
      found = false;
      result = NULL;
      lambda = l->val;
      mu = m->val;
      c = m->mult;
      if (l->mult == -1)
	{
	  conjgte (&mu);
	  c = c * minusoneto (wtfrm (&mu));
	}
      if (mu.A[1] < 2)
	if (mu.A[1] == 0)
	  {
	    snu (&result);
	    found = true;
	  }
	else if (mu.A[2] == 0)
	  {
	    snu (&result);
	    found = true;
	    mu = lambda;
	  }
      if (!found && (lambda.A[1] < 2))
	if (lambda.A[1] == 0)
	  if (mu.A[2] != 0)
	    found = true;
	  else
	    {
	      snu (&result);
	      found = true;
	      mu = lambda;
	    }
	else if (lambda.A[2] == 0)
	  {
	    snu (&result);
	    found = true;
	  }
      if (result != NULL)
	{
	  result->val = mu;
	  result->mult = 1;
	  result->next = NULL;
	}
      if (!found && (wtfrm (&mu) == 2))
	if (lambda.A[2] == 0)
	  {
	    n = lambda.A[1];
	    k = mu.A[2];
	    rho = nolls;
	    rho.A[1] = 2 * n - k;
	    rho.A[2] = k;
	    result = NULL;
	    for (termp = 1; termp <= (n + 2 - k) / 2; termp++)
	      {
		snu (&ptr);
		ptr->val = rho;
		ptr->mult = 1;
		if (result == NULL)
		  result = ptr;
		else
		  tail->next = ptr;
		tail = ptr;
		rho.A[1] = rho.A[1] - 2;
		rho.A[2] = rho.A[2] + 2;
	      }
	    found = true;
	    tail->next = NULL;
	  }
	else if (lambda.A[1] == 1)
	  {
	    n = len (&lambda);
	    k = mu.A[2];
	    rho = nolls;
	    twos = n - k;
	    ones = twos + 2 * k + 1;
	    result = NULL;
	    for (termp = 1; termp <= twos; termp++)
	      {
		rho.A[termp] = 2;
	      }
	    if (k == 1)
	      {
		rho.A[twos + 1] = 1;
		rho.A[twos + 2] = 1;
	      }
	    for (termp = 1; termp <= (n + 2 - k) / 2; termp++)
	      {
		snu (&ptr);
		ptr->val = rho;
		ptr->mult = 1;
		if (result == NULL)
		  result = ptr;
		else
		  tail->next = ptr;
		tail = ptr;
		if (twos > 0)
		  rho.A[twos] = 1;
		if (twos > 1)
		  rho.A[twos - 1] = 1;
		twos = twos - 2;
		rho.A[ones] = 1;
		rho.A[ones + 1] = 1;
		ones = ones + 2;
	      }
	    found = true;
	    tail->next = NULL;
	  }
      if (!found)
	{
	  lambdawt = 0;
	  lambdaln = 0;
	  k = wtfrm (&mu);
	  while (lambda.A[lambdaln + 1] != 0)
	    {
	      lambdaln = lambdaln + 1;
	      lambdawt = lambdawt + lambda.A[lambdaln];
	    }
	  conjug = (bool) (lambda.A[1] > lambdaln);
	  if (conjug)
	    {
	      conjgte (&lambda);
	      if ((bool) ((lambdawt) & 1))
		conjgte (&mu);
	    }
	  if ((lambdawt == 2) && (mu.A[2] == 0))
	    {
	      conjug = !conjug;
	      result = sameweight (2 * k, 2);
	    }
	  else
	    {
	      snu (&lambdap);
	      snu (&mup);
	      lambdap->val = lambda;
	      lambdap->mult = abs (l->mult);
	      lambdap->next = NULL;
	      mup->val = mu;
	      mup->mult = m->mult;
	      mup->next = NULL;
	      mser = seriesm (lambda.A[1], false);
	      skewed = lskew (lambdap, mser);
	      ldisp (&mser);
	      dispsfn (&lambdap);
	      expansion = listplethterm (skewed, mup, true);
	      ldisp (&skewed);
	      dispsfn (&mup);
	      result =
		deskewm (&expansion, lambdawt * k,
			 k * (lambdawt - lambda.A[1]));
	    }
	  if (conjug)
	    result = lconjgte (result);
	}
      if (c != 1)
	coeffset (&result, c, '*');
    }
  return result;
}

termptr
pleth (frame a, frame b)
{
  register termptr tmp;
  termptr l, m;

  snu (&l);
  {
    register termptr W41 = &(*l);

    W41->mult = 1;
    W41->val = a;
    W41->next = NULL;
  }
  snu (&m);
  {
    register termptr W42 = &(*m);

    W42->mult = 1;
    W42->val = b;
    W42->next = NULL;
  }
  tmp = plethysm (l, m);
  dispsfn (&l);
  dispsfn (&m);
  return tmp;
}

termptr
listplethterm (termptr list, termptr tterm, bool suppress)
{
  int coe;
  termptr res, dres, pres, ptr, sigi, sigmas, sigmal, li;

  coe = tterm->mult;
  tterm->mult = 1;
  if (list->next == NULL)
    dres = plethfactor (list, tterm);
  else
    {
      sigmas = allskews (tterm->val, suppress);
      dres = NULL;
      sigi = sigmas;
      while (sigi != NULL)
	{
	  sigmal = skewx (1, tterm->val, sigi->val, maxdim - 1);
	  li = sigmal;
	  res = NULL;
	  while (li != NULL)
	    {
	      ptr = plethfactor (list, li);
	      merge (&res, &ptr, true, true);
	      li = li->next;
	    }
	  ldisp (&sigmal);
	  ptr = listplethterm (list->next, sigi, false);
	  /* if (ptr != NULL) { */
	  pres = louter (res, ptr);
	  ldisp (&ptr);
	  ldisp (&res);
	  schur_restrict (&pres, plwt, 'w');
	  if (pres != NULL)
	    merge (&dres, &pres, true, true);
	  /* } */
	  sigi = sigi->next;
	}
      ldisp (&sigmas);
    }
  tterm->mult = coe;
  coeffset (&dres, coe, '*');
  return dres;
}

termptr
listplethlist (termptr list1, termptr list2)
{
  termptr result, subl, cleavedlist;

  result = NULL;
  if ((list1 != NULL) && (list2 != NULL))
    {
      cleavedlist = sfncopy (list1);
      if (cleavedlist == NULL) // sfncopy include a call to sort... 
	return result;

      cleave (&cleavedlist);
      while (list2 != NULL)
	{
	  subl = listplethterm (cleavedlist, list2, false);
	  merge (&result, &subl, true, true);
	  list2 = list2->next;
	}
      ldisp (&cleavedlist);
    }
  return result;
}

ocharptr
change (ocharptr list, int op, groop fromg, groop tog)
{
  double valr[4 - 1 + 1];
  register int i;
  double a, b, c, d;
  ocharptr ptr, result;
  register ocharptr W43;

  if (op != -1)
    {
      result = chrccopy (list);
      ptr = result;
      while (ptr != NULL)
	{
	  {
	    W43 = ptr;

	    progress ();
	    for (i = 1; i <= 4; i++)
	      {
		valr[i - 1] = W43->val.A[i] + 0.5 * (unsigned) (W43->spin);
	      }
	    if ((fromg.name == son) && !(bool) (fromg.rank & 1))
	      {
		i = fromg.rank / 2;
		if ((W43->lab == '-'))
		  valr[i - 1] = -valr[i - 1];
	      }
	    a = valr[1 - 1];
	    b = valr[2 - 1];
	    c = valr[3 - 1];
	    d = valr[4 - 1];
	    for (i = 1; i <= 4; i++)
	      {
		valr[i - 1] = 0;
	      }
	    switch ((int) (op))
	      {
	      case 0:
		valr[1 - 1] = a;
		break;
	      case 1:
		valr[1 - 1] = a + b;
		valr[2 - 1] = a - b;
		break;
	      case 2:
		valr[1 - 1] = a + b;
		valr[2 - 1] = a - c;
		valr[3 - 1] = b - c;
		break;
	      case 3:
		valr[1 - 1] = (a + b + c + d) / ((double) 2);
		valr[2 - 1] = (a + b - c - d) / ((double) 2);
		valr[3 - 1] = (a - b + c - d) / ((double) 2);
		valr[4 - 1] = (-a + b + c - d) / ((double) 2);
		break;
	      case 4:
		valr[1 - 1] = a / ((double) 2);
		break;
	      case 5:
		valr[1 - 1] = (a + b - c) / ((double) 2);
		valr[2 - 1] = (a - b + c) / ((double) 2);
		valr[3 - 1] = (a - b - c) / ((double) 2);
		break;
	      case 6:
		valr[1 - 1] = 2 * a;
		break;
	      case 7:
		valr[1 - 1] = (a + b) / ((double) 2);
		valr[2 - 1] = (a - b) / ((double) 2);
		break;
	      default:
		Caseerror (Line);
	      }
	    W43->spin = false;
	    for (i = 1; i <= 4; i++)
	      {
		W43->val.A[i] = Trunc (valr[i - 1]);
		a = W43->val.A[i];
		W43->spin = (bool) (W43->spin || ((valr[i - 1] - a) != 0));
	      }
	    W43->lab = ' ';
	    if (tog.name == son)
	      {
		i = tog.rank / 2;
		if (valr[i - 1] < 0)
		  {
		    W43->val.A[i] = abs (W43->val.A[i]);
		    W43->lab = '-';
		  }
		else if (valr[i - 1] != 0)
		  if (! (tog.rank & 1))
		    W43->lab = '+';
	      }
	  }
	  ptr = ptr->next;
	}
    }
  if (tog.rank > 1)
    ostndise (&result);
  osort (&result, false);
  return result;
}

void
cf (grptype grp, int r, int o)
{
  if ((G219_tog->name == grp) && (G219_tog->rank == r))
    (*G217_op) = o;
}

bool
readauto (string0 * buffx, int *p, int *op, groop * fromg,
	  groop * tog, int grno)
{
  register bool R215;
  bool bad;
  int oldn;
  groopArray oldg;
  groop altergp;
  int *F218;
  groop *F220;

  F220 = G219_tog;
  G219_tog = &(*tog);
  F218 = G217_op;
  G217_op = &(*op);
  oldg = currgrp;
  oldn = nprod;
  skipbl ((*buffx), &(*p));
  getgroup (&(*buffx), (*p), false);
  altergp = currgrp.A[1 - 1];
  currgrp = oldg;
  nprod = oldn;
  currgrp.A[grno - 1] = altergp;
  bad = erred;
  erred = false;
  if (!bad)
    {
      (*fromg) = oldg.A[grno - 1];
      (*G219_tog) = currgrp.A[grno - 1];
      (*G217_op) = -1;
      if (grno <= nprod)
	switch ((int) (fromg->name))
	  {
	  case son:
	    switch ((int) (fromg->rank))
	      {
	      case 2:
		cf (un, 1, 0);
		break;
	      case 3:
		cf (sung, 2, 6);
		cf (spn, 2, 6);
		break;
	      case 4:
		break;
	      case 5:
		cf (spn, 4, 1);
		break;
	      case 6:
		cf (sung, 4, 2);
		break;
	      case 8:
		cf (son, 8, 3);
		break;
	      default:
		error (MISTAKE, (*p));
	      }
	    break;
	  case sung:
	    switch ((int) (fromg->rank))
	      {
	      case 2:
		cf (son, 3, 4);
		cf (spn, 2, 0);
		break;
	      case 4:
		cf (son, 6, 5);
		break;
	      default:
		error (MISTAKE, (*p));
	      }
	    break;
	  case un:
	    switch ((int) (fromg->rank))
	      {
	      case 1:
		cf (son, 2, 4);
		break;
	      default:
		error (MISTAKE, (*p));
	      }
	    break;
	  case spn:
	    switch ((int) (fromg->rank))
	      {
	      case 2:
		cf (son, 3, 4);
		cf (sung, 2, 0);
		break;
	      case 4:
		cf (son, 5, 7);
		break;
	      default:
		error (MISTAKE, (*p));
	      }
	    break;
	  default:
	    error (MISTAKE, (*p));
	  }
      if ((*G217_op) == -1)
	{
	  bad = true;
	  error (NOT_IMPLEMENTED, (*p));
	}
    }
  else
    error (MISTAKE, (*p));
  if (bad)
    {
      currgrp = oldg;
      nprod = oldn;
    }
  ggroup = currgrp.A[grno - 1];
  do
    {
      (*p) = (*p) + 1;
    }
  while (!((buffx->A[(*p) - 1] == delim) || ((*p) == bcol)));
  if (buffx->A[(*p) - 1] == delim)
    {
      (*p) = (*p) + 1;
      skipbl ((*buffx), &(*p));
    }
  else
    error (MISSING_COMMA, (*p));
  R215 = (bool) (!bad);
  G217_op = F218;
  G219_tog = F220;
  return R215;
}

void
cancel (int brn)
{
  switch ((int) (brn))
    {
    case 43:
      if (f4load)
	{
	  free (f4tab);
	  f4load = false;
	  odisp (&f4index);
	}
      break;
    case 44:
      if (e6load)
	{
	  free (e6tab);
	  e6load = false;
	  odisp (&e6index);
	}
      break;
    case 47:
      if (e7load)
	{
	  free (e7tab);
	  e7load = false;
	  odisp (&e7index);
	}
      break;
    case 49:
      if (e8load)
	{
	  free (e8tab);
	  e8load = false;
	  odisp (&e8index);
	}
      break;
    case 50:
      if (e8soload)
	{
	  free (e8sotab);
	  e8soload = false;
	  odisp (&e8soindex);
	}
      break;
    case 51:
      if (e8suload)
	{
	  free (e8sutab);
	  e8suload = false;
	  odisp (&e8suindex);
	}
      break;
    case 52:
      if (e8e6load)
	{
	  free (e8e6tab);
	  e8e6load = false;
	  odisp (&e8e6index);
	}
      break;
    case 45:
      if (e6soload)
	{
	  free (e6sotab);
	  e6soload = false;
	  odisp (&e6soindex);
	}
      break;
    case 48:
      if (e7e6load)
	{
	  free (e7e6tab);
	  e7e6load = false;
	  odisp (&e7e6index);
	}
      break;
    case 46:
      if (e6g2load)
	{
	  free (e6g2tab);
	  e6g2load = false;
	  odisp (&e6g2index);
	}
      break;
    case 53:
      if (u27e6load)
	{
	  free (u27e6tab);
	  u27e6load = false;
	  odisp (&u27e6index);
	}
      break;
    case 54:
      if (su56e7load)
	{
	  free (su56e7tab);
	  su56e7load = false;
	  odisp (&su56e7index);
	}
      break;
    case 55:
      if (su248e8load)
	{
	  free (su248e8tab);
	  su248e8load = false;
	  odisp (&su248e8index);
	}
      break;
    case 56:
      if (e6f4load)
	{
	  free (e6f4tab);
	  e6f4load = false;
	  odisp (&e6f4index);
	}
      break;
    case 57:
      if (f4g2load)
	{
	  free (f4g2tab);
	  f4g2load = false;
	  odisp (&f4g2index);
	}
      break;
    case 58:
      if (e6su3g2load)
	{
	  free (e6su3g2tab);
	  e6su3g2load = false;
	  odisp (&e6su3g2index);
	}
      break;
    case 63:
      if (l168load)
	{
	  free (l168tab);
	  l168load = false;
	  odisp (&l168index);
	}
      break;
    default:
      Caseerror (Line);
    }
}

void
saveframe (text * t, frame fr)
{
  register int i, lastStep;

  lastStep = qqlen (fr);	/*28/6/97 */
  for (i = 1; i <= lastStep; i++)
    {
      (void) fprintf ((*t).fp, " %1d", fr.A[i]), Putl ((*t), 0);
    }
  Putchr ('\n', (*t));
}

void
loadframe (text * t, int size, frame * fr)
{
  register int i;

  (*fr) = nolls;
  for (i = 1; i <= size; i++)
    {
      Fscan ((*t)), Scan ("%ld", &Tmplng), fr->A[i] = Tmplng, Getx ((*t));
    }
}

void
savecom (string0 bufferx, int px)
{
  bool svs, justone;
  int theone, size;
  register int i;
  termptr ptr;
  ocharptr cptr;
  char wd[MAXSTRING], hit;
  string0 filen;
  text t;

  t.init = 0;t.fp=NULL;t.eoln=0;t.buf='\0';
  skipbl (bufferx, &px);
  readword (&bufferx, &px, wd);
  svs = interp ("svar", wd, 2);
  if (!svs)
    if (!interp ("rvar", wd, 2))
      error (MISTAKE, px);
  if (!erred)
    {
      hit = skipbl (bufferx, &px);
      justone = false;
      theone = 0;
      if (Member ((unsigned) (hit), numbers.S))
	{
	  readint (&bufferx, &px, &theone);
	  if ((theone <= 0) || ((theone > svarlimit) && svs)
	      || ((theone > rvarlimit) && !svs))
	    error (MISTAKE, px);
	  if (!erred)
	    if (svs)
	      erred = (bool) (svar.A[theone - 1] == NULL);
	    else
	      erred = (bool) (vari.A[theone - 1] == NULL);
	  justone = (bool) (!erred);
	}
      if (!erred)
	{
	  hit = skipbl (bufferx, &px);
	  readfilename (bufferx, &px, &filen);
	}
      if (!erred)
	{
	  if (overwrite (filen))
	    {
	      Rewritex (t, filen.A, sizeof (filen.A));
	      fprintf (output.fp, "save["), Putl (output, 0);
	      if (svs) // if svar
		{
		  fprintf (t.fp, "svar\n"), Putl (t, 1);
		  for (i = 1; i <= svarlimit; i++)
		    {
		      if ((svar.A[i - 1] != NULL)
			  && (!justone || (i == theone)))
			{
			  fprintf (t.fp, "%1d\n", i), Putl (t, 1);
			  ptr = svar.A[i - 1];
			  fprintf (output.fp, "%1d", i),
			    Putl (output, 0);
			  while (ptr != NULL)
			    {
			      size = len (&ptr->val);
			      fprintf (t.fp, "%1d\n", ptr->mult),
				Putl (t, 1);
			      fprintf (t.fp, "%1d", size), Putl (t, 0);
			      saveframe (&t, ptr->val);
			      ptr = ptr->next;
			    }
			  Putchr ('0', t), Putchr ('\n', t);
			}
		    }
		}
	      else // save a rvar
		{
		  fprintf (t.fp, "rvar\n"), Putl (t, 1);
		  for (i = 1; i <= rvarlimit; i++)
		    {
		      if ((vari.A[i - 1] != NULL)
			  && (!justone || (i == theone)))
			{
			  fprintf (t.fp, "%1d\n", i), Putl (t, 1);
			  cptr = vari.A[i - 1];
			  fprintf (output.fp, "%1d", i),
			    Putl (output, 0);
			  while (cptr != NULL)
			    {
			      {
				register ocharptr W58 = &(*cptr);

				size = qqlen (cptr->val);	/*26/6/97 */
				fprintf (t.fp, "%1d\n", W58->mult),
				  Putl (t, 1);
				fprintf (t.fp, "%1d", size),
				  Putl (t, 0);
				saveframe (&t, W58->val);
				size = qqlen (cptr->conval);	/*26/6/97 */
				fprintf (t.fp, "%1d", size),
				  Putl (t, 0);
				saveframe (&t, W58->conval);
				if ((W58->C6_double == false))
				  W58->conlab = spc;
				fprintf (t.fp, "%2d%2d%4d%4d\n",
						(unsigned) (W58->spin),
						(unsigned) (W58->
							    C6_double),
						(unsigned) (W58->lab),
						(unsigned) (W58->conlab)),
				  Putl (t, 1);
			      }
			      cptr = cptr->next;
			    }
			  Putchr ('0', t), Putchr ('\n', t);
			}
		    }
		}
	      Putchr ('$', t), Putchr ('\n', t);
	      fclose (t.fp);	
	      fprintf (output.fp, "]\n");
	    }
	}
    }
  //fclose (t.fp);
}

void
loadcom (string0 bufferx, int px)
{
  bool svs, gox;
  char ch;
  int k, i, n, c, sp, db, lb, cl, size, filesize;
  termptr ptr1, ptr2;
  ocharptr cptr1, cptr2;
  char wd[MAXSTRING];
  string0 filen;
  text t;

  t.init = 0;
  skipbl (bufferx, &px);
  readword (&bufferx, &px, wd);
  svs = interp ("svar", wd, 2);
  if (!svs)
    if (!interp ("rvar", wd, 2))
      error (MISTAKE, px);
  if (!erred)
    {
      skipbl (bufferx, &px);
      readfilename (bufferx, &px, &filen);
      if (!erred)
	{
	  if (!findfile (filen))
	    error (NO_SUCH_FILE, px);
	  else
	    {
	      Resetx2 (&t, filen.A, -1);
	      i = 0;
	      do
		{
		  ch = Getchr (t), Getl (&t);
		  i = i + 1;
		}
	      while (!(ch == '$'));
	      filesize = i - 1;
	      fclose (t.fp);
	      Resetx2 (&t, filen.A, -1);
	      ch = Getchr (t), Getl (&t);
	      k = 1;
	      if (svs)
		gox = (bool) (locase (ch) == 's');
	      else
		gox = (bool) (locase (ch) == 'r');
	      if (!gox)
		fprintf (output.fp, "error:file has wrong variable type.\n"),
		  Putl (output, 1);
	      else
		fprintf (output.fp, "load["), Putl (output, 0);
	      while (gox && (k <= filesize - 1))
		{
		  Fscan (t), Scan ("%ld", &Tmplng), n = Tmplng, Getx (t), Getl (&t);
		  fprintf (output.fp, "%1d", n), Putl (output, 0);
		  k = k + 1;
		  if (svs)
		    {
		      ldisp (&svar.A[n - 1]);
		      ptr1 = NULL;
		      do
			{
			  Fscan (t), Scan ("%ld", &Tmplng), c =
			    Tmplng, Getx (t), Getl (&t);
			  k = k + 1;
			  if (c != 0)
			    {
			      snu (&ptr2);
			      ptr2->mult = c;
			      ptr2->next = NULL;
			      Fscan (t), Scan ("%ld", &Tmplng), size =
				Tmplng, Getx (t);
			      loadframe (&t, size, &ptr2->val);
			      k = k + 1;
			      if (ptr1 == NULL)
				svar.A[n - 1] = ptr2;
			      else
				ptr1->next = ptr2;
			      ptr1 = ptr2;
			    }
			}
		      while (!(c == 0));
		    }
		  else
		    {
		      odisp (&vari.A[n - 1]);
		      cptr1 = NULL;
		      do
			{
			  Fscan (t), Scan ("%ld", &Tmplng), c =
			    Tmplng, Getx (t), Getl (&t);
			  k = k + 1;
			  if (c != 0)
			    {
			      cnu (&cptr2);
			      {
				register ocharptr W59 = &(*cptr2);

				W59->mult = c;
				W59->next = NULL;
				Fscan (t), Scan ("%ld", &Tmplng), size =
				  Tmplng, Getx (t);
				loadframe (&t, size, &W59->val);
				Getl (&t);
				k = k + 2;
				Fscan (t), Scan ("%ld", &Tmplng), size =
				  Tmplng, Getx (t);
				loadframe (&t, size, &W59->conval);
				k = k + 1;
				Fscan (t), Scan ("%ld", &Tmplng), sp =
				  Tmplng, Scan ("%ld", &Tmplng), db =
				  Tmplng, Scan ("%ld", &Tmplng), lb =
				  Tmplng, Scan ("%ld", &Tmplng), cl =
				  Tmplng, Getx (t);
				W59->spin = (bool) (sp == 1);
				W59->C6_double = (bool) (db == 1);
				W59->lab = lb;
				W59->conlab = cl;
			      }
			      if (cptr1 == NULL)
				vari.A[n - 1] = cptr2;
			      else
				cptr1->next = cptr2;
			      cptr1 = cptr2;
			    }
			}
		      while (!(c == 0));
		    }
		}
	      fclose (t.fp);
	      fprintf (output.fp, "]\n");
	    }
	}
    }
  //fclose (t.fp);			/*13/12/95 */
}

void
contragsplit (ocharptr list)
{
  ocharptr cov, tcov, selfc, tselfc, conc, tconc, ptr, temp, dummy;
  int w1, w2, l1, l2;

  cov = NULL;
  tcov = NULL;
  tselfc = NULL;
  selfc = NULL;
  conc = NULL;
  tconc = NULL;
  while (list != NULL)
    {
      dummy = list->next;
      list->next = NULL;
      ptr = chrccopy (list);
      list->next = dummy;
      temp = contrag (ptr, currgrp.A[1 - 1]);
      w1 = wtfrm (&temp->val);
      w2 = wtfrm (&ptr->val);
      l1 = len (&temp->val);
      l2 = len (&ptr->val);
      dispchr (&temp);
      if ((w1 != w2) || (l1 != l2))
	{
	  if ((w1 > w2))
	    {
	      if ((cov == NULL))
		{
		  cov = ptr;
		  tcov = cov;
		}
	      else
		{
		  tcov->next = ptr;
		  tcov = tcov->next;
		}
	    }
	  else
	    {
	      if ((conc == NULL))
		{
		  conc = ptr;
		  tconc = conc;
		}
	      else
		{
		  tconc->next = ptr;
		  tconc = tconc->next;
		}
	    }
	}
      else
	{
	  if ((selfc == NULL))
	    {
	      selfc = ptr;
	      tselfc = selfc;
	    }
	  else
	    {
	      tselfc->next = ptr;
	      tselfc = tselfc->next;
	    }
	}
      ptr = NULL;
      list = list->next;
    }
  osort (&cov, false);
  osort (&conc, false);
  osort (&selfc, false);
  odisp (&vari.A[18 - 1]);
  odisp (&vari.A[19 - 1]);
  odisp (&vari.A[20 - 1]);
  vari.A[18 - 1] = cov;
  vari.A[19 - 1] = selfc;
  vari.A[20 - 1] = conc;
}

termptr
signsnseq (int n, termptr list)
{
  termptr newlist, temp;
  int w, l;
  register int i;
  register int s;

  newlist = list;
  while (list != NULL)
    {
      register termptr W60 = list;

      w = wtfrm (&W60->val);
      l = len (&W60->val);
      for (s = 1; s <= setlimit; s++)
	{
	  snu (&temp);
	  if (s & 1)
	    temp->mult = -W60->mult;
	  else
	    temp->mult = W60->mult;
	  temp->val = nolls;
	  temp->val.A[1] = n - w + 1;
	  for (i = 2; i <= s; i++)
	    {
	      temp->val.A[i] = W60->val.A[i - 1] + 1;
	    }
	  if ((s < l))
	    for (i = s + 1; i <= maxdim - 1; i++)
	      {
		temp->val.A[i] = W60->val.A[i];
	      }
	  add (&newlist, &temp);
	  dispsfn (&temp);
	}
      list = W60->next;
    }
  return newlist;
}

void
conjsplit (termptr list)
{
  termptr cov, selfc, conc, ptr, temp;
  /*int       l1, l2; *//*12/12/95 */

  cov = NULL;
  selfc = NULL;
  conc = NULL;
  while (list != NULL)
    {
      register termptr W67 = &(*list);

      snu (&ptr);
      ptr->mult = W67->mult;
      ptr->val = W67->val;
      snu (&temp);
      temp->mult = W67->mult;
      temp->val = W67->val;
      conjgte (&temp->val);
      if ((testord (&temp->val, &ptr->val) == EQUAL))
	add (&selfc, &ptr);
      else if ((testord (&temp->val, &ptr->val) == GREATER))
	add (&conc, &ptr);
      else
	add (&cov, &ptr);
      dispsfn (&ptr);
      dispsfn (&temp);
      list = W67->next;
    }
  sort (&cov, true);
  sort (&selfc, true);
  sort (&conc, true);
  ldisp (&svar.A[18 - 1]);
  ldisp (&svar.A[19 - 1]);
  ldisp (&svar.A[20 - 1]);
  svar.A[18 - 1] = cov;
  svar.A[19 - 1] = selfc;
  svar.A[20 - 1] = conc;
}


ocharptr
upqseq (ocharptr list, int rank, int pp, int qq)
{
  register ocharptr R217;
  ocharptr newlist, temp;
  termptr slist, listv, listc, nlistv, nlistc, tlist;
  int p, q, x, w;
  register int i, lastStep;

  newlist = NULL;
  tlist = seriesx ('f', setlimit, full);
  slist = tlist;
  while (list != NULL)
    {
      register ocharptr W60 = &(*list);

      p = len (&W60->conval);
      q = len (&W60->val);
      x = rank - p - q;
      while (slist != NULL)
	{
	  w = wtfrm (&slist->val);
	  snu (&listv);
	  listv->mult = 1;
	  listv->val = slist->val;
	  snu (&listc);
	  listc->mult = 1;
	  listc->val = slist->val;
	  conjgte (&listc->val);
	  snu (&nlistv);
	  nlistv->mult = 1;
	  nlistv->val = W60->val;
	  snu (&nlistc);
	  nlistc->mult = 1;
	  nlistc->val = W60->conval;
	  for (i = 1; i <= x; i++)
	    {
	      nlistc->val.A[i + p] = 0;
	    }
	  lastStep = qlen (listc->val);
	  for (i = 1; i <= lastStep; i++)
	    {
	      nlistc->val.A[p + i + x] = listc->val.A[i];
	    }
	  lastStep = len (&listv->val);
	  for (i = 1; i <= lastStep; i++)
	    {
	      nlistv->val.A[q + i] = listv->val.A[i];
	    }
	  if (((len (&nlistv->val) <= qq) && (qlen (nlistc->val) <= pp)))
	    {
	      cnu (&temp);
	      if ((bool) ((w) & 1))
		temp->mult = -1;
	      else
		temp->mult = 1;
	      temp->C6_double = true;
	      temp->spin = false;
	      temp->lab = ' ';
	      temp->conlab = ' ';
	      temp->val = nlistv->val;
	      temp->conval = nlistc->val;
	      ostndise (&temp);
	      if (temp != NULL)
		oadd (&newlist, &temp);
	    }
	  dispsfn (&nlistv);
	  dispsfn (&nlistc);
	  dispsfn (&listv);
	  dispsfn (&listc);
	  slist = slist->next;
	}
      list = W60->next;
    }
  ldisp (&tlist);
  R217 = newlist;
  return R217;
}

termptr
expleth (termptr list)
{
  termptr newone, C30_pow, sf, pone;
  int k, l, p;
  register int i,j;

  j = 1;
  snu (&sf);
  sf->mult = 1;
  sf->val = nolls;
  l = len (&list->val);
  while ((j <= l))
    {
      k = list->val.A[j];
      p = 0;
      for (i = 1; i <= l; i++)
	{
	  if ((list->val.A[i] == k))
	    p = p + 1;
	}
      snu (&newone);
      newone->mult = 1;
      newone->val = nolls;
      newone->val.A[1] = k;
      newone->val.A[2] = 2;
      snu (&C30_pow);
      C30_pow->mult = 1;
      C30_pow->val = nolls;
      C30_pow->val.A[1] = p;
      pone = plethysm (newone, C30_pow);
      dispsfn (&newone);
      dispsfn (&C30_pow);
      sf = louter (sf, pone);
      ldisp (&pone);
      j = j + p + 1;
    }
  return sf;
}
