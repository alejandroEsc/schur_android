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
# include "Scanck.h"

/*
**	Start of program definitions
*/
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "s2.h"
# include "s3.h"
# include "s6.h"
# include "bignums.h"
# include "dimensions.h"
# include "m.h"
# include "s.h"
# include "branch.h"
# include "s7.h"
# include "gr.h"
# include "g.h"
# include "skew.h"
# include "write.h"
# include "s5.h"

static FILE *Tmpfil;
static long Tmplng;

void
distinct (termptr * list)
{
  termptr newlist, temp;
  int i, length;
  bool test;

  newlist = NULL;
  while ((*list) != NULL)
    {
      register termptr W1 = *list;

      i = 1;
      test = true;
      length = len (&W1->val);
      if (length > 1)
	do
	  {
	    if ((W1->val.A[i] == W1->val.A[i + 1]))
	      test = false;
	    else
	      i = i + 1;
	  }
	while (!((W1->val.A[i] == 0) || (i == length) || (test == false)));
      temp = W1->next;
      if ((test == false))
	dispsfn (&(*list));
      else
	{
	  W1->next = newlist;
	  newlist = (*list);
	}
      (*list) = temp;
    }
  sort (&newlist, true);
  *list = newlist;
}

int
nlambda (frame x)
{
  int n;
  register int i, lastStep;

  n = 0;
  lastStep = len (&x);
  for (i = 1; i <= lastStep; i++)
      n = n + (i - 1) * x.A[i];
  return n;
}

int
pwtfrm (frame * partition)
{
  int i, weights;

  i = 1;
  weights = 0;
  while (i <= maxdim && partition->A[i] != 0)	// Modified by FB (add i<= maxdim)
    {
      weights = weights + (i - 1) * partition->A[i];	// Modified by FB, was i * partition->A[i]
      i = i + 1;
    }
  return weights;
}

int
hookpart (int i, int j, frame a)
{
  frame ac;

  ac = a;
  conjgte (&ac);
  return (a.A[i] + ac.A[j] - i - j + 1);
}

void
prodhook (frame a, lbframe * polyhk)
{
  int hk;
  register int j;
  register int i, lastStep, lastStep2;
  lbframe result;
  char ph;

  ph = '-';
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyhk->A[i] = 0;
    }
  polyhk->A[0] = 1;
  lastStep = len (&a);
  for (j = 1; j <= lastStep; j++)
    {
      lastStep2 = a.A[j];
      for (i = 1; i <= lastStep2; i++)
	{
	  hk = hookpart (j, i, a);
	  poly1 (hk, ph, &result);
	  polyprod ((*polyhk), result, &(*polyhk));
	}
    }
}

void
qqprodhook (frame a, lbframe * polyhk)
{
  register int k;
  register int j;
  register int i, lastStep, lastStep2;
  lbframe result;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      result.A[i] = 0;
    }
  result.A[0] = 1;
  lastStep = len (&a);
  for (j = 1; j <= lastStep; j++)
    {
      lastStep2 = a.A[j];
      for (i = 1; i <= lastStep2; i++)
	{
	  for (k = 0; k <= POLYMAXDEGREE; k++)
	    {
	      polyhk->A[k] = 0;
	    }
	  polyhk->A[i] = 1;
	  if (i != j)
	    polyhk->A[j] = 1;
	  else
	    polyhk->A[i] = 2;
	  polyprod (result, (*polyhk), &result);
	}
    }
  (*polyhk) = result;
}

void
qqsfn (termptr sfn, lbframe * qpoly)
{
  register int i;
  lbframe qtpoly, qqpoly;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      qpoly->A[i] = 0;
    }
  prodhook (sfn->val, &qtpoly);
  qqprodhook (sfn->val, &qqpoly);
  polyprod (qtpoly, qqpoly, &(*qpoly));
  for (i = 1; i <= POLYMAXDEGREE; i++)  //Added by FB  20060308
	 qpoly->A[i] = sfn->mult * qpoly->A[i];

}

void
qqsfnlist (termptr list)
{
  lbframe tpoly, lpoly;
  register int i;
  char signp;
  termptr tlist;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      tpoly.A[i] = 0;
    }
  signp = '+';
  tlist = list;
  while (tlist != NULL)
    {
      qqsfn (tlist, &lpoly);
      polysum (signp, tpoly, lpoly, &tpoly);
      tlist = tlist->next;
    }
  ldisp (&tlist);
  writepoly (tpoly);
}

void
qsfn (termptr sfn, lbframe * qpoly)
{
  int power;
  register int i;
  lbframe qtpoly;

  for (i = 0; i < POLYMAXDEGREE; i++)
    qpoly->A[i] = 0;

  power = pwtfrm (&sfn->val);
  if (power < POLYMAXDEGREE)
    {
      prodhook (sfn->val, &qtpoly);
      for (i = power; i <= POLYMAXDEGREE; i++)	//modified by FB & FT
	qpoly->A[i] = sfn->mult * qtpoly.A[i - power];
      for (i = 0; i < power; i++)
	qpoly->A[i] = 0;
    }
}

void
qsfnlist (termptr list)
{
  lbframe tpoly, lpoly;
  register int i;
  char signp;
  termptr tlist;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      tpoly.A[i] = 0;
    }
  signp = '+';
  tlist = list;
  while (tlist != NULL)
    {
      qsfn (tlist, &lpoly);
      polysum (signp, tpoly, lpoly, &tpoly);
      tlist = tlist->next;
    }
  ldisp (&tlist);
  writepoly (tpoly);
}

void
polyprod (lbframe a, lbframe b, lbframe * ab)
{
  register int j;
  register int i;

  progress ();
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ab->A[i] = 0;
    }
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      for (j = 0; j <= POLYMAXDEGREE - i; j++)
	{
	  ab->A[i + j] = ab->A[i + j] + a.A[i] * b.A[j];
	}
    }
}

void
polysum (char signp, lbframe apoly, lbframe bpoly, lbframe * abpoly)
{
  register int i;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      if ((signp == '+'))
	abpoly->A[i] = apoly.A[i] + bpoly.A[i];
      else
	abpoly->A[i] = apoly.A[i] - bpoly.A[i];
    }
}

void
polyn1 (char ph, lbframe * mpx)
{
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      mpx->A[i] = 0;
      ap.A[i] = 0;
      bp.A[i] = 0;
    }
  bp.A[0] = 1;
  ap.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE; i++)
    {
      if ((ph == '+'))
	ap.A[i] = 1;
      else
	ap.A[i] = -1;
      polyprod (ap, bp, &(*mpx));
      bp = (*mpx);
      ap.A[i] = 0;
    }
}

void
partpoly (lbframe * pp)
{
  lbframe apoly, bpoly;
  int k;
  register int i;
  char ph;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      apoly.A[i] = 0;
      bpoly.A[i] = 0;
      pp->A[i] = 0;
    }
  (void) fprintf (output.fp, "k ="), Putl (output, 0);
  Fscan (input), Scan ("%ld", &Tmplng), k =
    Tmplng, Getx (input), Getl (&input);
  (void) fprintf (output.fp, "enter sign ="), Putl (output, 0);
  ph = Getchr (input), Getl (&input);
  bpoly.A[0] = 1;
  polyn1 (ph, &apoly);
  for (i = 1; i <= k; i++)
    {
      polyprod (apoly, bpoly, &(*pp));
      bpoly = (*pp);
    }
}

void
poly1 (int p, char ph, lbframe * mpp)
{
  int q, k;
  register int i;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      mpp->A[i] = 0;
    }
  mpp->A[0] = 1;
  q = (30 / p) + 1;
  for (i = 1; i <= q; i++)
    {
      k = i * p;
      if ((k <= POLYMAXDEGREE))
	{
	  mpp->A[k] = 1;
	  if (((ph == '+') && ((bool) ((i) & 1) == true)))
	    mpp->A[k] = -mpp->A[k];
	}
    }
}

void
poly2 (char ph, lbframe * ppoly)
{
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ppoly->A[i] = 0;
      ap.A[i] = 0;
      bp.A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE; i++)
    {
      poly1 (i, ph, &ap);
      polyprod (ap, bp, &(*ppoly));
      bp = (*ppoly);
    }
}

void
poly3 (char ph, lbframe * ppoly)
{
  int k;
  register int j;
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      ppoly->A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE - 1; i++)
    {
      for (j = i + 1; j <= POLYMAXDEGREE - i; j++)
	{
	  k = i + j;
	  poly1 (k, ph, &ap);
	  polyprod (ap, bp, &(*ppoly));
	  bp = (*ppoly);
	}
    }
}

void
polybp (char ph, lbframe * ppoly)
{
  int k;
  register int j;
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      ppoly->A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE-1; i++)
    {
      for (j = i + 1; j <= POLYMAXDEGREE - i; j++)
	{
	  k = i + j - 1;
	  poly1 (k, ph, &ap);
	  polyprod (ap, bp, &(*ppoly));
	  bp = (*ppoly);
	}
    }
}

void
polyn3 (char ph, lbframe * polyn)
{
  register int j;
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      polyn->A[i] = 0;
    }
  bp.A[0] = 1;
  ap.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE-1; i++)
    {
      for (j = i + 1; j <= POLYMAXDEGREE - i; j++)
	{
	  ap.A[i + j] = 1;
	  if ((ph == '-'))
	    ap.A[i + j] = -1;
	  polyprod (ap, bp, &(*polyn));
	  ap.A[i + j] = 0;
	  bp = (*polyn);
	}
    }
}

void
polyn4 (char ph, lbframe * polyn)
{
  register int j;
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      polyn->A[i] = 0;
    }
  bp.A[0] = 1;
  ap.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE; i++)
    {
      for (j = i; j <= POLYMAXDEGREE - i; j++)
	{
	  ap.A[i + j] = 1;
	  if ((ph == '-'))
	    ap.A[i + j] = -1;
	  polyprod (ap, bp, &(*polyn));
	  ap.A[i + j] = 0;
	  bp = (*polyn);
	}
    }
}

void
poly4 (char ph, lbframe * ppoly)
{
  int k;
  register int j;
  register int i;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      ppoly->A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= POLYMAXDEGREE; i++)
    {
      for (j = i; j <= POLYMAXDEGREE - i; j++)
	{
	  k = i + j;
	  poly1 (k, ph, &ap);
	  polyprod (ap, bp, &(*ppoly));
	  bp = (*ppoly);
	}
    }
}

void
poly5 (char ph, lbframe * ppoly)
{
  int k, range;
  register int i;
  lbframe ap, bp;

  range = POLYMAXDEGREE / 2;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      ppoly->A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= range; i++)
    {
      k = 2 * i;
      poly1 (k, ph, &ap);
      polyprod (ap, bp, &(*ppoly));
      bp = (*ppoly);
    }
}

void
poly6 (char ph, lbframe * ppoly)
{
  int k, range;
  register int j;
  register int i;
  lbframe ap, bp;

  range = POLYMAXDEGREE / 2;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      bp.A[i] = 0;
    }
  bp.A[0] = 1;
  for (i = 1; i <= range; i++)
    {
      k = 2 * i;
      for (j = 0; j <= POLYMAXDEGREE; j++)
	{
	  ap.A[j] = 0;
	}
      ap.A[0] = 1;
      ap.A[k] = 1;
      if ((ph == '-'))
	ap.A[k] = -1;
      polyprod (ap, bp, &(*ppoly));
      bp = (*ppoly);
    }
}

void
mseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  poly2 (ch, &(*polyn));
}

void
pseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '+';
  poly2 (ch, &(*polyn));
}

void
aseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  polyn3 (ch, &(*polyn));
}

void
cseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  polyn4 (ch, &(*polyn));
}

void
bseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  poly3 (ch, &(*polyn));
}

void
dseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  poly4 (ch, &(*polyn));
}

void
eseries (lbframe * polyn)
{
  register int i;
  lbframe ap, bp;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      polyn->A[i] = 0;
    }
  ch = '-';
  polyn1 (ch, &ap);
  polyn3 (ch, &bp);
  polyprod (ap, bp, &(*polyn));
}

void
gseries (lbframe * polyn)
{
  register int i;
  lbframe ap, bp;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      ap.A[i] = 0;
      bp.A[i] = 0;
      polyn->A[i] = 0;
    }
  ch = '+';
  polyn1 (ch, &ap);
  ch = '-';
  polyn3 (ch, &bp);
  polyprod (ap, bp, &(*polyn));
}

void
fseries (lbframe * polyn)
{
  register int i;
  char ch;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  poly2 (ch, &ap);
  poly3 (ch, &bp);
  polyprod (ap, bp, &(*polyn));
}

void
hseries (lbframe * polyn)
{
  register int i;
  char ch;
  lbframe ap, bp;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '+';
  poly2 (ch, &ap);
  ch = '-';
  poly3 (ch, &bp);
  polyprod (ap, bp, &(*polyn));
}

void
lseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  polyn1 (ch, &(*polyn));
}

void
qseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '+';
  polyn1 (ch, &(*polyn));
}

void
qbseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '+';
  polyn3 (ch, &polyup);
  ch = '-';
  poly3 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  polyprod ((*polyn), (*polyn), &(*polyn));
  poly5 (ch, &polyin);
  ch = '+';
  poly6 (ch, &polyup);
  polyprod (polyin, (*polyn), &(*polyn));
  polyprod (polyup, (*polyn), &(*polyn));
}

void
qaseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '-';
  polyn3 (ch, &polyup);
  ch = '+';
  poly3 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  polyprod ((*polyn), (*polyn), &(*polyn));
  poly5 (ch, &polyin);
  ch = '-';
  poly6 (ch, &polyup);
  polyprod (polyin, (*polyn), &(*polyn));
  polyprod (polyup, (*polyn), &(*polyn));
}

void
qeseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '-';
  polyn1 (ch, &polyup);
  polyprod (polyup, polyup, &polyup);
  ch = '+';
  poly5 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  ch = '-';
  polyn3 (ch, &polyup);
  ch = '+';
  poly3 (ch, &polyin);
  polyprod (polyup, polyin, &polyup);
  polyprod (polyup, polyup, &polyup);
  polyprod (polyup, (*polyn), &(*polyn));
}

void
qfseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '-';
  poly2 (ch, &polyup);
  polyprod (polyup, polyup, &polyup);
  ch = '+';
  poly6 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  ch = '-';
  poly3 (ch, &polyup);
  ch = '+';
  polyn3 (ch, &polyin);
  polyprod (polyup, polyin, &polyup);
  polyprod (polyup, polyup, &polyup);
  polyprod (polyup, (*polyn), &(*polyn));
}

void
qgseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '-';
  polyn3 (ch, &polyup);
  ch = '+';
  poly3 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  polyprod ((*polyn), (*polyn), &(*polyn));
}

void
qhseries (lbframe * polyn)
{
  char ch;
  lbframe polyup, polyin;

  ch = '+';
  polyn3 (ch, &polyup);
  ch = '-';
  poly3 (ch, &polyin);
  polyprod (polyup, polyin, &(*polyn));
  polyprod ((*polyn), (*polyn), &(*polyn));
}

void
qlseries (lbframe * polyn)
{
  lbframe polyup, polyin;

  lseries (&polyup);
  pseries (&polyin);
  polyprod (polyup, polyin, &(*polyn));
}

void
qmseries (lbframe * polyn)
{
  lbframe polyup, polyin;

  mseries (&polyup);
  qseries (&polyin);
  polyprod (polyup, polyin, &(*polyn));
}

void
qpseries (lbframe * polyn)
{
  lbframe polyup, polyin;

  pseries (&polyup);
  lseries (&polyin);
  polyprod (polyup, polyin, &(*polyn));
}

void
qqqseries (lbframe * polyn)
{
  lbframe polyup, polyin;

  qseries (&polyup);
  mseries (&polyin);
  polyprod (polyup, polyin, &(*polyn));
}

void
bpseries (lbframe * polyn)
{
  register int i;
  char ch;
  for (i = 0; i <= POLYMAXDEGREE; i++)
    {
      polyn->A[i] = 0;
    }
  ch = '-';
  polybp (ch, &(*polyn));
}

bool
qvalidser (char ser)
{
  register bool R228;
  bool qv;

  qv = false;
  locase (ser);
  switch ((int) (ser))
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
      qv = true;
      break;
    default:
      Caseerror (Line);
    }
  R228 = qv;
  return R228;
}



void
qallseries (char series, lbframe * apoly)
{
  series = locase (series);
  if (qvalidser (series))
    {
      if ((series == 'm'))
	mseries (&(*apoly));
      else if ((series == 'p'))
	pseries (&(*apoly));
      else if ((series == 'b'))
	bseries (&(*apoly));
      else if ((series == 'd'))
	dseries (&(*apoly));
      else if ((series == 'f'))
	fseries (&(*apoly));
      else if ((series == 'h'))
	hseries (&(*apoly));
      else if ((series == 'l'))
	lseries (&(*apoly));
      else if ((series == 'q'))
	qseries (&(*apoly));
      else if ((series == 'a'))
	aseries (&(*apoly));
      else if ((series == 'e'))
	eseries (&(*apoly));
      else if ((series == 'g'))
	gseries (&(*apoly));
      else if ((series == 'c'))
	cseries (&(*apoly));
      writepoly ((*apoly));
    }
  else
    print ( "qseries %cnot available\n", series);
}

void
qqallseries (char series, lbframe * apoly)
{
  series = locase (series);
  if (qvalidser (series))
    {
      if ((series == 'm'))
	qmseries (&(*apoly));
      else if ((series == 'p'))
	qpseries (&(*apoly));
      else if ((series == 'b'))
	qbseries (&(*apoly));
      else if ((series == 'd'))
	qbseries (&(*apoly));
      else if ((series == 'f'))
	qfseries (&(*apoly));
      else if ((series == 'h'))
	qhseries (&(*apoly));
      else if ((series == 'l'))
	qlseries (&(*apoly));
      else if ((series == 'q'))
	qqqseries (&(*apoly));
      else if ((series == 'a'))
	qaseries (&(*apoly));
      else if ((series == 'e'))
	qeseries (&(*apoly));
      else if ((series == 'g'))
	qgseries (&(*apoly));
      else if ((series == 'c'))
	qaseries (&(*apoly));
      writepoly ((*apoly));
    }
  else
    print ( "qqseries %cnot available\n", series);
}

termptr
ntensor (termptr list, int n)
{
  register termptr R229;
  termptr newlist = NULL,	// was uninitialized. FB
    temp;
  register int i;

  if (n <= 0)
    newlist = NULL;
  else if (n == 1)
    newlist = list;
  else if (n >= 2)
    {
      snu (&temp);
      temp->mult = 1;
      for (i = 1; i <= maxdim; i++)
	{
	  temp->val.A[i] = 0;
	}
      for (i = 1; i <= n; i++)
	{
	  newlist = louter (temp, list);
	  ldisp (&temp);
	  temp = newlist;
	}
    }
  R229 = newlist;
  return R229;
}

ocharptr
rntensor (ocharptr list, int n)
{
  register ocharptr R230;
  ocharptr newlist = NULL,	// was uninitialized. FB
    temp;
  register int i;

  if (n <= 0)
    newlist = NULL;
  else if (n == 1)
    newlist = list;
  else if (n >= 2)
    {
      cnu (&temp);
      temp->mult = 1;
      for (i = 1; i <= maxdim; i++)
	{
	  temp->val.A[i] = 0;
	}
      for (i = 1; i <= n; i++)
	{
	  newlist = kronk (temp, list, currgrp.A[1 - 1]);
	  odisp (&temp);
	  temp = newlist;
	}
    }
  R230 = newlist;
  return R230;
}

termptr
fullsa (int n)
{
  register termptr R231;
  termptr sfn, sfn1, sfn2, temp;

  sfn1 = NULL;
  if (!(bool) ((n) & 1))
    {
      temp = fullx (n);
      sfn = temp;
      while (sfn != NULL)
	{
	  sfn2 = fulling (sfn);
	  add (&sfn1, &sfn2);
	  sort (&sfn1, true);
	  sfn = sfn->next;
	}
      ldisp (&temp);
    }
  R231 = sfn1;
  return R231;
}

termptr
fullx (int n)
{
  register termptr R232;
  termptr sfn1, sfn;
  int p;
  register int i;

  p = 1;
  snu (&sfn);
  sfn->mult = 1;
  for (i = 1; i <= maxdim; i++)
    {
      sfn->val.A[i] = 0;
    }
  sfn->val.A[1] = n;
  sfn1 = useseries ('f', sfn, true, false, p);
  dispsfn (&sfn);
  schur_restrict (&sfn1, -n, 'w');
  limit (&sfn1, 1);
  R232 = sfn1;
  return R232;
}

termptr
fulling (termptr a)
{
  register termptr R233;
  int i, j, m, n;
  register int q;
  register int p;
  register int k;
  termptr  sfn1, sfn2, sfn3, sfn4;

  i = len (&a->val);
  m = 1;
  j = 1;
  snu (&sfn1);
  sfn1->mult = 1;
  for (k = 1; k <= maxdim; k++)
    {
      sfn1->val.A[k] = 0;
    }
  do
    {
      k = a->val.A[m];
      while (a->val.A[j + 1] == k)
	j = j + 1;
      n = j - m + 1;
      snu (&sfn2);
      sfn2->mult = 1;
      for (p = 1; p <= maxdim; p++)
	{
	  sfn2->val.A[p] = 0;
	}
      sfn2->val.A[1] = k;
      sfn2->val.A[2] = 2;
      if (n > 1)
	{
	  snu (&sfn3);
	  sfn3->mult = 1;
	  for (q = 1; q <= maxdim; q++)
	    {
	      sfn3->val.A[q] = 0;
	    }
	  sfn3->val.A[1] = n;
	  sfn4 = plethysm (sfn2, sfn3);
	  dispsfn (&sfn2);
	  dispsfn (&sfn3);
	  sfn2 = sfn4;
	}
      sfn3 = louter (sfn1, sfn2);
      ldisp (&sfn1);
      ldisp (&sfn2);
      sfn1 = sfn3;
      j = j + 1;
      m = j;
    }
  while (!(j > i));
  p = -1;
  limit (&sfn1, p);
  R233 = sfn1;
  return R233;
}

void
fixxcg (groop * grp, grptype g, int r1, int r2)
{
  {
    register groop *W140 = &(*grp);

    W140->name = g;
    W140->rank = r1;
    W140->rank2 = r2;
  }
}

ocharptr
sp_pleth (ocharptr list1, ocharptr list2, groop grp)
{
  register ocharptr R234;
  ocharptr tlist, newlist, result;	/*10/05/99 */
  termptr slist1, slist2, slist12;
  char taga = 'a', tagb = 'b';	// was uninitialized. FB


  group = grp.name;
  if ((group == spn))
    {
      taga = 'a';
      tagb = 'b';
    }
  else if (((group == son) || (group == on)))
    {
      taga = 'c';
      tagb = 'd';
    }
  newlist = so2k1so3brnch (taga, list1, grp.rank, false);
  slist1 = repsfn (newlist);
  slist2 = repsfn (list2);
  odisp (&newlist);
  slist12 = listplethlist (slist1, slist2);
  ldisp (&slist1);
  ldisp (&slist2);
  tlist = sfntochrc (slist12, false, ' ');
  ldisp (&slist12);
  result = unospnbrnch (tagb, tlist, grp.rank);
  R234 = gmodify (result, grp);	/*10/05/99 */
  odisp (&tlist);
  odisp (&result);
  return R234;
}

ocharptr
g2pleth (ocharptr list1, ocharptr list2)
{
  register ocharptr R235;
  ocharptr list, newlist, tlist;
  termptr slist1, slist2, slist12;

  newlist = g2_so7 (list1);
  tlist = newlist;
  fixxcg (&currgrp.A[1 - 1], son, 7, 0);
  list = so2k1so3brnch ('c', newlist, currgrp.A[1 - 1].rank, false);
  fixxcg (&currgrp.A[1 - 1], sung, 7, 0);
  odisp (&tlist);
  slist1 = repsfn (list);
  slist2 = repsfn (list2);
  odisp (&list);
  slist12 = listplethlist (slist1, slist2);
  ldisp (&slist1);
  ldisp (&slist2);
  list = sfntochrc (slist12, false, ' ');
  ldisp (&slist12);
  tlist = list;
  newlist = unospnbrnch ('d', list, currgrp.A[1 - 1].rank);
  odisp (&tlist);
  tlist = newlist;
  fixxcg (&currgrp.A[1 - 1], son, 7, 0);
  list = so7g2brnch (newlist);
  odisp (&tlist);
  tlist = list;
  fixxcg (&currgrp.A[1 - 1], g2, 2, 0);
  /* gumodify(&list,2); */
  R235 = gmodify (list, currgrp.A[1 - 1]);
  odisp (&tlist);
  return R235;
}

ocharptr
unpleth (ocharptr list1, ocharptr list2, groop grp)
{
  register ocharptr R236;
  ocharptr tlist;
  termptr slist1, slist2, slist12;

  slist1 = repsfn (list1);
  slist2 = repsfn (list2);
  slist12 = listplethlist (slist1, slist2);
  ldisp (&slist1);
  ldisp (&slist2);
  tlist = sfntochrc (slist12, false, ' ');
  ldisp (&slist12);
  R236 = gmodify (tlist, grp);
  odisp (&tlist);
  return R236;
}

ocharptr
gpleth (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr tempo = NULL;
  ocharptr chrc;

  echo = false;
  group = grp.name;
  qspecial = false;
  switch ((int) (group))
    {
    case un:
    case sung:
    case spn:
    case on:
    case son:
    case g2:
    case spnc:
    case sonc:
      chrc = gmodify (chrc1, grp); // alloc
    case nill:
    case sn:
    case an:;
      break;
    default:
      Caseerror (Line);
    }
  qspecial = true;
  {
    register groop *W141 = &grp;

    switch ((int) (W141->name))
      {
      case un:
      case sung:
	tempo = unpleth (chrc, chrc2, grp);
	odisp (&chrc);
	break;
      case spn:
      case on:
	tempo = sp_pleth (chrc, chrc2, grp);
	odisp (&chrc);
	break;
      case spnc:
	tempo = spncpleth (chrc1, chrc2, grp);
	break;
      case sonc:
	tempo = soncpleth (chrc1, chrc2, grp);
	break;
      case son:
	if (!(bool) ((grp.rank) & 1))
	  {
	    print ( "Even rank SO(n) not implemented\n");
	    tempo = NULL;
	  }
	else
	  {
	    tempo = sp_pleth (chrc, chrc2, grp);
	    odisp (&chrc);
	  }
	break;
      case g2:
	tempo = g2pleth (chrc, chrc2);
	odisp (&chrc);
	break;
      case f4:
      case e6:
      case e7:
      case e8:
      case unm:
      case sunm:
      case ospnm:
      case sn:
	inform ("not implemented", cr);
	tempo = NULL;
	break;
      case nill:
	tempo = NULL;
	break;
      default:
	Caseerror (Line);
      }
  }
  echo = true;
  return tempo;
}



void
sigmap (frame a, frame b, frame * sigma, int *r)
{
  int k,  x;
  register int p;
  register int i;
  frame aconj, bconj, mray, temp;

  x = len (&a);
  k = a.A[1];
  for (i = 1; i <= maxdim; i++)
    {
      mray.A[i] = 0;
      sigma->A[i] = 0;
      temp.A[i] = 0;
    }
  for (i = 1; i <= k; i++)
    {
      for (p = 1; p <= x; p++)
	{
	  if (a.A[p] == i)
	    mray.A[i] = mray.A[i] + 1;
	}
    }
  aconj = a;
  bconj = b;
  (*r) = 0;
  conjgte (&aconj);
  conjgte (&bconj);
  //h = len (&a);
  for (i = 1; i <= k; i++)
    {
      if ((aconj.A[i] - bconj.A[i] == 1))
	sigma->A[i] = 1;
    }
  for (i = 1; i <= k; i++)
    {
      if (sigma->A[i] > sigma->A[i + 1])
	{
	  (*r) = (*r) + 1;
	  temp.A[(*r)] = mray.A[i];
	}
    }
  (*sigma) = temp;
}

void
writefrme (frame x)
{
  int i, k, m, t;
  register int j;

  t = maxdim;
  while ((x.A[t] == 0) && (t > 0))
    t = t - 1;
  if (!pow_note)
    {
      i = 2;
      print ( "%1d", x.A[1]);
      if (x.A[1] > 9)
	print(" ");
      while ((i <= t))
	{
	  print ( "%1d", x.A[i]);
	  if (abs (x.A[i]) > 9)
	    print ( " ");
	  i = i + 1;
	}
    }
  else
    {
      i = 1;
      m = 1;
      while ((m <= t) && (i < maxdim))
	{
	  j = x.A[i];
	  k = 1;
	  do
	    {
	      i = i + 1;
	      if ((x.A[i] == j))
		k = k + 1;
	    }
	  while (!((x.A[i] != j) || (i == maxdim)));
	  if (k > cutoff)
	    {
	      if (((m > 1) && (abs (x.A[m]) > 9)))
		print ( " ");
	      print ( "%1d^%1d%1c", x.A[m], k, ' ');
	    }
	  else
	    for (j = 1; j <= k; j++)
	      {
		if (abs (x.A[m]) > 9)
		  {
		    print ("%1d%1c", x.A[m], ' ');
		    if (x.A[m] < 0)
		      print ( " ");
		  }
		else
		  print ( "%1d", x.A[m]);
	      }
	  m = i;
	}
      if (i == 1)
	print ( "0");
    }
}

void
hallp (frame a, frame b)
{
  termptr list;
  int r/*, h, wcount*/;
  register int i;
  frame sigma;
  termptr temp;
  char intertotal[500];
  char *s=intertotal;
  int q=0;

  if ((b.A[2] == 0))
    {
      list = outer (a, b);
      //wcount = 0;
      while (list != NULL)
	{
	  progress ();
	  //h = len (&list->val);
	  sigmap (list->val, a, &sigma, &r);
	  s +=sprintf (s, "(");
	  for (i = 1; i <= r; i++)
	    {
	      s +=sprintf (s, "%d", sigma.A[i]);
	    }
	  s += sprintf (s, "){");
	  s += wrtfrme2 (s,list->val);
	  s += sprintf (s, "}");
	  //wcount = wcount + 14 + r + h;
	  if ( (int)(q+strlen(intertotal)) >= (int)(tcol))
	    {
	      print ("\n");
	      q=0;
	    }
	  print ("%s",intertotal);
	  q +=strlen (intertotal);
	  s=intertotal;
	  if (list->next != NULL)
	    s += sprintf (s, " +");
	  temp=list;
	  list = list->next;
	  dispsfn (& temp);
	}
    }
  else
    print ( "The second partition must be restricted to one part only\n");
}

void
kostkamatrix (int n)
{
  termptr tsfn, ksfn, msfn, k, m;
  int msum, qq, ksum;
  bframe bigg;
  qq = 12;
  maxb = n;
  ksfn = (termptr) gensfnlist (n);
  msfn = ksfn;
  tsfn = ksfn;
  msum = 0;
  ksum = 0;
  while (msfn != NULL)
    {
      register termptr W157 = &(*msfn);

      msum = msum + 1;
      ksum = msum;
      snu (&m);
      m->mult = W157->mult;
      m->val = W157->val;
      m->next = NULL;
      while (ksfn != NULL)
	{
	  register termptr W158 = &(*ksfn);

	  if ((ksum >= msum))
	    {
	      snu (&k);
	      k->mult = 1;
	      k->val = W158->val;
	      k->next = NULL;
	      kostka (m, k, &bigg);
	      dispsfn (&k);
	      if (logging)
		{
		  wrtbigno (&logfile, &qq, &bigg);
		  Putchr (' ', logfile);
		}
	      wrtbigno (&output, &qq, &bigg);
	      Putchr (' ', output);
	    }
	  else
	      print ( "%10d", 0);
	  ksum = ksum + 1;
	  ksfn = W158->next;
	}
      ksfn = tsfn;
      print ("\n");
      dispsfn (&m);
      msfn = W157->next;
    }
  ldisp (&tsfn);
  maxb = maxl;
}

void
kostka (termptr a, termptr b, bframe * result)
{
  bool test;
  int n, w, wa, wb, i;
  termptr temp, tempa, dummy;
  ocharptr newlist;

  test = true;
  wa = wtfrm (&a->val);
  wb = wtfrm (&b->val);
  maxb = wa;
  if ((wa != wb) || (testord (&a->val, &b->val) == LESS))
    test = false;
  if ((test == true))
    {
      i = 1;
      snu (&tempa);
      tempa->mult = a->mult;
      tempa->val = a->val;
      w = 0;
      while ((b->val.A[i] >= 2))
	{
	  snu (&temp);
	  temp->mult = 1;
	  temp->val = nolls;
	  temp->val.A[1] = b->val.A[i];
	  temp->next = NULL;
	  dummy = lskew (tempa, temp);
	  ldisp (&tempa);
	  tempa = dummy;
	  w = w + temp->val.A[1];
	  dispsfn (&temp);
	  i = i + 1;
	}
      n = wa - w;
      if (n >= 1)
	{
	  newlist = sfntochrc (tempa, false, ' ');
	  ldisp (&tempa);
	  dimnsymm (n, newlist, &(*result));
	  odisp (&newlist);
	}
      else if (n == 0)
	{
	  if (tempa == NULL)
	    i = 0;
	  else
	    i = tempa->mult;
	  tobig (i, &(*result));
	  if (i != 0)
	    ldisp (&tempa);
	}
    }
  else
    {
      tobig (0, &(*result));

      if ((wa != wb))
	print ( "ERROR: Incompatible partitions?\n");
    }
  maxb = maxl;

}




ocharptr
spncpleth (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr R263;
  ocharptr tempun, temp, unreplist, sprlist, treplist, newlist, blist, tlist;
  int w, spnum, gnum, kk, k, kappa, wrep, pk;
  bool sspin;

  newlist = NULL;
  spnum = grp.rank;
  gnum = spnum / 2;
  kappa = chrc1->val.A[maxdim];
  k = wtfrm (&chrc2->val);
  w = lowestweight (chrc1);
  wrep = k * w;
  if ((k % 2 == 1) && (chrc1->spin == true))
    sspin = true;
  else
    sspin = false;
  kk = k / 2;
  pk = (k * kappa);
  if ((chrc1->spin == true))
    pk = pk + kk;

  treplist = genspnrlist (chrc1, kappa, grp);
  fixxcg (&grp, un, gnum, 0);
  unreplist = unpleth (treplist, chrc2, grp);
  fixxcg (&grp, spnc, spnum, 0);
  while (wrep <= (plwt - pk - 1))
    {
      temp = chrccopy (unreplist);
      sprlist = gensprunlist (unreplist, wrep, pk, sspin, grp);
      tempun = chrccopy (sprlist);
      oadd (&newlist, &sprlist);
      odisp (&sprlist);
      fixxcg (&grp, spnc, spnum, 0);
      plwt = plwt + gnum * pk;
      blist = genspnrlist (tempun, pk, grp);
      plwt = plwt - gnum * pk;
      odisp (&tempun);
      fixxcg (&grp, un, gnum, 0);
      tlist = chrcmult (-1, blist);
      odisp (&blist);
      unreplist = chrcadd (temp, tlist);
      odisp (&temp);
      odisp (&tlist);
      fixxcg (&grp, spnc, spnum, 0);
      wrep = wrep + 1;
    }
  odisp (&unreplist);
  osort (&newlist, true);
  R263 = newlist;
  return R263;
}

ocharptr
genunlist (ocharptr plrep, groop grp, bool zeroone)
{
  register ocharptr R264;
  termptr temp;
  ocharptr clist, newlist;

  temp = seriesx ('m', plwt, full);
  if (zeroone)
    schur_restrict (&temp, 1, 'e');
  else
    schur_restrict (&temp, 1, 'o');
  newlist = sfntochrc (temp, false, ' ');
  ldisp (&temp);
  clist = unpleth (newlist, plrep, grp);
  odisp (&newlist);
  R264 = clist;
  return R264;
}

ocharptr
genspnrlist (ocharptr repp, int kappa, groop grp)
{
  register ocharptr R265;
  ocharptr trep, urep, temp, newlist;
  int gnum;
  register int i;

  newlist = NULL;
  urep = (ocharptr) lsprun (&repp, currgrp.A[1 - 1].rank);
  trep = urep;
  gnum = grp.rank / 2;
  fixxcg (&grp, un, gnum, 0);
  while (urep != NULL)
    {
      register ocharptr W159 = &(*urep);

      cnu (&temp);
      temp->mult = W159->mult;
      temp->val = nolls;
      temp->C6_double = false;
      for (i = 1; i <= gnum; i++)
	{
	  temp->val.A[i] = W159->val.A[i] - kappa;
	}
      temp->spin = false;
      oadd (&newlist, &temp);
      dispchr (&temp);
      urep = W159->next;
    }
  odisp (&trep);
  osort (&newlist, true);
  R265 = newlist;
  return R265;
}

ocharptr
gensprunlist (ocharptr gunlist, int wrep, int kappa, bool sspin, groop grp)
{
  register ocharptr R266;
  ocharptr newlist, tlist, tgun;
  int spnum;

  newlist = NULL;
  spnum = 2 * grp.rank;
  rrestrict (&gunlist, -wrep, 'w');
  tgun = gunlist;
  while (gunlist != NULL)
    {
      register ocharptr W162 = &(*gunlist);

      cnu (&tlist);
      tlist->mult = W162->mult;
      tlist->val = W162->val;
      tlist->val.A[maxdim] = kappa;
      tlist->conval = nolls;
      tlist->C6_double = true;
      tlist->spin = sspin;
      oadd (&newlist, &tlist);
      dispchr (&tlist);
      gunlist = W162->next;
    }
  odisp (&tgun);
  fixxcg (&grp, spnc, spnum, 0);
  R266 = gmodify (newlist, grp);
  return R266;
}

void
bpower (int j, int k, bframe * result)
{
  register int i;

  if ((k < 0))
    print ( "ERROR:negative exponent\n");
  else if ((k == 0))
    factor (1, &(*result));
  else if ((k >= 1))
    factor (j, &(*result));
  if ((k >= 2))
    for (i = 1; i <= k - 1; i++)
      {
	kmult ((*result), j, &(*result));
      }
}

void
hclass (frame x, bframe * tresult)
{
  int i, j, k, m, t, w;
  bframe temp1, temp2, temp;

  w = wtfrm (&x);
  maxb = w;
  factorialn (w, &temp1);
  bigtofact (temp1, &temp1);
  t = len (&x);
  i = 1;
  m = 1;
  factor (1, &(*tresult));
  while ((m <= t) && (i <= t))
    {
      j = x.A[i];
      k = 1;
      if ((j != 0))
	do
	  {
	    i = i + 1;
	    if ((x.A[i] == j))
	      k = k + 1;
	  }
	while (!((x.A[i] != j) || (i > t)));
      factorialn (k, &temp2);
      bigtofact (temp2, &temp2);
      bpower (j, k, &temp);
      fmult ((*tresult), temp, &(*tresult));
      fmult ((*tresult), temp2, &(*tresult));
      m = i;
    }
  fdiv (temp1, (*tresult), &(*tresult));
  factobig ((*tresult), &(*tresult));
  maxb = maxl;
}

termptr
classlist (termptr list)
{
  register termptr R267;
  termptr newlist, temp, tlist;
  bframe btemp;

  newlist = NULL;
  while (list != NULL)
    {
      register termptr W163 = &(*list);

      snu (&tlist);
      tlist->mult = W163->mult;
      tlist->val = W163->val;
      hclass (tlist->val, &btemp);
      tlist->mult = W163->mult * frbig (btemp);
      temp = ladd (newlist, tlist);
      ldisp (&newlist);
      newlist = temp;
      dispsfn (&tlist);
      list = list->next;
    }
  sort (&newlist, true);
  R267 = newlist;
  return R267;
}

termptr
s_to_p (termptr x)
{
  register termptr R268;
  int w;
  termptr xlist, newlist, temp, clist, slist, plist, dlist;

  w = wtfrm (&x->val);
  dlist = (termptr) gensfnlist (w);
  temp = classlist (dlist);
  ldisp (&dlist);
  newlist = NULL;
  clist = temp;
  while (temp != NULL)
    {
      register termptr W164 = &(*temp);

      snu (&slist);
      slist->mult = x->mult;
      slist->val = x->val;
      snu (&plist);
      plist->mult = 1;
      plist->val = W164->val;
      dlist = (termptr) character (slist, plist);
      dispsfn (&slist);
      if (dlist != NULL)
	{
	  plist->mult = dlist->mult * W164->mult;
	  dispsfn (&dlist);
	}
      else
	plist->mult = 0;
      if (plist->mult != 0)
	{
	  xlist = ladd (newlist, plist);
	  ldisp (&newlist);
	  newlist = xlist;
	}
      dispsfn (&plist);
      temp = temp->next;
    }
  ldisp (&clist);
  R268 = newlist;
  return R268;
}

ocharptr
kinsert (ocharptr list, int k, groop grp)
{
  register ocharptr R269;
  ocharptr newlist, tlist;
  int t;

  newlist = NULL;
  if ((grp.name == spnc) || (grp.name == sonc))
    {
      while (list != NULL)
	{
	  register ocharptr W165 = &(*list);

	  cnu (&tlist);
	  tlist->mult = W165->mult;
	  tlist->C6_double = true;
	  tlist->spin = list->spin /*false */ ;
	  t = k;
	  /*if ((bool)((k) & 1) & (grp.name == spnc)) {
	     tlist->spin = true;
	     t = t - 1;
	     } */
	  tlist->val = W165->val;
	  tlist->conval = nolls;
	  tlist->val.A[maxdim] = t;
	  tlist->conval.A[1] = t;
	  oadd (&newlist, &tlist);
	  list = W165->next;
	}
    }
  R269 = newlist;
  return R269;
}

ocharptr
sponmodify (ocharptr list, groop grp)
{
  register ocharptr R270;
  ocharptr newlist, tlist;
  int r, k;

  newlist = NULL;
  list = gmodify (list, grp);
  if ((grp.name == spnc))
    {
      r = grp.rank;
      while (list != NULL)
	{
	  register ocharptr W166 = &(*list);
	  cnu (&tlist);
	  tlist->mult = W166->mult;
	  tlist->spin = false;
	  tlist->val = W166->val;
	  tlist->val.A[maxdim] = 0;
	  tlist->lab = ' ';	/*25/12/97 */
	  k = 2 * W166->val.A[maxdim];
	  if (W166->spin)
	    k = k + 1;
	  group = on;
	  omodify (&tlist, k);
	  tlist->conlab = ' ';
	  tlist->val.A[maxdim] = W166->val.A[maxdim];
	  tlist->spin = W166->spin;
	  tlist->conval = W166->conval;
	  tlist->C6_double = W166->C6_double;
	  group = spnc;
	  grp.rank = r;
	  oadd (&newlist, &tlist);
	  list = W166->next;
	}
      sposort (&newlist, true);
    }
  R270 = newlist;
  return R270;
}

ocharptr
associate (ocharptr list, groop grp)
{
  register ocharptr R271;
  ocharptr newlist, tlist;
  int k = 0;			// modified by FB, was uninitialized

  newlist = NULL;
  if ((grp.name == spnc) || (grp.name == on))
    {
      while (list != NULL)
	{
	  register ocharptr W167 = &(*list);

	  cnu (&tlist);
	  tlist->mult = W167->mult;
	  tlist->spin = W167->spin;
	  tlist->val = W167->val;
	  tlist->C6_double = W167->C6_double;
	  tlist->conval = W167->conval;
	  tlist->C6_double = W167->C6_double;
	  tlist->conlab = ' ';
	  tlist->lab = W167->lab;
	  if ((grp.name == spnc))
	    {
	      k = 2 * W167->val.A[maxdim];
	      if (W167->spin)
		k = k + 1;
	      omodify (&tlist, k);
	    }
	  if ((grp.name == on))
	    {
	      k = grp.rank;
	    }
	  if ((bool) ((k) & 1))
	    {
	      if (W167->lab == '#')
		tlist->lab = ' ';
	      else
		tlist->lab = '#';
	    }
	  if (k % 2 == 0)
	    {
	      if (((len (&W167->val) < k / 2) && (W167->spin == false)))
		{
		  if (W167->lab == '#')
		    tlist->lab = ' ';
		  else
		    tlist->lab = '#';
		}
	      else if (((len (&W167->val) == k / 2) || W167->spin))
		tlist->lab = ' ';
	    }
	  oadd (&newlist, &tlist);
	  list = W167->next;
	}
      osort (&newlist, true);
    }
  else
    print ( "Inappropriate group set\n");
  R271 = newlist;
  return R271;
}

ocharptr
sprextend (ocharptr list, groop grp)
{
  register ocharptr R272;
  ocharptr newlist, tlist;
  int r, k, h;
  register int i;

  newlist = NULL;
  if ((grp.name == spnc) || (grp.name == on))
    {
      while (list != NULL)
	{
	  register ocharptr W168 = &(*list);
	  if ((grp.name == spnc))
	    k = 2 * W168->val.A[maxdim];
	  else
	    k = grp.rank;
	  if (W168->spin)
	    k = k + 1;
	  r = len (&W168->val);
	  h = k - 2 * r;
	  cnu (&tlist);
	  tlist->mult = W168->mult;
	  tlist->val = W168->val;
	  tlist->lab = W168->lab;
	  tlist->conlab = W168->conlab;
	  tlist->spin = W168->spin;
	  tlist->conval = W168->conval;
	  tlist->C6_double = W168->C6_double;
	  if (h > 0)
	    for (i = (r + 1); i <= r + h; i++)
	      {
		tlist->val.A[i] = 1;
	      }
	  oadd (&newlist, &tlist);
	  list = W168->next;
	}
      osort (&newlist, true);
    }
  else
    print ( "sprextend: Inappropriate group set\n");
  R272 = newlist;
  return R272;
}

ocharptr
star (ocharptr list, groop grp)
{
  register ocharptr W169;
  ocharptr newlist, tlist, clist, savelist;
  ocharac oneterm;	// FB 20060323
  int k, r;

  newlist = NULL;
  oneterm.next=NULL;
  tlist=&oneterm;      // FB 20060323
  if ((grp.name == spnc))
    {
      list = gmodify (list, grp); // alloc
      savelist=list;
      while (list != NULL)
	{
	  W169 = list;
	  k = 2 * W169->val.A[maxdim];
	  if (W169->spin)
	    k = k + 1;
	  k = (k / 2);
	  r = len (&W169->val);
	  //cnu (&tlist);		// alloc
	  tlist->mult = W169->mult;
	  tlist->val = W169->val;
	  tlist->lab = W169->lab;
	  tlist->conlab = W169->conlab;
	  tlist->spin = W169->spin;
	  tlist->conval = W169->conval;
	  tlist->C6_double = W169->C6_double;
	  
	  if ((r <= k))
	    clist = sprextend (tlist, grp); //alloc
	  if ((r > k))
	    clist = sponmodify (tlist, grp); //alloc
	  clist->lab = ' ';
	  oadd (&newlist, &clist);
	  //odisp(&tlist);           // FB 20060323
	  list = W169->next;
	}
      osort (&newlist, true);
    }
  W169 = gmodify (newlist, grp); // alloc
  odisp(&newlist);// FB 20060323
  odisp(&savelist);// FB 20060323
  return W169;
}

ocharptr
soncpleth (ocharptr chrc1, ocharptr chrc2, groop grp)
{
  register ocharptr R274;
  ocharptr tempun, temp, unreplist, sprlist, treplist, newlist, blist, tlist;
  int w, spnum, gnum, /*kk,*/ k, kappa, wrep, pk;


  newlist = NULL;
  spnum = grp.rank;
  gnum = spnum / 2;
  kappa = chrc1->val.A[maxdim];
  k = wtfrm (&chrc2->val);
  w = lowestweight (chrc1);
  wrep = k * w;

  //kk = k / 2;
  pk = (k * kappa);
  treplist = gensonrlist (chrc1, kappa, grp);
  fixxcg (&grp, un, gnum, 0);
  unreplist = unpleth (treplist, chrc2, grp);
  fixxcg (&grp, sonc, spnum, 0);
  while (wrep <= (plwt - pk - 1))
    {
      temp = chrccopy (unreplist);
      sprlist = gensorunlist (unreplist, wrep, pk, grp);
      tempun = chrccopy (sprlist);
      oadd (&newlist, &sprlist);
      odisp (&sprlist);
      fixxcg (&grp, sonc, spnum, 0);
      plwt = plwt + gnum * pk;
      blist = gensonrlist (tempun, pk, grp);
      plwt = plwt - gnum * pk;
      odisp (&tempun);
      fixxcg (&grp, un, gnum, 0);
      tlist = chrcmult (-1, blist);
      odisp (&blist);
      unreplist = chrcadd (temp, tlist);
      odisp (&temp);
      odisp (&tlist);
      fixxcg (&grp, sonc, spnum, 0);
      wrep = wrep + 1;
    }
  odisp (&unreplist);
  osort (&newlist, true);
  /*rrestrict(&newlist, setlimit-2, 'w'); */
  R274 = newlist;
  return R274;
}

ocharptr
gensonrlist (ocharptr repp, int kappa, groop grp)
{
  register ocharptr R275;
  ocharptr trep, urep, temp, newlist;
  int gnum;
  register int i;

  newlist = NULL;
  urep = (ocharptr) lsoncbrun (&repp, currgrp.A[1 - 1].rank);
  trep = urep;
  gnum = grp.rank / 2;
  fixxcg (&grp, un, gnum, 0);
  while (urep != NULL)
    {
      register ocharptr W159 = &(*urep);

      cnu (&temp);
      temp->mult = W159->mult;
      temp->val = nolls;
      temp->C6_double = false;
      for (i = 1; i <= gnum; i++)
	{
	  temp->val.A[i] = W159->val.A[i] - kappa;
	}
      temp->spin = false;
      oadd (&newlist, &temp);
      dispchr (&temp);
      urep = W159->next;
    }
  odisp (&trep);
  osort (&newlist, true);
  R275 = newlist;
  return R275;
}

ocharptr
gensorunlist (ocharptr gunlist, int wrep, int kappa, groop grp)
{
  register ocharptr R276;
  ocharptr newlist, tlist, tgun;
  int spnum;

  newlist = NULL;
  spnum = 2 * grp.rank;
  rrestrict (&gunlist, -wrep, 'w');
  tgun = gunlist;
  while (gunlist != NULL)
    {
      register ocharptr W162 = &(*gunlist);

      cnu (&tlist);
      tlist->mult = W162->mult;
      tlist->val = W162->val;
      tlist->val.A[maxdim] = kappa;
      tlist->conval = nolls;
      tlist->C6_double = true;
      tlist->spin = false;
      tlist->conval.A[1] = kappa;
      oadd (&newlist, &tlist);
      dispchr (&tlist);
      gunlist = W162->next;
    }
  odisp (&tgun);
  fixxcg (&grp, sonc, spnum, 0);
  R276 = gmodify (newlist, grp);
  return R276;
}

int
lowestweight (ocharptr list)
{
  register int R277;
  int w, lw;

  w = maxdim;
  while (list != NULL)
    {
      register ocharptr W159 = &(*list);

      lw = wtfrm (&W159->val);
      if (lw < w)
	w = lw;
      list = W159->next;
    }
  R277 = w;
  return R277;
}

ocharptr
g2p (ocharptr list)
{
  register ocharptr R278;
  bool test;
  ocharptr temp, fun, tlist, plist;
  tlist = list;
  if ((ggroup.name == g2))
    {
      g2modify (&tlist, 0);
      if (tlist != NULL)
	{
	  racahg2 (&tlist);
	  test = false;
	  cnu (&fun);
	  fun->mult = 1;
	  fun->val = nolls;
	  fun->val.A[1] = 1;
	  fun->spin = false;
	  fun->C6_double = false;
	  fun->val.length=1;
	  do
	    {
	      plist = gpleth (fun, tlist, currgrp.A[1 - 1]);

	      if (plist == NULL)
		      fprintf(stderr,"plist==NULL 1\n");
	      osort (&plist, false);
	      if (plist == NULL)
		      fprintf(stderr,"plist==NULL 2\n");
	      temp = plist->next;
	      if (temp != NULL)
		{
		  temp->mult = -temp->mult;
		  racahg2 (&temp);
		  oadd (&tlist, &temp);
		  osort (&tlist, false);

		}
	      else
		test = true;
	    }
	  while (!((test == true)));
	  dispchr (&fun);
	}
      else
	print ("improper G(2) irrep entered \n");
    }
  else
    print ( "mistake! group not set as G(2) \n");
  R278 = tlist;
  return R278;
}

termptr
snredpleth (termptr * list)
{
  register termptr R279;
  bool test;
  termptr temp, tlist, plist;

  tlist = (*list);
  putsfn (&output, tlist, true);
  if (logging)
    putsfn (&logfile, tlist, true);
  if ((tlist != NULL))
    {
      test = false;
      do
	{
	  plist = plethonerinner (tlist);
	  sort (&plist, true);
	  putsfn (&output, plist, true);
	  if (logging)
	    putsfn (&logfile, plist, true);

	  temp = plist->next;
	  if (temp != NULL)
	    {
	      temp->next = NULL;
	      temp->mult = -temp->mult;
	      add (&tlist, &temp);
	      sort (&tlist, true);
	    }
	  else
	    test = true;
	}
      while (!((test == true)));
    }
  R279 = tlist;
  return R279;
}
