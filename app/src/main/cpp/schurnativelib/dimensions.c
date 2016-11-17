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

/** \file dimensions.c
 * Computation of dimensions
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
# include "s3.h"
# include "bignums.h"
# include "m.h"
# include "g.h"
# include "dimensions.h"

void
orthdim (int n, frame irrep, bframe * dimn, bool spinn)
{
  int s, p;
  register int i;
  register int j;
  bframe hk, subdimn1, subdimn2;
  frame orig;
  register int lastStep,lastStep2;

  p = n / 2;
  if (spinn)
    s = 1;
  else
    s = 0;
  if ((dimb && (n <= maxdim)))
    {
      factor (1, &subdimn1);
      factor (1, &subdimn2);
      for (j = 1; j <= p; j++)
	{
	  for (i = 1; i <= j; i++)
	    {
	      if (j > i)
		{
		  kmult (subdimn1, irrep.A[i] - irrep.A[j] - i + j,
			 &subdimn1);
		  kmult (subdimn2, -i + j, &subdimn2);
		}
	      if (((j > i) && (!(bool) ((n) & 1)))
		  || ((j >= i) && ((bool) ((n) & 1))))
		{
		  kmult (subdimn1,
			 irrep.A[i] + irrep.A[j] + n - i - j + s, &subdimn1);
		  kmult (subdimn2, n - i - j, &subdimn2);
		}
	    }
	}
      switch ((int) (group))
	{
	case on:
	case son:
	  if (((!(bool) ((n) & 1)) && ((irrep.A[p] > 0) || spinn)))
	    kmult (subdimn1, 2, &subdimn1);
	  break;
	case un:
	case sung:
	case mp:
	case spnc:
	case sonc:
	case unc:
	case spn:
	case g2:
	case f4:
	case e6:
	case e7:
	case e8:
	case sunm:
	case ospnm:
	case unm:
	case nill:
	case sn:
	case an:;
	  break;
	default:
	  Caseerror (Line);
	}
      fdiv (subdimn1, subdimn2, &(*dimn));
    }
  else
    {
      orig = irrep;
      conjgte (&irrep);
      factor (1, &(*dimn));
      {
	lastStep = irrep.A[1];
	for (j = 1; j <= lastStep; j++)
	  {
	    {
	      lastStep2 = orig.A[j];
	      for (i = 1; i <= lastStep2; i++)
		{
		  if ((!spinn))
		    {
		      if (i >= j)
			kmult ((*dimn),
			       n + orig.A[i] + orig.A[j] - i - j + s,
			       &(*dimn));
		      else
			kmult ((*dimn),
			       n - irrep.A[i] - irrep.A[j] + i + j - 2 +
			       s, &(*dimn));
		    }
		  else if (i > j)
		    kmult ((*dimn), n + orig.A[i] + orig.A[j] - i - j + s,
			   &(*dimn));
		  else
		    kmult ((*dimn),
			   n - irrep.A[i] - irrep.A[j] + i + j - 2 + s,
			   &(*dimn));
		}
	    }
	  }
      }
      hklth (irrep, &hk);
      fdiv ((*dimn), hk, &(*dimn));
      if (spinn)
	for (i = 1; i <= p; i++)
	  {
	    kmult ((*dimn), 2, &(*dimn));
	  }
    }
}

void
sympdim (int n, frame irrep, bframe * dimn)
{
  int p;
  register int j;
  register int k;
  register int i;
  register int lastStep,lastStep2;
  bframe hk, subdimn1, subdimn2;
  frame orig;

  if ((dimb && (n <= maxdim)))
    {
      factor (1, &subdimn1);
      factor (1, &subdimn2);
      p = n / 2;
      for (i = 1; i <= p; i++)
	{
	  kmult (subdimn1, irrep.A[i] + p - i + 1, &subdimn1);
	  kmult (subdimn2, p - i + 1, &subdimn2);
	}
      for (k = 2; k <= p; k++)
	{
	  for (i = 1; i <= k - 1; i++)
	    {
	      kmult (subdimn1, irrep.A[i] - irrep.A[k] + k - i, &subdimn1);
	      kmult (subdimn1,
		     irrep.A[i] + irrep.A[k] + n - i - k + 2, &subdimn1);
	      kmult (subdimn2, k - i, &subdimn2);
	      kmult (subdimn2, n + 2 - i - k, &subdimn2);
	    }
	}
      fdiv (subdimn1, subdimn2, &(*dimn));
    }
  else
    {
      orig = irrep;
      conjgte (&irrep);
      factor (1, &(*dimn));
      {
	lastStep = irrep.A[1];
	for (j = 1; j <= lastStep; j++)
	  {
	    {
	      lastStep2 = orig.A[j];
	      for (i = 1; i <= lastStep2; i++)
		{
		  if (i > j)
		    kmult ((*dimn), n + orig.A[i] + orig.A[j] - i - j + 2,
			   &(*dimn));
		  else
		    kmult ((*dimn), n - irrep.A[i] - irrep.A[j] + i + j,
			   &(*dimn));
		}
	    }
	  }
      }
      hklth (irrep, &hk);
      fdiv ((*dimn), hk, &(*dimn));
    }
}

void
dimnunit (int n, ocharptr charlist, bframe * dimn)
{
  register int j;
  register int i;
  register int lastStep,lastStep2;
  frame irrep, orig, pt1, pt2, cpt1, cpt2;
  bframe subdimn, subdimn2, cvdimn, cndimn, hk;

  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr W22 = &(*charlist);

      if (n > 1)
	{
	  if (W22->C6_double)
	    {
	      pt1 = W22->val;
	      cpt1 = pt1;
	      conjgte (&cpt1);
	      pt2 = W22->conval;
	      cpt2 = pt2;
	      conjgte (&cpt2);
	      factor (1, &cvdimn);
		lastStep = cpt1.A[1];
		for (i = 1; i <= lastStep; i++)
		  {
		      lastStep2 = pt1.A[i];
		      for (j = 1; j <= lastStep2; j++)
			{
			  kmult (cvdimn,
				 n - cpt1.A[j] - cpt2.A[i] + i + j - 1,
				 &cvdimn);
			}
		  }
	      factor (1, &cndimn);
		lastStep = cpt2.A[1];
		for (i = 1; i <= lastStep; i++)
		  {
		      lastStep2 = pt2.A[i];
		      for (j = 1; j <= lastStep2; j++)
			{
			  kmult (cndimn,
				 n + pt1.A[j] + pt2.A[i] - i - j + 1,
				 &cndimn);
			}
		  }
	      hklth (pt1, &hk);
	      fdiv (cvdimn, hk, &cvdimn);
	      hklth (pt2, &hk);
	      fdiv (cndimn, hk, &cndimn);
	      fmult (cvdimn, cndimn, &subdimn);
	    }
	  else
	    {
	      factor (1, &subdimn);
	      factor (1, &subdimn2);
	      if ((dimb && (n <= maxdim)))
		{
		  for (j = 2; j <= n; j++)
		    {
		      for (i = 1; i <= j - 1; i++)
			{
			  kmult (subdimn,
				 W22->val.A[i] - W22->val.A[j] - i +
				 j, &subdimn);
			  kmult (subdimn2, j - i, &subdimn2);
			}
		    }
		  fdiv (subdimn, subdimn2, &subdimn);
		}
	      else
		{
		  irrep = W22->val;
		  orig = irrep;
		  conjgte (&orig);
		    lastStep = irrep.A[1];
		    for (j = 1; j <= lastStep; j++)
		      {
			  lastStep2 = orig.A[j];
			  for (i = 1; i <= lastStep2; i++)
			    {
			      kmult (subdimn, n - i + j, &subdimn);
			      kmult (subdimn2,
				     irrep.A[i] + orig.A[j] - i - j + 1,
				     &subdimn2);
			    }
		      }
		  fdiv (subdimn, subdimn2, &subdimn);
		}
	    }
	}
      else if (n == 1)
	factor (1, &subdimn);
      if (W22->mult != 1)
	kmult (subdimn, W22->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W22->next;
    }
}

void
dimnsymp (int n, ocharptr charlist, bframe * dimn)
{
  bframe subdimn;

  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr W39 = &(*charlist);

      sympdim (n, W39->val, &subdimn);
      if (W39->mult != 1)
	kmult (subdimn, W39->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W39->next;
    }
}

void
dimnso3 (ocharptr charlist, bframe * dimn)
{
  int j;
  bframe subdim;

  tobig (0, &(*dimn));
  tobig (0, &subdim);
  while (charlist != NULL)
    {
      register ocharptr W40 = &(*charlist);

      if (W40->spin)
	j = 2 * W40->val.A[1] + 2;
      else
	j = 2 * W40->val.A[1] + 1;
      tobig (W40->mult * j, &subdim);
      bigadd ((*dimn), subdim, &(*dimn));
      charlist = W40->next;
    }
}
void
dimnorth (int n, ocharptr charlist, bframe * dimn)
{
  bframe subdimn;

  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr W41 = &(*charlist);

      orthdim (n, W41->val, &subdimn, W41->spin);
      if (W41->mult != 1)
	kmult (subdimn, W41->mult, &subdimn);
      if (((W41->lab == '+') || (W41->lab == '-')))
	kdiv (subdimn, 2, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W41->next;
    }
}


void
dimnalt (int n, ocharptr charlist, bframe * dimn)
{
  int j, k, m, pk;
  register int q;
  register int p;
  register int i;
  register int lastStep;
  bframe hk, subdimn;

  anselect (&charlist, n);
  pk = 1;
  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr W42 = &(*charlist);

      factor (1, &subdimn);
      if (W42->spin)
	{
	  m = wtfrm (&charlist->val);
	  k = len (&charlist->val) + 1;
	  j = (n - k) / 2;
	  for (i = 1; i <= j; i++)
	    {
	      kmult (subdimn, 2, &subdimn);
	    }
	  for (i = 0; i <= m - 1; i++)
	    {
	      kmult (subdimn, (n - i), &subdimn);
	    }
	  for (p = 1; p <= k - 1; p++)
	    {
	      for (q = p + 1; q <= k; q++)
		{
		  kmult (subdimn,
			 (charlist->val.A[p] - charlist->val.A[q]), &subdimn);
		}
	    }
	  for (p = 1; p <= k; p++)
	    {
	      kmult (subdimn, (n - m - charlist->val.A[p]), &subdimn);
	    }
	  for (i = 1; i <= k; i++)
	    {
		lastStep = (charlist->val.A[i]);
		for (p = 1; p <= lastStep; p++)
		  {
		    kdiv (subdimn, p, &subdimn);
		  }
	    }
	  for (p = 1; p <= k; p++)
	    {
	      kdiv (subdimn, (n - m + charlist->val.A[p]), &subdimn);
	    }
	  for (p = 1; p <= k - 1; p++)
	    {
	      for (q = p + 1; q <= k; q++)
		{
		  kdiv (subdimn,
			(charlist->val.A[p] + charlist->val.A[q]), &subdimn);
		}
	    }
	  if (W42->lab != ' ')
	    kdiv (subdimn, 2, &subdimn);
	}
      else
	{
	  if (wtfrm (&W42->val) != n)
	    error (BAD_IRREP, pk);
	  else
	    {
	      for (i = 2; i <= n; i++)
		{
		  kmult (subdimn, i, &subdimn);
		}
	      hklth (W42->val, &hk);
	      fdiv (subdimn, hk, &subdimn);
	      if (W42->lab != ' ')
		kdiv (subdimn, 2, &subdimn);
	    }
	}
      if (W42->mult != 1)
	kmult (subdimn, W42->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W42->next;
    }
}

void
dimnsymm (int n, ocharptr charlist, bframe * dimn)
{
  int j, k, m, pk;
  register int q;
  register int p;
  register int i;
  register int lastStep;
  bframe hk, subdimn;

  snselect (&charlist, n, qspecial);
  pk = 1;
  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr tmp_ptr = charlist;

      factor (1, &subdimn);  // will use subdimn as a bignum based on prime factors
      if (tmp_ptr->spin)
	{
	  m = wtfrm (&charlist->val);
	  k = len (&charlist->val);
	  j = (n - k) / 2;
	  for (i = 1; i <= j; i++)
	    {
	      kmult (subdimn, 2, &subdimn);
	    }
	  for (i = 0; i <= m - 1; i++)
	    {
	      kmult (subdimn, (n - i), &subdimn);
	    }
	  for (p = 1; p <= k - 1; p++)
	    {
	      for (q = p + 1; q <= k; q++)
		{
		  kmult (subdimn,
			 (charlist->val.A[p] - charlist->val.A[q]), &subdimn);
		}
	    }
	  for (p = 1; p <= k; p++)
	    {
	      kmult (subdimn, (n - m - charlist->val.A[p]), &subdimn);
	    }
	  for (i = 1; i <= k; i++)
	    {
		lastStep = (charlist->val.A[i]);
		for (p = 1; p <= lastStep; p++)
		  {
		    kdiv (subdimn, p, &subdimn);
		  }
	    }
	  for (p = 1; p <= k; p++)
	    {
	      kdiv (subdimn, (n - m + charlist->val.A[p]), &subdimn);
	    }
	  for (p = 1; p <= k - 1; p++)
	    {
	      for (q = p + 1; q <= k; q++)
		{
		  kdiv (subdimn,
			(charlist->val.A[p] + charlist->val.A[q]), &subdimn);
		}
	    }
	  if (tmp_ptr->lab != ' ')
	    kdiv (subdimn, 2, &subdimn);
	}
      else		// no spin
	{
	  if (wtfrm (&tmp_ptr->val) != n)
	    error (BAD_IRREP, pk);
	  else
	    {
	      for (i = 2; i <= n; i++)
		{
		  kmult (subdimn, i, &subdimn);
		}
	      hklth (tmp_ptr->val, &hk);
	      fdiv (subdimn, hk, &subdimn);
	    }
	}
      if (tmp_ptr->mult != 1)
	kmult (subdimn, tmp_ptr->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = tmp_ptr->next;
    }
}

void
dimne6 (ocharptr charlist, bframe * dimn)
{
  bframe subdimn, divi;
  int wt, s;
  register int r;
  register int q;
  register int p;
  register int i;
  register ocharptr W90;

  tobig (0, &(*dimn));
  for (i = 6; i <= maxl; i++)
    {
      divi.A[i] = 0;
    }
  divi.A[0] = 1;
  divi.A[1] = 45;
  divi.A[2] = 10;
  divi.A[3] = 5;
  divi.A[4] = 3;
  divi.A[5] = 1;
  while (charlist != NULL)
    {
      W90 = &(*charlist);

      wt = 0;
      for (i = 1; i <= 6; i++)
	{
	  wt = W90->val.A[i] + wt;
	}
      if ((wt % 2 == 0))
	{
	  W90->val.A[0] = -1;
	  s = W90->val.A[1];
	  for (i = 2; i <= 6; i++)
	    {
	      s = s - W90->val.A[i];
	    }
	  factor (W90->val.A[1] + 11, &subdimn);
	  for (p = 2; p <= 6; p++)
	    {
	      for (q = p + 1; q <= 7; q++)
		{
		  kmult (subdimn,
			 W90->val.A[p] - W90->val.A[q] + q - p, &subdimn);
		  for (r = q + 1; r <= 7; r++)
		    {
		      kmult (subdimn,
			     2 * W90->val.A[p] +
			     2 * W90->val.A[q] +
			     2 * W90->val.A[r] + s + 38 - 2 * (p +
							       q +
							       r), &subdimn);
		    }
		}
	    }
	  fdiv (subdimn, divi, &subdimn);
	  if (W90->mult != 1)
	    kmult (subdimn, W90->mult, &subdimn);
	  factobig (subdimn, &subdimn);
	  bigadd (subdimn, (*dimn), &(*dimn));
	}
      charlist = W90->next;
    }
}

void
dimne7 (ocharptr charlist, bframe * dimn)
{
  bframe subdimn, divi;
  int wt;
  register int s;
  register int r;
  register int q;
  register int p;
  register int t;

  tobig (0, &(*dimn));
  for (t = 8; t <= maxl; t++)
    {
      divi.A[t] = 0;
    }
  divi.A[0] = 1;
  divi.A[1] = 47;
  divi.A[2] = 22;
  divi.A[3] = 10;
  divi.A[4] = 6;
  divi.A[5] = 3;
  divi.A[6] = 2;
  divi.A[7] = 1;
  while (charlist != NULL)
    {
      register ocharptr W103 = &(*charlist);

      wt = 0;
      for (t = 1; t <= 7; t++)
	{
	  wt = wt + W103->val.A[t];
	}
      if ((wt % 2 == 0))
	{
	  factor (1, &subdimn);
	  for (p = 2; p <= 8; p++)
	    {
	      kmult (subdimn, W103->val.A[1] - W103->val.A[p] + 9 + p,
		     &subdimn);
	      for (q = p + 1; q <= 8; q++)
		{
		  kmult (subdimn,
			 W103->val.A[p] - W103->val.A[q] + q - p, &subdimn);
		  for (r = q + 1; r <= 7; r++)
		    {
		      for (s = r + 1; s <= 8; s++)
			{
			  kmult (subdimn,
				 (wt / 2) - W103->val.A[p] -
				 W103->val.A[q] -
				 W103->val.A[r] -
				 W103->val.A[s] + p + q + r +
				 s - 13, &subdimn);
			}
		    }
		}
	    }
	  fdiv (subdimn, divi, &subdimn);
	  if (W103->mult != 1)
	    kmult (subdimn, W103->mult, &subdimn);
	  factobig (subdimn, &subdimn);
	  bigadd (subdimn, (*dimn), &(*dimn));
	}
      charlist = W103->next;
    }
}

void
dimne8 (ocharptr charlist, bframe * dimn)
{
  bframe subdimn, divi;
  int wt;
  register int r;
  register int q;
  register int p;
  register int t;

  tobig (0, &(*dimn));
  for (t = 11; t <= maxl; t++)
    {
      divi.A[t] = 0;
    }
  divi.A[0] = 1;
  divi.A[1] = 97;
  divi.A[2] = 47;
  divi.A[3] = 21;
  divi.A[4] = 14;
  divi.A[5] = 8;
  divi.A[6] = 6;
  divi.A[7] = 4;
  divi.A[8] = 3;
  divi.A[9] = 2;
  divi.A[10] = 1;
  while (charlist != NULL)
    {
      register ocharptr W116 = &(*charlist);

      factor (1, &subdimn);
      wt = 0;
      for (t = 1; t <= 8; t++)
	{
	  wt = W116->val.A[t] + wt;
	}
      if ((wt % 3 == 0))
	{
	  for (p = 2; p <= 9; p++)
	    {
	      kmult (subdimn, W116->val.A[1] - W116->val.A[p] + 20 + p,
		     &subdimn);
	      for (q = p + 1; q <= 9; q++)
		{
		  kmult (subdimn,
			 W116->val.A[p] - W116->val.A[q] + q - p, &subdimn);
		  kmult (subdimn,
			 W116->val.A[1] + W116->val.A[p] +
			 W116->val.A[q] - (wt / 3) + 28 - p - q, &subdimn);
		  for (r = q + 1; r <= 9; r++)
		    {
		      kmult (subdimn,
			     (wt / 3) - W116->val.A[p] -
			     W116->val.A[q] - W116->val.A[r] + p +
			     q + r - 8, &subdimn);
		    }
		}
	    }
	  fdiv (subdimn, divi, &subdimn);
	  if (W116->mult != 1)
	    kmult (subdimn, W116->mult, &subdimn);
	  factobig (subdimn, &subdimn);
	  bigadd (subdimn, (*dimn), &(*dimn));
	}
      charlist = W116->next;
    }
}

void
dimng2 (ocharptr charlist, bframe * dimn)
{
  bframe subdimn;

  tobig (0, &(*dimn));
  while (charlist != NULL)
    {
      register ocharptr W125 = &(*charlist);

      factor (1, &subdimn);
      kmult (subdimn,
	     (W125->val.A[1] + W125->val.A[2] + 4) * (W125->val.A[1] -
						      2 * W125->val.A[2] + 1),
	     &subdimn);
      kmult (subdimn,
	     (2 * W125->val.A[1] - W125->val.A[2] + 5) * (W125->val.A[1] +
							  3) *
	     (W125->val.A[2] + 1), &subdimn);
      kmult (subdimn, W125->val.A[1] - W125->val.A[2] + 2, &subdimn);
      kdiv (subdimn, 120, &subdimn);
      if (W125->mult != 1)
	kmult (subdimn, W125->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W125->next;
    }
}

void
dimnf4 (ocharptr charlist, bframe * dimn)
{
  bframe subdimn, divi;
  int x;
  register int i;

  tobig (0, &(*dimn));
  for (i = 0; i <= maxl; i++)
    {
      divi.A[i] = 0;
    }
  divi.A[0] = 1;
  divi.A[1] = 15;
  divi.A[2] = 7;
  divi.A[3] = 4;
  divi.A[4] = 2;
  divi.A[5] = 1;
  while (charlist != NULL)
    {
      register ocharptr W128 = &(*charlist);

      if (W128->spin)
	x = 1;
      else
	x = 0;
      factor (1, &subdimn);
      kmult (subdimn,
	     (11 + x + 2 * W128->val.A[1]) * (5 + x + 2 * W128->val.A[2]),
	     &subdimn);
      kmult (subdimn,
	     (3 + x + 2 * W128->val.A[3]) * (1 + x + 2 * W128->val.A[4]),
	     &subdimn);
      kmult (subdimn,
	     (8 + x + W128->val.A[1] + W128->val.A[2]) * (7 + x +
							  W128->val.A[1] +
							  W128->val.A[3]),
	     &subdimn);
      kmult (subdimn,
	     (6 + x + W128->val.A[1] + W128->val.A[4]) * (4 + x +
							  W128->val.A[2] +
							  W128->val.A[3]),
	     &subdimn);
      kmult (subdimn,
	     (3 + x + W128->val.A[2] + W128->val.A[4]) * (2 + x +
							  W128->val.A[3] +
							  W128->val.A[4]),
	     &subdimn);
      kmult (subdimn,
	     (3 + W128->val.A[1] - W128->val.A[2]) * (4 + W128->val.A[1] -
						      W128->val.A[3]),
	     &subdimn);
      kmult (subdimn,
	     (5 + W128->val.A[1] - W128->val.A[4]) * (1 + W128->val.A[2] -
						      W128->val.A[3]),
	     &subdimn);
      kmult (subdimn,
	     (2 + W128->val.A[2] - W128->val.A[4]) * (1 + W128->val.A[3] -
						      W128->val.A[4]),
	     &subdimn);
      kmult (subdimn,
	     10 + 2 * x + W128->val.A[1] + W128->val.A[2] + W128->val.A[3] +
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     9 + x + W128->val.A[1] + W128->val.A[2] + W128->val.A[3] -
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     7 + x + W128->val.A[1] + W128->val.A[2] - W128->val.A[3] +
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     6 + W128->val.A[1] + W128->val.A[2] - W128->val.A[3] -
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     5 + x + W128->val.A[1] - W128->val.A[2] + W128->val.A[3] +
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     4 + W128->val.A[1] - W128->val.A[2] + W128->val.A[3] -
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     2 + W128->val.A[1] - W128->val.A[2] - W128->val.A[3] +
	     W128->val.A[4], &subdimn);
      kmult (subdimn,
	     1 - x + W128->val.A[1] - W128->val.A[2] - W128->val.A[3] -
	     W128->val.A[4], &subdimn);
      fdiv (subdimn, divi, &subdimn);
      if (W128->mult != 1)
	kmult (subdimn, W128->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, (*dimn), &(*dimn));
      charlist = W128->next;
    }
}

void
casmunit (int n, ocharptr list, bframe * cas)
{
  frame lam;
  int wt;
  register int i;

  {
    register ocharptr W129 = &(*list);

    wt = wtfrm (&W129->val);
    lam = W129->val;
    if (W129->C6_double)
      {
	i = 1;
	while (W129->conval.A[i] != 0)
	  {
	    lam.A[n + 1 - i] = -W129->conval.A[i];
	    i = i + 1;
	  }
	wt = wt - wtfrm (&W129->conval);
      }
    tobig (0, &(*cas));
    for (i = 1; i <= n; i++)
      {
	cadd ((*cas), n * lam.A[i] * (lam.A[i] - 2 * i + n + 1), &(*cas));
      }
    cadd ((*cas), -wt * wt, &(*cas));
  }
}

void
casmsymp (int n, ocharptr list, bframe * cas)
{
  register int i;

  {
    register ocharptr W132 = &(*list);

    tobig (0, &(*cas));
    for (i = 1; i <= n / 2; i++)
      {
	cadd ((*cas), W132->val.A[i] * (W132->val.A[i] - 2 * i + n + 2),
	      &(*cas));
      }
  }
}

void
casmorth (int n, ocharptr list, bframe * cas)
{
  int extra;
  register int i;

  {
    register ocharptr W135 = &(*list);

    if (W135->spin)
      extra = 1;
    else
      extra = 0;
    if ((bool) ((n) & 1))
      tobig (extra * n * (n / 2), &(*cas));
    else
      tobig (extra * (n - 1) * (n / 2), &(*cas));
    for (i = 1; i <= n / 2; i++)
      {
	cadd ((*cas),
	      4 * W135->val.A[i] * (W135->val.A[i] - 2 * i + n + extra),
	      &(*cas));
      }
  }
}

void
casmg2 (ocharptr list, bframe * cas)
{
  {
    register ocharptr W138 = &(*list);

    tobig ((W138->val.A[1] - W138->val.A[2]) * (W138->val.A[1] -
						W138->val.A[2]), &(*cas));
    cadd ((*cas),
	  W138->val.A[1] * W138->val.A[2] + 5 * W138->val.A[1] -
	  W138->val.A[2], &(*cas));
  }
}

void
casmf4 (ocharptr list, bframe * cas)
{
  int extra;

  {
    register ocharptr W139 = &(*list);

    if (W139->spin)
      extra = 1;
    else
      extra = 0;
    tobig (W139->val.A[1] * (W139->val.A[1] + 11 + extra), &(*cas));
    cadd ((*cas), W139->val.A[2] * (W139->val.A[2] + 5 + extra), &(*cas));
    cadd ((*cas), W139->val.A[3] * (W139->val.A[3] + 3 + extra), &(*cas));
    cadd ((*cas), W139->val.A[4] * (W139->val.A[4] + 1 + extra) + 11 * extra,
	  &(*cas));
  }
}

void
casme6 (ocharptr list, bframe * cas)
{
  int wt;
  register int i;

  {
    register ocharptr W140 = &(*list);

    tobig (3 * W140->val.A[1] * (W140->val.A[1] + 22), &(*cas));
    for (i = 2; i <= 6; i++)
      {
	cadd ((*cas), 6 * W140->val.A[i] * (W140->val.A[i] - 2 * i + 9),
	      &(*cas));
      }
    wt =
      W140->val.A[2] + W140->val.A[3] + W140->val.A[4] + W140->val.A[5] +
      W140->val.A[6];
    cadd ((*cas), -wt * wt, &(*cas));
  }
}

void
casme7 (ocharptr list, bframe * cas)
{
  int wt;
  register int i;

  {
    register ocharptr W143 = &(*list);

    wt = wtfrm (&W143->val) / 2;
    tobig (40 * W143->val.A[1], &(*cas));
    for (i = 1; i <= 7; i++)
      {
	cadd ((*cas), 2 * W143->val.A[i] * (W143->val.A[i] - 2 * i + 6),
	      &(*cas));
      }
    cadd ((*cas), wt * (2 - wt), &(*cas));
  }
}

void
casme8 (ocharptr list, bframe * cas)
{
  int wt;
  register int i;

  {
    register ocharptr W146 = &(*list);

    wt = wtfrm (&W146->val) / 3;
    tobig (42 * W146->val.A[1], &(*cas));
    for (i = 1; i <= 8; i++)
      {
	cadd ((*cas), W146->val.A[i] * (W146->val.A[i] - 2 * i + 5), &(*cas));
      }
    cadd ((*cas), wt * (1 - wt), &(*cas));
  }
}

void
dimngrp (int n, ocharptr charlist, bframe * dimn)
{
  register int i;
  for (i = 1; i <= maxl; i++)
    {
      dimn->A[i] = 0;
    }
  if ((charlist != NULL))
    {
      switch ((int) (group))
	{
	case un:
	case sung:
	  dimnunit (n, charlist, &(*dimn));
	  break;
	case spn:
	  dimnsymp (n, charlist, &(*dimn));
	  break;
	case son:
	  if (n != 3)
	    dimnorth (n, charlist, &(*dimn));
	  else
	    dimnso3 (charlist, &(*dimn));
	  break;
	case on:
	  dimnorth (n, charlist, &(*dimn));
	  break;
	case sn:
	  dimnsymm (n, charlist, &(*dimn));
	  break;
	case an:
	  dimnalt (n, charlist, &(*dimn));
	  break;
	case g2:
	  dimng2 (charlist, &(*dimn));
	  break;
	case f4:
	  dimnf4 (charlist, &(*dimn));
	  break;
	case e6:
	  dimne6 (charlist, &(*dimn));
	  break;
	case e7:
	  dimne7 (charlist, &(*dimn));
	  break;
	case e8:
	  dimne8 (charlist, &(*dimn));
	  break;
	case l168:
	  diml168 (charlist, &(*dimn));
	  break;
	case nill:
	case unm:
	case ospnm:
	case sunm:
	case spnc:
	case mp:
	  print ( "Inappropriate group \n");
	  break;
	default:
	  Caseerror (Line);
	}
    }
}
void
casmgrp (int n, ocharptr charlist, bframe * casm)
{
  switch ((int) (group))
    {
    /* case un: */ case sung:
      casmunit (n, charlist, &(*casm));
      break;
    case spn:
      casmsymp (n, charlist, &(*casm));
      break;
      /* case spnc:
         casmnsymp(n, charlist, &(*casm));
         break ; */
    /* case on: */ case son:
      casmorth (n, charlist, &(*casm));
      break;
    case g2:
      casmg2 (charlist, &(*casm));
      break;
    case f4:
      casmf4 (charlist, &(*casm));
      break;
    case e6:
      casme6 (charlist, &(*casm));
      break;
    case e7:
      casme7 (charlist, &(*casm));
      break;
    case e8:
      casme8 (charlist, &(*casm));
      break;
    case an:
    case sn:
    case unm:
    case sunm:
    case un:
    case on:
    case mp:
    case ospnm:
    case spnc:
    case unc:
    case nill:
      print ( "Inappropriate group \n");
      break;
    default:
      Caseerror (Line);
    }
}

void
d2index (int n, ocharptr list, bframe * dind)
{
  register int i;
  bframe bigcas, bigdim, bigsub, fcas, fdim, fsub, Base;
  ocharptr nextone;
  for (i = 0; i <= maxl; i++)
    {
      dind->A[i] = 0;
    }
  switch ((int) (group))
    {
    case sung:
    case un:
      factor (n * (n * n - 1), &Base);
      break;
    case son:
    case on:
      if (n == 2)
	factor (n * (n - 1) * 2, &Base);
      else
	factor (n * (n - 1) * 4, &Base);
      break;
    case spn:
      factor (n * (n + 1), &Base);
      break;
    case g2:
      factor (42, &Base);
      break;
    case f4:
      factor (312, &Base);
      break;
    case e6:
      factor (2808, &Base);
      break;
    case e7:
      factor (3192, &Base);
      break;
    case e8:
      factor (14880, &Base);
      break;
    case an:
    case sn:
    case spnc:
    case sonc:
    case unc:
    case mp:
    case unm:
    case ospnm:
    case sunm:
    case nill:;
      break;
    default:
      Caseerror (Line);
    }
  switch ((int) (group))
    {
    case sung:
    case son:
    case spn:			/*case on: *//*1/1/98 */
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
      while (list != NULL)
	{
	  register ocharptr W153 = &(*list);

	  nextone = W153->next;
	  W153->next = NULL;
	  casmgrp (n, list, &bigcas);
	  dimngrp (n, list, &bigdim);
	  W153->next = nextone;
	  list = nextone;
	  bigtofact (bigcas, &fcas);
	  bigtofact (bigdim, &fdim);
	  fmult (fcas, fdim, &fsub);
	  fdiv (fsub, Base, &fsub);
	  factobig (fsub, &bigsub);
	  bigadd (bigsub, (*dind), &(*dind));
	}
      break;
    case un:
    case on:
    case en:
    case sn:
    case nill:
    case an:
    case mp:
    case spnc:
    case unm:
    case sunm:
    case ospnm:;
      break;
    default:
      Caseerror (Line);
    }
}


/*termptr dynklbl();*/

void
todynk (ocharptr lambda, frame * a, int n)
{
  int k;
  register int i;

  k = n / 2;
  for (i = 1; i <= maxdim; i++)
    {
      a->A[i] = 0;
    }
  {
    register ocharptr W156 = &(*lambda);

    switch ((int) (group))
      {
      case un:
      case on:
      case en:
      case sn:
      case nill:
      case an:
      case mp:
      case spnc:
      case unm:
      case sunm:
      case ospnm:;
	break;
      case sung:
	for (i = 1; i <= n - 1; i++)
	  {
	    a->A[i] = W156->val.A[i] - W156->val.A[i + 1];
	  }
	break;
      case son:
	if ((bool) ((n) & 1))
	  {
	    for (i = 1; i <= k - 1; i++)
	      {
		a->A[i] = W156->val.A[i] - W156->val.A[i + 1];
	      }
	    if (W156->spin)
	      a->A[k] = 2 * W156->val.A[k] + 1;
	    else
	      a->A[k] = 2 * W156->val.A[k];
	  }
	else
	  {
	    for (i = 1; i <= k - 2; i++)
	      {
		a->A[i] = W156->val.A[i] - W156->val.A[i + 1];
	      }
	    if (W156->lab == '+')
	      if (W156->spin)
		{
		  a->A[k] = W156->val.A[k - 1] - W156->val.A[k];
		  a->A[k - 1] = W156->val.A[k - 1] + W156->val.A[k] + 1;
		}
	      else
		{
		  a->A[k] = W156->val.A[k - 1] - W156->val.A[k];
		  a->A[k - 1] = W156->val.A[k - 1] + W156->val.A[k];
		}
	    else if (W156->spin)
	      {
		a->A[k] = W156->val.A[k - 1] + W156->val.A[k] + 1;
		a->A[k - 1] = W156->val.A[k - 1] - W156->val.A[k];
	      }
	    else
	      {
		a->A[k] = W156->val.A[k - 1] + W156->val.A[k];
		a->A[k - 1] = W156->val.A[k - 1] - W156->val.A[k];
	      }
	  }
	break;
      case spn:
	for (i = 1; i <= k; i++)
	  {
	    a->A[i] = W156->val.A[i] - W156->val.A[i + 1];
	  }
	break;
      case g2:
	a->A[1] = W156->val.A[2];
	a->A[2] = W156->val.A[1] - 2 * W156->val.A[2];
	break;
      case f4:
	a->A[1] = W156->val.A[2] - W156->val.A[3];
	a->A[2] = W156->val.A[3] - W156->val.A[4];
	if (W156->spin)
	  {
	    a->A[3] = 2 * W156->val.A[4] + 1;
	    a->A[4] =
	      W156->val.A[1] - W156->val.A[2] - W156->val.A[3] -
	      W156->val.A[4] - 1;
	  }
	else
	  {
	    a->A[3] = 2 * W156->val.A[4];
	    a->A[4] =
	      W156->val.A[1] - W156->val.A[2] - W156->val.A[3] -
	      W156->val.A[4];
	  }
	break;
      case e6:
	for (i = 1; i <= 5; i++)
	  {
	    a->A[i] = W156->val.A[i + 1] - W156->val.A[i + 2];
	  }
	a->A[6] =
	  (W156->val.A[1] - W156->val.A[2] - W156->val.A[3] - W156->val.A[4] +
	   W156->val.A[5] + W156->val.A[6]) / 2;
	break;
      case e7:
	for (i = 1; i <= 6; i++)
	  {
	    a->A[i] = W156->val.A[8 - i] - W156->val.A[9 - i];
	  }
	a->A[7] =
	  (W156->val.A[1] - W156->val.A[2] - W156->val.A[3] - W156->val.A[4] -
	   W156->val.A[5] + W156->val.A[6] + W156->val.A[7]) / 2;
	break;
      case e8:
	for (i = 1; i <= 7; i++)
	  {
	    a->A[i] = W156->val.A[9 - i] - W156->val.A[10 - i];
	  }
	a->A[8] =
	  (W156->val.A[1] - 2 * W156->val.A[2] - 2 * W156->val.A[3] -
	   2 * W156->val.A[4] + W156->val.A[5] + W156->val.A[6] +
	   W156->val.A[7] + W156->val.A[8]) / 3;
	break;
      default:
	Caseerror (Line);
      }
  }
}

termptr
dynklbl (int n, ocharptr chrc)
{
  register termptr R145;
  termptr dum, list, sub;

  list = NULL;
  sub = NULL;
  switch ((int) (group))
    {
    case an:
    case sn:
    case un:
    case on:
    case unm:
    case sunm:
    case ospnm:
    case mp:
    case spnc:
    case unc:
    case sonc:
      print ( "Inappropriate group\n");
      break;
    case son:
    case sung:
    case spn:
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
      while (chrc != NULL)
	{
	  snu (&sub);
	  {
	    register termptr W171 = &(*sub);

	    todynk (chrc, &W171->val, n);
	    W171->mult = chrc->mult;
	    W171->next = list;
	    list = sub;
	  }
	  chrc = chrc->next;
	}
      sub = NULL;
      while (list != NULL)
	{
	  dum = list;
	  list = list->next;
	  dum->next = sub;
	  sub = dum;
	}
      break;
    default:
      Caseerror (Line);
    }
  R145 = sub;
  return R145;
}


/*ocharptr partlbl();*/

void
frdynk (ocharptr * lambda, frame a, int n)
{
  int k;
  register int i;

  k = n / 2;
  if ((*lambda) == NULL)
    cnu (&(*lambda));
  {
    register ocharptr W172 = &(*(*lambda));
    for (i = 1; i <= maxdim; i++)
      {
	W172->val.A[i] = 0;
      }
    W172->next = NULL;
    W172->spin = false;
    W172->lab = ' ';
    W172->C6_double = false;
    switch ((int) (group))
      {
      case un:
      case on:
      case sn:
      case en:
      case nill:
      case unc:
      case spnc:
      case sonc:
      case an:
      case unm:
      case sunm:
      case ospnm:;
	break;
      case sung:

	for (i = n - 1; i >= 1; i--)
	  {
	    W172->val.A[i] = W172->val.A[i + 1] + a.A[i];
	  }
	break;
      case son:
	if ((bool) ((n) & 1))
	  {
	    if ((bool) ((a.A[k]) & 1))
	      W172->spin = true;
	    W172->val.A[k] = a.A[k] / 2;

	    for (i = k - 1; i >= 1; i--)
	      {
		W172->val.A[i] = W172->val.A[i + 1] + a.A[i];
	      }
	  }
	else
	  {
	    if ((bool) ((a.A[k] + a.A[k - 1]) & 1))
	      W172->spin = true;
	    if (a.A[k - 1] > a.A[k])
	      {
		W172->lab = '+';
		W172->val.A[k] = (a.A[k - 1] - a.A[k]) / 2;
		W172->val.A[k - 1] = (a.A[k - 1] + a.A[k]) / 2;
	      }
	    else
	      {
		W172->val.A[k] = (a.A[k] - a.A[k - 1]) / 2;
		W172->val.A[k - 1] = (a.A[k] + a.A[k - 1]) / 2;
		W172->lab = '-';
	      }
	    if (!W172->spin && (W172->val.A[k] == 0))
	      W172->lab = ' ';

	    for (i = k - 2; i >= 1; i--)
	      {
		W172->val.A[i] = W172->val.A[i + 1] + a.A[i];
	      }
	  }
	break;
      case spn:

	for (i = k; i >= 1; i--)
	  {
	    W172->val.A[i] = W172->val.A[i + 1] + a.A[i];
	  }
	break;
      case g2:
	W172->val.A[1] = 2 * a.A[1] + a.A[2];
	W172->val.A[2] = a.A[1];
	break;
      case f4:
	if ((bool) ((a.A[3]) & 1))
	  W172->spin = true;
	W172->val.A[4] = a.A[3] / 2;
	W172->val.A[3] = W172->val.A[4] + a.A[2];
	W172->val.A[2] = W172->val.A[3] + a.A[1];
	W172->val.A[1] = a.A[1] + 2 * a.A[2] + 3 * a.A[3] / 2 + a.A[4];
	break;
      case e6:
	W172->val.A[6] = a.A[5];

	for (i = 5; i >= 2; i--)
	  {
	    W172->val.A[i] = W172->val.A[i + 1] + a.A[i - 1];
	  }
	W172->val.A[1] =
	  a.A[1] + 2 * a.A[2] + 3 * a.A[3] + 2 * a.A[4] + a.A[5] + 2 * a.A[6];
	break;
      case e7:
	W172->val.A[7] = a.A[1];

	for (i = 6; i >= 2; i--)
	  {
	    W172->val.A[i] = W172->val.A[i + 1] + a.A[8 - i];
	  }
	W172->val.A[1] =
	  2 * a.A[1] + 3 * a.A[2] + 4 * a.A[3] + 3 * a.A[4] + 2 * a.A[5] +
	  a.A[6] + 2 * a.A[7];
	break;
      case e8:
	W172->val.A[8] = a.A[1];

	for (i = 7; i >= 2; i--)
	  {
	    W172->val.A[i] = W172->val.A[i + 1] + a.A[9 - i];
	  }
	W172->val.A[1] =
	  2 * a.A[1] + 3 * a.A[2] + 4 * a.A[3] + 5 * a.A[4] + 6 * a.A[5] +
	  4 * a.A[6] + 2 * a.A[7] + 3 * a.A[8];
	break;
      default:
	Caseerror (Line);
      }
  }
}

ocharptr
partlbl (int n, termptr list)
{
  register ocharptr R146;
  ocharptr dum, chrc, sub;

  chrc = NULL;
  sub = NULL;
  switch ((int) (group))
    {
    case an:
    case sn:
    case mp:
    case spnc:
    case un:
    case on:
    case unm:
    case sunm:
    case ospnm:
    case nill:
    case unc:
    case sonc:
      print ( "Inappropriate group\n");
      break;
    case son:
    case spn:
    case sung:
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
      while (list != NULL)
	{
	  cnu (&sub);
	  {
	    register ocharptr W189 = &(*sub);

	    frdynk (&sub, list->val, n);
	    W189->mult = list->mult;
	    W189->next = chrc;
	    chrc = sub;
	  }
	  list = list->next;
	}
      sub = NULL;
      while (chrc != NULL)
	{
	  dum = chrc;
	  chrc = chrc->next;
	  dum->next = sub;
	  sub = dum;
	}
      break;
    default:
      Caseerror (Line);
    }
  R146 = sub;
  return R146;
}

void
putprop (ocharptr chrc, bool ptod)
{
  int norml, order, qq;
  bframe dimn, dindex, casimir;
  termptr dynky;

  if ((ggroup.name != nill) && (chrc != NULL))
    {
      register groop *W190 = &ggroup;

      group = W190->name;
      norml = 0;
      order = 0;
      switch ((int) (W190->name))
	{
	case sung:		/* case un: */
	  norml = 2 * W190->rank * W190->rank;
	  order = W190->rank - 1;
	  break;
	case spn:
	  norml = 4 * (W190->rank / 2 + 1);
	  order = W190->rank / 2;
	  break;
	  /* case spnc:
	     norml = 4;
	     order = W190->rank / 2;
	     break ; */
	/* case on: */ case son:
	  norml = 8 * (W190->rank - 2);
	  order = W190->rank / 2;
	  break;
	case g2:
	  norml = 12;
	  order = 2;
	  break;
	case f4:
	  norml = 18;
	  order = 4;
	  break;
	case e6:
	  norml = 144;
	  order = 6;
	  break;
	case e7:
	  norml = 72;
	  order = 7;
	  break;
	case e8:
	  norml = 60;
	  order = 8;
	  break;
	case nill:
	case sn:
	case unm:
	case sunm:
	case un:
	case on:
	case ospnm:
	case mp:
	case unc:
	case spnc:
	case an:;
	  break;
	default:
	  Caseerror (Line);
	}
      switch ((int) (W190->name))
	{
	/*case un: */ case sung:	/* case on: */
	case son:
	case spn:
	case g2:
	case f4:
	case e6:
	case e7:
	case e8:
	  dynky = dynklbl (W190->rank, chrc);
	  print ("<dynkin label> ");
	  qq = 16;
	  wrtdlst (&output, &qq, tcol, order, dynky);
	  if (logging)
	    {
	      qq = 16;
	      wrtdlst (&logfile, &qq, pcol, order, dynky);
	    }
	  print ("\n");
	  ldisp (&dynky);
	  break;
	case nill:
	case sn:
	case unm:
	case sunm:
	case un:
	case on:
	case ospnm:
	case mp:
	case spnc:
	case unc:
	case sonc:
	case an:;
	  break;
	default:
	  Caseerror (Line);
	}
      if ((ptod || (W190->name == spnc)))
	{
	  switch ((int) (W190->name))
	    {
	    case un:
	    case sung:
	    case on:
	    case son:
	    case spn:
	    case g2:
	    case f4:
	    case e6:
	    case e7:
	    case e8:
	    case an:
	    case sn:
	    case l168:
	      dimngrp (W190->rank, chrc, &dimn);
	      print ( "dimension = ");
	      qq = 12;
	      wrtbigno (&output, &qq, &dimn);
	      if (logging)
	      {
		qq = 12;
                wrtbigno (&logfile, &qq, &dimn);
	      }
	      break;
	    case nill:
	    case unm:
	    case sunm:
	    case ospnm:
	    case mp:
	    case spnc:
	    case sonc:
	    case unc:
	      print ( "Inappropriate group. No action taken\n");
	      print ( "All values set to zero\n");
	      break;
	    default:
	      Caseerror (Line);
	    }
	  switch ((int) (W190->name))
	    {
	    /* case un: */ case sung:	/*case on: */
	    case son:
	    case spn:
	    case g2:
	    case f4:
	    case e6:
	    case e7:
	    case e8:
	      print ( "   ");
	      qq = qq + 3;
	      if (chrc->next == NULL)
		{
		  casmgrp (W190->rank, chrc, &casimir);
		  print ( "%1d*2nd-casimir=", norml);
		  qq = qq + 18;
		  if (logging)
		  {int saveqq=qq;
		    wrtbigno (&logfile, &saveqq, &casimir);
		  }
		  wrtbigno (&output, &qq, &casimir);
		  print("\n");
		}
	      break;
	    case nill:
	    case sn:
	    case unm:
	    case sunm:
	    case spnc:
	    case ospnm:
	    case mp:
	    case an:
	    case un:
	    case on:;
	      break;
	    default:
	      Caseerror (Line);
	    }
	  switch ((int) (W190->name))
	    {
	    /* case un: */ case sung:	/*case on: */
	    case son:
	    case spn:
	    case g2:
	    case f4:
	    case e6:
	    case e7:
	    case e8:
	      d2index (W190->rank, chrc, &dindex);
	      print ( "2nd-dynkin = ");
	      if (logging)
               {int saveqq=qq+13;
                 wrtbigno (&logfile, &saveqq, &dindex);
	       }
	      qq = qq + 13;
	      wrtbigno (&output, &qq, &dindex);
	      print("\n");
	      break;
	    case nill:
	    case sn:
	    case unm:
	    case sunm:
	    case ospnm:
	    case mp:
	    case spnc:
	    case an:
	    case un:
	    case on:;
	      break;
	    default:
	      Caseerror (Line);
	    }
	}
	      print("\n");
    }
}


void
putdindex (ocharptr chrc)
{
  int /*norml, order, */ qq;
  bframe dindex;

  if ((chrc != NULL))
    {
      register groop *W191 = &ggroup;

      group = W191->name;
      //norml = 0;
      //order = 0;
      switch ((int) (W191->name))
	{
	case sung:		/* case un: */
	  //norml = 2 * W191->rank * W191->rank;
	  //order = W191->rank - 1;
	  break;
	case spn:
	  //norml = 4 * (W191->rank / 2 + 1);
	  //order = W191->rank / 2;
	  break;
	/* case on: */ case son:
	  //norml = 8 * (W191->rank - 2);
	  //order = W191->rank / 2;
	  break;
	case g2:
	  //norml = 12;
	  //order = 2;
	  break;
	case f4:
	  //norml = 18;
	  //order = 4;
	  break;
	case e6:
	  //norml = 144;
	  //order = 6;
	  break;
	case e7:
	  //norml = 72;
	  //order = 7;
	  break;
	case e8:
	  //norml = 60;
	  //order = 8;
	  break;
	  /*  case spnc:
	     norml = 4;
	     order = W191->rank / 2;
	     break ; */

	case an:
	case sn:
	case mp:
	case spnc:
	case sonc:
	case unc:
	case l168:
	case unm:
	case sunm:
	case ospnm:
	case nill:
	case un:
	case on:;
	  break;
	default:
	  Caseerror (Line);
	}
      switch ((int) (W191->name))
	{
	case sung:		/* case un: */
	case spn:		/* case on: */
	case son:
	case g2:
	case f4:
	case e6:
	case e7:
	case e8:
	  d2index (W191->rank, chrc, &dindex);
	  print ("2nd-dynkin = ");
	  qq = 13;
	  wrtbigno (&output, &qq, &dindex);
	  if (logging)
	  {
	    qq = 13;
	    wrtbigno (&logfile, &qq, &dindex);
	  }
	  print("\n");
	  break;
	  /*  case nill:
	     warn("group not set:No", cont);
	     inform(" action taken;", cr);
	     break ; */
	case an:
	case sn:
	case mp:
	case spnc:
	case l168:
	case unm:
	case sunm:
	case ospnm:
	case un:
	case on:
	  print ( "Inappropriate group. No action taken\n");
	  print ( "All values set to zero\n");
	  break;
	default:
	  Caseerror (Line);
	}
    }
}


void
dimnprop (ocharptr chrc)
{
  int qq;
  bframe dimn;

  if ((chrc != NULL))
    {
      register groop *W192 = &ggroup;

      group = W192->name;
      switch ((int) (W192->name))
	{
	case an:
	case sn:
	case un:
	case sung:
	case on:
	case son:
	case spn:
	case g2:
	case f4:
	case e6:
	case e7:
	case e8:
	case l168:
	  dimngrp (W192->rank, chrc, &dimn);
	  print ( "dimension = ");
	  qq = 12;
	  wrtbigno (&output, &qq, &dimn);
	  if (logging)
	    {
	      qq = 12;
	      wrtbigno (&logfile, &qq, &dimn);
	    }
	  print ("\n");
	  break;
	case mp:
	case spnc:
	case unm:
	case sunm:
	case ospnm:
	case nill:
	case sonc:
	case unc:
	  print ( "Inappropriate group. No action taken\n");
	  print ( "All values set to zero\n");
	  break;
	default:
	  Caseerror (Line);
	}
    }
}


void
utrace (prodtype pr)
{
  int j;
  bframe dimn, sdimn, subdimn;

  tobig (0, &dimn);
  while (pr != NULL)
    {
      register prodtype W193 = &(*pr);

      factor (1, &subdimn);
      for (j = 1; j <= nprod; j++)
	{
	  if ((currgrp.A[j - 1].rank > 1))
	    {
	      tobig (0, &sdimn);
	      group = currgrp.A[j - 1].name;
	      dimngrp (currgrp.A[j - 1].rank, W193->prods.A[j - 1], &sdimn);
	      bigtofact (sdimn, &sdimn);
	      fmult (subdimn, sdimn, &subdimn);
	    }
	  else if ((currgrp.A[j - 1].rank == 1))
	    kmult (subdimn, W193->prods.A[j - 1]->val.A[1], &subdimn);
	}
      kmult (subdimn, W193->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, dimn, &dimn);
      pr = pr->next;
    }
  inform ("Utrace = ;", cont);
  j = 12;
  wrtbigno (&output, &j, &dimn);
  Putchr ('\n', output);
  if (logging)
    {
      wrtbigno (&logfile, &j, &dimn);
      Putchr ('\n', logfile);
    }
}

void
pproperties (prodtype pr)
{
  int j;
  bframe dimn, sdimn, subdimn;

  tobig (0, &dimn);
  while (pr != NULL)
    {
      register prodtype W196 = &(*pr);

      factor (1, &subdimn);
      for (j = 1; j <= nprod; j++)
	{
	  tobig (0, &sdimn);
	  group = currgrp.A[j - 1].name;
	  dimngrp (currgrp.A[j - 1].rank, W196->prods.A[j - 1], &sdimn);
	  bigtofact (sdimn, &sdimn);
	  fmult (subdimn, sdimn, &subdimn);
	}
      kmult (subdimn, W196->mult, &subdimn);
      factobig (subdimn, &subdimn);
      bigadd (subdimn, dimn, &dimn);
      pr = pr->next;
    }
  inform ("Dimension = ;", cont);
  j = 12;
  wrtbigno (&output, &j, &dimn);
  Putchr ('\n', output);
  if (logging)
    {
      wrtbigno (&logfile, &j, &dimn);
      Putchr ('\n', logfile);
    }
}

void
adjoint (frame * adj)
{
  register int j;
  for (j = 0; j <= maxdim; j++)
    {
      adj->A[j] = 0;
    }
  {
    register groop *W201 = &ggroup;

    switch ((int) (W201->name))
      {
      /*case un: */ case sung:
	adj->A[1] = 2;
	for (j = 2; j <= W201->rank - 1; j++)
	  {
	    adj->A[j] = 1;
	  }
	break;
      /* case on: */ case son:
	adj->A[1] = 1;
	adj->A[2] = 1;
	break;
      case spn:
	adj->A[1] = 2;
	break;
      case g2:
	adj->A[1] = 2;
	adj->A[2] = 1;
	break;
      case f4:
	adj->A[1] = 1;
	adj->A[2] = 1;
	break;
      case e6:
	adj->A[1] = 2;
	break;
      case e7:
	adj->A[1] = 2;
	for (j = 2; j <= 7; j++)
	  {
	    adj->A[j] = 1;
	  }
	break;
      case e8:
	adj->A[1] = 2;
	for (j = 2; j <= 8; j++)
	  {
	    adj->A[j] = 1;
	  }
	break;
      case unc:
      case spnc:
      case mp:
      case sonc:
      case nill:
      case sn:
      case an:
      case unm:
      case sunm:
      case ospnm:
      case un:
      case on:;
	break;
      default:
	Caseerror (Line);
      }
  }
}

void
dynkinn (int rank, int n, ocharptr mu, ocharptr lam, bframe * dynkn)
{
  bframe tsums, sum, clam, cmu, clm, cn, ct, pct, dnu, dlam;
  ocharptr help, nu, top;
  bool neg;
  register int j;

  tobig (0, &tsums);
  cnu (&help);
  casmgrp (rank, mu, &cmu);
  dimngrp (rank, lam, &dlam);
  while (lam != NULL)
    {
      (*help) = (*lam);
      help->mult = 1;
      help->next = NULL;
      nu = kronk (mu, help, ggroup);
      casmgrp (rank, help, &clam);
      bigadd (clam, cmu, &clm);
      tobig (0, &sum);
      top = nu;
      while (top != NULL)
	{
	  (*help) = (*top);
	  help->next = NULL;
	  dimngrp (rank, help, &dnu);
	  casmgrp (rank, help, &cn);
	  bigsubtr (cn, clm, &ct, &neg);
	  if (n == 0)
	    tobig (1, &ct);
	  bigtofact (ct, &ct);
	  pct = ct;
	  if (n > 0)
	    for (j = 1; j <= n - 1; j++)
	      {
		fmult (ct, pct, &ct);
	      }
	  bigtofact (dnu, &dnu);
	  fmult (ct, dnu, &ct);
	  factobig (ct, &ct);
	  bigadd (sum, ct, &sum);
	  top = top->next;
	}
      cmult (sum, lam->mult, &sum);
      bigadd (tsums, sum, &tsums);
      lam = lam->next;
      odisp (&nu);
    }
  bigtofact (tsums, &tsums);
  bigtofact (dlam, &dlam);
  fdiv (tsums, dlam, &tsums);
  tsums.A[1] = tsums.A[1] - n;
  factobig (tsums, &tsums);
  odisp (&help);
  (*dynkn) = tsums;
}

void
casmnsymp (int n, ocharptr list, bframe * cas)
{
  int k, nn;
  register int i;

  {
    register ocharptr W210 = &(*list);

    k = W210->val.A[maxdim];
    nn = n / 2;
    if (W210->spin)
      k = 2 * k + 1;
    else
      k = 2 * k;
    tobig (0, &(*cas));
    for (i = 1; i <= n / 2; i++)
      {
	cadd ((*cas),
	      (2 * (W210->val.A[i] + nn + 1 - i) +
	       k) * (2 * (W210->val.A[i] + nn + 1 - i) + k) - 4 * i * i,
	      &(*cas));
      }
  }
}

void
diml168 (ocharptr charlist, bframe * dimn)
{
  int j = 0;			// j not initialized  corrected FB
  bframe subdim;

  tobig (0, &(*dimn));
  tobig (0, &subdim);
  while (charlist != NULL)
    {
      register ocharptr W211 = &(*charlist);
      if (W211->val.A[1] == 0)
	j = 1;
      if (W211->val.A[1] == 1)
	j = 6;

      if (W211->val.A[1] == 2)
	j = 7;

      if (W211->val.A[1] == 3)
	j = 8;

      if (W211->val.A[1] == 4)
	j = 3;

      if (W211->val.A[1] == 5)
	j = 3;
      if (j == 0)		// coherency test added by FB
	fprintf (stderr, "ERROR j=0 in diml168\n"), exit (1);
      j = j * W211->mult;
      tobig (j, &subdim);
      bigadd ((*dimn), subdim, &(*dimn));
      charlist = W211->next;
    }
}
