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
# include <stdlib.h>
# include <unistd.h>
# include <math.h>
# include "standard.h"
# include "define.h"
# include "ReadWrite.h"
# include "sets_mgmt.h"

/*
**	Start of program definitions
*/
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "s2.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "s6.h"
# include "g.h"
# include "skew.h"
# include "s4.h"

tabptr *G183_brtab;

ocharptr
formbb (termptr con, termptr cov, bool suppress, bool spinor, char lable)
{
  register ocharptr R177;
  ocharptr lastptr, list;
  termptr contr;

  lastptr = NULL;
  while (cov != NULL)
    {
      contr = con;
      while (contr != NULL)
	{
	  cnu (&list);
	  {
	    register ocharptr W2 = &(*list);

	    if ((contr->val.A[1] == 0) && suppress)
	      W2->C6_double = false;
	    else
	      {
		W2->C6_double = true;
		W2->conval = contr->val;
	      }
	    W2->val = cov->val;
	    W2->mult = contr->mult * cov->mult;
	    W2->spin = spinor;
	    W2->lab = lable;
	    W2->next = lastptr;
	    lastptr = list;
	    contr = contr->next;
	  }
	}
      cov = cov->next;
    }
  R177 = lastptr;
  return R177;
}

/*ocharptr spnunbr();*/
ocharptr
spnunbr1 (ocharptr lambda, int n)
{
  register ocharptr R179;
  ocharptr list, sublist;
  termptr slam, xi, xit, temp1, temp2, temp, bb;
  int u1val;
  register int i;

  n = n / 2;
  list = NULL;
  snu (&slam);
  {
    register termptr W3 = &(*slam);

    W3->val = lambda->val;
    W3->mult = lambda->mult;
    W3->next = NULL;
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
	  bb = seriesx ('d', -1, temp1->val);
	  temp = temp1;
	  temp1 = temp1->next;
	  temp->next = NULL;
	  temp2 = lskew (temp, bb);
	  ldisp (&bb);
	  dispsfn (&temp);
	  sublist = formbb (xit, temp2, false, false, ' ');
	  gumodify (&sublist, n);
	  oadd (&list, &sublist);
	  ldisp (&temp2);
	}
      dispsfn (&xit);
    }
  sublist = list;
  while (sublist != NULL)
    {
      register ocharptr W4 = &(*sublist);

      if (W4->C6_double)
	{
	  u1val = 2 * (wtfrm (&W4->val) - wtfrm (&W4->conval));
	  for (i = 1; i <= n; i++)
	    {
	      W4->val.A[i] =
		W4->val.A[i] + W4->conval.A[1] - W4->conval.A[n - i + 1];
	    }
	  W4->C6_double = false;
	}
      else
	{
	  u1val = 2 * wtfrm (&W4->val);
	  for (i = 1; i <= n; i++)
	    {
	      W4->val.A[i] = W4->val.A[i] - W4->val.A[n];
	    }
	}
      for (i = maxdim - 1; i >= 1; i--)
	{
	  W4->val.A[i + 1] = W4->val.A[i];
	}
      W4->val.A[1] = u1val;
      sublist = W4->next;
    }
  dispsfn (&slam);
  osort (&list, true);
  R179 = list;
  return R179;
}

ocharptr
spnunbr (ocharptr lambda, int n)
{
  register ocharptr R178;
  ocharptr brnch, subbrnch;

  brnch = NULL;
  while (lambda != NULL)
    {
      register ocharptr W11 = &(*lambda);

      subbrnch = spnunbr1 (lambda, n);
      oadd (&brnch, &subbrnch);
      lambda = W11->next;
    }
  osort (&brnch, true);
  R178 = brnch;
  return R178;
}

void
fixcg (groop * grp, grptype g, int r1, int r2)
{
  {
    register groop *W12 = &(*grp);

    W12->name = g;
    W12->rank = r1;
    W12->rank2 = r2;
  }
}

void
fixochar (ocharptr * help, int m, frame v, frame c, char l, char cl,
	  bool s, bool d)
{
  cnu (&(*help));
  {
    register ocharptr W13 = &(*(*help));

    W13->mult = m;
    W13->val = v;
    W13->conval = c;
    W13->lab = l;
    W13->conlab = cl;
    W13->spin = s;
    W13->C6_double = d;
    W13->next = NULL;
  }
}

void
brload1 (text * gfile, tabptr * gtab, ocharptr * gindex, bool * gload)
{
  int i, ppp;
  string0 buff2;
  ocharptr indexb, dummy;

  {
    register tabptr W14 = &(*(*gtab));

    (*gindex) = NULL;
    i = 0;
    do
      {
	i = i + 1;
	W14->tab.A[i - 1] = Getchr ((*gfile));
	if (i % 80 == 0)
	  Getl (&(*gfile));
      }
    while (!((W14->tab.A[i - 1] == cont) || (i == brtabmax)));
    //fprintf(stderr,"loading data i=%d\n",i);
    if (i == brtabmax)
      aaargh ("branchtable full", false);
    if (i % 80 != 0)
      Getl (&(*gfile));
    readacard (&(*gfile), &buff2, &ppp);
    readachrc (&(*gfile), &buff2, &ppp, &indexb);
    while (indexb != NULL)
      {
	dummy = indexb->next;
	indexb->next = (*gindex);
	(*gindex) = indexb;
	indexb = dummy;
      }
    (*gload) = true;
    fclose (gfile->fp);
     /**/			/*13/12/95 */
  }
}

void
brload (int bno)
{
/* --------------------------------------------------------- */

  char *tmpstring  ;
  char tmp[80];

  tmpstring = tmp;
  tmpstring[0] = '\0';

  sprintf (tmpstring, "%s/dat/", dataPath);
  if (debug_schur)
	  fprintf(stderr,"datapath=%s bno=%d\n",tmpstring,bno);
/* --------------------------------------------------------- */

  if (bno == 43)
    {
      strcat (tmpstring, "f4file.dat");
      Resetx2 (&f4file, tmpstring, -1);
      f4tab = (tabptr) malloc ((unsigned) (sizeof (*f4tab)));
      brload1 (&f4file, &f4tab, &f4index, &f4load);
    }
  else if (bno == 44)
    {
      strcat (tmpstring, "e6file.dat");
      Resetx2 (&e6file, tmpstring, -1);
      e6tab = (tabptr) malloc ((unsigned) (sizeof (*e6tab)));
      brload1 (&e6file, &e6tab, &e6index, &e6load);
    }
  else if (bno == 47)
    {
      strcat (tmpstring, "e7file.dat");
      Resetx2 (&e7file, tmpstring, -1);
      e7tab = (tabptr) malloc ((unsigned) (sizeof (*e7tab)));
      brload1 (&e7file, &e7tab, &e7index, &e7load);
    }
  else if (bno == 49)
    {
      strcat (tmpstring, "e8file.dat");
      Resetx2 (&e8file, tmpstring, -1);
      e8tab = (tabptr) malloc ((unsigned) (sizeof (*e8tab)));
      brload1 (&e8file, &e8tab, &e8index, &e8load);
    }
  else if (bno == 50)
    {
      strcat (tmpstring, "e8so16.dat");
      Resetx2 (&e8so16, tmpstring, -1);
      e8sotab = (tabptr) malloc ((unsigned) (sizeof (*e8sotab)));
      brload1 (&e8so16, &e8sotab, &e8soindex, &e8soload);
    }
  else if (bno == 51)
    {
      strcat (tmpstring, "e8su2e7.dat");
      Resetx2 (&e8su2e7, tmpstring, -1);
      e8sutab = (tabptr) malloc ((unsigned) (sizeof (*e8sutab)));
      brload1 (&e8su2e7, &e8sutab, &e8suindex, &e8suload);
    }
  else if (bno == 52)
    {
      strcat (tmpstring, "e8su3e6.dat");
      Resetx2 (&e8su3e6, tmpstring, -1);
      e8e6tab = (tabptr) malloc ((unsigned) (sizeof (*e8e6tab)));
      brload1 (&e8su3e6, &e8e6tab, &e8e6index, &e8e6load);
    }
  else if (bno == 45)
    {
      strcat (tmpstring, "e6so10.dat");
      Resetx2 (&e6so10, tmpstring, -1);
      e6sotab = (tabptr) malloc ((unsigned) (sizeof (*e6sotab)));
      brload1 (&e6so10, &e6sotab, &e6soindex, &e6soload);
    }
  else if (bno == 48)
    {
      strcat (tmpstring, "e7e6file.dat");
      Resetx2 (&e7e6file, tmpstring, -1);
      e7e6tab = (tabptr) malloc ((unsigned) (sizeof (*e7e6tab)));
      brload1 (&e7e6file, &e7e6tab, &e7e6index, &e7e6load);
    }
  else if (bno == 46)
    {
      strcat (tmpstring, "e6g2file.dat");
      Resetx2 (&e6g2file, tmpstring, -1);
      e6g2tab = (tabptr) malloc ((unsigned) (sizeof (*e6g2tab)));
      brload1 (&e6g2file, &e6g2tab, &e6g2index, &e6g2load);
    }
  else if (bno == 53)
    {
      strcat (tmpstring, "u27e6.dat");
      Resetx2 (&u27e6, tmpstring, -1);
      u27e6tab = (tabptr) malloc ((unsigned) (sizeof (*u27e6tab)));
      brload1 (&u27e6, &u27e6tab, &u27e6index, &u27e6load);
    }
  else if (bno == 54)
    {
      strcat (tmpstring, "su56e7.dat");
      Resetx2 (&su56e7, tmpstring, -1);
      su56e7tab = (tabptr) malloc ((unsigned) (sizeof (*su56e7tab)));
      brload1 (&su56e7, &su56e7tab, &su56e7index, &su56e7load);
    }
  else if (bno == 55)
    {
      strcat (tmpstring, "su248e8.dat");
      Resetx2 (&su248e8, tmpstring, -1);
      su248e8tab = (tabptr) malloc ((unsigned) (sizeof (*su248e8tab)));
      brload1 (&su248e8, &su248e8tab, &su248e8index, &su248e8load);
    }
  else if (bno == 56)
    {
      strcat (tmpstring, "e6f4.dat");
      Resetx2 (&e6f4, tmpstring, -1);
      e6f4tab = (tabptr) malloc ((unsigned) (sizeof (*e6f4tab)));
      brload1 (&e6f4, &e6f4tab, &e6f4index, &e6f4load);
    }
  else if (bno == 57)
    {
      strcat (tmpstring, "f4g2.dat");
      Resetx2 (&f4g2, tmpstring, -1);
      f4g2tab = (tabptr) malloc ((unsigned) (sizeof (*f4g2tab)));
      brload1 (&f4g2, &f4g2tab, &f4g2index, &f4g2load);
    }
  if (bno == 58)
    {
      strcat (tmpstring, "e6su3g2.dat");
      Resetx2 (&e6su3g2, tmpstring, -1);
      e6su3g2tab = (tabptr) malloc ((unsigned) (sizeof (*e6su3g2tab)));
      brload1 (&e6su3g2, &e6su3g2tab, &e6su3g2index, &e6su3g2load);

    }
  if (bno == 60)
    {
      strcat (tmpstring, "e8f4g2.dat");
      Resetx2 (&e8f4g2, tmpstring, -1);
      e8f4g2tab = (tabptr) malloc ((unsigned) (sizeof (*e8f4g2tab)));
      brload1 (&e8f4g2, &e8f4g2tab, &e8f4g2index, &e8f4g2load);

    }
  if (bno == 59)
    {
      strcat (tmpstring, "so26f4.dat");
      Resetx2 (&so26f4, tmpstring, -1);
      so26f4tab = (tabptr) malloc ((unsigned) (sizeof (*so26f4tab)));
      brload1 (&so26f4, &so26f4tab, &so26f4index, &so26f4load);
    }
  if (bno == 63)
    {
      strcat (tmpstring, "l168.dat");
      Resetx2 (&l168file, tmpstring, -1);
      l168tab = (tabptr) malloc ((unsigned) (sizeof (*l168tab)));
      brload1 (&l168file, &l168tab, &l168index, &l168load);
    }
}



/*ocharptr extract();*/

void
extrno (int *n, int *num)
{
  bool neg;

  neg = (bool) ((*G183_brtab)->tab.A[(*n) - 1] == '~');
  if (neg)
    (*n) = (*n) + 1;
  if ((*G183_brtab)->tab.A[(*n) - 1] == '*')
    {
      (*num) =
	((unsigned) ((*G183_brtab)->tab.A[(*n) + 1 - 1]) -
	 (unsigned) ('0')) * 10 + (unsigned) ((*G183_brtab)->tab.A[(*n) + 2 -
								   1]) -
	(unsigned) ('0');
      (*n) = (*n) + 3;
    }
  else
    {
      (*num) = (unsigned) ((*G183_brtab)->tab.A[(*n) - 1]) - (unsigned) ('0');
      (*n) = (*n) + 1;
    }
  if (neg)
    (*num) = -(*num);
}

ocharptr
extract (int n, int multy, tabptr brtab)
{
  register ocharptr R180;
  ocharptr charc, lastptr;
  register int i;
  int temp;
  tabptr *F184;

  F184 = G183_brtab;
  G183_brtab = &brtab;
  lastptr = NULL;
  {
    register tabptr W15 = &(*(*G183_brtab));

    do
      {
	cnu (&charc);
	{
	  register ocharptr W16 = &(*charc);

	  W16->spin = false;
	  W16->C6_double = false;
	  W16->val = nolls;
	  W16->conval = nolls;
	  W16->lab = ' ';
	  W16->conlab = ' ';
	  n = n + 1;
	  W16->next = lastptr;
	  extrno (&n, &W16->mult);
	  W16->mult = W16->mult * multy;
	  /*if (!(Member((unsigned)(W15->tab.A[n - 1]), Conset[0]))) */
	  if (((W15->tab.A[n - 1] != 's') && (W15->tab.A[n - 1] != 'S')))
	    W16->spin = false;
	  else
	    {
	      W16->spin = true;
	      n = n + 1;
	    }
	  for (i = 0; i <= maxdim; i++)
	    {
	      W16->val.A[i] = 0;
	      //bval.A[i] = 0;
	    }
	  i = 1;
	  do
	    {
	      extrno (&n, &temp);

	      if ((temp <= 127) && (temp >= -128))
	      {
		W16->val.A[i] = temp;
	        //bval.A[i] = temp;
	      }
	      /*else
		{
		  fprintf (stderr, "Error must be <= 127 and >=-128\n");
		  return (NULL);
		}*/
	      i = i + 1;
	    }
	  while (!
		 ((W15->tab.A[n - 1] != '*')
		  && !(Member ((unsigned) (W15->tab.A[n - 1]), numbers.S))));
	  if (!(W15->tab.A[n - 1] == ';'))
	    W16->C6_double = false;
	  else
	    {
	      n = n + 1;
	      W16->C6_double = true;
	      /*if (!(Member((unsigned)(W15->tab.A[n - 1]), Conset[1]))) */
	      if (((W15->tab.A[n - 1] != 's') && (W15->tab.A[n - 1] != 'S')))
		W16->spin = false;
	      else
		{
		  W16->spin = true;
		  n = n + 1;
		}
	      for (i = 0; i <= maxdim; i++)
		{
		  W16->conval.A[i] = 0;
		  //bconval.A[i] = 0;
		}
	      i = 1;
	      do
		{
		  extrno (&n, &temp);
		  if ((temp <= 127) && (temp >= -128))
		  {
		    W16->conval.A[i] = temp;
	            //bconval.A[i] = temp;
		  }
		  i = i + 1;
		}
	      while (!
		     ((W15->tab.A[n - 1] != '*')
		      &&
		      !(Member ((unsigned) (W15->tab.A[n - 1]), numbers.S))));
	    }
	  /*if ((Member((unsigned)(W15->tab.A[n - 1]), Conset[2]))) */
	  if (((W15->tab.A[n - 1] == ',') || (W15->tab.A[n - 1] == '/')))
	    W16->lab = ' ';
	  else
	    {
	      W16->lab = W15->tab.A[n - 1];
	      n = n + 1;
	    }
	  lastptr = charc;
	}
      }
    while (!(W15->tab.A[n - 1] == '/'));
  }
  R180 = charc;
  G183_brtab = F184;
  return R180;
}

void
exceptgroupload (int *brno)
{
  if (((*brno) == 43) && (!f4load))
    brload (43);
  if (((*brno) == 44) && (!e6load))
    brload (44);
  if (((*brno) == 45) && (!e6soload))
    brload (45);
  if (((*brno) == 46) && (!e6g2load))
    brload (46);
  if (((*brno) == 47) && (!e7load))
    brload (47);
  if (((*brno) == 48) && (!e7e6load))
    brload (48);
  if (((*brno) == 49) && (!e8load))
    brload (49);
  if (((*brno) == 50) && (!e8soload))
    brload (50);
  if (((*brno) == 51) && (!e8suload))
    brload (51);
  if (((*brno) == 52) && (!e8e6load))
    brload (52);
  if (((*brno) == 53) && (!u27e6load))
    brload (53);
  if (((*brno) == 54) && (!su56e7load))
    brload (54);
  if (((*brno) == 55) && (!su248e8load))
    brload (55);
  if (((*brno) == 56) && (!e6f4load))
    brload (56);
  if (((*brno) == 57) && (!f4g2load))
    brload (57);
  if (((*brno) == 58) && (!e6su3g2load))
    brload (58);
  if (((*brno) == 59) && (!so26f4load))
    brload (59);
  if (((*brno) == 60) && (!e8f4g2load))
    brload (60);
  if (((*brno) == 63) && (!l168load))
    brload (63);
}

void
fixg (groop * grp, grptype g, int r1, int r2)
{
  {
    register groop *W21 = &(*grp);

    W21->name = g;
    W21->rank = r1;
    W21->rank2 = r2;
  }
}

/*ocharptr unso3brnch();*/

/*void plethp();*/

void
expand (int startx, int stop, int a, barr * expansion)
{
  int length;
  register int j;
  register int i;
  barr helper;
  bool first;

  first = true;
  for (i = 1; i <= maxsize; i++)
    {
      expansion->A[i] = 0;
    }
  expansion->A[0] = 1;
  length = 0;
  for (i = startx; i <= stop; i++)
    {
      if (length + i + a > maxsize)
	{
	  if (first)
	    print ("possible internal error - check result\n");
	  length = length - i - a;
	  first = false;
	}
      for (j = 0; j <= maxsize; j++)
	{
	  helper.A[j] = 0;
	}
      for (j = 0; j <= length; j++)
	{
	  helper.A[j + a + i] = expansion->A[j];
	}
      length = length + i + a;
      for (j = 0; j <= length; j++)
	{
	  expansion->A[j] = expansion->A[j] - helper.A[j];
	}
    }
  progress ();
}

void
divide (int lim, barr n, barr d, barr * res)
{
  barr nc;
  register int k;
  register int j;

  progress ();
  for (j = 0; j <= lim; j++)
    {
      nc.A[j] = n.A[j];
    }
  for (j = 0; j <= lim; j++)
    {
      res->A[lim - j] = nc.A[j];
      for (k = j; k <= lim; k++)
	{
	  nc.A[k] = nc.A[k] - res->A[lim - j] * d.A[k - j];
	}
    }
}

void
evalpleth (int *numso3prts, int *k1, double p, barr * rese)
{
  register int j;
  barr a, b;
  for (j = 0; j <= maxsize; j++)
    {
      rese->A[j] = 0;
    }
  (*numso3prts) = Trunc (p * (*k1));
  expand (1, (*k1), Trunc (2 * p), &a);
  expand (2, (*k1), 0, &b);
  divide ((*numso3prts), a, b, &(*rese));
  progress ();
}

void
arrayprod (double hfy, int *m2y, int *m2prod, barr * resy, barr * prodx)
{
  int help = 0;			// modified by FB, was uninitialized
  register int i;
  register int k;
  register int j;
  barr helper;
  double startx, stop;

  progress ();
  for (j = 0; j <= maxsize; j++)
    {
      helper.A[j] = 0;
    }
  for (j = 0; j <= *m2prod; j++)
    {
      for (k = 0; k <= *m2y; k++)
	{
	  startx = fabs (k + hfy - j / ((double) 2));
	  stop = fabs (k + hfy + j / ((double) 2));
	  for (i = 0; i <= stop - startx; i++)
	    {
	      help = Trunc (2 * startx + 2 * i);
	      helper.A[help] = helper.A[help] + prodx->A[j] * resy->A[k];
	    }
	}
    }
  (*m2prod) = help;
  for (j = 0; j <= maxsize; j++)
    {
      prodx->A[j] = helper.A[j];
    }
}

void
pplethsfn (double p, frame sfn, int lensfn, int *m2prod, barr * so3parts)
{
  int numnon0hs, m2y, m2x, hfx;
  register int j;
  register int k;
  int pp[maxl];
  double hfy;
  barr resx, resy, prodx;
  termptr temp, list, top;

  progress ();
  snu (&list);
  list->val = sfn;
  list->mult = 1;
  temp = factorsof (list->mult, list->val);
  dispsfn (&list);
  list = temp;
  for (k = 0; k <= maxsize; k++)
    {
      so3parts->A[k] = 0;
    }
  top = list;
  while (list != NULL)
    {
      numnon0hs = 0;
      for (j = 1; j <= lensfn; j++)
	{
	  if (list->val.A[j] != 0)
	    {
	      numnon0hs = numnon0hs + 1;
	      pp[numnon0hs] = list->val.A[j];
	    }
	}
      evalpleth (&m2x, &pp[1], p, &resx);
      progress ();
      if ((pp[1] * p) != Trunc (pp[1] * p))
	hfx = 1;
      else
	hfx = 0;
      for (j = 0; j <= maxsize; j++)
	{
	  prodx.A[j] = 0;
	}
      for (j = 0; j <= m2x; j++)
	{
	  k = 2 * j + hfx;
	  prodx.A[k] = resx.A[j];
	}
      (*m2prod) = 2 * m2x + hfx;
      for (j = 2; j <= numnon0hs; j++)
	{
	  evalpleth (&m2y, &pp[j], p, &resy);
	  progress ();
	  if ((pp[j] * p) != Trunc (pp[j] * p))
	    hfy = 0.5;
	  else
	    hfy = 0;
	  arrayprod (hfy, &m2y, &(*m2prod), &resy, &prodx);
	}
      for (j = 0; j <= *m2prod; j++)
	{
	  so3parts->A[j] = so3parts->A[j] + prodx.A[j] * list->mult;
	}
      list = list->next;
    }
  ldisp (&top);
}

void
plethp (double p, int lensfn, frame sfn, barr * so3parts, int *numso3prts)
{
  pplethsfn (p, sfn, lensfn, &(*numso3prts), &(*so3parts));
}

void
gradtest (ocharptr list, int n)
{
  int lenlist, lentemp = 0, k;
  register int i;
  ocharptr temp;

  lenlist = len (&list->val);
  cnu (&temp);
  for (i = 1; i <= maxdim; i++)
    {
      temp->val.A[i] = 0;
    }
  for (i = 0; i <= n; i++)
    {
	k = n - i;
	temp->val.A[k] = list->val.A[1] - list->val.A[i + 1];
	lentemp = len (&temp->val);
    }
  if (lentemp < lenlist)
    list->val = temp->val;
  odisp (&temp);
}

ocharptr
unso3brnch (ocharptr list1, int n)
{
  int lenval, numso3irreps, i;
  register int l;
  register int k;
  register int j;
  double p;
  barr so3irreps;
  ocharptr top, entry;
  barr helper;

  progress ();
  for (j = 0; j <= maxsize; j++)
    {
      so3irreps.A[j] = 0;
    }
  while (list1 != (NULL))
    {
      register ocharptr W69 = &(*list1);

      //weightx = wtfrm (&W69->val);
      lenval = len (&W69->val);
      if (lenval <= n)
	{
	  if ((lenval == n))
	    for (i = 1; i <= n; i++)
	      {
		W69->val.A[i] = W69->val.A[i] - W69->val.A[n];
	      }
	  gradtest (list1, n);
	  p = (n - 1) / ((double) 2);
	  /*spinirrep = false;
	  if (p != Trunc (p))
	    spinirrep = (bool) (weightx != ((weightx / 2) * 2));*/
	  if ((W69->val.A[1] == 0))
	    so3irreps.A[0] = so3irreps.A[0] + W69->mult;
	  else
	    {
	      plethp (p, lenval, W69->val, &helper, &numso3irreps);
	      for (j = 0; j <= numso3irreps; j++)
		{
		  so3irreps.A[j] = so3irreps.A[j] + W69->mult * helper.A[j];
		}
	    }
	}
      list1 = W69->next;
    }
  top = NULL;
  for (k = maxsize; k >= 0; k--)
    {
      if (so3irreps.A[k] != 0)
	{
	  cnu (&entry);
	  {
	    register ocharptr W76 = &(*entry);
	    for (l = 0; l <= maxdim; l++)
	      {
		W76->val.A[l] = 0;
	      }
	    W76->C6_double = false;
	    W76->spin = (bool) (!(k == ((k / 2) * 2)));
	    W76->lab = ' ';
	    W76->mult = so3irreps.A[k];
	    W76->val.A[1] = k / 2;
	    W76->next = top;
	    progress ();
	  }
	  top = entry;
	}
    }
  return top;
}
