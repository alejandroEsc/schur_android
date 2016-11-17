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
# include <unistd.h>

# include "standard.h"
# include "define.h"
# include "ReadWrite.h"
# include "../config.h"

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
# include "s4.h"
# include "s5.h"
# include "s6.h"
# include "s8.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "q1.h"
# include "q2.h"
# include "bignums.h"
# include "branch.h"
# include "gr.h"
# include "g.h"
# include "dimensions.h"
# include "tableaux.h"
# include "s0.h"

prodtype discardedprds, pcurrent;

void
getsset1 (text * fyle, string0 * buff, int *p, int *n)
{
  termptr temp;

  if (((*n) >= 1) && ((*n) <= svarlimit))
    {
      temp = getsfn (&(*fyle), &(*buff), &(*p));
      ldisp (&svar.A[(*n) - 1]);
      svar.A[(*n) - 1] = temp;
      if (srjctindex > 0)
	if (svar.A[(*n) - 1] == sreject.A[srjctindex - 1])
	  {
	    sreject.A[srjctindex - 1] = NULL;
	    srjctindex = srjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
getset1 (text * fyle, string0 * buff, int *p, int *n)
{
  ocharptr temp;

  if (((*n) >= 1) && ((*n) <= rvarlimit))
    {
      temp = getchrc (&(*fyle), &(*buff), &(*p));
      odisp (&vari.A[(*n) - 1]);
      vari.A[(*n) - 1] = temp;
      if (rjctindex > 0)
	if (vari.A[(*n) - 1] == reject.A[rjctindex - 1])
	  {
	    reject.A[rjctindex - 1] = NULL;
	    rjctindex = rjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
getpset1 (text * fyle, string0 * buff, int *p, int *n)
{
  prodtype temp;

  if (((*n) >= 1) && ((*n) <= pvarlimit))
    {
      temp = getprod (&(*fyle), &(*buff), &(*p));
      pdisp (&pvar.A[(*n) - 1]);
      pvar.A[(*n) - 1] = temp;
      if (prjctindex > 0)
	if (pvar.A[(*n) - 1] == preject.A[prjctindex - 1])
	  {
	    preject.A[prjctindex - 1] = NULL;
	    prjctindex = prjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
getsset (text * fyle, string0 * buff, int *p)
{
  int n;
  termptr temp;

  readint (&(*buff), &(*p), &n);
  if ((n >= 1) && (n <= svarlimit))
    {
      temp = getsfn (&(*fyle), &(*buff), &(*p));
      ldisp (&svar.A[n - 1]);
      svar.A[n - 1] = temp;
      if (srjctindex > 0)
	if (svar.A[n - 1] == sreject.A[srjctindex - 1])
	  {
	    sreject.A[srjctindex - 1] = NULL;
	    srjctindex = srjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
getset (text * fyle, string0 * buff, int *p)
{
  int n;
  ocharptr temp;

  readint (&(*buff), &(*p), &n);
  if ((n >= 1) && (n <= rvarlimit))
    {
      temp = getchrc (&(*fyle), &(*buff), &(*p));
      odisp (&vari.A[n - 1]);
      vari.A[n - 1] = temp;
      if (rjctindex > 0)
	if (vari.A[n - 1] == reject.A[rjctindex - 1])
	  {
	    reject.A[rjctindex - 1] = NULL;
	    rjctindex = rjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
getpset (text * fyle, string0 * buff, int *p)
{
  int n;
  prodtype temp;

  readint (&(*buff), &(*p), &n);
  if ((n >= 1) && (n <= pvarlimit))
    {
      temp = getprod (&(*fyle), &(*buff), &(*p));
      pdisp (&pvar.A[n - 1]);
      pvar.A[n - 1] = temp;
      if (prjctindex > 0)
	if (pvar.A[n - 1] == preject.A[prjctindex - 1])
	  {
	    preject.A[prjctindex - 1] = NULL;
	    prjctindex = prjctindex - 1;
	  }
    }
  else
    {
      warn ("illegal variable", cont);
      inform (":No action taken", cr);
    }
}

void
ndynk (text * fyle, string0 * buff, int *p, int tipe)
{
  ocharptr lam, mu;
  bframe dynkn;
  int norml = 0,		//modified by FB, was uninitialized
     n;
  frame fval;

  tobig (0, &dynkn);
  switch ((int) (ggroup.name))
    {
    case sn:
    case unm:
    case nill:
    case sunm:
    case un:
    case on:
    case ospnm:
    case spnc:
    case mp:
    case sonc:
    case an:
    case unc:;
      break;
    /* case un: */ case sung:
    case son:			/* case on: */
    case spn:
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
      readint (&(*buff), &(*p), &n);
      if (tipe == 1)
	readchrc (&(*fyle), &(*buff), &(*p), &mu);
      else
	{
	  adjoint (&fval);
	  cnu (&mu);
	  {
	    register ocharptr W5 = &(*mu);

	    W5->lab = ' ';
	    W5->C6_double = false;
	    W5->next = NULL;
	    W5->spin = false;
	    W5->mult = 1;
	    W5->val = fval;
	  }
	}
      lam = getchrc (&(*fyle), &(*buff), &(*p));
      ppn = n;
      ppt = tipe;
      trace = false;
      tipe = ppt;
      n = ppn;
      group = ggroup.name;
      {
	register groop *W6 = &ggroup;

	switch ((int) (W6->name))
	  {
	  case sung:		/* case un: */
	    norml = 2 * W6->rank * W6->rank;
	    break;
	  case spn:
	    norml = 4 * (W6->rank / 2 + 1);
	    break;
	  /* case on: */ case son:
	    norml = 8 * (W6->rank - 2);
	    break;
	  case g2:
	    norml = 12;
	    break;
	  case f4:
	    norml = 18;
	    break;
	  case e6:
	    norml = 144;
	    break;
	  case e7:
	    norml = 72;
	    break;
	  case e8:
	    norml = 60;
	    break;
	  default:
	    Caseerror (Line);
	  }
      }
      dynkinn (ggroup.rank, n, mu, lam, &dynkn);
      if (n >= 1)
	{
	  bigtofact (dynkn, &dynkn);
	  kdiv (dynkn, norml, &dynkn);
	  factobig (dynkn, &dynkn);
	}
      odisp (&mu);
      if (logging)
	Putchr ('\n', logfile);
      if (n <= 1)
	print (" trace %1d\n", n);
      else
	print ("%1d**%1d*trace %1d", norml, n - 1, n);
      break;
    default:
      Caseerror (Line);
    }
  /*error(GROUP_NOT_SET, (*p)); */
  switch ((int) (ggroup.name))
    {
    case sung:
    case son:
    case spn:
    case g2:
    case f4:
    case e6:
    case e7:
    case e8:
      inform ("th order Casimir", cont);
      inform (" invariant is ", cont);
      wrtbigno (&output, &(*p), &dynkn);
      Putchr ('\n', output);
      if (logging)
	{
	  wrtbigno (&logfile, &(*p), &dynkn);
	  Putchr ('\n', logfile);
	}
      break;
    case sn:
    case unm:
    case nill:
    case sunm:
    case un:
    case on:
    case ospnm:
    case spnc:
    case mp:
    case sonc:
    case an:
    case unc:;
      break;
    default:
      Caseerror (Line);
    }
}


bool
ibranch (text * fyle, string0 * buff, int *p)
{
  register bool R288;
  char word[MAXSTRING];
  prodtype pr, top;
  ocharptr temp, help, help1;
  int nprd, nn, n, m, n2, m2,  brno;
  groop group1;

  inform ("Branch Mode", cr);
  mode = brm;
  do
    {
      fflush (stdout);
#ifdef HAVE_LIBREADLINE
      strcpy (prompt, "enter branching & rule numbers> ");
#else
      prompt[0] = '\0';
#endif
      readacard (&(*fyle), &(*buff), &(*p));
      prompt[0] = '\0';
      readword (&(*buff), &(*p), word);
      if (!
	  (interp ("stop", word, 4)
	   || interp ("tableofbranchingrules", word, 3)
	   || interp ("exitmode", word, 4)
	   || interp ("end", word, 3) || interp ("quit", word, 4)))
	{
	  (*p) = 1;
	  readint (&(*buff), &(*p), &brno);
	  if (((brno < 1) || (brno > 63)))
	    print ("Inappropriate Branching Rule number\n");
	  else
	    {
	      nprod = 1;
	      nprd = nprod;
	      //tprod = nprod;
	      /*if (!(Member((unsigned)(brno), Conset[2])))        
	         readint(&(*buff), &(*p), &n); */
	      if (((brno < 40) || (brno > 59)))
		readint (&(*buff), &(*p), &n);
	      switch ((int) (brno))
		{
		case 4:
		case 5:
		case 8:
		case 11:
		case 14:
		case 15:
		case 16:
		case 22:
		case 23:
		case 24:
		case 27:
		case 32:
		case 33:
		case 34:
		case 35:
		case 37:
		case 38:
		case 39:
		  readint (&(*buff), &(*p), &m);
		  break;
		case 1:
		case 2:
		case 3:
		case 6:
		case 7:
		case 9:
		case 10:
		case 12:
		case 13:
		case 17:
		case 18:
		case 19:
		case 20:
		case 21:
		case 25:
		case 26:
		case 28:
		case 29:
		case 30:
		case 31:
		case 40:
		case 41:
		case 42:
		case 43:
		case 44:
		case 45:
		case 46:
		case 47:
		case 48:
		case 49:
		case 50:
		case 51:
		case 52:
		case 53:
		case 54:
		case 55:
		case 56:
		case 36:
		case 60:
		case 61:
		case 62:
		case 57:
		case 58:
		case 59:
		case 63:
		  break;
		default:
		  Caseerror (Line);
		}
	      if (((brno == 33) || (brno == 34)))
		{
		  readint (&(*buff), &(*p), &n2);
		  readint (&(*buff), &(*p), &m2);
		}

	      if (((brno > 0) && (brno <= 63)))
		switch ((int) (brno))
		  {
		  case 1:
		  case 2:
		  case 3:
		  case 6:
		  case 7:
		  case 62:
		    fixcg (&group1, un, n, 0);
		    break;
		  case 4:
		    fixcg (&group1, un, (n + m), 0);
		    break;
		  case 5:
		    fixcg (&group1, un, (n * m), 0);
		    break;
		  case 8:
		    fixcg (&group1, sung, (n + m), 0);
		    break;
		  case 9:
		  case 10:
		  case 11:
		  case 12:
		  case 13:
		    fixcg (&group1, spn, n, 0);
		    break;
		  case 14:
		    fixcg (&group1, spn, (n + m), 0);
		    break;
		  case 15:
		    fixcg (&group1, spn, (n * m), 0);
		    break;
		  case 16:
		    fixcg (&group1, sn, (n + m), 0);
		    break;
		  case 17:
		    fixcg (&group1, sn, n, 0);
		    break;
		  case 18:
		  case 19:
		  case 20:
		  case 21:
		  case 25:
		    fixcg (&group1, on, n, 0);
		    break;
		  case 22:
		    fixcg (&group1, on, (n + m), 0);
		    break;
		  case 23:
		    fixcg (&group1, on, (n * m), 0);
		    break;
		  case 24:
		    fixcg (&group1, on, (n * m), 0);
		    break;
		  case 26:
		    fixcg (&group1, son, n, 0);
		    break;
		  case 27:
		    fixcg (&group1, son, (n + m), 0);
		    break;
		  case 28:
		  case 29:
		    fixcg (&group1, son, 4, 0);
		    break;
		  case 30:
		  case 31:
		    fixcg (&group1, son, 7, 0);
		    break;
		  case 32:
		    fixcg (&group1, sunm, n, m);
		    break;
		  case 33:
		    fixcg (&group1, sunm, (n + n2), (m + m2));
		    break;
		  case 34:
		    fixcg (&group1, unm, (n * n2 + m * m2),
			   (n * m2 + n2 * m));
		    break;
		  case 35:
		    fixcg (&group1, ospnm, n, m);
		    break;
		  case 36:
		    fixcg (&group1, spnc, n, 0);
		    break;
		  case 37:
		    fixcg (&group1, spnc, n, 0);
		    break;
		  case 38:
		    fixcg (&group1, spnc, n * m, 0);
		    break;
		  case 39:
		    fixcg (&group1, mp, n * m, 0);
		    break;
		  case 40:
		  case 41:
		  case 42:
		    fixcg (&group1, g2, 2, 0);
		    break;
		  case 43:
		  case 57:
		    fixcg (&group1, f4, 4, 0);
		    break;
		  case 44:
		  case 45:
		  case 46:
		  case 56:
		  case 58:
		    fixcg (&group1, e6, 6, 0);
		    break;
		  case 47:
		  case 48:
		    fixcg (&group1, e7, 7, 0);
		    break;
		  case 60:
		  case 49:
		  case 50:
		  case 51:
		  case 52:
		    fixcg (&group1, e8, 8, 0);
		    break;
		  case 53:
		    fixcg (&group1, sung, 27, 0);
		    break;
		  case 54:
		    fixcg (&group1, sung, 56, 0);
		    break;
		  case 55:
		    fixcg (&group1, sung, 248, 0);
		    break;
		  case 61:
		    fixcg (&group1, sonc, n, 0);
		    break;
		  case 59:
		    fixcg (&group1, son, 26, 0);
		    break;
		  case 63:
		    fixcg (&group1, sn, 8, 0);
		    break;
		  default:
		    Caseerror (Line);
		  }

	      ggroup = currgrp.A[1 - 1];
	      currgrp.A[1 - 1] = group1;
	      //group2 = currgrp.A[1 - 1];
	      putgroup1 (currgrp, 2);
	      inform (" to ;", cont);
	      fixsubgroup (&brno, &n, &m, &n2, &m2);
	      putgroup1 (currgrp, 2);
	      echo = false;
	      nprd = nprod;
	      inform (" ;", cr);
	      /*if (Member((unsigned)(brno), Conset[6])) */
	      if (((brno >= 43) && (brno <= 60)))
		{
		  echo = false;
		  exceptgroupload (&brno);
		  echo = true;
		}
	      if (brno == 63)
		{
		  echo = false;
		  exceptgroupload (&brno);
		  echo = true;
		}
	      do
		{
		  //inform ("BRM>", cr);
		  readacard (&(*fyle), &(*buff), &(*p));
		  readword (&(*buff), &(*p), word);
		  if (!
		      (interp ("stop", word, 4)
		       || interp ("exitmode", word, 4)
		       || interp ("end", word, 3)
		       || interp ("quit", word, 4)
		       || interp ("tableofbranchingrules", word, 3)))
		    {
		      (*p) = 1;
		      pnu (&pr);
		      fixcg (&group1, group1.name, group1.rank, group1.rank2);
		      pr->prods.A[1 - 1] =
			getchrc (&(*fyle), &(*buff), &(*p));
		      temp = pr->prods.A[1 - 1];
		      pr->prods.A[1 - 1] = gmodify (temp, group1);
		      odisp (&temp);
		      if ((pr->prods.A[1 - 1] != NULL))
			{
			  pr->mult = 1;
			  top = pr;
			  {
			    register prodtype W7 = &(*pr);

			    help = W7->prods.A[1 - 1];
			    help1 = W7->prods.A[1 - 1];
			    dobranch (&help, &help1, brno, n, m, n2, m2);
			    W7->prods.A[1 - 1] = help1;
			  }
			  pr = top;
			  nn = nprod;
			  nprod = 1;
			  pr = prodexpand (top);
			  nprod = nn;
			  nprd = 1;
			  top = pr;
			  fixgroups (&pr, brno, nprd);
			  pdisp (&pr);
			  pr = top;
			  /*su1crunch(&pr); */
			  pr = pmodify (pr);
			  schur_psort (&pr, true);
			  if (orderlist)
			    pr = preverselist (pr);
			  putprod (pr);
			  pdisp (&pr);
			  odisp (&help1);
			  if ((brno == 27) && ((n == 1) || (m == 1)))
			    {
			      fixcg (&group1, son, (n + m), 0);
			      fixsubgroup (&brno, &n, &m, &n2, &m2);
			    }
			  else if ((brno == 10) && (n == 2))
			    {
			      fixcg (&group1, spn, n, 0);
			      fixsubgroup (&brno, &n, &m, &n2, &m2);
			    }	/*15/2/98 */
			}
		      else
			print ("Irrep modifies to null\n");
		    }
		}
	      while (!
		     ((interp ("stop", word, 4)
		       || interp ("exitmode", word, 4)
		       || interp ("end", word, 3)
		       || interp ("quit", word, 4))));
	    }
	}
      strcpy (instr, word);
    }
  while (!
	 (interp ("end", word, 3)
	  || interp ("exitmode", word, 4) || interp ("quit", word, 4)));
  if (interp ("exitmode", word, 4))
    {
      nprod = 0;
      discardedchrcs = NULL;
      for (srjctindex = 1; srjctindex <= srjctlimit; srjctindex++)
	{
	  sreject.A[srjctindex - 1] = NULL;
	}
      srjctindex = 0;
      current = NULL;
      discardedprds = NULL;
      for (rjctindex = 1; rjctindex <= rjctlimit; rjctindex++)
	{
	  reject.A[rjctindex - 1] = NULL;
	}
      rjctindex = 0;
      for (prjctindex = 1; prjctindex <= prjctlimit; prjctindex++)
	{
	  preject.A[prjctindex - 1] = NULL;
	}
      prjctindex = 0;
      ggroup.name = nill;
      ggroup.rank = 0;
      ggroup.rank2 = 0;
      current = NULL;
      pcurrent = NULL;
    }
  R288 = interp ("exitmode", word, 4);
  return R288;
}

void
enter (void)
{
  readword (&buffz, &kkz, instr);
  readint (&buffz, &kkz, &jz);
//  Putchr ('\n', output);  //FB 20070311
//  for (kkz = 1; kkz <= 50; kkz++)
//    {
//      Putchr (buffz.A[kkz - 1], output);
//    }
//  Putchr ('\n', output); 
  memcpy (prompt, buffz.A, bcol);	
  for (kkz = bcol-1; kkz >= 0 && prompt[kkz] == ' '; kkz--)
    ;
  prompt[kkz + 2] = '\0';
  kkz = bcol;

  readacard (&input, &buffz, &ppz);
  if (interp ("svar", instr, 2))
    getsset1 (&input, &buffz, &ppz, &jz);
  else if (interp ("rvar", instr, 2))
    getset1 (&input, &buffz, &ppz, &jz);
  else if (interp ("varfordpreps", instr, 1))
    getpset1 (&input, &buffz, &ppz, &jz);
  else if (interp ("logfile", instr, 3))
    logcom (buffz, &ppz);
  else if (interp ("group", instr, 1))
    {
      getgroup (&buffz, ppz, true);
      ggroup = currgrp.A[1 - 1];
    }
  prompt[0] = '\0';
}

void
sfnmode (void)
{
  termptr p1, p2, p3;
  bool test;
  int n, qq;
  bframe bbig;

  inform ("Schur Function ;", cont);
  inform ("Mode    ;", cr);
  mode = sfnm;
  qq = 12;
  do
    {
      erred = false;
      //stack = false;
/*      if (fnex && fnptr!=NULL)
	{
	  buffz = fnptr->gbuff;
	  fnptr = fnptr->next;
	  ppz = 1;
	  fnex = (bool) (fnptr != NULL);
	  if (!fnex)
	    iosup = false;
	  if (debug_schur && fnex)
	    fprintf (stderr, "fnex_sfnmode %s\n",buffz.A);
	}
      else
	{
		fnex = false;*/
	  //inform ("SFN>", cr);
	  if (bell)
	    Putchr (bel, output);
	  readacard (&input, &buffz, &ppz);
//	}
      kkz = ppz;
      readword (&buffz, &kkz, instr);
      if (interp ("setsvar", instr, 5)) 
      {
	getsset (&input, &buffz, &kkz);
      }
      else if (interp ("dead", instr, 4))
	{
	  p1 = getsfn (&input, &buffz, &kkz);
	  p2 = getsfn (&input, &buffz, &kkz);
	  p3 = getsfn (&input, &buffz, &kkz);
	  test = dead (p1->val, p2, p3);
	  print ("dead =%s\n", boolstr(test));
	}
      else if (interp ("hclass", instr, 3))
	{
	  p1 = getsfn (&input, &buffz, &kkz);
	  hclass (p1->val, &bbig);
	  print ("Number of elements in the class = \n");
	  wrtbigno (&output, &qq, &bbig);
	  Putchr ('\n', output);
	  if (logging)
	    {
	      wrtbigno (&logfile, &qq, &bbig);
	      Putchr ('\n', logfile);
	    }
	}
      else if (interp ("return", instr, 3))
	{
	  Putchr (cr, output), Putchr ('\n', output);
	  if (logging)
	    Putchr ('\n', logfile);
	}
      else if (interp ("maxcoeff", instr, 4))
	maxscoeff (getsfn (&input, &buffz, &kkz));
      else if (interp ("countcoeffsinlist", instr, 6))
	ssummult (getsfn (&input, &buffz, &kkz));
      else if (interp ("mult_splitintotwolists", instr, 7))
	multsplit (getsfn (&input, &buffz, &kkz));
      else if (interp ("counttermsinlist", instr, 6))
	tssum (getsfn (&input, &buffz, &kkz));
      else if (interp ("squares", instr, 2))
	sumabssquares (getsfn (&input, &buffz, &kkz));
      else if (interp ("sfnmode", instr, 3))
	fprintf (stderr, "You are already in SFN Mode\n");
      else if (common (sfnm, instr, buffz, &kkz))
	;
      else if (interp ("content", instr, 7))
	{
	  bool displaySum = false;
	  char word[MAXSTRING];
	  readint (&buffz, &kkz, &n);
	  p1 = getsfn (&input, &buffz, &kkz);
	  readword (&buffz, &kkz, word);
	  if (word != NULL && (UpperCase (word[0]) == 'T'))
	    displaySum = true;
	  htablx (p1, output.fp, n, false, true, displaySum);
	  if (logging)
	    htablx (p1, logfile.fp, n, false, false, displaySum);
	}
      else if (interp ("mult_list", instr, 6))
	{
	  multlist (getsfn (&input, &buffz, &kkz));
	  //srjctindex = srjctindex - 1;
	}
      else if (interp ("nlambda", instr, 2))
	{
	  p1 = getsfn (&input, &buffz, &kkz);
	  n = nlambda (p1->val);
	  print ("nlambda = %10d\n", n);
	}
      else if (interp ("kmatrix", instr, 2))
	{
	  readint (&buffz, &kkz, &n);
	  kostkamatrix (n);
	}
      else if (interp ("kostka", instr, 1))
	{
	  p1 = getsfn (&input, &buffz, &kkz);
	  p2 = getsfn (&input, &buffz, &kkz);
	  kostka (p1, p2, &bbig);
	  print ("Kostka matrix element = \n"),
	    wrtbigno (&output, &qq, &bbig);
	  Putchr ('\n', output);
	  if (logging)
	    {
	      wrtbigno (&output, &qq, &bbig);
	      Putchr ('\n', output);
	    }
	}
      else if (interp ("HIVESLRCoefficient", instr, 8))
	{
	  unsigned long res;
	  unsigned short l;
	  bool displayResults = false, secho;
	  // lambda mu nu dilatation
	  p1 = getsfn (&input, &buffz, &kkz);
	  p2 = getsfn (&input, &buffz, &kkz);
	  p3 = getsfn (&input, &buffz, &kkz);
	  if (p1->next != NULL || p2->next != NULL || p3->next != NULL)
	    print ("Erreur only simple partitions are allowed here");

	  l = MAX (MAX (p1->val.length, p2->val.length), p3->val.length);

	  p1->val.length = l;
	  p2->val.length = l;
	  p3->val.length = l;
	  readint (&buffz, &kkz, &n);
	  if (n == 0)
	    n = 1;

	  secho = echo;
	  echo = false;
	  boole (buffz, &kkz, &displayResults);

	  res = hivesLRcoef (p1->val, p2->val, p3->val, l, n, displayResults);
	  echo = secho;
	  print ("Littlewood Richardson Coefficient = %lu\n", res);
	  srjctindex = srjctindex - 1;
	}
      else if (interp ("paritysequence", instr, 6))
	{
	  int ppp;
	  p1 = getsfn (&input, &buffz, &kkz);
	  ppp = paritysequence (p1);
	  if (ppp)
	    print ("parity of sequence is even\n");
	  else
	    print ("parity of sequence is odd\n");
	  srjctindex = srjctindex - 1;
	}
      else if (interp ("yhooks", instr, 2))
	{
	  bool displaySum = false;
	  char word[MAXSTRING];

	  n = 0;
	  p1 = getsfn (&input, &buffz, &kkz);
	  readword (&buffz, &kkz, word);
	  if (word != NULL && (UpperCase (word[0]) == 'T'))
	    displaySum = true;

	  htablx (p1, output.fp, n, true, true, displaySum);
	  if (logging)
	    htablx (p1, logfile.fp, n, true, false, displaySum);
	}
      else if (interp ("Youngdiagrams", instr, 2))
	{
	  p1 = getsfn (&input, &buffz, &kkz);
	  dtablx (p1, false, output.fp, true);
	  if (logging)
	    dtablx (p1, false, logfile.fp, false);
	}
      else
	if (!
	    (interp ("exitmode", instr, 4)
	     || interp ("end", instr, 3)
	     || interp ("quit", instr, 4)
	     || interp ("dpmode", instr, 3)
	     || interp ("repmode", instr, 3)) && (buffz.A[1 - 1] != ' '))
	{
	  p1 = getsfn (&input, &buffz, &ppz);
	  putsfn (&output, p1, true);
	  /*if (logging)
	     putsfn (&logfile, p1, true); */
	}

      allfalse ();		// FB 20060402
      scollctgarbage ();
    }
  while (!
	 (interp ("exitmode", instr, 4)
	  || interp ("end", instr, 3)
	  || interp ("quit", instr, 4)
	  || interp ("dpmode", instr, 3) || interp ("repmode", instr, 3)));
} //sfnmode

void
repmode (void)
{
  termptr dynks;
  ocharptr chrc, chrc1, chrc2;
  bool test;
  inform ("REP mode       ;", cr);
  mode = repm;
  ggroup = currgrp.A[1 - 1];
  if (nprod > 0)
    {
      nprod = 1;
      if (echo)
	putgroup (currgrp);
    }
  do
    {
      erred = false;
      if (bell)
	Putchr (bel, output);
      readacard (&input, &buffz, &ppz);
      kkz = ppz;
      readword (&buffz, &kkz, instr);
      if (interp ("setrvar", instr, 5))
	getset (&input, &buffz, &kkz);
      else if (interp ("generic", instr, 7))
	{
	  chrc1 = getchrc (&input, &buffz, &kkz);
	  chrc2 = getchrc (&input, &buffz, &kkz);
	  if (chrc1 != NULL && chrc2 != NULL)
	    {
	      test = generic (chrc1, chrc2, currgrp.A[1 - 1]);
	      print ("generic=%s\n", boolstr(test));
	    }
	  else
	    error (MISSING_PARAMETER, 1);
	}
      else if (interp ("hstdlist", instr, 5))
	{
	  chrc1 = getchrc (&input, &buffz, &kkz);
	  test = hstd (chrc1, currgrp.A[1 - 1]);
	  print ("highlystandard list =%s\n", boolstr(test));
	}
      else if (interp ("hstd", instr, 4))
	{
	  chrc1 = getchrc (&input, &buffz, &kkz);

	  test = highlystandard (chrc1, currgrp.A[1 - 1]);
	  print ("highlystandard =%s\n", boolstr(test));
	}
      else if (interp ("group", instr, 1))
	{
	  getgroup (&buffz, kkz, true);
	  ggroup = currgrp.A[1 - 1];
	}
      else if (interp ("dimension", instr, 3))
	dimnprop (getchrc (&input, &buffz, &kkz));
      else if (interp ("dynkinindex", instr, 7))
	putdindex (getchrc (&input, &buffz, &kkz));
      else if (interp ("return", instr, 3))
	{
	  Putchr (cr, output), Putchr ('\n', output);
	  if (logging)
	    Putchr ('\n', logfile);
	}
      else if (interp ("dpmode", instr, 3))
	return;
      else if (interp ("d_to_plabel", instr, 6))
	{
	  readlist (&input, &buffz, &kkz, &dynks, false);
	  group = ggroup.name;
	  inform ("<<partition>> ;", cont);
	  chrc = partlbl (ggroup.rank, dynks);
	  putchrc (&output, chrc, true);
	  odisp (&chrc);
	  ldisp (&dynks);
	}
      else if (interp ("maxcoeff", instr, 4))
	maxrcoeff (getchrc (&input, &buffz, &kkz));
      else if (interp ("sfnmode", instr, 3))
	sfnmode ();
      else if (interp ("repmode", instr, 3))
	;
      else if (interp ("propertyofreplist", instr, 4))
	{
	  chrc = getchrc (&input, &buffz, &kkz);
	  chrc1 = gmodify (chrc, currgrp.A[1 - 1]);
	  putprop (chrc1, true);
	  odisp (&chrc1);
	}
      else if (interp ("schar", instr, 5))
	{
	  chrc1 = getchrc (&input, &buffz, &kkz);
	  chrc2 = getchrc (&input, &buffz, &kkz);
	  schar (chrc1, chrc2);
	}
      else if (interp ("p_to_dlabel", instr, 6))
	{
	  chrc = getchrc (&input, &buffz, &kkz);
	  chrc1 = gmodify (chrc, currgrp.A[1 - 1]);
	  putprop (chrc1, false);
	  odisp (&chrc1);
	}
      else if (interp ("splitintospinandtensor", instr, 3))
	spinsplit (getchrc (&input, &buffz, &kkz));
      else if (interp ("squares", instr, 2))
	sumabssquaresrep (getchrc (&input, &buffz, &kkz));
      else if (interp ("countcoeffsinlist", instr, 6))
	summult (getchrc (&input, &buffz, &kkz));
      else if (interp ("onscalar", instr, 4))
	{
	  chrc1 = getchrc (&input, &buffz, &kkz);
	  chrc2 = getchrc (&input, &buffz, &kkz);
	  onttscalar (chrc1, chrc2, currgrp.A[1 - 1]);
	}
      else if (interp ("counttermsinlist", instr, 6))
	tsum (getchrc (&input, &buffz, &kkz));
      else if (interp ("casimirgeneralnthtrace", instr, 8))
	ndynk (&input, &buffz, &kkz, 1);
      else if (interp ("casimirnthordertrace", instr, 3))
	ndynk (&input, &buffz, &kkz, 2);
      else if (interp ("whatgroup", instr, 5))
	putgroup (currgrp);
      else if (interp ("consplit", instr, 8))	// added by FB
	contragsplit (getchrc (&input, &buffz, &kkz));
      else if (common (repm, instr, buffz, &kkz))
	;
      else
	if (!(interp ("exitmode", instr, 4)
	      || interp ("end", instr, 3)
	      || interp ("quit", instr, 4)
	      || interp ("dpmode", instr, 3) || interp ("sfnmode", instr, 3)))
	putchrc (&output, getchrc (&input, &buffz, &ppz), true);
      collectgarbage ();
    }
  while (!(interp ("exitmode", instr, 4)
	   || interp ("end", instr, 3)
	   || interp ("quit", instr, 4)
	   || interp ("dpmode", instr, 3) || interp ("sfnmode", instr, 3)));
}

void
dpmode (void)
{
  int n;
  buffz = buffi;

  inform ("DPrep Mode (with", cont);
  inform (" function)", cr);
  mode = dpm;
  ppz = 1;			//added by FB
  if (nprod > 0)
    if (echo)
      putgroup (currgrp);
  do
    {
      erred = false;
      if (fnex)
	{
	  buffz = fnptr->gbuff;
	  fnptr = fnptr->next;
	  fnex = (bool) (fnptr != NULL);
	  if (!fnex)
	    iosup = false;
	  if (debug_schur && fnex)
	    fprintf (stderr, "fnex %s\n", buffz.A);
	  ppz = 1;
	}
      else
	{
	  //inform ("DP>;", cr);
	  if (bell)
	    Putchr (bel, output);
	  readacard (&input, &buffz, &ppz);
	}
      kkz = ppz;
      readword (&buffz, &kkz, instr);
      if (interp ("setrvar", instr, 5))
	getset (&input, &buffz, &kkz);
      else if (interp ("setvarindpmode", instr, 4))
	getpset (&input, &buffz, &kkz);
      else if (interp ("return", instr, 3))
	{
	  Putchr (cr, output), Putchr ('\n', output);
	  if (logging)
	    Putchr ('\n', logfile);
	}
      else if (interp ("brmode", instr, 3))
	{
	  
	  if (ibranch (&input, &buffz, &kkz))
	    {
	      inform ("DPrep Mode (with", cont);
	      inform (" function)", cr);
	    }
	  mode = dpm;
	}
      else if (interp ("dpmode", instr,3))
      {
	      fprintf(stderr,"You are already in DPrep mode\n");
      }
      else if (interp ("group", instr, 1))
	{
	  getgroup (&buffz, kkz, true);
	  ggroup = currgrp.A[1 - 1];
	}
      else if (interp ("whatgroup", instr, 5))
	putgroup (currgrp);
      else if (interp ("rvar", instr, 2))
	{
	  readint (&buffz, &kkz, &ppz);
	  putchrc (&output, vari.A[ppz - 1], true);
	}
      else if (interp ("dimension", instr, 2))
	pproperties (getprod (&input, &buffz, &kkz));
      else if (interp ("uonetrace", instr, 5))
	utrace (getprod (&input, &buffz, &kkz));
      else if (interp ("countcoeffsinlist", instr, 6))
	psummult (getprod (&input, &buffz, &kkz));
      else if (interp ("counttermsinlist", instr, 6))
	tpsum (getprod (&input, &buffz, &kkz));
      else if (interp ("maxcoeff", instr, 4))
	maxdcoeff (getprod (&input, &buffz, &kkz));

      else if (common (dpm, instr, buffz, &kkz))
	;
      else if (interp ("canceldatfile", instr, 3))
	{
	  readint (&buffz, &kkz, &n);
	  if (!erred)
	    cancel (n);
	}

      else
	if (!(interp ("end", instr, 3)
	      || interp ("quit", instr, 4)
	      || interp ("repmode", instr, 3)
	      || interp ("sfnmode", instr, 3)))
	putprod (getprod (&input, &buffz, &ppz));
      pcollectgarbage ();
      scollctgarbage ();
    }
  while (!(interp ("end", instr, 3)
	   || interp ("quit", instr, 4)
	   || interp ("repmode", instr, 3) || interp ("sfnmode", instr, 3)));
}
