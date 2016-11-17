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

/** \file g.c
 */

/*
**	Definitions for i/o
*/
# include <stdio.h>
# include <string.h>
# include "standard.h"
# include "define.h"
# include "sets_mgmt.h"

/*
**	Start of program definitions
*/
# include "dim.h"
# include "type.h"

# include "var.h"
prodtype pcurrent;

# include "utils.h"
# include "s1.h"
# include "s2.h"
# include "s3.h"
# include "s4.h"
# include "s5.h"
# include "s6.h"
# include "s7.h"
# include "m.h"
# include "r.h"
# include "s.h"
# include "gr.h"
# include "q1.h"
# include "q2.h"
# include "bignums.h"
# include "branch.h"
# include "skew.h"
# include "label.h"
# include "rib_to_s.h"
# include "g.h"

/** analyse sfn expression recursively*/
termptr
getsfn (text * fyle, string0 * buffg, int *p)
{
  termptr sfn, sfn1, sfn2, sfn3;
  int n, ppg, qq, s, kk, m;
  register int ir;
  char ch;
  char word[MAXSTRING], word1[MAXSTRING];
  bool stack, series, ppp, song;
  lbframe polyg;

  qq = 20;
  ppg = (*p);
  stack = true;
  readword (&(*buffg), &(*p), word);
  sfn = NULL;
  qtest = false;
  if (interp ("absolutevalue", word, 3))
    {
      sfn = absval (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("rsamewtsfns", word, 5))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = rsameweight (sfn1->val);
    }
  else if (interp ("qsamewtsfns", word, 5))
    {
      allfalse ();
      qfn = true;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = qsameweight (&sfn1, &sfn2);
    }
  else if (interp ("wgen", word, 4))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn3 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = wgenerate (sfn1, sfn2, sfn3);
    }
  else if (interp ("add", word, 3))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = ladd (sfn1, sfn2);
    }
  else if (interp ("classlist", word, 5))
    {
      allfalse ();
      sfn = classlist (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("character", word, 4))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = character (sfn1, sfn2);
    }
  else if (interp ("snredpleth", word, 5))
    {
      allfalse ();
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = snredpleth (&sfn);
    }
  else if (interp ("snchar", word, 2))
    {
      allfalse ();
      readint (&(*buffg), &(*p), &n);
      sfn = (termptr) snchar (getsfn (&(*fyle), &(*buffg), &(*p)), n);
    }
  else if (interp ("complement", word, 5))
    {
      int m;
      allfalse ();
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      sfn = complement (getsfn (&(*fyle), &(*buffg), &(*p)), n, m);
    }
  else if (interp ("zraise", word, 2))
    {
      readint (&(*buffg), &(*p), &kk);
      readint (&(*buffg), &(*p), &qq);
      s = 1;
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      zraise (sfn, kk, qq, s);
      srjctindex = srjctindex - 1;
    }
 /* else if (interp ("paritysequence", word, 6))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      ppp = paritysequence (sfn);
      if (ppp)
	print ("parity of sequence is even\n");
      else
	print ("parity of sequence is odd\n");
      srjctindex = srjctindex - 1;
    }*/
  else if (interp ("wseq", word, 4))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      kk = wtfrm (&sfn->val);
      print ( "weight = %10d\n", kk);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("lseq", word, 4))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      kk = len (& sfn->val);
      print ( "length = %10d\n", kk);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("index", word, 5))
    sfn = (termptr) indexx (getsfn (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("allskewsfn", word, 3))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = lskcmpt (sfn1);
      chkskcmpt (sfn1, &sfn);
      allfalse ();
    }
  else if (interp ("attachpartitiontosfn", word, 2))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = attach (sfn1, sfn2, true);
    }
  else if (interp ("vmult", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = multmono (&sfn1, &sfn2, n);
    }
  else if (interp ("onmod", word, 5))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &kk);
      sfn = osnmod (getsfn (&(*fyle), &(*buffg), &(*p)), n, kk);
      Putchr (sfn->slab, output), Putchr ('\n', output);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("signseq", word, 4))
    {
      char word2[MAXSTRING];
      readint (&(*buffg), &(*p), &kk);
      readint (&(*buffg), &(*p), &n);
      readword (&(*buffg), &(*p), word2);
      readword (&(*buffg), &(*p), word1);
      ch = word1[0];
      if (((ch == 'c') || (ch == 'a')))
	song = false;
      else
	song = true;

      if (((ch == 't') || (ch == 'f')))
	{
	  if (ch == 't')
	    ppp = true;
	  else
	    ppp = false;

	  sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
	  sfn = signseq (&sfn1, kk, n, word2[0], ppp, song);
	  if (sfn != NULL)
	    srjctindex = srjctindex - 1;
	}
      else
	{
	  error (MISTAKE, (*p));
	  sfn = NULL;
	}
    }
  else if (interp ("genprod", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = genprod (&n);
    }
  else if (interp ("gwt", word, 3))	// added by FB to conform help file GWT
    {
      readint (&(*buffg), &(*p), &n);
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn1 = sfncopy (sfn2);
      schur_restrict (&sfn2, n, 'w');
      srjctindex = srjctindex - 1;
      sfn3 = sfnmult (-1, sfn2);
      sfn = ladd (sfn1, sfn3);
    }
  else if (interp ("smon", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = sfnmon (&sfn1, &sfn2, n);
    }
  else if (interp ("o_restrict", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = louter2 (sfn1, sfn2, n);
    }
  else if (interp ("nstdise", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      nstndise (&sfn, n);
      sort (&sfn, false);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("ch_coeffstooneforsfns", word, 4))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      smult (sfn);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("ch_phaseofsfns", word, 4))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      coeffset (&sfn, 1, '^');
      stack = false;
    }
  else if (interp ("compare", word, 4))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = compare (sfn2, sfn1);
    }
  else if (interp ("expandsfnlist", word, 3))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      cleave (&sfn);
      stack = false;
    }
  else if (interp ("conjadd", word, 7))
    {
      sfn = conjadd (getsfn (&(*fyle), &(*buffg), &(*p)));
      srjctindex = srjctindex - 1;
    }
  else if (interp ("conjugatesfnlist", word, 4))
    {
      sfn = lconjgte (getsfn (&(*fyle), &(*buffg), &(*p)));
      srjctindex = srjctindex - 1;
    }
  else if (interp ("conv_d_to_s", word, 11))
    sfn = prodsfn (getprod (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("conv_r_to_s", word, 11))
    sfn = repsfn (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("rm_partsequaln", word, 8))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      limit (&sfn, n);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_oddparts", word, 8))
    {
      n = -1;
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      limit (&sfn, n);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_evenparts", word, 8))
    {
      n = -2;
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      limit (&sfn, n);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_nmparts", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      partssfn (&sfn, n, m);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_partitionfromsfn", word, 4))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = attach (sfn1, sfn2, false);
    }
  else if (interp ("rm_repeatedpartssfns", word, 4))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      n = 0;
      distinct (&sfn);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("e_to_hsymmfn", word, 6))
    {
      sfn = elemtohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      homo = true;
    }
  else if (interp ("e_to_fsymmfn", word, 6))
    {
      sfn = homotomono (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      forg = true;
    }
  else if (interp ("e_to_msymmfn", word, 6))
    {
      sfn = elemtomono (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      mmono = true;
    }
  else if (interp ("e_to_ssymmfn", word, 6))
    {
      allfalse ();
      sfn = elemtosfn (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("f_to_hsymmfn", word, 6))
    {
      sfn = monotoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      homo = true;
    }
  else if (interp ("f_to_esymmfn", word, 6))
    {
      sfn = monotohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      elem = true;
    }
  else if (interp ("f_to_msymmfn", word, 6))
    {
      sfn1 = monotohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      sfn = elemtomono (sfn1);
      allfalse ();
      mmono = true;
    }
  else if (interp ("firstpart", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, n, 'f');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("frobenius", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, n, 'x');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("m_to_fsymmfn", word, 6))
    {
      sfn1 = monotoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      sfn = homotomono (sfn1);
      allfalse ();
      forg = true;
    }
  else if (interp ("s_to_fsymmfn", word, 6))
    {
      sfn1 = sfntohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      sfn = elemtomono (sfn1);
      allfalse ();
      forg = true;
    }
  else if (interp ("f_to_ssymmfn", word, 6))
    {
      sfn1 = monotoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      sfn = homotosfn (sfn1);
      allfalse ();
    }
  else if (interp ("rm_oddrksfnsonly", word, 7))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, 1, 'p');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_oddwtinlist", word, 7))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, 1, 'e');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("exitmode", word, 4))
    sfn = NULL;
  else if (interp ("hallpolygnomialproduct", word, 4))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      hallp (sfn1->val, sfn2->val);
      sfn = NULL;
      qtest = true;
    }
  else if (interp ("h_to_esymmfn", word, 6))
    {
      sfn = homotoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      elem = true;
    }
  else if (interp ("h_to_fsymmfn", word, 6))
    {
      sfn = elemtomono (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      forg = true;
    }
  else if (interp ("h_to_msymmfn", word, 6))
    {
      sfn = homotomono (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      mmono = true;
    }
  else if (interp ("h_to_ssymmfn", word, 6))
    {
      sfn = homotosfn (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
    }
  else if (interp ("i_sfnqfnproduct", word, 6))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      qfn = true;
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = sqinner (sfn1, sfn2);
    }
  else if (interp ("i_plethysmrd", word, 4))
    {
      sfn = plethonerinner (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      nreduce = true;
    }
  else if (interp ("i_qfnproduct", word, 3))
    {
      allfalse ();
      qfn = true;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = lqinner (sfn1, sfn2);
      qfn = false;
    }
  else if (interp ("insertpartitionintosfn", word, 3))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = insert (sfn1, sfn2);
    }
  else if (interp ("intdividecoeffs", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      if (n == 0)
	print ("division by zero not allowed\n");
      else
	coeffset (&sfn, n, '\\');
      /*stack = false; */
      srjctindex = srjctindex - 1;	/*26/8/98 */
    }
  else if (interp ("i_sfnproduct", word, 1))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      allfalse ();
      sfn = linner (sfn1, sfn2);
    }
  else if (interp ("latticetest", word, 3))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      ppp = platticetest (sfn);
      print ("lattice test = %s\n", boolstr (ppp));
      srjctindex = srjctindex - 1;
      allfalse ();
    }
  else if (interp ("lastresult", word, 4))
    sfn = sfncopy (scurrent);
  else if (interp ("lengthofpartitionsselect", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, n, 'l');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("makewtofsfnton", word, 4))
    {
      nreduce = false;
      readint (&(*buffg), &(*p), &n);
      sfn = makeweight (n, getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("m_timessfnproduct", word, 4))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      allfalse ();
      sfn = gordan (sfn1, sfn2);
    }
  else if (interp ("m_to_esymmfn", word, 6))
    {
      sfn = monotoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      elem = true;
    }
  else if (interp ("m_to_hsymmfn", word, 6))
    {
      sfn = monotohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      homo = true;
    }
  else if (interp ("m_to_sfnsymmfn", word, 6))
    {
      sfn = monotosfn (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
    }
  else if (interp ("mult_partsbyanint", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      multpart (sfn, n);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("mult_ntimes", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = ntensor (sfn1, n);
      allfalse ();
    }
  else if (interp ("nskew", word, 3))
  {
    readint(&(*buffg), &(*p), &n);
   sfn1 = getsfn(&(*fyle), &(*buffg), &(*p));
        sfn2 = getsfn(&(*fyle), &(*buffg), &(*p));
        sfn = nskew(sfn1,sfn2,n);
        //stack = false;
        allfalse();
  }
  else if (interp ("rm_evenrksfnsonly", word, 8))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      n = 0;
      schur_restrict (&sfn, 1, 'q');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("rm_evenwtinlist", word, 8))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      n = 0;
      schur_restrict (&sfn, 1, 'o');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("o_pfnproduct", word, 3))
    {
      allfalse ();
      qfn = true;
      pfn = true;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      qstndise (&sfn1);
      qstndise (&sfn2);
      redu = false;
      sfn = lqouter (&sfn1, &sfn2);
    }
  else if (interp ("o_qfnproduct", word, 3))
    {
      allfalse ();
      qfn = true;
      pfn = false;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      qstndise (&sfn1);
      qstndise (&sfn2);
      redu = false;
      sfn = lqouter (&sfn1, &sfn2);
      if (sfn1 == NULL) srjctindex = srjctindex - 1; // added by FB 20060810
      if (sfn2 == NULL) srjctindex = srjctindex - 1; // added by FB 20060810
    }
  else if (interp ("o_sfnproduct", word, 1))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      stndise (&sfn1);          // added by FB 20060809 
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      stndise (&sfn2);          // added by FB 20060809 
      sfn = louter (sfn1, sfn2);
      if (sfn1 == NULL) srjctindex = srjctindex - 1;// added by FB 20060810
      if (sfn2 == NULL) srjctindex = srjctindex - 1;// added by FB 20060810
      allfalse ();
    }
  else if (interp ("plethysm", word, 2))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = listplethlist (sfn1, sfn2);
    }
  else if (interp ("p_to_ssymmfn", word, 6))
    {
      allfalse ();
      sfn = powersumtosfn (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("s_to_psymmfn", word, 6))
    {
      allfalse ();
      sfn = s_to_p (getsfn (&(*fyle), &(*buffg), &(*p)));
      psum = true;
    }
  else if (interp ("qexpandspecialseries", word, 4))
    {
      allfalse ();
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      qsfnlist (sfn);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("qqexpandspecialseries", word, 5))
    {
      allfalse ();
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      qqsfnlist (sfn);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("qqseries", word, 4))
    {
      allfalse ();
      readword (&(*buffg), &(*p), word1);
      series = (bool) ((!(Member ((unsigned) (word1[0]),
				  Union (numbers.S, plus_minus_etc.S))))
		       && (word1[1] == ' '));
      if (series && UpperCase(word1[0])<'R' && UpperCase(word1[0])!='I' && UpperCase(word1[0])!='J') {
        Claimset ();
	qqallseries (word1[0], &polyg);
      }
      else
	print ("qqseries %c is not available\n", word1[0]);
      qtest = true;
    }
  else if (interp ("qseries", word, 4))
    {
      allfalse ();
      readword (&(*buffg), &(*p), word1);
      series = (bool) ((!  (Member ((unsigned) (word1[0]),
		   Union (numbers.S, plus_minus_etc.S))))
		&& (word1[1] == ' '));
      if (series && UpperCase(word1[0])<'R' && UpperCase(word1[0])!='I' && UpperCase(word1[0])!='J')
      {
        Claimset ();
	qallseries (word1[0], &polyg);
      }
      else
	print ( "qseries %c is not available\n", word1[0]);
      qtest = true;
    }
  else if (interp ("raiseinverseop", word, 6))
    {
      allfalse ();
      readint (&(*buffg), &(*p), &s);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = lraise (sfn, -s);
    }
  else if (interp ("raiseop", word, 5))
    {
      allfalse ();
      readint (&(*buffg), &(*p), &s);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = craise (sfn, s);
    }
  else if (interp ("rib_to_s", word, 8))
  {
	  allfalse ();
	  sfn = getsfn (&(*fyle), &(*buffg), &(*p));
	  if (sfn->next != NULL)
		  fprintf(stderr,"Incorrect input \n");
	  else
		  sfn = rib_to_s(sfn);
  }
  else if (interp ("rd_raiseinverseop", word, 9))
    {
      redu = true;
      allfalse ();
      nreduce = true;
      readint (&(*buffg), &(*p), &s);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = lraise (sfn, s);
      redu = false;
    }
  else if (interp ("rd_raiseop", word, 4))
    {
      redu = true;
      allfalse ();
      nreduce = true;
      readint (&(*buffg), &(*p), &s);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = craise (sfn, s);
      redu = false;
    }
  else if (interp ("rd_i_qfnproduct", word, 6))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      allfalse ();
      nreduce = true;
      qfn = true;
      sfn = rqinner (sfn1, sfn2);
    }
  else if (interp ("rd_i_sfnproduct", word, 4))
    {
      allfalse ();
      nreduce = true;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = rinner (sfn1, sfn2);
    }
  else if (interp ("rm_firstpartofsfn", word, 4))
    {
      allfalse ();
      nreduce = true;
      sfn = reducedq (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("rp_reporsfnbywt", word, 4))
    sfn = swt (getsfn (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("rp_sfncoeffbyint", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      coeffset (&sfn, n, '=');
      stack = false;
    }
  else if (interp ("riemannlist", word, 8))
    {
      allfalse ();
      readint (&(*buffg), &(*p), &n);
      if (n > 0)
	sfn = fullx (n);
      else
	sfn = NULL;
    }
  else if (interp ("riemannplethlist", word, 8))
    {
      allfalse ();
      readint (&(*buffg), &(*p), &n);
      if (n > 0)
	sfn = fullsa (n);
      else
	sfn = NULL;
    }
  else if (interp ("riemannscalarsordern", word, 8))
    {
      allfalse ();
      sfn = fulling (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("samewtsfns", word, 4))
    {
      allfalse ();
      sfn = leqwt (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("s_to_esymmfn", word, 6))
    {
      sfn = sfntoelem (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      elem = true;
    }
  else if (interp ("s_to_hsymmfn", word, 6))
    {
      sfn = sfntohomo (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      homo = true;
    }
  else if (interp ("s_to_msymmfn", word, 6))
    {
      sfn = lsfntomono (getsfn (&(*fyle), &(*buffg), &(*p)));
      allfalse ();
      mmono = true;
      srjctindex = srjctindex - 1;
    }
  else if (interp ("s_to_qsymmfn", word, 6))
    {
      sfn = s_to_q (getsfn (&(*fyle), &(*buffg), &(*p)));
      qfn = true;
      qstndise (&sfn);
      sort (&sfn, false);	/*7/6/99 */
      allfalse ();
      qfn = true;
    }
  else if (interp ("sumsquares", word, 5))
    {
      sfn = sumsquares (getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("q_to_sdual", word, 7))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      qstndise (&sfn1);
      sfn=q_to_s(sfn1);
      allfalse ();
      ldisp(&sfn1);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("sk_pfn", word, 4))
    {
      allfalse ();
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      n = (*p);
      readword (&(*buffg), &n, word1);
      series =
	(bool) ((!
		 (Member
		  ((unsigned) (word1[0]),
		   Union (numbers.S, plus_minus_etc.S))))
		&& (word1[1] == ' '));
      Claimset ();
      if (series)
	{
	  sfn2 = useseries (locase (word1[0]), sfn1, false, false, n);
	  (*p) = n;
	}
      else
	sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      redu = true;
      sfn = lqouter (&sfn1, &sfn2);
      pfn = true;
      qfn = true;
      if (series)
	ldisp (&sfn2);
    }
  else if (interp ("sk_qfn", word, 4))
    {
      allfalse ();
      qfn = true;
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      n = (*p);
      readword (&(*buffg), &n, word1);
      series =
	(bool) ((!
		 (Member
		  ((unsigned) (word1[0]),
		   Union (numbers.S, plus_minus_etc.S))))
		&& (word1[1] == ' '));
      Claimset ();
      if (series)
	{
	  sfn2 = useseries (locase (word1[0]), sfn1, false, false, n);
	  (*p) = n;
	}
      else
	sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      redu = true;
      sfn = lqouter (&sfn1, &sfn2);
      if (series)
	ldisp (&sfn2);
    }
  else if (interp ("sk_sfn", word, 2))
    {
      allfalse ();
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      n = (*p);
      readword (&(*buffg), &n, word1);
      series =
	(bool) ((!
		 (Member
		  ((unsigned) (word1[0]),
		   Union (numbers.S, plus_minus_etc.S))))
		&& (word1[1] == ' '));
      Claimset ();
      if (series)
	{
	  sfn = useseries (locase (word1[0]), sfn, false, true, n);
	  (*p) = n;
	}
      else
	sfn = lskew (sfn, getsfn (&(*fyle), &(*buffg), &(*p)));
    }
  else if (interp ("std_q", word, 5))
    {
      qfn = true;
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      qstndise (&sfn);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("std", word, 3))
    {
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      allfalse ();
      stndise (&sfn);
      sort (&sfn, false);
      srjctindex = srjctindex - 1;
    }
  else if (interp ("subtract", word, 3))
    {
      sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn2 = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn3 = sfnmult (-1, sfn2);
      sfn = ladd (sfn1, sfn3);
      ldisp (&sfn3);
    }
  else if (interp ("svar", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      //stack=false;
      if ((n >= 1) && (n <= svarlimit)) {
	sfn = sfncopy (svar.A[n - 1]);
      }
      else
	{
	  warn ("illegal variable", cr);
	  sfn = NULL;
	}
    }
  else
    if (interp
	("seriestermsthatskew", word, 8) || interp ("seriestointwt", word, 3))
    {
      allfalse ();
      if (!interp ("seriestermsthatskew", word, 8))	// modified by FB was interp ("seriestointwt", word, 3)+!
	{
	  readint (&(*buffg), &(*p), &n);
	  if (n >= 0)
	    {
	      snu (&sfn1);
	      sfn1->mult = 1;
	      for (ir = 1; ir <= maxdim; ir++)
		{
		  sfn1->val.A[ir] = 0;
		}
	      sfn1->val.A[1] = n;
	    }
	  else
	    sfn1 = NULL;
	}
      else
	sfn1 = getsfn (&(*fyle), &(*buffg), &(*p));
      readword (&(*buffg), &(*p), word1);
      if (validser (word1[0]) && (word1[1] == ' '))
	{
	  sfn = useseries (locase (word1[0]), sfn1, !interp ("seriestermsthatskew", word, 8), false, (*p));	// modified by FB was interp ("seriestointwt", word, 3)+!
	  if (!interp ("seriestermsthatskew", word, 8))	// modified by FB was interp ("seriestointwt", word, 3)+!
	    dispsfn (&sfn1);
	}
      else
	inform ("Bad series.", cr);
    }
  else if (interp ("wt", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      schur_restrict (&sfn, n, 'w');
      srjctindex = srjctindex - 1;
    }
  else if (interp ("mult_select", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      readword(&(*buffg), &(*p), word1);
      if (word1[0]!=' ' && (word1[0]=='t'|| word1[0]=='T'))
	restrict2 (&sfn, n, 'm', true);
      else {
        restrict2 (&sfn, n, 'm', false);
	}
      srjctindex = srjctindex - 1;
    }
  else if (interp ("mult_coeffsbyanint", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      sfn = getsfn (&(*fyle), &(*buffg), &(*p));
      sfn = sfnmult (n, sfn);
    }
  else if (interp ("yshapeselect", word, 2))
    {
      readword (&(*buffg), &(*p), word1);
      if (((word1[0] == 'c') || (word1[0] == 's')
	    || (word1[0] == 'd') || (word1[0] == 'r')
	    || (word1[0] == 't')) && (word1[1] == ' '))
	{
	  sfn = getsfn (&(*fyle), &(*buffg), &(*p));
	  n = -1;
	  schur_restrict (&sfn, n, word1[0]);
	  srjctindex = srjctindex - 1;
	}
      else
	{
	  sfn = NULL;
	  print ( "ERROR:invalid letter %c%.24s%c entered\n", qt, word, qt);
	}
    }
  else if (interp ("zero", word, 1))
    sfn = NULL;
  else
    {
      (*p) = ppg;
      while ((Member ((unsigned) (buffg->A[(*p) - 1]), brackets.S))
	     && ((*p) < bcol))
	(*p) = (*p) + 1;
      if (buffg->A[(*p) - 1] == cont)
	{
	  readacard (&(*fyle), &(*buffg), &(*p));
	  sfn = getsfn (&(*fyle), &(*buffg), &(*p));
	  stack = false;
	}
      else
	if (((buffg->A[(*p) - 1] == '0') || (buffg->A[(*p) - 1] == '1')
	     || (buffg->A[(*p) - 1] == '2') || (buffg->A[(*p) - 1] == '3')
	     || (buffg->A[(*p) - 1] == '4') || (buffg->A[(*p) - 1] == '5')
	     || (buffg->A[(*p) - 1] == '6') || (buffg->A[(*p) - 1] == '7')
	     || (buffg->A[(*p) - 1] == '8') || (buffg->A[(*p) - 1] == '9')
	     || (buffg->A[(*p) - 1] == '!') || (buffg->A[(*p) - 1] == '~')
	     || (buffg->A[(*p) - 1] == '-') || (buffg->A[(*p) - 1] == '+')))
	readlist (&(*fyle), &(*buffg), &(*p), &sfn, false);
	else //{
	  //srjctindex = srjctindex - 1;
	  error (MISTAKE, *p);
	//}
    }
  if (orderlist)
    sfn = reverselist (sfn);
  if (stack)
    sstkhandle (sfn);
  return sfn;
} // getsfn


void
allfalse (void)
{
  homo = false;
  elem = false;
  mmono = false;
  qfn = false;
  pfn = false;
  nreduce = false;
  forg = false;
  psum = false;
}

ocharptr
getchrc (text * fyle, string0 * buffg, int *p)
{
  ocharptr chrc, chrc1, chrc2, chrc3;
  int n, ppg, iauto, m;
  groop fg, tg;
  char word[MAXSTRING];
  bool stack;

  ppg = (*p);
  stack = true;
  readword (&(*buffg), &(*p), word);
  chrc = NULL;
  if (interp ("lastresult", word, 4))
    chrc = gchrccopy (current);
  else if (interp ("kinsert", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = kinsert (getchrc (&(*fyle), &(*buffg), &(*p)), n, currgrp.A[0]);
    }
  else if (interp ("plg", word, 3))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = g2p (chrc1);
      rjctindex--;        // FB 20060323
    }
  else if (interp ("associate", word, 5))
    chrc = associate (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
  else if (interp ("sponmodify", word, 5))
    chrc = sponmodify (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
  else if (interp ("sprextend", word, 5))
    chrc = sprextend (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
  else if (interp ("spstar", word, 6))
    chrc = star (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
  else if (interp ("compare", word, 4))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = rcompare (chrc2, chrc1);
    }
  else if (interp ("attachpartitiontosfn", word, 2))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = rattach (chrc1, chrc2, true);
    }
  else if (interp ("plethysmforgroup", word, 2))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = gpleth (chrc1, chrc2, currgrp.A[0]);
    }
  else if (interp ("product", word, 1))
    {
      chrc1 = getchrc (fyle, buffg, p);
      chrc2 = getchrc (fyle, buffg, p);
      chrc = kronk (chrc1, chrc2, currgrp.A[0]);
    }
  else if (interp ("signseq", word, 7))
    {
      readint (&(*buffg), &(*p), &n);
      m = 0;
      if (currgrp.A[0].name == un)
	readint (&(*buffg), &(*p), &m);
      chrc = signseqgr (getchrc (&(*fyle), &(*buffg), &(*p)), n, m);
    }
  else if (interp ("mult_select", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, n, 'm');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("mult_ntimes", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = rntensor (chrc1, n);
    }
  else if (interp ("firstpart", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, n, 'f');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("fprod", word, 5))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = fprod (chrc1, chrc2, currgrp.A[0]);
    }
  else if (interp ("ffprod", word, 6))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = ffprod (chrc1, chrc2, currgrp.A[0]);
    }
  else if (interp ("rm_oddwtinlist", word, 7))
    {
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, 1, 'e');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("rm_evenwtinlist", word, 8))
    {
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, 1, 'o');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("rm_partsequaln", word, 8))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      unlimit (&chrc, n);
      rjctindex = rjctindex - 1;
    }
  else if (interp ("rm_nmparts", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      partsrep (&chrc, n, m);
      rjctindex = rjctindex - 1;
    }
  else if (interp ("ch_labelforon", word, 4))
    chrc = clabel (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("racahnotation", word, 5))
    {
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      racahg2 (&chrc);
      //stack = false;
      rjctindex = rjctindex - 1;      // added by FB 20060321
    }
  else if (interp ("fracahnotation", word, 6))
    {
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      fracahg2 (&chrc);
      stack = false;
    }
  else if (interp ("rp_firstpartbyspin", word, 4))
    chrc = qspin (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("rm_sonevenlabels", word, 5))
    {
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      dlabel (chrc);
      stack = false;
    }
  else if (interp ("rp_reporsfnbywt", word, 4))
    chrc = rwt (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("scalarinner", word, 7))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      chrc = scalarinner (n, m);
    }
  else if (interp ("ch_spinindex", word, 4))
    chrc = cspin (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("intdividecoeffs", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rdiv (&chrc, n);
      stack = false;
    }
  else if (interp ("contragradientrep", word, 7))
    chrc = contrag (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
  else if (interp ("ch_uonereps", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = muc (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0], n);
    }
  else if (interp ("std_reps", word, 3))
    {
      if (currgrp.A[0].name == unc)
	readint (&(*buffg), &(*p), &kkz);
      chrc = gmodify (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
    }
  else if (interp ("sprch", word, 5))
    chrc = sprch (getchrc (&(*fyle), &(*buffg), &(*p)));
  else if (interp ("fusion", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = fusion (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0], n);
    }
  else if (interp ("conv_s_to_r", word, 11))
    {
      chrc = gettrans (&(*fyle), &(*buffg), &(*p));
      scollctgarbage ();
    }
  else if (interp ("rm_uonewtovermax", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = umax (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0], n);
    }
  else if (interp ("conv_d_to_r", word, 11))
    {
      chrc = prodcon (getprod (&(*fyle), &(*buffg), &(*p)));
      pcollectgarbage ();
    }
  else if (interp ("mult_coeffsbyanint", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = chrcmult (n, chrc);
    }
  else if (interp ("autoorisomorphism", word, 2))
    {
      if (readauto (&(*buffg), &(*p), &iauto, &fg, &tg, 1))
	{
	  chrc = change (getchrc (&(*fyle), &(*buffg), &(*p)), iauto, fg, tg);
	  putgroup (currgrp);
	}
    }
  else if (interp ("lengthofpartitionsselect", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, n, 'l');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("wt", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = getchrc (&(*fyle), &(*buffg), &(*p));
      rrestrict (&chrc, n, 'w');
      rjctindex = rjctindex - 1;
    }
  else if (interp ("add", word, 3))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc = chrcadd (chrc1, chrc2);
    }
  else if (interp ("subtract", word, 3))
    {
      chrc1 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc2 = getchrc (&(*fyle), &(*buffg), &(*p));
      chrc3 = chrcmult (-1, chrc2);
      chrc = chrcadd (chrc1, chrc3);
      odisp (&chrc3);
    }
  else if (interp ("covariant", word, 3))
    {
      chrc = covariant (getchrc (&(*fyle), &(*buffg), &(*p)), currgrp.A[0]);
    }
  else if (interp ("uonedivint", word, 5))
    {
      readint (&(*buffg), &(*p), &n);
      chrc = (ocharptr) uadd (getchrc (&(*fyle), &(*buffg), &(*p)), n, 2);
    }
  else if (interp ("zero", word, 1))
    chrc = NULL;
  else if (interp ("rvariable", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      if (currgrp.A[0].name==nill) { // added by FB20060303
	error (GROUP_NOT_SET,0);
      }
      else {
      if ((n >= 1) && (n <= rvarlimit))
	chrc = gchrccopy (vari.A[n - 1]);
      else
	{
	  warn ("illegal variable", cr);
	  chrc = NULL;
	}
      }
    }
  else
    {
      (*p) = ppg;
      if (word[0] == cont)
	{
	  readacard (&(*fyle), &(*buffg), &(*p));
	  chrc = getchrc (&(*fyle), &(*buffg), &(*p));
	  stack = false;
	}
      else
	{
	  while ((Member ((unsigned) (buffg->A[(*p) - 1]), brackets.S))
		 && ((*p) < bcol))
	    (*p) = (*p) + 1;
	  if (((buffg->A[(*p) - 1] == '0') || (buffg->A[(*p) - 1] == '1')
	       || (buffg->A[(*p) - 1] == '2') || (buffg->A[(*p) - 1] == '3')
	       || (buffg->A[(*p) - 1] == '4') || (buffg->A[(*p) - 1] == '5')
	       || (buffg->A[(*p) - 1] == '6') || (buffg->A[(*p) - 1] == '7')
	       || (buffg->A[(*p) - 1] == '8') || (buffg->A[(*p) - 1] == '9')
	       || (buffg->A[(*p) - 1] == '!') || (buffg->A[(*p) - 1] == '~')
	       || (buffg->A[(*p) - 1] == '-') || (buffg->A[(*p) - 1] == '+')
	       || (buffg->A[(*p) - 1] == 's') || (buffg->A[(*p) - 1] == 'S')))
	    readchrc (&(*fyle), &(*buffg), &(*p), &chrc);
	  else
            error (MISTAKE, (*p));
	}
    }
  if (orderlist)
    chrc = rreverselist (chrc);
  stackhandle (chrc,stack);
  return chrc;
} // getchrc


prodtype
getprod (text * fyle, string0 * buffg, int *p)
{
  prodtype prd1, prd2, prd3, prd;
  int n, ppg, k, m;
  char word1[MAXSTRING], word[MAXSTRING];
  bool signg, r, d, stack;
  frame aa;
  char ser;
  leftb = 0;
  rightb = 0;
  ppg = (*p);
  stack = true;
  readword (&(*buffg), &(*p), word);
  prd = NULL;
  if (interp ("rule", word, 4))
    {
      prd = rule (&(*fyle), &(*buffg), &(*p));
      prjctindex = prjctindex - 1;
    }
  else if (interp ("product", word, 1))
    {
      prd1 = getprod (&(*fyle), &(*buffg), &(*p));
      prd2 = getprod (&(*fyle), &(*buffg), &(*p));
      prd = pkronk (prd1, prd2);
    }
  else if (interp ("intdividecoeffs", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      prd = getprod (&(*fyle), &(*buffg), &(*p));
      pdiv (&prd, n);
      stack = false;
    }
  else if (interp ("swapgroup", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      prd = swapgroups (getprod (&(*fyle), &(*buffg), &(*p)), n, m);
      prjctindex = prjctindex - 1;
    }
  else if (interp ("unui", word, 4))
    {
      readint (&(*buffg), &(*p), &m);
      prd = unu1 (m);
    }
  else if (interp ("hecke", word, 5))
    {
      readpart (&(*buffg), &(*p), &aa);
      prd = hecke (aa);
    }
  else if (interp ("branch", word, 2))
    {
      prd2 = pbranch ( &(*buffg), &(*p));
      if (prd2 !=NULL) {
         prd = pmodify (prd2);
         prjctindex--;
      }
      else 
      {
	 if (debug_schur)
	      fprintf(stderr,"pbranch returned NULL\n");
	 fprintf(stderr,"Warning: this branching rule number may be false\n");
	 fprintf(stderr,"Type ?Table for the table of branching rules.\n");
	 if (prjctindex >0) prjctindex--;
      }
    }
  else if (interp ("ch_labelforon", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 6);
  else if (interp ("std_onedprep", word, 5))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 9);
  else if (interp ("mult_select", word, 6))
    {
      readint (&(*buffg), &(*p), &n);
      prd = getprod (&(*fyle), &(*buffg), &(*p));
      prestrict (&prd, n, 'm');
      prjctindex = prjctindex - 1;
    }
  else if (interp ("wt", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      prd = pwrestrict (getprod (&(*fyle), &(*buffg), &(*p)), m, n, 'w');
      prjctindex = prjctindex - 1;
    }
  else if (interp ("len", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      prd = pwrestrict (getprod (&(*fyle), &(*buffg), &(*p)), m, n, 'l');
      prjctindex = prjctindex - 1;
    }
  else if (interp ("uoneaddint", word, 5))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 8);
  else if (interp ("uonedivint", word, 5))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 10);
  else if (interp ("rp_firstpartbyspin", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 7);
  else if (interp ("ch_uonereps", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 1);
  else if (interp ("rm_uonewtovermax", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 2);
  else if (interp ("ch_spinindex", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 4);
  else if (interp ("rp_reporsfnbywt", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 3);
  else if (interp ("contragradientrep", word, 7))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 5);
  else if (interp ("autoorisomorphism", word, 2))
    prd = pautom (&(*fyle), &(*buffg), &(*p));
  else if (interp ("kinsert", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 14);
  else if (interp ("spstar", word, 6))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 11);
  else if (interp ("associate", word, 5))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 12);
  else if (interp ("starequivalent", word, 4))
    prd = pmu (&(*fyle), &(*buffg), &(*p), 13);
  else if (interp ("rm_group", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      if ((n <= nprod) && (n > 0))
	prd = crunchup (getprod (&(*fyle), &(*buffg), &(*p)), n);
      else
	prd = NULL;
      if (prd != NULL)
	prjctindex = prjctindex - 1;
    }
  else if (interp ("contractgroups", word, 4))
    {
      readint (&(*buffg), &(*p), &n);
      readint (&(*buffg), &(*p), &m);
      if ((n != m) && (n > 0) && (m > 0))
	{
	  k = (*p);
	  readword (&(*buffg), &k, word1);
	  if ((interp ("i_sfnproduct", word1, 1)
	       || interp ("o_sfnproduct", word1, 1)
	       || interp ("sk_sfn", word1, 2)
	       || interp ("plethysmsfnouter", word1, 2)
	       || interp ("rd_i_sfnproduct", word1, 4)
	       || interp ("mixedtensorreps", word1, 3)))
	    (*p) = k;
	  prd = contract (getprod (&(*fyle), &(*buffg), &(*p)), n, m, word1);
	  prjctindex = prjctindex - 1;
	}
      else
	error (MISTAKE, (*p));
    }
  else if (interp ("lastresult", word, 4))
    prd = prodcopy (pcurrent);
  else if (interp ("std_reps", word, 3))
    {
      prd = getprod (&(*fyle), &(*buffg), &(*p));
      prd = pmodify (prd);
      prjctindex = prjctindex - 1;
    }
  else if (interp ("mult_coeffsbyanint", word, 2))
    {
      readint (&(*buffg), &(*p), &n);
      prd = getprod (&(*fyle), &(*buffg), &(*p));
      prd = prodmult (n, prd);
    }
  else if (interp ("add", word, 3))
    {
      prd1 = getprod (&(*fyle), &(*buffg), &(*p));
      prd2 = getprod (&(*fyle), &(*buffg), &(*p));
      prd = prodadd (prd1, prd2);
    }
  else if (interp ("subtract", word, 3))
    {
      prd1 = getprod (&(*fyle), &(*buffg), &(*p));
      prd2 = getprod (&(*fyle), &(*buffg), &(*p));
      prd3 = prodmult (-1, prd2);
      prd = prodadd (prd1, prd3);
      /*pdisp(&prd3); */ /**/
    }
  else if (interp ("zero", word, 1))
    prd = NULL;
  else if (interp ("macmixedseries", word, 4))
  {
      readint (&(*buffg), &(*p), &n);
      readword (&(*buffg), &(*p), word1);
      ser = word1[0];
      signg =false;
      r= false;
      d= false;
      prd = macseries (ser, true /*w*/, signg, r, d, true/*m*/, n);
  }
  else if (interp ("macseries", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      readword (&(*buffg), &(*p), word1);
      ser = word1[0];
      signg =false;
      r= false;
      d= false;
      if (ser=='a' || ser=='b' || ser=='c' 
		      || ser=='d' )
	      d=true;
      if (ser=='e')
	      signg=true;
      
      prd = macseries (ser, true /*w*/, signg, r, d, false /*m*/, n);
    }
  else if (interp ("inverseseries", word, 3))
    {
      readint (&(*buffg), &(*p), &n);
      prd = inverseseries (getprod (&(*fyle), &(*buffg), &(*p)), n);
    }
  else if (interp ("variablefordpreps", word, 1))
    {
      readint (&(*buffg), &(*p), &n);
      if ((n >= 1) && (n <= pvarlimit))
	prd = prodcopy (pvar.A[n - 1]);
      else
	{
	  warn ("illegal variable", cr);
	  prd = NULL;
	}
    }
  else
    {
      (*p) = ppg;
      if (word[0] == cont)
	{
	  readacard (&(*fyle), &(*buffg), &(*p));
	  prd = getprod (&(*fyle), &(*buffg), &(*p));
	  stack = false;
	}
      else
	{
	  skipbl ((*buffg), &(*p));
	  if (((buffg->A[(*p) - 1] == '0') || (buffg->A[(*p) - 1] == '1')
	       || (buffg->A[(*p) - 1] == '2') || (buffg->A[(*p) - 1] == '3')
	       || (buffg->A[(*p) - 1] == '4') || (buffg->A[(*p) - 1] == '5')
	       || (buffg->A[(*p) - 1] == '6') || (buffg->A[(*p) - 1] == '7')
	       || (buffg->A[(*p) - 1] == '8') || (buffg->A[(*p) - 1] == '9')
	       || (buffg->A[(*p) - 1] == '[') || (buffg->A[(*p) - 1] == ']')
	       || (buffg->A[(*p) - 1] == '-') || (buffg->A[(*p) - 1] == '+')))
	    {
	      readprod (&(*fyle), &(*buffg), &(*p), &prd);
	      collectgarbage ();
	    }
	  else
	    error (MISTAKE, (*p));
	}
    }
  if (orderlist)
    prd = preverselist (prd);
  if (stack)
    pstackhandle (prd);
  return prd;
} // getprod

prodtype
pmu (text * fyle, string0 * buffg, int *p, int ss)
{
  register prodtype R323;
  prodtype list, newlist;
  int grno, m;
  char word[MAXSTRING];
  ocharptr temp, temp1;

  if (((ss == 1) || (ss == 2) || (ss == 8) || (ss == 10)))
    readint (&(*buffg), &(*p), &m);
  if ((ss == 14))
    {				// modified by FB to conform to examples
      readint (&(*buffg), &(*p), &m);
      readword (&(*buffg), &(*p), word);
      grno = 0;
      if (interp ("group", word, 1))
	readint (&(*buffg), &(*p), &grno);
    }
  if (((ss >= 1) && (ss <= 13)))
    {
      readword (&(*buffg), &(*p), word);
      if (interp ("group", word, 1))
	readint (&(*buffg), &(*p), &grno);
      else
	{
	  error (GROUP_NOT_SET, 0);
	  return (NULL);
	}
    }
  list = getprod (&(*fyle), &(*buffg), &(*p));
  newlist = list;
  while (newlist != NULL)
    {
      temp = newlist->prods.A[grno - 1];
      if ((ss == 1))
	newlist->prods.A[grno - 1] = muc (temp, currgrp.A[grno - 1], m);
      if ((ss == 2))
	newlist->prods.A[grno - 1] = umax (temp, currgrp.A[grno - 1], m);
      if ((ss == 3))
	newlist->prods.A[grno - 1] = rwt (temp);
      if ((ss == 4))
	newlist->prods.A[grno - 1] = cspin (temp);
      if ((ss == 5))
	newlist->prods.A[grno - 1] = contrag (temp, currgrp.A[grno - 1]);
      if ((ss == 6))
	newlist->prods.A[grno - 1] = clabel (temp);
      if ((ss == 7))
	newlist->prods.A[grno - 1] = qspin (temp);
      if ((ss == 8))
	newlist->prods.A[grno - 1] = uadd (temp, m, 1);
      if ((ss == 9))
	{
	  if (currgrp.A[grno - 1].name == unc)
	    readint (&(*buffg), &(*p), &kkz);
	  newlist->prods.A[grno - 1] = gmodify (temp, currgrp.A[grno - 1]);
	}
      if ((ss == 10))
	newlist->prods.A[grno - 1] = uadd (temp, m, 2);
      if ((ss == 11))
	newlist->prods.A[grno - 1] = star (temp, currgrp.A[grno - 1]);
      if ((ss == 12))
	newlist->prods.A[grno - 1] = associate (temp, currgrp.A[grno - 1]);
      if ((ss == 13))
	{
	  temp1 = newlist->prods.A[grno + 0];
	  newlist->prods.A[grno - 1] = star (temp, currgrp.A[grno - 1]);
	  newlist->prods.A[grno + 0] = associate (temp1, currgrp.A[grno + 0]);
	}
      if ((ss == 14))
	{
	  rrestrict (&temp, m, 'w');
	  //newlist->prods.A[grno - 1] = temp;
	  newlist->prods.A[grno - 1] = kinsert (temp, m, currgrp.A[grno - 1]);	//added by FB
	}
      odisp (&temp);
      newlist = newlist->next;
    }
  newlist = list;
  newlist = prodexpand (newlist);
  schur_psort (&newlist, true);
  R323 = newlist;
  return R323;
} // prodtype

ocharptr
gettrans (text * fyle, string0 * buffg, int *p)
{
  register ocharptr R324;
  ocharptr listchrc, lastchrc;
  termptr scon, scov, tempcon;
  int ppg;
  char tag;
  char word[MAXSTRING];
  bool spinor, mixed;

  ppg = (*p);
  spinor = false;
  mixed = false;
  scon = NULL;
  scov = NULL;
  listchrc = NULL;
  lastchrc = NULL;
  readword (&(*buffg), &(*p), word);
  if (interp ("spin", word, 4))
    {
      ppg = (*p);
      spinor = true;
      readword (&(*buffg), &(*p), word);
    }
  if (interp ("mixedtensorreps", word, 3))
    {
      ppg = (*p);
      mixed = true;
      readword (&(*buffg), &(*p), word);
    }
  if (interp ("label", word, 3))
    {
      while ((Member ((unsigned) (buffg->A[(*p) - 1]), brackets.S))
	     && ((*p) < bcol))
	(*p) = (*p) + 1;
      tag = buffg->A[(*p) - 1];
      ppg = (*p) + 1;
    }
  else
    tag = ' ';
  (*p) = ppg;
  if (mixed)
    scon = getsfn (&(*fyle), &(*buffg), &(*p));
  scov = getsfn (&(*fyle), &(*buffg), &(*p));
  if (!mixed)
    R324 = sfntochrc (scov, spinor, tag);
  else
    {
      while (scov != NULL)
	{
	  tempcon = scon;
	  while (tempcon != NULL)
	    {
	      cnu (&listchrc);
	      {
		register ocharptr W6 = &(*listchrc);

		W6->mult = tempcon->mult * scov->mult;
		W6->val = scov->val;
		W6->lab = tag;
		W6->spin = spinor;
		W6->next = lastchrc;
		W6->C6_double = true;
		W6->conval = tempcon->val;
		W6->conlab = tag;
		lastchrc = listchrc;
		tempcon = tempcon->next;
	      }
	    }
	  scov = scov->next;
	}
      osort (&listchrc, false);
      R324 = listchrc;
    }
  scollctgarbage ();
  srjctindex = 1;
  scollctgarbage ();
  return R324;
} // gettrans

prodtype
pbranch (string0 * buffg, int *p)
{
  register prodtype R325;
  char word[MAXSTRING];
  prodtype pr, top;
  ocharptr help, help1;
  int nprd, n, m, n2, m2, grno, brno;
  register int jr;

  pr = NULL;
  if (nprod >= 1)
    {
      readint (&(*buffg), &(*p), &brno);
      nprd = nprod;

      if ((brno >= 1) && (brno <= 63))
	{
	  if (((brno < 43) || (brno > 59)))
	    readint (&(*buffg), &(*p), &n);
	  if (((brno >= 4) && (brno <= 39)))
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
		readint (&(*buffg), &(*p), &m);
		break;
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
	      case 36:
		break;
	      default:
		Caseerror (Line);
	      }

	  if (((brno == 33) || (brno == 34)))
	    {
	      readint (&(*buffg), &(*p), &n2);
	      readint (&(*buffg), &(*p), &m2);
	    }
	  readword (&(*buffg), &(*p), word);
	  if (interp ("group", word, 1))
	    {
	      readint (&(*buffg), &(*p), &grno);
	      pr = pmodify (getprod (&input, &(*buffg), &(*p)));
	      for (jr = grno; jr <= nprod - 1; jr++)
		{
		  currgrp.A[jr - 1] = currgrp.A[jr + 0];
		}
	      fixsubgroup (&brno, &n, &m, &n2, &m2);

	      if (((brno >= 43) && (brno <= 60)))
		exceptgroupload (&brno);
	      if (brno == 63)
		exceptgroupload (&brno);
	      top = pr;
	      /*echo = true; *//*21/3/98 */
	      while (pr != NULL)
		{
		  register prodtype W9 = &(*pr);

		  help = W9->prods.A[grno - 1];
		  help1 = W9->prods.A[nprod - 1];
		  for (jr = grno; jr <= nprd - 1; jr++)
		    {
		      W9->prods.A[jr - 1] = W9->prods.A[jr + 0];
		    }
		  dobranch (&help, &help1, brno, n, m, n2, m2);
		  W9->prods.A[nprd - 1] = help1;
		  pr = pr->next;
		}
	      pr = top;
	      n = nprod;
	      nprod = nprd;
	      pr = prodexpand (top);
	      nprod = n;
	      pdisp (&top);
	      top = pr;
	      fixgroups (&pr, brno, nprd);
	      pr = top;
	      /*su1crunch(&pr); */
	      if (echo)
		putgroup1 (currgrp, 1);
	      /*echo = false; *//*1/1/98 */
	      schur_psort (&pr, true);
	    }
	  else
	    inform ("Sorry -try again", cr);
	}
      else
	print ( "Inappropriate Branching Rule Number\n");
    }
  else
    error (GROUP_NOT_SET, 2);
  R325 = pr;
  return R325;
} // pbranch

prodtype
rule (text * fyle, string0 * buffg, int *p)
{
  int op[maxprod + 1];
  register prodtype R326;
  int n2, n1, gn, k, w, m, rnk;
  register int jr, ir;
  frame v;
  bool sum, mylist, series;
  char word[MAXSTRING];
  termptr zeta, z, x, y, u, list, plist;
  ocharptr ch, help;
  prodtype sublist, prd, pr, top;

  prd = getprod (&(*fyle), &(*buffg), &(*p));
  sublist = NULL;
  for (jr = 0; jr <= maxprod; jr++)
    {
      op[jr] = 0;
    }
  k = (*p);
  n2 = 0;
  n1 = 0;
  sum = false;
  mylist = false;
  readword (&(*buffg), &k, word);
  if (interp ("sum", word, 3)) {
      (*p) = k;
      readword (&(*buffg), &k, word);
      sum = true;
      mylist = interp ("mylistofsfns", word, 2);
      if (mylist)
	(*p) = k;
    }
  do {
      m = (*p);
      readword (&(*buffg), &(*p), word);
      k = (*p);
      gn = 0;
      if (strlen(stripwhite(word))!=0
	  && !interp ("with", word, 4)) {
	      int oldp=*p;
	  readint (&(*buffg), &(*p), &gn);
	  if (gn > nprod) {
	      gn = 0;
	      error (MISTAKE, oldp);
	      print ("(bad product number.)\n");
	    }
	}
      if ((gn != 0) && (!interp ("with", word, 4)))
	if (interp ("i_plethysmrd", word, 4))
	  op[gn] = 8;
	else if (interp ("i_sfnproduct", word, 1)) {
	    op[gn] = 1;
	    n2 = gn;
	  }
	else if (interp ("o_sfnproduct", word, 1))
	  op[gn] = 2;
	else if (interp ("sk_sfn", word, 2)) {
	    op[gn] = 3;
	    n1 = gn;
	  }
	else if (interp ("conjugatesfnlist", word, 4))
	  op[gn] = 4;
	else if (interp ("equalsfnlist", word, 2))
	  op[gn] = 5;
	else if (interp ("ch_phaseofsfns", word, 4))
	  op[gn] = 6;
	else if (interp ("plethysmsfnouter", word, 2))
	  op[gn] = 7;
	else if (interp ("rp_reporsfnbywt", word, 4))
	  op[gn] = 9;
	else if (interp ("makewtsfn", word, 4))
	  op[gn] = 10;
    }
  while (!(interp ("with", word, 4) || ((*p) >= bcol)));
  if (!interp ("with", word, 4))
    (*p) = m;
  if (sum) {
      if (mylist)
	zeta = getsfn (&(*fyle), &(*buffg), &k);
      (*p) = k;
      top = prd;
      while (top != NULL) {
	  register prodtype W14 = &(*top);

	  if (!mylist) {
	      if (n1 != 0)
		zeta = skcompat (W14->prods.A[n1 - 1]->val);
	      else
		zeta = eqwt (wtfrm (&W14->prods.A[n2 - 1]->val));
	      if ((n1 > 0) && (n2 > 0)) {
		  w = wtfrm (&W14->prods.A[n2 - 1]->val);
		  x = NULL;
		  y = zeta;
		  u = NULL;
		  while (y != NULL) {
		      if (wtfrm (&y->val) == w) {
			  snu (&u);
			  (*u) = (*y);
			  u->next = x;
			  x = u;
			}
		      y = y->next;
		    }
		  ldisp (&zeta);
		  zeta = x;
		}
	    }
	  z = zeta;
	  while (z != NULL) {
	      pnu (&pr);
	      (*pr) = (*top);
	      pr->next = NULL;
	      for (jr = 1; jr <= nprod; jr++) {
		  cnu (&ch);
		  (*ch) = (*W14->prods.A[jr - 1]);
		  ch->next = NULL;
		  pr->prods.A[jr - 1] = ch;
		  ch = NULL;
		}

	      for (jr = 1; jr <= nprod; jr++) {
		  if (((op[jr] <= 3) || (op[jr] == 7))) {
		      help = pr->prods.A[jr - 1];
		      if ((op[jr] == 1))
			list = inner (help->val, z->val);
		      if ((op[jr] == 2))
			list = outer (help->val, z->val);
		      if ((op[jr] == 3))
			list = skew (help->val, z->val);
		      if ((op[jr] == 7))
			list = pleth (help->val, z->val);
		      if ((op[jr] == 0)) {//  keep the group unchanged
			snu(&list);
			list->val =help->val;
			list->mult=1;
		      }
		      pr->prods.A[jr - 1] = sfntochrc (list, help->spin, help->lab);
		      ldisp (&list);
		      odisp (&help);
		    }
		  if (op[jr] == 6) {
		      pr->prods.A[jr - 1]->mult =
			z->mult * minusoneto (wtfrm (&z->val));
		    }
		  if (op[jr] == 4) {
		      v = z->val;
		      conjgte (&v);
		      pr->prods.A[jr - 1]->val = v;
		    }
		  if (op[jr] == 5)
		    pr->prods.A[jr - 1]->val = z->val;
		}
	      pr->next = sublist;
	      sublist = pr;
	      z = z->next;
	    }
	  ldisp (&zeta);
	  top = top->next;
	  if ((op[jr] == 5) && mylist) {
	      ldisp (&zeta);
	      srjctindex = srjctindex - 1;
	    }
	}
      if ((op[jr] == 5) && (!mylist))
	ldisp (&zeta);
      pstackhandle (sublist);
      schur_psort (&sublist, true);
      R326 = sublist;
    } else { 			// if not sum
      series = false;
      zeta = NULL;
      (*p) = k;

      if (!((op[gn] == 4) || (op[gn] == 8) || (op[gn] == 10))) {
	  if ((n1 == 0) && (op[gn]==0)) {   // keep the group unchanged
		  list =NULL;
	  }
	  else if (n1 == 0)
	    list = getsfn (&(*fyle), &(*buffg), &(*p));
	  else {
	      readword (&(*buffg), &k, word);
	      k = (*p);
	      if ((!(Member ((unsigned) (word[0]), numbers.S)))
		  && (word[1] == ' ') && (word[0] != '!'))
		series = true;
	      else
		list = getsfn (&(*fyle), &(*buffg), &(*p));
	    }
	}
      top = prd;
      snu (&plist);
      while (top != NULL) {
	  register prodtype W19 = &(*top);
	  for (jr = 1; jr <= nprod; jr++) {
	      if (op[jr] != 0) {
		  plist->val = W19->prods.A[jr - 1]->val;
		  plist->mult = W19->prods.A[jr - 1]->mult;
		  switch ((int) (op[jr])) {
		    case 1:
		      zeta = linner (plist, list);
		      ldisp (&list);
		      srjctindex = srjctindex - 1;
		      break;
		    case 2:
		      zeta = louter (plist, list);
		      ldisp (&list);
		      srjctindex = srjctindex - 1;
		      break;
		    case 3:
		      if (!series) {
			  zeta = lskew (plist, list);
			  /*ldisp(&list);
			     srjctindex = srjctindex - 1; */
		      } else {
			  zeta =
			    useseries (locase (word[0]), plist,
				       false, true, (*p));
			  if (top->next != NULL)
			    (*p) = k;
		      }
		      break;
		    case 4:
		      conjgte (&W19->prods.A[jr - 1]->val);
		      break;
		    case 6:
		      W19->prods.A[jr - 1]->mult =
			W19->prods.A[jr -
				     1]->mult *
			minusoneto (wtfrm (&W19->prods.A[jr - 1]->val));
		      break;
		    case 7:
		      zeta = listplethlist (plist, list);
		      ldisp (&list);
		      srjctindex = srjctindex - 1;
		      break;
		    case 8:
		      nreduce = true;
		      zeta = plethonerinner (plist);
		      break;
		    case 9:
		      W19->prods.A[jr - 1]->val.A[1] =
			wtfrm (&W19->prods.A[jr - 1]->val);
		      for (ir = 2; ir <= maxdim; ir++) {
			  W19->prods.A[jr - 1]->val.A[ir] = 0;
			}
		      break;
		    case 10:
		      rnk = currgrp.A[jr - 1].rank;
		      zeta = makeweight (rnk, plist);
		      nreduce = false;
		      break;
		    default:
		      Caseerror (Line);
		    }

		  switch ((int) (op[jr])) {
		    case 1:
		    case 2:
		    case 3:
		    case 7:
		    case 8:
		    case 10:
		      help = W19->prods.A[jr - 1];
		      W19->prods.A[jr - 1] =
			sfntochrc (zeta, help->spin, help->lab);
		      odisp (&help);
		      ldisp (&zeta);
		      break;
		    case 4:
		    case 5:
		    case 6:
		    case 9:
		      break;
		    default:
		      Caseerror (Line);
		    }
		}
	    }
	  top = top->next;
	}
      ldisp (&plist);
      schur_psort (&prd, true);
      R326 = prd;
    }
  return R326;
} // rule

prodtype
pautom (text * fyle, string0 * buffg, int *p)
{
  register prodtype R327;
  prodtype newlist, ptr, list;
  int iauto, grno;
  char word[MAXSTRING];
  groop tg, fg;
  bool test;
  ocharptr temp;

  test = false;
  newlist = NULL;
  readword (&(*buffg), &(*p), word);
  if (interp ("group", word, 1))
    {
      readint (&(*buffg), &(*p), &grno);
      fg = currgrp.A[grno - 1];
      test = readauto (&(*buffg), &(*p), &iauto, &fg, &tg, grno);
      list = getprod (&(*fyle), &(*buffg), &(*p));
      newlist = prodcopy (list);
      ptr = newlist;
      if ((test == true))
	{
	  while (newlist != NULL)
	    {
	      temp = newlist->prods.A[grno - 1];
	      newlist->prods.A[grno - 1] = change (temp, iauto, fg, tg);
	      odisp (&temp);
	      newlist = newlist->next;
	    }
	  newlist = ptr;
	  newlist = prodexpand (newlist);
	  pdisp (&ptr);
	  schur_psort (&newlist, true);
	}
    }
  if (echo)
    putgroup1 (currgrp, 1);
  R327 = newlist;
  return R327;
}

void
readprod (text * fyle, string0 * buffg, int *p, prodtype * plist)
{
  prodtype sublist, pr;
  int ssign, jg, k, entry;
  /*register */ int pl;
  string0 buffg2;
  sublist = NULL;
  skipbl ((*buffg), &(*p));
  do
    {
      skipbl ((*buffg), &(*p));
      if ((*p) < bcol)
	{
	  pnu (&pr);
	  {
	    register prodtype W24 = &(*pr);

	    ssign = 1;
	    if (((buffg->A[(*p) - 1] == '+') || (buffg->A[(*p) - 1] == '-')))
	      {
		if (buffg->A[(*p) - 1] == '-')
		  ssign = -1;
		(*p) = (*p) + 1;
	      }
	    if (buffg->A[(*p) - 1] == cont)
	      readacard (&(*fyle), &(*buffg), &(*p));
	    skipbl ((*buffg), &(*p));
	    if (buffg->A[(*p) - 1] != '[')
	      readint (&(*buffg), &(*p), &W24->mult);
	    else
	      W24->mult = 1;
	    if (buffg->A[(*p) - 1] == '[')
	      leftb = leftb + 1;
	    W24->mult = W24->mult * ssign;
	    if ((*p) < bcol)
	      (*p) = (*p) + 1;
	    jg = 1;
	    k = (*p);
	    while (jg <= nprod)
	      {
		entry = 1;
		memset(buffg2.A,' ',bcol);
		if (buffg->A[k - 1] == ']')
		  rightb = rightb + 1;
		while ((k < bcol) && !((buffg->A[k - 1] == qt)
			    || (buffg->A[k - 1] == '*')
			    || (buffg->A[k - 1] == ']')))
		  {
		    buffg2.A[entry - 1] = buffg->A[k - 1];
		    entry = entry + 1;
		    k = k + 1;
		    if (buffg->A[k - 1] == cont)
		      readacard (&(*fyle), &(*buffg), &k);
		  }
		if (buffg->A[k - 1] == ']')
		  rightb = rightb + 1;
		if (k < bcol)
		  k = k + 1;
		pl = 1;
		W24->prods.A[jg - 1] = getchrc (&(*fyle), &buffg2, &pl);

		if (W24->prods.A[jg - 1] != NULL)
		  W24->prods.A[jg - 1]->val.A[maxdim] =
		    W24->prods.A[jg - 1]->conval.A[1];
		rjctindex = rjctindex - 1;
		jg = jg + 1;
	      }
	    W24->next = sublist;
	    sublist = pr;
	  }
	  (*p) = k;
	  skipbl ((*buffg), &(*p));
	}
    }
  while (! (!((buffg->A[(*p) - 1] == '+') || (buffg->A[(*p) - 1] == '-'))
	  || ((*p) >= bcol)));
  skipbl ((*buffg), &(*p));
  if ((buffg->A[(*p) - 1] == delim) && ((*p) < bcol))
    (*p) = (*p) + 1;
  if (((leftb == rightb) && (leftb != 0)) && (k != bcol))
    {
      schur_psort (&sublist, true);
      (*plist) = sublist;
      sublist = NULL;
    }
  else
    {
      (*plist) = NULL;
      //error (1, (*p));
      pdisp (&pr);
      print ( "missing brackets or group not set?\n");
    }
} //readprod
