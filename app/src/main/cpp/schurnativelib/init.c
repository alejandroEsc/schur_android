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

/*\file init.c
 * This file initialises variables and generates a list of prime numbers used in big number routines.*/

/*
**	Definitions for i/o
*/
# include <stdio.h>
# include <math.h>
# include <sys/types.h>
# include <dirent.h>
# include <limits.h>
# include <strings.h>

# include "standard.h"
# include "define.h"
# include "sets_mgmt.h"
# include "mymalloc.h"
# include "dim.h"
# include "type.h"
# include "var.h"
# include "utils.h"
# include "s1.h"
# include "g.h"
# include "init.h"

/*
**	Start of program definitions
*/

prodtype discardedprds, pcurrent;

/*Generates a list of prime numbers. Called by b.c. It will be used by bignums*/
void
genprime (bframe * primesi)
{
  int ji, k;
  register int l;
  bool prime;
  for (l = 0; l <= maxl; l++)
    {
      primesi->A[l] = 0;
    }
  primesi->A[0] = 1;
  primesi->A[1] = 2;
  primesi->A[2] = 3;
  ji = 0;
  k = 2;
  while (k < maxl)
    {
      ji = ji + 6;
      prime = true;
      for (l = 2; l <= Trunc (sqrt ((double) ji - 1)); l++)
	{
	  if (((ji - 1) % l) == 0)
	    prime = false;
	}
      if (prime)
	{
	  k = k + 1;
	  primesi->A[k] = ji - 1;
	}
      prime = true;
      for (l = 2; l <= Trunc (sqrt ((double) ji + 1)); l++)
	{
	  if (((ji + 1) % l) == 0)
	    prime = false;
	}
      if (prime)
	{
	  k = k + 1;
	  if (k <= maxl)
	    primesi->A[k] = ji + 1;
	}
    }
}

/*Initializes globally defined variables from var.h.*/
void
initialise (void)
{
  register int j1;
  register int i1;
  char ostype[100] = OSTYPE;

  displayModes[sfnm] = "SFN> ";
  displayModes[brm] = "BRM> ";
  displayModes[repm] = "REP> ";
  displayModes[dpm] = "DP> ";

  // strcpy(redLine, "\x1B[31;40m"); red on black
  strcpy (colorsBegin, "\x1B[32;40m");	//green on black
  strcpy (colorsEnd, "\x1B[00m");	//back to normal
  strcpy (ostype, OSTYPE);	//OSTYPE is a defined value created by ./configure
  if (ostype == NULL || strncasecmp (ostype, "solaris", 7) == 0)	// pbs with sun terminals
    {
      colorsBegin[0] = '\0';
      colorsEnd[0] = '\0';
    }

  instr = mymalloc (MAXSTRING);
  genprime (&primes);

  strtoSet ("?_abcdefghijklmnopqrstuvwxyz", letterset.S);
  strtoSet ("0123456789", numbers.S);
  strtoSet (" ()<>{}", brackets.S);
  strtoSet ("!+-~", plus_minus_etc.S);
  strtoSet ("0123456789!~", numbersEtc.S);

/*
  printf("0:");
  displayCharsOfSet(Conset[0]);
  printf("1:");
  displayCharsOfSet(Conset[1]);
  printf("5:");
  displayCharsOfSet(Conset[5]);
  printf("8:");
  displayCharsOfSet(Conset[8]);
  printf("9:");
  displayCharsOfSet(Conset[9]);
  printf("10:");
  displayCharsOfSet(Conset[10]);
  printf("11:");
  displayCharsOfSet(Conset[11]);
 */ 

  start = 1;
  fin = 2;
  bs = 8;
  cr = 13;
  esc = 27;
  bel = 7;
  backspace = bs;
  lt = 11;
  rt = 12;
  go = 15;
  nth = 21;
  sth = 10;
  del = 19;
  //tlines = 40;
  columns = 5;
  logging = false;
  sb_conj = true;
  more = false;
  digits = true;
  rep = false;
  logopen = false;
  echo = true;
  bell = false;
  nreduce = false;
  gotline = false;
  iosup = false;
  fnex = false;
  poly = false;
  xreduce = false;
  mmono = false;
  homo = false;
  elem = false;
  forg = false;
  dimb = true;
  ibm = false;
  tex = false;
  sslab = false;
  qtest = false;
  color = false;
  psum = false;
  fnptr = NULL;
  pfn = false;
  wprod = false;
  pow_note = true;
  sb_ListOutput = false;
  sb_prog = false ; /* if true dislay progression of computation by a rotating bar*/
  /* plen = maxdim;  was used for sound duration */ 
  maxb = maxl;
  plwt = maxdim;
  cutoff = 1;
  discardedsfns = NULL;
  scurrent = NULL;
  terms = 0;
  prompt[0]='\0';
  for (i1 = 1; i1 <= bcol; i1++)
      logname.A[i1 - 1] = ' ';
  qfn = false;
  for (i1 = 0; i1 <= maxdim; i1++)
    {
      nolls.A[i1] = 0;
      full.A[i1] = INT_MAX;
    }
  for (i1 = 0; i1 <= maxl; i1++)
      nulls.A[i1] = 0;
  for (i1 = 1; i1 <= bcol; i1++)
      buffz.A[i1 - 1] = ' ';
  buffer = buffz;
  buffi = buffz;

  full.A[maxdim] = 0;
  full.length = maxdim - 1;
  nolls.length = 0;
  for (i1 = 1; i1 <= svarlimit; i1++)
      svar.A[i1 - 1] = NULL;
  for (i1 = 1; i1 <= maxfn; i1++)
      fn.A[i1 - 1] = NULL;
  discardedchrcs = NULL;
  f4load = false;
  e6load = false;
  e7load = false;
  e8load = false;
  e6f4load = false;
  e8soload = false;
  e8suload = false;
  e8e6load = false;
  e6soload = false;
  e7e6load = false;
  e6g2load = false;
  u27e6load = false;
  su56e7load = false;
  su248e8load = false;
  e6f4load = false;
  f4g2load = false;
  e6su3g2load = false;
  e8f4g2load = false;
  racah = false;
  redu = false;
  qspecial = true;
  for (svarindex = 1; svarindex <= svarlimit; svarindex++)
      svar.A[svarindex - 1] = NULL;
  for (srjctindex = 1; srjctindex <= srjctlimit; srjctindex++)
      sreject.A[srjctindex - 1] = NULL;
  srjctindex = 0;
  scurrent = NULL;
  orderlist = false;
  discardedprds = NULL;
  for (varindex = 1; varindex <= rvarlimit; varindex++) 
      vari.A[varindex - 1] = NULL;
  for (rjctindex = 1; rjctindex <= rjctlimit; rjctindex++) {
      reject.A[rjctindex - 1] = NULL;
      reject2.A[rjctindex - 1] = NULL;
    }
  rjctindex = 0;
  rjctindex2 = 0;
  ggroup.name = nill;
  ggroup.rank = 0;
  ggroup.rank2 = 0;
  current = NULL;
  pcurrent = NULL;
  count = 0;
  seccount = 0;
  rcount = 0;
  pcount = 0;
  setlimit = 12;
  weight = false;
  for (pvarindex = 1; pvarindex <= pvarlimit; pvarindex++)
      pvar.A[pvarindex - 1] = NULL;
  for (prjctindex = 1; prjctindex <= prjctlimit; prjctindex++)
      preject.A[prjctindex - 1] = NULL;
  prjctindex = 0;
  for (j1 = 1; j1 <= maxprod; j1++) {
	currgrp.A[j1-1].name = nill;
	currgrp.A[j1-1].rank = 0;
	currgrp.A[j1-1].rank2 = 0;
    }

  discardedfns = NULL;
  for (j1 = 1; j1 <= maxfn; j1++)
      fn.A[j1 - 1] = NULL;

  nprod = 0;
  qsn = 0;

  dataPath = getenv ("SCHURLIB");
  if (dataPath == NULL || strlen(dataPath)==0) {
      dataPath = (char *) malloc (sizeof (char) * 600);
      if (opendir (DATADIR) != NULL)
	strcpy (dataPath, DATADIR);
      else {
	  if (strcmp (ostype, "cygwin") == 0) {	//on windows try to access registry
	      FILE *f;
	      f =
		fopen
		("/proc/registry/HKEY_LOCAL_MACHINE/SOFTWARE/SCHUR/Settings/SCHURLIB",
		 "r");
	      if (f != NULL) {
		  char tempo[500];
		  fgets (tempo, 500, f);
		  //sprintf (dataPath, "%s/share/schur", tempo);
		  sprintf (dataPath, "/share/schur");
		}
	      else {
		fprintf (stderr,
			 "Sorry I can't find any data directory (windows problem ?)\n");
		  exit (EXIT_FAILURE);
	      }
	    }
	  else if (opendir (DEFAULTDATADIR) !=NULL)
	  	strcpy(dataPath,DEFAULTDATADIR);
	  else{
	    fprintf (stderr,
		     "Sorry, I can't find any valid data directory (%s or %s) - try to use SCHURLIB env. variable\n", dataPath, DEFAULTDATADIR);
	    exit (EXIT_FAILURE);
          }
	}
    }
    else if (debug_schur)
	fprintf(stderr,"SCHURLIB environment variable found ");
  if (debug_schur)
    fprintf (stderr, "dataPath=%s\n", dataPath);
}
