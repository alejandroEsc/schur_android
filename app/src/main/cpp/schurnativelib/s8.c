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
# include <string.h>
# include <stdlib.h>
# include <errno.h>

# include "standard.h"
# include "define.h"
# include "ReadWrite.h"
# include "mymalloc.h"
# include "Scanck.h"
# include "sets_mgmt.h"

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
# include "s6.h"
//# include "loadcomList.h"
# include "savecomList.h"
# include "apropos.h"
# include "FrontPage.h"
# include "s8.h"

//# define  Cmpstr(x, y) strncmp((x), (y), sizeof(x))
# define mydisplay(m)   printf ("%s%s%s",colorsBegin,m,colorsEnd);

//static FILE *Tmpfil;

bool
getfunc (string0 buff, int *p, int *n)
{
  readint (&buff, &(*p), &(*n));
  /*if (!(Member((unsigned)((*n)), Conset[7])) && !erred) */
  if (((*n < 1) || (*n > maxfn)) && !erred)
    error (BAD_FUNCTION, (*p));
  return (bool) (!erred);
}

/**
 * Display function 
 */
void
writefn (fnptrs fnx)
{
  fnptrs f;
  char s[MAXSTRING], *p;

  f = fnx;
  while (f != NULL)
    {
      memcpy(s,f->gbuff.A,bcol);
      p=findLastSpace(s);
      *p='\0';
      print ("%s\n",s);
      f = f->next;
    }
}

/* Functions common to all modes */
bool
common (modes md, char *instr8c, string0 buffer8, int *px8)
{
  register bool found;
  bool secho = false;
  char  hit;
  string0 fnname;
  text fnfile;

  found = true;
  if (interp ("apropos", instr8c, 4))
    {
      char wd[MAXSTRING];

      readword (&buffer8, &(*px8), wd);
      apropos (wd);
    }
  else if (interp ("columns", instr8c, 3)) {
	  int savecolumns=columns;
    readint (&buffer8, &(*px8), &columns);
    if (columns==0) {
	    columns=savecolumns;
	    fprintf(output.fp, "columns=%d\n",columns);
    }
  }
  else if (interp ("entervar", instr8c, 3))
    enter ();
  else if (interp ("fn", instr8c, 2))
    {
      if ((getfunc (buffer8, &(*px8), &fnn) && (fn.A[fnn - 1] != NULL)))
	{
	  fnex = true;
	  iosup = true;
	  fnptr = fn.A[fnn - 1];
	  return(true);
	}
      else
	(void) fprintf (output.fp, "Error:function not set\n"),
	  Putl (output, 1);
    }
  else if (interp ("heapstatus", instr8c, 4))
    print ("Count  = %u Rcount = %u Pcount = %u\n", count, rcount, pcount);
  else if (interp ("logfile", instr8c, 3))
    logcom (buffer8, &(*px8));
  else if (interp ("loadfile", instr8c, 2)) {
/*	if (sb_ListOutput) 
	   loadcomList(buffer8, *px8);
    	else */
	   loadcom (buffer8, *px8);
  }
  else if (interp ("remark", instr8c, 3))
    ;
  else if (interp ("readfnfromdisk", instr8c, 5)) {
      if (getfunc (buffer8, &(*px8), &fnn)) {
	  hit = skipbl (buffz, &(*px8));
	  if (hit == qt) {
	      readfilename (buffz, &(*px8), &fnname);
	      if (!erred && (strlen(stripwhite(fnname.A))!= 0)) {
		  if (findfile (fnname)) {
		      Resetx2 (&fnfile, fnname.A, -1);
		      readfnset (&fnfile, &fn.A[fnn - 1]);
		      fclose (fnfile.fp);
		      Putchr ('\n', output);
		    }
		  else
		    error (NO_SUCH_FILE, *px8);
		}
	    }
	  else
	    error (MISSING_QUOTE, *px8);
	}
    }
  else if (interp ("savesetvar", instr8c, 4)) {
      if (sb_ListOutput)
	savecomList (buffer8, *px8);
      else
	savecom (buffer8, *px8);
    }
  else if (interp ("sb_allfalse", instr8c, 6))
    allfalse ();
  else if (interp ("sb_weight", instr8c, 4))
    boole (buffer8, &(*px8), &weight);
  else if (interp ("sb_bell", instr8c, 4))
    boole (buffer8, &(*px8), &bell);
  else if (interp ("sb_conj", instr8c, 6))
    boole (buffer8, &(*px8), &sb_conj);
  else if (interp ("sb_cut", instr8c, 5)) {
	  int savecutoff=cutoff;
    readint (&buffer8, &(*px8), &cutoff);
    if (cutoff==0) {
	    cutoff=savecutoff;
	    fprintf(output.fp, "cutoff=%d\n",cutoff);
    }
  }
  else if (interp ("sb_debug", instr8c, 8))
    boole (buffer8, px8, &debug_schur);
  else if (interp ("sb_dimension", instr8c, 6))
    boole (buffer8, &(*px8), &dimb);
  else if (interp ("sb_digits", instr8c, 4))
    boole (buffer8, px8, &digits);
  else if (interp ("sb_echo", instr8c, 4))
    {
      boole (buffer8, &(*px8), &secho);
      if (secho)
	echo = true;
      else
	echo = false;
    }
  else if (interp ("sb_listoutput", instr8c, 7))	// added by FB 20060315
    boole (buffer8, &(*px8), &sb_ListOutput);
  else if (interp ("sb_more", instr8c, 4))
    boole (buffer8, &(*px8), &more);
  else if (interp ("sb_progress", instr8c, 7))
    boole (buffer8, &(*px8), &sb_prog);
  else if (interp ("sb_powernotation", instr8c, 4))
    boole (buffer8, &(*px8), &pow_note);
  else if (interp ("sb_qfn", instr8c, 4))
    boole (buffer8, &(*px8), &qfn);
  else if (interp ("sb_rd_notation", instr8c, 5))
    boole (buffer8, &(*px8), &nreduce);
  else if (interp ("sb_reverseorderbool", instr8c, 4))
    boole (buffer8, &(*px8), &orderlist);
  else if (interp ("sb_texoutput", instr8c, 4))
    {
      boole (buffer8, &(*px8), &tex);
      pow_note = true;
    }
  else if (interp ("sb_wprod", instr8c, 5))
    boole (buffer8, &(*px8), &wprod);
  else if (interp ("setfn", instr8c, 4))
    {
      if (getfunc (buffer8, &(*px8), &fnn))
	readfnset (&input, &fn.A[fnn - 1]);
    }
/* changing the size of bignums... too dangerous for the moment */
/*  else if (interp ("set_maxb", instr8c, 8))
    readint (&buffer8, &(*px8), &maxb);*/
/* was used to parameter sound duration 
  else if (interp ("set_plen", instr8c, 8)) 
    readint (&buffer8, &(*px8), &plen); */
  else if (interp ("set_pwt", instr8c, 6)) {
	  int saveplwt=plwt;
    readint (&buffer8, &(*px8), &plwt);
    if (plwt==0) {
	    plwt=saveplwt;
	    fprintf (output.fp, "pwt=%d\n",plwt);
    }
  }
  else if (interp ("setlimit", instr8c, 6)) {
	  int savelimit=setlimit;
    readint (&buffer8, &(*px8), &setlimit);
    if (setlimit==0) {
	    setlimit=savelimit;
	    fprintf (output.fp, "setlimit=%d\n",setlimit);
    }
  }
  else if (interp ("statusofschur", instr8c, 4))
    statuscom (md);
  /*else
     if (interp("sound",instr8c,3)){readint(&buffer8, &(*px8), &plen);
     sound(plen);} */
  else if (interp ("?", instr8c, 1))
    helpcom (buffer8, 1);
  else if (interp ("help", instr8c, 4))
    helpcom (buffer8, 0);
  else if (interp ("pause", instr8c, 3))
    {
      int key;

      fprintf (output.fp, "to continue strike any key"), Putl (output,
								      0);
      key= Getchr (input); UNUSED(key);
      Getl (&input);
    }
  else if (interp ("supressoutputtoscreen", instr8c, 3))
    boole (buffer8, &(*px8), &iosup);
  else if (interp ("wrfntodisk", instr, 7))
    {
      if (getfunc (buffer8, &(*px8), &fnn))
	{
	  hit = skipbl (buffer8, &(*px8));
	  if (hit == qt)
	    {
	      readfilename (buffer8, &(*px8), &fnname);
	      if (!erred)
		{
		  if (overwrite (fnname))
		    {
		      fnptrs fnpt;

		      Rewritex (fnfile, fnname.A, sizeof (fnname.A));
		      fnpt = fn.A[fnn - 1];
		      while (fnpt != NULL)
			{
			  for (jz = 1; jz <= bcol; jz++)
			    {
			      Putchr (fnpt->gbuff.A[jz - 1], fnfile);
			    }
			  Putchr ('\n', fnfile);
			  fnpt = fnpt->next;
			}
		      (void) fprintf (fnfile.fp, "stop\n");
		      fclose (fnfile.fp);
		      Putchr ('\n', output);
		    }
		}
	    }
	  else
	    error (MISSING_QUOTE, *px8);
	}
    }
  else if (interp ("wrfntoscreen", instr8c, 2))
    {
      if (getfunc (buffer8, &(*px8), &fnn))
	writefn (fn.A[fnn - 1]);
    }
  else
    found = false;
  return found;
}

void
boole (string0 buffer8, int *px8, bool * b)
{
  char wd[MAXSTRING];

  skipbl (buffer8, &(*px8));
  readword (&buffer8, &(*px8), wd);
  if (interp ("true", wd, 1))
    (*b) = true;
  else if (interp ("false", wd, 1))
    (*b) = false;
  else if (echo)
    print ("%s\n", (*b ? "True" : "False"));
}

void
logcom (string0 buffer8, int *px8)
{
  char hit;
  bool b = false;
  string0 l;

  hit = skipbl (buffer8, &(*px8));
  if (hit == qt)
    {
      readfilename (buffer8, &(*px8), &l);
      if (!erred)
	if (l.A[1 - 1] == 0)
	  {
	    if (!logopen)
	      error (NO_LOG_FILE, (*px8));
	    else
	      {
		fclose (logfile.fp);
		logopen = false;
		logging = false;
	      }
	  }
	else if (overwrite (l))
	  {
	    if (logopen)
	      fclose (logfile.fp);
	    logname = l;
	    Rewritex (logfile, l.A, sizeof (l.A));
	    logging = true;
	    logopen = true;
	  }
    }
  else
    {
      boole (buffer8, &(*px8), &b);
      if (!erred)
	if (!logopen)
	  error (NO_LOG_FILE, (*px8));
	else
	  logging = b;
    }
}

void
statuscom (modes m)
{
  register int i8, i;

  print
    ("digits:%c reverse:%c more:%c TeX:%c debug:%c setlimit:%d pwt:%d logging:",
     boolchar(digits), boolchar(orderlist), boolchar(more), boolchar(tex),
     boolchar(debug_schur), setlimit, plwt);
  if (logging)
    fprintf (output.fp, "'%s'\n", logname.A);
  else
    fprintf (output.fp, "%c\n", boolchar(logging));

  print ("     ");
  for (i = 1; i <= maxfn; i++)
    {
      if (i % 10 == 0)
	print ("%d", i);
      else if (i % 10 != 1 || i == 1)
	print (".");
    }
  print ("\n");
  print (" fns:");
  for (i8 = 1; i8 <= maxfn; i8++)
    {
      if (fn.A[i8 - 1] != NULL)
	print ("%1d", i8 % 10);
      else if (i8 % 10 == 0)
	print (".");
      else
	print ("_");
    }
  print ("\n");
  switch (m)
    {
    case sfnm:
      print ("svar:");
      for (i8 = 1; i8 <= svarlimit; i8++)
	{
	  if (svar.A[i8 - 1] != NULL)
	    print ("%1d", i8 % 10);
	  else if (i8 % 10 == 0)
	    print (".");
	  else
	    print ("_");
	}
      print (" Schur Function mode");
      break;
    case repm:
      print ("rvar:");
      for (i8 = 1; i8 <= pvarlimit; i8++)
	{
	  if (vari.A[i8 - 1] != NULL)
	    print ("%1d", i8 % 10);
	  else if (i8 % 10 == 0)
	    print (".");
	  else
	    print ("_");
	}
      print (" REPresentation mode");
      break;
    case dpm:
      print (" var:");
      for (i8 = 1; i8 <= rvarlimit; i8++)
	{
	  if (pvar.A[i8 - 1] != NULL)
	    print ("%1d", i8 % 10);
	  else if (i8 % 10 == 0)
	    print (".");
	  else
	    print ("_");
	}
      print (" Direct Product mode");
      break;
    default:
      aaargh ("bad mode.", true);
    }
  print ("\n");
}

void
readhelpname (string0 buffer8, int *px8, char *f)
{
  register int i8r;

  //(*px8) = 1;
  skipbl (buffer8, &(*px8));
  (*px8) = (*px8) + 1;
  skipbl (buffer8, &(*px8));
  i8r = 1;
  if ((Member ((unsigned) (locase (buffer8.A[(*px8) - 1])), letterset.S)))
    {
      while ((Member
	      ((unsigned) (locase (buffer8.A[(*px8) - 1])), letterset.S))
	     && (i8r < MAXSTRING - 1))
	{
	  f[i8r - 1] = buffer8.A[(*px8) - 1];
	  (*px8) = (*px8) + 1;
	  i8r = i8r + 1;
	}
      f[i8r - 1] = '\0';
    }
  else
    {
      f[0] = 'x';
      f[1] = '\0';
    }
}

/* ----------------------------------------------------------------- */

void
helpcom (string0 buffer8, int px8)
{
  FILE *hlp;
  register int j8h;
  register int line;
  string0 file; 
  char tfile[MAXSTRING], nfile[MAXSTRING];
  char  qch;
  bool test, try;
/* ----------------------------------------------------------- */
  char  *tmpstr3;
  char pathvartmp[80];
  size_t plength;

/* ----------------------------------------------------------- */

  test = more;
  more = true;
  try = false;
  qch = ' ';
  memset (file.A, '\0', bcol);
  readhelpname (buffer8, &px8, tfile);

/* --------------------------------------------------------------- */
  tmpstr3 = pathvartmp;
  sprintf (tmpstr3, "%s/help/", dataPath);

  j8h = 0;
  do
    file.A[j8h++] = *tmpstr3++;
  while (*tmpstr3 != 0);
  plength = strlen (file.A);
  nfile[0] = '\0';
/* --------------------------------------------------------------- */
  if (interp (tfile, "remark", 3))
    (void) strcpy (nfile, "REMark");
  else if (interp (tfile, "abs", 3))
    (void) strcpy (nfile, "ABSoluteValue");
  else if (interp (tfile, "add", 3))
    (void) strcpy (nfile, "ADD");
  else if (interp (tfile, "all", 3))
    (void) strcpy (nfile, "ALLskewSfn");
  else if (interp (tfile, "associate", 5))
    (void) strcpy (nfile, "ASSOCiate");
  else if (interp (tfile, "at", 2))
    (void) strcpy (nfile, "ATtachPartitionToSfn");
  else if (interp (tfile, "au", 2))
    (void) strcpy (nfile, "AUtoOrIsoMorphism");

  else if (interp (tfile, "brm", 3))
    (void) strcpy (nfile, "BRMode");
  else if (interp (tfile, "br", 2))
    (void) strcpy (nfile, "BRanch");

  else if (interp (tfile, "can", 3))
    (void) strcpy (nfile, "CANcelDatFile");
  else if (interp (tfile, "casimirg", 8))
    (void) strcpy (nfile, "CASIMIRGnthTrace");
  else if (interp (tfile, "cas", 3))
    (void) strcpy (nfile, "CASimirNthordertrace");
  else if (interp (tfile, "ch_c", 4))
    (void) strcpy (nfile, "CH_CoeffsToOneForSfns");
  else if (interp (tfile, "ch_l", 4))
    (void) strcpy (nfile, "CH_LabelForOn");
  else if (interp (tfile, "ch_p", 4))
    (void) strcpy (nfile, "CH_PhaseOfSfns");
  else if (interp (tfile, "ch_s", 4))
    (void) strcpy (nfile, "CH_SpinIndex");
  else if (interp (tfile, "ch_u", 4))
    (void) strcpy (nfile, "CH_UoneReps");
  else if (interp (tfile, "character", 4))
    (void) strcpy (nfile, "CHARacter");
  else if (interp (tfile, "class", 5))
    (void) strcpy (nfile, "CLASS");
  else if (interp (tfile, "col", 3))
    (void) strcpy (nfile, "COLumns");
  else if (interp (tfile, "complement", 5))
    (void) strcpy (nfile, "COMPLement");
  else if (interp (tfile, "comp", 4))
    (void) strcpy (nfile, "COMPare");
  else if (interp (tfile, "contragredientrep", 7))
    (void) strcpy (nfile, "CONTRAGredientRep");
  else if (interp (tfile, "covariant", 3))
    (void) strcpy (nfile, "COVariant");
  else if (interp (tfile, "conjadd", 7))
    (void) strcpy (nfile, "CONJADD");
  else if (interp (tfile, "conj", 4))
    (void) strcpy (nfile, "CONJugateSfnList");
  else if (interp (tfile, "consplit", 8))
    (void) strcpy (nfile, "CONSPLIT");
  else if (interp (tfile, "content", 7))
    (void) strcpy (nfile, "CONTENT");
  else if (interp (tfile, "contractgroups", 4))
    (void) strcpy (nfile, "CONTractGroups");
  else if (interp (tfile, "conv_d_to_s", 11))
    (void) strcpy (nfile, "CONV_D_TO_Sfn");
  else if (interp (tfile, "conv_d_to_r", 11))
    (void) strcpy (nfile, "CONV_D_TO_Rep");
  else if (interp (tfile, "conv_r_to_s", 11))
    (void) strcpy (nfile, "CONV_R_TO_Sfn");
  else if (interp (tfile, "conv_s_to_r", 11))
    (void) strcpy (nfile, "CONV_S_TO_Rep");
  else if (interp (tfile, "countc", 6))
    (void) strcpy (nfile, "COUNTCoeffsInList");
  else if (interp (tfile, "countt", 6))
    (void) strcpy (nfile, "COUNTTermsInList");

  else if (interp (tfile, "DEAD", 4))
    (void) strcpy (nfile, "DEAD");
  else if (interp (tfile, "dim", 3))
    (void) strcpy (nfile, "DIMension");
  else if (interp (tfile, "dpm", 3))
    (void) strcpy (nfile, "DPMode");
  else if (interp (tfile, "dynkini", 7))
    (void) strcpy (nfile, "DYNKINIndex");
  else if (interp (tfile, "d_to_p", 6))
    (void) strcpy (nfile, "D_TO_Plabel");

  else if (interp (tfile, "end", 3))
    (void) strcpy (nfile, "END");
  else if (interp (tfile, "e_to_f", 6))
    (void) strcpy (nfile, "E_TO_FSymmFn");
  else if (interp (tfile, "e_to_h", 6))
    (void) strcpy (nfile, "E_TO_HSymmFn");
  else if (interp (tfile, "e_to_m", 6))
    (void) strcpy (nfile, "E_TO_MSymmFn");
  else if (interp (tfile, "e_to_s", 6))
    (void) strcpy (nfile, "E_TO_SSymmFn");
  else if (interp (tfile, "ent", 3))
    (void) strcpy (nfile, "ENTerVar");
  else if (interp (tfile, "eq", 2))
    (void) strcpy (nfile, "EQualSfnList");
  else if (interp (tfile, "exp", 3))
    (void) strcpy (nfile, "EXPandSfnList");
  else if (interp (tfile, "exit", 4))
    (void) strcpy (nfile, "EXITmode");

  else if (interp (tfile, "firstp", 6))
    (void) strcpy (nfile, "FIRSTPart");
  else if (interp (tfile, "fprod", 5))
    (void) strcpy (nfile, "FPROD");
  else if (interp (tfile, "ffprod", 6))
    (void) strcpy (nfile, "FFPROD");
  else if (interp (tfile, "fn", 2))
    (void) strcpy (nfile, "FN");
  else if (interp (tfile, "frobenius", 4))
    (void) strcpy (nfile, "FROBenius");
  if (interp (tfile, "f_to_e", 6))
    (void) strcpy (nfile, "F_TO_ESymmFn");
  else if (interp (tfile, "f_to_h", 6))
    (void) strcpy (nfile, "F_TO_HSymmFn");
  else if (interp (tfile, "f_to_m", 6))
    (void) strcpy (nfile, "F_TO_MSymmFn");
  else if (interp (tfile, "f_to_s", 6))
    (void) strcpy (nfile, "F_TO_SSymmFn");
  else if (interp (tfile, "fus", 3))
    (void) strcpy (nfile, "FUSion");

  else if (interp (tfile, "generic", 7))
    (void) strcpy (nfile, "GENERIC");
  else if (interp (tfile, "genprod", 3))
    (void) strcpy (nfile, "GENprod");
  else if (interp (tfile, "gr", 2))
    (void) strcpy (nfile, "GRoup");
  else if (interp (tfile, "gwt", 3))
    (void) strcpy (nfile, "GWT");

  else if (interp (tfile, "hall", 4))
    (void) strcpy (nfile, "HALLpolynomialProduct");
  else if (interp (tfile, "heap", 4))
    (void) strcpy (nfile, "HEAPstatus");
  else if (interp (tfile, "help", 4))
    (void) strcpy (nfile, "HELP");
  else if (interp (tfile, "h_to_e", 6))
    (void) strcpy (nfile, "H_TO_ESymmFn");
  else if (interp (tfile, "h_to_f", 6))
    (void) strcpy (nfile, "H_TO_FSymmFn");
  else if (interp (tfile, "h_to_m", 6))
    (void) strcpy (nfile, "H_TO_MSymmFn");
  else if (interp (tfile, "h_to_s", 6))
    (void) strcpy (nfile, "H_TO_SSymmFn");
  else if (interp (tfile, "hecke", 5))
    (void) strcpy (nfile, "HECKE");
  else if (interp (tfile, "hiveslrcoef", 11))
    (void) strcpy (nfile, "HIVESLRCOEFficient");
  else if (interp (tfile, "hstdlist", 5))
    (void) strcpy (nfile, "HSTDList");
  else if (interp (tfile, "hstd", 4))
    (void) strcpy (nfile, "HSTD");

  else if (interp (tfile, "inverseseries", 3))
    (void) strcpy (nfile, "INVerseseries");
  else if (interp (tfile, "i_sfnq", 6))
    (void) strcpy (nfile, "I_SFNQfnProduct");
  else if (interp (tfile, "i_pl", 4))
    (void) strcpy (nfile, "I_PLethysmRd");
  else if (interp (tfile, "i_q", 3))
    (void) strcpy (nfile, "I_QfnProduct");
  else if (interp (tfile, "index", 5))
    (void) strcpy (nfile, "INDEXsequence");
  else if (interp (tfile, "ins", 3))
    (void) strcpy (nfile, "INSertPartitionIntoSfn");
  else if (interp (tfile, "int", 3))
    (void) strcpy (nfile, "INTegerDivideCoeffs");
  else if (interp (tfile, "i", 1))
    (void) strcpy (nfile, "I_sfnProduct");

  else if (interp (tfile, "kinsert", 4))
    (void) strcpy (nfile, "KINSert");
  else if (interp (tfile, "kmatrix", 2))
    (void) strcpy (nfile, "KMatrix");
  else if (interp (tfile, "kostka", 1))
    (void) strcpy (nfile, "Kostka");

  else if (interp (tfile, "lab", 3))
    (void) strcpy (nfile, "LABel");
  else if (interp (tfile, "latticetest", 3))
    (void) strcpy (nfile, "LATticetest");
  else if (interp (tfile, "last", 4))
    (void) strcpy (nfile, "LASTresult");
  else if (interp (tfile, "lsequence", 4))
    (void) strcpy (nfile, "LSEQuence");
  else if (interp (tfile, "lines", 5))
    (void) strcpy (nfile, "LINES");
  else if (interp (tfile, "len", 3))
    (void) strcpy (nfile, "LENgthOfPartitionsSelect");
  else if (interp (tfile, "log", 3))
    (void) strcpy (nfile, "LOGfile");
  else if (interp (tfile, "lo", 2))
    (void) strcpy (nfile, "LOadFile");

  else if (interp (tfile, "macmixedseries", 4))
    (void) strcpy (nfile, "MACMixedSeries");
  else if (interp (tfile, "macseries", 3))
    (void) strcpy (nfile, "MACseries");
  else if (interp (tfile, "make", 4))
    (void) strcpy (nfile, "MAKEwtOfSfnToN");
  else if (interp (tfile, "maxc", 4))
    (void) strcpy (nfile, "MAXCoeffInList");
  else if (interp (tfile, "mix", 3))
    (void) strcpy (nfile, "MIXedTensorReps");
  else if (interp (tfile, "m_ti", 4))
    (void) strcpy (nfile, "M_TImesSfnProduct");
  else if (interp (tfile, "m_to_e", 6))
    (void) strcpy (nfile, "M_TO_ESymmFn");
  else if (interp (tfile, "m_to_f", 6))
    (void) strcpy (nfile, "M_TO_FSymmFn");
  else if (interp (tfile, "m_to_h", 6))
    (void) strcpy (nfile, "M_TO_HSymmFn");
  else if (interp (tfile, "m_to_s", 6))
    (void) strcpy (nfile, "M_TO_SSymmFn");
  else if (interp (tfile, "modsfn", 6))
    (void) strcpy (nfile, "MODSFN");
  else if (interp (tfile, "modspnr", 7))
    (void) strcpy (nfile, "MODSPNR");
  else if (interp (tfile, "mult_list", 6))
    (void) strcpy (nfile, "MULT_List");
  else if (interp (tfile, "mult_n", 6))
    (void) strcpy (nfile, "MULT_Ntimes");
  else if (interp (tfile, "mult_p", 6))
    (void) strcpy (nfile, "MULT_PartsByAnInt");
  else if (interp (tfile, "mult_sp", 7))
    (void) strcpy (nfile, "MULT_SPlitIntoTwoLists");
  else if (interp (tfile, "mult_s", 6))
    (void) strcpy (nfile, "MULT_SelectInList");
  else if (interp (tfile, "mult_u", 6))
    (void) strcpy (nfile, "MULT_UnitSfnIntTimes");
  else if (interp (tfile, "mu", 2))
    (void) strcpy (nfile, "MUlt_CoeffsByAnInt");
  else if (interp (tfile, "my", 2))
    (void) strcpy (nfile, "MYlistOfSfns");

  else if (interp (tfile, "nl", 2))
    (void) strcpy (nfile, "NLambda");
  else if (interp (tfile, "NSKew", 3))
    strcpy (nfile, "NSKew");
  else if (interp (tfile, "nstdise", 4))
    (void) strcpy (nfile, "NSTDise");

  else if (interp (tfile, "onsc", 4))
    (void) strcpy (nfile, "ONSCalar");
  else if (interp (tfile, "o_p", 3))
    (void) strcpy (nfile, "O_PfnProduct");
  else if (interp (tfile, "o_r", 3))
    (void) strcpy (nfile, "O_Restrict");
  else if (interp (tfile, "o_q", 3))
    (void) strcpy (nfile, "O_QfnProduct");
  else if (interp (tfile, "o", 1))
    (void) strcpy (nfile, "O_sfnProduct");

  else if (interp (tfile, "paritysequence", 6))
    (void) strcpy (nfile, "PARITYsequence");
  else if (interp (tfile, "pause", 5))
    (void) strcpy (nfile, "PAUSE");
  else if (interp (tfile, "plg", 3))
    (void) strcpy (nfile, "PLG");
  else if (interp (tfile, "p_to_d", 6))
    (void) strcpy (nfile, "P_TO_Dlabel");
  else if (interp (tfile, "pl", 2))
    (void) strcpy (nfile, "PLethysm");
  else if (interp (tfile, "p_to_s", 6))
    (void) strcpy (nfile, "P_TO_SSymmFn");
  else if (interp (tfile, "prop", 4))
    (void) strcpy (nfile, "PROPertyOfRepList");
  else if (interp (tfile, "p", 1))
    (void) strcpy (nfile, "ProductKronecker");

  else if (interp (tfile, "q_to_sd", 7))
    (void) strcpy (nfile, "Q_TO_SDual");
  else if (interp (tfile, "qexp", 4))
    (void) strcpy (nfile, "QEXPandSpecialSeries");
  else if (interp (tfile, "qqexp", 5))
    (void) strcpy (nfile, "QQEXPandSpecialSeries");
  else if (interp (tfile, "qqse", 4))
    (void) strcpy (nfile, "QQSEries");
  else if (interp (tfile, "qser", 4))
    (void) strcpy (nfile, "QSERies");
  else if (interp (tfile, "QUIT", 4))
    (void) strcpy (nfile, "QUIT");
  else if (interp (tfile, "racah", 5))
    (void) strcpy (nfile, "RACAHnotation");
  else if (interp (tfile, "fracah", 6))
    (void) strcpy (nfile, "FRACAHnotation");
  else if (interp (tfile, "raisei", 6))
    (void) strcpy (nfile, "RAISEInverseOp");
  else if (interp (tfile, "raise", 5))
    (void) strcpy (nfile, "RAISEop");
  else if (interp (tfile, "rd_raisei", 9))
    (void) strcpy (nfile, "RD_RAISEInverseOp");
  else if (interp (tfile, "rd_r", 4))
    (void) strcpy (nfile, "RD_RaiseOp");
  else if (interp (tfile, "rd_i_q", 6))
    (void) strcpy (nfile, "RD_I_QfnProduct");
  else if (interp (tfile, "rd_i", 4))
    (void) strcpy (nfile, "RD_I_sfnProduct");
  else if (interp (tfile, "readf", 5))
    (void) strcpy (nfile, "READFnFromDisk");
  else if (interp (tfile, "remark", 3))
    (void) strcpy (nfile, "REMark");
  else if (interp (tfile, "rep", 3))
    (void) strcpy (nfile, "REPmode");
  else if (interp (tfile, "return", 3))
    (void) strcpy (nfile, "RETurn");
  else if (interp (tfile, "rib_to_s", 8))
    (void) (strcpy (nfile, "RIB_TO_S"));
  else if (interp (tfile, "riemannl", 8))
    (void) strcpy (nfile, "RIEMANNList");
  else if (interp (tfile, "riemannp", 8))
    (void) strcpy (nfile, "RIEMANNPlethList");
  else if (interp (tfile, "riemanns", 8))
    (void) strcpy (nfile, "RIEMANNScalarsOrderN");
  else if (interp (tfile, "rm_evenparts", 8))
    (void) strcpy (nfile, "RM_EVENPARTS");
  else if (interp (tfile, "rm_evenr", 8))
    (void) strcpy (nfile, "RM_EVENRkSfnsOnly");
  else if (interp (tfile, "rm_evenw", 8))
    (void) strcpy (nfile, "RM_EVENWtInList");
  else if (interp (tfile, "rm_f", 4))
    (void) strcpy (nfile, "RM_FirstPartOfSfn");
  else if (interp (tfile, "rm_g", 4))
    (void) strcpy (nfile, "RM_Group");
  else if (interp (tfile, "rm_oddparts", 8))
    (void) strcpy (nfile, "RM_ODDPARTS");
  else if (interp (tfile, "rm_oddr", 7))
    (void) strcpy (nfile, "RM_ODDRkSfnsOnly");
  else if (interp (tfile, "rm_oddw", 7))
    (void) strcpy (nfile, "RM_ODDWtInList");
  else if (interp (tfile, "rm_parts", 8))
    (void) strcpy (nfile, "RM_PARTSequalN");
  else if (interp (tfile, "rm_nmparts", 6))
    (void) strcpy (nfile, "RM_NMParts");
  else if (interp (tfile, "rm_p", 4))
    (void) strcpy (nfile, "RM_PartitionFromSfn");
  else if (interp (tfile, "rm_r", 4))
    (void) strcpy (nfile, "RM_RepeatedPartsSfns");
  else if (interp (tfile, "rm_so", 5))
    (void) strcpy (nfile, "RM_SOnEvenLabels");
  else if (interp (tfile, "rm_u", 4))
    (void) strcpy (nfile, "RM_UoneWtOverMax");
  else
    /*if (interp(tfile, "router", 6))
       (void)strncpy(tfile.A, "ROUTER" + 0, sizeof(tfile.A));
       else */
  if (interp (tfile, "rp_f", 4))
    (void) strcpy (nfile, "RP_FirstPartBySpin");
  else if (interp (tfile, "rp_r", 4))
    (void) strcpy (nfile, "RP_RepOrSfnByWt");
  else if (interp (tfile, "rp_s", 4))
    (void) strcpy (nfile, "RP_SfnCoeffByInt");
  else if (interp (tfile, "rsame", 5))
    (void) strcpy (nfile, "RSAMEwtSfnList");
  else if (interp (tfile, "rule", 4))
    (void) strcpy (nfile, "RULE");
  else if (interp (tfile, "rv", 2))
    (void) strcpy (nfile, "RVar");
  else if (interp (tfile, "scalari", 7))
    (void) strcpy (nfile, "SCALARInner");
  else if (interp (tfile, "swapgroups", 4))
    (void) strcpy (nfile, "SWAPgroups");
  else if (interp (tfile, "sprch", 5))
    (void) strcpy (nfile, "SPRCH");
  else if (interp (tfile, "same", 4))
    (void) strcpy (nfile, "SAMEwtSfns");
  else if (interp (tfile, "save", 4))
    (void) strcpy (nfile, "SAVEsetVar");
  else if (interp (tfile, "sb_b", 4))
    (void) strcpy (nfile, "SB_Bell");
  else if (interp (tfile, "sb_con", 6))
    (void) strcpy (nfile, "SB_CONjecture");
  else if (interp (tfile, "sb_cut", 6))
    (void) strcpy (nfile, "SB_CUT");
  else if (interp (tfile, "sb_debug", 8))
    (void) strcpy (nfile, "SB_DEBUG");
  else if (interp (tfile, "sb_dim", 6))
    (void) strcpy (nfile, "SB_DIMension");
  else if (interp (tfile, "sb_d", 4))
    (void) strcpy (nfile, "SB_Digits");
  else if (interp (tfile, "sb_e", 4))
    (void) strcpy (nfile, "SB_Echo");
  else if (interp (tfile, "sb_list", 7))
    (void) strcpy (nfile, "SB_LISToutput");
  else if (interp (tfile, "sb_m", 4))
    (void) strcpy (nfile, "SB_More");
  else if (interp (tfile, "sb_progress", 7))
    (void) strcpy (nfile, "SB_PROGress");
  else if (interp (tfile, "sb_pow", 6))
    (void) strcpy (nfile, "SB_POWerNotation");
  else if (interp (tfile, "sb_q", 4))
    (void) strcpy (nfile, "SB_Qfn");
  else if (interp (tfile, "sb_rd", 5))
    (void) strcpy (nfile, "SB_RD_notation");
  else if (interp (tfile, "sb_psum", 7))
    (void) strcpy (nfile, "SB_PSUM");
  else if (interp (tfile, "sb_wprod", 4))
    (void) strcpy (nfile, "SB_WPROD");
  else if (interp (tfile, "sb_rev", 6))
    (void) strcpy (nfile, "SB_REVerseOrder");
  else if (interp (tfile, "sb_t", 4))
    (void) strcpy (nfile, "SB_TexOutPut");
  else if (interp (tfile, "serieste", 8))
    (void) strcpy (nfile, "SERIESTErmsThatSkew");
  else if (interp (tfile, "ser", 3))
    (void) strcpy (nfile, "SERiesToIntWt");
  else if (interp (tfile, "setf", 4))
    (void) strcpy (nfile, "SETFn");
  else if (interp (tfile, "setlimit", 6))
    (void) strcpy (nfile, "SETLIMit");
  else if (interp (tfile, "set_pwt", 5))
    (void) strcpy (nfile, "SET_PWT");
  else if (interp (tfile, "setv", 4))
    (void) strcpy (nfile, "SETVarInDPmode");
  else if (interp (tfile, "setrv", 5))
    (void) strcpy (nfile, "SETRVar");
  else if (interp (tfile, "setsv", 5))
    (void) strcpy (nfile, "SETSVar");
  else if (interp (tfile, "sfn", 3))
    (void) strcpy (nfile, "SFNmode");
  else if (interp (tfile, "signseq", 7))
    (void) strcpy (nfile, "SIGNSEQuence");
  else if (interp (tfile, "s_to_e", 6))
    (void) strcpy (nfile, "S_TO_ESymmFn");
  else if (interp (tfile, "s_to_f", 6))
    (void) strcpy (nfile, "S_TO_FSymmFn");
  else if (interp (tfile, "s_to_h", 6))
    (void) strcpy (nfile, "S_TO_HSymmFn");
  else if (interp (tfile, "s_to_m", 6))
    (void) strcpy (nfile, "S_TO_MSymmFn");
  else if (interp (tfile, "s_to_q", 6))
    (void) strcpy (nfile, "S_TO_QsymmFn");
  else if (interp (tfile, "s_to_p", 6))
    (void) strcpy (nfile, "S_TO_PsymmFn");
  else if (interp (tfile, "sk_p", 4))
    (void) strcpy (nfile, "SK_Pfn");
  else if (interp (tfile, "sk_q", 4))
    (void) strcpy (nfile, "SK_Qfn");
  else if (interp (tfile, "sk", 2))
    (void) strcpy (nfile, "SK_sfn");
  else if (interp (tfile, "schar", 5))
    (void) strcpy (nfile, "SCHAR");
  else if (interp (tfile, "snchar", 2))
    (void) strcpy (nfile, "SNchar");
  else if (interp (tfile, "smon", 4))
    (void) strcpy (nfile, "SMON");
  else if (interp (tfile, "sponmodify", 5))
    (void) strcpy (nfile, "SPONModify");
  else if (interp (tfile, "sprextend", 5))
    (void) strcpy (nfile, "SPREXtend");
  else if (interp (tfile, "spl", 3))
    (void) strcpy (nfile, "SPLitIntoSpinAndTensor");
  else if (interp (tfile, "spin", 4))
    (void) strcpy (nfile, "SPIN");
  else if (interp (tfile, "status", 4))
    (void) strcpy (nfile, "STATusOfSchur");
  else if (interp (tfile, "spstar", 6))
    (void) strcpy (nfile, "SPSTAR");
  else if (interp (tfile, "squares", 2))
    (void) strcpy (nfile, "SQuares");
  else if (interp (tfile, "std_o", 5))
    (void) strcpy (nfile, "STD_OneDprep");
  else if (interp (tfile, "std_q", 5))
    (void) strcpy (nfile, "STD_Qfn");
  else if (interp (tfile, "std", 3))
    (void) strcpy (nfile, "STD");
  else if (interp (tfile, "STARequivalent", 4))
    (void) strcpy (nfile, "STARequivalent");
  else if (interp (tfile, "stop", 4))
    (void) strcpy (nfile, "STOP");
  else if (interp (tfile, "SUB", 3))
    (void) strcpy (nfile, "SUBtract");
  else if (interp (tfile, "sumsquares", 5))
    (void) strcpy (nfile, "SUMSQuares");
  else if (interp (tfile, "sum", 3))
    (void) strcpy (nfile, "SUM");
  else if (interp (tfile, "sup", 3))
    (void) strcpy (nfile, "SUPpressOutputToScreen");
  else if (interp (tfile, "sv", 2))
    (void) strcpy (nfile, "SVar");
  else if (interp (tfile, "tab", 3))
    (void) strcpy (nfile, "TABleOfBranchingRules");
  else if (interp (tfile, "uonea", 5))
    (void) strcpy (nfile, "UONEAddInteger");
  else if (interp (tfile, "uoned", 5))
    (void) strcpy (nfile, "UONEDivInteger");
  else if (interp (tfile, "uonet", 5))
    (void) strcpy (nfile, "UONETrace");
  else if (interp (tfile, "vmult", 2))
    (void) strcpy (nfile, "VMult");
  else if (interp (tfile, "v", 1))
    (void) strcpy (nfile, "VarForDpreps");
  else if (interp (tfile, "wsequence", 4))
    (void) strcpy (nfile, "WSEQuence");
  else if (interp (tfile, "wtofreporsfnselect", 2))
    (void) strcpy (nfile, "WTofRepOrSfnSelect");
  else if (interp (tfile, "whatg", 5))
    (void) strcpy (nfile, "WHATGroup");
  else if (interp (tfile, "with", 4))
    (void) strcpy (nfile, "WITH");
  else if (interp (tfile, "wrfntod", 7))
    (void) strcpy (nfile, "WRFNTODisk");
  else if (interp (tfile, "wr", 2))
    (void) strcpy (nfile, "WRfnToScreen");
  else if (interp (tfile, "yhooklengths", 2))
    (void) strcpy (nfile, "YHooklengths");
  else if (interp (tfile, "youngdiagrams", 2))
    (void) strcpy (nfile, "YOungDiagrams");
  else if (interp (tfile, "yshapeselect", 2))
    (void) strcpy (nfile, "YShapeSelect");
  else if (interp (tfile, "z", 1))
    (void) strcpy (nfile, "Zero");
  else if (interp (tfile, "x", 1))
    strcpy (nfile, "HELP");
  //try = true;

  if (nfile[0] == '\0')
    {
      strcpyuppercase (nfile, tfile);	/* try the file with identical name but in uppercase letters */
    }
  strncat (file.A, nfile, MAXSTRING - plength);

  if (try)
    /*(void)fprintf(output.fp, "ERROR ? must be followed by letters\n"), Putl(output, 1); */
    (void) fprintf (output.fp, "`%.80s` not found.\n", file.A);
  if (findfile (file) && !try)
    {
      hlp = fopen (file.A, "r");
      if (hlp == NULL)
	{
	  strerror (errno);
	}
      else
	{
	  qch = ' ';
	  fflush (stdin);
	  line = 1;
	  while (!feof (hlp) && (qch != 'q'))
	    {
	      if (fgets (tfile, MAXSTRING, hlp) != NULL)
		{
		  fputs (tfile, output.fp);
		  if (logging)
		    fputs (tfile, logfile.fp);
		  line = line + 1;
		  if ((line == tlines) && more)
		    {
		      qch = more_Question();
		      /*while (qqch.A[0] == '\n' && !feof (hlp))	// Return: one more line to display
			{
			  if (fgets (tfile, MAXSTRING, hlp) != NULL)
			    {
			      tfile[strlen (tfile) - 1] = '\0';	// no \n
			      fprintf (output.fp, "%s", tfile);
			    }
			  qqch.A[0] = getchar ();
			}
		      qch = locase (qqch.A[1 - 1]);*/
		      if (qch == '\n')
			 line = tlines-1;
		      else
		         line = 1;
		    }
		}
	    }
	  fprintf (output.fp, "\n");
	}
      fclose (hlp);
    }
  else
    (void) fprintf (output.fp, "`%.80s` not found.\n", file.A);
  more = test;			/*11/12/95 */
}				/*11/12/95 */
