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

// export a variable or a list in a file, List format.
// this function is called if sb_ListOutput is true
# include <stdio.h>

# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
# include "var.h"

// for len and qlen
# include "s1.h"

// for skipbl, readword, interp
# include "r.h"

// for error, Caserror,...
# include "utils.h"

# include "sets_mgmt.h"
# include "ReadWrite.h"
# include "savecomList.h"

void
saveframeList (text * t, frame fr)
{
  register int i, lastStep;

  //lastStep = MAX(len (fr),fr.length);
  //if (lastStep>maxdim)
  //  lastStep=maxdim;
  lastStep = len (&fr);
  if (lastStep==0)
	  lastStep=1;
  fprintf ((*t).fp, "[");
  for (i = 1; i <= lastStep; i++)
    {
      if (i > 1)
	fprintf ((*t).fp, ",");
      fprintf ((*t).fp, "%d", fr.A[i]);
    }
  fprintf ((*t).fp, "]");
}

char
charac_code (void)
{
  if ((qfn == true) && (pfn == false))
    return ('Q');
  else if ((pfn == true) && (qfn == true))
    return ('P');
  else if ((mmono == true))
    return ('m');
  else if ((forg == true))
    return ('f');
  else if ((homo == true))
    return ('h');
  else if ((psum == true))
    return ('p');
  else if ((elem == true))
    return ('e');
  else
    return (' ');
}

void
outputGroupCode (FILE * f, groopArray grp)
{
  int j;
  for (j = 1; j <= nprod; j++)
    {
      if ((nprod > 1) && (j <= nprod) && (j > 1))
	fprintf (f, "*");
      switch ((int) (grp.A[j - 1].name))
	{
	case sung:
	case sunm:
	  fprintf (f, "SU,");
	  break;
	case un:
	case unm:
	case unc:
	  fprintf (f, "U,");
	  break;
	case son:
	  fprintf (f, "SO,");
	  break;
	case on:
	  fprintf (f, "O,");
	  break;
	case spn:
	case spnc:
	  fprintf (f, "Sp,");
	  break;
	case ospnm:
	  fprintf (f, "OSp,");
	  break;
	case sonc:
	  fprintf (f, "SO*,");
	  break;
	case mp:
	  fprintf (f, "Mp,");
	  break;
	case l168:
	  fprintf (f, "L,");
	  break;
	case sn:
	  fprintf (f, "S,");
	  if (((bool) ((grp.A[j - 1].rank) & 1)))
	    qsn = 1;
	  else
	    qsn = 0;
	  break;
	case an:
	  fprintf (f, "A,");
	  if (((bool) ((grp.A[j - 1].rank) & 1)))
	    qsn = 1;
	  else
	    qsn = 0;
	  break;
	case e6:
	case e7:
	case e8:
	case en:
	  fprintf (f, "E,");
	  break;
	case f4:
	  fprintf (f, "F,");
	  break;
	case g2:
	  fprintf (f, "G,");
	  break;
	case nill:
	  inform ("No group set;", cr);
	  break;
	default:
	  Caseerror (Line);
	}
      if (grp.A[j - 1].name != nill)
	{
	  fprintf (f, "%d", grp.A[j - 1].rank);
	  switch ((int) (grp.A[j - 1].name))
	    {
	    case ospnm:
	    case unm:
	    case sunm:
	      fprintf (f, ",%d", grp.A[j - 1].rank2);
	      break;
	    case unc:
	      fprintf (f, ",%d", grp.A[j - 1].rank2);
	      break;
	    case spnc:
	      fprintf (f, ",R");
	      break;
	    case un:
	    case on:
	    case son:
	    case spn:
	    case sung:
	    case mp:
	    case an:
	    case sn:
	    case g2:
	    case f4:
	    case e6:
	    case e7:
	    case en:
	    case e8:
	    case nill:
	    case sonc:
	    case l168:
	      break;
	    default:
	      Caseerror (Line);
	    }
	}
    }
}				// end outputGroupCode

void
savecomList (string0 bufferx, int px)
{
  bool svs, justone;
  int theone;
  register int i;
  termptr ptr;
  ocharptr cptr;
  char wd[MAXSTRING], c, hit;
  string0 filen;
  text t;

  t.init = 0;
  t.fp = NULL;
  t.eoln = 0;
  t.buf = '\0';
  hit = skipbl (bufferx, &px);
  readword (&bufferx, &px, wd);
  svs = interp ("svariable", wd, 2);
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
	      fprintf (output.fp, "save[");
	      if (svs)
		{
		  fprintf (t.fp, "schur:=[[sfn]");
		  if ((c = charac_code ()) != ' ')
		    fprintf (t.fp, ",%c", c);	// format : [[sfn], optional_code, var1, var2,...]
		  else
	 	    fprintf (t.fp, ",s");
		  for (i = 1; i <= svarlimit; i++)
		    {
		      if ((!justone && (svar.A[i - 1] != NULL))
			  || (justone && i == theone))
			{
			  fprintf (t.fp, ",[");	// format of var_i : [term1,term2,...]
			  ptr = svar.A[i - 1];
			  fprintf (output.fp, "%d", i);	// display on screen the variable currently saved
			  while (ptr != NULL)	// format of term_i : [mult,[partition],optional_#]
			    {
			      fprintf (t.fp, "[%d,", ptr->mult);
			      saveframeList (&t, ptr->val);
			      if (ptr->slab == '#')
				fprintf (t.fp, ",\"#\"");
			      fprintf (t.fp, "]");	// end of term_i
			      ptr = ptr->next;
			      if (ptr != NULL)	// at least one more term
				fprintf (t.fp, ",");
			    }
			  fprintf (t.fp, "]");	// end of var_i
			}
		    }
		  fprintf (t.fp, "];\n");	// end of list of vars.
		}
	      else
		{		// rvar.. format : [[rep, groupe], var1, var2, ...]
		  fprintf (t.fp, "schur:=[[rep,");
		  outputGroupCode (t.fp, currgrp);
		  fprintf (t.fp, "]");
		  for (i = 1; i <= rvarlimit; i++)
		    {
		      if ((!justone && (vari.A[i - 1] != NULL))
			  || (justone && i == theone))
			{
			  cptr = vari.A[i - 1];
			  fprintf (output.fp, "%d", i);		// on screen
			  // format var_i : [term1,term2,...]
			  fprintf (t.fp, ",[");
			  // format term_i [mult, spin, conval_list, conlab, val_list, lab]
			  while (cptr != NULL)
			    {
			      fprintf (t.fp, "[%d,%d", cptr->mult, (unsigned) (cptr->spin));
			      if ( ! cptr->C6_double)
				      cptr->conlab = spc;

			  fprintf (t.fp, ",");
			  saveframeList (&t, cptr->conval);
			  fprintf (t.fp, ",\"%c\"", cptr->conlab);
			      fprintf (t.fp, ",");
			      saveframeList (&t, cptr->val);
			      fprintf (t.fp, ",\"%c\"]", cptr->lab);
			      cptr = cptr->next;
			      if (cptr != NULL)
				fprintf (t.fp, ",");
			    }
			  fprintf (t.fp, "]");	// end of var_i
			}
		    }
		  fprintf (t.fp, "];\n");	// end of vars.
		}
	    }
	  fclose (t.fp);
	  fprintf (output.fp, "]\n");	// on screen
	}
    }
}
