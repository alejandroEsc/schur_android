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

#include <stdbool.h>
#include "dim.h"
#include "define.h"
#include "type.h"
#include "var.h"
#include "r.h"
#include "utils.h"
#include "write.h"

#define RESTORE_CURSOR_POSITION()   printf("[17D[1A                     [21D")

/** return the frame in s in display format (TeX or not) and the length of the
 * resulting string */
int
wrtfrme2 (char *s, frame x)
{
  int i, k, m, t;
  register int j;
  char *saveptr = s;

  t = maxdim - 1;
  while ((x.A[t] == 0) && (t > 0))
    t = t - 1;
  x.length = t;
  if (!pow_note) {
    i = 2;
    s += sprintf (s, "%1d", x.A[1]);	// kind of strcat
    if (x.A[1] > 9) {
      if (tex)
	s += sprintf (s, "\\");
      s += sprintf (s, " ");
    }
    while ((i <= t)) {
      s += sprintf (s, "%1d", x.A[i]);
      if (abs (x.A[i]) > 9) {
	if (tex)
	  s += sprintf (s, "\\");
	s += sprintf (s, " ");
      }
      i = i + 1;
    }
  } else {			//pow_note
    i = 1;
    m = 1;
    while ((m <= t) && (i < maxdim - 1)) {
      j = x.A[i];
      k = 1;
      do {
	i = i + 1;
	if (x.A[i] == j)
	  k = k + 1;
      } while ((x.A[i] == j) && (i < maxdim - 1));
      if (k > cutoff) {
	s += sprintf (s, "%d", x.A[m]);
	if (abs (x.A[m]) > 9) {
	  if (tex)
	    s += sprintf (s, "\\");
	  s += sprintf (s, " ");
	}
	if (tex) {
	  s += sprintf (s, "^{%d}", k);
	  if (i <= t)
	    s += sprintf (s, "\\ ");
	} else {
	  s += sprintf (s, "^%1d", k);
	  if (i <= t)
	    s += sprintf (s, " ");
	}
      } else {
	for (j = 1; j <= k; j++) {
	  if (abs (x.A[m]) > 9) {
	    if (tex)
	      s += sprintf (s, "%d\\ ", x.A[m]);
	    else
	      s += sprintf (s, "%d ", x.A[m]);
	    if ((x.A[m] < 0))
	      if (tex)
		s += sprintf (s, "\\ ");
	      else
		s += sprintf (s, " ");
	  } else
	    s += sprintf (s, "%1d", x.A[m]);
	}
      }
      m = i;
    }				// while
    if (i == 1)
      s += sprintf (s, "0");
  }
  return (s - saveptr);
}

int
multToString (char *s, int mult, bool start)
{
  int i = 0, c = 0;
  if (!start) {
    if (tex)
      c = sprintf (s, "\\ ");
    else
      c = sprintf (s, " ");
    i = c;
    s += c;
  }

  if (mult <= -1) {
    if (tex)
      c = sprintf (s, "-\\ ");
    else
      c = sprintf (s, "- ");
    s += c;
    i += c;
  }

  if (!start && mult >= 1) {
    if (tex)
      c = sprintf (s, "+\\ ");
    else
      c = sprintf (s, "+ ");
    s += c;
    i += c;
  }
  if (abs (mult) > 1) {
    c = sprintf (s, "%d", abs (mult));
    s += c;
    i += c;
  }

  return (i);
}

void
wrttlst2 (text * fyle, int *qq, char a, char b, termptr list, bool vdu)
{
  bool startx, protecttex = tex;
  int line = 1;
  int curcol = 0;
  char key = ' ';
  termptr ptr;
  char intertotal[500];
  char *inter;

  if (list == NULL && !qtest) {
    fprintf (fyle->fp, "zero\n");
    return;
  }
  if (qfn || pfn || mmono || forg || homo || psum || elem)
    protecttex = false;		/* do not protect {} from TeX interpretation */

  startx = true;
  while ((list != NULL) && (key != 'q')) {
    ptr = list;
    curcol++;
    inter = intertotal;
    if (tex)
      inter += sprintf (inter, "$");
    inter += multToString (inter, ptr->mult, startx);

    if (startx)
      startx = false;
    if ((qfn == true) && (pfn == false))
      inter += sprintf (inter, "Q_");
    else if ((pfn == true) && (qfn == true))
      inter += sprintf (inter, "P_");
    else if ((mmono == true))
      inter += sprintf (inter, "m_");
    else if ((forg == true))
      inter += sprintf (inter, "f_");
    else if ((homo == true))
      inter += sprintf (inter, "h_");
    else if ((psum == true))
      inter += sprintf (inter, "p_");
    else if ((elem == true))
      inter += sprintf (inter, "e_");

    if (protecttex && a == '{')
      inter += sprintf (inter, "\\");
    inter += sprintf (inter, "%c", a);
    inter += wrtfrme2 (inter, ptr->val);
    startx = false;
    if (protecttex && b == '}')
      inter += sprintf (inter, "\\");
    inter += sprintf (inter, "%c", b);
    if ((ptr->slab == '#')) {
      if (tex)
	inter += sprintf (inter, "_");
      inter += sprintf (inter, "#");
    }
    if (tex) {
      inter += sprintf (inter, "$");
      if (curcol % columns == 0)
	inter += sprintf (inter, "\\\\\n");
      else if (ptr->next != NULL)
	inter += sprintf (inter, "&");
    }

    if (!tex && (int) (*qq + strlen (intertotal)) > (int) (tcol)) {
      fprintf (fyle->fp, "\n");
      *qq = 0;
      if (vdu && more) {
	line = line + 1;

	if (line == tlines - 1) {
	  if ((key = more_Question ()) == 'q')
	    return;
	  else {
	    if (key == '\n')	// just one more line
	      line = tlines - 2;
	    else
	      line = 1;
	  }
	}

      }
    }
    *qq += fprintf (fyle->fp, "%s", intertotal);
    list = ptr->next;
  }
}


void
writer2 (text * fyle, int *qq, char a, char b, ocharptr list, bool vdu)
{
  char key, spinchar;
  bool startx, protecttex;
  int i, line, k, m, z = 0, curcol = 0;	/* curcol is used for tex output */
  register int l;
  ocharptr ptrlist;
  caystype caystmp;
  char intertotal[500];
  char *s;

  line = 1;
  key = '?';
  spinchar = 's';
  if (list == NULL) {
    *qq += fprintf ((*fyle).fp, "zero");
    return;
  }
  startx = true;
  protecttex = false;
  ptrlist = list;

  while ((ptrlist != NULL) && (key != esc)) {
    s = intertotal;
    curcol++;

    if (tex)
      fprintf ((*fyle).fp, "$");
    s += multToString (s, ptrlist->mult, startx);

    if (startx)
      startx = false;
    switch ((int) (currgrp.A[1 - 1].name)) {
    case sn:
      if (nreduce) {
	a = '<';
	b = '>';
      } else {
	a = '{';
	b = '}';
	protecttex = true;
      }
      break;
    case an:
      if (nreduce) {
	a = '<';
	b = '>';
      } else {
	a = '{';
	b = '}';
	protecttex = true;
      }
      break;
    case sung:
    case un:
    case unm:
    case sunm:
    case unc:
      a = '{';
      b = '}';
      protecttex = true;
      break;
    case son:
    case on:
    case sonc:
      a = '[';
      b = ']';
      break;
    case spn:
      a = '<';
      b = '>';
      break;
    case spnc:
    case mp:
      a = '<';
      b = '>';
      break;
    case ospnm:
      a = '[';
      b = '>';
      break;
    case e6:
    case e7:
    case e8:
    case f4:
    case g2:
    case en:
    case l168:
      a = '(';
      b = ')';
      break;
    case nill:
      warn ("no group set;", cr);
      break;
    default:
      Caseerror (Line);
    }
    if (currgrp.A[1 - 1].name == e6)
      caystmp = pe6;
    else
      caystmp = norm;

    if (tex && protecttex)
      s += sprintf (s, "\\");
    s += sprintf (s, "%c", a);
    if ((currgrp.A[1 - 1].name == unc))
      s += sprintf (s, "(");
    if (ptrlist->spin) {
      s += sprintf (s, "%c", spinchar);
      if (((currgrp.A[1 - 1].name != spnc)
	   && (currgrp.A[1 - 1].name != mp)
	   && (currgrp.A[1 - 1].name != sonc)))
	s += sprintf (s, ";");
    }
    if ((ptrlist->C6_double && (currgrp.A[1 - 1].name != spnc)
	 && (currgrp.A[1 - 1].name != mp)
	 && (currgrp.A[1 - 1].name != sonc))) {
      s += wrtfrme2 (s, ptrlist->conval);
      s += sprintf (s, ";");
    }
    if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == mp)
	|| (currgrp.A[1 - 1].name == sonc))
      s += sprintf (s, "%d(", ptrlist->val.A[maxdim]);

    if (caystmp == pe6) {
      s += sprintf (s, "%d:", ptrlist->val.A[1]);
      if (ptrlist->val.A[2] == 0) {
	s += sprintf (s, "0");
	i = 3;
      } else
	i = 2;
      if (!pow_note) {
	while (ptrlist->val.A[i] != 0) {
	  s += sprintf (s, "%1d", ptrlist->val.A[i]);
	  i = i + 1;
	}
      } else {
	if ((group == g2))
	  if ((racah == true)) {
	    z = ptrlist->val.A[1];
	    ptrlist->val.A[1] = z - ptrlist->val.A[2];
	  }
	m = i;
	while (ptrlist->val.A[m] != 0) {
	  l = ptrlist->val.A[i];
	  i = i + 1;
	  k = 1;
	  while (ptrlist->val.A[i] == l) {
	    i = i + 1;
	    k = k + 1;
	  }
	  if (k > cutoff) {
	    if (ptrlist->val.A[m] > 9) {
	      if (tex)
		s += sprintf (s, "\\");
	      s += sprintf (s, " ");
	    }
	    s += sprintf (s, "%1d^%1d", ptrlist->val.A[m], k);
	  } else
	    for (l = 1; l <= k; l++) {
	      if (ptrlist->val.A[m] > 9)
		s +=
		  sprintf (s, "%d%s", ptrlist->val.A[m], tex ? "\\ " : " ");
	      else
		s += sprintf (s, "%d", ptrlist->val.A[m]);
	    }
	  m = i;
	}
	if ((group == g2))
	  if ((racah == true))
	    ptrlist->val.A[1] = z;
      }
    } else if ((currgrp.A[1 - 1].name != spnc)
	       && (currgrp.A[1 - 1].name != mp)
	       && (currgrp.A[1 - 1].name != sonc))
      s += wrtfrme2 (s, ptrlist->val);
    if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == mp)) {
      s += wrtfrme2 (s, ptrlist->val);
      s += sprintf (s, ")>");
    }
    if ((currgrp.A[1 - 1].name == spnc)) {
      if (ptrlist->lab != ' ') {
	if (tex)
	  s += sprintf (s, "_");
	s += sprintf (s, "#");
      }
    }
    if ((currgrp.A[1 - 1].name == sonc)) {
      s += wrtfrme2 (s, ptrlist->val);
      s += sprintf (s, ")]");
    }
    if ((currgrp.A[1 - 1].name == unc))
      s += sprintf (s, ")");
    if ((currgrp.A[1 - 1].name != spnc) && (currgrp.A[1 - 1].name != mp)
	&& (currgrp.A[1 - 1].name != sonc)) {
      if (tex && protecttex)
	s += sprintf (s, "\\");
      s += sprintf (s, "%c", b);
      if (ptrlist->lab != ' ') {
	if (tex)
	  s += sprintf (s, "_");
	s += sprintf (s, "%c", ptrlist->lab);
      }
    }

    if (tex) {
      s += sprintf (s, "$");
      if (curcol % columns == 0)
	s += sprintf (s, "\\\\\n");
      else if (ptrlist->next != NULL)
	s += sprintf (s, "&");
    } else if ((int) (*qq + strlen (intertotal)) > (int) (tcol)) {
      fprintf (fyle->fp, "\n");
      *qq = 0;
      if (vdu && more) {
	line = line + 1;

	if (line == tlines - 1) {
	  if ((key = more_Question ()) == 'q')
	    return;
	  else {
	    if (key == '\n')	// just one more line
	      line = tlines - 2;
	    else
	      line = 1;
	  }
	}

      }
    }
    *qq += fprintf (fyle->fp, "%s", intertotal);

    ptrlist = ptrlist->next;
  }
}

void
writeprod2 (text * fyle, int *qq, prodtype pr, bool vdu)
{
  bool startx, protecttex;
  int z = 0, i, line, k, m, curcol = 0;
  register int l;
  register int j;
  caystype cays;
  prodtype ptr;
  char a = '{', b = '}', spinchar, key;
  char intertotal[500];
  char *s;

  line = 1;
  key = '?';
  spinchar = 's';
  if (pr == NULL) {
    *qq += fprintf ((*fyle).fp, "zero");
    return;
  }
  startx = true;
  protecttex = false;
  ptr = pr;
  while ((ptr != NULL) && (key != esc)) {
    s = intertotal;
    curcol++;

    if (tex)
      fprintf ((*fyle).fp, "$");
    s += multToString (s, ptr->mult, startx);
    if (startx)
      startx = false;

    for (j = 1; j <= nprod; j++) {
      switch ((int) (currgrp.A[j - 1].name)) {
      case sn:
	if (nreduce) {
	  a = '<';
	  b = '>';
	} else {
	  a = '{';
	  b = '}';
	  protecttex = true;
	}
	break;
      case an:
	if (nreduce) {
	  a = '<';
	  b = '>';
	} else {
	  a = '{';
	  b = '}';
	  protecttex = true;
	}
	break;
      case sung:
      case un:
      case unm:
      case sunm:
      case unc:
	a = '{';
	b = '}';
	protecttex = true;
	break;
      case son:
      case on:
      case sonc:
	a = '[';
	b = ']';
	break;
      case spn:
	a = '<';
	b = '>';
	break;
      case spnc:
      case mp:
	a = '<';
	b = '>';
	//spinc = ptr->prods.A[j - 1]->spin;
	break;
      case ospnm:
	a = '[';
	b = '>';
	break;
      case e6:
      case e7:
      case e8:
      case f4:
      case g2:
      case en:
      case l168:
	a = '(';
	b = ')';
	break;
      case nill:
	warn ("no group set;", cr);
	break;
      default:
	Caseerror (Line);
      }
      if (currgrp.A[j - 1].name == e6)
	cays = pe6;
      else
	cays = norm;
      if (tex && protecttex)
	s += sprintf (s, "\\%c", a);
      else
	s += sprintf (s, "%c", a);
      if (currgrp.A[j - 1].name == unc)
	s += sprintf (s, "(");
      if (ptr->prods.A[j - 1]->spin) {
	s += sprintf (s, "%c", spinchar);
	if (((currgrp.A[j - 1].name != spnc)
	     && (currgrp.A[1 - 1].name != mp)
	     && (currgrp.A[1 - 1].name != sonc)))
	  s += sprintf (s, ";");
      }
      if (ptr->prods.A[j - 1]->C6_double && (currgrp.A[j - 1].name != spnc)
	  && (currgrp.A[1 - 1].name != mp)
	  && (currgrp.A[1 - 1].name != sonc)) {
	s += wrtfrme2 (s, ptr->prods.A[j - 1]->conval);
	s += sprintf (s, ";");
      }
      if (((currgrp.A[j - 1].name == spnc)
	   || (currgrp.A[1 - 1].name == mp)
	   || (currgrp.A[1 - 1].name == sonc))) {
	s += sprintf (s, "%d(", ptr->prods.A[j - 1]->val.A[maxdim]);
      }
      if (cays == pe6) {
	s += sprintf (s, "%d:", ptr->prods.A[j - 1]->val.A[1]);
	if (ptr->prods.A[j - 1]->val.A[2] == 0) {
	  s += sprintf (s, "0");
	  i = 3;
	} else
	  i = 2;
	if (!pow_note) {
	  while (ptr->prods.A[j - 1]->val.A[i] != 0) {
	    s += sprintf (s, "%d", ptr->prods.A[j - 1]->val.A[i]);
	    i = i + 1;
	  }
	} else {
	  register ocharptr W33 = ptr->prods.A[j - 1];

	  if ((group == g2))
	    if ((racah == true)) {
	      z = W33->val.A[1];
	      W33->val.A[1] = z - W33->val.A[2];
	    }
	  m = i;
	  while ((W33->val.A[m] != 0)) {
	    l = W33->val.A[i];
	    i = i + 1;
	    k = 1;
	    while (W33->val.A[i] == l) {
	      i = i + 1;
	      k = k + 1;
	    }
	    if (k > cutoff) {
	      s += sprintf (s, "%d", W33->val.A[m]);
	      if (W33->val.A[m] > 9) {
		if (tex)
		  s += sprintf (s, "\\ ");
		else
		  s += sprintf (s, " ");
	      }
	      s += sprintf (s, "^%d", k);
	    } else
	      for (l = 1; l <= k; l++) {
		if (W33->val.A[m] > 9) {
		  if (tex)
		    s += sprintf (s, "%d\\ ", W33->val.A[m]);
		  else
		    s += sprintf (s, "%d ", W33->val.A[m]);

		} else
		  s += sprintf (s, "%d", W33->val.A[m]);
	      }
	    m = i;
	  }
	  if ((group == g2))
	    if ((racah == true))
	      W33->val.A[1] = z;
	}
      } else if (((currgrp.A[j - 1].name != spnc)
		  && (currgrp.A[j - 1].name != mp)
		  && (currgrp.A[1 - 1].name != sonc)))
	s += wrtfrme2 (s, ptr->prods.A[j - 1]->val);
      if ((currgrp.A[j - 1].name == spnc)
	  || (currgrp.A[j - 1].name == mp)) {
	s += wrtfrme2 (s, ptr->prods.A[j - 1]->val);
	s += sprintf (s, ")>");
      }
      if ((currgrp.A[1 - 1].name == sonc)) {
	s += wrtfrme2 (s, ptr->prods.A[j - 1]->val);
	s += sprintf (s, ")]");
      }
      if (currgrp.A[j - 1].name == unc)
	s += sprintf (s, ")");
      if (((currgrp.A[j - 1].name != spnc)
	   && (currgrp.A[j - 1].name != mp)
	   && (currgrp.A[1 - 1].name != sonc))) {
	if (tex && protecttex)
	  s += sprintf (s, "\\");
	s += sprintf (s, "%c", b);
	if (ptr->prods.A[j - 1]->lab != ' ') {
	  if (tex)
	    s += sprintf (s, "_");
	  s += sprintf (s, "%c", ptr->prods.A[j - 1]->lab);
	}
      }
    }

    if (tex) {
      s += sprintf (s, "$");
      if (curcol % columns == 0)
	s += sprintf (s, "\\\\\n");
      else if (ptr->next != NULL)
	s += sprintf (s, "&");
    } else if ((int) (*qq + strlen (intertotal)) > (int) (tcol)) {
      fprintf (fyle->fp, "\n");
      *qq = 0;
      if (vdu && more) {
	line = line + 1;

	if (line == tlines - 1) {
	  if ((key = more_Question ()) == 'q')
	    return;
	  else {
	    if (key == '\n')	// just one more line
	      line = tlines - 2;
	    else
	      line = 1;
	  }
	}

      }
    }
    *qq += fprintf (fyle->fp, "%s", intertotal);

    ptr = ptr->next;
  }
}

void
writepoly (lbframe apoly)
{
  register int p;
  bool notfirst = false;	// if the poly is equals to 0, notfirst will stay at false.
  char intertotal[maxdim], *s;	// used for the display width mgmt.
  char var = 'q';
  int index;

  //Putchr ('\n', output);

  s = intertotal;
  if (apoly.A[0] != 0) {	// degree 0
    if (apoly.A[0] == 1)
      s += sprintf (s, "1");
    else
      s += sprintf (s, "%d", apoly.A[0]);
    notfirst = true;
  }

  if (apoly.A[1] != 0) {	// degree 1
    if (notfirst) {
      if (tex)
	s += sprintf (s, "\\");
      if (apoly.A[1] == 1) {
	s += sprintf (s, " +%c", var);
      } else if (apoly.A[1] == -1) {
	s += sprintf (s, " -%c", var);
      } else {
	s += sprintf (s, " %+d%c", apoly.A[1], var);
      }
    } else {			// first term
      if (apoly.A[1] == 1)
	s += sprintf (s, "%c", var);
      else if (apoly.A[1] == -1)
	s += sprintf (s, "-%c", var);
      else
	s += sprintf (s, "%d%c", apoly.A[1], var);
    }
    notfirst = true;
  }

  print ("%s", intertotal);
  index = strlen (intertotal);
  for (p = 2; p <= setlimit; p++)	// degrees>1
  {
    s = intertotal;

    if (apoly.A[p] != 0) {
      if (notfirst) {
	if (tex)
	  s += sprintf (s, "\\");
	if (apoly.A[p] == 1)
	  s += sprintf (s, " +");
	else if (apoly.A[p] == -1)
	  s += sprintf (s, " -");
	else
	  s += sprintf (s, " %+d", apoly.A[p]);
      } else {
	if (apoly.A[p] == 1);
	else if (apoly.A[p] == -1)
	  s += sprintf (s, "-");
	else
	  s += sprintf (s, "%d", apoly.A[p]);
      }
      if (tex)
	s += sprintf (s, "%c^{%d}", var, p);
      else
	s += sprintf (s, "%c^%d", var, p);
      notfirst = true;
      index += strlen (intertotal);
      if (index > tcol) {
	print ("\n%s", intertotal);
	index = strlen (intertotal);
      } else
	print ("%s", intertotal);
    }
  }
  if (!notfirst)		// the poly is equals to 0
    print ("0");

  print ("\n");
}
