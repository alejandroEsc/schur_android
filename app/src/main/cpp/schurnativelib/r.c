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
/** \file r.c
 * Definitions for i/o
*/
#include <stdio.h>
#include <limits.h>

#include <unistd.h>
#include <stdbool.h>

#include "define.h"
#include "mymalloc.h"
#include "sets_mgmt.h"
#include "Scanck.h"
#include "../config.h"

#ifdef HAVE_LIBREADLINE
#include <readline/readline.h>
#include <readline/history.h>

bool use_readline = true;

#else

#define readline(m)   m
#define add_history(m)
bool use_readline = false;
#endif

#define RESTORE_CURSOR_POSITION()   printf("[17D[1A    ")

static FILE *Tmpfil;

/*
**	Start of program definitions
*/
#include "dim.h"
#include "type.h"
#include "var.h"
#include "utils.h"
#include "m.h"
#include "init.h"
#include "s1.h"
#include "s2.h"
#include "r.h"

caystype *G144_cays;
bool *G142_texx;
caystype *G140_cays;

/** Ask for more output or not */
char
more_Question (void)
{
  char key;

  printf ("%s--More [q,m,Spacebar,Return]--%s", colorsBegin, colorsEnd);
  fflush (stdout);
  key = getchar ();
  if (key != '\n')
    getchar ();			//drop the final \n
  fflush (stdin);

//up one line, spaces to clean line and back cursor to beginning of the line
  printf ("[1A[2K");
  return (key);
}

/** skip brackets and parenthesis () {} <> */
char
skipbl (string0 b, int *p)
{
  while (b.A[*p - 1] != '\0' && Member ((unsigned) b.A[*p - 1], brackets.S)
	 && *p < bcol)
    *p += 1;
  if (b.A[*p - 1] == '\0' || (*p == bcol && b.A[*p - 1] == ' '))
    return (cr);
  else
    return (locase (b.A[*p - 1]));
}

void
readfilename (string0 bufferx, int *px, string0 * f)
{
  char hit;
  register int i;
  for (i = 1; i <= bcol; i++) {
    f->A[i - 1] = ' ';
  }
  hit = skipbl (bufferx, &(*px));
  if (hit != qt)
    error (MISSING_QUOTE, (*px));
  else {
    (*px) = (*px) + 1;
    skipbl (bufferx, &(*px));
    i = 1;
    while ((bufferx.A[(*px) - 1] != qt) && ((*px) < pcol)) {
      f->A[i - 1] = locase (bufferx.A[(*px) - 1]);
      (*px) = (*px) + 1;
      i = i + 1;
    }
    if ((f->A[1 - 1] == 'p') && (f->A[2 - 1] == 'r')
	&& (f->A[3 - 1] == 'n'))
      f->A[4 - 1] = '.';
    if (bufferx.A[(*px) - 1] != qt)
      error (MISSING_QUOTE, (*px));
    else
      f->A[i - 1] = 0;
  }
}

/** Ask y/n question  */
bool
dothis (void)
{
  string0 ch;
  bool yn;
  register int i;
  for (i = 1; i <= bcol; i++) {
    ch.A[i - 1] = ' ';
  }
  fprintf (output.fp, "(y/n)?"), Putl (output, 0);
  do {
    Fscan (input), Scan ("%s", ch.A), Getx (input), Getl (&input);
    ch.A[1 - 1] = locase (ch.A[1 - 1]);
  }
  while (ch.A[0] != '+' && ch.A[0] != '-' && ch.A[0] != 'n'
	 && ch.A[0] != 'y');

  yn = (ch.A[0] == '+') || (ch.A[0] == 'y');
  if (yn)
    Putchr ('y', output), Putchr ('\n', output);
  else
    Putchr ('n', output), Putchr ('\n', output);
  return yn;
}


/** Ask if file must be overwrited or not  */
bool
overwrite (string0 l)
{
  if (findfile (l)) {
    fprintf (output.fp, "file already exists, ok to over-write"),
      Putl (output, 0);
    return dothis ();
  } else
    return (true);
}

bool
findfile (string0 name)
{
  FILE *f;
  if ((f = fopen (name.A, "r")) == NULL)
    return (false);
  else {
    fclose (f);
    return (true);
  }
}

termptr
reverselist (termptr l)
{
  register termptr R134;
  termptr q, r;

  q = NULL;
  while (l != NULL) {
    r = l->next;
    l->next = q;
    q = l;
    l = r;
  }
  R134 = q;
  return R134;
}

ocharptr
rreverselist (ocharptr l)
{
  register ocharptr R135;
  ocharptr q, r;

  q = NULL;
  while (l != NULL) {
    r = l->next;
    l->next = q;
    q = l;
    l = r;
  }
  R135 = q;
  return R135;
}

prodtype
preverselist (prodtype l)
{
  register prodtype R136;
  prodtype q, r;

  q = NULL;
  while (l != NULL) {
    r = l->next;
    l->next = q;
    q = l;
    l = r;
  }
  R136 = q;
  return R136;
}

// called in sfnmode to read a line into buffx
// from fyle (stdin or a file on disk), p is the pointer inside buffx
// if fnex=true, the line comes from a user-defined function
void
readacard (text * fyle, string0 * buffx, int *p)
{
  char a, *ptr;
  int col;
  register int i;

  *p = 1;

  memset (buffx->A, ' ', bcol);
  //use_readline = (use_readline && isatty (fileno (stdin)));
  use_readline = (use_readline && isatty (STDIN_FILENO));

  i = 1;
  if (fnex && prompt[0] == '\0')	//we are reading a user-defined function and not enter
  {
    if (fnptr != NULL) {
      memcpy (buffx->A, fnptr->gbuff.A, bcol);
      fnptr = fnptr->next;
      if (fnptr == NULL)	// bring back normal output for the last line of the function
	iosup = false;
      if (debug_schur && fnex) {
	char word[bcol];
	strncpy (word, buffx->A, bcol);
	word[bcol - 1] = '\0';
	fprintf (stderr, "fnex %s\n", word);
      }
    } else			// fnptr==NULL
    {
      //fprintf (stderr, "Error in function : no continuation\n");
      fnex = false;
      iosup = false;
      return;
    }
  } else if ((fnex && prompt[0] != '\0') || fyle->fp == stdin) {
    char *resu;
    if (use_readline) {
      if (prompt[0] != '\0')
	ptr = readline (prompt);
      else
	ptr = readline (displayModes[mode]);
      if (ptr != NULL && strlen (ptr) > bcol) {
	add_history (ptr);
	error (LINE_TOO_LONG, bcol);
	return;
      }
    } else {
      if (prompt != NULL)
	printf ("%s", prompt);
      else
	printf ("%s", displayModes[mode]);
      ptr = (char *) malloc (sizeof (char) * bcol + 1);
      resu = fgets (ptr, bcol + 1, stdin);
      if (strlen (ptr) > bcol) {
	error (LINE_TOO_LONG, bcol);
	ptr[bcol - 1] = '\0';
      }
      if (resu != NULL && ptr[strlen (ptr) - 1] == '\n') {
	ptr[strlen (ptr) - 1] = '\0';
      }
      if (resu == NULL)
	ptr = NULL;
    }
    if (ptr != NULL && *ptr != '\0') {
      if (use_readline)
	add_history (ptr);
      if (strlen (ptr) > bcol) {
	error (LINE_TOO_LONG, bcol);
      }
      strcpy (buffx->A, stripwhite (ptr));
      buffx->A[strlen (buffx->A)] = ' ';	// strange but needed
    }

    if (ptr == NULL)		// Ctrl-D = end.
      strcpy (buffx->A, "end ");
    if (ptr != NULL)
      free (ptr);

  } else {			// not stdin, must be a user-defined function on disk
    do {
      a = Getchr (*fyle);
      if ((a == backspace) || ((i == 1) && (a == ' ')) || (a == 10)
	  || (a == 13))
	if (i > 2)
	  i = i - 1;
	else
	  i = 1;
      else {
	buffx->A[i - 1] = a;
	i = i + 1;
      }
    }
    while (!feof (fyle->fp) && !(Eoln (*fyle) || (i >= bcol)));
    if (i == bcol && !Eoln (*fyle) && a != ' ')
      error (LINE_TOO_LONG, bcol);
    buffx->A[i] = '\0';
    if (feof (fyle->fp))
      strcpy (buffx->A, "end ");
    //buffx->A[0] = ' ';
  }

  if (echo && logging) {
    Putchr ('\n', logfile);
    fprintf (logfile.fp, "->"), Putl (logfile, 0);
    if (bcol < pcol - 2)
      col = bcol;
    else
      col = pcol - 2;
    for (i = 1; i <= col; i++) {
      Putchr (buffx->A[i - 1], logfile);
    }
    Putchr ('\n', logfile);
  }
}				// readacard

void
readint (string0 * buffx, int *p, int *res)
{
  bool neg, bad;

  *res = 0;
  neg = false;
  bad = false;
  while (Member (buffx->A[(*p) - 1], brackets.S) && (*p < bcol))
    *p += 1;
  if (*p < bcol) {
    if ((buffx->A[*p - 1] == '-') || (buffx->A[*p - 1] == '+')) {
      neg = (buffx->A[*p - 1] == '-');
      *p = *p + 1;
      while (Member (buffx->A[*p - 1], brackets.S) && (*p < bcol))
	*p = *p + 1;
    }
    if (*p < bcol) {
      if (Member (buffx->A[*p - 1], numbers.S)) {
	do {
	  bad = *res > (INT_MAX / 10);
	  if (!bad)
	    *res = *res * 10 + charval (buffx->A[*p - 1]);
	  else
	    error (MISSING_INTEGER, *p);
	  *p = *p + 1;
	}
	while (!(!(Member ((unsigned) (buffx->A[(*p) - 1]), numbers.S))
		 || ((*p) >= bcol) || bad));
	if (neg)
	  *res = -*res;
      }
    }
  }
  if (buffx->A[*p - 1] == ',')
    *p += 1;
}


/** read an elt of a partition */
void
readelt (string0 * buffx, int *p, int *elt)
{
  int s;

  while ((Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))
	 && ((*p) < bcol))
    (*p) = (*p) + 1;
  if (buffx->A[(*p) - 1] == '~') {
    s = -1;
    (*p) = (*p) + 1;
  } else
    s = 1;
  if ((buffx->A[(*p) - 1] == '!') || !digits) {
    if (buffx->A[(*p) - 1] == '!')
      (*p) = (*p) + 1;
    readint (&(*buffx), &(*p), &(*elt));
  } else {
    (*elt) = charval (buffx->A[(*p) - 1]);
    (*p) = (*p) + 1;
  }
  (*elt) = (*elt) * s;
}

void
readachrc (text * fyle, string0 * buffx, int *p, ocharptr * list)
{
  char hit, ch;
  int q, coeff;
  ocharptr lastptr, dum;
  frame tframe;

  lastptr = NULL;
  do {
    hit = skipbl ((*buffx), &(*p));
    if ((*p) < bcol) {
      cnu (&(*list));
      {
	register ocharptr tempo = *list;

	tempo->C6_double = false;
	tempo->spin = false;
	tempo->conlab = ' ';
	tempo->lab = ' ';
	tempo->mult = 1;
	tempo->conval = nolls;
	tempo->val = nolls;
	if (((buffx->A[(*p) - 1] == '+') || (buffx->A[(*p) - 1] == '-'))) {
	  if (buffx->A[(*p) - 1] == '-')
	    tempo->mult = -1;
	  (*p) = (*p) + 1;
	  hit = skipbl ((*buffx), &(*p));
	}
	if (buffx->A[(*p) - 1] == cont)
	  readacard (&(*fyle), &(*buffx), &(*p));
	if (locase (buffx->A[(*p) - 1]) != 's') {
	  q = (*p);
	  readelt (&(*buffx), &q, &coeff);
	  if (((buffx->A[q - 1] == '.') || (buffx->A[q - 1] == 's')
	       || (buffx->A[q - 1] == 'S'))) {
	    tempo->mult = tempo->mult * coeff;
	    if (buffx->A[q - 1] == '.') {
	      (*p) = q + 1;
	      hit = skipbl ((*buffx), &(*p));
	    } else
	      (*p) = q;
	  }
	}
	if (((buffx->A[(*p) - 1] == 's') || (buffx->A[(*p) - 1] == 'S'))) {
	  tempo->spin = true;
	  (*p) = (*p) + 1;
	} else
	  tempo->spin = false;
	readpart (&(*buffx), &(*p), &tempo->val);
	if (buffx->A[(*p) - 1] == ';') {
	  tempo->C6_double = true;
	  (*p) = (*p) + 1;
	  tframe = tempo->val;
	  readpart (&(*buffx), &(*p), &tempo->val);
	  tempo->conval = tframe;
	  tempo->val.A[maxdim] = tempo->conval.A[1];
	}
	ch = buffx->A[*p - 1];
	tempo->lab = ' ';
	//if (Member ((unsigned) (ch), Conset[5]))
	if (ch == '#' || ch == '+' || ch == '-') {
	  q = *p + 1;
	  hit = skipbl (*buffx, &q);
	  if ((buffx->A[q - 1] == delim) || (q == bcol)
	      || ((hit == '+') || (hit == '-') || (hit == '#'))) {
	    tempo->lab = ch;
	    *p = q;
	  }
	}
	tempo->next = lastptr;
	lastptr = *list;
	hit = skipbl (*buffx, p);
      }
    }
  }
  while (!(((buffx->A[(*p) - 1] != '+') && (buffx->A[(*p) - 1] != '-')
	    && (buffx->A[(*p) - 1] != '#')) || ((*p) > bcol - 1)));
  hit = skipbl ((*buffx), &(*p));
  if (buffx->A[(*p) - 1] == delim)
    (*p) = (*p) + 1;
  (*list) = NULL;
  while (lastptr != NULL) {
    dum = lastptr->next;
    lastptr->next = (*list);
    (*list) = lastptr;
    lastptr = dum;
  }
}

void
readpart (string0 * buffx, int *p, frame * partr)
{
  int j;
  register int k;
  register int i;
  bool test;
  *partr = nolls;
  for (i = 1; i <= maxdim; i++)
    partr->A[i] = 0;
  while ((Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))
	 && ((*p) < bcol))
    (*p) = (*p) + 1;
  /*if ((Member
     ((unsigned) (buffx->A[(*p) - 1]), (Union (numbers.S, Conset[8]))))
     && ((*p) < bcol)) */
  if (Member (buffx->A[*p - 1], numbersEtc.S) && (*p < bcol)) {
    i = 0;
    test = false;
    do {
      i = i + 1;
      readelt (&(*buffx), &(*p), &j);
      if ((j <= maxdim) && (j >= -maxdim))
	partr->A[i] = j;
      else {			// added by FB 
	/*fprintf (stderr,
	   "Parts greater than %d or lower than %d are ignored (%d)\n",
	   maxdim, -maxdim, j); */
	j = 1;
	test = true;
      }
      while ((Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))
	     && ((*p) < bcol))
	(*p) = (*p) + 1;
      if (buffx->A[(*p) - 1] == '^') {
	(*p) = (*p) + 1;
	readelt (&(*buffx), &(*p), &j);
	for (k = 1; k <= j - 1; k++) {
	  if ((i + k) < maxdim)
	    partr->A[i + k] = partr->A[i];
	  else
	    test = true;
	}
	i = i + j - 1;
	while ((Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))
	       && ((*p) < bcol))
	  (*p) = (*p) + 1;
      }
    }
    /*while (!
       ((!(Member
       ((unsigned) (buffx->A[(*p) - 1]),
       (Union (numbers.S, Conset[9]))))) || (i == maxdim)
       || test)); */
    while (Member (buffx->A[*p - 1], numbersEtc.S) && (i < maxdim) && !test);

    Claimset ();
    /*if (test
       ||
       (Member
       ((unsigned) (buffx->A[(*p) - 1]),
       (Union (numbers.S, Conset[10]))))) */
    if (test || Member (buffx->A[*p - 1], numbersEtc.S)) {
      fprintf (output.fp, "ERROR maxdim = %d exceeded\n",
	       maxdim), Putl (output, 1);
      while (buffx->A[*p - 1] != '+' && buffx->A[*p - 1] != '-' && *p < bcol)
	*p = *p + 1;
    }
    Claimset ();
  }
  Claimset ();
  partr->length = i;
}

void
readchrc (text * fyle, string0 * buffx, int *p, ocharptr * list)
{
  char ch;
  char hit;
  int q, coeff, savep, initp = *p - 1;
  ocharptr lastptr, dum;
  frame tframe;

  lastptr = NULL;
  do {
    hit = skipbl ((*buffx), &(*p));
    if ((*p) < bcol) {
      cnu (&(*list));
      {
	register ocharptr tempo = *list;

	tempo->C6_double = false;
	tempo->spin = false;
	tempo->conlab = ' ';
	tempo->lab = ' ';
	tempo->mult = 1;
	tempo->conval = nolls;
	tempo->val = nolls;
	if (((buffx->A[(*p) - 1] == '+') || (buffx->A[(*p) - 1] == '-'))) {
	  if (buffx->A[(*p) - 1] == '-')
	    tempo->mult = -1;
	  (*p) = (*p) + 1;
	  hit = skipbl ((*buffx), &(*p));
	}
	if (hit != cr) {
	  if (buffx->A[(*p) - 1] == cont)
	    readacard (&(*fyle), &(*buffx), &(*p));
	  if (locase (buffx->A[(*p) - 1]) != 's') {
	    q = (*p);
	    readelt (&(*buffx), &q, &coeff);
	    if (((buffx->A[q - 1] == '.') || (buffx->A[q - 1] == 'S')
		 || (buffx->A[q - 1] == 's'))) {
	      tempo->mult = tempo->mult * coeff;
	      if (buffx->A[q - 1] == '.') {
		(*p) = q + 1;
		hit = skipbl ((*buffx), &(*p));
	      } else
		(*p) = q;
	    }
	  }
	  if (((buffx->A[(*p) - 1] == 's')
	       || (buffx->A[(*p) - 1] == 'S'))) {
	    tempo->spin = true;
	    (*p) = (*p) + 1;
	  } else
	    tempo->spin = false;
	  savep = *p - 1;
	  readpart (&(*buffx), &(*p), &(tempo->val));
	  if (buffx->A[(*p) - 1] == ';') {
	    tempo->C6_double = true;
	    (*p) = (*p) + 1;
	    tframe = tempo->val;
	    savep = *p - 1;
	    readpart (buffx, p, &tempo->val);
	    tempo->conval = tframe;
	    tempo->val.A[maxdim] = tempo->conval.A[1];
	  } else if (currgrp.A[1 - 1].name == spnc) {
	    char str[bcol];
	    int i;
	    for (i = 0; i <= savep - initp && i < bcol - 1; i++)
	      str[i] = buffx->A[i + initp];
	    str[i] = '\0';
	    printf
	      ("Warning: Ambiguous, a part of the input needs a ';' : %s\n",
	       str);
	  }

	  ch = buffx->A[(*p) - 1];
	  tempo->lab = ' ';
	  if (((ch == '+') || (ch == '-') || (ch == '#'))) {
	    q = (*p) + 1;
	    hit = skipbl ((*buffx), &q);
	    if ((buffx->A[q - 1] == delim) || (q == bcol)
		|| ((hit == '+') || (hit == '-') || (hit == '#'))) {
	      tempo->lab = ch;
	      (*p) = q;
	    }
	  }
	  tempo->next = lastptr;
	  lastptr = *list;
	  hit = skipbl ((*buffx), &(*p));
	} else {
	  error (PARENTHESIS_NOT_NEEDED, 1);
	  dispchr (list);
	}
      }
    }
  }
  while (!(((buffx->A[(*p) - 1] != '+') && (buffx->A[(*p) - 1] != '-')
	    && (buffx->A[(*p) - 1] != '#')) || ((*p) > bcol - 1)));

  hit = skipbl (*buffx, p);
  if (buffx->A[(*p) - 1] == delim)
    (*p) = (*p) + 1;
  (*list) = NULL;
  while (lastptr != NULL) {
    dum = lastptr->next;
    lastptr->next = *list;
    *list = lastptr;
    lastptr = dum;
  }
}

void
readword (string0 * buffx, int *p, char *word)
{
  char hit;
  int i;
  bool test;

  while ((Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))
	 && ((*p) < bcol))
    (*p) = (*p) + 1;

  memset (word, ' ', bcol);	// initialize
  word[bcol - 1] = '\0';

  i = 1;
  test = false;
  if ((*p) < bcol) {
    do {
      if (!(Member ((unsigned) (buffx->A[(*p) - 1]), brackets.S))) {
	if (buffx->A[(*p) - 1] == '?') {
	  word[i - 1] = '?';
	  test = true;
	} else if (((test == true) && (buffx->A[(*p) - 1] != '?')))
	  word[i - 1] = buffx->A[(*p) - 1];
	else if (!test)
	  word[i - 1] = locase (buffx->A[(*p) - 1]);
	i = i + 1;
      }
      (*p) = (*p) + 1;
    }
    while (!(((*p) >= bcol) || (i > 24)
	     ||
	     !(Member
	       ((unsigned) (locase (buffx->A[(*p) - 1])), letterset.S))));
    hit = skipbl ((*buffx), &(*p));
    if (hit == delim)
      (*p) = (*p) + 1;
  }
}

void
readlist (text * fyle, string0 * buffx, int *p, termptr * list, bool tidy)
{
  termptr lastptr, dum;
  int coeff;
  int q;
  char ch, hit;

  lastptr = NULL;
  do {
    hit = skipbl ((*buffx), &(*p));
    if ((*p) < bcol) {
      snu (&(*list));
      {
	register termptr tempo = &(*(*list));
	for (q = 1; q <= maxdim; q++) {
	  tempo->val.A[q] = 0;
	}
	tempo->mult = 1;
	if (((buffx->A[(*p) - 1] == '+') || (buffx->A[(*p) - 1] == '-'))) {
	  if (buffx->A[(*p) - 1] == '-')
	    tempo->mult = -1;
	  (*p) = (*p) + 1;
	  hit = skipbl ((*buffx), &(*p));
	}
	if (buffx->A[(*p) - 1] == cont)
	  readacard (&(*fyle), &(*buffx), &(*p));

	q = (*p);
	readelt (&(*buffx), &q, &coeff);
	if (buffx->A[q - 1] == '.') {
	  tempo->mult = tempo->mult * coeff;
	  (*p) = q + 1;
	}
	readpart (&(*buffx), &(*p), &tempo->val);
	ch = buffx->A[(*p) - 1];
	tempo->slab = ' ';
	if (ch == '#') {
	  q = (*p) + 1;
	  hit = skipbl ((*buffx), &q);
	  if (buffx->A[q - 1] == delim || q == bcol
	      || (hit == '+' || hit == '-')) {
	    tempo->slab = ch;
	    (*p) = q;
	  }
	}
	tempo->next = lastptr;
	lastptr = (*list);
      }
      hit = skipbl ((*buffx), &(*p));
    }
  }
  while (!((buffx->A[*p - 1] != '+' && buffx->A[(*p) - 1] != '-')
	   || *p >= bcol - 1));
  if (buffx->A[(*p) - 1] == delim)
    *p += 1;
  *list = NULL;
  while (lastptr != NULL) {
    dum = lastptr->next;
    lastptr->next = *list;
    *list = lastptr;
    lastptr = dum;
  }
  if (tidy) {
    if (qfn == true)
      qstndise (list);
    else
      stndise (list);
    sort (list, false);
  }
}

/*
void
wrtfrme (text * fyle, frame x)
{
  int i, k, m, t;
  register int j;

  t = maxdim - 2;
  while ((x.A[t] == 0) && (t > 0))
    t = t - 1;
  if (!pow_note) {
    i = 2;
    (void) fprintf ((*fyle).fp, "%1d", x.A[1]), Putl ((*fyle), 0);
    if (((x.A[1] > 9) && tex))
      (void) fprintf ((*fyle).fp, "\\"), Putl ((*fyle), 0);
    else
      Putchr (' ', (*fyle));
    while ((i <= t)) {
      (void) fprintf ((*fyle).fp, "%1d", x.A[i]), Putl ((*fyle), 0);
      if (((abs (x.A[i]) > 9) && tex))
	(void) fprintf ((*fyle).fp, "\\"), Putl ((*fyle), 0);
      else
	Putchr (' ', (*fyle));
      i = i + 1;
    }
  } else {
    i = 1;
    m = 1;
    while ((m <= t) && (i < maxdim - 2)) {
      j = x.A[i];
      k = 1;
      do {
	i = i + 1;
	if ((x.A[i] == j))
	  k = k + 1;
      }
      while (!((x.A[i] != j) || (i == maxdim - 2)));
      if (k > cutoff) {
	if (((m > 1) && (abs (x.A[m]) > 9)))
	  if (tex)
	    (void) fprintf ((*fyle).fp, "\\"), Putl ((*fyle), 0);
	  else
	    Putchr (' ', (*fyle));
	if (tex) {
	  (void) fprintf ((*fyle).fp, "%1d^", x.A[m]), Putl ((*fyle), 0);
	  if ((k > 9))
	    (void) fprintf ((*fyle).fp, "{%1d}\\", k), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, "%1d", k), Putl ((*fyle), 0);
	} else {
	  (void) fprintf ((*fyle).fp, "%1d^%1d", x.A[m], k);
	  Putl ((*fyle), 0);
	  if (i <= t)
	    fprintf ((*fyle).fp, " ");	// added by FB 
	}
      } else
	for (j = 1; j <= k; j++) {
	  if (abs (x.A[m]) > 9) {
	    if (tex)
	      (void) fprintf ((*fyle).fp, "%1d\\", x.A[m]), Putl ((*fyle), 0);
	    else
	      (void) fprintf ((*fyle).fp, "%2d ", x.A[m]), Putl ((*fyle), 0);
	    if ((x.A[m] < 0))
	      if (tex)
		(void) fprintf ((*fyle).fp, "\\"), Putl ((*fyle), 0);
	      else
		Putchr (' ', (*fyle));
	  } else
	    (void) fprintf ((*fyle).fp, "%1d", x.A[m]), Putl ((*fyle), 0);
	}
      m = i;
    }
    if (i == 1)
      Putchr ('0', (*fyle));
  }
}
*/

int
csize (ocharptr list)
{
  register int R137;
  int i, j, k, l, m, n, jm, abmult;

  (*G142_texx) = false;
  switch ((int) (currgrp.A[1 - 1].name)) {
  case sn:
  case sung:
  case un:
  case unm:
  case sunm:
  case unc:
    (*G142_texx) = true;
    break;
  case on:
  case an:
  case ospnm:
  case son:
  case spn:
  case mp:
  case spnc:
  case sonc:
  case g2:
  case f4:
  case e6:
  case e7:
  case e8:
  case l168:
    break;
  default:
    Caseerror (Line);
  }

  {
    register ocharptr W21 = &(*list);

    if (!pow_note) {
      i = 1;
      if (abs (W21->val.A[1]) > 9)
	if (abs (W21->val.A[1]) > 99)
	  j = 3;
	else
	  j = 2;
      else
	j = 0;
      if (W21->val.A[1] < 0)
	j = j + 1;
      if (((*G140_cays) == e6bra) || ((*G140_cays) == pe6)
	  || ((*G140_cays) == triple)) {
	if (W21->val.A[2] > 9)
	  if (W21->val.A[2] > 99)
	    j = j + 3;
	  else
	    j = j + 2;
	i = 2;
      }
      while (W21->val.A[i + 1] != 0) {
	i = i + 1;
	if (W21->val.A[i] > 9)
	  if (W21->val.A[i] > 99)
	    j = j + 3;
	  else
	    j = j + 2;
      }
    } else {
      if (((*G140_cays) == e6bra) || ((*G140_cays) == pe6)
	  || ((*G140_cays) == triple)) {
	if (abs (W21->val.A[1]) > 9)
	  if (abs (W21->val.A[1]) > 99)
	    j = 3;
	  else
	    j = 2;
	else
	  j = 0;
	if (W21->val.A[1] < 0)
	  j = j + 1;
	if (W21->val.A[2] > 9)
	  if (W21->val.A[2] > 99)
	    j = j + 3;
	  else
	    j = j + 2;
	i = 3;
	m = 3;
      } else {
	i = 1;
	m = 1;
	j = 0;
      }
      while ((W21->val.A[m] != 0) && (m < maxdim)) {
	l = W21->val.A[i];
	i = i + 1;
	k = 1;
	while ((W21->val.A[i] == l) && (i < maxdim)) {
	  i = i + 1;
	  k = k + 1;
	}
	if (W21->val.A[m] > 9)
	  if (W21->val.A[m] > 99)
	    n = 3;
	  else
	    n = 2;
	else
	  n = 1;
	if (k > cutoff) {
	  if (n > 1)
	    j = n + j + 3;
	  else
	    j = n + j + 2;
	  if (k > 9)
	    j = j + 2;
	  else
	    j = j + 1;
	} else
	  j = j + k * n;
	m = i;
      }
      if (i != 1)
	i = 0;
    }
    k = i + j;
    switch ((int) ((*G140_cays))) {
    case norm:
      break;
    case pe6:
      k = k + 1;
      break;
    case e6bra:
      k = k + 2;
      break;
    case prod:
      k = k + 2;
      break;
    case triple:
      k = k + 4;
      break;
    default:
      Caseerror (Line);
    }
    if (W21->spin) {
      k = k + 2;
      if ((*G140_cays) == prod)
	k = k + 2;
    }
    if ((tex && W21->spin))
      k = k + 5;
    if ((tex && (*G142_texx)))
      k = k + 2;
    if (W21->lab != ' ')
      k = k + 1;
    if (((W21->lab != ' ') && tex))
      k = k + 1;
    if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == mp)
	|| (currgrp.A[1 - 1].name == sonc))
      k = k + 3;
    if (W21->C6_double) {
      if (((*G140_cays) == prod) && (W21->conlab != ' '))
	k = k + 1;
      i = 1;
      if (W21->conval.A[1] > 9)
	if (W21->conval.A[1] > 99)
	  j = 4;
	else
	  j = 3;
      else
	j = 1;
      while (W21->conval.A[i + 1] != 0) {
	i = i + 1;
	if (W21->conval.A[i] > 9)
	  if (W21->conval.A[i] > 99)
	    j = j + 3;
	  else
	    j = j + 2;
      }
      k = k + i + j;
    }
    abmult = abs (W21->mult);
    jm = 1;
    do {
      abmult = abmult / 10;
      if (abmult > 0)
	jm = jm + 1;
    }
    while (!(abmult == 0));
    if (abs (W21->mult) == 1)
      jm = jm - 1;
    R137 = k + jm + 5;
  }
  return R137;
}

/*
void
writer (text * fyle, int *qq, int col, char a, char b, caystype cays,
	ocharptr list, bool vdu)
{
  char key, spinchar;
  bool startx, texx;
  int i, cs, line, k, m, z = 0,	// modified by FB, was unitialized
    tcolumns, sline, abmult;
  register int l;
  ocharptr top;
  caystype *F141;
  bool *F143;

  F143 = G142_texx;
  G142_texx = &texx;
  F141 = G140_cays;
  G140_cays = &cays;
  line = 1;
  key = '?';
  spinchar = 's';
  sline = 1;
  tcolumns = 0;
  (void) fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
  (*qq) = 4;
  top = list;
  if (list == NULL)
    {
      (void) fprintf ((*fyle).fp, "zero"), Putl ((*fyle), 0);
      (*qq) = (*qq) + 4;
    }
  startx = true;
  while ((list != NULL) && (key != esc))
    {
      register ocharptr W22 = &(*list);

      tcolumns = tcolumns + 1;
      cs = 0;
      if (currgrp.A[1 - 1].name == e6)
	(*G140_cays) = pe6;
      else
	(*G140_cays) = norm;
      cs = cs + csize (list);
      abmult = abs (W22->mult);
      cs = cs + 4;
      do
	{
	  abmult = abmult / 10;
	  if (abmult > 0)
	    cs = cs + 1;
	}
      while (!(abmult == 0));
      if (abs (W22->mult) == 1)
	cs = cs - 1;
      if ((((*qq) + cs < col) && (!tex)) || (tex && (tcolumns <= columns)))
	{
	  if (startx)
	    {
	      if (tex)
		(void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0);
	      if (((W22->mult < 0) && !tex))
		(void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
	      else
		cs = cs - 3;
	      if (((W22->mult < 0) && tex))
		(void) fprintf ((*fyle).fp, " - \\"), Putl ((*fyle), 0);
	      startx = false;
	    }
	  else
	    {
	      if (tex)
		if (tcolumns == 1)
		  (void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0);
		else
		  (void) fprintf ((*fyle).fp, "&$"), Putl ((*fyle), 0);
	      if (tex)
		if (W22->mult < 0)
		  (void) fprintf ((*fyle).fp, " - \\"), Putl ((*fyle), 0);
		else
		  (void) fprintf ((*fyle).fp, " + \\"), Putl ((*fyle), 0);
	      if (((!tex) && (W22->mult < 0)))
		(void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
	      if (((!tex) && (W22->mult > 0)))
		(void) fprintf ((*fyle).fp, " +"), Putl ((*fyle), 0);
	    }
	  if (abs (W22->mult) != 1)
	    (void) fprintf ((*fyle).fp, "%1d", abs (W22->mult)),
	      Putl ((*fyle), 0);
	  (*G142_texx) = false;
	  switch ((int) (currgrp.A[1 - 1].name))
	    {
	    case sn:
	      if (nreduce)
		{
		  a = '<';
		  b = '>';
		}
	      else
		{
		  a = '{';
		  b = '}';
		  (*G142_texx) = true;
		}
	      break;
	    case an:
	      if (nreduce)
		{
		  a = '<';
		  b = '>';
		}
	      else
		{
		  a = '{';
		  b = '}';
		  (*G142_texx) = true;
		}
	      break;
	    case sung:
	    case un:
	    case unm:
	    case sunm:
	    case unc:
	      a = '{';
	      b = '}';
	      (*G142_texx) = true;
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
	    (*G140_cays) = pe6;
	  else
	    (*G140_cays) = norm;
	  if ((tex && (*G142_texx)))
	    (void) fprintf ((*fyle).fp, "\\%c", a), Putl ((*fyle), 0);
	  else
	    Putchr (a, (*fyle));
	  if ((currgrp.A[1 - 1].name == unc))
	    Putchr ('(', (*fyle));
	  if ((W22->spin && (tex == false)))
	    {
	      Putchr (spinchar, (*fyle));
	      if (((currgrp.A[1 - 1].name != spnc)
		   && (currgrp.A[1 - 1].name != mp)
		   && (currgrp.A[1 - 1].name != sonc)))
		Putchr (';', (*fyle));
	    }
	  if ((W22->spin && tex && (currgrp.A[1 - 1].name != spnc)
	       && (currgrp.A[1 - 1].name != mp)
	       && (currgrp.A[1 - 1].name != sonc)))
	    (void) fprintf ((*fyle).fp, "s;"), Putl ((*fyle), 0);
	  if (((W22->spin && tex && (currgrp.A[1 - 1].name == spnc))
	       || ((currgrp.A[1 - 1].name == mp))))
	    (void) fprintf ((*fyle).fp, "s"), Putl ((*fyle), 0);
	  if ((W22->C6_double && (currgrp.A[1 - 1].name != spnc)
	       && (currgrp.A[1 - 1].name != mp)
	       && (currgrp.A[1 - 1].name != sonc)))
	    {
	      wrtfrme (&(*fyle), W22->conval);
	      Putchr (';', (*fyle));
	    }
	  if ((currgrp.A[1 - 1].name == spnc) || (currgrp.A[1 - 1].name == mp)
	      || (currgrp.A[1 - 1].name == sonc))
	    (void) fprintf ((*fyle).fp, "%1d(", W22->val.A[maxdim]),
	      Putl ((*fyle), 0);
	  if (((*G140_cays) == pe6))
	    {
	      (void) fprintf ((*fyle).fp, "%1d", W22->val.A[1]),
		Putl ((*fyle), 0);
	      Putchr (':', (*fyle));
	      if (W22->val.A[2] == 0)
		{
		  Putchr ('0', (*fyle));
		  i = 3;
		}
	      else
		i = 2;
	      if (!pow_note)
		{
		  while (W22->val.A[i] != 0)
		    {
		      (void) fprintf ((*fyle).fp, "%1d", W22->val.A[i]),
			Putl ((*fyle), 0);
		      i = i + 1;
		    }
		}
	      else
		{
		  if ((group == g2))
		    if ((racah == true))
		      {
			z = W22->val.A[1];
			W22->val.A[1] = z - W22->val.A[2];
		      }
		  m = i;
		  while ((W22->val.A[m] != 0))
		    {
		      l = W22->val.A[i];
		      i = i + 1;
		      k = 1;
		      while (W22->val.A[i] == l)
			{
			  i = i + 1;
			  k = k + 1;
			}
		      if (k > cutoff)
			{
			  if (W22->val.A[m] > 9)
			    Putchr (' ', (*fyle));
			  (void) fprintf ((*fyle).fp, "%1d^%1d",
					  W22->val.A[m], k), Putl ((*fyle),
								   0);
			}
		      else
			for (l = 1; l <= k; l++)
			  {
			    if (W22->val.A[m] > 9)
			      (void) fprintf ((*fyle).fp, "%1d",
					      W22->val.A[m]),
				Putl ((*fyle), 0);
			    else
			      (void) fprintf ((*fyle).fp, "%1d",
					      W22->val.A[m]),
				Putl ((*fyle), 0);
			  }
		      m = i;
		    }
		  if ((group == g2))
		    if ((racah == true))
		      W22->val.A[1] = z;
		}
	    }
	  else
	    if ((currgrp.A[1 - 1].name != spnc)
		&& (currgrp.A[1 - 1].name != mp)
		&& (currgrp.A[1 - 1].name != sonc))
	    wrtfrme (&(*fyle), W22->val);
	  if ((currgrp.A[1 - 1].name == spnc)
	      || (currgrp.A[1 - 1].name == mp))
	    {
	      wrtfrme (&(*fyle), W22->val);
	      (void) fprintf ((*fyle).fp, ")>"), Putl ((*fyle), 0);
	    }
	  if ((currgrp.A[1 - 1].name == spnc))
	    {
	      if ((tex && (W22->lab == '#')))
		(void) fprintf ((*fyle).fp, "\\#"), Putl ((*fyle), 0);
	      else if (W22->lab != ' ')
		Putchr ('#', (*fyle));
	    }
	  if ((currgrp.A[1 - 1].name == sonc))
	    {			//25/12/97
	      wrtfrme (&(*fyle), W22->val);
	      (void) fprintf ((*fyle).fp, ")]"), Putl ((*fyle), 0);
	    }
	  if ((currgrp.A[1 - 1].name == unc))
	    Putchr (')', (*fyle));
	  if ((currgrp.A[1 - 1].name != spnc) && (currgrp.A[1 - 1].name != mp)
	      && (currgrp.A[1 - 1].name != sonc))
	    {
	      if ((tex && (*G142_texx)))
		(void) fprintf ((*fyle).fp, "\\%c", b), Putl ((*fyle), 0);
	      else
		Putchr (b, (*fyle));
	      if (((W22->lab != ' ') && (tex == false)))
		Putchr (W22->lab, (*fyle));
	      if (((W22->lab != ' ') && tex))
		if (W22->lab == '#')
		  (void) fprintf ((*fyle).fp, "\\#"), Putl ((*fyle), 0);
		else
		  (void) fprintf ((*fyle).fp, "_%c", W22->lab), Putl ((*fyle),
								      0);
	    }

	  if (tex)
	    Putchr ('$', (*fyle));
	  (*qq) = (*qq) + cs + 1;
	  list = list->next;
	}
      else
	{
	  tcolumns = 0;
	  if ((vdu && (!tex)))
	    Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
	  else if (tex)
	    (void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
	  else
	    Putchr ('\n', (*fyle));
	  if (tex)
	    {
	      sline = sline + 1;
	      if (sline == tlines)
		{
		  sline = 1;
		  Putchr ('\n', (*fyle));
		  Putchr ('\n', (*fyle));
		}
	    }
	  if (vdu && more)
	    {
	      line = line + 1;
	      if (line == tlines)
		{
		  line = 1;
		  if ((key = more_Question ()) == 'q')
		    return;
		  else
		    {
		      RESTORE_CURSOR_POSITION ();
		      if (key != 'm')	// just one more line
			line = tlines - 1;
		    }
		}
	      else
		(void) fprintf (output.fp, "    "), Putl (output, 0);
	    }
	  else
	    (void) fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
	  (*qq) = 5;
	}
    }
  list = top;
  if ((vdu && (!tex)))
    Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
  else if (tex)
    (void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
  else
    Putchr ('\n', (*fyle));
  G140_cays = F141;
  G142_texx = F143;
}
*/

int
tsize (termptr list)
{
  register int R138;
  int i, j, k, l, m, n, tz, leng;

  {
    register termptr W25 = &(*list);

    leng = qlen (W25->val);
    if (!pow_note) {
      i = 1;
      if (W25->val.A[1] > 9)
	j = 2;
      else
	j = 0;
      while ((i <= leng)) {
	i = i + 1;
	if (W25->val.A[i] > 9)
	  j = j + 2;
      }
    } else {
      i = 1;
      m = 1;
      j = 1;
      while ((m <= leng)) {
	l = W25->val.A[i];
	i = i + 1;
	k = 1;
	while ((W25->val.A[i] == l) && (i < maxdim)) {
	  i = i + 1;
	  k = k + 1;
	}
	if (W25->val.A[m] > 9)
	  if (W25->val.A[m] > 99)
	    n = 3;
	  else
	    n = 2;
	else
	  n = 1;
	if (k > cutoff) {
	  if (n > 1)
	    j = n + j + 3;
	  else
	    j = n + j + 2;
	  if (k > 9)
	    j = j + 2;
	  else
	    j = j + 1;
	} else
	  j = j + k * n;
	m = i;
      }
      if (i != 1)
	i = 0;
    }
    if (abs (W25->mult) != 1)
      if (abs (W25->mult) > 9)
	tz = i + j + 7;
      else
	tz = i + j + 6;
    else
      tz = i + j + 5;
    if (((qfn == true) || (pfn == true)))
      tz = tz + 2;
    if ((mmono == true) || (homo == true) || (elem == true) || (forg == true)
	|| (psum == true))
      tz = tz + 2;
    if (abs (W25->mult) > 99)
      if (abs (W25->mult) > 999)
	tz = tz + 2;
      else
	tz = tz + 1;
  }
  R138 = tz;
  return R138;
}

/*
void
wrttlst (text * fyle, int *qq, int col, char a, char b, termptr list,
	 bool vdu)
{
  bool startx;
  int ts, line, tcolumns, sline;
  char key;

  line = 1;
  key = '?';
  tcolumns = 0;
  sline = 1;
  if ((*qq) == 1) {
    (void) fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
    (*qq) = 5;
  }
  if (((list == NULL) && (qtest == false))) {
    (void) fprintf ((*fyle).fp, "zero"), Putl ((*fyle), 0);
    (*qq) = (*qq) + 4;
  }
  startx = true;
  while ((list != NULL) && (key != esc)) {
    register termptr W26 = &(*list);

    ts = tsize (list);
    tcolumns = tcolumns + 1;
    if ((((*qq) + ts < col) && (!tex)) || (tex && (tcolumns <= columns))) {
      if (startx) {
	if (tex)
	  (void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0);
	if (W26->mult < 0)
	  if (tex)
	    (void) fprintf ((*fyle).fp, " - \\"), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
	else
	  ts = ts - 3;
	startx = false;
      } else {
	if (tex)
	  if (tcolumns == 1)
	    (void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, "&$"), Putl ((*fyle), 0);
	if (W26->mult < 0)
	  if (tex)
	    (void) fprintf ((*fyle).fp, " - \\"), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
	if (W26->mult > 0)
	  if (tex)
	    (void) fprintf ((*fyle).fp, " + \\"), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, " +"), Putl ((*fyle), 0);
      }
      if (abs (W26->mult) != 1)
	(void) fprintf ((*fyle).fp, "%1d", abs (W26->mult)),
	  Putl ((*fyle), 0);
      if ((qfn == true) && (pfn == false))
	(void) fprintf ((*fyle).fp, "Q_"), Putl ((*fyle), 0);
      else if ((pfn == true) && (qfn == true))
	(void) fprintf ((*fyle).fp, "P_"), Putl ((*fyle), 0);
      else if ((mmono == true))
	(void) fprintf ((*fyle).fp, "m_"), Putl ((*fyle), 0);
      else if ((forg == true))
	(void) fprintf ((*fyle).fp, "f_"), Putl ((*fyle), 0);
      else if ((homo == true))
	(void) fprintf ((*fyle).fp, "h_"), Putl ((*fyle), 0);
      else if ((psum == true))
	(void) fprintf ((*fyle).fp, "p_"), Putl ((*fyle), 0);
      else if ((elem == true))
	(void) fprintf ((*fyle).fp, "e_"), Putl ((*fyle), 0);
      if (tex
	  && !(qfn || nreduce || pfn || mmono || homo || elem || forg
	       || psum))
	(void) fprintf ((*fyle).fp, "\\%c", a), Putl ((*fyle), 0);
      else
	Putchr (a, (*fyle));
      wrtfrme (&(*fyle), W26->val);
      if (tex
	  && !(qfn || nreduce || pfn || mmono || homo || elem || forg
	       || psum))
	(void) fprintf ((*fyle).fp, "\\%c$", b), Putl ((*fyle), 0);
      else
	if (tex
	    && (qfn || nreduce || pfn || mmono || homo || elem || forg
		|| psum))
	(void) fprintf ((*fyle).fp, "%c$", b), Putl ((*fyle), 0);
      else
	Putchr (b, (*fyle));
      if ((W26->slab == '#'))
	(void) fprintf ((*fyle).fp, "\\#"), Putl ((*fyle), 0);
      if (tex
	  && !(qfn || mmono || homo || elem || forg || pfn || nreduce
	       || psum))
	(*qq) = (*qq) + ts + 3;
      else
	(*qq) = (*qq) + ts + 1;
      list = W26->next;
    } else {
      tcolumns = 0;
      if ((vdu && (!tex)))
	fprintf (fyle->fp, " \n");
      //  Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
      else if ((vdu && tex))
	(void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
      else if (tex)
	(void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
      else
	Putchr ('\n', (*fyle));
      if (tex) {
	sline = sline + 1;
	if (sline == tlines) {
	  sline = 1;
	  Putchr ('\n', (*fyle));
	  Putchr ('\n', (*fyle));
	}
      }
      if (vdu && more) {
	line = line + 1;
	if (line == tlines) {
	  line = 1;
	  if ((key = more_Question ()) == 'q')
	    return;
	  else {
	    RESTORE_CURSOR_POSITION ();
	    if (key != 'm')	// just one more line
	      line = tlines - 1;
	  }
	} else
	  (void) fprintf (output.fp, "    "), Putl (output, 0);
      } else
	(void) fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
      (*qq) = 5;
    }
  }
  if ((vdu && (!tex)))
    Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
  else if ((vdu && tex))
    (void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
  else if (tex)
    (void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
  else
    Putchr ('\n', (*fyle));
}
*/

int
psize (ocharptr list)
{
  register int R139;
  int z, i, j, k, l, m;

  {
    register ocharptr W27 = &(*list);

    z = W27->val.A[1];
    if ((group == g2))
      if ((racah == true))
	W27->val.A[1] = z - W27->val.A[2];
    if (!pow_note) {
      i = 1;
      if (W27->val.A[i] == 0)
	j = 1;
      else
	j = 0;
      while (W27->val.A[i] != 0) {
	k = abs (W27->val.A[i]);
	if (W27->val.A[i] < 0)
	  j = j + 1;
	if (k > 9)
	  if (k > 99)
	    j = j + 4;
	  else
	    j = j + 3;
	else
	  j = j + 1;
	i = i + 1;
      }
    } else {
      i = 1;
      if (W27->val.A[i] == 0)
	j = 1;
      else
	j = 0;
      while (W27->val.A[i] != 0) {
	k = abs (W27->val.A[i]);
	l = 1;
	m = W27->val.A[i];
	i = i + 1;
	while ((W27->val.A[i] == m) && (i < maxdim - 1)) {
	  l = l + 1;
	  i = i + 1;
	}
	if (m < 0)
	  m = 1;
	else
	  m = 0;
	if (k > 9)
	  if (k > 99)
	    m = m + 4;
	  else
	    m = m + 3;
	else
	  m = m + 1;
	if (l > cutoff) {
	  j = j + m + 2;
	  if (l > 9)
	    j = j + 2;
	  else
	    j = j + 1;
	} else
	  j = j + l * m;
      }
    }
    if ((*G144_cays) == pe6) {
      j = j + 1;
      if (W27->val.A[2] == 0)
	j = j + 1;
    }
    if (((W27->spin && tex) == false))
      j = j + 2;
    if ((W27->spin && tex))
      j = j + 7;
    if (W27->lab != ' ')
      j = j + 1;
    if (((W27->lab != ' ') && tex))
      j = j + 1;
    if (W27->C6_double) {
      i = 1;
      if (W27->conval.A[i] == 0)
	j = j + 1;
      while (W27->conval.A[i] != 0) {
	k = abs (W27->conval.A[i]);
	if (W27->conval.A[i] < 0)
	  j = j + 1;
	if (k > 9)
	  if (k > 99)
	    j = j + 4;
	  else
	    j = j + 3;
	else
	  j = j + 1;
	i = i + 1;
      }
    }
    W27->val.A[1] = z;
    R139 = j + 2;
  }
  return R139;
}

/*
void
writeprod (text * fyle, prodtype pr, int col, bool vdu)
{
  bool startx, texx;
  int z = 0,			// modified by FB, was uninitialized
    i, cs, qq, abmult, line, k, m, sline, tcolumns;
  register int l;
  register int j;
  caystype cays;
  prodtype top;
  char a = '{', b = '}',	// modified by FB, was uninitialized
    spinchar, key;
  caystype *F145;

  F145 = G144_cays;
  G144_cays = &cays;
  line = 1;
  key = '?';
  spinchar = 's';
  sline = 1;
  tcolumns = 0;
  qq = fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
  top = pr;
  if (pr == NULL) {
    (void) fprintf ((*fyle).fp, "zero"), Putl ((*fyle), 0);
    qq = qq + 4;
  }
  startx = true;
  if (tex) {
    fprintf ((*fyle).fp, "\\begin{tabular}{*{%d}{l}}\n$", columns),
      Putl ((*fyle), 1);
  }
  while ((pr != NULL) && (key != esc)) {
    register prodtype W28 = &(*pr);

    tcolumns = tcolumns + 1;
    cs = 0;
    for (j = 1; j <= nprod; j++) {
      if (currgrp.A[j - 1].name == e6)
	(*G144_cays) = pe6;
      else
	(*G144_cays) = norm;
      cs = cs + psize (W28->prods.A[j - 1]);
    }
    abmult = abs (W28->mult);
    cs = cs + 4;
    do {
      abmult = abmult / 10;
      if (abmult > 0)
	cs = cs + 1;
    }
    while (!(abmult == 0));
    if (abs (W28->mult) == 1)
      cs = cs - 1;
    if (((qq + cs < col) && (!tex)) || (tex && (tcolumns <= columns))) {
      if (startx) {
	if (tex)
	   (void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0); 
	if (((W28->mult < 0) && !tex))
	  qq +=fprintf ((*fyle).fp, " -"), Putl ((*fyle), 0);
	else
	  cs = cs - 3;
	if (((W28->mult < 0) && tex))
	  (void) fprintf ((*fyle).fp, " -\\ "), Putl ((*fyle), 0);
	startx = false;
      } else {
	if (tex)
	  if (tcolumns == 1)
	    (void) fprintf ((*fyle).fp, "\\+$"), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, "&$"), Putl ((*fyle), 0);
	if (tex)
	  if (W28->mult < 0)
	    (void) fprintf ((*fyle).fp, " -\\ "), Putl ((*fyle), 0);
	  else
	    (void) fprintf ((*fyle).fp, " +\\ "), Putl ((*fyle), 0);
	if (((!tex) && (W28->mult < 0)))
	  qq +=fprintf ((*fyle).fp, " - "), Putl ((*fyle), 0);
	if (((!tex) && (W28->mult > 0)))
	  qq += fprintf ((*fyle).fp, " + "), Putl ((*fyle), 0);
      }
      if (abs (W28->mult) != 1)
	fprintf ((*fyle).fp, "%1d", abs (W28->mult)), Putl ((*fyle), 0);

      for (j = 1; j <= nprod; j++) {
	texx = false;
	switch ((int) (currgrp.A[j - 1].name)) {
	case sn:
	  if (nreduce) {
	    a = '<';
	    b = '>';
	  } else {
	    a = '{';
	    b = '}';
	    texx = true;
	  }
	  break;
	case an:
	  if (nreduce) {
	    a = '<';
	    b = '>';
	  } else {
	    a = '{';
	    b = '}';
	    texx = true;
	  }
	  break;
	case sung:
	case un:
	case unm:
	case sunm:
	case unc:
	  a = '{';
	  b = '}';
	  texx = true;
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
	  //spinc = W28->prods.A[j - 1]->spin;
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
	  (*G144_cays) = pe6;
	else
	  (*G144_cays) = norm;
	if ((tex && texx))
	  (void) fprintf ((*fyle).fp, "\\%c", a), Putl ((*fyle), 0);
	else
	  Putchr (a, (*fyle));
	if (currgrp.A[j - 1].name == unc)
	  Putchr ('(', (*fyle));
	if ((W28->prods.A[j - 1]->spin && (tex == false))) {
	  Putchr (spinchar, (*fyle));
	  if (((currgrp.A[j - 1].name != spnc)
	       && (currgrp.A[1 - 1].name != mp)
	       && (currgrp.A[1 - 1].name != sonc)))
	    Putchr (';', (*fyle));
	}
	if ((W28->prods.A[j - 1]->spin && tex
	     && (currgrp.A[j - 1].name != spnc)
	     && (currgrp.A[1 - 1].name != mp)
	     && (currgrp.A[1 - 1].name != sonc)))
	  (void) fprintf ((*fyle).fp, "s;"), Putl ((*fyle), 0);
	if ((W28->prods.A[j - 1]->C6_double && (currgrp.A[j - 1].name != spnc)
	     && (currgrp.A[1 - 1].name != mp)
	     && (currgrp.A[1 - 1].name != sonc))) {
	  wrtfrme (&(*fyle), W28->prods.A[j - 1]->conval);
	  Putchr (';', (*fyle));
	}
	if (((currgrp.A[j - 1].name == spnc)
	     || (currgrp.A[1 - 1].name == mp)
	     || (currgrp.A[1 - 1].name == sonc))) {
	  if (((W28->prods.A[j - 1]->spin) && tex))
	    Putchr ('s', (*fyle));
	  (void) fprintf ((*fyle).fp, "%1d;(",	// modified by FB was %1d;(
			  W28->prods.A[j - 1]->val.A[maxdim]),
	    Putl ((*fyle), 0);
	}
	if (((*G144_cays) == pe6)) {
	  (void) fprintf ((*fyle).fp, "%1d",
			  W28->prods.A[j - 1]->val.A[1]), Putl ((*fyle), 0);
	  Putchr (':', (*fyle));
	  if (W28->prods.A[j - 1]->val.A[2] == 0) {
	    Putchr ('0', (*fyle));
	    i = 3;
	  } else
	    i = 2;
	  if (!pow_note) {
	    while (W28->prods.A[j - 1]->val.A[i] != 0) {
	      (void) fprintf ((*fyle).fp, "%1d",
			      W28->prods.A[j - 1]->val.A[i]),
		Putl ((*fyle), 0);
	      i = i + 1;
	    }
	  } else {
	    register ocharptr W33 = &(*W28->prods.A[j - 1]);

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
		if (W33->val.A[m] > 9)
		  Putchr (' ', (*fyle));
		(void) fprintf ((*fyle).fp, "%1d^%1d",
				W33->val.A[m], k), Putl ((*fyle), 0);
	      } else
		for (l = 1; l <= k; l++) {
		  if (W33->val.A[m] > 9)
		    (void) fprintf ((*fyle).fp, "%1d",
				    W33->val.A[m]), Putl ((*fyle), 0);
		  else
		    (void) fprintf ((*fyle).fp, "%1d",
				    W33->val.A[m]), Putl ((*fyle), 0);
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
	  wrtfrme (&(*fyle), W28->prods.A[j - 1]->val);
	if ((currgrp.A[j - 1].name == spnc)
	    || (currgrp.A[j - 1].name == mp)) {
	  wrtfrme (&(*fyle), W28->prods.A[j - 1]->val);
	  (void) fprintf ((*fyle).fp, ")>"), Putl ((*fyle), 0);
	}
	if ((currgrp.A[1 - 1].name == sonc)) {
	  wrtfrme (&(*fyle), W28->prods.A[j - 1]->val);
	  (void) fprintf ((*fyle).fp, ")]"), Putl ((*fyle), 0);
	}
	if (currgrp.A[j - 1].name == unc)
	  Putchr (')', (*fyle));
	if (((currgrp.A[j - 1].name != spnc)
	     && (currgrp.A[j - 1].name != mp)
	     && (currgrp.A[1 - 1].name != sonc))) {
	  if ((tex && texx))
	    (void) fprintf ((*fyle).fp, "\\%c", b), Putl ((*fyle), 0);
	  else
	    Putchr (b, (*fyle));
	  if ((W28->prods.A[j - 1]->lab != ' ') && !tex)
	    Putchr (W28->prods.A[j - 1]->lab, (*fyle));
	  if ((W28->prods.A[j - 1]->lab != ' ') && tex)
	    if ((W28->prods.A[j - 1]->lab == '#'))
	      (void) fprintf ((*fyle).fp, "\\#"), Putl ((*fyle), 0);
	    else
	      (void) fprintf ((*fyle).fp, "_%c",
			      W28->prods.A[j - 1]->lab), Putl ((*fyle), 0);
	}
      }
      if (tex)
	Putchr ('$', (*fyle));
      qq = qq + cs + 1;
      pr = pr->next;
    } else {
      tcolumns = 0;
      if ((vdu && (!tex)))
	Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
      else if (tex)
	(void) fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
      else
	Putchr ('\n', (*fyle));
      if (tex) {
	sline = sline + 1;
	if (sline == tlines) {
	  sline = 1;
	  Putchr ('\n', (*fyle));
	  Putchr ('\n', (*fyle));
	}
      }
      if (vdu && more) {
	line = line + 1;
	if (line == tlines) {
	  line = 1;
	  if ((key = more_Question ()) == 'q')
	    return;
	  else {
	    RESTORE_CURSOR_POSITION ();
	    if (key != 'm')	// just one more line
	      line = tlines - 1;
	  }
	} else
	  (void) fprintf (output.fp, "    "), Putl (output, 0);
      } else
	(void) fprintf ((*fyle).fp, "    "), Putl ((*fyle), 0);
      qq = 5;
    }
  }
  pr = top;
  if ((vdu && (!tex)))
    Putchr (cr, (*fyle)), Putchr ('\n', (*fyle));
  else if (tex) {
    fprintf ((*fyle).fp, "\\cr\n"), Putl ((*fyle), 1);
    fprintf ((*fyle).fp, "\\end{tabular}\n"), Putl ((*fyle), 1);
  } else
    Putchr ('\n', (*fyle));
  G144_cays = F145;
}*/

/* read a user-defined function and store it in fu */
void
readfnset (text * fyle, fnptrs * fu)
{
  fnptrs help = NULL, f;
  char word[bcol];
  string0 buffx;
  int k;

  if (*fu != NULL)
    fdisp (fu);
  do {
    fprintf (output.fp, "=-"), Putl (output, 1);
    readacard (&(*fyle), &buffx, &k);
    if (debug_schur) {
      memcpy (word, buffx.A, bcol);
      word[bcol - 1] = '\0';
      fprintf (output.fp, "%s", word);
    }
    fprintf (output.fp, "\n");
    readword (&buffx, &k, word);
    if (!feof (fyle->fp) && !interp ("stop", word, 4)
	&& (buffx.A[0] != ' ')) {
      fnu (&f);
      f->gbuff = buffx;
      if (*fu == NULL) {
	(*fu) = f;
	help = f;
      } else {
	if (help == NULL) {	// FB coherence test
	  aaargh ("error in readfnset", true);
	  exit (1);
	}
	help->next = f;
	help = f;
      }
    }
  }
  while (!feof (fyle->fp) && !interp ("stop", word, 4));
}
