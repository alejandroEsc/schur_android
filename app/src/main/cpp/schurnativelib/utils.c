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

/* This file contains Utilities */

# include <stdio.h>
# include <strings.h>
# include <string.h>
# include <stdarg.h>
# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
# include "var.h"
# include "mymalloc.h"
# include "sets_mgmt.h"
# include "../config.h"
# include "utils.h"

char
locase (char ch)
{
  if ((ch >= 'A') && (ch <= 'Z'))
    return ((unsigned) (ch) + 'a' - 'A');
  else
    return ch;
}

#ifndef __GNUC__
char *
strdup(char *old)
{
  size_t len = strlen(old) + 1; 
  char *new = malloc(len); 
  return new ? memcpy(new, old, len) : NULL;
}
#endif

void
strcpylowercase (char *dest, char *src)
{
  char *ptrd, *ptrs;

  for (ptrs = src, ptrd = dest; *ptrs != '\0'; ptrs++, ptrd++)
    if ((*ptrs >= 'A') && (*ptrs <= 'Z'))
      *ptrd = *ptrs - 'A' + 'a';
    else
      *ptrd = *ptrs;
  *ptrd = '\0';
}

void
strcpyuppercase (char *dest, char *src)
{
  char *ptrd, *ptrs;

  for (ptrs = src, ptrd = dest; *ptrs != '\0'; ptrs++, ptrd++)
    if ((*ptrs >= 'a') && (*ptrs <= 'z'))
      *ptrd = *ptrs + 'A' - 'a';
    else
      *ptrd = *ptrs;
  *ptrd = '\0';
}

bool
interp (char *instrx, char *cand, int n)
{
  return (strncasecmp (instrx, cand, n) == 0);
}

void
inform (char *message, char comm)
{
  char *ptr = strdup (message), *ptr2;

  if ((ptr2 = rindex (ptr, ';')) != NULL)
    *ptr2 = '\0';
  if (! quiet) {
  	print (ptr);
  	fflush (output.fp);
  	if (comm == cr)
    		print ("\n");
  }
  echo = true;
  free (ptr);
}

void
aaargh (char *msg, bool intern)
{
  fprintf (output.fp, "%cerror: %.24s\n", bel, msg), Putl (output, 1);
  if (intern)
    fprintf (output.fp, "( internal! )\n"), Putl (output, 1);
  fprintf (output.fp, "Please contact %s.\n", PACKAGE_BUGREPORT),
    Putl (output, 1);
  exit (1);
}

int
charval (char k)
{

  if (Member ((unsigned) (k), numbers.S))
    return (k - '0');
  else
    {
      fprintf (stderr, "char=%c(%x)", k, k);
      aaargh ("bad charval", true);
      return 0;
    }
}


char
boolchar (bool b)
{
  return (b ? 'T' : 'F');
}

int
sgn (int s)
{
  sign = (unsigned) (s > 0) - (unsigned) (s < 0);
  return sign;
}

void
progress (void)
{
  static char codes[] = "|/-\\";
  static short beat = 0;
  static short count = 0;

  if (sb_prog)
    {
      if (count == 0)
	{
	  printf ("%c\x1B[1D", codes[beat]);
	  fflush (stdout);
	  //putchar(codes[beat]);fflush(stdout);
	  beat = (beat + 1) % strlen (codes);
	}
      if (++count > 30)
	count = 0;		// to slow down !
    }
}

int
minusoneto (int n)
{
  if (n & 1)			/*if odd */
    return -1;
  else
    return 1;
}

void
warn (char *warning, char comm) {
  print (">>>>>");
  inform (warning, comm);
}

void
error (int errn, int posn) {
  char *errmsg;
  register int i;

  if (posn < bcol)
    posn += strlen (displayModes[mode]);
  errmsg = mymalloc (MAXSTRING);
  if (!erred) {
      if (errn < WARNING_LIMIT)
	erred = true;
      ern = errn;
      errpos = posn;
      if (logging)
	(void) fprintf (logfile.fp, "  "), Putl (logfile, 0);
      switch ((int) (errn))
	{
	case MISTAKE:
	  (void) strcpy (errmsg, "mistake.");
	  break;
	case COMMAND_NOT_FOUND:
	  (void) strcpy (errmsg, "command not found");
	  break;
	case MISSING_COMMA:
	  (void) strcpy (errmsg, "comma?");
	  break;
	case MISSING_PARAMETER:
	  (void) strcpy (errmsg, "parameter?");
	  break;
	case MISSING_INTEGER:
	  (void) strcpy (errmsg, "int?");
	  break;
	case TOO_MANY_PARTS:
	  (void) strcpy (errmsg, "too many parts.");
	  break;
	case WRONG_SERIES:
	  (void) strcpy (errmsg, "series?");
	  break;
	case WRONG_DIGIT:
	  (void) strcpy (errmsg, "digit?");
	  break;
	case NOT_IMPLEMENTED:
	  (void) strcpy (errmsg, "not implemented.");
	  break;
	case CANNOT_BE_ZERO:
	  (void) strcpy (errmsg, "cannot be 0.");
	  break;
	case INVALID_SERIES:
	  (void) strcpy (errmsg, "invalid series.");
	  break;
	case FILENAME_ERR:
	  (void) strcpy (errmsg, "file name?");
	  break;
	case NO_LOG_FILE:
	  (void) strcpy (errmsg, "no log file");
	  break;
	case BAD_CONJUGATE:
	  (void) sprintf(errmsg, "part too big >%d.",maxdim);
	  break;
	case BAD_OUTER_SKEW:
	  (void) strcpy (errmsg, "bad outer/skew.");
	  break;
	case NO_SUCH_FILE:
	  (void) strcpy (errmsg, "no such file.");
	  break;
	case BAD_FUNCTION:
	  (void) strcpy (errmsg, "bad function.");
	  break;
	case NOT_IN_FUNCTION:
	  (void) strcpy (errmsg, "not in function.");
	  break;
	case MISSING_QUOTE:
	  (void) strcpy (errmsg, "quote?");
	  break;
	case GROUP_NOT_SET:
	  (void) strcpy (errmsg, "group not set");
	  break;
	case ZERO_DIVISOR:
	  (void) strcpy (errmsg, "zero divisor?");
	  break;
	case BAD_IRREP:
	  (void) strcpy (errmsg, "bad irrep?");
	  break;
	case PART_TOO_BIG:
	  (void) sprintf (errmsg, "Part too big (max=%d).", maxdim-1);
	  break;
	case LINE_TOO_LONG:
	  (void) sprintf (errmsg, "line too long (limit is %d char.)", bcol);
	  break;
	case PARENTHESIS_NOT_NEEDED:
	  (void) strcpy (errmsg, "( not needed )");
	  break;
	case ARE_YOU_SURE:
	  (void) strcpy (errmsg, "( are you sure?)");
	  break;
	default:
	  aaargh ("bad error #", true);
	}
      // try to locate the error with a ^ sign under it
      if (posn < MAXSTRING) {
	  for (i = 1; i <= posn - 1; i++)
	    print (" ");
	  print ("^\n");
	  inform (errmsg, cr);
      } else
	inform (errmsg, cr);
    }
  free (errmsg);
//  print(" Hit return to continue... but this program may crash...");
//  c=getchar();
}

/* print to stdout and to logfile if needed */
void
print (char *format, ...)
{
  va_list ap;

  va_start (ap, format);
  vprintf (format, ap);
  va_end (ap);
  if (logging) {
    va_start (ap, format); // due to pb on 64 bits mach., the va list have to be "restarted". 
    vfprintf (logfile.fp, format, ap);
    va_end (ap);
  }
}

void
Caseerror (int n)
{
  if (mode == repm || mode == dpm)
    {
      print (">>>>>ERROR: group not set or another mistake\n");
    }
  else
    {
      fprintf (stderr, "Missing case limb: line %d\n", n);
      //exit (1);
    }
}


/**  drop spaces at the begining of a string */
char *
stripwhite (char *str)
{
  char *start;

  for (start = str;
       *start != '\0' && (*start == ' ' || *start == '\n' || *start == '\r'
			  || *start == '\t'); start++);

  return start;
}

char *
findLastSpace (char *str)
{
  char *start;
  start = str;
  str += bcol - 1;
  while (*str == ' ' && str != start)
    str--;
  if (str != start + bcol - 1)
    str++;
  return (str);
}


/** return length of a partition */
unsigned short
len (frame * f)
{
  register int l;
  l = 0;
  while (l < maxdim && f->A[l + 1] != 0)
    l = l + 1;
  f->length = l;		// refresh frame length if needed
  return l;
}

/** sum elements of a partition but returns 0 if invalid,
 * used by hiveslrcoeff */
unsigned
sumOfPartition (frame * part, unsigned n)
{
  unsigned i, res = 0;

  for (i = 1; i <= n; i++)
    {
      res += part->A[i];
      if (i < n && part->A[i + 1] > part->A[i])	// not a valid partition
	return (0);
    }
  return (res);
}

/** returns max element of a partition */
unsigned
maxOfPartition (frame * part)
{
  int i = 0, m = 0;

  if (part->length > maxdim)
    part->length = len (part);
  for (i = 1; i <= part->length; i++)
    if (part->A[i] > m)
      m = part->A[i];
  return (m);
}

void
free_all (void)
{
  free (instr);
}
