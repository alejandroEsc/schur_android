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

/** \file FillCommandTab.c
 */

# include "standard.h"

#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <strings.h>
#include <string.h>

#include "dim.h"
#include "type.h"
#include "var.h"
#include "FillCommandTab.h"

#include "../config.h"
#ifdef HAVE_LIBREADLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

static int
cmpstringp (const void *p1, const void *p2)
{
  /* The actual arguments to this function are "pointers to
   * pointers to char", but strcasecmp() arguments are "pointers
   * to char", hence the following cast plus dereference */

  //return strcasecmp(* (char * const *) p1, * (char * const *) p2);
  return strcmp (*(char *const *) p1, *(char *const *) p2);
}

int
fillCommandTab (void)
{
  //FILE *f;
  struct dirent *d;
  DIR *dirstream;
  char helpdir[MAXSTRING];
  //char s[MAXSTRING], *p;

  sprintf (helpdir, "%s/help", dataPath);

  dirstream = opendir (helpdir);

  int i = 0;

  if (dirstream == NULL)
    {
      perror ("Error opening help files directory\n");
      exit (1);
    }

  while ( (d=readdir(dirstream)) != NULL) {
       if (d->d_name[0] >= 'A' && d->d_name[0]<='Z'
          && d->d_name[strlen(d->d_name)-1] !='~') {
             if (i >= MAXNBCOMMANDS) {
               fprintf(stderr, "ERROR: Table of commands is full !\n");
               exit(1);
             }
             CommandTab[i++] = strdup(d->d_name);
       }
 }

 qsort((void *)CommandTab, i, sizeof (char *), cmpstringp);

 closedir(dirstream);

 return(i);
}

/* Generator function for command completion.  STATE lets us
 * know whether to start from scratch; without any state
 * (i.e. STATE == 0), then we start at the top of the list. */
char *
command_generator (const char *text, int state)
{
  static int list_index, len;
  char *name, *t;

  t = (char *) text;
  if (text[0] == '?')
    t = (char *) text + 1;

  /* If this is a new word to complete, initialize now.  This
   * includes saving the length of TEXT for efficiency, and
   * initializing the index variable to 0. */
  if (!state)
    {
      list_index = 0;
      len = strlen (t);
    }

  /* Return the next name which partially matches from the
   * command list. */
  while ((name = CommandTab[list_index]))
    {
      list_index++;

      if (strncasecmp (name, t, len) == 0)
	if (text[0] == '?')
	  {
	    char *c = (char *) malloc (strlen (name) + 2);
	    sprintf (c, "?%s", name);
	    return (c);
	  }
	else
	  return (strdup (name));
    }

  /* If no names matched, then return NULL. */
  return ((char *) NULL);
}

/* Attempt to complete on the contents of TEXT.  START and END
 * bound the region of rl_line_buffer that contains the word to
 * complete.  TEXT is the word to complete.  We can use the entire
 * contents of rl_line_buffer in case we want to do some simple
 * parsing.  Return the array of matches, or NULL if there aren't any. */
char **
Command_completion (const char *text, int start, int end)
{
  UNUSED(start),UNUSED(end);
  char **matches;

  matches = (char **) NULL;
#ifdef HAVE_LIBREADLINE
  /*if (text[0]=='?') //forget the leading ? for completion
     {
     matches = rl_completion_matches (text+1, command_generator);
     }
     else */
  matches = rl_completion_matches (text, command_generator);
#endif
  return (matches);
}

void
initialise_readline (void)
{
#ifdef HAVE_LIBREADLINE
  /* Allow conditional parsing of the ~/.inputrc file. */
  rl_readline_name = "fillCommand";

  /* Tell the completer that we want a crack first. */
  rl_attempted_completion_function = Command_completion;
#endif
}
