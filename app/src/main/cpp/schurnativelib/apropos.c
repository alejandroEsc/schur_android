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

/** \file apropos.c
 */
# include "apropos.h"

/** Search in the help files if a string match
 *  @param s : the string to find */
void
apropos (char *s)
{
  int i;
  char *ptr;
  char search[MAXSTRING], cmd[MAXSTRING], fil[MAXSTRING];
  FILE *f;

  i = strlen (s) - 1;
  while (i > 0 && s[i] == ' ')
    i--;
  s[i + 1] = '\0';
  strcpyuppercase (search, s);

  if (debug_schur)
    fprintf (stderr, "looking for %s in %d commands\n", search, nbCommands);

  for (i = 0; i < nbCommands; i++)
    {
      strcpyuppercase (cmd, CommandTab[i]);
      if (debug_schur)
	fprintf (stderr, "<%s> ", cmd);

      if ((ptr = strstr (cmd, search)) != NULL)
	{
	  printf ("%s:", CommandTab[i]);
	  sprintf (fil, "%s/help/%s", dataPath, CommandTab[i]);
	  f = fopen (fil, "r");
	  if (f != NULL)
	    {
	      while (fgets (fil, MAXSTRING, f) != NULL)
		     if (strncmp (stripwhite(fil), "Description:",
				 strlen ("Description:")) == 0)
			     break;
		;
	      if (fil != NULL)
		{
		  char *ptr = index (fil, '-');
		  if (ptr == NULL)
		    printf ("%s", fil);
		  else
		    {
		      ptr++;
		      if (ptr[strlen (ptr) - 1] == '\n')
			ptr[strlen (ptr) - 1] = '\0';
		      /* display only the first line of the help file */
		      ptr[tcol - strlen (CommandTab[i]) - 4] = '\0';
		      printf ("%s...\n", ptr);
		    }
		}
	    }
	}
    }
  printf ("\n");
}
