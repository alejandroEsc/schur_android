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

/** \file FrontPage.c
 * display screen of schur
 */

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

#include "standard.h"
#include "dim.h"
// For dpm
#include "type.h"
#include "var.h"
// For statuscom
#include "s8.h"
#include "FrontPage.h"

#define mydisplay(m)   fprintf (out,"%s%s%s\n",colorsBegin,m,colorsEnd);

void
FrontPage (FILE * out)
{

  mydisplay
    ("+-----------------------------------------------------------------------------+");

  fprintf (out, "%s|                            SCHUR %-43s|%s\n",
	   colorsBegin, VERSION, colorsEnd);
  mydisplay
    ("| Copyright (C) 1996 Brian G. Wybourne,                                       |");
  mydisplay
    ("|               2006 Franck BUTELLE, Steven M. Christensen,                   |");
  mydisplay
    ("|                    Ronald C. KING & Frederic TOUMAZET                       |");
  mydisplay
    ("| SCHUR comes with ABSOLUTELY NO WARRANTY. This is free software, and you are |");
  mydisplay
    ("| welcome to redistribute it under certain conditions; type ?LICENSE for      |");
  mydisplay
    ("| details.  - If you wish to exit, enter END or QUIT                          |");
  mydisplay
    ("|           - If you wish to obtain help: HELP or ?cmd or APROPOS search      |");
  mydisplay
    ("| Please report bugs to                                                       |");
  mydisplay
    ("|                                Franck.Butelle@lipn.fr & toumazet@univ-mlv.fr|");
  mydisplay
    ("+-----------------------------------------------------------------------------+");

  fprintf (out, "\n");

  fprintf (out, "                              *** Status ***\n");
  statuscom (dpm);
  fprintf (out, "                              *------------*\n\n");
}
