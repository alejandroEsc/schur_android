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

/** \file schur.c
 */

/*
**	Definitions for i/o
*/
# include <stdio.h>
# include <stdbool.h>
# include <getopt.h>

# include "standard.h"
# include "define.h"
# include "dim.h"
# include "type.h"
//# include "var.h"
# include "utils.h"
# include "init.h"
# include "FrontPage.h"
# include "../config.h"

//for dpmode(), repmode(), sfnmode()
# include "s0.h"

//for interp()
# include "r.h"

//for initialise_readline()
# include "FillCommandTab.h"

void DisplayHelp(void);
void DisplayVersion(void);

text input, output;

void 
DisplayHelp(void) 
{
  printf("\n SCHUR \n\n");
  printf("An Interactive Program For \n");
  printf("Calculating Properties \n");
  printf("Of Lie Groups and \n");
  printf("Symmmetric\n");
  printf("Functions\n");
  printf("\nGPL version %s\n", VERSION);
  printf("\nSCHUR comes with ABSOLUTELY NO WARRANTY. This is free software, and you are\n");
  printf("welcome to redistribute it under certain conditions; type ? LICENSE for details.\n");

  printf("\nthe HELP command give the command list (over 200 commands)\n");
  printf("\n ? COMMAND         give help on COMMAND\n");
  printf("\n*** Command-line options:\n");

  printf("-d : start in debug mode\n");
  printf("-h : display this help\n");
  printf("-q : a bit more quiet, may be useful for shell scripts\n");
  printf("-f filename : read filename as a script for SCHUR\n");
  printf("-v : display version and stop\n");

  printf("\nPlease report bugs to\n");
  printf("Franck.Butelle@lipn.univ-paris13.fr and Frederic.Toumazet@lipn.univ-paris13.fr (you may use -d option when running schur to help us to find the error)\n");
}

void 
DisplayVersion(void) 
{
  printf("SCHUR GPL version %s\n", VERSION);
}

/*
**	Start of program code
*/
int
main (int argc, char *argv[])
{  
  char c;
  int option_index = 0;

  const char* script_filename = NULL;

  const char* const short_options = "dhqvf:"; // f option need an argument
  static struct option long_options[] =
  {
     {"debug", 0, NULL, 'd'},
     {"help", 0, NULL, 'h'},
     {"quiet",0,NULL, 'q'},
     {"file",1,NULL, 'f'},
     {"version", 0, NULL, 'v'},
     {0, 0, 0, 0}
  };


  input.fp = stdin;
  input.eoln = 0;
  input.eof = 0;
  input.out = 0;
  input.init = 0;
  input.buf = '\0';
  output.fp = stdout;
  output.eoln = 0;
  output.eof = 0;
  output.out = 0;
  output.init = 0;
  output.buf = '\0';

# ifdef STDINIT
  (void) (Getx (input));
# endif

  initialise ();
  debug_schur=false; quiet=false; 
  strcpy (instr, "dpm" ); // default initial mode
   
   while ((c = getopt_long (argc, argv, short_options, long_options, &option_index))!= -1)
   {
        switch(c) {
	  case 'd': debug_schur=true; break;
          case 'h': DisplayHelp(); return EXIT_SUCCESS; 
          case 'q': quiet=true; break;
          case 'f': 
		    script_filename = optarg; 
		    input.fp = fopen(script_filename, "r");
		    if (input.fp == NULL)
			    fprintf(stderr,"File not found\n");
		    break;
          case 'v': DisplayVersion(); return EXIT_SUCCESS;
	  default: fprintf(stderr, "Unknown parameter %d",c); DisplayHelp(); return EXIT_FAILURE;
	}
   }

   if (optind < argc)
   {
     fprintf (stderr, "The following arguments are unknown :");
     while (optind < argc)
       fprintf (stderr, "%s ", argv[optind++]);
     fprintf (stderr, "\n");
   }

  initialise_readline();
  nbCommands=fillCommandTab();	// fill the list of commands for command completion

  if (! quiet)
  	FrontPage(output.fp); // licence, status etc.

  do
    {
      if (interp ("dpmode", instr, 3) || interp ("exitmode", instr, 4))
	dpmode ();
      else if (interp ("repmode", instr, 3))
	repmode ();
      else if (interp ("sfnmode", instr, 3))
	sfnmode ();
    }
  while (!(interp ("end", instr, 3) || interp ("quit", instr, 4)));
  if (logging)
    {
      fprintf (output.fp, "%s has been closed.\n",logname.A);
    }
  free_all();
  return (0);
}
