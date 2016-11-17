#include "ReadWrite.h"

void 
Resetx2(text *f, char *name, int l) 
{
	f->eof=0;
	f->out=0;
	/*if (f->init !=1) {*/
		f->fp = Fopen(name, l, Rmode);
		f->init = 1;
/*	}
	else
	    fseek(f->fp, 0L, SEEK_SET);*/

	if ((f->buf = fgetc(f->fp)) == '\n')
	{
		f->buf = ' ';
		f->eoln = 1;
	}
	else
		f->eoln = 0;
	f->eof=feof(f->fp);
}

// open filename n (null for a temporary file) of size l (-1 if unknown)
// with mode m ("w"rite or "r"ead);
FILE *
Fopen (char *n, int l, char *m)
{
  FILE *f;
  register char *s;
  static char ch = 'A';
  static char tmp[MAXFILENAME];


  if (n == NULL) // create a temporary file
    sprintf (tmp, "/tmp/ptc%d%c", getpid (), ch++);
  else
    {
      if (l < 0)
	l = strlen (n);
      strncpy (tmp, n, MAXFILENAME);
      for (s = &tmp[MAXFILENAME - 1];
	   *s == ' ' || *s == '\0' || s - tmp > l;)
	*s-- = '\0';
      if (tmp[MAXFILENAME - 1] != '\0')
	{
	  (void) fprintf (stderr, "Too long filename '%s'\n", n);
	  exit (1);
	}
    }
  s = tmp;
  if ((f = fopen (s, m)) == NULL)
    {
      (void) fprintf (stderr, "Cannot open: %s\n", s);
      exit (1);
    }
  if (n == NULL)
    unlink (tmp);
  return (f);
}
