# include <stdio.h>
# include <unistd.h>
# include "define.h"
# include "standard.h"

void
Getl (text * f)
{
  while (f->eoln == 0)
    Getx (*f);
  Getx (*f);
}
