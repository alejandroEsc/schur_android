# include <stdio.h>
# include <sys/types.h>
# include <unistd.h>

# include "define.h"

//# ifdef READONLY
#define Rmode  "r"
//# else
//#define Rmode  "r+"
//# endif

# define  Finish(f) ((f).out && !(f).eoln) ? (Putchr('\n', f), 0) : 0, !fseek((f).fp, 0L, 0)
# define  Closex(f) (f).init = ((f).init ? (Finish(f), fclose((f).fp), 0) : 0), (f).fp = NULL

/*
 * This has been replaced by function Resetx2
# define  Resetx(f, n, l) (f).init = (f).init ? (Finish(f)) : (((f).fp = Fopen(n, l, Rmode)), 1), (f).eof = (f).out = 0, Getx(f)
*/

# ifdef WRITEONLY
#define Wmode 	"w"
# else
#define Wmode 	"w+"
# endif

//# define  Rewrite(f, n, l) (f).init = (f).init ? !fseek((f).fp, 0L, 0) : (((f).fp = Fopen(n, l, Wmode)), 1), (f).out = (f).eof = 1
# define  Rewritex(f, n, l) (f).init = (f).init ? (Finish(f)) : (((f).fp = Fopen(n, l, Wmode)), 1), (f).out = (f).eof = (f).eoln = 1

FILE *Fopen (char *, int, char *);
void Resetx2(text *, char *, int );


