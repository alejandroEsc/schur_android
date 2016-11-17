# ifndef DEFINE_H
# define DEFINE_H 1
# include "standard.h"
extern text output;
extern text input;

# define  Fread(x, f) fread((void *)&x, sizeof(x), 1, f)
# define  Get(f) Fread((f).buf, (f).fp)
# define  Getx(f) (f).eoln ? (((f).eoln = 0, (f).buf = '\0'), (f).init = 0) : ((f).init = 1, (f).eoln = (((f).buf = fgetc((f).fp)) == '\n') ? (((f).buf = ' '), 1) : 0)
# define  Getchr(f) (f).init ? (f).buf : (Getx(f), (f).buf) , Getx(f)
# define  Fscan(f) (f).init ? ungetc((f).buf, (f).fp) : 0, Tmpfil = (f).fp
# define  Scan(p, a) Scanck(fscanf(Tmpfil, p, a))
  
# define  Eoln(f) ((f).eoln ? true : false)
# define  Eof(f) ((((f).init == 0) ? (Get(f), (f).init = ((f).buf == '\n') ? (((f).buf = ' ', (f).eoln = 1), 1) : 1) : 0, ((f).eof ? 1 : feof((f).fp))) ? true : false)
# define  Fwrite(x, f) fwrite((void *)&x, sizeof(x), 1, f)
  
# define  Put(f) Fwrite((f).buf, (f).fp)
# define  Putx(f) (f).eoln = ((f).buf == '\n'), (void)fputc((f).buf, (f).fp)
# define  Putchr(c, f) (f).buf = (c), Putx(f)
# define  Putl(f, v) (f).eoln = v
  
# define  Close(f) (f).init = ((f).init ? (fclose((f).fp), 0) : 0), (f).fp = NULL

#define min(x,y)  (x<y ? x : y)
#define MIN(x,y)  (x<y ? x : y)
#define MAX(x,y)  (x>y ? x : y)


# ifndef MAXFILENAME
# define  MAXFILENAME 128
# endif
  
# define  Line __LINE__

#endif /* ! DEFINE_H */
