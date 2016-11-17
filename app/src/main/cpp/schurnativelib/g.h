# ifndef G_H
# define G_H 1


termptr getsfn(text *,string0 *, int *);

ocharptr getchrc(text *,string0 *, int *);

prodtype getprod(text *,string0 *, int *);

ocharptr gettrans(text *,string0 *, int *);

void readprod(text *, string0 *,int *,prodtype *);

void allfalse(void);

prodtype pmu(text *,string0 *, int *,int);

prodtype pbranch(string0 *, int *);

prodtype rule(text *,string0 *, int *);

prodtype pautom(text *,string0 *, int *);
# endif  /* ! G_H */
