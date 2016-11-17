void wrtbigno(text *,int *, bframe * );

void wrtdlst(text *, int *, int, int, termptr);

int frbig(bframe);

void bigadd(bframe, bframe, bframe *);

void bigsubtr(bframe, bframe, bframe *, bool *);

void cmult(bframe, int, bframe *);

void cadd(bframe, int, bframe *);

void factor(int, bframe *);

void fmult(bframe, bframe, bframe *);

void factobig(bframe, bframe *);

void cdiv(bframe, int, bframe *, int *);

void bigtofact(bframe, bframe *);

void fdiv(bframe, bframe, bframe *);

void kmult(bframe, int, bframe *);

void kdiv(bframe, int, bframe *);

void maxscoeff(termptr);

void maxrcoeff(ocharptr);

void maxdcoeff(prodtype);

void summult(ocharptr);

void ssummult(termptr);

void psummult(prodtype);

void tsum(ocharptr);

void tssum(termptr);

void tpsum(prodtype);

void hklth(frame, bframe *);

void tobig(int, bframe *);

void factorialn(int, bframe *);

void sumabssquares(termptr);

void multlist(termptr);

int dsize(termptr);

void sumabssquaresrep(ocharptr);
