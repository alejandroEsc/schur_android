void dsecond(termptr *,frame *,int, int,int, int);
void dfirst(termptr *,frame *,int, int,int, int);

termptr dgenerate(termptr *,frame *);

int dcoeffq(termptr *,frame *);
int icoeffq(termptr *,termptr *);

termptr difchrq(termptr *);

termptr reducedq(termptr);

termptr qouterskew(termptr *,termptr *);

termptr lqouter(termptr *,termptr *);

termptr qinner(termptr, termptr);
termptr lqinner(termptr, termptr);

ocharptr qspin(ocharptr);

termptr qmult(termptr);

termptr nraise(termptr, int, int, int);
termptr yraise(termptr, int, int, int);
termptr ayraise(termptr, int);
void zraise(termptr, int, int,int);
termptr xraise(termptr, int);
termptr craise(termptr, int);
termptr lraise(termptr, int);

termptr oqinner(termptr, termptr);


termptr sqinner(termptr, termptr);
termptr rqinner(termptr, termptr);
termptr lsqinner(frame, frame, int);
termptr lsinner(frame *,frame *,int, int, int);

termptr sladd(termptr, termptr);

void signsort(termptr);
void dsort(termptr *,termptr *,frame);
void nsort(termptr, int);
void mmsort(termptr *, int);
