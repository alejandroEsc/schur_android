termptr jstndise(termptr *);

void nqup(termptr *);

termptr joinmn(termptr *,termptr *);

void jnum(termptr *,int *,int *,int *);

termptr jmod(termptr *);
termptr rmod(termptr *);

int repmult(termptr *,termptr *,termptr *);

bool dead(frame ,termptr ,termptr );
bool deadtest(frame *,termptr *, termptr *);

termptr qsameweight(termptr *,termptr *);
termptr qparts(termptr *);

int rest(frame *, int);
int qmaxin(frame *);

void maker(termptr *, frame *);

void makeray(frame *,qframe *,int *,int *);

bool latticetest(termptr *);

void salam(termptr *, int);
void salam1(termptr *,frame *, int, int, int, int, int, bool *);

void second(termptr *,termptr *, int, int, int, int, int, int, int);
void rsecond(termptr *, frame *, int, int, int, int, int);

void first(termptr *,termptr *, int, int, int, int, int, int);
void rfirst(termptr *,frame *, int, int, int, int);

termptr wgenerate(termptr,termptr,termptr);
termptr rgenerate(termptr *,frame *);

int coeff(termptr,termptr,termptr);
int rcoeff(termptr *,frame *);

int qfactor(termptr,termptr,termptr);

bool platticetest(termptr);

bool paritysequence(termptr);

termptr s_to_q(termptr);
termptr q_to_s(termptr);

termptr indexx(termptr);


