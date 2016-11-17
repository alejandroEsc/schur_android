void qup(termptr *);

int qlen(frame);
int qqlen(frame);

void schur_psort(prodtype *,bool);
void Qsort(termptr *,bool);
void sort(termptr *,bool);
void osort(ocharptr *,bool);
void sposort(ocharptr *,bool);
void gsort(ocharptr *,bool, caystype);
void insort(int, frame, termptr *,termptr *);

void merge(termptr *,termptr *,bool ,bool);
void cmerge(ocharptr *,ocharptr *,bool, bool);

void pdiv(prodtype *,int);
void rdiv(ocharptr *,int);

void pnu(prodtype *);

void pdisp(prodtype *);
void dispprod(prodtype *);

int ptestord(prodtype, prodtype);
int qtestord(frame *, frame *);

int qcffrmfrm(frame, frame);

void sstkhandle(termptr);
void stackhandle(ocharptr, bool);
void pstackhandle(prodtype);

void cnu(ocharptr *);
void snu(termptr *);

void fdisp(fnptrs *);
void ldisp(termptr *);
void odisp(ocharptr *);

void dispchr(ocharptr *);
void dispsfn(termptr *);

void fnu(fnptrs *);

int blen(string0 *);

int wtfrm(frame *);

int cffrmfrm(frame, frame);
int cfptrfrm(termptr, frame);
int cfptrptr(termptr, termptr);
int ccfptrptr(ocharptr, ocharptr);

void coeffset(termptr *,int ,char);

termptr sameweight(int, int);

void limit(termptr *, int);


int testord(frame *,frame *);
int otestord(ocharptr, ocharptr);
int sotestord(ocharptr, ocharptr);

int stmod(int, int);

termptr swt(termptr);
ocharptr rwt(ocharptr);

void conjgte(frame *);
termptr lconjgte(termptr);
bool sconjgte(frame *);
bool conjtest(frame *);
termptr conjadd(termptr );

void add(termptr *, termptr *);
void oadd(ocharptr *, ocharptr *);
termptr ladd(termptr, termptr);

termptr absval(termptr);

termptr sfncopy(termptr);
ocharptr chrccopy(ocharptr);

prodtype prodcopy(prodtype);

void mergemin(frame, frame, frame *);
void mergemax(frame, frame, frame *);

termptr skcompat(frame);
termptr lskcmpt(termptr);

termptr eqwt(int);

termptr sfnmult(int, termptr);

ocharptr chrcadd(ocharptr, ocharptr);
ocharptr chrcmult(int, ocharptr);

ocharptr sfntochrc(termptr, bool, char);


void pexpand(prodtype *, prodtype *, int);

prodtype prodmult(int, prodtype);
prodtype prodadd(prodtype, prodtype);
prodtype prodexpand(prodtype);

prodtype pcheck(prodtype);

ocharptr gchrccopy(ocharptr);

void coeffsetch(ocharptr *,int ,char);

ocharptr absvalch(ocharptr);
