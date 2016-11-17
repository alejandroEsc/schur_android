void swapgroupname(int, int);
prodtype swapgroups(prodtype, int, int);

ocharptr sp2nk(ocharptr *,int,int);
termptr sp2spun(termptr *,int ,int, int);
termptr spun(termptr *,int ,int);
ocharptr sprun(ocharptr *,int);
ocharptr lsprun(ocharptr *,int);

ocharptr metaplet(int ,int);

ocharptr lspnrsp2on(ocharptr *,int);
ocharptr spnrsp2on(ocharptr *,int);

ocharptr so4so3(ocharptr);
ocharptr so4su2su2(ocharptr);
ocharptr sp2nsu2son(ocharptr);

void fixgroups(prodtype *, int ,int);
void fixsubgroup(int *,int *,int *,int *, int *);

void dobranch(ocharptr *, ocharptr *,int ,int, int,int ,int);

void su1crunch(prodtype *);

void putgroup1(groopArray, int);

prodtype hecke(frame);

prodtype unu1(int);

void schar(ocharptr, ocharptr);

ocharptr removek(ocharptr, int);

termptr character(termptr, termptr);

termptr sremovek(termptr, int);

termptr gensfnlist(int);

termptr snchar(termptr, int);

prodtype macseries(char,bool, bool, bool, bool, bool, int);
prodtype inverseseries(prodtype, int);

ocharptr soncbrun(ocharptr *, int);
ocharptr lsoncbrun(ocharptr *, int);
termptr soncun(termptr *, int, int);
