termptr osnmod(termptr, int, int);/*12/7/94*/   

termptr signseq(termptr *,int, int, char, bool, bool);/*12/7/94*/
termptr signsnseq(int, termptr);
ocharptr signseqgr(ocharptr, int, int);
ocharptr signseqson(ocharptr, int, int);
ocharptr signseqsp2k(ocharptr, int, int);
ocharptr signseqon(ocharptr, int, int);
ocharptr signseqso2k(ocharptr, int, int);

void multsplit(termptr);
void spinsplit(ocharptr);

void chkskcmpt(termptr, termptr *);

termptr attach(termptr, termptr, bool);
ocharptr rattach(ocharptr, ocharptr, bool);

termptr insert(termptr, termptr);

/*termptr linsert(termptr, termptr);*/

termptr deskewm(termptr *,int ,int);

termptr plethysm(termptr, termptr);
termptr pleth(frame, frame);
termptr listplethlist(termptr, termptr);
termptr listplethterm(termptr, termptr, bool);

ocharptr change(ocharptr, int, groop, groop);

bool readauto(string0 *,int *,int *, groop *,groop *, int);

void cancel(int);

void saveframe(text *, frame);
void loadframe(text*, int, frame *);

void savecom(string0, int);
void loadcom(string0, int);

void contragsplit(ocharptr);
void conjsplit (termptr list);

void cf(grptype, int, int);

termptr plethfactor(termptr, termptr);

ocharptr upqseq(ocharptr,int,int, int);

termptr expleth(termptr);
