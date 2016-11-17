#ifndef S2_H
#define S2_H 1
void partssfn(termptr *, int, int);

void partsrep(ocharptr *, int, int);

void categrze(termptr *, termarray *);

void schur_restrict(termptr *, int, char);
void restrict2(termptr *, int, char, bool);
void rrestrict(ocharptr *, int, char);
void prestrict(prodtype *, int, char);
prodtype pwrestrict(prodtype, int, int, char);

void smult(termptr);

void scollctgarbage(void);
void collectgarbage(void);
void pcollectgarbage(void);

termptr ql(char, int);

termptr xv(char, int);

prodtype crunchup(prodtype, int);

ocharptr prodcon(prodtype);
termptr prodsfn(prodtype);

termptr repsfn(ocharptr);

void unlimit(ocharptr *, int);

void onexp(ocharptr *, int);

ocharptr cspin(ocharptr);

ocharptr muc(ocharptr, groop, int);

ocharptr umax(ocharptr, groop, int);

ocharptr contrag(ocharptr, groop);

termptr leqwt(termptr);

int frank(frame);

termptr outerx(int, frame ,frame , int);
termptr outer(frame, frame);
termptr outer2( frame ,frame , int);
termptr louter(termptr, termptr);
termptr louter2(termptr, termptr, int);
//termptr router(termptr *, termptr *, int);

termptr factorsof(int, frame);

termptr sfntohomo(termptr);
termptr homotosfn(termptr);

termptr elemtosfn(termptr);
termptr sfntoelem(termptr);

termptr elemtohomo(termptr);
termptr homotoelem(termptr);

termptr reduinnprd(frame, frame);

termptr plethonerinner(termptr);

termptr makeweight(int, termptr);

termptr inner(frame, frame);
termptr linner(termptr, termptr);
termptr rinner(termptr, termptr);
ocharptr scalarinner(int, int);

void putchrc(text *, ocharptr, bool);
void putsfn(text *, termptr, bool);
void putprod(prodtype);
void putgroup(groopArray);
void getgroup(string0 *,int ,bool);


ocharptr oformbb(ocharptr, ocharptr, bool, bool);

termptr multmono(termptr *, termptr *, int);

termptr genprod(int *);

termptr sfnmon(termptr *, termptr *, int);

termptr complement(termptr, int, int);

termptr sumsquares(termptr);

#define SEGPTRMAX 64
#endif // S2_H
