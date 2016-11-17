#ifndef R_H
#define R_H 1

char more_Question(void);

bool dothis(void);

bool overwrite(string0);

bool findfile(string0);

char skipbl(string0, int *);


termptr reverselist(termptr);
ocharptr rreverselist(ocharptr);
prodtype preverselist(prodtype);

void readacard(text *,string0 *,int *);
void readelt(string0 *,int *,int *);
void readachrc(text *,string0 *,int *, ocharptr *);
void readfilename(string0, int *, string0 *);
void readint(string0 *,int *,int *);
void readpart(string0 *, int *, frame *);
void readchrc(text *,string0 *,int *, ocharptr *);
void readword(string0 *, int *, char *);
void readfnset(text *, fnptrs *);
void readlist(text *,string0 *,int *,termptr *, bool);

//void wrtfrme(text *, frame);
//void writer(text *, int *,int, char, char, caystype , ocharptr, bool);
//void wrttlst(text *, int *,int, char, char, termptr, bool);
//void writeprod(text *, prodtype, int, bool);

int csize(ocharptr);
int tsize(termptr);
int psize(ocharptr);

#endif /* ! R_H */
