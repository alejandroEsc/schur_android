void distinct(termptr *);

int pwtfrm(frame *);

int nlambda(frame);

int hookpart(int, int,frame );

void prodhook(frame, lbframe *);
void qqprodhook(frame, lbframe *);

void qsfn(termptr , lbframe *);
void qsfnlist(termptr);
void qqsfn(termptr, lbframe *);
void qqsfnlist(termptr);

void polyprod(lbframe, lbframe, lbframe *);
void polysum(char,lbframe, lbframe, lbframe *);
void polyn1(char, lbframe *);
void partpoly(lbframe *);
void poly1(int, char, lbframe *);
void poly2(char, lbframe *);
void poly3(char, lbframe *);
void polybp(char, lbframe *);
void polyn3(char, lbframe *);
void polyn4(char, lbframe *);
void poly4(char, lbframe *);
void poly5(char, lbframe *);
void poly6(char, lbframe *);

bool qvalidser(char);

void qallseries(char, lbframe *);
void qqallseries(char, lbframe *);

termptr ntensor(termptr, int);
ocharptr rntensor(ocharptr, int);

termptr fullx(int);
termptr fulling(termptr);
termptr fullsa(int);

ocharptr sp_pleth(ocharptr, ocharptr, groop);
ocharptr unpleth(ocharptr, ocharptr, groop);
ocharptr gpleth(ocharptr, ocharptr, groop);
ocharptr g2pleth(ocharptr, ocharptr);
ocharptr g2p(ocharptr);

void sigmap(frame, frame, frame *, int *);

void hallp(frame, frame);

void kostka(termptr, termptr, bframe *);
void kostkamatrix(int);

ocharptr spncpleth(ocharptr, ocharptr, groop);
ocharptr sprpleth(ocharptr, ocharptr, groop);

ocharptr genunlist(ocharptr, groop, bool);
ocharptr genspnrlist(ocharptr, int, groop);
ocharptr gensprunlist(ocharptr, int, int, bool, groop);

termptr  classlist(termptr);

termptr  s_to_p(termptr);

ocharptr kinsert(ocharptr, int,  groop);

void bpower(int, int, bframe *);

void hclass(frame, bframe *);

ocharptr associate(ocharptr, groop);

ocharptr sponmodify(ocharptr, groop);

ocharptr sprextend(ocharptr, groop);

ocharptr star(ocharptr, groop);

ocharptr soncpleth(ocharptr,ocharptr, groop);

ocharptr gensonrlist(ocharptr, int, groop);
ocharptr gensorunlist(ocharptr, int, int, groop);

void mseries(lbframe *);
void pseries(lbframe *);
void aseries(lbframe *);
void cseries(lbframe *);
void bseries(lbframe *);
void dseries(lbframe *);
void eseries(lbframe *);
void gseries(lbframe *);
void fseries(lbframe *);
void hseries(lbframe *);
void lseries(lbframe *);
void qseries(lbframe *);
void qbseries(lbframe *);
void qaseries(lbframe *);
void qeseries(lbframe *);
void qfseries(lbframe *);
void qgseries(lbframe *);
void qhseries(lbframe *);
void qlseries(lbframe *);
void qmseries(lbframe *);
void qpseries(lbframe *);
void qqqseries(lbframe *);
void bpseries(lbframe *);

void fixxcg(groop *, grptype, int, int);

void writefrme(frame);

int lowestweight(ocharptr);

termptr snredpleth(termptr *);
