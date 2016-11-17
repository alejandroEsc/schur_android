void orthdim(int, frame ,bframe *, bool);
void sympdim(int, frame, bframe *);

void dimnunit(int, ocharptr, bframe *);
void dimnsymp(int, ocharptr, bframe *);
void dimnso3(ocharptr, bframe *);
void dimnorth(int, ocharptr, bframe *);
void dimnalt(int, ocharptr, bframe *);
void dimnsymm(int, ocharptr, bframe *);
void dimne6(ocharptr, bframe *);
void dimne7(ocharptr, bframe *);
void dimne8(ocharptr, bframe *);
void dimng2(ocharptr, bframe *);
void dimnf4(ocharptr, bframe *);
void dimngrp(int, ocharptr, bframe *);
void dimnprop(ocharptr);

void casmunit(int, ocharptr, bframe *);
void casmsymp(int, ocharptr, bframe *);
void casmorth(int, ocharptr, bframe *);
void casmg2(ocharptr, bframe *);
void casmf4(ocharptr, bframe *);
void casme6(ocharptr, bframe *);
void casme7(ocharptr, bframe *);
void casme8(ocharptr, bframe *);
void casmgrp(int, ocharptr, bframe *);
void casmnsymp(int, ocharptr, bframe *);

void d2index(int, ocharptr, bframe *);

termptr dynklbl(int, ocharptr);
void dynkinn(int, int, ocharptr, ocharptr, bframe *);

ocharptr partlbl(int, termptr);

void putprop(ocharptr, bool);
void putdindex(ocharptr);

void utrace(prodtype);

void pproperties(prodtype);

void adjoint(frame *);

void todynk(ocharptr, frame *, int);
void frdynk(ocharptr *, frame, int);

void diml168(ocharptr, bframe *);
