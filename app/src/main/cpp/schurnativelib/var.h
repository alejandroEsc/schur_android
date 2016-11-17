/** \file var.h
 * definition of global variables
 */
#ifndef VAR_H
#define VAR_H 1

char    *dataPath;
char	start, fin, del, bs, cr, esc,
	bel, lt, rt, go, nth, sth,
	backspace;
string0	abuff, buffz, buffer,buffi;
string0	logname;
int	nprod, iii, ppn, ppt,jz, fnn,
	kkz, ppz, comnum, ern, errpos,
	sign, seccount,setlimit, maxb,
	count, rcount, pcount, /*!< for garbage collection, see heapstatus */
	qsn, cutoff, liveprod, liverep, livesfn, plwt, /* plen,*/
	terms, leftb, rightb , columns, 
  	nbCommands;  /*!< number of commands known by SCHUR, computed from help files */
char    *instr;
bframe	primes, nulls;
frame	full,  /*!< partition initialized to be fully filled */
	nolls ;/*!< partition initialized to zeros */
grptype	group;
modes	mode;
char    *displayModes[4],
	prompt[500], 
        colorsBegin[20], /*!< display text in color (depends on OSTYPE)*/
        colorsEnd[20] ,  /*!< back on normal color */
        ostype[100] ;	 /*!< to copy OSTYPE value defined by configure*/

bool	logging, f4load, e6g2load, e6load, e7load, e7e6load,e6f4load,
	e8load, e8soload, e8suload, e8e6load, e6soload, u27e6load,f4g2load,
        e6su3g2load,e8f4g2load,weight,wprod,so26f4load,l168load,
	su56e7load, su248e8load, echo, 
	fnex, /*!< set to true when interpreting a user-defined function, see fnptr*/
	iosup, notdone,
	trace, erred, logopen, rep, racah, more,
	digits, gotline,  bell, qfn, orderlist,
	pow_note, redu, tex, nreduce, xreduce, sslab,
	qspecial, poly, qtest, color, pfn, dimb, sb_conj, sb_ListOutput,
	ibm, mmono, homo, elem, forg, psum, debug_schur, sb_prog, quiet;

text	f4file, e6file, e6g2file, u27e6, su56e7, e7file,e6f4,f4g2,so26f4,
	e7e6file, su248e8, e8file, file1, file2, e8so16,e6su3g2,e8f4g2,
	e8su2e7, e8su3e6, e6so10,l168file, logfile, archive;

ocharptr	f4index, e6index, e6g2index, u27e6index, su56e7index, e7index,
	e7e6index, e8index, discardedchrcs, e8soindex, e8suindex, e8e6index,
        e8f4g2index,so26f4index, lambda1,l168index,
	e6su3g2index,su248e8index, e6soindex,e6f4index,f4g2index, current;
termptr	discardedsfns, scurrent, dynk, spnil;
tabptr	f4tab, e6tab, e6g2tab, u27e6tab, su56e7tab, e7tab,
	e7e6tab, e8tab, e8sotab, e8sutab, e8e6tab, e6sotab,e6f4tab,f4g2tab,
	su248e8tab,e6su3g2tab,e8f4g2tab,so26f4tab,l168tab;
unsigned short	svarindex;
unsigned short	srjctindex;
unsigned short	varindex;
unsigned short	rjctindex;
unsigned short  rjctindex2;
unsigned short	pvarindex;
unsigned short	prjctindex;
groop	ggroup;
groopArray	currgrp;
fnptrs	fnptr, /*!< user defined function currently interpreted if fnex=true*/
	discardedfns;
#endif /* ! VAR_H */
