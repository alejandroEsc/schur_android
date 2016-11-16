2016-11-15 Alejandro Escobar <jaescobar.cell@gmail.com>
        * Android version started.

2014-05-15 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
        * version 6.10 
	* output of big partitions corrected (an old standing
	* bug)
	* be more strict on non-regression checks
	* reintroduce SB_LISToutput that save vars in a file in
	* a list format which may be
	easier to parse for external applications

2013-10-17 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* version 6.09 
	* some code cleaning and corrections
	* compile on mac Os X (thanks to
	* teake.nutma@aei.mpg.de)
	* option -q (quiet) a bit more quiet ! 

2013-01-11 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* version 6.08 
	* some memory leaks correction
	* some bugs (unitialized array) correction
	* make check correction
	* new option -f filename to use filename as as script

2012-02-12 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* version 6.07
	* rule command correction +some little corrections
	* ok with 64 bits systems since version 6.06a

2008-06-22 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* long lines bug correction +some others (functions
	* etc)

2007-03-14 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add APROPOS command a kind of man -k... PDF Manual
	* use clickable
	(hyperefs) links

2007-03-11 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Sometimes branch command crashes schur. Corrected. 
	* Also corrected some typos in the Example part of the
	* manual. User-defined functions
	(macros) are useable in all modes now (see SETFn, FN,
etc).

2007-02-14 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add RIB_TO_S

2007-02-05 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* prompt without CR !

2007-02-01 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Q_TO_SsymFns-> Q_TO_SDual; Ctrl D ok; output
	* cleaning.

2007-01-28 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add HIVESLRCOEFficient and SB_DEBUG

2006-10-19 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* correction on plethysm computation in some case schur
	* crashes
	* an effort for make check to work on other platforms
	* (Solaris for example!)

2006-10-10 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Steven M. Christensen joined us.
	* an effort to take into account solaris terminals -
	* work on progress!

2006-09-19 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* command completion

2006-08-10 Ron C. King and Franck Butelle
<Franck.Butelle@lipn.univ-paris13.fr>
	* RD_RAISE RD_RAISEI correction  

2006-07-04 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* correct some warnings

2006-05-28 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Conversion of the manual from Plain TeX to LaTeX
	* (bibliography to be
	    improved) relative references, tables, equations.

2006-04-30 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Correction/cleaning of output width of results but
	* TeX or LaTeX
	output still to do. help files limited to 67 char.
width. 
	* Correction of tests (extract_example)

2006-04-22
	* Latex version of the manual (just the begining of a
	* long job!),
	corrections to mac, macm
	* use environment variables LINES and COLUMNS to
	* automatically adapt to
	terminal window size changes (if not presents,
remplaced by 37 and 79 
	    resp.) - for help files the limit is fixed at 79
columns.

2006-04-05
	* Much work on logging capabilities, help files
	* corrected...

2006-03-30 
	* Correction of CONTENT etc display of tableaux, help
	* files corrected,
	variables syntax...

2006-03-23 
	* Correction of racah/fracah, inverseries,
	* SB_LISToutput improved

2006-03-17 Ron C. King and Franck Butelle
<Franck.Butelle@lipn.univ-paris13.fr>
	* much work on help file correction of descriptions and
	* examples.
	* some bugs corrections. Ron officially joined the
	* team.
	* add SB_LISToutput to have the ability to save
	* variables in a Maple
	compatible format - work in progress.

2006-03-11 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr> 
	* change MAXINT to INT_MAX and use limits.h instead of
	* values.h to make
	it compile on cygwin - First windows version (from
schur6 point of
	    view).

2006-03-02 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr> 
	* bug correction, rep mode and group not set, rv
	* command breaks the
	program.

2006-03-02 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
and Frederic
Toumazet <Frederic.Toumazet@lipn.univ-paris13.fr> 
	* with the help of Ron King we find a bug in signseq.
	* change in structure term to have a length, needed
	* when zeros are 
	significant, for example in insert command. This change
implies many 
	small changes and more to come !
	* some more code cleaning

2006-02-25 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* some more code cleaning
	* add NSKEW code,doc and example from an old mail of
	* Brian Weybourne
	* add documentation and examples on HEAPstatus

2006-02-23 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* code cleaning on for loops, may improve execution
	* time

2006-02-04 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* code cleaning on for loops, may improve execution
	* time

2006-02-02 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* some code cleaning.

2006-01-30 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* Appendix A of the manual is now constructed from the
	* help files 

2006-01-25 Frederic Toumazet
<Frederic.Toumazet@lipn.univ-paris13.fr> and Franck Butelle
<Franck.Butelle@lipn.univ-paris13.fr>
	* add STARequivalent, INVerseseries..
	* TABlesOfBranchingRules updated.

2006-01-24 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add SPSTAR, SPLitIntoSpinAndTensor to "valid"
	* examples.

2006-01-18 Frederic Toumazet
<Frederic.Toumazet@lipn.univ-paris13.fr> and Franck Butelle
<Franck.Butelle@lipn.univ-paris13.fr>
	* add F_TO_S, SERIESTErmsThatSkew to "valid" examples.

2006-01-17 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add FPROD, FFPROD, RP_RepOrSfnByWt, GENERIC,
	* HALLpolynomialProduct to
	"valid" examples
	* add readline call to have an history of commands

2006-01-14 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* spec file for rpm package construction

2006-01-10 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add CANcelDatFile, SPREXtend, SPRCH, STD_Qfn,
	* UONETrace, YHooklengths
	to "valid" examples
	* QEXPandSpecialSeries corrected and added to "valid"
	* examples

2006-01-08 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add LSEQ, QSAME, GWT to "valid" examples
	* some corrections in help file reading
	* some output cleaning

2005-12-27 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* add KINSert, QQSERies and QSERies to "valid" examples

2005-12-23 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* fix some memory leaks - but some remains
	* no more warnings

2005-12-15 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* wseq, parity corrected and added to "valid" examples  

2005-12-14 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* autotools corrections:
		* ./configure && make && make install finally
		* works  
		* make distcheck finally works 

2005-12-13 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* purge headers and multiples definitions of some
	* functions
	* boolean -> bool
	* still some warnings (uninitialized variables)

2005-12-01 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* 145 tests validated
	* no more string8, string16, string24

2005-11-06 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr> 
	* first try of autoconf/automake   

2005-11-01 Franck Butelle <Franck.Butelle@lipn.univ-paris13.fr>
	* tiny changes to decrease number of warnings 
  	* use of protoize and indent to clean the code (now
  	* ansi C)	

