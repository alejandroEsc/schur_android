bin_PROGRAMS = schur
schur_SOURCES = define.h dim.h Getl.h g.h gr.h m.h q1.h q2.h r.h s0.h s1.h s2.h s3.h s4.h s5.h s6.h s7.h s8.h s.h standard.h type.h var.h \
		g.c Getl.c gr.c \
		dimensions.h dimensions.c \
		bignums.h bignums.c \
		init.h init.c \
		m.c q1.c q2.c r.c s0.c s1.c s2.c s3.c s4.c s5.c s6.c s7.c s8.c s.c \
		apropos.h apropos.c \
		branch.h branch.c \
		utils.h utils.c \
		mymalloc.h mymalloc.c \
		sets_mgmt.h sets_mgmt.c \
		ReadWrite.h ReadWrite.c\
		Scanck.h Scanck.c \
		skew.h skew.c \
		outerskew.h outerskew.c outerskew_inc.c\
		label.h label.c \
		FrontPage.h FrontPage.c\
		savecomList.h savecomList.c\
		tableaux.h tableaux.c\
		write.h write.c\
		FillCommandTab.h FillCommandTab.c\
		hivesLRcoef.h hivesLRcoef.c\
		rib_to_s.h rib_to_s.c\
		schur.c
		
schur_LDADD = $(LIBREADLINE)

# -fomit-frame-pointer make executable smaller and faster but forbid debugging.
AM_CFLAGS =  -W -Wall -Wstrict-prototypes -Wmissing-prototypes -Wno-parentheses -Wredundant-decls $(AST_NESTED_FUNCTIONS) -DDATADIR=\"$(pkgdatadir)\" -DOSTYPE=\"$(build_os)\"
