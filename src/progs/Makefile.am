AM_CPPFLAGS = -I$(top_srcdir)/src $(all_includes)

AM_CXXFLAGS = $(BOOST_CPPFLAGS)

#LINKEDITLIB= -lmiraEdIt -lmiraExamine -lmiraKNN_abi373
#LINKEDITDIR= -L../examine -L../EdIt -L../knn_abi373 -L../knn_alf

LINKEDITLIB=
LINKEDITDIR=


SUFFIXES= .itxt

ciclean:
	-@rm -f compileinfo.itxt

compileinfo.itxt:
	@/bin/echo -n "Compiled by: " >compileinfo.itxt.xxd
	@whoami >>compileinfo.itxt.xxd
	@date >>compileinfo.itxt.xxd
	@/bin/echo -n "On: " >>compileinfo.itxt.xxd
	@uname -a >>compileinfo.itxt.xxd
	@echo "Making .H out of .xxd" $@
	@$(SHELL) -ec 'grep -v "^#" compileinfo.itxt.xxd | xxd -i >compileinfo.itxt.xxd.H'
	@rm compileinfo.itxt.xxd




bin_PROGRAMS= mira
noinst_PROGRAMS= miratest miramer

noinst_HEADERS= quirks.H

if BUILDMIRASTATIC
AM_LDFLAGS=  $(BLOODYDARWINSTATLIBDIRHACK) -L../modules -L../debuggersupport -L../memorc -L../io/ -L../util/ -L../errorhandling/ -L../mira -L../caf $(LINKEDITDIR) $(BOOST_LDFLAGS) -all-static
else
AM_LDFLAGS=  -L../modules -L../debuggersupport -L../memorc -L../io/ -L../util/ -L../errorhandling/ -L../mira -L../caf $(LINKEDITDIR) $(BOOST_LDFLAGS)
endif

if BUILDMIRAMEMORC
MEMORC= -lmemorc
endif

MIRALIBS= -lmiramodules -lmiradebuggersupport $(MEMORC) -lmira -lmiraestass -lmiracaf -lmirautil -lmiradptools -lmirafio -lmiraerrorhandling -lmiradebuggersupport -lmirafmttext $(LINKEDITLIB)

install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	rm -f miramem$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) miramem$(EXEEXT) && \
	rm -f mirabait$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) mirabait$(EXEEXT) &&\
	rm -f miraconvert$(EXEEXT) && \
	$(LN_S) mira$(EXEEXT) miraconvert$(EXEEXT)

#	$(LN_S) miramer$(EXEEXT) miradiff$(EXEEXT) && \
#	$(LN_S) mira$(EXEEXT) miraSearchESTSNPs$(EXEEXT) && \
#	$(LN_S) mira$(EXEEXT) miraclip$(EXEEXT) && \
#	rm -f miraSearchESTSNPs$(EXEEXT) && \
#	rm -f miraclip$(EXEEXT) && \
#	rm -f miradiff$(EXEEXT) && \
#

uninstall-hook:
	cd $(DESTDIR)$(bindir) && \
	rm -f miraSearchESTSNPs$(EXEEXT) && \
	rm -f mirabait$(EXEEXT) && \
	rm -f miramem$(EXEEXT) && \
	rm -f miradiff$(EXEEXT) && \
	rm -f miraclip$(EXEEXT)


strip:
	du -sk $(bin_PROGRAMS)
	strip $(bin_PROGRAMS)
	du -sk $(bin_PROGRAMS)

privinst:
	cp $(bin_PROGRAMS) /home/bach/pub/bin

binclean:
	rm -rf $(bin_PROGRAMS)


mira.C: ciclean compileinfo.itxt
mira_SOURCES= mira.C quirks.C
mira_LDADD= $(MIRALIBS)

miramer.C:
miramer_SOURCES= miramer.C
miramer_LDADD= $(MIRALIBS)

miratest.C:
miratest_SOURCES= miratest.C
miratest_LDADD= $(MIRALIBS)

