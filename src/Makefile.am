
GITBUILD=`$(top_srcdir)/build-aux/git-version-gen $(top_srcdir)/.tarball-version`

noinst_HEADERS = version.H
BUILT_SOURCES = version.H
EXTRA_DIST = version.stub

target_not_existing_always_run:
	@echo >/dev/null

version.H: version.stub target_not_existing_always_run
	@cp version.stub version.H
	@echo \#define GITBUILD \"${GITBUILD}\" >>version.H
	@echo \#define MIRALIBVERSION GITBUILD >>version.H
	@echo \#define MIRAVERSION MIRALIBVERSION \" \" MLV_QUIET MLV_TIME\"\" >>version.H


SUBDIRS = \
	debuggersupport\
	io\
	stdinc\
	errorhandling\
	util\
	mira\
	caf\
	memorc\
	modules\
	progs\
	support\
	scripts

#	EdIt\
#	examine\
#	knn_abi373\
#	knn_alf\
#
#	examples_programming\
#