GITSHORT=`$(top_srcdir)/build-aux/git-version-gen isshort $(top_srcdir)/.tarball-version`
MIRAVERSION=${GITSHORT}_${build_os}_${BITPACKAGE}_${MIRASTATIC}
MIRADIR=${PACKAGE_TARNAME}_${MIRAVERSION}

BUNDLETARGET=darwin

3rdparty:
	rm -rf 3rdparty
	cp -R -L ../src/3rdparty .
	tar cvjf mira_3rdparty_`date +%d-%m-%Y`.tar.bz2 3rdparty

libbundle-nobundle:
	@echo "Do not need to bundle on this platform"

libbundle-darwin:
	@echo "Using Macdylibbundler to create OSX compatible dylib references"
	(cd ${MIRADIR}; dylibbundler -b -od -d lib -p @executable_path/../lib/ -x bin/mira; cd ..)

copybin:
	@echo "Making distribution for ${BITPACKAGE}"
	-rm -rf ${MIRADIR} *.tar *.gz
	mkdir ${MIRADIR}
	mkdir ${MIRADIR}/bin
	cp ../src/progs/mira ${MIRADIR}/bin
	ln -s -f mira ${MIRADIR}/bin/mirabait
	ln -s -f mira ${MIRADIR}/bin/miraconvert
	ln -s -f mira ${MIRADIR}/bin/miramem
#	ln -s -f mira ${MIRADIR}/bin/miraSearchESTSNPs
#	ln -s -f mira ${MIRADIR}/bin/miraclip
	mkdir ${MIRADIR}/scripts
	cp ../src/scripts/fixACE4consed.tcl ${MIRADIR}/scripts
	cp ../src/scripts/fasta2frag.tcl ${MIRADIR}/scripts

copydocs:
	mkdir ${MIRADIR}/docs
	cp -R -L ../doc/docbook/Defi*.html ${MIRADIR}/docs
	( cd ${MIRADIR}/docs; ln -s DefinitiveGuideToMIRA.html index.html)
	cp -R -L ../doc/docbook/*.pdf ${MIRADIR}/docs
	cp -R -L ../doc/docbook/bookfigures ${MIRADIR}/docs
	cp -R -L ../doc/docbook/images ${MIRADIR}/docs
	cp -R -L ../doc/docbook/doccss ${MIRADIR}/docs

distribin: copybin libbundle-${BUNDLETARGET}
	strip ${MIRADIR}/bin/*
	tar cvjf ${PACKAGE_TARNAME}_binonly_${MIRAVERSION}.tar.bz2 ${MIRADIR}

distridbg: copybin libbundle-${BUNDLETARGET}
	rm -rf ${MIRADIR}/scripts
	tar cvjf ${PACKAGE_TARNAME}_dbgonly_${MIRAVERSION}.tar.bz2 ${MIRADIR}

distrib: copybin libbundle-${BUNDLETARGET} copydocs
	cp ../THANKS ${MIRADIR}
	cp README ${MIRADIR}
	cp LICENCE ${MIRADIR}
	cp -R -L ../src/support ${MIRADIR}
	rm -rf ${MIRADIR}/support/Makefile*
#	cp ../src/mira/demoparam.prm ${MIRADIR}/paramtemplates.par
	cp ../src/mira/CHANGES.txt ${MIRADIR}
	cp ../src/mira/CHANGES_old.txt ${MIRADIR}
#	cp -R -L ../minidemo ${MIRADIR}
#	( cd ${MIRADIR}/minidemo/demo1; lndir ../data/scf; lndir ../data/exp_set1; ls -1 *exp >mira_in.fofn; rm fofn)
	tar cvjf ${PACKAGE_TARNAME}_${MIRAVERSION}_fullsymbols.tar.bz2 ${MIRADIR}
	strip ${MIRADIR}/bin/*
	tar cvjf ${PACKAGE_TARNAME}_${MIRAVERSION}.tar.bz2 ${MIRADIR}


clean:
	rm -rf mira_*
	rm -rf 3rdparty*
