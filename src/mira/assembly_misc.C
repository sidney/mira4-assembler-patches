/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 * Copyright (C) 2000 and later by Bastien Chevreux
 *
 * All rights reserved.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 */

// functions to process reads
// currently in namespace and object assembly


#include "boost/unordered_map.hpp"
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "mira/assembly.H"
#include "mira/dataprocessing.H"
#include "mira/hashstats.H"


using namespace std;


//#define CEBUG(bla)   {if(CEBUGFLAG) {cout << bla; cout.flush();}}
#define CEBUG(bla)




//#define CEBUG(bla)   {if(id1==2282 && id2==342) {cout << bla; cout.flush();}}
//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla;}






/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::performHashAnalysis(bool useminkmer, bool dorarekmermask, int32 version, const string prefix, const string postfix, const string logname)
{
  FUNCSTART("void Assembly::performHashAnalysis()");
  //CEBUG("BEFORE\n");
  //for(uint32 actid=0; actid<AS_readpool.size(); actid++){
  //  Read & r=AS_readpool.getRead(actid);
  //  r.integrityCheck();
  //  Read::setCoutType(Read::AS_TEXT);
  //  cout << r;
  //}

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
  hashstatistics_parameters const & hs_params= AS_miraparams[0].getHashStatisticsParams();

  uint8 basesperhash=skim_params.sk_basesperhash;
  if(sizeof(uint64) < 8 && basesperhash > 15) basesperhash=15;

  uint32 nastyrepeatratio=0;
  if(AS_needsskimfornastyrepeats && hs_params.hs_masknastyrepeats){
    // BaCh 22.04.2013: MNR tags are now always deleted by assignReadBaseStatistics(), so must be also re-set
    // TODO: check whether AS_needsskimfornastyrepeats is still needed
    //AS_needsskimfornastyrepeats=false;
    AS_needsskimfornastyrepeats=true;
    nastyrepeatratio=hs_params.hs_nastyrepeatratio;
  }

  HashStatistics s3;
  string filenameforhs;
  string stathsfn(buildDefaultCheckpointFileName("static_hashstat.bin"));
  {
    s3.setHashFrequencyRatios(hs_params.hs_freqest_minnormal,
			      hs_params.hs_freqest_maxnormal,
			      hs_params.hs_freqest_repeat,
			      hs_params.hs_freqest_heavyrepeat,
			      hs_params.hs_freqest_crazyrepeat,
			      hs_params.hs_nastyrepeatratio,
			      hs_params.hs_nastyrepeatcoverage);

    s3.setAvgHashFreqMinimum(hs_params.hs_freq_covestmin);

    vector<uint32> minkmer;
    if(useminkmer){
      for(auto & mp : AS_miraparams){
	minkmer.push_back(mp.getAssemblyParams().as_clipmask_rarekmers);
      }
    }

    if(fileExists(stathsfn)){
      s3.loadHashStatistics(AS_readpool,stathsfn,basesperhash);
    }else{
      s3.prepareHashStatistics(AS_miraparams[0].getDirectoryParams().dir_tmp,
			       AS_readpool,
			       true,
			       false,
			       false,
			       true,
			       1,
			       basesperhash,
			       hs_params.hs_million_hashes_per_buffer,
			       hs_params.hs_rare_kmer_early_kill,
			       filenameforhs
	);
    }
    s3.showHashStatisticsInfo();

    // very rough estimator of coverage for contigs
    if(basesperhash>=15){
      AS_hashstat_avghashfreq=s3.getAvgHashFreqRaw();
      cout << "Estimator of average coverage: " << AS_hashstat_avghashfreq << endl;
      AS_assemblyinfo.setLargeContigCovForStats(AS_hashstat_avghashfreq/3);
    }

    cout << "Assigning statistics values:\n";
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    {
      uint32 mincountkmerforks=2;
      if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	mincountkmerforks=AS_hashstat_avghashfreq/6; // /6 to conservative??
	if(mincountkmerforks<2) mincountkmerforks=2;
      }
      s3.assignReadBaseStatistics_MultiThread(skim_params.sk_numthreads, nastyrepeatratio>0,minkmer,mincountkmerforks);
    }
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

    if(AS_miraparams[0].getAssemblyParams().as_buntify_reads){
      AS_dataprocessing.buntifyReadsByHashFreq_Pool(AS_readpool,basesperhash);
      AS_dataprocessing.addKMerForkTags_Pool(AS_readpool,basesperhash);
    }
  }

  //CEBUG("AFTER\n");
  //for(uint32 actid=0; actid<AS_readpool.size(); actid++){
  //  Read & r=AS_readpool.getRead(actid);
  //
  //  Read::setCoutType(Read::AS_TEXT);
  //
  //  if(r.getName()=="FF5UQ0101A62BE.fn"
  //     || r.getName()=="FFPHEER01DATWZ"
  //     || r.getName()=="FFPHEER01AK3C0"){
  //    cout << r;
  //  }
  //}

  //if(nastyrepeatratio){
  if(hs_params.hs_repeatlevel_in_infofile){
    string filename;

    if(logname.size()){
      filename=buildFileName(version, prefix, postfix, logname, "");
    }else{
      //filename=buildFileName(version, prefix, postfix,
      //			     as_fixparams.as_outfile_stats_readrepeats,
      //			     ".lst");

      //filename=buildDefaultInfoFileName(version, prefix, postfix,
      filename=buildDefaultInfoFileName(-1, "", "",
					"",
					as_fixparams.as_outfile_stats_readrepeats,
					".lst",
					true);
    }

    cout << "Writing read repeat info to: " << filename << " ... ";
    cout.flush();

    uint32 howmanys=0;
    uint32 howmanyt=0;
    uint32 repanalysislevel=hs_params.hs_repeatlevel_in_infofile;
    if(repanalysislevel<5) repanalysislevel=5;
    if(repanalysislevel>8) repanalysislevel=8;

    ofstream fout;
    fout.open(filename.c_str(), ios::out|ios::trunc);
    for(uint32 rpi=0; rpi<AS_readpool.size(); rpi++){
      Read & actread= AS_readpool.getRead(rpi);
      if(!actread.hasValidData()
	 || !actread.isUsedInAssembly()) continue;
      bool mustshow=false;
      if(actread.hasTag(Read::REA_tagentry_idHAF5,-1)) {
	if(repanalysislevel==5) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idHAF6,-1)) {
	if(repanalysislevel<=6) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idHAF7,-1)) {
	if(repanalysislevel<=7) mustshow=true;
      }else if(actread.hasTag(Read::REA_tagentry_idMNRr,-1)) {
	if(repanalysislevel<=8) mustshow=true;
      }
      if(mustshow){
	bool countedthisseq=false;
	for(uint32 tn=0; tn<actread.getNumOfTags(); tn++){
	  const multitag_t & acttag=actread.getTag(tn);
	  if(acttag.to-acttag.from +1 >= basesperhash){
	    mustshow=false;
	    if(acttag.identifier==Read::REA_tagentry_idHAF5) {
	      if(repanalysislevel==5) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idHAF6) {
	      if(repanalysislevel<=6) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idHAF7) {
	      if(repanalysislevel<=7) mustshow=true;
	    }else if(acttag.identifier==Read::REA_tagentry_idMNRr) {
	      if(repanalysislevel<=8) mustshow=true;
	    }
	    if(mustshow){
	      if(!countedthisseq){
		countedthisseq++;
		++howmanys;
	      }
	      ++howmanyt;
	      fout << actread.getName() << '\t'
		   << acttag.getIdentifierStr() << '\t';
	      for(uint32 readpos=acttag.from; readpos<=acttag.to; readpos++){
		fout << static_cast<char>(toupper(actread.getBaseInSequence(readpos)));
	      }
	      fout << '\n';
	    }
	  }
	}
      }
    }

    cout << howmanys << " sequences with " << howmanyt << " masked stretches." << endl;
  }

  // quick check for estimated coverage in genome data
  if(!AS_donequickdenovocoveragecheck){
    // atm only for de-novo, think about doing it for mapping
    // (though there may be good reasons for high coverage in mappings)
    if(basesperhash>=17
       && warnAtHighCoverages(AS_hashstat_avghashfreq)
       && AS_miraparams[0].getNagAndWarnParams().nw_check_coverage==NWSTOP){
      MIRANOTIFY(Notify::FATAL,"High average coverage detected, see output log above respectively the 'WARNING' files in the info directory for more information. In case you wish to force MIRA to disregard this safety check, consider using '-NW:cac=warn' or '-NW:cac=no'");
    }
    AS_donequickdenovocoveragecheck=true;
  }

  if(dorarekmermask) {
    string tmpfname;
    tmpfname=buildFileName(version,"","",
			   as_fixparams.as_tmpf_clippings,
			   ".txt","",false);
    AS_dataprocessing.startLogging(tmpfname,true);

    AS_dataprocessing.performRareKMERMasking_Pool(AS_readpool,"pre-assembly");
    AS_dataprocessing.stopLogging();
  }

  if(0){
    AS_dataprocessing.performKMERRepeatTagging_Pool(AS_readpool,basesperhash);
  }

  if(0){
    Read::setCoutType(Read::AS_TEXT);
    cout << "AFTER2\n";
    for(uint32 actid=0; actid<AS_readpool.size(); ++actid){
      Read & r=AS_readpool.getRead(actid);
      cout << r;
    }
    exit(999);
  }

  if(hs_params.hs_masknastyrepeats && hs_params.hs_apply_digitalnormalisation){
    if(!fileExists(stathsfn)){
      //cout << "Renaming / moving static hash statistics " << filenameforhs << " to " << stathsfn << " ... "; cout.flush();
      cout << "Renaming / moving static hash statistics ... "; cout.flush();
      fileRename(filenameforhs,stathsfn);
      cout << "done." << endl;
    }
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    cout << "Performing digital normalisation: "; cout.flush();
    AS_dataprocessing.performDigitalNormalisation_Pool(AS_readpool,s3,&AS_debrisreason);
    cout << "done\n";
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  }


  if(AS_logflag_dumphashanalysis){
    string logfilename=buildFileName(version, "", "",
				     "elog.dp.hashanalysis_pass",
				     ".lst");
    //string logfilename=AS_miraparams[0].getDirectoryParams().dir_tmp+"/elog.dp.hashanalysis.lst";

    cout << "elog hashan: " << logfilename << endl;
    ofstream logfout;
    logfout.open(logfilename.c_str(), ios::out|ios::trunc);

    for(uint32 rpi=0; rpi<AS_readpool.size(); rpi++){
      Read::setCoutType(Read::AS_TEXT);
      logfout << AS_readpool[rpi];
    }
  }

  FUNCEND();
  return;
}









/*************************************************************************
 *
 * expects reads to have baseflags set  (by performHashAnalysis())
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint64 Assembly::performNewProposedCutbackClips(const string & logname, const string & logprefix)
{
  FUNCSTART("void Assembly::performProposedCutbackClips(const string & logname, const string & logprefix)");

  bool doit=false;
  for(auto mpi=0; mpi < AS_miraparams.size(); ++mpi){
    doit|=AS_seqtypespresent[mpi] && AS_miraparams[mpi].getAssemblyParams().as_clip_proposeendclips;
  }

  if(!doit) return 0;

  cout << "Hash analysis for proposed cutbacks:";

  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  hashstatistics_parameters const & hs_params= AS_miraparams[0].getHashStatisticsParams();

  {
    uint8 basesperhash=as_fixparams.as_clip_pec_basesperhash;
    if(sizeof(uint64) < 8 && basesperhash > 15) basesperhash=15;

    HashStatistics s3;

    s3.setHashFrequencyRatios(hs_params.hs_freqest_minnormal,
			      hs_params.hs_freqest_maxnormal,
			      hs_params.hs_freqest_repeat,
			      hs_params.hs_freqest_heavyrepeat,
			      hs_params.hs_freqest_crazyrepeat,
			      hs_params.hs_nastyrepeatratio,
			      hs_params.hs_nastyrepeatcoverage);

    vector<uint32> dummy;

    string filenameforhs;
    s3.prepareHashStatistics(AS_miraparams[0].getDirectoryParams().dir_tmp,
			     AS_readpool,
			     true,
			     false,
			     false,
			     true,
			     AS_miraparams[0].getAssemblyParams().as_clip_pec_mkfr,
			     basesperhash,
			     hs_params.hs_million_hashes_per_buffer,
			     hs_params.hs_rare_kmer_early_kill,
			     filenameforhs
      );
    s3.showHashStatisticsInfo();
    cout << "Assigning statistics values:\n";
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    s3.assignReadBaseStatistics_MultiThread(skim_params.sk_numthreads, false, dummy,0);
    if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);


    auto avgcov=s3.getAvgHashFreqRaw();

    if(basesperhash>=17
       && AS_miraparams[0].getAssemblyParams().as_clip_pec_mkfr <2){
      if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
	 && AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA] && avgcov >=50){
	cout << "Detected probable higher coverage in Illumina genome project, setting: -CL:pmkfr=2\n";
	const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams()).as_clip_pec_mkfr=2;
      }
    }

    // quick check for estimated coverage in genome data
    if(!AS_donequickdenovocoveragecheck){
      // atm only for de-novo, think about doing it for mapping
      // (though there may be good reasons for high coverage in mappings)
      if(basesperhash>=17
	 && warnAtHighCoverages(avgcov)
	 && AS_miraparams[0].getNagAndWarnParams().nw_check_coverage==NWSTOP){
	MIRANOTIFY(Notify::FATAL,"High average coverage detected, see output log above respectively the 'WARNING' files in the info directory for more information. In case you wish to force MIRA to disregard this safety check, consider using '-NW:cac=warn' or '-NW:cac=no'");
      }
      // Nope, do not set that here. The hash statistics of the main loop have the final word!
      // AS_donequickdenovocoveragecheck=true;
    }

  }

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname.c_str(), ios::out|ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  if(as_fixparams.as_dateoutput) dateStamp(cout);
  cout << '\n';

  static string ggcproblem="ggc";

  cout << "Looking for proposed cutbacks ... "; cout.flush();
  Read::setCoutType(Read::AS_TEXT);

  uint32 cbleft=0;
  uint32 cbright=0;
  uint32 killed=0;
  uint64 numbasesclipped=0;
  for(uint32 actid=0; actid<AS_readpool.size(); actid++){
    Read & r=AS_readpool.getRead(actid);

    if(r.hasValidData()
       && AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_proposeendclips
       && r.hasBaseHashStats()
       && !(r.isBackbone()
	    || r.isRail())){

      CEBUG("lar " << r.getName() << " ");

      bool hasbeenclipped=false;

      uint32 oldlen=r.getLenClippedSeq();

      {
	int32 lpos=r.getLeftClipoff();
	vector<Read::bposhashstat_t>::const_iterator bhsI=r.getBPosHashStats().begin();
	advance(bhsI,lpos);
	for(int32 lend=static_cast<int32>(r.getLenClippedSeq()); lpos<lend; ++lpos, ++bhsI) {
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq >0
	     && (bhsI->fwd.getFrequency() > AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq
		 || bhsI->rev.getFrequency() > AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_ffreq)) {
	    CEBUG("ffreq stop at " << lpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_ffr
	     && ( bhsI->fwd.hasConfirmedFwdRev()
		  || bhsI->rev.hasConfirmedFwdRev())) {
	    CEBUG("ffore stop at " << lpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_fcmst
	     && ( bhsI->fwd.hasConfirmedMultipleSeqType()
		  || bhsI->rev.hasConfirmedMultipleSeqType())) {
	    CEBUG("fcmst stop at " << lpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_fsalp
	     && ( bhsI->fwd.hasSeenAtLowPos()
		  || bhsI->rev.hasSeenAtLowPos())) {
	    CEBUG("fsalp stop at " << lpos << "\n");
	    break;
	  }
	}

	if(lpos != r.getLeftClipoff()){
	  hasbeenclipped=true;

	  if(lpos>0 && lpos>r.getLenSeq()) lpos=r.getLenSeq();
	  CEBUG("pcb l: " << r.getName() << " " << r.getLeftClipoff()
		<< " " << lpos << endl);
	  logfout << logprefix << " left "
		  << r.getName() << '\t'
		  << r.getLeftClipoff() << " -> ";
	  if(lpos==r.getLenSeq()){
	    r.setRQClipoff(r.getLeftClipoff());
	    logfout << "killed\n";
	  }else{
	    r.setLQClipoff(lpos);
	    logfout << r.getLeftClipoff() << '\n';
	  }
	  cbleft++;
	}
      }

      {
      	int32 rpos=r.getRightClipoff();
	vector<Read::bposhashstat_t>::const_iterator bhsI=r.getBPosHashStats().begin();
	advance(bhsI,rpos);

	for(int32 rend=r.getLeftClipoff(); rpos >rend; --rpos){
	  --bhsI;
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq
	     && (bhsI->fwd.getFrequency() > AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq
		 || bhsI->rev.getFrequency() > AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bfreq)) {
	    CEBUG("bfreq stop at " << rpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bfr
	     && ( bhsI->fwd.hasConfirmedFwdRev()
		  || bhsI->rev.hasConfirmedFwdRev())) {
	    CEBUG("bfore stop at " << rpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bcmst
	     && ( bhsI->fwd.hasConfirmedMultipleSeqType()
		  || bhsI->rev.hasConfirmedMultipleSeqType())) {
	    CEBUG("bcmst stop at " << rpos << "\n");
	    break;
	  }
	  if(AS_miraparams[r.getSequencingType()].getAssemblyParams().as_clip_pec_bsalp
	     && ( bhsI->fwd.hasSeenAtLowPos()
		  || bhsI->rev.hasSeenAtLowPos())) {
	    CEBUG("fsalp stop at " << rpos << "\n");
	    break;
	  }
      	}

      	if(rpos != r.getRightClipoff()){
	  hasbeenclipped=true;

      	  CEBUG("pcb r: " << r.getName() << " " << r.getRightClipoff()
      		<< " " << rpos << endl);
      	  logfout << logprefix << " right "
      		  << r.getName() << '\t'
      		  << r.getRightClipoff() << " -> ";
      	  r.setRQClipoff(rpos);
      	  cbright++;
      	  logfout << r.getRightClipoff() << '\n';

	  // special handling of Solexa GGC.G error
	  // from point of right clip, 15 bases backwards:
	  //  search for first ggc.g and clip there
	  if(r.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
	     && AS_miraparams[0].getAssemblyParams().as_clip_pec_sxaggcxg
	     && r.getLenClippedSeq() >=15){
	    //Read::setCoutType(Read::AS_TEXTSHORT);
	    //cout << r;
	    string searchstr=r.getSeqAsChar();
	    boost::to_lower(searchstr);
	    int64 searchstart=r.getRightClipoff()-15;
	    if(searchstart<0) searchstart=0;
	    size_t found;
	    do{
	      found=searchstr.find(ggcproblem,searchstart);
	      if (found!=string::npos){
		searchstart=found+1;
		if(found < r.getRightClipoff()
		   && found+4<r.getRightClipoff()
		   && searchstr[found+4]=='g'){
		  logfout << logprefix << "possible Solexa GGC.G problem "
			  << r.getName() << '\t' << r.getRQClipoff() << " -> ";
		  r.setRQClipoff(static_cast<int32>(found+4));
		  logfout << r.getRQClipoff() << '\n';
		  found=string::npos; // stop the loop
		}
	      }
	    }while(found!=string::npos);
	  }
      	}
      }

      if(hasbeenclipped){
	CEBUG("clipstat yes\n" << r << endl);
	numbasesclipped+=oldlen-r.getLenClippedSeq();
	if(oldlen
	   && (r.getLenClippedSeq() < AS_miraparams[r.getSequencingType()].getAssemblyParams().as_minimum_readlength )){
	  killed++;
	  logfout << logprefix << " "
		  << r.getName() << " killed, remaining length ("
		  << r.getLenClippedSeq() << ")\n";
	}
      }else{
	CEBUG("clipstat no\n" << r << endl);
      }
    }
  }

  logfout.close();

  cout << "done.\nPerformed clips:"
       << "\n\tNum reads cliped left: " << cbleft
       << "\n\tNum reads cliped right: " << cbright
       << "\n\tNum reads completely killed: " << killed
       << "\n\tTotal bases clipped         : " << numbasesclipped
       << "\n\n";

  // now, set the align parameters to enforce clean ends
  for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
    align_parameters & alpar=const_cast<align_parameters &>(AS_miraparams[st].getAlignParams());
    alpar.ads_enforce_clean_ends=true;
    alpar.ads_clean_end_distance=AS_miraparams[0].getSkimParams().sk_basesperhash;
  }


  AS_dataprocessing.clipPolyBaseAtEnd_Pool(AS_readpool,logprefix);

  FUNCEND();

  return numbasesclipped;
}
//#define CEBUG(bla)






/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::cutBackPossibleChimeras(const string & logname, const string & logprefix, const vector<int32> & chuntleftcut, const vector<int32> & chuntrightcut, vector<bool> & chimeracutflag)
{
  FUNCSTART("void Assembly::cutBackPossibleChimeras(const string & logname, const string & logprefix, const vector<int32> & chuntleftcut, const vector<int32> & chuntrightcut)");

  BUGIFTHROW(chuntleftcut.size()!=chuntrightcut.size() && chuntleftcut.size() != AS_readpool.size(),"Arrays mismatch? chuntleftcut.size()!=chuntrightcut.size && chuntleftcut.size() != AS_readpool.size()");

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname.c_str(), ios::out|ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  cout << "Cutting back possible chimeras ... "; cout.flush();

  if(!chimeracutflag.empty()){
    chimeracutflag.clear();
    chimeracutflag.resize(chuntleftcut.size(),false);
  }

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  for(uint32 actreadid=0;actreadid<chuntleftcut.size();actreadid++){
    Read & actread=AS_readpool.getRead(actreadid);
    if(actread.hasValidData()
       && !(actread.isBackbone()
	    || actread.isRail())){
      bool didcut=false;
      if(as_fixparams.as_clip_skimchimeradetection
	 && (chuntleftcut[actreadid]>0
	     || chuntrightcut[actreadid]>0)){
	logfout << logprefix << " possible chimera: " << actread.getName()
		<< "\t["
		<< actread.getLeftClipoff()
		<< ","
		<< actread.getRightClipoff()
		<< "[ using cfrag " << chuntleftcut[actreadid] << ":" << chuntrightcut[actreadid]
		<< " cut back to ";

	actread.setLSClipoff(actread.getLeftClipoff()+chuntleftcut[actreadid]);
	actread.setRSClipoff(actread.getLeftClipoff()+(chuntrightcut[actreadid]-chuntleftcut[actreadid])+1);
	didcut=true;
	if(!chimeracutflag.empty()){
	  chimeracutflag[actreadid]=true;
	}

	logfout << '['
		<< actread.getLeftClipoff()
		<< ","
		<< actread.getRightClipoff()
		<< "[\n";
      }

      if(!didcut
	 && (chuntleftcut[actreadid]<0
	     || chuntrightcut[actreadid]<0)){
	if(as_fixparams.as_clip_skimjunkdetection){
	  logfout << logprefix << " removed possible junk: " ;
	}else{
	  logfout << logprefix << " untouched possible junk: " ;
	}
	logfout << actread.getName()
		<< "\t["
		<< -chuntleftcut[actreadid]
		<< ","
		<< -chuntrightcut[actreadid]
		<< '\n';
	if(as_fixparams.as_clip_skimjunkdetection){
	  actread.setLSClipoff(actread.getLeftClipoff()-chuntleftcut[actreadid]);
	  actread.setRSClipoff(actread.getRightClipoff()+chuntrightcut[actreadid]);
	  if(!chimeracutflag.empty()){
	    chimeracutflag[actreadid]=true;
	  }
	}
      }
    }
  }

  cout << "done.\n";
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
/*
void Assembly::performPool_AdaptorRightClip(const string & logname, const string & logprefix, const uint8 seqtype)
{
  FUNCSTART("void Assembly::performPool_SolexaAdaptorRightClip(const string & logname, const string & logprefix, const uint8 seqtype);)");

// BOOST: regex not compatible with _GLIBCXX_DEBUG
#ifdef _GLIBCXX_DEBUG
  cout << "_GLIBCXX_DEBUG not compatible with BOOST::regex :-(\n";
  return;
#endif

  BUGIFTHROW(seqtype>=ReadGroupLib::SEQTYPE_END,"Unknown seqtype " << static_cast<uint16>(seqtype) << "given.");

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  struct masterslavere_t {
    boost::regex masterre;
    std::list<boost::regex> slaveres;
    bool hasmaster;

    masterslavere_t(): hasmaster(false) {};
  };

  // prepare regular expressions
  list<masterslavere_t> adapres;
  {
    istringstream tmpis;
    if(seqtype==ReadGroupLib::SEQTYPE_SOLEXA){
      static const char regexfile[] = {
#include "adaptorsregex.solexa.xxd.H"
	,0
      };
      tmpis.str(regexfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_IONTORRENT){
      static const char regexfile[] = {
#include "adaptorsregex.iontor.xxd.H"
	,0
      };
      tmpis.str(regexfile);
    }

    masterslavere_t tmpmsre;
    string line;

    while(true){
      getline(tmpis,line);
      if(tmpis.eof()) break;
      if(line[0]=='>'){
	adapres.push_back(tmpmsre);
	line.erase(0,1);         // get away the ">"
	boost::trim(line);
	if(!line.empty()){
	  boost::to_upper(line);
	  adapres.back().masterre=boost::regex(line);
	  adapres.back().hasmaster=true;
	}
      }else{
	BUGIFTHROW(adapres.empty(),"Oooops, found no master expression?");
	boost::to_upper(line);
	adapres.back().slaveres.push_back(boost::regex(line));
      }
    }
  }

  ReadPool adappool(&AS_miraparams);
  {
    istringstream tmpis;

    if(seqtype==ReadGroupLib::SEQTYPE_SOLEXA){
      static const char adapfile[] = {
#include "adaptorsforclip.solexa.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_IONTORRENT){
      static const char adapfile[] = {
#include "adaptorsforclip.iontor.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }else if(seqtype==ReadGroupLib::SEQTYPE_454GS20){
      static const char adapfile[] = {
#include "adaptorsforclip.454.xxd.H"
	,0
      };
      tmpis.str(adapfile);
    }

    string line;
    while(true){
      getline(tmpis,line);
      if(tmpis.eof()) break;
      line.erase(0,1);         // get away the ">"
      if(!line.empty()){
	size_t ereadidx=adappool.provideEmptyRead();
	Read & actread=adappool[ereadidx];
	actread.disallowAdjustments();
	actread.setName(line);
	getline(tmpis,line);
	if(tmpis.eof()) break;
	actread.setSequenceFromString(line);
      }
    }
  }

  //adappool.dumpPoolInfo(cout);

  // Go back if nothing to be searched
  if(adappool.size()==0 && adapres.size()==0) return;

  cout << "Starting " << ReadGroupLib::getNameOfSequencingType(seqtype) << " known adaptor right clip ... "; cout.flush();

  Skim adapskim;
  adapskim.skimStreamPrepare(adappool,7,1);

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "Searching multithread now ... \n"; cout.flush();

  cout << static_cast<int16>(AS_miraparams[0].getSkimParams().sk_numthreads) << endl;

  vector<int32> clipres;
  adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,AS_miraparams[0].getSkimParams().sk_numthreads);
  //adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,1);
  //adapskim.findAdaptorRightClip(AS_readpool,clipres,seqtype,9,8);

  BUGIFTHROW(clipres.size()!=AS_readpool.size(),"clipres.size()!=AS_readpool.size()???");

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname.c_str(), ios::out|ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  uint32 numclipped=0;

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "Searching for " <<  ReadGroupLib::getNameOfSequencingType(seqtype) << " partial end adaptors ... \n"; cout.flush();
  ProgressIndicator<int64> P(0, AS_readpool.size());
  for(uint32 actid=0; actid < AS_readpool.size(); actid++){
    P.progress(actid);
    Read & actread = AS_readpool.getRead(actid);
    if(actread.hasValidData()
       && actread.getSequencingType()==seqtype
       && !(actread.isBackbone() || actread.isRail())){

      auto oldrsclip=actread.getRSClipoff();
      if(clipres[actid]>=0){
	if(clipres[actid] < oldrsclip){
	  ++numclipped;
	  actread.setRSClipoff(clipres[actid]);
	  logfout << logprefix << " "
		  << ReadGroupLib::getNameOfSequencingType(seqtype)
		  << " adaptor: " << actread.getName()
		  << " changed right clip from " << oldrsclip << " to " << clipres[actid] << "\n";
	}
      }else if(!adapres.empty()){
	string seq(actread.getSeqAsChar());
	boost::to_upper(seq);

	boost::match_results<std::string::const_iterator> what;
	boost::match_flag_type flags = boost::match_default;
	std::string::const_iterator start, end;

	for(auto & msre : adapres){
	  bool dosearch=true;
	  if(msre.hasmaster){
	    if(!regex_search(start, end, what, msre.masterre, flags)) {
	      dosearch=false;
	    }
	  }
	  bool breakit=false;
	  if(dosearch){
	    for(auto & thisre : msre.slaveres){
	      start = seq.begin();
	      end = seq.end();
	      if(regex_search(start, end, what, thisre, flags)) {
		if(what.position()< oldrsclip){
		  actread.setRSClipoff(what.position());
		  logfout << logprefix << " "
			  << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
			  << " partial end adaptor: " << actread.getName()
			  << " changed right clip from " << oldrsclip << " to " << what.position() << "\n";
		  breakit=true;
		  break;
		}
	      }
	    }
	  }
	  if(breakit) break;
	}

      }
    }
  }

  P.finishAtOnce();

  if(!logname.empty()){
    logfout.close();
  }

  cout << "done. Clipped " << numclipped << " reads.\n";

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);

  FUNCEND();
  return;
}

*/


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::correctContigs()
{
#ifdef MIRA_HAS_EDIT
  FUNCSTART("void Assembly::correctContigs()");

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\nEditing contigs:" << endl;

  EDITParameters eparams;

  //  eparams.setDoEval();
  eparams.setStrictEvaluation(false);
  eparams.setConfirmationThreshold(0.5);
  eparams.setShowProgress(true);
  eparams.setVerbose(0);
  eparams.setShowProgress(true);


  list<Contig>::iterator I = AS_contigs.begin();
  int32 ccounter=0;
  ProgressIndicator<int64> P(0, AS_contigs.size());

  while(I!=AS_contigs.end()){
    P.progress(ccounter);
    try {
      //	CEBUG("Editing contig:" << ccounter << endl);
      //	CEBUG(*I);
      cout << "Editing contig:" << ccounter << endl;
      editContigBack(*I, eparams);
      ScfBuffer::discard();
      cout << "deleting star columns" << ccounter << endl;
      I->deleteStarOnlyColumns(0, I->getContigLength()-1);
      cout << "marking repeats" << ccounter << endl;

      Contig::repeatmarker_stats_t repstats;
      vector<bool> readsmarkedsrm;
      I->newMarkPossibleRepeats(repstats, readsmarkedsrm);

      //	CEBUG("Corrected contig:" << endl);
      //	CEBUG(*I);
    }
    catch(Notify n){
      n.handleError("Error while examining fault-region");
    }

    I++;ccounter++;
  }

  P.finishAtOnce();

  cout << endl;

  FUNCEND();
#endif
  return;
}






/*************************************************************************
 *
 * Calculates possible sequence vector leftovers at the left side of a read
 * Reads that get a clip must be of Sanger type
 *
 * Does not clip backbone reads, rail reads, multicopyreads
 *  AND not areas protected by Staden GenBank Feature tags
 *
 * Clipping itself must be done afterwards in the performSeqVectorClippings()
 *  function. This was split in two parts to allow releasing of the
 *  big memory chunks AS_readhmcovered, AS_readhitmiss, etc.
 *
 *************************************************************************/


void Assembly::calcPossibleSeqVectorClipoffs(int32 version, const string prefix, const string postfix, const string logname)
{
  FUNCSTART("void Assembly::calcPossibleSeqVectorClipoffs(int32 version, const string prefix, const string postfix, const string logname)");

  if(AS_readhmcovered.size()==0 || AS_readhitmiss.size()==0) {
    cout << "\nNo vector clipping information available, aborting vector clip.\n";
    FUNCEND();
    return;
  }

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\nCalculating possible vector leftovers ... ";
  cout.flush();
  //ProgressIndicator P (0, AS_readhmcovered.size()-1);

  AS_clipleft.clear();
  AS_clipright.clear();
  AS_clipleft.resize(AS_readhmcovered.size(),-1);
  AS_clipright.resize(AS_readhmcovered.size(),-1);

  string filename;
  if(logname.size()){
    filename=buildFileName(version, prefix, postfix, logname, ".txt");
  }else{
    filename=buildFileName(version, prefix, postfix,
			   AS_miraparams[0].getAssemblyParams().as_tmpf_vectorclip,
			   ".txt");
  }

  ofstream logout(filename.c_str(), ios::out | ios::trunc);

  for(uint32 id=0; id<AS_readhmcovered.size(); id++) {
    if(AS_readpool.getRead(id).getSequencingType() != ReadGroupLib::SEQTYPE_SANGER
       || AS_readpool.getRead(id).isBackbone()
       || AS_readpool.getRead(id).isRail()
       || AS_multicopies[id]>0
      ) continue;


    //P.progress(id);

    uint32 clippos=0;
    bool mustclip=false;
    for(uint32 actpos=0; actpos<AS_readhmcovered[id].size(); actpos++) {
      if(actpos-clippos > 5) break;
      if(AS_readhmcovered[id][actpos]>=4) {
	if(AS_readhitmiss[id][actpos]) {
	  if(100.0/static_cast<double>(AS_readhmcovered[id][actpos])*static_cast<double>(AS_readhitmiss[id][actpos]) >= 30.0) {
	    clippos=actpos;
	    mustclip=true;
	  }
	}
      }
    }
    clippos++;

    // check that no GenBank Feature tags protect the area, else clip less
    {

      // FIXME: put all checks for that into read.C (*sigh*)

      for(uint32 i=0; i<AS_readpool.getRead(id).getNumOfTags(); i++){
	const multitag_t & acttag=AS_readpool.getRead(id).getTag(i);
	if(!acttag.isSourceMIRA()){
	  if(acttag.from<clippos) clippos=acttag.from;
	  if(acttag.to<=clippos) clippos=0;
	}
      }
    }

    // auf clip verzichten wenn nur 1 base betroffen (sieht zu doof aus)
    if(mustclip && clippos>1) {
      uint32 maxcliplenallowed=AS_miraparams[AS_readpool.getRead(id).getSequencingType()].getAssemblyParams().as_clip_vector_maxlenallowed;
      if(maxcliplenallowed == 0 || clippos <= maxcliplenallowed) {
	//AS_readpool.getRead(id).setClipoffs(AS_readpool.getRead(id).getLeftClipoff()+clippos,
	//				    AS_readpool.getRead(id).getRightClipoff(),
	//				    false);

	//AS_clipleft[id]=AS_readpool.getRead(id).getLeftClipoff()+clippos;

	AS_clipleft[id]=clippos;

	logout << "Clipped " << clippos << " bases on the left of " << AS_readpool.getRead(id).getName() << "\n";

      } else {
	if(clippos > maxcliplenallowed) {
	  logout << "Not clipped " << clippos << " bases on the left of " << AS_readpool.getRead(id).getName() << " , too long.\n";
	}
      }
    }
  }

  logout.close();

  //P.progress(AS_readhmcovered.size());
  cout << "done.\n";

  AS_steps[ASVECTORSCLIPPED]=1;
  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
}




/*************************************************************************
 *
 * Reads must be Sanger type
 *
 *
 *************************************************************************/

void Assembly::performSeqVectorClippings()
{
  FUNCSTART("void Assembly::performSeqVectorClippings()");

  cout << "\nPerforming vector clipping ... ";
  cout.flush();

  for(uint32 id=0; id<AS_clipleft.size(); id++) {
    if(AS_clipleft[id]>=0
       && AS_readpool.getRead(id).isSequencingType(ReadGroupLib::SEQTYPE_SANGER)) {
      AS_readpool.getRead(id).setClipoffs(AS_readpool.getRead(id).getLeftClipoff()+AS_clipleft[id],
					  AS_readpool.getRead(id).getRightClipoff(),
					  false);
    }
  }
  FUNCEND();

  AS_clipleft.clear();

  cout << "done." << endl;

  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

struct cliplen_t{
  int32 len;
  bool changed;
};


//#define CEBUGFLAG 1
void Assembly::extendADS(int32 version, const string prefix, const string postfix, const string logname)
{
  FUNCSTART("void Assembly::extendADS(int32 version, const string prefix, const string postfix, const string logname)");

//  if(AS_steps[ASADSLISTOK]==0){
//    makeAlignments();
//  }


#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  // TODO: change to use different Aligns / MIRAparams depending
  //   on Sanger / 454 (/ PacBio ???) reads

  // TODO: what about PacBio? currently not used, but should it?

  MIRAParameters tmpparams = AS_miraparams[0];

  const_cast<align_parameters &>(tmpparams.getAlignParams()).al_min_relscore=5;

  assembly_parameters const & as_params= tmpparams.getAssemblyParams();

  string filename;
  if(logname.size()){
    filename=buildFileName(version, prefix, postfix, logname, ".txt");
  }else{
    filename=buildFileName(version, prefix, postfix,
			   as_params.as_tmpf_adsextend,
			   ".txt");
  }

  ofstream logout(filename.c_str(), ios::out | ios::trunc);


  vector<cliplen_t> clips(AS_readpool.size());
  for(uint32 i=0; i<clips.size(); i++){
    clips[i].len=0;
    clips[i].changed=false;
  }

  list<AlignedDualSeq> madsl;

  try{
    // true for using memcache
    Align bla(&tmpparams);

    cout << "\n";
    if(as_params.as_dateoutput) dateStamp(cout);
    cout << "\nSearching possible read extensions (for Sanger and/or 454):\n";

    ProgressIndicator<int32> P(0, static_cast<int32>(AS_adsfacts.size())-1);
    uint32 pindic=0;

    vector<AlignedDualSeqFacts>::const_iterator I = AS_adsfacts.begin();
    for(;I!=AS_adsfacts.end();I++){
      P.progress(pindic++);
      // first try: prolongate to end.
      int32 id1=I->getID1();
      int32 id2=I->getID2();

      // no sense to calc read extensions for reads where both seqtypes are said
      //  not to use extensions
      if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension == false
	 && AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension == false) continue;

      if(AS_permanent_overlap_bans.checkIfBanned(id1,id2)) {
	CEBUG("PermBan for: " << id1 << " " << id2 <<"\tskipping\n");
	continue;
      }

      CEBUG("\n\nid1: " << id1 << "\t" << AS_readpool.getRead(id1).getName() <<endl);
      CEBUG("id2: " << id2 << "\t" << AS_readpool.getRead(id2).getName() <<endl);

      // normally the sequences should have a length >0
      // but due to some clipping being done after SKIM (chimera etc.), it
      //  may happen they are 0 now. If that's the case, don't bother
      //  looking at.
      if(AS_readpool[id1].getLenClippedSeq() == 0
	 || AS_readpool[id2].getLenClippedSeq() == 0) continue;

      // check for sequencing types
      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)) continue;
      // let's allow PacBio HQ

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_IONTORRENT)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_TEXT)) continue;

      if( AS_readpool.getRead(id1).isSequencingType(ReadGroupLib::SEQTYPE_ABISOLID)
	  || AS_readpool.getRead(id2).isSequencingType(ReadGroupLib::SEQTYPE_ABISOLID)) continue;

      //if(clips[id1].changed && clips[id2].changed){
      //	CEBUG(id1 << " and " << id2 <<" already changed.\n");
      //	continue;
      //}

      madsl.clear();

#if CEBUGFLAG > 0
      //Read::setCoutType(Read::AS_TEXT);
      Read::setCoutType(Read::AS_TEXTCLIPS);
      CEBUG(AS_readpool.getRead(id1));
      CEBUG(AS_readpool.getRead(id2));
#endif

      if(I->getSequenceDirection(id1) * I->getSequenceDirection(id2) > 0){

	CEBUG("doalign\n");

	// evil hack warning
	// the &(* ...) construction is needed for gcc3 as it cannot convert
	//  a vector<char> iterator to char *   (*sigh*)

	int32 extendlen1=AS_readpool.getRead(id1).getRightExtend();
	int32 extendlen2=AS_readpool.getRead(id2).getRightExtend();

	if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension == false) {
	  extendlen1=0;
	}
	if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension == false){
	  extendlen2=0;
	}

	CEBUG("l1: " <<AS_readpool.getRead(id1).getLenClippedSeq() << endl);
	CEBUG("e1: " <<extendlen1 << endl);
	CEBUG("l2: " <<AS_readpool.getRead(id2).getLenClippedSeq() << endl);
	CEBUG("e2: " <<extendlen2 << endl);

	if(extendlen1 >= 10 || extendlen2 >= 10){
	  bla.acquireSequences(
	    &(*AS_readpool.getRead(id1).getActualSequence().begin())
	    +AS_readpool.getRead(id1).getLeftClipoff(),
	    AS_readpool.getRead(id1).getLenClippedSeq()+extendlen1,
	    &(*AS_readpool.getRead(id2).getActualSequence().begin())
	    +AS_readpool.getRead(id2).getLeftClipoff(),
	    AS_readpool.getRead(id2).getLenClippedSeq()+extendlen2,
	    id1, id2, 1, 1, true, I->getOffsetInAlignment(id2));
	  bla.fullAlign(&madsl,false,false);

	  if(madsl.size()==0){
	    CEBUG("No results, less radical try.\n");

	    int32 tryseqlen1=0;
	    if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension) {
	      if(clips[id1].changed){
		extendlen1-=clips[id1].len;
	      }
	      extendlen1/=2;
	      tryseqlen1=AS_readpool.getRead(id1).getLenClippedSeq()+extendlen1;
	      if(clips[id1].changed){
		tryseqlen1+=clips[id1].len;
	      }
	      if(AS_readpool.getRead(id1).getLeftClipoff()+tryseqlen1 >= static_cast<int32>(AS_readpool.getRead(id1).getLenSeq())) {
		CEBUG("t1o: " <<tryseqlen1 << endl);
		tryseqlen1=AS_readpool.getRead(id1).getLenClippedSeq()+AS_readpool.getRead(id1).getRightExtend();
		CEBUG("t1n: " <<tryseqlen1 << endl);
	      }
	    }

	    int32 tryseqlen2=0;
	    if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension) {
	      if(clips[id2].changed){
		extendlen2-=clips[id2].len;
	      }
	      extendlen2/=2;
	      tryseqlen2=AS_readpool.getRead(id2).getLenClippedSeq()+extendlen2;
	      if(clips[id2].changed){
		tryseqlen2+=clips[id2].len;
	      }
	      if(AS_readpool.getRead(id2).getLeftClipoff()+tryseqlen2 >= static_cast<int32>(AS_readpool.getRead(id2).getLenSeq())) {
		CEBUG("t2o: " <<tryseqlen2 << endl);
		tryseqlen2=AS_readpool.getRead(id2).getLenClippedSeq()+AS_readpool.getRead(id2).getRightExtend();
		CEBUG("t2n: " <<tryseqlen2 << endl);
	      }
	    }

	    CEBUG("cc1: " <<clips[id1].changed << endl);
	    CEBUG("cl1: " <<clips[id1].len << endl);
	    CEBUG("l1: " <<AS_readpool.getRead(id1).getLenClippedSeq() << endl);
	    CEBUG("t1: " <<tryseqlen1 << endl);
	    CEBUG("cc2: " <<clips[id2].changed << endl);
	    CEBUG("cl2: " <<clips[id2].len << endl);
	    CEBUG("l2: " <<AS_readpool.getRead(id2).getLenClippedSeq() << endl);
	    CEBUG("t2: " <<tryseqlen2 << endl);
	    if(extendlen1 < 5 && extendlen2 < 5) {
	      CEBUG("skip" << endl);
	      continue;
	    }

	    if(tryseqlen1>0 && tryseqlen2>0){
	      bla.acquireSequences(
		&(*AS_readpool.getRead(id1).getActualSequence().begin())
		+AS_readpool.getRead(id1).getLeftClipoff(),
		tryseqlen1,
		&(*AS_readpool.getRead(id2).getActualSequence().begin())
		+AS_readpool.getRead(id2).getLeftClipoff(),
		tryseqlen2,
		id1, id2, 1, 1, true, I->getOffsetInAlignment(id2));
	    }
	  }
	}
      }else{
	if(I->getSequenceDirection(id2)>0){
	}else{
	}
      }

      if(madsl.size()==0){
	CEBUG("No results\n");
      }else{
	int32 bestweight=0;
	list<AlignedDualSeq>::iterator J;
	for(J= madsl.begin(); J!=madsl.end(); ){
	  if(J->isValid()==false){
	    J=madsl.erase(J);
	  }else{
	    if(J->getWeight()>bestweight) bestweight=J->getWeight();
	    J++;
	  }
	}
	// take only the best
	for(J= madsl.begin(); J!=madsl.end();){
	  if(J->getWeight() != bestweight){
	    J=madsl.erase(J);
	  } else {
	    J++;
	  }
	}
//	  cout << "Ext. 1st success: " << id1 << "\t" << id2 << "\n";
//	  cout << *I;
//	  cout << *(madsl.begin());

	int32 lens1=0;
	int32 lens2=0;
	if(madsl.begin()->clipper(as_params.as_readextension_window_len,
				  as_params.as_readextension_window_maxerrors,
				  lens1, lens2)){
//	    cout << "Lalala\n";

	  lens1-=AS_readpool.getRead(id1).getLenClippedSeq();
	  lens2-=AS_readpool.getRead(id2).getLenClippedSeq();
	  CEBUG("o1: " << AS_readpool.getRead(id1).getLenClippedSeq() << "\tn: " << lens1);
	  CEBUG("\no2: " << AS_readpool.getRead(id2).getLenClippedSeq() << "\tn: " << lens2<<endl);


	  if(AS_miraparams[AS_readpool.getRead(id1).getSequencingType()].getAssemblyParams().as_use_read_extension){
	    if(lens1>5 && lens1>clips[id1].len){
	      clips[id1].len=lens1;
	      clips[id1].changed=true;
	    }
	  }

	  if(AS_miraparams[AS_readpool.getRead(id2).getSequencingType()].getAssemblyParams().as_use_read_extension){
	    if(lens2>5 && lens2>clips[id2].len){
	      clips[id2].len=lens2;
	      clips[id2].changed=true;
	    }
	  }
	}
      }
    }
    P.finishAtOnce();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  int32 lenplus=0;
  int32 numchanged=0;
  for(uint32 rid=0; rid<clips.size(); rid++){
    if(AS_readpool.getRead(rid).isBackbone()
       || AS_readpool.getRead(rid).isRail()) continue;
    // contig join spoiler! do not extend back again!
    if(AS_readpool.getRead(rid).hasTag(Read::REA_defaulttag_CJSP.identifier)) continue;
    if(AS_miraparams[AS_readpool.getRead(rid).getSequencingType()].getAssemblyParams().as_use_read_extension) continue;

    if(clips[rid].changed){
      CEBUG("ID: " << rid << "\t" << AS_readpool.getRead(rid).getName() << "\toldlen: " << AS_readpool.getRead(rid).getLenClippedSeq());
      CEBUG("\tgained: " << clips[rid].len << endl);
      numchanged++;
      lenplus+=clips[rid].len;

      logout << AS_readpool.getRead(rid).getName() << "\t" << clips[rid].len << "\n";

      AS_readpool.getRead(rid).setClipoffs(AS_readpool.getRead(rid).getLeftClipoff(),
					 AS_readpool.getRead(rid).getLeftClipoff()+AS_readpool.getRead(rid).getLenClippedSeq()+clips[rid].len-1,
					 false);

      if(AS_readpool.getRead(rid).checkRead()){
	cout << AS_readpool.getRead(rid);
	MIRANOTIFY(Notify::INTERNAL, AS_readpool.getRead(rid).checkRead()) ;
      }
    }
  }

  cout << "\nChanged length of " << numchanged << " sequences."<< endl;
  if(numchanged!=0){
    cout << "Mean length gained in these sequences: " << static_cast<double>(lenplus)/ static_cast<double>(numchanged) << " bases." << endl;
  }

  logout.close();

  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
  return;
}
//#define CEBUGFLAG 0




#define CEBUG(bla)   {cout << bla; cout.flush();}

void Assembly::analyseOverlapHashProfile(vector<uint8> & profile, vector<skimedges_t>::const_iterator seI, ADSEstimator & adse)
{
  vector<uint32> longeststretch(7,0);
  vector<uint32> currentstretch(7,0);

  for(size_t pi=0; pi<profile.size(); pi++){
    //CEBUG(pi << '\t' << static_cast<uint16>(profile[pi]) << '\n');
    for(size_t si=0; si<7; si++){
      if(si==profile[pi]){
	currentstretch[si]++;
	if(currentstretch[si]>longeststretch[si]) longeststretch[si]=currentstretch[si];
      }else{
	currentstretch[si]=0;
      }
    }
  }

  if(longeststretch[3]<5){
    if(AS_skimstaken[seI->skimindex]==true){
      cout << "Remove seI: " << *seI;
      cout << "stretches:\n";
      for(size_t si=0; si<7; si++){
	cout << si << ' ' << longeststretch[si] << endl;
      }

      AS_skimstaken[seI->skimindex]=false;
      AS_numskimoverlaps[seI->rid1]--;
      AS_numskimoverlaps[seI->linked_with]--;
    }
  }
}

#define CEBUG(bla)



/*************************************************************************
 *
 * sorter to sort Contig::templateguessinfo_t from low to high
 *  first on readgroup id
 *  then segment_placement
 *  then on template size
 *
 *************************************************************************/

inline bool Assembly__sortTemplateGuessInfo_(const Contig::templateguessinfo_t & a,
					     const Contig::templateguessinfo_t & b);
inline bool Assembly__sortTemplateGuessInfo_(const Contig::templateguessinfo_t & a,
					     const Contig::templateguessinfo_t & b)
{
  if(a.rgid==b.rgid){
    if(a.splace_seen==b.splace_seen){
      return a.tsize_seen<b.tsize_seen;
    }
    return a.splace_seen<b.splace_seen;
  }
  return a.rgid<b.rgid;
}

/*************************************************************************
 *
 * Destroys AS_templateguesses by sorting it, therefore cleared at end of func
 *
 *
 *************************************************************************/

struct atg_tmp_t {
  uint32 count;
  Contig::templateguessinfo_t tg;
  double mean;
  double stdev;
  double skewness;
  int32 deduced_min;
  int32 deduced_max;

  friend std::ostream & operator<<(std::ostream &ostr, const atg_tmp_t & atgt){
    ostr << "rgid: " << atgt.tg.rgid.getLibId()
	 << "\tc: " << atgt.count
	 << "\tsp: " << static_cast<int16>(atgt.tg.splace_seen)
	 << "\tm: " << atgt.mean
	 << "\td: " << atgt.stdev
	 << "\ts: " << atgt.skewness
	 << "\t-: " << atgt.deduced_min
	 << "\t+: " << atgt.deduced_max;
    return ostr;
  }
};


//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Assembly::analyseTemplateGuesses()
{
//  cout << "tgs: " << AS_templateguesses.size() << endl;
  if(AS_templateguesses.empty()) return;
  sort(AS_templateguesses.begin(),AS_templateguesses.end(),Assembly__sortTemplateGuessInfo_);


  {
    uint32 ind=0;
    for(auto & tge : AS_templateguesses){
      CEBUG(ind
	    << "\trgid: " << tge.rgid.getLibId()
	    << "\tspl: " << static_cast<int16>(tge.splace_seen)
	    << "\tts: " << tge.tsize_seen
	    << '\n'; ++ind);
    }
  }


  vector<atg_tmp_t> atgpred;

  // collapse AS_templateguesses into predictions grouped by rgid and segment_placement
  auto tgI=AS_templateguesses.cbegin();
  for(; tgI!=AS_templateguesses.cend() && tgI->rgid.isDefaultNonValidReadGroupID(); ++tgI) {}
  if(tgI!=AS_templateguesses.cend()){
    while(tgI!=AS_templateguesses.cend()){
      auto tgS=tgI;
      // search end iterator for that grouping
      for(; tgI!=AS_templateguesses.cend() && tgI->rgid==tgS->rgid && tgI->splace_seen == tgS->splace_seen; ++tgI) {};
      auto tgE=tgI;
      auto numvals=tgI-tgS;
      if(numvals>=100){
	atgpred.resize(atgpred.size()+1);
	atgpred.back().count=numvals;
	atgpred.back().tg=*tgS;

	uint64 sum=0;
	for(tgI=tgS; tgI!=tgE; ++tgI) sum+=tgI->tsize_seen;
	double mean=static_cast<double>(sum)/static_cast<double>(numvals);
	double sp2sum=0.0d;
	for(tgI=tgS; tgI!=tgE; ++tgI) {
	  auto m1=static_cast<double>(tgI->tsize_seen - mean);
	  sp2sum+=m1*m1;    // for stdev
	}
	double stdev=sqrt(sp2sum/(numvals-1));
	// skewness, see
	//   http://www.itl.nist.gov/div898/handbook/eda/section3/eda35b.htm
	//   http://www.tc3.edu/instruct/sbrown/stat/shape.htm
	//
	// skewness >0 left skew
	// skewness <0 right skew
	// abs skewness >= 1 very skewed
	//     skewness >= .5 mildly skewed
	//     skewness >= .25 slightly skewed (own definition as we're working with higher number of measurements
	//
	// but for large data sets with tens or hundreds of thousands of measurement, even a few outliers (<1%) on one side
	//  (like happens in misassembled contigs) can induce a skew >1
	// Solution: recalc skew on a subset of values which takes only values at 3*stdev
	// also recalc new stdev on this subset

	sp2sum=0.0d;
	double sp3sum=0.0d;
	double sxstdev=stdev*3;
	uint64 numvals2=0;
	for(tgI=tgS; tgI!=tgE; ++tgI) {
	  auto m1=static_cast<double>(tgI->tsize_seen - mean);
	  if(abs(m1)<=sxstdev){
	    sp2sum+=m1*m1;    // for stdev
	    sp3sum+=m1*m1*m1; // for skewness
	    ++numvals2;
	  }
	}
	stdev=sqrt(sp2sum/(numvals2-1));
	double skewness=0;
	if(stdev>0) skewness=sp3sum/((numvals2-1)*stdev*stdev*stdev);

	// homebrew
	double lefttailfactor=2.0d;
	double righttailfactor=2.0d;
	if(skewness<=-1.0d){
	  lefttailfactor=2.4d;
	  righttailfactor=1.6d;
	}else if(skewness<=-0.5d){
	  lefttailfactor=2.2d;
	  righttailfactor=1.8d;
	}else if(skewness<=-0.25d){
	  lefttailfactor=2.1d;
	  righttailfactor=1.9d;
	}else if(skewness>=1.0d){
	  lefttailfactor=1.6d;
	  righttailfactor=2.4d;
	}else if(skewness>=0.5d){
	  lefttailfactor=1.8d;
	  righttailfactor=2.2d;
	}else if(skewness>=0.25d){
	  lefttailfactor=1.9dd;
	  righttailfactor=2.1d;
	}

	atgpred.back().count=numvals2;
	atgpred.back().mean=mean;
	atgpred.back().stdev=stdev;
	atgpred.back().skewness=skewness;
	if(stdev>0){
	  atgpred.back().deduced_min=(mean-lefttailfactor*stdev > 0) ? (mean-lefttailfactor*stdev) : 0;
	  atgpred.back().deduced_max=mean+righttailfactor*stdev;
	}else{
	  atgpred.back().deduced_min=mean-mean/10;
	  atgpred.back().deduced_max=mean+mean/10;
	}
      }
    }
  }

  cout << "ATG PREDICTIONS\n";
  cout << std::fixed << std::setprecision(10);
  for(auto & ae : atgpred){
    cout << ae << endl;
  }

  // per readgroup, search the prediction with the highest count: that's going to be our final prediction
  auto atgpI=atgpred.cbegin();
  while(atgpI!=atgpred.cend()){
    auto best=*atgpI;
    for(; atgpI!=atgpred.cend() && atgpI->tg.rgid==best.tg.rgid; ++atgpI){
      if(atgpI->count >best.count) best=*atgpI;
    }
    cout << "Final prediction: " << best << endl;
    if(best.tg.rgid.wantSegmentPlacementEstimate()){
      best.tg.rgid.setSegmentPlacementCode(best.tg.splace_seen);
      best.tg.rgid.setWantSegmentPlacementEstimate(false);
      cout << "Set segment placement code.\n";
      AS_guessedtemplatevalues=true;
    } else if(best.tg.rgid.wantTemplateInfoEstimate()){
      best.tg.rgid.setInsizeFrom(best.deduced_min);
      best.tg.rgid.setInsizeTo(best.deduced_max);
      best.tg.rgid.setWantTemplateSizeEstimate(false);
      cout << "Set template size.\n";
      AS_guessedtemplatevalues=true;
    }
  }

  AS_templateguesses.empty();
}
//#define CEBUG(bla)


bool Assembly::warnAtSmileCoverage()
{
  string wstr;
  if(AS_assemblyinfo.huntForSmileCoverage(wstr)){
    AS_warnings.setWarning("CONCOV_SUSPICIOUS_DISTRIBUTION",0,"Suspicious distribution of contig coverages",wstr);
  }
  return !wstr.empty();
}

bool Assembly::warnAtHighCoverages(uint32 measuredcov)
{
  bool retval=false;
  if(AS_miraparams[0].getNagAndWarnParams().nw_check_coverage!=NWNONE
     && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
     && measuredcov>AS_miraparams[0].getNagAndWarnParams().nw_check_covvalue){
    string wstr("You are running a genome ");
    if(AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
      wstr+="mapping";
    }else{
      wstr+="de-novo";
    }
    wstr+=" assembly and the current best estimation for average coverage is "
      +boost::lexical_cast<string>(measuredcov)
      +"x (note that this number can be +/- 20% off the real value). This is ";
    if(measuredcov>=80) wstr+="a pretty high coverage,";
    wstr+="higher than the current warning threshold of "
      +boost::lexical_cast<string>(AS_miraparams[0].getNagAndWarnParams().nw_check_covvalue)
      +"x."
      "\n\nYou should try to get the average coverage not higher than, say, 60x to 100x for Illumina data or 40x to 60x for 454 and Ion Torrent data. Hybrid assemblies should target a total coverage of 80x to 100x as upper bound. For that, please downsample your input data."
      "\n\nThis warning has two major reasons:"
      "\n- for MIRA and other overlap based assemblers, the runtime and memory requirements for ultra-high coverage projects grow exponentially, so reducing the data helps you there"
      "\n- for all assemblers, the contiguity of an assembly can also suffer if the coverage is too high, i.e. you get more contigs than you would otherwise. Causes for this effect can be non-random sequencing errors or low frequency sub-populations with SNPs which become strong enough to be mistaken for possible repeats.";
    if(measuredcov>=150){
      if(measuredcov>=300){
	wstr+="\nA coverage of >300x ... no really, are you kidding me? *sigh*";
      }
      wstr+="\nWith the coverage you currently have, you *really* should downsample your data. You. Have. Been. Warned!";
    }
    AS_warnings.setWarning("ASCOV_VERY_HIGH",1,"Very high average coverage",wstr);
    retval=true;
  }
  return retval;
}



/*************************************************************************
 *
 * mean function, with side-effect and definetely to be used only where
 *  it is currently called from
 *
 * Used to guess which inserts in a contigs are wrong after a first
 *  mapping round and makes sure they get removed
 *
 * E.g.:
 *        bb    cagtcatga***ctgcatgca
 *        r1    cag*catga***ctgcatgca
 *        r2    cag*catga***ctgcatgca
 *        r3    cag*catga***ctgcatgca
 *        ...
 *        rX    cag*catgaTTTctgcatgca
 *
 * with rX being some weird read (maybe low frequency variant, or sequencing
 *  error, or ...)
 *
 * Normally, the two stage mapping would calculate an intermediate new backbone
 *  to be
 *
 *        bbi   cag*catga***ctgcatgca
 *
 * but that is obviously not really the best guess and currently leads to
 *  misalignments with reads ending in this area. E.g.:
 *
 *        bbi   cag*catga***ctgcatgca
 *        rY               actgcatgca
 *
 * instead of
 *
 *        bbi   cag*catga***ctgcatgca
 *        rY            a***ctgcatgca
 *
 * This function will fake the contig CON_counts structure to enable
 *  Contig::deleteStarOnlyColumns() to remove the probably wrong inserts, so
 *  that the new intermediate contig looks like this:
 *
 *        bbi   cag*catgactgcatgca
 *
 * Side-effects:
 *  - CON_counts of contig is probably not reflecting truth anymore
 *  - reads of the contig also get edited, but this routine
 *    is currently used just before all the reads are discarded anyway.
 *
 *************************************************************************/

void Assembly::priv_removePotentiallyWrongBaseInserts(Contig & con)
{
  FUNCSTART("void Assembly::priv_removePotentiallyWrongBaseInserts(Contig & con)");

  // try to look for strains which are only not in backbone
  vector<string> straincons(ReadGroupLib::getNumOfStrains());
  vector<base_quality_t> dummyqual;

  for(uint32 sid=0; sid<ReadGroupLib::getNumOfStrains(); ++sid){
    bool takesid=false;
    for(uint32 rgid=1; rgid<ReadGroupLib::getNumReadGroups();++rgid){
      auto rg=ReadGroupLib::getReadGroupID(rgid);
      if(rg.getStrainID()==sid
	 && !rg.isBackbone()
	 && !rg.isRail()){
	takesid=true;
	break;
      }
    }
    if(takesid){
      con.newConsensusGet(straincons[sid],dummyqual,sid);
    }
  }

  // look if we have *any* consensus made (meaning we'd have same strain backbone AND mapped reads).
  //  If not, recalc consensi for all strains
  {
    bool needrecalc=true;
    for(auto & s : straincons){
      if(!s.empty()) {
	needrecalc=false;
	break;
      }
    }
    if(needrecalc){
      for(uint32 sid=0; sid<ReadGroupLib::getNumOfStrains(); ++sid){
	con.newConsensusGet(straincons[sid],dummyqual,sid);
      }
    }
  }

  {
    bool stillnocons=true;
    for(auto & s : straincons){
      if(!s.empty()) {
	stillnocons=false;
	break;
      }
    }
    BUGIFTHROW(stillnocons,"Ummmm ... no cons built???");
  }

  // fill up empty consensi
  for(auto & s : straincons){
    if(s.empty()) s.resize(con.getContigLength(),'*');
  }

  Contig::cccontainer_t & cc = const_cast<Contig::cccontainer_t &>(con.getConsensusCounts());

  auto ccI=cc.begin();
  for(uint32 actcontigpos=0; actcontigpos<straincons[0].size(); ++actcontigpos, ++ccI){
    bool maydelete=true;
    for(auto & s : straincons){
      if(s[actcontigpos]!='*'){
	maydelete=false;
	break;
      }
    }
    if(maydelete){
      ccI->total_cov=65535;
      ccI->star=65535;
    }
  }
  con.deleteStarOnlyColumns(0,con.getContigLength(),false,65535);
}
