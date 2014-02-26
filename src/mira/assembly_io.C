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


#include "util/fileanddisk.H"

#include "mira/assembly.H"
#include "mira/maf_parse.H"

#include "caf/caf.H"

#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>


using namespace std;


#define CEBUG(bla)
#define CEBUGF(bla)





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpContigs()
{
  cout << "The assembled project has " << AS_contigs.size() << " objects.\n";

  Contig::setCoutType(Contig::AS_TEXT);
  list<Contig>::iterator I=AS_contigs.begin();
  while(I!=AS_contigs.end()){
    cout << *I++ << "\n";
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

void Assembly::loadSequenceData()
{
  FUNCSTART("void Assembly::loadSequenceData()");

#ifdef TIMERESTRICTED
  if(AS_timesup) {
    cerr << "\n\nThis version of MIRA is old, please get a newer version of the assembler.\n";
    cerr << "\nCanonical download page: http://www.chevreux.org/mira_downloads.html\n";
    exit(0);
  }
#endif

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif

  discard();

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif

  if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
  cout << "\n";

  if(AS_resumeasembly){
    loadSequenceData_resume();
  }else{
    //loadSequenceData_new();
    loadSequenceDataFromManifest();
  }

#ifdef MIRAMEMORC
  MemORC::statistics();
#endif
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadSequenceData_resume()
{
  FUNCSTART("void Assembly::loadSequenceData_resume()");

//  cout << "Assembly::loadSequenceData_resume() must be adapted to readgroups\n"; exit(0);

  ReadGroupLib::discard();

  MAFParse mafp(&AS_readpool, nullptr, &AS_miraparams);
  vector<uint32> dummy;
  mafp.load(buildDefaultCheckpointFileName(AS_miraparams[0].getAssemblyParams().as_infile_chkptMAF),
	    ReadGroupLib::SEQTYPE_SANGER,
	    1,
	    dummy,
	    false,
	    nullptr
	  );

  bool templatesusable=AS_readpool.makeTemplateIDs();
  if(!templatesusable) {
    cout << "No useful template information found.\n";
  }

//  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
//
//  vector<uint32> lrperseqtype(ReadGroupLib::SEQTYPE_END,0);
//
//  size_t seqsloaded=loadMAF(buildDefaultCheckpointFileName(as_fixparams.as_infile_chkptMAF),
//			    ReadGroupLib::SEQTYPE_SANGER,
//			    1,                    // load directly
//			    lrperseqtype);
//
//  dumpSomeStatistics();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::loadSequenceDataFromManifest()
{
  FUNCSTART("void Assembly::loadSequenceDataFromManifest()");

  vector<Manifest::manifestloadentry_t> & manifestdata2load=AS_manifest.MAN_manifestdata2load;

  // quick consistency check for mapping
  {
    bool seenbb=false;
    bool seennotbb=false;
    for(const auto & mle : manifestdata2load){
      if(mle.loadasbackbone) {
	seenbb=true;
      }else{
	seennotbb=true;
      }
    }
    if(AS_miraparams[0].getAssemblyParams().as_assemblyjob_mapping){
      if(!seenbb){
	MIRANOTIFY(Notify::FATAL,"The \"job=...\" definition of the manifest says you want a mapping assembly, but no backbone sequence was given in the readgroups.\n\nDid you forget an 'is_reference' in one of the readgroups?");
      }
    }else{
      if(seenbb){
	MIRANOTIFY(Notify::FATAL,"The \"job=...\" definition of the manifest says you want a de-novo assembly, but there is a backbone/reference sequence given in the readgroups.\n\nThis is a slight contradiction, make up your mind please.");
      }
    }
  }

  // for delayed transformation reads to contigs
  // needed to simplify life when loading .fna and .gff3 with annotations
  // this gives the Readpool loader the opportunity to load the sequence,
  //  then the annotation, map the annotation to the sequence
  // then only make the contigs
  //
  // if this were not done this way, one would have to annotate both the readpool read
  //   AND search contigs for this read and annotate that one too
  // TODO: might be an idea for "re-annotation"
  vector<readid_t> readsasbackbonecontigs;

  // this loads the data completely (no contig or read callbacks given)
  streamSequenceDataFromManifest(AS_miraparams,
				 AS_manifest,
				 AS_readpool,
				 &AS_bbcontigs,
				 &readsasbackbonecontigs);

  checkForReadNameLength(AS_miraparams[0].getNagAndWarnParams().nw_check_mrnlvalue,
			 AS_miraparams[0].getNagAndWarnParams().nw_check_maxreadnamelength==NWSTOP);

  bool hassomeerror=false;
  {
    string logname(buildFileName(0,"","",
				 AS_miraparams[0].getAssemblyParams().as_tmpf_clippings,
				 ".txt","",false));
    string logprefix("load ancillary: ");

  for(const auto & mle : manifestdata2load){
    for(auto & fnfte : mle.ancillaryfilesfoundfordata){
      cout << "Merging ancillary data " << fnfte.fn << " type " << fnfte.ft << endl;
	if(fnfte.ft=="xml"){
	  AS_readpool.mergeXMLTraceInfo(fnfte.fn);
	}else if(fnfte.ft=="ssaha2"){
	  AS_readpool.mergeSSAHA2SMALTVecScreenData(fnfte.fn,
						    false,
						    logname,
						    logprefix);
	}else if(fnfte.ft=="smalt"){
	  AS_readpool.mergeSSAHA2SMALTVecScreenData(fnfte.fn,
						    true,
						    logname,
						    logprefix);
	}else{
	  BUGIFTHROW(true, "Unknown ancillary file type '" << fnfte.fn << "' for data '" << fnfte.ft << "'\n");
	}
      }
    }
  }

  if(hassomeerror){
    MIRANOTIFY(Notify::FATAL,"While looking at or loading ancillary files named in the manifest, some errors occured. Please check the log.");
  }

  // delayed creation of backbone contigs from single reads
  if(!readsasbackbonecontigs.empty()){
    Contig con(&AS_miraparams, AS_readpool);
    for(auto & rpi : readsasbackbonecontigs){
      AS_bbcontigs.push_back(con);
      AS_bbcontigs.back().addFirstRead(rpi,1);
    }
  }

  if(!AS_bbcontigs.empty()){
    postLoadBackbone();
  }

  // special handling for 454 data (20.09.2008)
  //  as there's no separate seqvec clip in the XML files (until now)
  //  and MIRA now uses a "clip back, extend later" strategy for 454
  //  reads, the 454 adaptor must be protected from extension as it
  //  happens often enough that two reads start the sequencing adaptor
  //  at the same time ... WHICH WOULD THEN BE UNCOVERED!
  //
  // E.g.
  //  ACCGTCAGTCAGTCAGTGTTGACGTGTCAccctgagacacgcaacaggggatagacaaggca
  //  ACCGTACGTCAG*CAGTGTTGACGTGTCAccctgagacacgcaacaggggatagacaaggca
  //
  // Two possible worarounds
  // 1) instruct extendADS() not to extend into lower case (bad, relies
  //    on case information)
  // 2) transform the right qual clip into a vec clip if there is no
  //    vec clip
  //
  // we'll do number 2 here

  // find out what we have in the pool
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 rgi=0; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(!rgid.isRail() && !rgid.isBackbone()){
      AS_seqtypespresent[rgid.getSequencingType()]=true;
    }
  }

#if CPP_READ_SEQTYPE_END != 8
#error "Check if new seqtype needs same workaround."
#endif
  if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_454GS20]){
    uint32 changecount=0;
    for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
      if(!AS_readpool.getRead(rnr).isBackbone() && !AS_readpool.getRead(rnr).isRail()){
	if(AS_readpool.getRead(rnr).hasValidData()){
	  // if no right seq vec but a right clip
	  if(AS_readpool.getRead(rnr).getRSClipoff() == static_cast<int32>(AS_readpool.getRead(rnr).getLenSeq())
	     && AS_readpool.getRead(rnr).getRQClipoff() != static_cast<int32>(AS_readpool.getRead(rnr).getLenSeq())){
	    // make right seq vec = right clip
	    AS_readpool.getRead(rnr).setRSClipoff(AS_readpool.getRead(rnr).getRQClipoff());
	    changecount++;
	  }
	}
      }
    }
    if(changecount){
      cout << "Note: " << changecount << " reads with 454 data had quality clips given, but no sequencing vector clip.\n"
	   << "For MIRA to run properly with read extension, those quality clips have been\n"
	   << "changed to sequencing vector clips.\n\n";
    }
  }

  postLoad();
  cout << "Have read pool with " << AS_readpool.size() << " reads.\n";

  dumpSomeStatistics();

  basicDataChecks();
  basicReadGroupChecks();

  clipsAfterLoad();

  dumpSomeStatistics();

  vector<uint32> dummy;
  if(AS_bbcontigs.empty()) {
    sortReadPool(dummy);
  }else{
    addRailsToBackbones();
    sortReadPool(dummy);

    list<Contig>::iterator bbcI=AS_bbcontigs.begin();
    for(; bbcI!=AS_bbcontigs.end(); ++bbcI){
      bbcI->exchangeReadIDs(dummy);
      // a few things to be redone
      //bbcI->setupAsBackBoneContig();
    }
    dumpSomeStatistics();
  }

  ofstream fout;
  fout.open((AS_miraparams[0].getDirectoryParams().dir_tmp+'/'+AS_miraparams[0].getAssemblyParams().as_tmpf_poolinfo+".lst").c_str(), ios::out| ios::trunc);
  AS_readpool.dumpPoolInfo(fout);
  fout.close();

  FUNCEND();
  return;
}


void Assembly::sortReadPool(vector<uint32> & dummy)
{
  FUNCSTART("void Assembly::sortReadPool()");

  AS_readpool.checkTemplateIDs("Before sorting: template partners have template mismatch!");

  dummy.clear();
  cout << "Sorting reads ... "; cout.flush();
  AS_readpool.sortPoolToMIRAStandard(dummy);
  cout << "done."<<endl;

  if(AS_debrisreason.size()){
    BUGIFTHROW(AS_debrisreason.size() != dummy.size(),"AS_debrisreason.size() != dummy.size()");
    vector<uint8> tmp;
    tmp.resize(AS_debrisreason.size());
    auto dstI=tmp.begin();
    auto srcI=dummy.begin();
    for(; dstI!=tmp.end(); ++srcI, ++dstI){
      *dstI=AS_debrisreason[*srcI];
    }
    AS_debrisreason.swap(tmp);
  }

  // The reads have Template Partner IDs internally which are now wrong
  // Get them corrected
  {
    vector<uint32> backmap(dummy.size(),-1);
    for(uint32 ri=0;ri<backmap.size(); ++ri){
      backmap[dummy[ri]]=ri;
    }

    for(uint32 ri=0;ri<AS_readpool.size(); ++ri){
      CEBUG(cout << AS_readpool[ri].getName() << "\t" << AS_readpool[ri].getTemplatePartnerID());
      if(AS_readpool[ri].getTemplatePartnerID()>=0){
	CEBUG("\t" << AS_readpool[AS_readpool[ri].getTemplatePartnerID()].getName());
	AS_readpool[ri].setTemplatePartnerID(backmap[AS_readpool[ri].getTemplatePartnerID()]);
	CEBUG("\t" << AS_readpool[AS_readpool[ri].getTemplatePartnerID()].getName());
      }else{
	CEBUG("\t-");
      }
      CEBUG("\n");
    }

    AS_readpool.checkTemplateIDs("After sorting: template partners have template mismatch!");
  }
}

/*************************************************************************
 *
 * All pointer parameters may be nullptr
 * - if contigsptr is nullptr, then the CAF and MAF loading routines will
 *   know they do not need to build contig structures, saves time and
 *   memory
 * - readsasbackbonecontigs will contain a list of readpool ids of those
 *   *single reads!* which were not embedded in a contig but should be
 *   seen as backbone sequence (like e.g., loaded from FASTA, GFF3, FASTQ,
 *   etc.pp
 * - if the callbacks for contigs and/or reads are nullptr, then obviously
 *   no callback will occur.
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::streamSequenceDataFromManifest(std::vector<MIRAParameters> & miraparams, Manifest & man, ReadPool & readpool, list<Contig> * contigsptr, vector<readid_t> * readsasbackbonecontigs, void (*ccallback)(list<Contig> &, ReadPool &), void (*rcallback)(ReadPool &))

{
  FUNCSTART("streamSequenceDataFromManifest()");

  auto & manifestdata2load=man.MAN_manifestdata2load;

  bool hassomeerror=false;
  for(const auto & mle : manifestdata2load){
    for(auto & fnfte : mle.mainfilesfoundfordata){
      string fn2;

      if(fnfte.ft=="fasta"){
	fn2=fnfte.fn+".qual";
      }

      size_t oldrpsize=readpool.size();
      size_t oldbbsize=0;
      if(contigsptr!=nullptr) oldbbsize=contigsptr->size();

      if((mle.loadasbackbone || contigsptr!=nullptr)
	 && (fnfte.ft=="maf" || fnfte.ft=="caf")){
	// Note: this branch really only if we want contigs or load as backbone (may be
	//  two different usages from the caller function)
	// MAF & CAF without need for this will still load, using the readpool.loadData_rgid()
	//  branch further below
	cout << "Loading reads or assembled contigs ";
	if(mle.loadasbackbone) cout << "as reference backbone ";
	cout << "from " << fnfte.fn << " type " << fnfte.ft << endl;
	if(fnfte.ft=="caf") {
	  CAF tcaf(&readpool, contigsptr, &miraparams);
	  vector<uint32> dummy;
	  tcaf.load(fnfte.fn.c_str(),
		    mle.rgid.getSequencingType(),
		    1,
		    dummy,
		    false,
		    ccallback,
		    rcallback);
	}else if(fnfte.ft=="maf") {
	  MAFParse mafp(&readpool, contigsptr, &miraparams);
	  vector<uint32> dummy;
	  mafp.load(fnfte.fn.c_str(),
		    mle.rgid.getSequencingType(),
		    1,
		    dummy,
		    false,
		    ccallback,
		    rcallback);
	}else{
	  BUGIFTHROW(true,"Should never arrive here???");
	}
      }else{
	if(mle.loadasbackbone){
	  cout << "Loading reference backbone from " << fnfte.fn << " type " << fnfte.ft << endl;
	}else{
	  cout << "Loading reads from " << fnfte.fn << " type " << fnfte.ft << endl;
	}
	readpool.loadData_rgid(fnfte.ft,fnfte.fn,fn2,mle.rgid,
				   false,
				   rcallback);
      }

      // should have already been set by the manifest parser, just in case
      BUGIFTHROW(mle.loadasbackbone && !mle.rgid.isBackbone(),"manifest entry says to load as backbone, but readgroup is not backbone?");

      if(contigsptr!=nullptr && contigsptr->size() != oldbbsize){
	cout << "contained " << contigsptr->size()-oldbbsize << " contigs. Only the contigs will be added as backbone.\n";
      }else if(mle.loadasbackbone){
	// there were no contigs in this file: set qualities (if needed) and mark the reads to be added as contigs later
	for(size_t rpi=oldrpsize; rpi<readpool.size(); ++rpi) {
	  if(readpool[rpi].hasValidData()){
	    if(!readpool[rpi].hasQuality()){
	      readpool[rpi].setQualities(readpool[rpi].getReadGroupID().getDefaultQual());
	      readpool[rpi].setQualityFlag(false);
	    }
	    if(readsasbackbonecontigs!=nullptr) readsasbackbonecontigs->push_back(static_cast<readid_t>(rpi));
	  }
	}
      }
    }
  }

  if(hassomeerror){
    MIRANOTIFY(Notify::FATAL,"While looking at or loading data files named in the manifest, some errors occured. Please check the log.");
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::addRailsToBackbones()
{
  FUNCSTART("void Assembly::addRailsToBackbones()");
  vector<uint32> lrperseqtype(ReadGroupLib::SEQTYPE_END,0);
  uint32 longestread=0;

  for(uint32 i=0; i<AS_readpool.size(); ++i){
    if(AS_readpool[i].isRail()
       || AS_readpool[i].isBackbone()
       || AS_readpool[i].isCoverageEquivalentRead()) continue;
    if(AS_readpool[i].getLenClippedSeq()>longestread) longestread=AS_readpool[i].getLenClippedSeq();
    if(AS_readpool[i].getLenClippedSeq()>lrperseqtype[AS_readpool[i].getSequencingType()]) lrperseqtype[AS_readpool[i].getSequencingType()]=AS_readpool[i].getLenClippedSeq();
  }

  if(longestread==0){
    MIRANOTIFY(Notify::FATAL, "No read with sequence length >0 present? Did you provide data to load?");
  }

  // if wanted, determine solexahack parameter automatically
  if(AS_miraparams[0].getAlignParams().al_solexahack_maxerrors<0
     && lrperseqtype[ReadGroupLib::SEQTYPE_SOLEXA]){
    cout << "-AL:shme is < 0, automatically determining optimal value.\n";
    AS_miraparams[0].getNonConstAlignParams().al_solexahack_maxerrors=
      lrperseqtype[ReadGroupLib::SEQTYPE_SOLEXA]*15/100;
    cout << "set -AL:shme "
	 << AS_miraparams[0].getAlignParams().al_solexahack_maxerrors
	 << '\n';
  }

  // if wanted, determine rail length and overlap automatically
  if(AS_miraparams[0].getAssemblyParams().as_backbone_raillength == 0){
    cout << "-SB:brl is 0, automatically determining optimal value.\n";

    // add 15% to longest read (so accomodate insertion), then times 2
    uint32 newraillength=(longestread*115/100) * 2;
    if(newraillength > 32760){
      cout << "Optimal rail would be longer than 32k, adjusting down to 32k.\n";
      newraillength=32760;
    }
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength=newraillength;
    cout << "brl: "
	 << AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength
	 << '\n';
  }
  if(AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap == 0){
    cout << "-SB:bro is 0, automatically determining optimal value.\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap=
      AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength/2;
    cout << "bro: "
	 << AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap
	 << '\n';
  }
  if(AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap >=
     AS_miraparams[0].getAssemblyParams().as_backbone_raillength){
    cout << "-SB:bro is >= -SB:brl ... adjusting -SB:bro to (-SB:brl)-1\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_backbone_railoverlap=
      AS_miraparams[0].getNonConstAssemblyParams().as_backbone_raillength-1;
  }

  // TODO: see whether rails need strains!!!
  // Answer: yes, it helps. E.g.: markFeaturesByConsensus() and other routines do not
  //  need to work around a strain which says "I have 3 strains" (backbone, rails (==empty) and
  //  reads).
  // Set the strain name to be equal to the strain name of the first backbone read encountered

  size_t numrailscreated=0;

  for(auto & bbc : AS_bbcontigs){
    bbc.recalcTemplateIDsAndStrainPresent();

    string strainname;
    bool foundbbread=false;
    for(auto & cr : bbc.getContigReads()){
      if(cr.isBackbone()){
	strainname=cr.getStrainName();
	foundbbread=true;
	break;
      }
    }
    BUGIFTHROW(foundbbread==false,"no backbone read found in backbone???");

    bool bbvalue=true;

    // add the rails
    if(bbvalue) {
      numrailscreated+=bbc.addRails(
	AS_miraparams[0].getAssemblyParams().as_backbone_raillength,
	AS_miraparams[0].getAssemblyParams().as_backbone_railoverlap,
	strainname,
	AS_miraparams[0].getAssemblyParams().as_backbone_strainname_forceforall,
	AS_miraparams[0].getAssemblyParams().as_backbone_rail_fromstrain,
	false);
    }
  }

  AS_debrisreason.resize(AS_readpool.size(),DEBRIS_NOTDEBRIS);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::clipsAfterLoad()
{
  AS_debrisreason.clear();
  AS_debrisreason.resize(AS_readpool.size(),DEBRIS_NOTDEBRIS);

  vector<unique_ptr<DataProcessing>> dpv(AS_miraparams[0].getAssemblyParams().as_numthreads);

  for(uint32 ti=0; ti<dpv.size(); ++ti){
    dpv[ti]=unique_ptr<DataProcessing>(new DataProcessing(&AS_miraparams));
    string logname(buildFileName(0,"","",
				 AS_miraparams[0].getAssemblyParams().as_tmpf_clippings
				 + "_t" + boost::lexical_cast<string>(ti),
				 ".txt","",false));
    cout << logname << endl;
    dpv[ti]->startLogging(logname,false);
  }

  string logprefix("loadclip: ");
  cout << "Post-load clips:\n";

  DataProcessing::stdTreatmentPool_MultiThread(AS_dataprocessing,dpv,AS_readpool,&AS_debrisreason,logprefix,true);

  cout << endl;
  if(AS_dataprocessing.DP_stats.cphix174){
    cout << "SEARCH MSG: PhiX 174 found: " << AS_dataprocessing.DP_stats.cphix174 << endl;
  }
  if(AS_dataprocessing.DP_stats.cadapright){
    cout << "CLIP MSG: Adaptor right found: " << AS_dataprocessing.DP_stats.cadapright << endl;
  }
  if(AS_dataprocessing.DP_stats.cadaprightpartial){
    cout << "CLIP MSG: Partial adaptor right found: " << AS_dataprocessing.DP_stats.cadaprightpartial << endl;
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpRailReads(ofstream & fout)
{
  Read::setCoutType(Read::AS_FASTA);
  for(uint32 i=0; i<AS_readpool.size(); ++i){
    if(AS_readpool[i].isRail()){
      fout << AS_readpool[i].isRail();
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::basicDataChecks()
{
  FUNCSTART("void Assembly::basicDataChecks()");

  uint32 nummsglr=0;
  uint32 nummsgqual=0;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    if(AS_readpool[ri].getLenClippedSeq() > 29990
       && !AS_readpool[ri].isBackbone()){
      if(++nummsglr<=10){
	cout << "Read " << AS_readpool[ri].getName() << " has " << AS_readpool[ri].getLenClippedSeq() << " bases (clipped). Too long (>29990)\n";
	if(nummsglr==10) cout << "More long reads may exist, but stopping output here.\n";
      }
    }
    bool bahqual=false;
    for(auto q : AS_readpool[ri].getQualities()){
      if(q>=100){
	bahqual=true;
	break;
      }
    }
    if(bahqual){
      if(++nummsgqual<=10){
	cout << "Read " << AS_readpool[ri].getName() << " has quality values >100, this is illegal.\n";
	if(nummsgqual==10) cout << "More reads with bad quals may exist, but stopping output here.\n";
      }
    }
  }

  if(nummsglr){
    MIRANOTIFY(Notify::FATAL,"MIRA found " << nummsglr << " very long reads (too long for normal reads), see log above.");
  }
  if(nummsgqual){
    MIRANOTIFY(Notify::FATAL,"MIRA found " << nummsgqual << " reads with illegal qualities, see log above.");
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::basicReadGroupChecks()
{
  FUNCSTART("void Assembly::basicReadGroupChecks()");

  bool haserror_ep=false; // expected pair
  bool haserror_nsp=false; // no segment placement
  for(uint32 rglid=1; rglid<ReadGroupLib::getNumReadGroups(); ++rglid){
    auto rgid = ReadGroupLib::getReadGroupID(rglid);
    cout << "Checking pairs of readgroup " << rglid << " (named: '" << rgid.getGroupName() << "'): ";
    cout.flush();
    uint32 numreadswithtpids=0;
    uint32 numreads=0;
    for(uint32 rid=0; rid<AS_readpool.size(); ++rid){
      if(AS_readpool[rid].getReadGroupID()==rgid){
	++numreads;
	if(AS_readpool[rid].getTemplatePartnerID()>=0) ++numreadswithtpids;
      }
    }
    cout << " found " << numreadswithtpids<<endl;
    if(rgid.expectsReadPairs()){
      if(numreads && numreadswithtpids==0) {
	cout << "WARNING: in the above readgroup, no read is paired although the manifest says there should be pairs. This is fishy!\n";
	haserror_ep=true;
      }
    }else{
      if(numreadswithtpids>0) {
	string wstr("In readgroup ");
	wstr+=boost::lexical_cast<string>(rgid.getLibId())+" (named: '"+rgid.getGroupName()+"') paired reads were found but no pairing information given in the manifest. MIRA will estimate 'template_size' and 'segment_placement'.\nYou can suppress this warning by using the keyword 'autopairing' in the readgroup definition of the manifest file.";
	AS_warnings.setWarning("READGROUP_UNEXPECTED_PAIRS",2,"Unexpected pairs in readgroup",wstr);
      }
    }
  }

  if(haserror_ep){
    MIRANOTIFY(Notify::FATAL,"MIRA found readgroups where pairs are expected but no read has a partner. See log above and then check your input please (either manifest file or data files loaded or segment_naming scheme).");
  }

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void Assembly::checkForReadNameLength(uint32 stoplength, bool stop)
{
  FUNCSTART("void Assembly::checkForReadNameLength(uint32 stoplength)");

  if(stoplength==0) return;

  // if names are too long, print out only the first 20 of each read group
  vector<uint32> countperrg(ReadGroupLib::getNumReadGroups(),20);

  uint32 count=0;
  bool maybemore=false;
  for(uint32 ri=0; ri<AS_readpool.size(); ++ri){
    if(AS_readpool[ri].getName().size()>stoplength){
      if(countperrg[AS_readpool[ri].getReadGroupID().getLibId()]){
	if(--countperrg[AS_readpool[ri].getReadGroupID().getLibId()] == 0) maybemore=true;
	if(count==0) {
	  cout << "List of read names which have problems with name length:\n";
	}
	cout << "Name too long: " << AS_readpool[ri].getName() << '\n';
      }
      ++count;
    }
  }
  if(maybemore){
    cout << count << " reads had a long name length, for brevity's sake not all were listed.\n";
  }
  if(count>0){
    string emsg=boost::lexical_cast<string>(count)+
      " reads were detected with names longer than "
      +boost::lexical_cast<string>(stoplength)+
      " characters.\n\n"
      "While MIRA and many other programs have no problem with that, some older "
      "programs have restrictions concerning the length of the read name.\n"
      "\nExample given: the pipeline\n"
      "     CAF -> caf2gap -> gap2caf\n"
      "will stop working at the gap2caf stage if there are read names having > 40 characters "
      "where the names differ only at >40 characters.\n"
      "\nThis is a warning only, but as a couple of people were bitten by this, the default "
      "behaviour of MIRA is to stop when it sees that potential problem.\n"
      "\nYou might want to rename your reads to have <= "
      +boost::lexical_cast<string>(stoplength) +
      " characters.\n"
      "\nOn the other hand, you also can ignore this potential problem and force MIRA to "
      "continue by using the parameter: '-NW:cmrnl=warn'  or  '-NW:cmrnl=no'\n";
    if(stop){
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << endl;
    }
  }
  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::postLoadBackbone()
{
  FUNCSTART("void Assembly::postLoadBackbone()");

  AS_hasbackbones=true;

  cout << "Deleting gap columns in backbones ... "; cout.flush();
  for(auto I=AS_bbcontigs.begin(); I!=AS_bbcontigs.end(); I++){
    I->deleteStarOnlyColumns(0,I->getContigLength());
  }

  cout << "Postprocessing backbone(s) ... this may take a while."<< endl;

  // mark all reads loaded in backbone as backbone
  // check that they are not named "ContigX"
  // set the strain to "backbone"
  // Backbones will not be included is Skim, makeAlignments etc.


  static const boost::regex badseqnameexp("^Contig[0-9]+$");
  //return regex_match(s, e);

  // set MFSM tags
  if(0){
    vector<multitag_t::mte_id_t> idstoreplace;
    {
      string tmp="FLTR";
      idstoreplace.push_back(Read::getTagID(tmp));
      tmp="FrRN";
      idstoreplace.push_back(Read::getTagID(tmp));
    }
    string mfsm="MFSM";
    multitag_t::mte_id_t mtid=Read::getTagID(mfsm);
    for(auto & ce : AS_bbcontigs){
      auto & conreads=ce.getContigReads();
      for(auto & cre : conreads){
	for(uint32 it=0; it<cre.getNumOfTags(); ++it){
	  multitag_t tmp=cre.getTag(it);
	  if(find(idstoreplace.begin(),idstoreplace.end(),tmp.identifier)!=idstoreplace.end()){
	    tmp.identifier=mtid;
	    tmp.source=multitag_t::MT_tagsrcentry_idMIRA;
	    Read & ncr=const_cast<Read &>(cre);
	    ncr.addTagO(tmp);
	  }
	}
      }
    }
  }


  bool contignamesok=true;
  {
    uint32 bbnum=0;
    cout << AS_bbcontigs.size() << " to process\n";
    for(auto & ce : AS_bbcontigs){
      ++bbnum;
      // first, find a name for that contig

      auto & conreads=ce.getContigReads();

      // if it is a single read contig,
      //  set the name for that contig to be the name of the read
      if(conreads.size()==1) {
	if(conreads.begin()->getName().size()){
	  if(conreads.begin()->getName()[0]=='C'){
// BOOST: regex not compatible with _GLIBCXX_DEBUG
#ifndef _GLIBCXX_DEBUG
	    if(regex_match(conreads.begin()->getName(), badseqnameexp)){
	      cout << "Bad name for backbone sequence " << bbnum << ": " << conreads.begin()->getName() << '\n';
	      cout << "Backbone sequences may NOT be name 'ContigX' with 'X' being any number.\n";
	      contignamesok=false;
	    }
#endif
	  }
	}else{
	  cout << "There's a backbone sequence (number " << bbnum << ") without a name? Too bad, not allowed.\n";
	  contignamesok=false;
	}

	/* BaCh 14.09.2010: why did I add _bb? Let's remove and see whether it's better.
	   BaCh 26.10.2010: I remembered. A contig may not be named like a sequence, CAF will
	                    dump an error.
			    Solution: ...?
	*/
	ce.setContigName(conreads.begin()->getName()+"_bb");

      }

      cout << ce.getContigName() << "\t" << ce.getContigLength() << endl;

      //bool bbvalue=true;

      ////  except singlets?!
      //if(conreads.size()==1 && I->getContigLength()<4000) {
      //	bbvalue=false;
      //}

      // now let the contig do the rest of the setup
      ce.setupAsBackBoneContig();
    }
  }

  if(!contignamesok){
    MIRANOTIFY(Notify::FATAL,"Some backbones had either no names or a bad name (see log above). Stopping here, fix your sequence names.\n")
  }

  ReadGroupLib::dumpStrainIDSummary();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::postLoad()
{
  FUNCSTART("void Assembly::postLoad()");

  if(!AS_readpool.checkForDuplicateReadNames()){
    MIRANOTIFY(Notify::FATAL,"MIRA found duplicate read names in your data (see log above for more info).\n\nThis should never, never be!\n\nYou may have loaded a file more than once in the manifest or\nreads may be present multiple times across your input file(s).\nEither way: fix that!\n");
  }

  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  // how many sequences in assembly


  // count how many have valid data, clean up stars in reads
  // count how many have SCF data: if none, switch off editing.

  {
    cout << "Checking reads for trace data (loading qualities if needed):\n";

    ProgressIndicator<int32> P (0, AS_readpool.size());

    ofstream fout;
    fout.open((dir_params.dir_tmp+'/'+as_fixparams.as_outfile_stats_reads_invalid).c_str(), ios::out| ios::trunc);

    bool can_edit=false;
    AS_num_reads_valid=0;
    for(uint32 i=0;i<AS_readpool.size();i++){
      P.progress(i);
      if(AS_readpool.getRead(i).hasValidData()){
	AS_num_reads_valid++;
	if(AS_readpool[i].isBackbone()){
	  AS_hasbackbones=true;
	}
	AS_readpool.getRead(i).removeGapsFromRead();
	if(AS_readpool.getRead(i).hasSCFData(true)){
	  can_edit=true;
	}
      } else {
	if(!AS_readpool.getRead(i).getName().empty()) {
	  fout << AS_readpool.getRead(i).getName() << "\n";
	//} else if (!AS_readpool.getRead(i).getEXPName().empty()) {
	//  fout << AS_readpool.getRead(i).getEXPName() << "\n";
	} else {
	  fout << "Unknown read (loaded as number: " << i << ")\n";
	}
      }
    }
    P.finishAtOnce();
    cout << endl;

    if(!can_edit) {
      cout << "No SCF data present in any read, EdIt automatic contig editing for Sanger data is now switched off.\n";
      const_cast<edit_parameters &>(AS_miraparams[0].getEditParams()).ed_edit_automatic_contic_editing=false;
    }
    fout.close();
  }

  cout << AS_num_reads_valid << " reads with valid data for assembly.\n";

  if(!AS_num_reads_valid){
    MIRANOTIFY(Notify::FATAL, "No valid read in assembly?");
  }


  bool templatesusable=AS_readpool.makeTemplateIDs();
  if(!templatesusable) {
    cout << "No useful template information found.\n";
  }


  //re-adjust bbcontigs template and strain ids in contig reads
  if(AS_hasbackbones){
    list<Contig>::iterator I=AS_bbcontigs.begin();
    for(; I!=AS_bbcontigs.end(); I++){
      I->recalcTemplateIDsAndStrainPresent();
    }
  }

  ReadGroupLib::dumpStrainIDSummary();

  // look for quality values in reads
  {
    bool stopall=false;
    for(uint32 i=0;i<AS_readpool.size();i++){
      Read & actread=AS_readpool.getRead(i);
      if(actread.isBackbone()
	 || actread.isRail()) continue;
      if(AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_enforce_qualsinreads
	 && actread.hasValidData()
	 && actread.hasQuality()==false
	 && !actread.hasUserDefaultQuality()){
	cout << "No quality data found: (" << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType()) << ") " << actread.getName() << '\n';
	stopall=true;
      }
    }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
    if(stopall) {
      MIRANOTIFY(Notify::FATAL,"Some reads had no quality values given (see log above),\nplease check your input data.\nIf sure that this is ok for your data, switch off this check with -AS:epoq=no\nfor any sequencing type you wish (Sanger, 454, IonTorrent, PacbioLQ, PacbioHQ, Text, Solexa, ...).\nAlso consider the '--noqualities' parameter setting.\nAlternatively, you can force to switch off this check for specific readgroups by using the 'default_qual' setting in the manifest.")
    }

  }


  FUNCEND();
  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpSomeStatistics()
{
  FUNCSTART("void Assembly::dumpSomeStatistics()");
  // initialise the assembly_structure and do some statistics

  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  uint32 backbonereads=0;
  uint32 railreads=0;
  {
    ofstream fout;
    fout.open((dir_params.dir_tmp+'/'+as_fixparams.as_outfile_stats_reads_tooshort).c_str(), ios::out| ios::trunc);

    AS_num_reads_too_small=0;
    for(uint32 i=0;i<AS_readpool.size();i++){
      Read & actread=AS_readpool.getRead(i);

      //AS_ok_for_assembly[i]=1;
      actread.setUsedInAssembly(true);

      if(actread.hasValidData()==false){
	cout << actread.getName() << ": unable to load or other reason for invalid data.\n";
	//AS_ok_for_assembly[i]=0;
	actread.setUsedInAssembly(false);
      }else{
	if(actread.isBackbone()){
	  actread.setUsedInAssembly(false);
	  backbonereads++;
	} else if(actread.isRail()) {
	  railreads++;
	}else{
	  // throw out on minumum length if no template partner is present
	  if(actread.getLenClippedSeq() < AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength){
	    //cout << "Short length: "
	    //	 << actread.getName() << " ("
	    //	 << actread.getShortNameOfSequencingType(actread.getSequencingType())
	    //	 << "): only " << actread.getLenClippedSeq()
	    //	 << " good bases, ";
	    fout << actread.getName();
	    if(actread.getTemplatePartnerID() == -1){
	      //cout << "need: "
	      //	   << AS_miraparams[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength
	      //	   << ". No paired end partner, rejected.\n";
	      fout << " too small and no paired end\n";
	      AS_num_reads_too_small++;
	      //AS_ok_for_assembly[i]=0;
	      actread.setUsedInAssembly(false);
	    }else{
	      if(actread.getLenClippedSeq() < 20){
		//cout << "really too small, rejected.\n";
		fout << " too small even with paired end\n";
		AS_num_reads_too_small++;
		//AS_ok_for_assembly[i]=0;
		actread.setUsedInAssembly(false);
	      }else{
		fout << " saved by paired-end\n";
		//cout << "accepted as paired-end partner is present.\n";
	      }
	    }
	  }
	}
      }
    }
    fout.close();
  }

  if(AS_logflag_dumpusedids){
    ofstream fout;
    fout.open((dir_params.dir_tmp+"/elog.usedids.lst").c_str(), ios::out | ios::trunc);
    for(uint32 i=0; i<AS_used_ids.size(); i++){
      if(AS_readpool[i].isUsedInAssembly()){
	fout << AS_readpool[i].getName()<< '\n';
      }
    }
    fout.close();
  }

  // TODO: also take reads too short into statistics
  //

  vector<int32> withoutclipsbytype(ReadGroupLib::SEQTYPE_END,0);
  vector<int32> readcountbytype(ReadGroupLib::SEQTYPE_END,0);
  vector<int32> readusedbytype(ReadGroupLib::SEQTYPE_END,0);
  vector<int32> numnoqualbytype(ReadGroupLib::SEQTYPE_END,0);

  vector<uint64_t> meanlengthbytype(ReadGroupLib::SEQTYPE_END,0);
  vector<uint64_t> meantlengthbytype(ReadGroupLib::SEQTYPE_END,0);

  //remove("log.noqualities");
  ofstream fout;
  fout.open((dir_params.dir_tmp+"/miralog.noqualities").c_str(), ios::out| ios::trunc);
  for(uint32 i=0;i<AS_readpool.size();i++){
    Read & actread=AS_readpool.getRead(i);
    if(actread.isBackbone() == false
       && actread.isRail() == false){
      if(actread.hasQuality()==false){
	numnoqualbytype[actread.getSequencingType()]++;
	fout << actread.getName() << endl;
      }
      readcountbytype[actread.getSequencingType()]++;
      meantlengthbytype[actread.getSequencingType()]+=actread.getLenSeq();
      if(actread.isUsedInAssembly()){
	readusedbytype[actread.getSequencingType()]++;
	meanlengthbytype[actread.getSequencingType()]+=actread.getLenClippedSeq();
      }
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
      switch(actread.getSequencingType()){
      case ReadGroupLib::SEQTYPE_TEXT :
      case ReadGroupLib::SEQTYPE_SANGER :
      case ReadGroupLib::SEQTYPE_PACBIOHQ :
      case ReadGroupLib::SEQTYPE_PACBIOLQ :
      case ReadGroupLib::SEQTYPE_IONTORRENT :
      case ReadGroupLib::SEQTYPE_454GS20 :
      case ReadGroupLib::SEQTYPE_SOLEXA : {
	if(actread.getLenClippedSeq() == actread.getLenSeq()){
	  withoutclipsbytype[actread.getSequencingType()]++;
	}
	break;
      }
      case ReadGroupLib::SEQTYPE_ABISOLID : {
	MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 20a.");
	break;
      }
      default : {
	cerr << "Sequencing type " << actread.getSequencingType() << " unknown?\n";
	MIRANOTIFY(Notify::FATAL, "Found unknown sequencing type in read.");
      }
      }
    }
  }
  fout.close();

  vector<uint64_t> totallengthbytype=meanlengthbytype;
  vector<uint64_t> totaltlengthbytype=meantlengthbytype;

  for(uint32 i=0; i<meantlengthbytype.size(); i++){
    if(readcountbytype[i]) meantlengthbytype[i]/=readcountbytype[i];
    if(readusedbytype[i]) meanlengthbytype[i]/=readusedbytype[i];
  }

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
  cout << "\n===========================================================================\n";
  cout << "Pool statistics:\n";
  cout << "Backbones: " << backbonereads << "\tBackbone rails: " << railreads << "\n";

  cout << "\n\t\tSanger\t454\tIonTor\tPcBioHQ\tPcBioLQ\tText\tSolexa\tSOLiD\n";
  cout << "\t\t------------------------------------------------------------\n";
  cout << "Total reads\t" << readcountbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << readcountbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';
  cout << "Reads wo qual\t" << numnoqualbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << numnoqualbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';
  cout << "Used reads\t" << readusedbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << readusedbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';
  cout << "Avg tot rlen\t" << meantlengthbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << meantlengthbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';
  cout << "Avg rlen used\t" << meanlengthbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << meanlengthbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';
  cout << "W/o clips\t" << withoutclipsbytype[ReadGroupLib::SEQTYPE_SANGER];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_454GS20];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_IONTORRENT];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_PACBIOHQ];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_PACBIOLQ];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_TEXT];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_SOLEXA];
  cout << '\t' << withoutclipsbytype[ReadGroupLib::SEQTYPE_ABISOLID] << '\n';

  cout << '\n';
  for(uint32 i=0; i<totallengthbytype.size(); i++){
    if(totaltlengthbytype[i]>0){
      cout << ReadGroupLib::getNameOfSequencingType(i) << "\ttotal bases: " << totaltlengthbytype[i]
	   << "\tused bases in used reads: " << totallengthbytype[i] << '\n';
    }
  }

  cout << "===========================================================================\n\n";


  if(AS_readpool.size()-AS_num_reads_too_small-backbonereads-railreads<=0) {
    MIRANOTIFY(Notify::FATAL, "No read can be used for assembly.");
  }

  cout << endl;

  FUNCEND();
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveExtTmpContig(Contig & con, string basename)
{
  FUNCSTART("void Assembly::saveExtTmpContig(string prepost)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  //directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();

  if(con.getNumReadsInContig() > 1
     ||   as_fixparams.as_output_exttmp_alsosinglets){

    if (as_fixparams.as_output_exttmp_caf) {
      string filename=basename+".caf";

      cout << "Logging this contig to file: " << filename << endl;

      ofstream cafout(filename.c_str(), ios::out | ios::trunc);
      Contig::setCoutType(Contig::AS_CAF);
      cafout << con;
      cafout.close();
    }

    if (as_fixparams.as_output_exttmp_ace) {
      string filename=basename+".ace";
      cout << "Logging this contig to file: " << filename << endl;

      ofstream aceout(filename.c_str(), ios::out | ios::trunc);
      Contig::setCoutType(Contig::AS_ACE);
      aceout << con;
      aceout.close();
    }


    if (as_fixparams.as_output_exttmp_fasta) {
      string filename=basename+".fasta";
      string qualname=filename+".qual";

      cout << "Logging this contig to files: " << filename << "  and  " << qualname << endl;

      ofstream fastaout(filename.c_str(), ios::out | ios::trunc);
      Contig::setCoutType(Contig::AS_FASTA);
      fastaout << con;
      fastaout.close();
      ofstream qualout(qualname.c_str(), ios::out | ios::trunc);
      Contig::setCoutType(Contig::AS_FASTAQUAL);
      qualout << con;
      qualout.close();
    }


    if (as_fixparams.as_output_exttmp_gap4da) {
      string dirname=basename+".gap4da";

      cout << "Logging this contig to directory: " << dirname << endl;
      if(ensureDirectory(dirname,true)){
	MIRANOTIFY(Notify::FATAL, "Cannot make sure the directory exist? Aborting.");
      }

      Contig::setCoutType(Contig::AS_GAP4DA);
      ofstream fofnout((dirname+"/fofn").c_str(), ios::out | ios::trunc);
      con.saveAsGAP4DA(dirname, fofnout);
      fofnout.close();
    }
  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::buildDefaultCheckpointFileName(const string & filename)
{
  return AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/"+filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::buildDefaultInfoFileName(int32 version, const string & prefix, const string & postfix, const string & basename, const string & defaultname, const string & defaultextension, bool removeold)
{
  string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_info;
  }

  string filename;
  if(basename.size()){
    filename=buildFileName(version, prefix, postfix,
			   basename, defaultextension,
			   "",
			   removeold);
  }else{
    filename=buildFileName(version, prefix, postfix,
			   defaultname, defaultextension,
			   dirname,
			   removeold);
  }

  return filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::buildDefaultResultsFileName(int32 version, const string & prefix, const string & postfix, const string & basename, const string & defaultname, const string & defaultextension, bool removeold)
{
  string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_results;
  }

  string filename;
  if(basename.size()){
    filename=buildFileName(version, prefix, postfix,
			   basename, defaultextension,
			   "",
			   removeold);
  }else{
    filename=buildFileName(version, prefix, postfix,
			   defaultname, defaultextension,
			   dirname,
			   removeold);
  }

  return filename;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getContigReadListFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_crlist,
    ".txt");
}

void Assembly::saveContigReadList(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getContigReadListFilename(version, prefix, postfix, basename));
  assout::saveContigReadList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getStatisticsFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_contigstats,
    ".txt");
}

void Assembly::saveStatistics(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getStatisticsFilename(version, prefix, postfix, basename));
  assout::saveStatistics(AS_contigs,filename, deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getAssemblyInfoFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_info,
    ".txt");
}

void Assembly::saveAssemblyInfo(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getAssemblyInfoFilename(version, prefix, postfix, basename));
  assout::saveAssemblyInfo(AS_assemblyinfo,filename, deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getLargeContigsInfoFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_largecontigs,
    ".txt");
}

void Assembly::saveLargeContigsInfo(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getLargeContigsInfoFilename(version, prefix, postfix, basename));
  assout::saveLargeContigsInfo(AS_assemblyinfo,filename, deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveDebrisList(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  FUNCSTART("void Assembly::saveDebrisInfo(int32 version, const string & prefix, const string & postfix, const string & debrisfilename)");

  string filename(buildDefaultInfoFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_stats_debrislist,
		    ".txt"));

  cout << "Saving debris list to file: " << filename << endl;
  ofstream fout(filename.c_str(), ios::out | ios::trunc);

  for(uint32 i=0; i< AS_readpool.size(); i++){
    if(AS_isdebris[i]) {
      fout << AS_readpool.getRead(i).getName();
      switch(AS_isdebris[i]) {
      case DEBRIS_NOTDEBRIS : {
	// NOP, can never happen but quieten compiler warning
	break;
      }
      case DEBRIS_SHORTONLOAD : {
	fout << "\tSHORTONLOAD\n";
	break;
      }
      case DEBRIS_UNSPECIFIED : {
	fout << "\tUNSPECIFIED\n";
	break;
      }
      case DEBRIS_NOOVERLAP : {
	fout << "\tNO_OVERLAP\n";
	break;
      }
      case DEBRIS_NOTMAPPED : {
	fout << "\tNOT_MAPPED\n";
	break;
      }
      case DEBRIS_ABORTEDCONTIGCREATION : {
	fout << "\tABORTED_CONTIG_CREATION\n";
	break;
      }
      case DEBRIS_TINYCONTIG : {
	fout << "\tTINY_CONTIG\n";
	break;
      }
      case DEBRIS_TINYCLUSTER : {
	fout << "\tTINY_CLUSTER\n";
	break;
      }
      case DEBRIS_TINYCLUSTERORPHAN : {
	fout << "\tTINY_CLUSTER_ORPHAN\n";
	break;
      }
      case DEBRIS_UNSAVEDSINGLET : {
	fout << "\tUNSAVED_SINGLET\n";
	break;
      }
      case DEBRIS_DIGITAL_NORMALISATION : {
	fout << "\tDIGITAL_NORMALISATION\n";
	break;
      }
      case DEBRIS_CLIP_BADSOLEXAEND : {
	fout << "\tCLIP_BAD_SOLEXA_END\n";
	break;
      }
      case DEBRIS_CLIP_KNOWNADAPTORRIGHT : {
	fout << "\tCLIP_KNOWNADAPTORRIGHT\n";
	break;
      }
      case DEBRIS_CLIP_QUALMINTHRESHOLD : {
	fout << "\tCLIP_QUALMINTHRESHOLD\n";
	break;
      }
      case DEBRIS_CLIP_LOWERCASEFRONT : {
	fout << "\tCLIP_LOWERCASEFRONT\n";
	break;
      }
      case DEBRIS_CLIP_LOWERCASEBACK : {
	fout << "\tCLIP_LOWERCASEBACK\n";
	break;
      }
      case DEBRIS_CLIP_QUALCLIPS : {
	fout << "\tCLIP_QUALCLIPS\n";
	break;
      }
      case DEBRIS_CLIP_MASKEDBASES : {
	fout << "\tCLIP_MASKEDBASES\n";
	break;
      }
      case DEBRIS_CLIP_BADSEQUENCESERACH : {
	fout << "\tCLIP_BADSEQUENCESERACH\n";
	break;
      }
      case DEBRIS_CLIP_POLYBASEATEND : {
	fout << "\tCLIP_POLYBASEATEND\n";
	break;
      }
      case DEBRIS_CLIP_POLYAT : {
	fout << "\tCLIP_POLYAT\n";
	break;
      }
      case DEBRIS_CLIP_MINLEFTCLIP : {
	fout << "\tCLIP_MINLEFTCLIP\n";
	break;
      }
      case DEBRIS_CLIP_MINRIGHTCLIP : {
	fout << "\tCLIP_MINRIGHTCLIP\n";
	break;
      }
      case DEBRIS_CLIP_PHIX174 : {
	fout << "\tCLIP_PHIX174\n";
	break;
      }
      default : {
	fout << "\tNO_CODE_YET?_" << static_cast<uint16>(AS_isdebris[i]) << "\n";
      }
      }
    }
  }
  fout.close();

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getReadTagListFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultInfoFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_stats_readtags,
    ".txt");
}

void Assembly::saveReadTagList(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getReadTagListFilename(version, prefix, postfix, basename));
  assout::saveReadTagList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getConsensusTagListFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
   return buildDefaultInfoFileName(
     version, prefix, postfix, basename,
     AS_miraparams[0].getAssemblyParams().as_outfile_stats_contigtags,
     ".txt");
}

void Assembly::saveConsensusTagList(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
  {
  string filename(getConsensusTagListFilename(version, prefix, postfix, basename));
  assout::saveConsensusTagList(AS_contigs,filename,deleteoldfile);
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveSNPList(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(buildDefaultInfoFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_stats_snpanalysis,
		    ".txt"));
  assout::saveSNPList(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveFeatureAnalysis(int32 version, const string & prefix, const string & postfix, const string & faname, const string & fsname, const string & fcname, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveFeatureAnalysis(int32 version, const string & prefix, const string & postfix, const string & faname, const string & faname, bool deleteoldfile)");

  string dirname;
  if(version>=0){
    dirname=AS_miraparams[0].getDirectoryParams().dir_tmp;
  }else{
    dirname=AS_miraparams[0].getDirectoryParams().dir_info;
  }

  string filenamea;
  if(faname.size()){
    filenamea=buildFileName(version, prefix, postfix, faname, ".txt");
  }else{
    filenamea=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featureanalysis,
			    ".txt",
			   dirname);
  }

  string filenames;
  if(fsname.size()){
    filenames=buildFileName(version, prefix, postfix, fsname, ".txt");
  }else{
    filenames=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featuresummary,
			    ".txt",
			   dirname);
  }

  string filenamec;
  if(fcname.size()){
    filenamec=buildFileName(version, prefix, postfix, fcname, ".txt");
  }else{
    filenamec=buildFileName(version, prefix, postfix,
			    AS_miraparams[0].getAssemblyParams().as_outfile_stats_featuresequences,
			    ".txt",
			   dirname);
  }

  assout::saveFeatureAnalysis(AS_contigs,AS_readpool,
			      filenamea,filenames,filenamec,
			      deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getFASTAFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAUNPADDED,
    ".fasta");
}
string Assembly::getFASTAPaddedFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAPADDED,
    ".fasta");
}

void Assembly::saveAsFASTA(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getFASTAFilename(version, prefix, postfix, basename));
  string paddedfilename(getFASTAPaddedFilename(version, prefix, postfix, basename));
  assout::saveAsFASTA(AS_contigs,filename,paddedfilename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveStrainsAsFASTAQUAL(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveStrainsAsFASTAQUAL(int32 version, const string & prefix, const string & postfix, const string & fastaname)");

  string filename(buildDefaultResultsFileName(
		    version, prefix, postfix, basename,
		    AS_miraparams[0].getAssemblyParams().as_outfile_FASTAPADDED,
		    ""));
  assout::saveStrainsAsFASTAQ(AS_contigs,AS_readpool,
			      filename,
			      false,0,0,
			      deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getTCSFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_TCS,
    ".tcs");
}
void Assembly::saveAsTCS(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsTCS(int32 version, const string & prefix, const string & postfix, const string & tcsname)");

  string filename(getTCSFilename(version, prefix, postfix, basename));
  assout::saveAsTCS(AS_contigs,filename,deleteoldfile);

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getCAFFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_CAF,
    ".caf");
}

void Assembly::saveAsCAF(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getCAFFilename(version, prefix, postfix, basename));
  assout::saveAsCAF(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getMAFFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_MAF,
    ".maf");
}

void Assembly::saveAsMAF(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getMAFFilename(version, prefix, postfix, basename));
  assout::saveAsMAF(AS_contigs,filename,deleteoldfile);
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getTXTFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_TXT,
    ".txt");
}

void Assembly::saveAsTXT(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string filename(getTXTFilename(version, prefix, postfix, basename));
  assout::saveAsTXT(AS_contigs,filename,deleteoldfile);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getACEFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_ACE,
    ".ace");
}
void Assembly::saveAsACE(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsACE(int32 version, const string & prefix, const string & postfix, const string & acename)");

  string filename(getACEFilename(version, prefix, postfix, basename));
  assout::saveAsACE(AS_contigs,filename,deleteoldfile);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getWiggleFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_WIGGLE,
    ".wig");
}
void Assembly::saveAsWiggle(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsWiggle(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)");

  string filename(getWiggleFilename(version, prefix, postfix, basename));
  assout::saveAsWiggle(AS_contigs,filename,deleteoldfile,false);

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getGAP4DAFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outdir_GAP4DA,
    ".gap4da");
}

void Assembly::saveAsGAP4DA(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  string subdirname(getGAP4DAFilename(version, prefix, postfix, basename));
  assout::saveAsGAP4DA(AS_contigs,subdirname,deleteoldfile);
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::getHTMLFilename(int32 version, const string & prefix, const string & postfix, const string & basename)
{
  return buildDefaultResultsFileName(
    version, prefix, postfix, basename,
    AS_miraparams[0].getAssemblyParams().as_outfile_HTML,
    ".html");
}

void Assembly::saveAsHTML(int32 version, const string & prefix, const string & postfix, const string & basename, bool deleteoldfile)
{
  FUNCSTART("void Assembly::saveAsHTML(int32 version, const string & prefix, const string & postfix, const string & htmlname, bool deleteoldfile)");

  string filename(getHTMLFilename(version, prefix, postfix, basename));

  string projectname(AS_miraparams[0].getAssemblyParams().as_projectname_out);

  cout << "Saving contigs to file: " << filename << endl;

  //ofstream htmlout(filename.c_str(), ios::out | ios::trunc);
  assout::dumpContigListAsHTML(AS_contigs, filename, deleteoldfile, projectname);

  FUNCEND();
}
