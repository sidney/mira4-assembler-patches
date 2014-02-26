/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2006 and later by Bastien Chevreux
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

#include "contig.H"
#include "util/progressindic.H"

using namespace std;


//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif


#ifdef PARANOIABUGTRACKFLAG
#define paranoiaBUGSTAT(statement) { statement;}
#define paranoiaBUGIF(ifcond, statement) { if(ifcond) {statement;}}
#else
#define paranoiaBUGSTAT(statement)
#define paranoiaBUGIF(ifcond, statement)
#endif




bool Contig::hasEditableOvercallData() const
{
  FUNCSTART("bool Contig::hasEditableOvercallData() const");

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  for(uint32 rgid=0; rgid<CON_readsperreadgroup.size(); ++rgid){
    if(CON_readsperreadgroup[rgid]>0
       && ((*CON_miraparams)[ReadGroupLib::getReadGroupID(rgid).getSequencingType()].getEditParams().ed_homopolymer_overcalls)) return true;
  }

  return false;
}

bool Contig::hasSeqTypeData(uint8 seqtype) const
{
  FUNCSTART("bool Contig::hasSeqTypeData(uint8 seqtype) const");

  for(uint32 rgid=0; rgid<CON_readsperreadgroup.size(); ++rgid){
    if(CON_readsperreadgroup[rgid]>0
       && ReadGroupLib::getReadGroupID(rgid).getSequencingType() == seqtype) return true;
  }

  return false;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)
//#define CEBUG(bla)   {cout << bla; cout.flush();}
uint32 Contig::createPacBioDarkStrobeEdits(list<pbdse_t> & pbdsev)
{
  FUNCSTART("uint32 Contig::createPacBioDarkStrobeEdits(list<pbdse_t> & pbdse)");

  BUGIFTHROW(true,"need redo ce1 for PlacedContigReads");

////++++////  pbdsev.clear();
////++++////
////++++////  int32 maxce=0;
////++++////
////++++////  string consseq;
////++++////  vector<base_quality_t> consqual;
////++++////
////++++////  newConsensusGet(consseq, consqual);
////++++////
////++++////  vector<contigread_t>::const_iterator crI=CON_reads.begin();
////++++////
////++++////    int32 dbg_goodestim=0;
////++++////    int32 dbg_medestim=0;
////++++////    int32 dbg_badestim=0;
////++++////    int32 dbg_sumgoodestim=0;
////++++////    int32 dbg_summedestim=0;
////++++////    int32 dbg_sumbadestim=0;
////++++////
////++++////  for(; crI != CON_reads.end(); crI++){
////++++////    if(crI->orpid<0
////++++////       // TODO: what aboutPacBioLQ ?
////++++////       || !crI->read.isSequencingType(ReadGroupLib::SEQTYPE_PACBIOHQ)) continue;
////++++////
////++++////    const Read & actread = crI->read;
////++++////
////++++////    BUGIFTHROW(actread.getLenClippedSeq() == 0, "Read: " << actread.getName() << " has 0 length?");
////++++////
////++++////    int32 cpos=crI->offset;
////++++////    int32 cposincr=1;
////++++////    if(crI->direction < 0) {
////++++////      cposincr=-1;
////++++////      cpos+=actread.getLenClippedSeq()-1;
////++++////    }
////++++////    int32 arpos=crI->contigPos2RawReadPos(cpos);
////++++////
////++++////    auto ccI=CON_counts.begin();
////++++////    advance(ccI,cpos);
////++++////
////++++////    char rbase;
////++++////    char cbase;
////++++////    int32 runcounter=0;
////++++////    int32 ncounter=0;
////++++////    int32 gapcounter=0;
////++++////    int32 changeestim=0;
////++++////
////++++////    bool runend=false;
////++++////
////++++////    // start/stop of run in read (unclipped read pos)
////++++////    int32 arposs=0;
////++++////    int32 arpose=0;
////++++////    for(uint32 ri=0; ri<actread.getLenClippedSeq(); ri++, arpos++, cpos+=cposincr, ccI+=cposincr){
////++++////      if(crI->direction < 0) {
////++++////	//rbase=actread.getBaseInComplementSequence(arpos);
////++++////	rbase=actread.getBaseInSequence(arpos);
////++++////      }else{
////++++////	rbase=actread.getBaseInSequence(arpos);
////++++////      }
////++++////      //cout << actread.getName() << " " << arpos << " " << rbase << "\t" << cpos << '\t' << consseq[cpos] << '\t' << changeestim << endl;
////++++////
////++++////      cbase=consseq[cpos];
////++++////      if(cbase == '*'){
////++++////	if(ccI->N > ccI->star) cbase ='N';
////++++////      }
////++++////
////++++////      if(runcounter>0){
////++++////	if(rbase=='*'){
////++++////	  runcounter++;
////++++////	  gapcounter++;
////++++////	  if(cbase!='*'){
////++++////	    changeestim++;
////++++////	  }
////++++////	}else if(rbase=='N' ||rbase=='n'){
////++++////	  runcounter++;
////++++////	  ncounter++;
////++++////	  if(cbase=='*'){
////++++////	    changeestim--;
////++++////	  }
////++++////	}else{
////++++////	  runend=true;
////++++////	  arpose=arpos-1;
////++++////	}
////++++////      }else{
////++++////	if(rbase=='*'){
////++++////	  arposs=arpos;
////++++////	  runcounter++;
////++++////	  gapcounter++;
////++++////	  if(cbase!='*'){
////++++////	    changeestim++;
////++++////	  }
////++++////	}else if(rbase=='N' ||rbase=='n'){
////++++////	  arposs=arpos;
////++++////	  runcounter++;
////++++////	  ncounter++;
////++++////	  if(cbase=='*'){
////++++////	    changeestim--;
////++++////	  }
////++++////	}
////++++////      }
////++++////      if(runend){
////++++////	CEBUG("runend. rc: " << runcounter << "\tgc: " << gapcounter << "\tnc: " << ncounter << "\tce: " << changeestim << "\ts: " << arposs << "\te: " << arpose << endl);
////++++////	if(ncounter>=5 && changeestim !=0){
////++++////	  // on 100/100 data
////++++////	  //direction good in 80% of the time when using 2/3 of change estimate
////++++////	  // on 100/100 data
////++++////	  // looks like 3/5 is sweet spot
////++++////
////++++////	  //if(abs(changeestim)>2) {
////++++////	  //  changeestim=changeestim*3/5;
////++++////	  //} else if(abs(changeestim)==2) {
////++++////	  //  changeestim/=2;
////++++////	  //}
////++++////	  if(abs(changeestim)>3) {
////++++////	    changeestim=changeestim*3/5;
////++++////	  } else if(abs(changeestim)!=1) {
////++++////	    changeestim/=2;
////++++////	  }
////++++////
////++++////	  {
////++++////	    pbdse_t ptmp;
////++++////	    ptmp.rid=crI->orpid;
////++++////	    ptmp.rposs=arposs;
////++++////	    ptmp.rpose=arpose;
////++++////	    ptmp.changeestim=changeestim;
////++++////	    pbdsev.push_back(ptmp);
////++++////	  }
////++++////	  CEBUG("Stored " << actread.getName() << "\tce: " << changeestim << endl);
////++++////
////++++////	  maxce=max(abs(changeestim),maxce);
////++++////	  int32 newestim=ncounter+changeestim;
////++++////	  if(abs(100-ncounter) > abs(100-newestim)){
////++++////	    dbg_goodestim++;
////++++////	    dbg_sumgoodestim+=changeestim;
////++++////	  }else if(abs(100-ncounter) == abs(100-newestim)){
////++++////	    dbg_medestim++;
////++++////	    dbg_summedestim+=changeestim;
////++++////	  }else{
////++++////	    dbg_badestim++;
////++++////	    dbg_sumbadestim+=changeestim;
////++++////	  }
////++++////	}
////++++////	runcounter=0;
////++++////	gapcounter=0;
////++++////	ncounter=0;
////++++////	changeestim=0;
////++++////	arposs=0;
////++++////	arpose=0;
////++++////	runend=false;
////++++////      }
////++++////    }
////++++////  }
////++++////
////++++////  CEBUG("Good: " << dbg_goodestim << "\tSum: " << dbg_sumgoodestim
////++++////	<< "\nMed : " << dbg_medestim << "\tSum: " << dbg_summedestim
////++++////	<< "\nBad : " << dbg_badestim << "\tSum: " << dbg_sumbadestim
////++++////	<< "\nMax : " << maxce
////++++////	<< endl);
////++++////
////++++////  //pbdsev.sort(Contig::pbdse_t_comparator);
////++++////
////++++////  FUNCEND();
////++++////  return static_cast<uint32>(maxce);
  return 0;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * check stars only / N only columns and remove them
 * readjust the offsets if needed
 * from = including, to= excluding! (normal C like interval)
 *
 * returns number of deleted columns
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
uint32 Contig::deleteStarOnlyColumns(int32 from, int32 to, bool alsononly, uint32 mincov)
{
  FUNCSTART("void Contig::deleteStarOnlyColumns(int32 from, int32 to, bool alsononly, uint32 mincov, bool updateconcounts)");

  // Problem: mapping assemblies with lots of gap columns (Ion!)
  // For "normal" mapping assemblies, one reference read coveres everything
  // Which means: the for loop over the reads will always start at the first read
  //  and sift over a lot of non-interesting reads before hitting the ones which
  //  are covering a gap column. Lots of time wastage, especially in high
  //  coverage contigs.
  //
  // Resolution:
  // For contigs with just one backbone read, cheat. Find the backbone read beforehand,
  //  and in the for loop, do not use getFirstPCRIForReadsCoveringPosition() but
  //  getFirstPCRIForNonBBReadsCoveringPosition().
  //  Handle the backbone read separately.
  //
  // Contigs with >1 bb read are rare (would have needed to use CAF/MAF alignment mapping
  //  with one long grand-parent reference as reference, not many users do that), but these
  //  have to bite the bullet and go through reads normally.

  int32 checklen=to-from;
  if(checklen==0) return 0;
  uint32 numdel=0;

  paranoiaBUGSTAT(checkContig());

  BUGIFTHROW(from>=CON_counts.size(), "from>=CON_counts.size() ???");
  if(to>CON_counts.size()) to=CON_counts.size();
  BUGIFTHROW(to<from, "to " << to << " < from " << from << " ?");

  if(CON_fixedconsseq.size()){
    nukeSTLContainer(CON_fixedconsseq);
    nukeSTLContainer(CON_fixedconsqual);
  }
  definalise();

  CEBUG("Starcheck from: " << from << " to " << to << " (excluding)\nAlso N only" << alsononly <<endl);

  PlacedContigReads::const_iterator firstbbI(CON_reads.end());
  PlacedContigReads::const_iterator pcrJ=firstbbI;

  auto numbbreads=getNumBackbones();
  if(numbbreads>0){
    for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end(); ++pcrI){
      if(pcrI->isBackbone()) {
	firstbbI=pcrI;
	break;
      }
    }
  }

  CEBUG("Num bb reads: " << numbbreads << endl);

  CEBUG("Going in with:\n"; dumpAsText(cout,from,to));

  auto ccI= CON_counts.begin();
  BOUNDCHECK(from, 0, CON_counts.size());
  advance(ccI, from);

  // take care if inserting "continue" statements into this loop
  // ++ccI and ++actread are at end of loop!!!
  for(int32 actcontigpos=from; actcontigpos < to; ){
    bool todelete=false;
    CEBUG("dsoc acp: " << actcontigpos << endl);// << *ccI << endl);
    if(ccI->total_cov >= mincov
       && (ccI->getBBChar() == '*'
	   || ccI->getBBChar() == '@')){
      if(ccI->star==ccI->total_cov){
	CEBUG("Complete star row. Deleting star in consensus and in other reads.\n");
	todelete=true;
      } else if(alsononly){
	if(ccI->A == ccI->N
	   && ccI->C == ccI->N
	   && ccI->G == ccI->N
	   && ccI->T == ccI->N
	   && ccI->X == 0){
	  CEBUG("Complete N row. Deleting in consensus and in other reads.\n");
	  todelete=true;
	}
      }
    }

    if(todelete){
      ++numdel;
      if(numbbreads==1){
	pcrJ=getFirstPCRIForNonBBReadsCoveringPosition(actcontigpos);

	// delete gaps / N in backbone read
	if(actcontigpos >= firstbbI.getReadStartOffset()
	   && actcontigpos < (firstbbI.getReadStartOffset() + firstbbI->getLenClippedSeq())){
	  CEBUG("BB in range, deleting.");
	  // internaloffset is the offset of the * / N in this read
	  int32 internaloffset= actcontigpos-firstbbI.getReadStartOffset();
	  CEBUG("\tInternal offset: " << internaloffset);
	  if(firstbbI.getReadDirection() > 0){
	    const_cast<Read &>(*firstbbI).deleteBaseFromClippedSequence(internaloffset);
	  }else{
	    const_cast<Read &>(*firstbbI).deleteBaseFromClippedComplementSequence(internaloffset);
	  }
	}
      }else{
	pcrJ=getFirstPCRIForReadsCoveringPosition(actcontigpos);
      }

      CEBUG("TODELETE! " << actcontigpos << "\npcrJ points at " << pcrJ->getName() << endl);
      // Deleting stars / Ns in other reads
      for(; pcrJ != CON_reads.end() && pcrJ.getReadStartOffset() <= actcontigpos; ++pcrJ){
	if(!(numbbreads == 1 && pcrJ->isBackbone()) // remember: 1 bb read -> already handled
	   && actcontigpos >= pcrJ.getReadStartOffset()
	   && actcontigpos < (pcrJ.getReadStartOffset() + pcrJ->getLenClippedSeq())){
	  CEBUG("In range " << pcrJ->getName() << " at " << pcrJ.getReadStartOffset() << " to " << pcrJ.getReadStartOffset() + pcrJ->getLenClippedSeq());
	  // internaloffset is the offset of the * / N in this read
	  int32 internaloffset= actcontigpos-pcrJ.getReadStartOffset();
	  CEBUG("\tInternal offset: " << internaloffset << '\n');
	  if(pcrJ.getReadDirection() > 0){
	    const_cast<Read &>(*pcrJ).deleteBaseFromClippedSequence(internaloffset);
	  }else{
	    const_cast<Read &>(*pcrJ).deleteBaseFromClippedComplementSequence(internaloffset);
	  }
	}
      }

      CEBUG("shiftReads(" << actcontigpos << ",-1)\n");
      CON_reads.shiftReads(actcontigpos,-1);
      shiftMarkerPositions(actcontigpos,-1);

      // push down consensus tags
      updateTagBaseDeleted(actcontigpos);

      // erase the place in CON_counts
      CEBUG("counts_offset: " << (ccI-CON_counts.begin()) << endl);
      CEBUG("ccI before: " << ccI << endl);
      ccI=CON_counts.erase(ccI);
      CEBUG("ccI after: " << ccI << endl);
      CEBUG("After DEL " << actcontigpos << ":\n"; dumpAsText(cout,from,to));

      // do not increase actcontigpos and ccI here, we stay at the same position
      //  as the old gap column disappeared
      // But: "to" needs to be decreased
      --to;
    }else{
      ++ccI;
      ++actcontigpos;
    }
  }

  paranoiaBUGSTAT(checkContig());

  FUNCEND();
  return numdel;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/


struct simplebasecounter_t {
  char base;
  uint32 counter;
};


//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint32 Contig::editTrickyOvercalls(const bool onlysetPSHPtags, const bool noSRMreadsallowed, vector<bool> & readsmarkedsrm)
{
  FUNCSTART("uint32 Contig::editTrickyOvercalls(const bool onlysetPSHPtags, const bool noSRMreadsallowed, vector<bool> & readsmarkedsrm)");

  if(onlysetPSHPtags){
    cout << "Marking";
  }else{
    cout << "Editing";
    definalise();
  }
  cout << " tricky overcalls (";
  if(noSRMreadsallowed){
    cout << "no ";
  }
  cout << "SRMs allowed):\n";

  bool needcheckSRM=false;
  if(readsmarkedsrm.size() != CON_readpool->size()){
    // if we were called with an incorrectly sized readsmarkedsrm vector
    //  (empty or whatever), assume no SRM freshly set
    readsmarkedsrm.clear();
    readsmarkedsrm.resize(CON_readpool->size(),false);
  }else{
    for(uint32 i=0; i<readsmarkedsrm.size(); ++i){
      if(readsmarkedsrm[i]){
	needcheckSRM=true;
	break;
      }
    }
  }

  // reserve enough space for one edit every 10th base in the consensus
  vector<edit454command_t> all454editcommands;
  all454editcommands.reserve(CON_counts.size()/10);

  // initialise the iterator for getting through the contig
  rcci_t rcci(this);
  {
    vector<int32> allowedstrainids; // empty would be all ids
    vector<uint8> allowedreadtypes;
    //allowedreadtypes.push_back(ReadGroupLib::SEQTYPE_454GS20);
    rcci.init(allowedstrainids,
	      allowedreadtypes,
	      false,            // don't take rails
	      false,           // nor backbones
	      false);   // nor reads without readpool-reads
  }

  vector<simplebasecounter_t> basecounter(5);
  {
    static string bases("ACGT*");
    for(uint32 i=0;i<basecounter.size(); i++){
      basecounter[i].base=bases[i];
      basecounter[i].counter=0;
    }
  }

  uint32 editlock=0;

  auto ccI=CON_counts.begin();
  ProgressIndicator<int32> P(0, static_cast<int32>(CON_counts.size())-1);

  for(uint32 actcontigpos=0; actcontigpos<CON_counts.size(); ++ccI, rcci.advance(), ++actcontigpos){

    P.progress(actcontigpos);

    if(editlock>0) {
      --editlock;
      continue;
    }

    // check for disagreement in this column
    uint32 setcount=0;

    if(ccI->A > 0) ++setcount;
    if(ccI->C > 0) ++setcount;
    if(ccI->G > 0) ++setcount;
    if(ccI->T > 0) ++setcount;
    if(ccI->star > 0) ++setcount;

    // if more than 1, then there are bases in this column which disagree
    //  (caution: this might be also just an "N", or some other IUPACs)
    // tricky == those with * in column
    if(setcount>1 && ccI->star > 0) {
      if(needcheckSRM){
	// have a look at all reads at this position: if one of them was
	//  marked in this pass with a SRM, we don't do the 454 edit here
	//  (because some "corrections" might be due to misassemblies and
	//  would be plain wrong)
	bool readfreshlymarkedsrm=false;
	for(auto & pcrI : rcci.getPCRIsInCol()){
	  if(readsmarkedsrm[pcrI.getORPID()]){
	    readfreshlymarkedsrm=true;
	    break;
	  }
	}
	if(readfreshlymarkedsrm) {
	  continue;
	}
      }

      CEBUG("cpos: " << actcontigpos << "\t" << ccI->A << " " << ccI->C << " " << ccI->G << " " << ccI->T << " " << ccI->star << endl);

      for(auto & bce : basecounter){
	bce.counter=0;
      }
      for(auto & pcrI : rcci.getPCRIsInCol()){
	CEBUG("read: " << pcrI->getName());
	char  base=pcrI.getBase(actcontigpos);
	CEBUG("\t" << base << endl);

	for(auto & bce : basecounter){
	  if(bce.base==base) ++bce.counter;
	}
      }


      {
	size_t maxsize=0;
	int32 maxsize_i=-1;
	size_t runnerup=0;
	int32 runnerup_i=-1;

	for(uint32 bci=0; bci < basecounter.size(); ++bci){
	  CEBUG("ACGT*"[bci] << "\t" << basecounter[bci].counter << " ");
	  if(basecounter[bci].counter>=maxsize){
	    runnerup=maxsize;
	    runnerup_i=maxsize_i;
	    maxsize=basecounter[bci].counter;
	    maxsize_i=bci;
	  }else if(basecounter[bci].counter>=runnerup){
	    runnerup=basecounter[bci].counter;
	    runnerup_i=bci;
	  }
	}

	CEBUG("\nmaxsize: " << maxsize << "\ti: " << maxsize_i << "\n");
	CEBUG("runnerup: " << runnerup << "\ti: " << runnerup_i << "\n");


	// index==4 is gap
	if((maxsize_i==4 && runnerup_i>=0)
	   || (maxsize_i>=0 && runnerup_i==4)){
	  if(maxsize+runnerup == 0){
	    //addTagToConsensus(actcontigpos, actcontigpos,'=',"T454","DOH?");
	  }else{
	    if(100*runnerup/(maxsize+runnerup) >= 40){
	      ostringstream ostr;
	      ostr << static_cast<char>(basecounter[maxsize_i].base) << ": " << maxsize;
	      ostr << " " << static_cast<char>(basecounter[runnerup_i].base) << ": " << runnerup;
	      ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
	      addTagToConsensus(actcontigpos,
				actcontigpos,
				'=',
				multitag_t::getIdentifierStr(CON_tagentry_idDGPc).c_str(),
				ostr.str().c_str(),
				true);

	      CEBUG(ostr.str() << '\n');
	    }

	    if(basecounter[maxsize_i].base=='*'){
	      editlock=edit454checkTrickies(basecounter[runnerup_i].base,
					    actcontigpos,
					    rcci.getPCRIsInCol(),
					    all454editcommands,
					    onlysetPSHPtags,
					    noSRMreadsallowed);
	    }else{
	      editlock=edit454checkTrickies(basecounter[maxsize_i].base,
					    actcontigpos,
					    rcci.getPCRIsInCol(),
					    all454editcommands,
					    onlysetPSHPtags,
					    noSRMreadsallowed);
	    }
	    if(editlock>0) editlock+=2;
	  }
	}
      }
    }
  }
  P.finishAtOnce();

  CEBUG("Generated " << all454editcommands.size() << " tricky mark/edit commands.\n");

  uint32 retval=all454editcommands.size();
  if(all454editcommands.size()) {
    uint32 numwedits=0;
    static multitag_t tagPSHP;
    tagPSHP.setIdentifierStr("PSHP");
    tagPSHP.source=multitag_t::MT_tagsrcentry_idMIRA;
    static multitag_t tagR454;
    tagR454.setIdentifierStr("R454");
    tagR454.source=multitag_t::MT_tagsrcentry_idMIRA;

    if(onlysetPSHPtags){
      cout << "Marking tricky overcall runs in " << all454editcommands.size() << " cases.\n";
    }else{
      cout << "Checking for allowed edits in " << all454editcommands.size() << " cases.\n";
    }
    sort(all454editcommands.begin(),
	 all454editcommands.end(), Contig::edit454command_t_comparator);

    for(uint32 aeci=0; aeci<all454editcommands.size(); ++aeci){
      auto & pcrI =all454editcommands[aeci].pcrI;

      if(onlysetPSHPtags){
	CEBUG("Mark read: " << pcrI->getName());
      }else{
	if(((*CON_miraparams)[pcrI->getSequencingType()].getEditParams().ed_homopolymer_overcalls)){
	  CEBUG("EDIT read: " << pcrI->getName());
	}else{
	  // user does not want edits in this type of reads
	  continue;
	}
      }
      CEBUG("\tbase: " << all454editcommands[aeci].base);
      CEBUG("\tpos: " << all454editcommands[aeci].readpos << endl);

      if(onlysetPSHPtags){
	uint32 zeroqualcounts=0;
	uint32 runfrom=0;
	uint32 runto=0;
	try{
	  uint32 runlength=getBaseRunLength(*pcrI,
					    all454editcommands[aeci].readpos,
					    all454editcommands[aeci].base,
					    runfrom,
					    runto,
					    zeroqualcounts,
					    true);
	  if(runlength){
	    if(runfrom > 0) --runfrom;
	    if(runto < pcrI->getLenSeq()-1) ++runto;
	    tagPSHP.from=runfrom;
	    tagPSHP.to=runto;
	    const_cast<Read &>(*pcrI).addTagO(tagPSHP);
	  }
	}
	catch(...){
	}
      }else{
	const_cast<Read &>(*pcrI).deleteWeakestBaseInRun(all454editcommands[aeci].base,
							 all454editcommands[aeci].readpos,
							 true);

	++numwedits;

	// just in case addTag decided to throw ...
	// too lazy to do it right at the moment
	try{
	  tagR454.from=all454editcommands[aeci].readpos-1;
	  tagR454.to=all454editcommands[aeci].readpos+1;
	  const_cast<Read &>(*pcrI).addTagO(tagR454);
	}
	catch(...){
	}
      }

    }

    if(!onlysetPSHPtags){
      //cout << "Performed " << numwedits << " edits.\n";
      retval=numwedits;
    }
  }

  cout << "\n";


  FUNCEND();

  return retval;
}



struct trickyrunshistogram_t {
  vector<PlacedContigReads::const_iterator> seen_pcrIs;
  vector<uint32> realreadpos_of_reads;
  vector<char> realbasehypo_of_reads;
  uint32 count_zeroqualadjusted;
  uint32 count_havezeroqual;
};



// returns length of base run that is looked at if edited
// else 0
uint32 Contig::edit454checkTrickies(const char basehypo, const uint32 actcontigpos, const vector<PlacedContigReads::const_iterator> & pcrIs_in_col, vector<edit454command_t> & editcommands, const bool onlysetPSHPtags, const bool noSRMreadsallowed)
{
  FUNCSTART("uint32 Contig::edit454checkTrickies(const char basehypo, const uint32 actcontigpos, const vector<PlacedContigReads::const_iterator> & pcrIs_in_col, vector<edit454command_t> & editcommands, const bool onlysetPSHPtags, const bool noSRMreadsallowed)");

  CEBUG("Checking tricky hypothesis: " << basehypo << " " << actcontigpos << endl);

  if(noSRMreadsallowed){
    for(auto & pcrI : pcrIs_in_col){
      if(pcrI->hasTag(Read::REA_tagentry_idSRMr)){
	CEBUG("Read with SRM found, but this is presently not allowed at this stage.\n");
	FUNCEND();
	return 0;
      }
    }
  }


  vector<trickyrunshistogram_t> trh;
  trh.reserve(100);

  // After stepping through all reads, the following variable shows
  //  whether additional bases in this column are all equal to "basehypo"
  bool cleanbasehypocolumn=true;

  uint32 maxspan=0;

  auto & ed_params = (*CON_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getContigParams();
  bool lookedatread=false;
  for(auto & pcrI : pcrIs_in_col){
    // do not analyse backbones and rails
    // look also at reads the user does not want to have edited ... they can tilt the balance toward the
    //  right number for those which need an edit!
    if(pcrI->isBackbone()
       || pcrI->isRail()) continue;
    //|| !((*CON_miraparams)[pcrI->getSequencingType()].getEditParams().ed_homopolymer_overcalls)) continue;

    lookedatread=true;

    int32 readpos=pcrI.contigPos2UnclippedReadPos(actcontigpos);

    CEBUG("read: " << pcrI->getName()) ; cout.flush();

    // to simplify life, we´re going to have this routine
    //  only work with the forward direction of reads
    // for reads in reverse direction in contig, calculate the "real" values
    int32 realreadpos=readpos;
    char realbasehypo=basehypo;
    if(pcrI.getReadDirection() < 0){
      realbasehypo=dptools::getComplementBase(realbasehypo);
      realreadpos=pcrI->calcComplPos(readpos);
    }


    CEBUG(" " << pcrI->getBaseInSequence(realreadpos) << " rbh: " << realbasehypo << " rrp: " << realreadpos << endl);

    if(pcrI->getBaseInSequence(realreadpos) != '*'
       && pcrI->getBaseInSequence(realreadpos) != realbasehypo){
      cleanbasehypocolumn=false;
    }

    uint32 zeroqualcounts=0;
    uint32 runfrom=0;
    uint32 runto=0;
    CEBUG("search runlength" << endl);
    uint32 runlength=getBaseRunLength(*pcrI,
				      realreadpos,
				      realbasehypo,
				      runfrom,
				      runto,
				      zeroqualcounts,
				      true);
    CEBUG("RL " << pcrI->getName() << ": " << runlength << " " << zeroqualcounts << endl);

    uint32 span=runto-runfrom;
    if(span>maxspan) maxspan=span;

    if(runlength>=trh.size()){
      CEBUG("extend histogramm" << endl);
      size_t initstart=trh.size();
      trh.resize(runlength+1);
      for(;initstart<trh.size(); initstart++){
	trh[initstart].count_zeroqualadjusted=0;
	trh[initstart].count_havezeroqual=0;
      }
    }
    CEBUG("rl in histogramm" << endl);
    trh[runlength].seen_pcrIs.push_back(pcrI);
    trh[runlength].realreadpos_of_reads.push_back(realreadpos);
    trh[runlength].realbasehypo_of_reads.push_back(realbasehypo);
    if(zeroqualcounts>0) trh[runlength].count_havezeroqual++;;

    // we will delete only one zeroqual base
    if(zeroqualcounts>0 && runlength>0){
      trh[runlength-1].count_zeroqualadjusted++;
    } else {
      trh[runlength].count_zeroqualadjusted++;
    }
    CEBUG("done" << endl);
  }

  bool edit=false;
  if(lookedatread){
    CEBUG("Cleanbasehypocolumn: " << cleanbasehypocolumn);
    CEBUG("\nsearch zqa\n");

    // search for maximum count_zeroqualadjusted in trh
    //  (but not in the zero-length counts)
    uint32 maxcountzqa=0;
    uint32 maxcountzqa_index=0;
    for(uint32 i=0; i<trh.size(); i++){
      CEBUG(i << "\t" << trh[i].seen_pcrIs.size() << "\t" << trh[i].count_zeroqualadjusted<< "\t" << trh[i].count_havezeroqual << endl);
      // changing > to >= favours longer runs with the exact same prediction
      //  of adjusted reads
      if(trh[i].count_zeroqualadjusted>=maxcountzqa){
	maxcountzqa=trh[i].count_zeroqualadjusted;
	maxcountzqa_index=i;
      }
    }

    CEBUG("search higher" << endl);

    uint32 runsatmczqa=static_cast<uint32>(trh[maxcountzqa_index].seen_pcrIs.size());
    uint32 runslargermczqa=0;
    uint32 runswillingtoshorten=0;
    uint32 expectedlarger=0;
    for(uint32 i=maxcountzqa_index+1; i<trh.size(); i++){
      runslargermczqa+=static_cast<uint32>(trh[i].seen_pcrIs.size());
      expectedlarger+=trh[i].count_zeroqualadjusted;
      runswillingtoshorten+=trh[i].count_havezeroqual;
    }

    ostringstream ostr;

    const uint32 ratiomultiplier=10;
    //cout << "TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO \n";
    //cout << "TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO \n";
    //cout << "TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO \n";
    //cout << "Change multiplier back to 10!\n";

    // General routines destined to get out most annoying things
    // start them if edit not triggered previously
    if(!edit){
      CEBUG("General rules.\n");
      CEBUG("trh.size(): " << trh.size());
      CEBUG("\nmaxcountzqa_index: " << maxcountzqa_index);
      CEBUG("\ncleanbasehypocolumn: " << cleanbasehypocolumn);
      CEBUG("\npcrIs_in_col.size(): " << pcrIs_in_col.size());
      CEBUG("\nrunsatmczqa: " << runsatmczqa);
      CEBUG('\n');

      if(trh.size()>maxcountzqa_index+2) {
	CEBUG("Funny looking histogram, no edit!\n");
	ostr << "Funny histogram ";
      }else{
	if(runslargermczqa>0){
	  // we have some kind of overcalls
	  if(runsatmczqa >= runslargermczqa*5){
	    /* if number of runs of probable length > 10x of the higher runs,
	       edit. E.g.
	       .....aAAA.....
	       .....*AAA.....
	       .....*AAA.....
	       .....*AAA.....
	       .....*AAA.....
	       etc.
	    */
	    edit=true;
	    ostr << "Smoking gun  ";
	  }
	} else if(maxcountzqa_index==1 && cleanbasehypocolumn
		  && pcrIs_in_col.size() >= runsatmczqa*5){
	  /* if it's a some unmotivated bases hanging around, being all
	     the same and less than 10% of the column. E.g.

	     Edit             Edit               No edit
	     ....TaGGG.....   ....TaGGG..... 	....TaGGG.....
	     ....T*GGG.....   ....TaGGG..... 	....TcGGG.....
	     ....T*GGG.....   ....T*GGG..... 	....T*GGG.....
	     ....T*GGG.....   ....T*GGG..... 	....T*GGG.....
	     ....T*GGG.....   ....T*GGG..... 	....T*GGG.....
	     etc.             etc.               etc.
	  */
	  edit=true;
	  ostr << "Lone ranger  ";
	  maxcountzqa_index=0;
	}
      }
    }
    if(edit) {
      CEBUG("Editing!\n");
      ostr << "Edit: " << basehypo;
      if(trh.size()>maxcountzqa_index+2) {
	/* ??? strategie leicht ändern in: wenn auch edits bei
	   maxcountzqa_index+2, dann erst die */
	ostr << " Caution!";
	CEBUG("CAUTION!\n");
	maxcountzqa_index++;
      }
    }else{
      ostr << "No edit.";
    }
    ostr << " PL: " << maxcountzqa_index;      // probable length
    ostr << " RS: " << runsatmczqa;            // reads that are this length
    ostr << " RL: " << runslargermczqa;        // reads that are larger
    ostr << " RWTS: " << runswillingtoshorten; // reads willing to shorten
    ostr << " EL: " << expectedlarger;         // reads expected to be larger
    //  after editing

    addTagToConsensus(actcontigpos, actcontigpos,'=',"H454",ostr.str().c_str(),true);

    CEBUG(ostr.str() << '\n');

    if(edit){
      // current strategy:
      // now that we know that we need to shorten a few reads at this
      //  position, shorten *all* reads larger than the length we determined
      //  (maxcountzqa_index) by one base (the weakest base). However, the
      //  base quality of the weakest base now does not need to be 0.

      uint32 acttrh=maxcountzqa_index+1;

      // addition: if only marking PSHP runs, mark all reads at that position
      if(onlysetPSHPtags){
	acttrh=0;
      }

      for(; acttrh<trh.size(); ++acttrh){
	uint32 readnr=0;
	for(auto & pcrI : trh[acttrh].seen_pcrIs){
	  CEBUG("Saving EDIT: read: " << pcrI->getName() << " " << pcrI->getBaseInSequence(trh[acttrh].realreadpos_of_reads[readnr]) << " (" << trh[acttrh].realbasehypo_of_reads[readnr] << ") at " << trh[acttrh].realreadpos_of_reads[readnr] << endl);
	  editcommands.push_back(edit454command_t(pcrI,
						  pcrI.getURDID(),
						  trh[acttrh].realreadpos_of_reads[readnr], //realreadpos
						  trh[acttrh].realbasehypo_of_reads[readnr]) //realbasehypo
	    );
	  ++readnr;
	}
      }
    }
  }

  FUNCEND();
  if(edit) return maxspan;
  return 0;
}







/*************************************************************************
 *
 * Performs edits on bases in reads which lack a HAF2-7 (or HAF3-7) tag
 *
 * editmode:
 *  0 : very conservative, only columns with a single discrepancy
 *  1 : conservative, only columns where non-consensus bases occur once
 *  2 : slightly less conservative: edits where consensus base 8 times more
 *      frequent than non-cons base
 *  3 : like 2, but also edits bases with count <=2 and a HAF2 tag which
 *      is directly bordered by HAF3-7 tags
 *      NOT SUITED FOR EST/RNASEQ assemblies!!!
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint32 Contig::editSingleDiscrepancyNoHAFTag(vector<bool> & readsmarkedsrm, uint8 editmode)
{
  FUNCSTART("uint32 Contig::editSingleDiscrepancyNoHAFTag(vector<bool> & readsmarkedsrm)");

  if(editmode>3) editmode=3;

  cout << "\nSearching for single discr. without HAF tags:\n";

  uint32 numedits=0;

  // see whether we have some HAF tags at all (if not -> no edit!)
  {
    bool foundhaf=false;
    for(auto pcrI=CON_reads.begin(); pcrI!=CON_reads.end() && !foundhaf; ++pcrI){
      if(pcrI.getORPID() < 0
	 || pcrI->isBackbone()
	 || pcrI->isRail()) continue;

      foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF3);
      if(!foundhaf) foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF5);
      if(!foundhaf) foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF4);
      if(!foundhaf) foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF6);
      if(!foundhaf) foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF7);
      if(!foundhaf) foundhaf=pcrI->hasTag(Read::REA_tagentry_idHAF2);
    }
    if(!foundhaf){
      cout << "No read with HAF tags found.\n";
      FUNCEND();
      return 0;
    }
  }

  definalise();

  bool needcheckSRM=false;
  if(readsmarkedsrm.size() != CON_readpool->size()){
    // if we were called with an incorrectly sized readsmarkedsrm vector
    //  (empty or whatever), assume no SRM freshly set
    readsmarkedsrm.clear();
    readsmarkedsrm.resize(CON_readpool->size(),false);
  }else{
    for(uint32 i=0; i<readsmarkedsrm.size(); ++i){
      if(readsmarkedsrm[i]){
	needcheckSRM=true;
	break;
      }
    }
  }

  // initialise the iterator for getting through the contig
  rcci_t rcci(this);
  {
    vector<int32> allowedstrainids; // empty would be all ids
    vector<uint8> allowedreadtypes;
    //allowedreadtypes.push_back(ReadGroupLib::SEQTYPE_454GS20);
    rcci.init(allowedstrainids,
	      allowedreadtypes,
	      false,            // don't take rails
	      false,           // nor backbones
	      false);   // nor reads without readpool-reads
  }


  vector<uint8> allowededits(6);
  vector<simplebasecounter_t> basecounter(6);
  {
    static const string bases="ACGT*N";
    for(uint32 i=0;i<basecounter.size(); i++){
      basecounter[i].base=bases[i];
      basecounter[i].counter=0;
    }
  }

  ProgressIndicator<int32> P(0, static_cast<int32>(CON_counts.size())-1);
  auto ccI=CON_counts.begin();

  for(uint32 actcontigpos=0; actcontigpos<CON_counts.size(); ++actcontigpos, ++ccI, rcci.advance()){

    P.progress(actcontigpos);

    // must be at least at coverage 5
    if(ccI->total_cov<5) continue;

    // check for disagreement in this column
    uint32 setcount=0;

    if(ccI->A > 0) ++setcount;
    if(ccI->C > 0) ++setcount;
    if(ccI->G > 0) ++setcount;
    if(ccI->T > 0) ++setcount;
    if(ccI->N > 0) ++setcount;
    if(ccI->star > 0) ++setcount;

    // if more than 1, then there are bases in this column which disagree
    //  (caution: this might be also just an "N", or some other IUPACs)
    if(setcount>1) {
      if(needcheckSRM){
	// have a look at all reads at this position: if one of them was
	//  marked in this pass with a SRM, we don't do the edit here
	//  (because some "corrections" might be due to misassemblies and
	//  would be plain wrong)
	bool readfreshlymarkedsrm=false;
	for(auto & tpcrI : rcci.getPCRIsInCol()){
	  if(tpcrI.getORPID() >=0 && readsmarkedsrm[tpcrI.getORPID()]){
	    readfreshlymarkedsrm=true;
	    break;
	  }
	}
	if(readfreshlymarkedsrm) {
	  continue;
	}
      }


      CEBUG("cpos: " << actcontigpos << "\t" << ccI->A << " " << ccI->C << " " << ccI->G << " " << ccI->T << " " << ccI->N << " " << ccI->star << endl);

      for(auto & bce : basecounter){
	bce.counter=0;
      }
      for(auto & tpcrI : rcci.getPCRIsInCol()){
	char base;
	int32 readpos=tpcrI.contigPos2UnclippedReadPos(actcontigpos);

	//CEBUG("read: " << tpcrI->getName() << '\n');

	//int32 realreadpos;
	if(tpcrI.getReadDirection() > 0){
	  base=static_cast<char>(toupper(tpcrI->getBaseInSequence(readpos)));
	  //realreadpos=readpos;
	}else{
	  base=static_cast<char>(toupper(tpcrI->getBaseInComplementSequence(readpos)));
	  //realreadpos=tpcrI->calcComplPos(readpos);
	}

	for(auto & bce : basecounter){
	  if(bce.base==base) ++bce.counter;
	}
      }

      uint32 maxcount=0;
      for(uint32 bci=0; bci<basecounter.size(); ++bci){
	maxcount=max(maxcount,basecounter[bci].counter);
	CEBUG("bc["<<bci<<"]: " << basecounter[bci].base << " " << basecounter[bci].counter << endl);
      }

      if(maxcount<3) continue;

      for(auto & aee : allowededits) aee=0;

      // maxi = index of consensus base if >=0
      int32 maxi=-1;

      uint32 nddb=0; // num different discrepancy bases
      for(uint32 bci=0; bci<basecounter.size(); ++bci){
	if(basecounter[bci].counter==maxcount) {
	  if(maxi>=0){
	    maxi=-1;
	    break;
	  }
	  maxi=bci;
	}
	if(editmode <=1){
	  if(basecounter[bci].counter==1) {
	    allowededits[bci]=1;
	    ++nddb;
	  }
	}else if(editmode==2){
	  if(basecounter[bci].counter>0 && basecounter[bci].counter*8 < maxcount) {
	    allowededits[bci]=1;
	    ++nddb;
	  }
	}else if(editmode==3){
	  if(basecounter[bci].counter>0 && basecounter[bci].counter<=2 && basecounter[bci].counter*8 < maxcount) {
	    allowededits[bci]=1;
	    ++nddb;
	  }
	}else{
	  MIRANOTIFY(Notify::INTERNAL,"unknown editmode " << static_cast<uint16>(editmode));
	}
      }

      CEBUG("mc: " << maxcount << "\tmi: " << maxi << "\tnddb: " << nddb << '\t');
      for(uint32 bci=0;bci<basecounter.size(); ++bci){
	CEBUG(basecounter[bci].base << static_cast<uint16>(allowededits[bci]) << ' ');
      }
      CEBUG(endl);

      if(maxi<0
	 || (editmode==0 && nddb>1)) continue;

      char repbase=static_cast<char>(tolower(basecounter[maxi].base));

      for(auto & tpcrI : rcci.getPCRIsInCol()){
	if(tpcrI->isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)) continue;

	int32 readpos=tpcrI.contigPos2UnclippedReadPos(actcontigpos);

	CEBUG("read: " << tpcrI->getName());

	char base;
	int32 realreadpos;
	if(tpcrI.getReadDirection() > 0){
	  base=static_cast<char>(toupper(tpcrI->getBaseInSequence(readpos)));
	  realreadpos=readpos;
	}else{
	  base=static_cast<char>(toupper(tpcrI->getBaseInComplementSequence(readpos)));
	  realreadpos=tpcrI->calcComplPos(readpos);
	}

	CEBUG('\t' << base);

	// not optimal regarding speed, but let's not care atm
	bool docheck=false;
	for(uint32 bci=0;bci<basecounter.size(); ++bci){
	  if(basecounter[bci].base==base){
	    if(allowededits[bci]) docheck=true;
	  }
	}

	if(docheck){
	  CEBUG("\tcheck");

	  bool doedit=false;
	  // not optimal regarding speed, but let's not care atm
	  bool foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF3,realreadpos);
	  if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF4,realreadpos);
	  if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF5,realreadpos);
	  if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF6,realreadpos);
	  if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF7,realreadpos);

	  bool foundtag2=tpcrI->hasTag(Read::REA_tagentry_idHAF2,realreadpos);

	  if(!(foundtag||foundtag2)){
	    // no HAF tag at all
	    doedit=true;
	  }else if(editmode==3 && !foundtag && foundtag2){
	    // HAF2 tag, check left and right for HAF3-7
	    CEBUG("\tcheckh2");
	    bool h2ok=false;
	    if(realreadpos>0){
	      auto rrp2=realreadpos-1;
	      foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF3,rrp2);
	      if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF4,rrp2);
	      if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF5,rrp2);
	      if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF6,rrp2);
	      if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF7,rrp2);

	      rrp2+=2;
	      if(foundtag && rrp2 < tpcrI->getLenSeq()){
		foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF3,rrp2);
		if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF4,rrp2);
		if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF5,rrp2);
		if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF6,rrp2);
		if(!foundtag) foundtag=tpcrI->hasTag(Read::REA_tagentry_idHAF7,rrp2);

		if(foundtag){
		  doedit=true;
		}
	      }
	    }
	  }
	  if(doedit){
	    CEBUG("\tdoedit");
	    if(tpcrI.getReadDirection() > 0){
	      const_cast<Read &>(*tpcrI).changeBaseInSequence(repbase,0,realreadpos);
	    }else{
	      const_cast<Read &>(*tpcrI).changeBaseInSequence(dptools::getComplementIUPACBase(repbase),0,realreadpos);
	    }
	    // replacing with a gap? technically, it's a deletion followed by an insertion
	    // therefore, the adjustment must be ... errr ... adjusted
	    if(repbase=='*' && tpcrI->usesAdjustments()){
	      const_cast<Read &>(*tpcrI).changeAdjustment(realreadpos,-1);
	    }
	    addTagToConsensus(actcontigpos, actcontigpos,'=',"ESDN","",true);
	    ++numedits;
	  }
	}
	CEBUG(endl);
      }
    }
  }

  P.finishAtOnce();

  if(numedits) rebuildConCounts();

  FUNCEND();
  return numedits;
}
//#define CEBUG(bla)






#define CEBUGF2(bla)  {cout << bla; cout.flush();}
//#define CEBUGF2(bla)

// looks very much like newMarkPossibleRepeats()


void Contig::editPBSledgeHammer(vector<bool> & readsmarkedsrm, uint32 & numcoledits, uint32 & numreadedits)
{
  FUNCSTART("void Contig::editPBSledgeHammer(vector<bool> & readsmarkedsrm)");

  BUGIFTHROW(true,"need redo ce2 for PlacedContigReads");

  definalise();

////++++////  CEBUGF2("SledgeHammer" << endl;);
////++++////
////++++////  numcoledits=0;
////++++////  numreadedits=0;
////++++////
////++++////  bool needcheckSRM=false;
////++++////  if(readsmarkedsrm.size() != CON_readpool->size()){
////++++////    // if we were called with an incorrectly sized readsmarkedsrm vector
////++++////    //  (empty or whatever), assume no SRM freshly set
////++++////    readsmarkedsrm.clear();
////++++////    readsmarkedsrm.resize(CON_readpool->size(),false);
////++++////  }else{
////++++////    for(uint32 i=0; i<readsmarkedsrm.size(); ++i){
////++++////      if(readsmarkedsrm[i]){
////++++////	needcheckSRM=true;
////++++////	break;
////++++////      }
////++++////    }
////++++////  }
////++++////
////++++////  checkContig();
////++++////
////++++////  // this is "the truth"
////++++////  string concons;
////++++////  {
////++++////    vector<base_quality_t> dummy;
////++++////    newConsensusGet(concons,dummy,-1);
////++++////  }
////++++////
////++++////  // get number of strains in this contig
////++++////  uint32 numstrains=0;
////++++////  for(uint32 i=0; i<CON_reads.size(); i++){
////++++////    if(static_cast<uint32>(CON_reads[i].read.getStrainID()) > numstrains){
////++++////      numstrains=CON_reads[i].read.getStrainID();
////++++////    }
////++++////  }
////++++////  // remember: default strain has ID =0, so we have +1 strains
////++++////  numstrains++;
////++++////
////++++////  //contig_parameters const & con_params= CON_miraparams->getContigParams();
////++++////
////++++////  static const string groupbases("ACGT*");
////++++////
////++++////  nngroups_t emptygroup;
////++++////  emptygroup.base='!';
////++++////  emptygroup.valid=false;
////++++////  emptygroup.forwarddircounter=0;
////++++////  emptygroup.complementdircounter=0;
////++++////  emptygroup.groupquality=0;
////++++////
////++++////  // groups per seqtype per strain, empty template
////++++////  vector<vector<vector <nngroups_t> > > emptygroups_st_st;
////++++////
////++++////  emptygroups_st_st.resize(ReadGroupLib::getNumSequencingTypes());
////++++////  for(uint32 seqtype=0; seqtype<ReadGroupLib::getNumSequencingTypes(); seqtype++){
////++++////    emptygroups_st_st[seqtype].resize(numstrains);
////++++////    for(uint32 strainid=0; strainid<numstrains; strainid++){
////++++////      for(uint32 actgroup=0; actgroup<groupbases.size(); actgroup++){
////++++////	emptygroups_st_st[seqtype][strainid].push_back(emptygroup);
////++++////	emptygroups_st_st[seqtype][strainid].back().base=groupbases[actgroup];
////++++////      }
////++++////    }
////++++////  }
////++++////
////++++////  // groups per seqtype per strain, the real thing
////++++////  vector<vector<vector <nngroups_t> > > groups_st_st;
////++++////
////++++////
////++++////  nnpos_rep_col_t emptyprc;
////++++////  emptyprc.ids.clear();
////++++////  emptyprc.groupbases.clear();
////++++////  emptyprc.groupquals.clear();
////++++////  emptyprc.type=Read::REA_tagentry_idEmpty;
////++++////  emptyprc.contigpos=0;
////++++////  emptyprc.is_dangerous=false;
////++++////  emptyprc.tagged=false;
////++++////
////++++////
////++++////
////++++////  vector<int8> maskshadow;
////++++////  vector<multitag_t::mte_id_t> masktagtypes;
////++++////  buildMaskShadow(maskshadow,masktagtypes,false);
////++++////
////++++////  CEBUGF2("Start." << endl);
////++++////
////++++////  // the ercci is for the analysis of the bases themselves
////++++////  //
////++++////  ercci_t ercci;
////++++////  extendedReadColContigIteratorInit(ercci,
////++++////				    false,        // don't take rails
////++++////				    true,        // take backbone
////++++////				    numstrains);
////++++////
////++++////  // initialise the iterator for getting through the contig
////++++////  rcci_t rcci;
////++++////  {
////++++////    vector<int32> allowedstrainids; // empty would be all ids
////++++////    vector<uint8> allowedreadtypes;
////++++////    allowedreadtypes.push_back(ReadGroupLib::SEQTYPE_PACBIOLQ);
////++++////    readColContigIteratorInit(rcci,
////++++////			      allowedstrainids,
////++++////			      allowedreadtypes,
////++++////			      false,            // don't take rails
////++++////			      false,           // nor backbones
////++++////			      false);   // nor reads without readpool-reads
////++++////  }
////++++////
////++++////
////++++////
////++++////  // evil hack:
////++++////  // we use nmpr_firstfillin (or maybe csbrm_fillin_groups_stst()?) and
////++++////  //  nmpr_rategroups, but need some slightly different parameters for PacBio.
////++++////  // Therefore, save the original values and swap them
////++++////  //  back at the end of the function
////++++////
////++++////  contig_parameters orig_pbl_con_params = (*CON_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getContigParams();
////++++////  contig_parameters & nc_pbl_con_params = const_cast<contig_parameters &>((*CON_miraparams)[ReadGroupLib::SEQTYPE_PACBIOLQ].getContigParams());
////++++////  // for nmpr_firstfillin
////++++////  nc_pbl_con_params.con_endreadmarkexclusionarea=10;
////++++////  nc_pbl_con_params.con_minrmbneighbourqual=0;
////++++////  // for nmpr_rategroups
////++++////  nc_pbl_con_params.con_mingroupqualforrmbtagging=10;
////++++////  nc_pbl_con_params.con_minreadspergroup=4;
////++++////
////++++////  ProgressIndicator<int32> P(0, static_cast<int32>(CON_counts.size()));
////++++////  auto ccI=CON_counts.begin();
////++++////  for(uint32 actcontigpos=0; actcontigpos<CON_counts.size() ;actcontigpos++, ccI++, extendedReadColContigIteratorAdvance(ercci), readColContigIteratorAdvance(rcci)){
////++++////    P.progress(actcontigpos);
////++++////    //if(actcontigpos>6000) continue;
////++++////    CEBUGF2("acp: " << actcontigpos << endl);
////++++////
////++++////    // check for disagreement in this column
////++++////    if((ccI->A > 0)+(ccI->C > 0)+(ccI->G > 0)+(ccI->T > 0)+(ccI->N > 0)+(ccI->star > 0) <= 1) continue;
////++++////
////++++////    // ok, there are some disagreements
////++++////    CEBUGF2("Disagreement pos " << actcontigpos << ' ' << *ccI << endl);
////++++////
////++++////    if(needcheckSRM){
////++++////      // have a look at all reads at this position: if one of them was
////++++////      //  marked in this pass with a SRM, we don't do the edit here
////++++////      //  (because some "corrections" might be due to misassemblies and
////++++////      //  would be plain wrong)
////++++////      vector<int32>::const_iterator I=rcci.read_ids_in_col.begin();
////++++////      bool readfreshlymarkedsrm=false;
////++++////      for(;I!=rcci.read_ids_in_col.end();I++){
////++++////	if(readsmarkedsrm[CON_reads[*I].orpid]){
////++++////	  readfreshlymarkedsrm=true;
////++++////	  break;
////++++////	}
////++++////      }
////++++////      if(readfreshlymarkedsrm) {
////++++////	continue;
////++++////      }
////++++////    }
////++++////
////++++////
////++++////    // TODO: continue here
////++++////    groups_st_st=emptygroups_st_st;
////++++////
////++++////    // put the bases of the different reads into groups
////++++////    nmpr_firstfillin(ercci, maskshadow, masktagtypes, groups_st_st);
////++++////    nmpr_rategroups(groups_st_st, ccI);
////++++////
////++++////    // nmpr_rategroups will always set gap base groups to non valid
////++++////    //  therefore, need to "correct that" for PacBioLQ
////++++////    //
////++++////    // at the same time, count how many different groups are set altogether
////++++////    uint32 numvalidgroups=0;
////++++////    vector<bool> validgroupmask(groupbases.size(),false);
////++++////    {
////++++////      for(uint32 seqtype=0; seqtype<ReadGroupLib::getNumSequencingTypes(); seqtype++){
////++++////	for(uint32 strainid=0; strainid<numstrains; strainid++){
////++++////	  CEBUGF2("seqt: " << seqtype << "\tstrid: " << strainid << '\n');
////++++////	  for(uint32 actgroup=0; actgroup<groupbases.size(); actgroup++){
////++++////	    if(seqtype==ReadGroupLib::SEQTYPE_PACBIOLQ){
////++++////	      CEBUGF2("Group " << actgroup << '\n');
////++++////	      CEBUGF2(groups_st_st[seqtype][strainid][actgroup]);
////++++////	    }
////++++////	    if(seqtype==ReadGroupLib::SEQTYPE_PACBIOLQ && actgroup==4){
////++++////	      if(groups_st_st[seqtype][strainid][actgroup].forwarddircounter>2
////++++////		 && groups_st_st[seqtype][strainid][actgroup].complementdircounter>2
////++++////		 && groups_st_st[seqtype][strainid][actgroup].groupquality>=nc_pbl_con_params.con_mingroupqualforrmbtagging
////++++////		 && static_cast<int32>(groups_st_st[seqtype][strainid][actgroup].ids.size()) >= nc_pbl_con_params.con_minreadspergroup){
////++++////	        groups_st_st[seqtype][strainid][actgroup].valid=true;
////++++////	      }
////++++////	    }
////++++////	    if(groups_st_st[seqtype][strainid][actgroup].valid){
////++++////	      CEBUGF2("Valid possible group " << actgroup << '\n');
////++++////	      CEBUGF2(groups_st_st[seqtype][strainid][actgroup]);
////++++////	      if(!validgroupmask[actgroup]){
////++++////		validgroupmask[actgroup]=true;
////++++////		numvalidgroups++;
////++++////	      }
////++++////	    }
////++++////	  }
////++++////	}
////++++////      }
////++++////    }
////++++////
////++++////    CEBUGF2("numvalidgroups: " << numvalidgroups <<endl);
////++++////
////++++////    // OK, make it only for clear cases
////++++////    if(numvalidgroups==1){
////++++////      char validbase='?';
////++++////      for(uint32 vgbcounter=0; vgbcounter<validgroupmask.size(); ++vgbcounter){
////++++////	if(validgroupmask[vgbcounter]){
////++++////	  validbase=groupbases[vgbcounter];
////++++////	  break;
////++++////	}
////++++////      }
////++++////      char consbase=static_cast<char>(toupper(concons[actcontigpos]));
////++++////      CEBUGF2("vb: " << validbase << "\tcb: " << consbase <<endl);
////++++////      if(consbase==validbase){
////++++////	// Good, good. Consensus and only valid base group agree
////++++////	// That is ... edit away everything else (within boundaries)
////++++////
////++++////	char lvalidbase=static_cast<char>(tolower(validbase));
////++++////	bool hasreadedit=false;
////++++////	for(uint32 strainid=0; strainid<groups_st_st[0].size(); strainid++){
////++++////	  for(uint32 readnr=0; readnr < rcci.read_ids_in_col.size(); readnr++){
////++++////	    int32    actreadid =rcci.read_ids_in_col[readnr];
////++++////	    if(CON_reads[actreadid].read.isSequencingType(ReadGroupLib::SEQTYPE_PACBIOLQ)){
////++++////	      contigread_t & ric =CON_reads[actreadid];
////++++////	      char          base =toupper(ric.getBase(rcci.actcontigpos));
////++++////	      if(base!=validbase){
////++++////		//ric.setBase(rcci.actcontigpos,'n',0);
////++++////		ric.setBase(rcci.actcontigpos,lvalidbase,0);
////++++////		hasreadedit=true;
////++++////		++numreadedits;
////++++////	      }
////++++////	    }
////++++////	  }
////++++////	}
////++++////	if(hasreadedit) ++numcoledits;
////++++////      }
////++++////    }
////++++////  }
////++++////
////++++////  // put back the original PacBio LQ contig params
////++++////  nc_pbl_con_params=orig_pbl_con_params;
////++++////
////++++////if(numreadedits) rebuildConCounts();

  FUNCEND();
}

// check ec55989_70815+ near start   aaagccgcgXcac


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::blindContig()
{
  FUNCSTART("void Contig::blindContig()");

  definalise();

  auto cI=CON_reads.begin();
  for(; cI != CON_reads.end(); ++cI){
    if(cI.getReadDirection() > 0){
      const_cast<Read &>(*cI).blindSeqData('c');
    }else{
      const_cast<Read &>(*cI).blindSeqData('g');
    }
  }

  FUNCEND();
}
/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::upDownCase(base_quality_t threshold)
{
  FUNCSTART("void Contig::upDownCase(base_quality_t threshold)");

  definalise();

  for(auto & pcre : CON_reads){
    const_cast<Read &>(pcre).upDownCase(threshold);
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)

//#define WIGGLETRACK
//void Contig::reduceReadsAtCoveragePeaks(ccctype_t avgcov, vector<uint8> & peakindicator, deque<readid_t> & readsremoved)
void Contig::reduceReadsAtCoveragePeaks(ccctype_t avgcov, vector<uint8> & peakindicator, unordered_set<readid_t> & readsremoved)
{
  FUNCSTART("void Contig::reduceReadsAtCoveragePeaks(ccctype_t avgcov, vector<uint8> & peakindicator, unordered_set<readid_t> & readsremoved)");

  CEBUG("rRATCO: " << avgcov << endl);

  BUGIFTHROW(peakindicator.size()!=CON_counts.size(),"peakindicator.size()!=CON_counts.size()");
  readsremoved.clear();

#ifdef WIGGLETRACK
  {
    ccctype_t maxc=0;
    for(auto & cce : CON_counts) maxc=max(maxc,cce.total_cov);
    ofstream ofs;
    dbgOpenWiggle(ofs,"cov.wig",getContigName(),"cov",maxc);
    for(auto & cce : CON_counts) ofs << cce.total_cov << '\n';
  }
#endif

  // First, make a copy of the total coverage, we'll work on that copy to track
  //  coverage when first virtually removing reads
  vector<ccctype_t> virtcoverage(CON_counts.size());
  {
    auto vI=virtcoverage.begin();
    for(auto & cce : CON_counts){
      *vI=cce.total_cov;
      ++vI;
    }
  }

  // collect single mates in peaks which *should* have a mate in this contig
  // virtually remove them down to a coverage of 100% of avg coverage
  //  do not remove those which would make the coverage fall below that value
  priv_rratcp_collectReadsToDelete(avgcov,peakindicator,virtcoverage,false,false,readsremoved);
#ifdef WIGGLETRACK
  dbgContainerToWiggle(virtcoverage,getContigName(),"vcovms");
#endif

  // collect paired mates in peaks of this contig
  // virtually remove them down to a coverage of 90% of avg coverage
  //  do not remove those which would make the coverage fall below that value
  priv_rratcp_collectReadsToDelete((avgcov*9)/10,peakindicator,virtcoverage,false,true,readsremoved);
#ifdef WIGGLETRACK
  dbgContainerToWiggle(virtcoverage,getContigName(),"vcovmp");
#endif

  // now collect unpaired reads which contribute to overcoverage (of still existing)
  // but only down to 66% of avg coverage
  priv_rratcp_collectReadsToDelete((avgcov*2)/3,peakindicator,virtcoverage,true,false,readsremoved);
#ifdef WIGGLETRACK
  dbgContainerToWiggle(virtcoverage,getContigName(),"vcovsr");
#endif

  // if reads were virtually removed, remove them now from the contig
  if(!readsremoved.empty()){
    uint32 numr=0;
    auto pcrI=CON_reads.begin();
    while(pcrI!=CON_reads.end()){
      if(readsremoved.find(pcrI.getORPID()) != readsremoved.end()) {
	pcrI=deleteRead(pcrI);
	++numr;
      }else{
	++pcrI;
      }
    }
    CEBUG("Removed " << numr << " reads, " << readsremoved.size() << " planned.\n");
  }

  //    check for contig without coverage -> error!
  {
    for(auto ccI=CON_counts.begin(); ccI!=CON_counts.end(); ++ccI){
      BUGIFTHROW(ccI->total_cov==0,"Ouch: found 0 coverage in " << getContigName() << " at position " << ccI-CON_counts.begin());
    }
  }
}
#undef CEBUG
#undef WIGGLETRACK

//#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUG(bla)
void Contig::priv_rratcp_collectReadsToDelete(ccctype_t minthresh, vector<uint8> & peakindicator, vector<ccctype_t> & virtcoverage, bool remunpaired, bool pairmustbeincontig, unordered_set<readid_t> & readsremoved)
{
  FUNCSTART("void Contig::priv_rratcp_helper(ccctype_t minthresh, vector<uint8> & peakindicator, vector<ccctype_t> & virtcoverage, bool checktprid, unordered_set<readid_t> & readsremoved)");

  BUGIFTHROW(remunpaired&&pairmustbeincontig,"remunpaired&&pairmustbeincontig ???");

  // careful with minthresh==1 as we might encounter the following situation: overlapping
  //  paired-end at a coverage of 2. These would be both removed and then ... bang.
  // Two possible solutions:
  //  1) force minthreshold of 2
  //  2) make priv_rratcp_checkDeletable() check what it does (small time penalty)
  // let's choose 1) here
  if(minthresh<2) return;  // coverage too low anyway

  CEBUG("Collecting at minthresh " << minthresh << "\t" << remunpaired << "\t" << pairmustbeincontig << endl);

  auto pcrI=CON_reads.begin();
  auto opcrI=pcrI;
  auto crE=CON_reads.end();
  for(; pcrI != crE; ++pcrI){
    CEBUG("Looking at " << pcrI.getORPID() << "\t" << pcrI->getName() << endl);
    opcrI=CON_reads.getIteratorOfReadpoolID(pcrI->getTemplatePartnerID());
    bool tpartnerincontig=opcrI!=crE;
    CEBUG(pcrI->isBackbone() << '\t' << pcrI->isRail() << '\t' << pcrI->isCoverageEquivalentRead() << '\t' << remunpaired << '\t' << pairmustbeincontig << '\t' << pcrI->hasTemplateInfo()  << '\t' << pcrI->getTemplatePartnerID() << '\t' << tpartnerincontig << endl);
    if(!pcrI->isBackbone()
       && !pcrI->isRail()
       && !pcrI->isCoverageEquivalentRead()){
      bool checkthis=false;
      bool checktp=false;
      if(!pcrI->hasTemplateInfo()){
	if(remunpaired) checkthis=true;
      }else{
	if(pairmustbeincontig){
	  if(tpartnerincontig){
	    checkthis=true;
	    checktp=true;
	  }
	}else{
	  if(!tpartnerincontig && shouldHaveTPartnerInContig(pcrI,opcrI)){
	    checkthis=true;
	  }
	}
      }

      bool delthis=false;
      if(checkthis && priv_rratcp_checkDeletable(pcrI,minthresh,peakindicator,virtcoverage)){
	delthis=true;
      }
      bool deltp=false;
      if(delthis && checktp){
	if(opcrI!=CON_reads.end() && priv_rratcp_checkDeletable(opcrI,minthresh,peakindicator,virtcoverage)){
	  deltp=true;
	}else{
	  delthis=false;
	}
      }
      // make sure that mates are not deleted twice
      // pairs of the same length with the same start site will
      //  also not be deleted, but when will this happen on real data? Almost never,
      //  so don't care
      if(delthis && deltp){
	if(pcrI.getReadStartOffset()>opcrI.getReadStartOffset()
	   || (pcrI.getReadStartOffset()==opcrI.getReadStartOffset()
	       && pcrI->getLenClippedSeq() >= opcrI->getLenClippedSeq())){
	  delthis=false;
	  deltp=false;
	}
      }

      // check whether previously removed
      if(delthis && readsremoved.find(pcrI.getORPID())!=readsremoved.end()) delthis=false;
      if(deltp && readsremoved.find(opcrI.getORPID())!=readsremoved.end()) deltp=false;

      if(delthis){
	CEBUG("PR removing " << pcrI->getName() << '\t' << pcrI.getReadStartOffset() << endl);
	auto vcI=virtcoverage.begin()+pcrI.getReadStartOffset();
	auto vcE=vcI+pcrI->getLenClippedSeq();
	for(; vcI!=vcE; ++vcI) *vcI-=1;;
	readsremoved.insert(pcrI.getORPID());
	if(deltp){
	  CEBUG("PR removing TP " << opcrI->getName() << '\t' << opcrI.getReadStartOffset() << endl);
	  vcI=virtcoverage.begin()+opcrI.getReadStartOffset();
	  vcE=vcI+opcrI->getLenClippedSeq();
	  for(; vcI!=vcE; ++vcI) *vcI-=1;;
	  readsremoved.insert(opcrI.getORPID());
	}
      }
    }
  }
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUG(bla)
bool Contig::priv_rratcp_checkDeletable(PlacedContigReads::const_iterator pcrI, ccctype_t minthresh, vector<uint8> & peakindicator, vector<ccctype_t> & virtcoverage)
{
  FUNCSTART("bool Contig::priv_rratcp_checkDeletable(PlacedContigReads::const_iterator pcrI, ccctype_t minthresh, vector<uint8> & peakindicator, vector<ccctype_t> & virtcoverage)");

  // careful here because of overlapping paired-end at coverage 2
  // instead of complicated checking, just make sure we bail out if we encounter
  //  a situation which should not be (and should have been handled earlier)
  BUGIFTHROW(minthresh<2,"minthresh<2 ???");

  CEBUG("Checking " << pcrI.getORPID() << "\t" << pcrI->getName() << endl);
  {
    auto piI=peakindicator.begin()+pcrI.getReadStartOffset();
    auto piE=piI+pcrI->getLenClippedSeq();
    for(; piI!=piE && *piI; ++piI);
    // not completely covered by peak? forget it, next read
    if(piI != piE){
      CEBUG("No peak coverage " << pcrI->getName() << '\t' << pcrI.getReadStartOffset() << endl);
      return false;
    }
  }
  {
    auto vcI=virtcoverage.begin()+pcrI.getReadStartOffset();
    auto vcE=vcI+pcrI->getLenClippedSeq();
    for(; vcI!=vcE && *vcI>minthresh; ++vcI);
    // coverage somewhere below threshold? forget it, next read
    if(vcI != vcE) {
      CEBUG("Would like to remove, but coverage drop " << pcrI->getName() << '\t' << pcrI.getReadStartOffset() << endl);
      return false;
    }
  }
  return true;
}
#undef CEBUG
