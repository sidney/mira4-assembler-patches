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

#include "contig.H"


using namespace std;


#define CEBUG(bla)

//#define PARANOIABUGTRACKFLAG
#ifdef PARANOIABUGTRACKFLAG
#define paranoiaBUGSTAT(statement) { statement;}
#define paranoiaBUGIF(ifcond, statement) { if(ifcond) {statement;}}
#else
#define paranoiaBUGSTAT(statement)
#define paranoiaBUGIF(ifcond, statement)
#endif




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla) {if(CON_cebugflag) cout << bla;}


// new ... calculate main consensus and cache it
// consensus for strains: moved as "on-demand" calulation to getConsensus()

//#define CEBUG(bla) {cout << bla;}
void Contig::calcConsensi(int32 mincoverage, base_quality_t minqual, char missingcoveragechar)
{
  FUNCSTART("void Contig::calcConsensi(int32 mincoverage, base_quality_t minqual, char missingcoveragechar)");

  //CON_cebugflag=true;

  CEBUG("calcConsensi(). mincov: " << mincoverage << "\tminqual: " << static_cast<uint16>(minqual) << "\tmissingcchar: " << missingcoveragechar << endl);

  // loading from different files (backbones reads etc) make the number of strain change over time
  // i.e., contigs loaded earlier may have a smaller CON_readsperstrain.size() than is good
  if(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains()){
    CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);
  }

  CEBUG("Rebuild for " << CON_readsperstrain.size() << " strain in readpool\n");

  if(CON_conscalc_mincov!=mincoverage
     || CON_conscalc_minqual!=minqual
     || CON_conscalc_missingchar!=missingcoveragechar
     || CON_strainconsseq.size()==0
     || CON_readsperstrain.size() > CON_strainconsseq.size()){    // BaCh 17.09.2013: bugfix ">=" into ">", else always needless recalcs

    CEBUG("need recalc.:\n");
    CEBUG("old values: " << CON_conscalc_mincov << " " << static_cast<uint16>(CON_conscalc_minqual) << " " << CON_conscalc_missingchar << " " << CON_strainconsseq.size() << " " << CON_readsperstrain.size() << endl);
    CEBUG("new values: " << mincoverage << " " << static_cast<uint16>(minqual) << " " << missingcoveragechar << " " << CON_strainconsseq.size() << " " << CON_readsperstrain.size() << endl);

    CON_conscalc_mincov=mincoverage;
    CON_conscalc_minqual=minqual;
    CON_conscalc_missingchar=missingcoveragechar;

    // make a consensus for every strain
    // cache them in CON_allconsseq and *qual
    makeIntelligentConsensus(CON_allconsseq,
			     CON_allconsqual,
			     &CON_alladjustments,
			     nullptr,
			     0,
			     CON_counts.size(),
			     0,
			     0,
			     -1,
			     missingcoveragechar);

    CEBUG("CON_readsperstrain.size(): " << CON_readsperstrain.size() << endl);

    CON_strainconsseq.clear();
    CON_strainconsqual.clear();
    CON_strainadjustments.clear();
    CON_strainconsseq.resize(CON_readsperstrain.size());
    CON_strainconsqual.resize(CON_readsperstrain.size());
    CON_strainadjustments.resize(CON_readsperstrain.size());

    uint32 numstrains=0;
    for(uint32 si=0; si<CON_readsperstrain.size(); si++){
      if(CON_readsperstrain[si]>0) numstrains++;
    }

    // now for all strains.
    // strains not present in contig will have sequence, quality and adjustments
    //  pre-filled with default values (@, 0 and -1)
    // strains present will have empty seq+qual+adj ... to be calculated on demand

    // Change: no precalculated strains, only pre-filled for clear cases!
    // Should getConsensus() ask for them later, they
    //  will be calculated on demand (using CON_conscalc_* values)
    for(uint32 si=0; si<CON_readsperstrain.size(); si++){
      if(CON_readsperstrain[si]==0) {
	CON_strainconsseq[si].resize(CON_allconsseq.size(),missingcoveragechar);
	CON_strainconsqual[si].resize(CON_allconsqual.size(),0);
	CON_strainadjustments[si].resize(CON_allconsqual.size(),-1);
      }else{
	if(numstrains==1 && mincoverage==0 && minqual==0){
	  // take over from allstrains in this very special case
	  CON_strainconsseq[si]=CON_allconsseq;
	  CON_strainconsqual[si]=CON_allconsqual;
	  CON_strainadjustments[si]=CON_alladjustments;
	}
      }
    }
  }else{
  }

  FUNCEND();

  return;
}
//#define CEBUG(bla)


// makes sure consensus and all adjoining structures are calculated and valid
// works like newConsensusGet() below, but does not return the consensus nor qualities
//
// trick: calling newConsensusGet() with the src==target will make sure things are not
//  unnecessarily copied
void Contig::ensureConsensus(int32 strainidtotake)
{
  FUNCSTART("void Contig::ensureConsensus(int32 strainidtotake)");

  CEBUG("ensureConsensus(): " << strainidtotake << endl);
  if(CON_allconsseq.empty() || CON_strainconsseq.empty()){
    CEBUG("ensureConsensus(): cons empty, need recalc" << endl);
    calcConsensi();
  }
  if(strainidtotake<0){
    newConsensusGet(CON_allconsseq, CON_allconsqual, strainidtotake);
  }else{
    newConsensusGet(CON_strainconsseq[strainidtotake], CON_strainconsqual[strainidtotake], strainidtotake);
  }
  FUNCEND();
}

void Contig::newConsensusGet(string & target, vector<base_quality_t> & qual, int32 strainidtotake)
{
  FUNCSTART("void Contig::newConsensusGet(string & target, vector<base_quality_t> & qual, int32 strainidtotake)");

//  if(CON_abortflag) {
//    uint16 * bombme=nullptr;
//    *bombme=0xdead;
//    cout << bombme;
//  }

  CEBUG("newConsensusGet(): gimme strain " << strainidtotake << endl);

  // loading from different files (backbones reads etc) make the number of strain change over time
  // i.e., contigs loaded earlier may have a smaller CON_readsperstrain.size() than is good
  if(CON_readsperstrain.size() < ReadGroupLib::getNumOfStrains()){
    CON_readsperstrain.resize(ReadGroupLib::getNumOfStrains(),0);
  }

  BUGIFTHROW(strainidtotake>=static_cast<int32>(CON_readsperstrain.size()),"strainidtotake>=CON_readsperstrain.size() ?");

  if(CON_fixedconsseq.size() && CON_fixedconsqual.size() && strainidtotake<0){
    CEBUG("newConsensusGet(): get fixed" << endl);
    target=CON_fixedconsseq;
    qual=CON_fixedconsqual;
  }else{
    CEBUG("CON_allconsseq.size(): " << CON_allconsseq.size() << endl);
    CEBUG("CON_strainconsseq.size(): " << CON_strainconsseq.size() << endl);
    if(CON_allconsseq.empty() || CON_strainconsseq.empty()
       || (strainidtotake>=0 && strainidtotake>=CON_strainconsseq.size())){
      CEBUG("something's empty, need recalc" << endl);
      calcConsensi();
    }
    if(strainidtotake<0){
      CEBUG("newConsensusGet(): get allcons" << endl);
      CEBUG("CON_allconsseq.size(): " << CON_allconsseq.size() << endl);
      CEBUG("CON_allconsqual.size(): " << CON_allconsqual.size() << endl);
      target=CON_allconsseq;
      qual=CON_allconsqual;
    }else{
      CEBUG("newConsensusGet(): get strain " << strainidtotake << endl);
      BUGIFTHROW(strainidtotake>=CON_strainconsseq.size(),"something's utterly wrong: strainidtotake>=CON_strainconsseq.size() ???");
      // on demand calculation
      if(CON_strainconsseq[strainidtotake].empty()){
	CEBUG("check on demand calculation\n");
	uint32 numstrains=0;
	for(uint32 si=0; si<CON_readsperstrain.size(); si++){
	  if(CON_readsperstrain[si]>0) numstrains++;
	}
	if(numstrains==1){
	  CEBUG("only 1 strain, can take main consensus\n");
	  CON_strainconsseq[strainidtotake]=CON_allconsseq;
	  CON_strainconsqual[strainidtotake]=CON_allconsqual;
	  CON_strainadjustments[strainidtotake]=CON_alladjustments;
	}else{
	  CEBUG("must do calculation\n");
	  makeIntelligentConsensus(CON_strainconsseq[strainidtotake],
				   CON_strainconsqual[strainidtotake],
				   &CON_strainadjustments[strainidtotake],
				   nullptr,
				   0,
				   CON_counts.size(),
				   CON_conscalc_mincov,
				   CON_conscalc_minqual,
				   strainidtotake,
				   CON_conscalc_missingchar);
	}
      }else{
	CEBUG("take cached\n");
      }
      target=CON_strainconsseq[strainidtotake];
      qual=CON_strainconsqual[strainidtotake];
    }
  }

  FUNCEND();
}
//#define CEBUG(bla)


//void Contig::OLDgetConsensus1(string & target, vector<base_quality_t> & qual, bool markspecials, int32 mincoverage, base_quality_t minqual, int32 strainidtotake, char missingcoveragechar, ostream * ostr, bool contagsintcs)
//{
//  FUNCSTART("void Contig::getConsensus(string & target, vector<base_quality_t> & qual, bool markspecials, int32 mincoverage, int32 strainidtotake, ostream * ostr, bool contagsintcs)");
//
//  CEBUG("getCons()\n");
//
//  CON_cebugflag=true;
//
//  if(CON_cheat_intelcons.empty()
//     || CON_cheat_intelcons_markspecials!=markspecials
//     || CON_cheat_intelcons_mincov!=mincoverage
//     || CON_cheat_intelcons_minqual!=minqual
//     || CON_cheat_intelcons_strainidtotake!=strainidtotake
//     || ostr != nullptr) {
//    bool mustcompute=true;
//    if(!CON_cheat_intelcons.empty()
//       && CON_cheat_intelcons_mincov==mincoverage
//       && CON_cheat_intelcons_minqual==minqual
//       && CON_cheat_intelcons_strainidtotake == strainidtotake) {
//      if(markspecials) {
//	if(CON_cheat_intelcons_markspecials) {
//	  mustcompute=false;
//	}
//      } else {
//	mustcompute=false;
//	if(CON_cheat_intelcons_markspecials) {
//	  // just make the sequence make uppercase
//	  for(uint32 i=0; i<CON_cheat_intelcons.size(); i++) {
//	    CON_cheat_intelcons[i]=toupper(CON_cheat_intelcons[i]);
//	  }
//	}
//      }
//    }
//    if(mustcompute || ostr != nullptr) {
//      makeIntelligentConsensus(CON_cheat_intelcons,
//			       CON_cheat_intelconsqual,
//			       0,
//			       CON_counts.size(),
//			       markspecials,
//			       mincoverage,
//			       minqual,
//			       strainidtotake,
//			       missingcoveragechar,
//			       ostr,
//			       contagsintcs);
//    }
//    CON_cheat_intelcons_markspecials=markspecials;
//    CON_cheat_intelcons_mincov=mincoverage;
//    CON_cheat_intelcons_minqual=minqual;
//    CON_cheat_intelcons_strainidtotake=strainidtotake;
//  }
//
//  target=CON_cheat_intelcons;
//  qual=CON_cheat_intelconsqual;
//
//  FUNCEND();
//
//  return;
//}






/*************************************************************************
 *
 * if the routines to decide for a base (helper2 routines) could not get
 *  clear base but the user wants one, this is doing a shootout based on
 *  majority vote.
 *
 * groups with forward and reverse count double in the read count
 *
 * if majority vote still does not work (all the same), then the last
 *  one wins. I.e., * takes precedence over T, this over G, over C, over A
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper3(char & thisbase, base_quality_t & thisqual, const vector<nngroups_t> & groups, const vector<char> & IUPACbasegroups)
{
  uint32 maxcount=0;

  for(uint32 actgroup=0; actgroup<groups.size(); actgroup++){
    for(uint32 actbase=0; actbase<IUPACbasegroups.size(); actbase++){
      if(IUPACbasegroups[actbase]==groups[actgroup].base){
	uint32 groupcount=groups[actgroup].urdids.size();
	if(groups[actgroup].forwarddircounter>0
	   && groups[actgroup].complementdircounter>0){
	  groupcount*=2;
	}
	if(groupcount>=maxcount){
	  maxcount=groupcount;
	  thisbase=groups[actgroup].base;
	  thisqual=groups[actgroup].groupquality;
	}
      }
    }
  }

}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
//#define FUNCSTART(bla)  static const char * THISFUNC = bla"  "; {cout << "enter " << THISFUNC << "\n"; cout.flush();}
//#define FUNCEND() {cout << "exit " << THISFUNC << "\n"; cout.flush();}


//#define CEBUG(bla) {cout << bla;}
void Contig::makeIntelligentConsensus_helper2_calcSOLEXA(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups, const base_quality_t maxqual, const uint32 maxcount, int32 strainidtotake)
{
  FUNCSTART("void Contig::makeIntelligentConsensus_helper2_calcSOLEXA(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups, const base_quality_t maxqual, const uint32 maxcount, int32 strainidtotake)");

  //makeIntelligentConsensus_helper3(thisbase,
  //				   thisqual,
  //				   groups,
  //				   IUPACbasegroups);
  //return;

  IUPACbasegroups.clear();


  // TODO: clarify whether to fill IUPACbasegroups also with merged base??

  bool hasmergedbases=(ccI->getOriginalBBChar()!='@') & (ccI->bbcounts[0]>0);
  bool mergedalsogroup=false;
  bool mergedisvalid=false;

  CEBUG("bbchar: " << ccI->getOriginalBBChar() << "\tbbcounts: " << ccI->bbcounts[0] << "\tbbbestquals: " << static_cast<uint16>(ccI->bbbestquals[0]) << "\tbbstrains: " << hex << static_cast<uint16>(ccI->bbstrains[0]) << dec << "\tHasmergedb: " << hasmergedbases <<'\n');

  // Ok, check whether the merged bases belong to this strain
  //  if not, well then no merged bases exist
  uint8 strainmask=255;

  CEBUG("Strainidtotake: " << strainidtotake << '\n');
  if(strainidtotake>=0) strainmask=getBBStrainMask(strainidtotake);
  if(!(ccI->bbstrains[0] & strainmask)) hasmergedbases=false;

  uint32 numgroups=0;
  base_quality_t lastseengroupqual=0;
  char lastseenbase='&';

  bool groupschosen[groups.size()];

  // search for excellent groups
  for(uint32 i=0; i<groups.size(); i++){
    groupschosen[i]=false;
    CEBUG("Lookatgroup: " << groups[i] << '\n');
    if(hasmergedbases && ccI->getOriginalBBChar() == groups[i].base){
      if((groups[i].urdids.size() + ccI->bbcounts[0]) >= 10
	 && (ccI->bbbestquals[0] >=35
	     || (groups[i].groupquality >40
		 && groups[i].forwarddircounter>0
		 && groups[i].complementdircounter>0))){
	// >= 6 solexa reads with qual >=30
	CEBUG("Excellent1: " << groups[i].base << '\n');
	numgroups++;
	groupschosen[i]=true;
	lastseengroupqual=groups[i].groupquality;
	lastseenbase=groups[i].base;
	mergedisvalid=true;
	if(!groups[i].urdids.empty()) mergedalsogroup=true;
      }
    }else{
      // or >=10 normal reads, gqual>=40 and forward/reverse directions
      //  are at least 25% of all reads
      if(groups[i].urdids.size()>=10
	 && groups[i].groupquality >=40
	 && groups[i].forwarddircounter>=groups[i].urdids.size()/4
	 && groups[i].complementdircounter>=groups[i].urdids.size()/4) {
	CEBUG("Excellent2: " << groups[i].base << '\n');
	numgroups++;
	groupschosen[i]=true;
	lastseengroupqual=groups[i].groupquality;
	lastseenbase=groups[i].base;
      }
    }
  }
  CEBUG("Numgroups excel: " << numgroups << '\n');
  CEBUG("Hasmergedbases: " << hasmergedbases << "\tMergedalsogroup: " << mergedalsogroup << '\n');


  // search for high groups
  if(numgroups==0){
    for(uint32 i=0; i<groups.size(); i++){
      CEBUG("Lookatgroup: " << groups[i] << '\n');
      if(hasmergedbases && ccI->getOriginalBBChar() == groups[i].base){
	if((groups[i].urdids.size() + ccI->bbcounts[0]) >= 6
	   && (ccI->bbbestquals[0] >=30
	       || (groups[i].groupquality >=30
		   && groups[i].forwarddircounter>0
		   && groups[i].complementdircounter>0))){
	  // >= 6 solexa reads with qual >=30
	  CEBUG("High1: " << groups[i].base << '\n');
	  numgroups++;
	  groupschosen[i]=true;
	  lastseengroupqual=groups[i].groupquality;
	  lastseenbase=groups[i].base;
	  mergedisvalid=true;
	  if(!groups[i].urdids.empty()) mergedalsogroup=true;
	}
      }else{
	// or >=10 normal reads, gqual>=30 and forward/reverse
	if(groups[i].urdids.size()>=6
	   && groups[i].groupquality >=30
	   && groups[i].forwarddircounter>0
	   && groups[i].complementdircounter>0) {
	  CEBUG("High2: " << groups[i].base << '\n');
	  numgroups++;
	  groupschosen[i]=true;
	  lastseengroupqual=groups[i].groupquality;
	  lastseenbase=groups[i].base;
	}
      }
    }
    CEBUG("Numgroups high: " << numgroups << '\n');
    CEBUG("Hasmergedbases: " << hasmergedbases << "\tMergedalsogroup: " << mergedalsogroup << '\n');
  }

  // TODO: insert good. fwd/rev and qual vs non-fwd/rev and qual

  // search for good groups. Same as high, but without mincount
  if(numgroups==0){
    for(uint32 i=0; i<groups.size(); i++){
      CEBUG("Lookatgroup: " << groups[i] << '\n');
      if(hasmergedbases && ccI->getOriginalBBChar() == groups[i].base){
	if(ccI->bbbestquals[0] >=30
	    || (groups[i].groupquality >=30
		&& groups[i].forwarddircounter>0
		&& groups[i].complementdircounter>0)){
	  // >= 6 solexa reads with qual >=30
	  CEBUG("Good1: " << groups[i].base << '\n');
	  numgroups++;
	  groupschosen[i]=true;
	  lastseengroupqual=groups[i].groupquality;
	  lastseenbase=groups[i].base;
	  mergedisvalid=true;
	  if(!groups[i].urdids.empty()) mergedalsogroup=true;
	}
      }else{
	// or gqual>=30 and forward/reverse
	if(groups[i].groupquality >=30
	   && groups[i].forwarddircounter>0
	   && groups[i].complementdircounter>0) {
	  CEBUG("Good2: " << groups[i].base << '\n');
	  numgroups++;
	  groupschosen[i]=true;
	  lastseengroupqual=groups[i].groupquality;
	  lastseenbase=groups[i].base;
	}
      }
    }
    CEBUG("Numgroups good: " << numgroups << '\n');
    CEBUG("Hasmergedbases: " << hasmergedbases << "\tMergedalsogroup: " << mergedalsogroup << '\n');
  }

  // search for medium groups
  if(numgroups==0){
    for(uint32 i=0; i<groups.size(); i++){
      if(!groups[i].urdids.empty() && groups[i].groupquality >=30) {
	numgroups++;
	groupschosen[i]=true;
	lastseengroupqual=groups[i].groupquality;
	lastseenbase=groups[i].base;
	if(hasmergedbases
	   && ccI->getOriginalBBChar() == groups[i].base
	   && ccI->bbbestquals[0] >=30) {
	  mergedisvalid=true;
	  mergedalsogroup=true;
	}
      }else{
	if(hasmergedbases && ccI->getOriginalBBChar() == groups[i].base
	   && ccI->bbbestquals[0] >=30){
	  mergedisvalid=true;
	  numgroups++;
	  groupschosen[i]=true;
	}
      }
    }
    CEBUG("Numgroups medium: " << numgroups << '\n');
    CEBUG("Hasmergedbases: " << hasmergedbases << "\tMergedalsogroup: " << mergedalsogroup << '\n');
  }

  bool badquals=false;
  if(numgroups==0){
    // no read in the groups? still might be ... but with bad quals.
    CEBUG("No group first time, redo without qual.\n");

    for(uint32 i=0; i<groups.size(); i++){
      if(!groups[i].urdids.empty()) {
	numgroups++;
	groupschosen[i]=true;
	lastseengroupqual=groups[i].groupquality;
	lastseenbase=groups[i].base;
	if(hasmergedbases
	   && ccI->getOriginalBBChar() == groups[i].base) {
	  mergedisvalid=true;
	  mergedalsogroup=true;
	}
      }else{
	if(hasmergedbases && ccI->getOriginalBBChar() == groups[i].base) {
	  mergedisvalid=true;
	  numgroups++;
	  groupschosen[i]=true;
	}
      }
    }
    badquals=true;

    CEBUG("Numgroups: " << numgroups << "\tbadquals: " << badquals << '\n');
    CEBUG("Hasmergedbases: " << hasmergedbases << "\tMergedalsogroup: " << mergedalsogroup << '\n');
  }

  char maxcount_base='%';
  base_quality_t maxcount_qual=0;

  if(numgroups==1){
    // just one group
    // see whether it is a valid merged base
    if(hasmergedbases && mergedisvalid) {
      // if merged, is it also covered by a readgroup?
      if(mergedalsogroup){
	// yes
	thisbase=ccI->getOriginalBBChar();
	thisqual=max(lastseengroupqual, ccI->bbbestquals[0]);
      }else{
	// no, only a merged base
	thisbase=ccI->getOriginalBBChar();
	thisqual=ccI->bbbestquals[0];
      }
    }else{
      // no, it's only one normal read group
      thisbase=lastseenbase;
      IUPACbasegroups.push_back(thisbase);
      thisqual=lastseengroupqual;
    }
  } else if(numgroups==0){
    thisbase='N';
    thisqual=0;
  } else {
    // let the fun begin ... *sigh*
    // at least two groups were detected, which means any combination of
    //  normal read groups (1-5) and eventually also a merged base
    // the groups were either detected all at good qual or all at bad qual

    // first decision taken on counts of bases
    vector<uint32> counts(groups.size());
    vector<base_quality_t> quals(groups.size());
    uint32 maxcounts=0;
    base_quality_t maxquals=0;      // carefull, this is not maxqual
    // from the parameter (rename that)
    uint32 totalcounts=0;
    bool hasgapchosen=false;

    for(uint32 i=0; i<groups.size(); i++){
      counts[i]=groups[i].urdids.size();
      quals[i]=groups[i].groupquality;
      if(groupschosen[i] && groups[i].base == '*') hasgapchosen=true;
      if(hasmergedbases
	 && ccI->getOriginalBBChar() == groups[i].base) {
	counts[i]+=ccI->bbcounts[0];
	quals[i]=max(quals[i],ccI->bbbestquals[0]);
      }
      totalcounts+=counts[i];
      maxcounts=max(maxcounts,counts[i]);
      maxquals=max(maxquals,quals[i]);
      if(maxcounts==counts[i]){
	maxcount_base=groups[i].base;
	maxcount_qual=quals[i];
      }
    }

    CEBUG("Maxcounts: " << maxcounts << "\ttotalcounts: " << totalcounts << '\n');

    thisbase=' ';
    thisqual=0;
    uint32 numtaken=0;
    for(uint32 i=0; i<groups.size(); i++){
      if(groupschosen[i]){
	numtaken++;
	if(thisbase==' ') {
	  thisbase=groups[i].base;
	  thisqual=quals[i];
	}else{
	  // take different approaches depending on whether user forces
	  //  a non-IUPAC or not
	  //  or whether a gap is part of the best groups
	  if(counts[i] == maxcounts
	     && (hasgapchosen
		 || (*CON_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getContigParams().con_force_nonIUPACconsensus_perseqtype)) {
	    thisbase=groups[i].base;
	    thisqual=quals[i];
	    numtaken=1;
	  }else{
	    thisbase=dptools::calcIUPACConsensus(thisbase, groups[i].base);
	    thisqual+=quals[i];
	  }
	}
	if(!groups[i].urdids.empty()) IUPACbasegroups.push_back(thisbase);
      }
    }
    if(numtaken) thisqual/=numtaken;
  }

  CEBUG("Settled on: '" << thisbase << "' " << static_cast<uint16>(thisqual) << '\n');

  if(!dptools::isValidACGTStarBase(thisbase)
     && maxcount_base!='%'
     && (*CON_miraparams)[ReadGroupLib::SEQTYPE_SOLEXA].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
    thisbase=maxcount_base;
    thisqual=maxcount_qual;
    CEBUG("User forced non-IUPAC per seqtype: '" << thisbase << "' " << static_cast<uint16>(thisqual) << '\n');
  }

  FUNCEND();
  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcPACBIOHQ(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups)
{

  // no info. atm, same as 454 without 40:60 rule

  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?
  // TODO: write PacBioLQ/HQ routine

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_PACBIOHQ].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  return;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calc454GS20(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups)
{
  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  // if max count is gap, but there are other bases
  if(maxsize_i==4 && runnerup>0){
    // apply 40:60 rule
    if(100*runnerup/(maxsize+runnerup) >= 40){
      swap(maxsize,runnerup);
      swap(maxsize_i,runnerup_i);
      swap(maxqual,runnerupqual);
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  //// TODO: testing
  //if(maxsize_i>=0 && runnerup_i>=0){
  //  if(maxsize+runnerup == 0){
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454","DOH?");
  //  }else if(100*runnerup/(maxsize+runnerup) >= 30){
  //    //ostringstream ostr;
  //    //ostr << static_cast<char>(groups[maxsize_i].base) << ": " << maxsize;
  //    //ostr << " " << static_cast<char>(groups[runnerup_i].base) << ": " << runnerup;
  //    //ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
  //    //
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454",ostr.str().c_str());
  //  }
  //}

  return;
}



/*************************************************************************
 *
 * atm a pure copy of 454
 *
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcIonTorrent(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups)
{
  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

  // if max count is gap, but there are other bases
  // Test, maybe IonTorrent is a bit different.
  // Actually, it's almost the best rule there is ...
  //  35 perhaps better, but needs to be verified.
  // Astonishing
  if(maxsize_i==4 && runnerup>0){
    // apply 40:60 rule
    if(100*runnerup/(maxsize+runnerup) >= 40){
      swap(maxsize,runnerup);
      swap(maxsize_i,runnerup_i);
      swap(maxqual,runnerupqual);
    }
  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_IONTORRENT].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  //// TODO: testing
  //if(maxsize_i>=0 && runnerup_i>=0){
  //  if(maxsize+runnerup == 0){
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454","DOH?");
  //  }else if(100*runnerup/(maxsize+runnerup) >= 30){
  //    //ostringstream ostr;
  //    //ostr << static_cast<char>(groups[maxsize_i].base) << ": " << maxsize;
  //    //ostr << " " << static_cast<char>(groups[runnerup_i].base) << ": " << runnerup;
  //    //ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
  //    //
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454",ostr.str().c_str());
  //  }
  //}

  return;
}




/*************************************************************************
 *
 * atm a pure copy of 454 & IonTorrent, except the gap override rule
 * is disabled
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcText(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups)
{
  // Idea: pure coverage, maximum wins.
  //       if two or more with same maximum:
  //          IUPAC (star goes under, sorry) or
  //          if non-IUPAC is wished, take the last
  // TODO: is there a better way?

  IUPACbasegroups.clear();

  int32 maxsize=0;
  int32 maxsize_i=-1;
  base_quality_t maxqual=0;
  int32 runnerup=0;
  int32 runnerup_i=-1;
  base_quality_t runnerupqual=0;

  size_t totalsize=0;

  for(uint32 i=0; i<groups.size(); i++){
    totalsize+=groups[i].urdids.size();
    if(static_cast<int32>(groups[i].urdids.size())>=maxsize){
      runnerup=maxsize;
      runnerup_i=maxsize_i;
      runnerupqual=maxqual;
      maxsize=static_cast<int32>(groups[i].urdids.size());
      maxsize_i=i;
      maxqual=groups[i].groupquality;
    }else if(static_cast<int32>(groups[i].urdids.size())>=runnerup){
      runnerup=static_cast<int32>(groups[i].urdids.size());
      runnerup_i=i;
      runnerupqual=groups[i].groupquality;
    }
  }

//  // if max count is gap, but there are other bases
//  // Test, maybe IonTorrent is a bit different.
//  // Actually, it's almost the best rule there is ...
//  //  35 perhaps better, but needs to be verified.
//  // Astonishing
//  if(maxsize_i==4 && runnerup>0){
//    // apply 40:60 rule
//    if(100*runnerup/(maxsize+runnerup) >= 40){
//      swap(maxsize,runnerup);
//      swap(maxsize_i,runnerup_i);
//      swap(maxqual,runnerupqual);
//    }
//  }

  int32 avgqual=0;
  for(uint32 i=0; i<groups.size(); i++){
    if(maxsize==static_cast<int32>(groups[i].urdids.size())){
      IUPACbasegroups.push_back(groups[i].base);
      avgqual+=groups[i].groupquality;
    }
  }

  if(totalsize==0 || IUPACbasegroups.empty()){
    /// Oooops? just an N???
    thisbase='N';
    thisqual=0;
    return;
  }

  if(IUPACbasegroups.size()==1) {
    thisbase=IUPACbasegroups[0];
    // reduce quality if there are doubts
    if(runnerup>0 && maxsize-runnerup < 10){
      avgqual-=runnerupqual;
      if(avgqual<0) avgqual=max(abs(avgqual),10);
    }
    thisqual=static_cast<base_quality_t>(avgqual);
  }else{
    if((*CON_miraparams)[ReadGroupLib::SEQTYPE_454GS20].getContigParams().con_force_nonIUPACconsensus_perseqtype) {
      thisbase=IUPACbasegroups[IUPACbasegroups.size()-1];
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }else{
      thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
      thisqual=static_cast<base_quality_t>(avgqual/IUPACbasegroups.size());
    }
  }

  //// TODO: testing
  //if(maxsize_i>=0 && runnerup_i>=0){
  //  if(maxsize+runnerup == 0){
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454","DOH?");
  //  }else if(100*runnerup/(maxsize+runnerup) >= 30){
  //    //ostringstream ostr;
  //    //ostr << static_cast<char>(groups[maxsize_i].base) << ": " << maxsize;
  //    //ostr << " " << static_cast<char>(groups[runnerup_i].base) << ": " << runnerup;
  //    //ostr << "  -  " << 100*runnerup/(maxsize+runnerup) << "%";
  //    //
  //    //addTagToConsensus(actcontigpos, actcontigpos,'=',"P454",ostr.str().c_str());
  //  }
  //}

  return;
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Contig::makeIntelligentConsensus_helper2_calcSangerQual(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, const vector<nngroups_t> & groups, vector<char> & IUPACbasegroups, const vector<char> & columnbases, const base_quality_t maxqual, const uint32 maxcount)
{
  // Idea: all groups with qual >= 30
  //    or groups with qual within (X) (with X==5 now) of maxqual
  //    or groups with count >= 3 and quality >=20 and forward/reverse

  IUPACbasegroups.clear();
  int32 avgqual=0;
  int32 groupstaken=0;
  int32 verygoodgroups=0;

  // New: if force non-IUPAC consensus is wished
  // decrease the level with which we look at the groups
  // gradually until we have to find something (qual == 0)
  base_quality_t goodgroupqual=35;
  base_quality_t lessergroupqual=goodgroupqual;
  if(lessergroupqual>=10) lessergroupqual-=10;

  do{
    if(goodgroupqual>=5) goodgroupqual-=5;
    if(lessergroupqual>=5) lessergroupqual-=5;
    for(uint32 i=0; i<groups.size(); i++){
      if(groups[i].urdids.size()
	 &&(groups[i].groupquality >= goodgroupqual
	    || groups[i].groupquality+5 >= maxqual
	    || (groups[i].urdids.size() >= 3          // TODO: check 2 or 3 (was 3, testing with 2)
		&& groups[i].groupquality >= lessergroupqual
		&& groups[i].forwarddircounter>0
		&& groups[i].complementdircounter>0))) {
	avgqual+=groups[i].groupquality;
	IUPACbasegroups.push_back(groups[i].base);
	groupstaken++;

	if(groups[i].groupquality >= goodgroupqual
	   && groups[i].forwarddircounter>0
	   && groups[i].complementdircounter>0) {
	  verygoodgroups++;
	}
      }
    }
  } while((*CON_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getContigParams().con_force_nonIUPACconsensus_perseqtype
	  && groupstaken==0 && goodgroupqual>0);

  CEBUG("groups taken: " << groupstaken << "\t");
#ifdef CEBUGFLAG
  {
    for(uint32 i=0; i< IUPACbasegroups.size(); i++){
      CEBUG(IUPACbasegroups[i]);
    }
  }
#endif
  CEBUG(endl);
  if(groupstaken == 0) {
    // this can still happen if we do not force a non-IUPAC
    // or if it's a column entirely made of N (or X) (or both)

    // well, calculate the base from all columnbases
    thisbase=dptools::calcIUPACConsensus(columnbases);
    avgqual=0;
    int32 n=0;
    uint32 totalids=0;
    for(uint32 i=0; i<groups.size()-1; i++){
      if(groups[i].urdids.size()) {
	totalids+=groups[i].urdids.size();
	avgqual+=groups[i].groupquality;
	n++;
      }
    }
    // test if there are more stars than other bases
    if(groups[groups.size()-1].urdids.size() > totalids) {
      // oh well, more stars than bases, make it a star
      thisbase='*';
      thisqual=groups[groups.size()-1].groupquality;
    } else if(n) {
      thisqual=avgqual/n;
    } else {
      thisqual=0;
    }
  } else if(groupstaken == 1) {
    // only one group taken (perfect)
    thisbase=IUPACbasegroups.front();
    int32 bahqual=maxqual;
    for(uint32 i=0; i<groups.size(); i++){
      if(groups[i].base!=thisbase) {
	bahqual-=groups[i].groupquality;
      }
    }
    if(bahqual<0) bahqual=0;
    thisqual=bahqual;
  } else {
    // ouch, more than one group was taken as good

    bool needIUPACconsensus=true;

    CEBUG("cpos: " << actcontigpos << "\tpossible IUPAC\tvgg: " << verygoodgroups << endl);
    //for(uint32 i=0; i< groups.size(); i++){
    //  cout << groups[i] << endl;
    //}

    // First, check whether we have a very good group
    //  as searched for above
    if(verygoodgroups==1){
      // if yes, rebuild IUPACbasegroups with very good groups
      // then, if only one group remains, perfect. Else it'll be
      //  a IUPAC

      CEBUG("recheck group: ");

      vector<char> newIUPACbasegroups;
      int32 newavgqual=0;
      int32 newgroupstaken=0;
      for(uint32 i=0; i<groups.size(); i++){
	if(groups[i].urdids.size()
	   && (groups[i].groupquality >= goodgroupqual
	       && groups[i].forwarddircounter>0
	       && groups[i].complementdircounter>0)){
	  newavgqual+=groups[i].groupquality;
	  newIUPACbasegroups.push_back(groups[i].base);
	  newgroupstaken++;
	}
      }

      if(newgroupstaken == 1) {
	CEBUG("success!");
	needIUPACconsensus=false;

	IUPACbasegroups=newIUPACbasegroups;
	avgqual=newavgqual;
	groupstaken=newgroupstaken;

	// only one group taken (perfect)
	thisbase=IUPACbasegroups.front();
	int32 bahqual=maxqual;
	for(uint32 i=0; i<groups.size(); i++){
	  if(groups[i].base!=thisbase) {
	    bahqual-=groups[i].groupquality;
	  }
	}
	if(bahqual<0) bahqual=0;
	thisqual=bahqual;
      }
    }
    CEBUG(endl);

    if(needIUPACconsensus){
      if((*CON_miraparams)[ReadGroupLib::SEQTYPE_SANGER].getContigParams().con_force_nonIUPACconsensus_perseqtype){
	// well, can't get a good decision, but user wants a non-IUPAC
	// make a majority-vote shootout among the chosen
	//   groups
	makeIntelligentConsensus_helper3(thisbase,
					 thisqual,
					 groups,
					 IUPACbasegroups);
      }else{
	// we'll have to make
	//  a IUPAC consensus here. Alas, if a star (*) is part of this group,
	//  it will 'go under' (the good bases always win). So check if the star
	//  has max quality.
	bool isstar=false;
	if(groups[groups.size()-1].groupquality == maxqual) {
	  // yep, now check if one of the bases has the same qual
	  isstar=true;
	  for(uint32 i=0; i<groups.size()-1; i++){
	    if(groups[i].groupquality == maxqual) {
	      // in dubio pro reo: it could be a base
	      isstar=false;
	    }
	  }
	}
	CEBUG("isstar: " << isstar << endl);
	if(isstar) {
	  thisbase='*';
	  thisqual=groups[groups.size()-1].groupquality;
	} else {
	  // really several different good bases.
	  thisqual=avgqual/groupstaken;
	  thisbase=dptools::calcIUPACConsensus(IUPACbasegroups);
	}
      }
    }
  }
}



/*************************************************************************
 *
 * For a given readtype, return thisbase and thisqual as result for this
 *  position in the contig
 * Also return the potential number of solutions (bases), the one which
 *  was chosen plus the ones which were not
 *
 * all other arguments are passovers from the main function (saving STL
 *  setup and memory allocation time)
 *
 * Side effects: manifold. Most notably: all the passovers from the
 *  main function which are not const will be overwritten with values
 *  reflecting the calculation of base and quality of this read type
 *
 * E.g.: the bases that were considered to make the consensus are in
 *  IUPACbasegroups, the base-groups in groups etc.
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}

void Contig::makeIntelligentConsensus_helper1(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const int32 mincoverage, vector<nngroups_t> & groups, vector<nngroups_t> & maskedshadowgroups, vector<char> & IUPACbasegroups, vector<char> & columnbases, const vector<PlacedContigReads::const_iterator> & read_pcrIs_in_col, vector<int8> & maskshadow, uint8 actreadtype, int32 strainidtotake, char missingcoveragechar)
{
  FUNCSTART("void Contig::makeIntelligentConsensus_helper1(char & thisbase, base_quality_t & thisqual, const uint32 actcontigpos, cccontainer_t::const_iterator ccI, const int32 mincoverage, vector<nngroups_t> & groups, vector<nngroups_t> & maskedshadowgroups, vector<char> & IUPACbasegroups, vector<char> & columnbases, const vector<int32> & read_pcrIs_in_col, vector<int8> & maskshadow, uint8 actreadtype, char missingcoveragechar)");

#if 0
  if(read_pcrIs_in_col.empty()){
    thisbase='!';
    FUNCEND();
    return;
  }
  if(CON_counts[actcontigpos].total_cov < mincoverage){
    thisbase='N';
    FUNCEND();
    return;
  }

#else

  size_t nummapped=0;
  if(actreadtype == ReadGroupLib::SEQTYPE_SOLEXA){
    // Ok, check whether the merged bases belong to this strain
    //  if not, well then no merged bases exist
    if(ccI->bbcounts[0]){
      uint8 strainmask=255;
      if(strainidtotake>=0) strainmask=getBBStrainMask(strainidtotake);
      if(ccI->bbstrains[0] & strainmask){
	nummapped=ccI->bbcounts[0];
//      cout << "looked at: " << *ccI << endl
//	   << "nummapped: " << nummapped << endl;
	//nummapped=0;
      }
    }
  }
  if(read_pcrIs_in_col.size()+nummapped == 0){
    thisbase='!';
    FUNCEND();
    return;
  }
  if(read_pcrIs_in_col.size()+nummapped < mincoverage){
    thisbase='N';
    FUNCEND();
    return;
  }
#endif

  const contig_parameters & con_rt_params = (*CON_miraparams)[actreadtype].getContigParams();

  CEBUG("conpos: " << actcontigpos << "\tnumreads: " << read_pcrIs_in_col.size() << endl);
  columnbases.clear();

  /*
     slowwwwwwwwwwwwwwwwwwwwwww

     on a 170MB CAF file, caf2fasta goes from 1:46 to 1:40 when
     replacing these two lines with the loop below (from 49s for
     output down to 43s)

     groups=emptygroups;
     maskedshadowgroups=emptygroups;

  */

#ifdef CLOCK_STEPS_CONS
  timeval tv,tvtotal,tvsub;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  for(uint32 groupi=0; groupi < groups.size(); groupi++){
    groups[groupi].groupquality=0;
    groups[groupi].valid=false;
    groups[groupi].forwarddircounter=0;
    groups[groupi].complementdircounter=0;
    groups[groupi].urdids.clear();
    groups[groupi].quals.clear();
    groups[groupi].directions.clear();
    //groups[groupi].strainids.clear();
    //groups[groupi].urdids_sorted.clear();
    //groups[groupi].uniqstrainids.clear();
    //groups[groupi].countstrainids.clear();

    maskedshadowgroups[groupi].groupquality=0;
    maskedshadowgroups[groupi].valid=false;
    maskedshadowgroups[groupi].forwarddircounter=0;
    maskedshadowgroups[groupi].complementdircounter=0;
    maskedshadowgroups[groupi].urdids.clear();
    maskedshadowgroups[groupi].quals.clear();
    maskedshadowgroups[groupi].directions.clear();
    //maskedshadowgroups[groupi].strainids.clear();
    //maskedshadowgroups[groupi].urdids_sorted.clear();
    //maskedshadowgroups[groupi].uniqstrainids.clear();
    //maskedshadowgroups[groupi].countstrainids.clear();
  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_MGROUPS]+=diffsuseconds(tv);
  gettimeofday(&tv,nullptr);
#endif

  bool xonly=true;
  bool hasmaskedset=false;

  if(nummapped>0 && ccI->getOriginalBBChar() != 'X') xonly=false;

  for(auto & pcrI : read_pcrIs_in_col){
    char           base;
    base_quality_t qual;

    //if(ric.read.isShortReadMapping()) continue;

#ifdef CLOCK_STEPS_CONS
    // way too expensive for this loop!
    //gettimeofday(&tvsub,nullptr);
#endif
    int32 readpos=pcrI.contigPos2UnclippedReadPos(actcontigpos);
#ifdef CLOCK_STEPS_CONS
    //CON_us_steps_cons[USCLOCONS_H1_CP2URP]+=diffsuseconds(tvsub);
#endif

    CEBUG("cc pcrI: " << pcrI << "\treadpos " << readpos);

    // yup, looks like this is better than using "pcrI->" all the time
    //  it is at least faster
    // maybe too hard for GCC to optimise across whole function

    auto & actread=*pcrI;
    auto rdir=pcrI.getReadDirection();

    int32 realreadpos;
    if(rdir>0){
      base=toupper(actread.getBaseInSequence(readpos));
      qual=actread.getQualityInSequence(readpos);
      realreadpos=readpos;
    }else{
      base=toupper(actread.getBaseInComplementSequence(readpos));
      qual=actread.getQualityInComplementSequence(readpos);
      realreadpos=actread.calcComplPos(readpos);
    }

    // Bad idea, this is soooo slow compared to the above
    //pcrI.getBaseAndQuality(actcontigpos,base,qual);
    //base=toupper(base);
    //if(rdir<0){
    //  realreadpos=actread.calcComplPos(readpos);
    //}

    columnbases.push_back(base);

    // this one gives a problem
    //   if(base != 'X' && base != '*') xonly=false;
    // when one retrieves the consensus of a single strain in multiple
    //  strain assemblies: it may well be only a star!
    // preliminary fix: back to checking of X only

    if(base != 'X') xonly=false;

    CEBUG("\t" << base << "\t" << (uint16) qual << endl);

    // bases near start/end of reads might be dangerous because
    //  of possible vector leftovers or still bad quality
    // so, if there is more than one read in this column,
    // ...
    // - if within 10 bases of vector clip, max quality = 2*distance to svec
    // - if not, lower the quality values of the bases that are near ends
    //   (con_endreadmarkexclusionarea bases) of (eventually clipped) reads
    // must be regular read (i.e. not backbone, rail or shortreadmapping

    if(likely(!actread.isRail()
	      && !actread.isBackbone()
	      && !actread.isCoverageEquivalentRead())){
      if(likely(read_pcrIs_in_col.size() >1)) {
	if((actread.getLSClipoff() >0
	    && realreadpos < actread.getLSClipoff()+10)
	   || (actread.getRSClipoff() < actread.getLenSeq()
	       && realreadpos > actread.getRSClipoff()-10)){
	  CEBUG(actread.getName()<< ": near seq vec in read, lowering the quality.\n");
	  int32 distance;
	  if(realreadpos < actread.getLSClipoff()+10){
	    distance=realreadpos-actread.getLSClipoff();
	  }else{
	    distance=actread.getRSClipoff()-realreadpos;
	  }
	  if(qual > distance*2){
	    qual=distance*2;
	  }
	} else if(realreadpos <= actread.getLeftClipoff()+con_rt_params.con_endreadmarkexclusionarea
		  || realreadpos >= actread.getRightClipoff()-con_rt_params.con_endreadmarkexclusionarea) {
	  CEBUG(actread.getName()<< ": near end of read, lowering the quality.\n");

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

	  if(actread.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
	    // decrease by one
	    if(likely(qual)) --qual;
	  }else{
	    // all other atm
	    if(qual>=10) {
	      qual-=10;
	    } else {
	      qual=5;
	    }
	    if(unlikely(qual<5)) qual=5;
	  }
	}
      }
    }

    // No railreads at all, have been taken out earlier!
    //// Quality of bases from railreads are set to 0 so as not
    ////  to be counted twice (as bases are the same as in backbone)
    //if(actread.isRail()) qual=0;



    bool maskedset=false;


    // TODO: rework the if(...hasTag(Read::REA_tagFpas,realreadpos))
    //  to check for all mask taks in a masktagstring vector

//    if(maskshadow[actcontigpos]) {
//      // remember that the readpos computing routine doesn't take care
//      //  of direction, so we have to complement that position in reverse cases
//      int32 realreadpos=readpos;
//      if(rdir<0){
//	realreadpos=actread.calcComplPos(readpos);
//      }
//
//      CEBUG("MASKED: " << actcontigpos << "\t");
//      CEBUG(actread.getName() << "\t" << realreadpos << "\t");
//      if(actread.hasTag(Read::REA_tagFpas,realreadpos)) {
//	CEBUG("in" << endl);
//	for(uint32 bindex=0; bindex<maskedshadowgroups.size(); bindex++) {
//	  if(maskedshadowgroups[bindex].base==base) {
//	    maskedshadowgroups[bindex].urdids.push_back(actreadid);
//	    maskedshadowgroups[bindex].quals.push_back(qual);
//	    maskedshadowgroups[bindex].directions.push_back(rdir);
//	    if(rdir>0){
//	      maskedshadowgroups[bindex].hasforwarddir=true;
//	    }else{
//	      maskedshadowgroups[bindex].hascomplementdir=true;
//	    }
//	    // special case: treat short read mapping as both forward and reverse
//	    if(actread.isCoverageEquivalentRead()){
//	      maskedshadowgroups[bindex].hasforwarddir=true;
//	      maskedshadowgroups[bindex].hascomplementdir=true;
//	    }
//	    maskedset=true;
//	    hasmaskedset=true;
//	    break;
//	  }
//	}
//      }
//    }


    if(likely(maskedset==false)) {
      CEBUG("out" << endl);
      //for(uint32 bindex=0; bindex<groups.size(); bindex++) {
      for(auto & groupe : groups){
	if(groupe.base==base) {
	  groupe.urdids.push_back(pcrI.getURDID());
	  groupe.quals.push_back(qual);
	  groupe.directions.push_back(rdir);
	  if(rdir>0){
	    groupe.forwarddircounter++;
	  }else{
	    groupe.complementdircounter++;
	  }

	  // special case: treat short read mapping as both forward and reverse
	  if(actread.isCoverageEquivalentRead()){
	    // we already counted forward
	    //groupe.forwarddircounter++;
	    groupe.complementdircounter++;
	  }
	  break;
	}
      }
    }
  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_PCRI]+=diffsuseconds(tv);
  gettimeofday(&tv,nullptr);
#endif

  if(hasmaskedset){
    // check whether there is anything in groups
    //  and maskedshadowgroups selected
    // if groups empty and masked not, copy masked to groups
    bool groupsempty=true;
    bool maskedempty=true;
    for(uint32 bindex=0; bindex<groups.size(); ++bindex) {
      if( groups[bindex].urdids.size()) groupsempty=false;
      if( maskedshadowgroups[bindex].urdids.size()) maskedempty=false;
    }
    if(groupsempty && maskedempty==false) {
      CEBUG("EMPTY: " << actcontigpos << endl);
      groups=maskedshadowgroups;
    }
  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_EGROUP]+=diffsuseconds(tv);
  gettimeofday(&tv,nullptr);
#endif

  //{
  //  // Cheat!
  //  bool hasmergedbases=(ccI->getOriginalBBChar()!='@') & (ccI->bbcounts[0]>0);
  //  bool mergedalsogroup=false;
  //
  //  // Ok, check whether the merged bases belong to this strain
  //  //  if not, well then no merged bases exist
  //  uint8 strainmask=255;
  //
  //  if(strainidtotake>=0) strainmask=getBBStrainMask(strainidtotake);
  //  if(!(ccI->bbstrains[0] & strainmask)) hasmergedbases=false;
  //
  //
  //}





  if(xonly) {
    thisbase='X';
  } else {
    base_quality_t maxqual=0;
    uint32 maxcount=0;
    for(uint32 i=0; i<groups.size(); i++){
      calcGroupQual(groups[i]);
      if(groups[i].groupquality>maxqual){
	maxqual=groups[i].groupquality;
      }
      if(groups[i].urdids.size()>maxcount){
	maxcount=groups[i].urdids.size();
      }
      //if(actcontigpos>830 && actcontigpos <920) {
      CEBUG(actcontigpos << "\tb: " << groups[i].base << "\tgq: " << (uint16) groups[i].groupquality << "\ts: " << groups[i].urdids.size() << endl);
      //}
    }

#ifdef CLOCK_STEPS_CONS
    CON_us_steps_cons[USCLOCONS_H1_GQUAL]+=diffsuseconds(tv);
    gettimeofday(&tv,nullptr);
#endif

    //if(actcontigpos>830 && actcontigpos <920) {
    CEBUG(actcontigpos << "maxqual: " << (uint16) maxqual << "\t" << "maxcount: " << maxcount << endl);
    //}

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

    switch(actreadtype) {
    case ReadGroupLib::SEQTYPE_SANGER : {
      makeIntelligentConsensus_helper2_calcSangerQual(thisbase,
						      thisqual,
						      actcontigpos,
						      groups,
						      IUPACbasegroups,
						      columnbases,
						      maxqual,
						      maxcount);
      break;
    }
    case ReadGroupLib::SEQTYPE_454GS20 : {
      makeIntelligentConsensus_helper2_calc454GS20(thisbase,
						   thisqual,
						   actcontigpos,
						   groups,
						   IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_IONTORRENT : {
      makeIntelligentConsensus_helper2_calcIonTorrent(thisbase,
						      thisqual,
						      actcontigpos,
						      groups,
						      IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_PACBIOHQ :
    case ReadGroupLib::SEQTYPE_PACBIOLQ : {
      // TODO: eventually two different for hq and lq ?
      makeIntelligentConsensus_helper2_calcPACBIOHQ(thisbase,
						    thisqual,
						    actcontigpos,
						    groups,
						    IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_TEXT : {
      makeIntelligentConsensus_helper2_calcText(thisbase,
						thisqual,
						actcontigpos,
						groups,
						IUPACbasegroups);
      break;
    }
    case ReadGroupLib::SEQTYPE_SOLEXA : {
      // the static cast for the strain id
//TODO: hier weiter
      makeIntelligentConsensus_helper2_calcSOLEXA(
	thisbase,
	thisqual,
	actcontigpos,
	ccI,
	groups,
	IUPACbasegroups,
	maxqual,
	maxcount,
	strainidtotake
	);
      break;
    }
    case ReadGroupLib::SEQTYPE_ABISOLID : {
      cout << "Actreadtype: " << static_cast<uint16>(actreadtype) << endl;
      MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support.");
      break;
    }
    default: {
      cout << "Actreadtype: " << static_cast<uint16>(actreadtype) << endl;
      MIRANOTIFY(Notify::INTERNAL, "Unknown read type.");
    }
    }

#ifdef CLOCK_STEPS_CONS
    CON_us_steps_cons[USCLOCONS_H1_CALLH2]+=diffsuseconds(tv);
    gettimeofday(&tv,nullptr);
#endif

    if(!dptools::isValidIUPACStarBase(thisbase)){
      CEBUG("ALERT! This is not a valid base: '" << thisbase << "'\t: " << static_cast<uint16>(thisbase) << '\n');
    }

  }

#ifdef CLOCK_STEPS_CONS
  CON_us_steps_cons[USCLOCONS_H1_TOTAL]+=diffsuseconds(tvtotal);
#endif


  FUNCEND();
}






/*************************************************************************
 *
 * Calculate the 'true' consensus and gives back a string with consensus
 *  and a vector with the base quality of each base
 *
 * strainidtotake: <0 means "all", >=0 means "exactly those reads with that id"
 *
 * If the ostream parameter is != nullptr, also writes a .tcs live file
 *  to it
 *
 *************************************************************************/

//#define CEBUG(bla) {cout << bla;}

void Contig::makeIntelligentConsensus(string & target, vector<base_quality_t> & qual, vector<int32> * targetadjustments, vector<int32> * targetadjustments_bla, int32 from, int32 to, int32 mincoverage, base_quality_t minqual, int32 strainidtotake, char missingcoveragechar)//, ostream * ostr, bool contagsintcs)
{
  FUNCSTART("void Contig::makeIntelligentConsensus(string & target, vector<base_quality_t> & qual, int32 from, int32 to, int32 mincoverage, base_quality_t minqual, int32 strainidtotake)");//, ostream * ostr, bool contagsintcs)");

  cout << "makeIntelligentConsensus() complete calc .. "; cout.flush();

  //CON_cebugflag=true;

  CEBUG("MIC\n");
  CEBUG("from " << from << endl);
  CEBUG("to " << to << endl);
  CEBUG("mincoverage " << mincoverage << endl);
  CEBUG("minqual " << static_cast<uint16>(minqual) << endl);
  CEBUG("missingcovchar " << missingcoveragechar << endl);
  CEBUG("strainidtotake " << strainidtotake << endl);


  BUGIFTHROW(from>to,"from>to?");
  BUGIFTHROW(from<0, "from < 0 ?");

  suseconds_t mict_fin=0;
  suseconds_t mict_pre=0;
  suseconds_t mict_shadow=0;
  suseconds_t mict_fallout=0;
  suseconds_t mict_newin=0;
  suseconds_t mict_helper1=0;
  suseconds_t mict_restofloop=0;
  suseconds_t mict_totalloop=0;

  timeval us_start;
  gettimeofday(&us_start,nullptr);

  // Finalising the contig initialises the output order structure vector
  finalise();
  mict_fin+=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  //const contig_parameters & con_params = CON_miraparams->getContigParams();

  if( to > static_cast<int32>(CON_counts.size())) to=CON_counts.size();
  int32 len_target=to-from;

  //target.resize(len_target);
  target.clear();
  target.reserve(len_target+10);
  qual.clear();
  qual.reserve(len_target+10);

  // for calculating the adjustments, only do this when whole
  //  consensus is calculated
  if(targetadjustments != nullptr && from==0 && to==CON_counts.size()) {
    targetadjustments->clear();
    targetadjustments->reserve(len_target+10);
  }
  int32 unpaddedposcounter=0;


  nngroups_t emptygroup;
  emptygroup.base='!';
  emptygroup.valid=false;
  // TODO: make also use of bothdirectionspresent in this routine?
  emptygroup.forwarddircounter=0;
  emptygroup.complementdircounter=0;
  emptygroup.groupquality=0;
  vector<nngroups_t> emptygroups;
  for(uint32 i=0; i<5; i++) {
    emptygroups.push_back(emptygroup);
    emptygroups[i].base= "ACGT*"[i];
  }

/*
  // if there's a stream, we're dumping TCS
  // initialise a quick lookup vector to point at the positions in the
  //  consensus that have a tag
  // a for any tag
  // d for dangerous tag
  vector<bool> tcs_aconstagpositions;
  vector<bool> tcs_dconstagpositions;
  if(ostr != nullptr && contagsintcs){
    tcs_aconstagpositions.resize(CON_counts.size(),false);
    tcs_dconstagpositions.resize(CON_counts.size(),false);
    vector<consensustag_t>::const_iterator I=CON_consensus_tags.begin();
    for(; I!=CON_consensus_tags.end(); I++) {
      for(uint32 i=I->from; i<=I->to; i++) tcs_aconstagpositions[i]=true;
      if(I->identifier == CON_tagentry_idSRMc
	 || I->identifier == CON_tagentry_idWRMc
	 || I->identifier == CON_tagentry_idDGPc
	 || I->identifier == CON_tagentry_idUNSc
	 || I->identifier == CON_tagentry_idIUPc){
	for(uint32 i=I->from; i<=I->to; i++) tcs_dconstagpositions[i]=true;
      }
    }
  }
  // temporary vectors for TCS
  vector<int32> tcs_totalgroupcount;
  tcs_totalgroupcount.reserve(emptygroups.size());
  vector<int32> tcs_totalgroupqual;
  tcs_totalgroupqual.reserve(emptygroups.size());
*/


  vector<nngroups_t> maskedshadowgroups=emptygroups;

  vector<char> IUPACbasegroups;
  IUPACbasegroups.reserve(10);

  vector<char> columnbases;
  columnbases.reserve(1000);

  // for each read type, we will predict a base and quality and such
  vector<char> predictedbases(ReadGroupLib::getNumSequencingTypes());
  vector<base_quality_t> predictedquality(ReadGroupLib::getNumSequencingTypes());

  // as well as the possible alternative bases
  vector<vector<char> > possiblebases(ReadGroupLib::getNumSequencingTypes());

  // as well as a goodness level estimate of the prediction
  vector<int8> predictlevel(ReadGroupLib::getNumSequencingTypes());
  vector<int8> predictlevelused(10);


  // for each read type a vector<nngroups_t> to hold all estimates
  //  regarding bases and qualities
  vector<vector<nngroups_t> > groupsvec;
  groupsvec.reserve(ReadGroupLib::getNumSequencingTypes());

  // a vector (read_ids_in_col) keeps the ids of the reads that are
  //  covering a specific position of the contig
  // reserving 1000 position should be enough 99.9% of all cases,
  //  is automatically extended by STL if needed.

  // now with different read types (sanger, 454 etc), we need a vector
  //  of iterators to the PCR. Was formerly a vector of read-ids, but
  //  not possible anymore

  vector<vector<PlacedContigReads::const_iterator> > read_pcrIs_in_col;
  read_pcrIs_in_col.resize(ReadGroupLib::getNumSequencingTypes());

  // ok, fill in some
  for(auto & rpice : read_pcrIs_in_col){
    rpice.reserve(1000);
    groupsvec.push_back(emptygroups);
  }

  // fill the vector in case we are starting within the contig
  if(from>0) {
    // TODO: this is wrong when going only after specific strainids!
    // correct that ASAN (as soon as needed)

    MIRANOTIFY(Notify::INTERNAL, "starting with from > 0 not available yet.");

    //getReadIDsAtContigPosition(read_pcrIs_in_col,from);
  }


  mict_pre=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  vector<int8> maskshadow;
  vector<multitag_t::mte_id_t> masktagstrings;

  // Bach: 17.08.2008
  // the new strategy of tagging poly-AT sites and keeping them in
  //  the read (no clipping) makes it necessary to keep sequence
  //  under Fpas tags as full valid member of consensus somputation
  // Therefore, Fpas may NOT be put into the masktagstrings anymore!
  // masktagstrings.push_back(Read::REA_tagFpas);

  buildMaskShadow(maskshadow,masktagstrings,false);

  mict_shadow=diffsuseconds(us_start);
  gettimeofday(&us_start,nullptr);

  vector<uint8> cached_contigseqtypes(ReadGroupLib::getNumSequencingTypes(),0);
  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); ++actseqtype){
    cached_contigseqtypes[actseqtype]=hasSeqTypeData(actseqtype);
  }

  auto pcrI=CON_reads.begin();
  auto ccI=CON_counts.cbegin();
  advance(ccI,from);

  // this is the loop that updates the vector that
  //  keeps track of the reads that are
  //  covering a specific position of the contig
  // works quite simple: reads in the vector that are ending get
  //  thrown out, new reads starting are entered.

  CEBUG("CON_counts.size(): " << CON_counts.size() << endl);

  timeval us_loop;

  for(uint32 actcontigpos=from; actcontigpos<to; ++ccI, ++actcontigpos){
    gettimeofday(&us_loop,nullptr);
    // updating the pcrIs of the reads at that position
    CEBUG("cc acp: " << actcontigpos << endl);
    // first delete those who fall out at this new position
    for(uint32 actseqtype=0; actseqtype<cached_contigseqtypes.size(); ++actseqtype){
      if(cached_contigseqtypes[actseqtype]){
	CEBUG("cc art " << actseqtype << " rpic " << read_pcrIs_in_col[actseqtype].size() << endl);
	auto Ifrom=read_pcrIs_in_col[actseqtype].begin();
	auto Ito=Ifrom;
	// pcrIs are expensive on copy: by not blindly copying every iterator but only when needed,
	//  we save between 30% and 40% time (maybe 5% total in a total cons. calc, but still)
	bool needcopy=false;
	for(;Ifrom != read_pcrIs_in_col[actseqtype].end(); ++Ifrom){
	  if(needcopy) *Ito=*Ifrom;
	  if(Ifrom->getReadStartOffset()+(*Ifrom)->getLenClippedSeq() > actcontigpos) {
	    ++Ito;
	  }else{
	    needcopy=true;
	    CEBUG("cc thrown out " << *Ito << endl);
	  }
	}
	//if(Ito != Ifrom) {
	//  read_pcrIs_in_col[actseqtype].resize(Ito-read_pcrIs_in_col[actseqtype].begin());
	//}
	// PlacedContigReads::const_iterator has no default constructor,
	//  i.e., resize cannot be used because it wants a default constructor
	//  though we always reduce the vector here
	// but we can use erase()
	CEBUG("loopend" << endl);
	if(needcopy){
	  CEBUG("cc popping " << Ito - Ifrom << " elements\n");
	  read_pcrIs_in_col[actseqtype].erase(Ito,Ifrom);
	}
      }
    }
    mict_fallout+=diffsuseconds(us_loop);
    gettimeofday(&us_loop,nullptr);

    // now insert ids of reads that have newly started at this position
    // Don't take railreads, backbone is there.
    for(;pcrI != CON_reads.end() && pcrI.getReadStartOffset() == actcontigpos; ++pcrI){
      if( ! pcrI->isRail()) {

	// new! one can also make a choice which strains to take into account
	// either all (strainidtotake<0 or exactly the one chosen)
	if(strainidtotake < 0
	   || pcrI->getStrainID() == strainidtotake) {
	  read_pcrIs_in_col[pcrI->getSequencingType()].push_back(pcrI);
	  CEBUG("cc taken " << pcrI << endl);
	}
      }
    }
    mict_newin+=diffsuseconds(us_loop);
    gettimeofday(&us_loop,nullptr);

    // for each read type, predict a base and a quality,
    //  but also get back the potential solutions and number thereof

    //vector<uint32> numsolutionsvec(ReadGroupLib::getNumSequencingTypes(),0);

    uint32 numreadtypeswithsolution=0;
    for(uint32 actseqtype=0; actseqtype<cached_contigseqtypes.size(); ++actseqtype){
      predictedbases[actseqtype]='!';
      predictedquality[actseqtype]=0;

      if(cached_contigseqtypes[actseqtype]){
	makeIntelligentConsensus_helper1(predictedbases[actseqtype],
					 predictedquality[actseqtype],
					 actcontigpos,
					 ccI,
					 mincoverage,
					 groupsvec[actseqtype],
					 maskedshadowgroups,
					 IUPACbasegroups,
					 columnbases,
					 read_pcrIs_in_col[actseqtype],
					 maskshadow,
					 actseqtype,
					 strainidtotake,
					 missingcoveragechar
	  );
	if(!IUPACbasegroups.empty()) possiblebases[actseqtype]=IUPACbasegroups;
      }

      if(predictedbases[actseqtype] != '!') numreadtypeswithsolution++;
    }

#ifndef CLOCK_STEPS_CONS
    mict_helper1+=diffsuseconds(us_loop);
#endif
    gettimeofday(&us_loop,nullptr);


    /* after the previous loop, we have the following important things
       per read type (other things too, less important at the moment):
       - the picked base and quality in predictedbases[] and
         predictedquality[]
       - possible alternatives (including picked one) in possiblebases[]
       - all group calculations for all read types in groupsvec[]
     */

    CEBUG("preda: ");
    for(uint32 actreadtype=0; actreadtype<ReadGroupLib::getNumSequencingTypes(); actreadtype++){
      CEBUG(predictedbases[actreadtype]);
    }
    CEBUG('\n');

    bool mismatchesintyperesults=false;
    if(numreadtypeswithsolution>1) {
      CEBUG("numreadtypes: " << numreadtypeswithsolution << endl);
      char basepick='!';
      for(uint32 actreadtype=0; actreadtype<ReadGroupLib::getNumSequencingTypes(); actreadtype++){
	CEBUG("pred[" << actreadtype << "]: " << predictedbases[actreadtype] << endl);
	if(predictedbases[actreadtype] != '!'){
	  if(basepick == '!'){
	    basepick=toupper(predictedbases[actreadtype]);
	  }else if(toupper(predictedbases[actreadtype]) != basepick){
	    if(basepick=='X'){
	      if(dptools::isValidIUPACStarBase(predictedbases[actreadtype])){
		basepick=toupper(predictedbases[actreadtype]);
	      }
	    }else if(basepick=='N'){
	      if(toupper(predictedbases[actreadtype]) != 'X'){
		basepick=toupper(predictedbases[actreadtype]);
	      }
	    }else{
	      if(toupper(predictedbases[actreadtype]) != 'X'
		 && toupper(predictedbases[actreadtype]) != 'N'){
		mismatchesintyperesults=true;
	      }
	    }
	  }
	}
      }
      CEBUG("Mismatch1: " << mismatchesintyperesults << endl);
    }


    // Party time, let's find out which base we take
    char thisbase='@';
    base_quality_t thisqual=0;

    bool hasSOLiD=false;
    if(hasSOLiD){
      MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 2.");
      break;
    }

    // if we have a mismatch: shootout
    if(mismatchesintyperesults){

      CEBUG("Mismatch. ");

      // look first at at uniqueness and quality

      // set default level 100 (not used) for all seqtypes
      predictlevel.clear();
      predictlevel.resize(ReadGroupLib::getNumSequencingTypes(),100);

      // count how many predictions of each level are present
      predictlevelused.clear();
      predictlevelused.resize(10,0);

      int8 bestpred=100;
      for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	// just do this for sanger, 454, PacBio, solexa, solid equally (might change later)

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

	switch(actseqtype){
	case ReadGroupLib::SEQTYPE_SANGER :
	case ReadGroupLib::SEQTYPE_454GS20 :
	case ReadGroupLib::SEQTYPE_IONTORRENT :
	case ReadGroupLib::SEQTYPE_PACBIOHQ :
	case ReadGroupLib::SEQTYPE_PACBIOLQ :
	case ReadGroupLib::SEQTYPE_TEXT :
	case ReadGroupLib::SEQTYPE_SOLEXA :
	case ReadGroupLib::SEQTYPE_ABISOLID : {
	  if(predictedbases[actseqtype] != '!') {
	    // TODO: still not optimal, should perhaps test all
	    //  goups of IUPAC *sigh*
	    if(dptools::isValidACGTStarBase(predictedbases[actseqtype])){
	      predictlevel[actseqtype]=rateGoodnessLevelOfConsensus(
		ccI,
		groupsvec[actseqtype][dptools::getIndexOfBase(predictedbases[actseqtype])],
		possiblebases[actseqtype].size(),
		actseqtype
		);
	      BUGIFTHROW(predictlevel[actseqtype]<0, "problem in hybrid cons calc: level < 0");
	      BUGIFTHROW(predictlevel[actseqtype]>=predictlevelused.size(), "level > expected maxiumum of 9");
	    }else{
	      predictlevel[actseqtype]=5;
	    }
	    predictlevelused[predictlevel[actseqtype]]++;
	    if(predictlevel[actseqtype] < bestpred) bestpred=predictlevel[actseqtype];
	    CEBUG("plevel " << predictedbases[actseqtype] << " [" << actseqtype << "]: " << static_cast<int16>(predictlevel[actseqtype]) << '\n');
	  }
	  break;
	}
	default: {
	  //throw Notify(Notify::INTERNAL, THISFUNC, "unexpected sequencing type 1");
	}
	}
      }

      CEBUG("bestpred: " << static_cast<int16>(bestpred) << '\n');

      ostringstream tagstr;
      ostringstream tmdstr; // tag method details

      // if there's more than one prediction with quality level of best
      //  prediction or one after, it's unresolved
      // else we'll take the best prediction
      if(predictlevelused[bestpred]+predictlevelused[bestpred+1] >= 2){
	/* if this happens, we have different base predictions for each
	   sequencing type. That's bad. Could be due to sequencing or
	   alignment or assembly error. Or when the strain sequenced
	   with one sequencing type was not 100% identical to the sequenced
	   strain with the other sequencing type(s).

	   But there's one saving attempt: if one is a normal base and the other(s)
	   an IUPAC and the normal base is contained in the IUPAC, then we'll
	   take this base as result and not the IUPAC. (gaps kill the saving
	   process, sorry)
	*/


	// saving attempt, pre-check
	thisbase = '!';
	bool canbesaved=true;
	for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	  if(predictlevel[actseqtype]==bestpred
	     || predictlevel[actseqtype]==bestpred-1) {
	    // make sure that gaps kill the whole saving process
	    if(predictedbases[actseqtype]=='*') canbesaved=false;

	    if(dptools::isValidBase(predictedbases[actseqtype])){
	      tmdstr << " / " << ReadGroupLib::getShortNameOfSequencingType(actseqtype)
		     << " (" << static_cast<uint16>(predictlevel[actseqtype])
		     << ") " << predictedbases[actseqtype]
		     << " " << static_cast<uint16>(predictedquality[actseqtype]);

	      if(thisbase=='!'){
		thisbase=predictedbases[actseqtype];
	      }else{
		thisbase=dptools::calcIUPACConsensus(thisbase,
						     predictedbases[actseqtype]);
	      }
	    }
	  }
	}

	canbesaved&=dptools::isValidBase(thisbase);

	// saving attempt, post-check (does the base fit in all IUPACs)
	if(canbesaved) {
	  tmdstr << "  is in  ";
	  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	    if(predictlevel[actseqtype]==bestpred
	       || predictlevel[actseqtype]==bestpred-1) {
	      if(!dptools::isValidBase(predictedbases[actseqtype])){
		if(!dptools::hasNucleicAcidInIUPAC(thisbase,predictedbases[actseqtype])){
		  canbesaved=false;
		} else {
		  tmdstr << " / " << ReadGroupLib::getShortNameOfSequencingType(actseqtype)
			 << " (" << static_cast<uint16>(predictlevel[actseqtype])
			 << ") " << predictedbases[actseqtype]
			 << " " << static_cast<uint16>(predictedquality[actseqtype]);
		}
	      }
	    }
	  }
	}

	CEBUG("TMDstr: " << tmdstr.str() << '\n');

	multitag_t::mte_id_t tagtype=CON_tagentry_idSTMU;
	if(canbesaved){
	  tagtype=CON_tagentry_idSTMS;
	  tagstr << "Solved vote: " << tmdstr.str();
	  // TODO: refine this: other quality and tell which sequencing types with which
	  //  results were taken into account
	  thisqual=10;
	}else{
	  tagstr << "Unresolved mismatch in seq. type consensus";

	  thisbase = '!';
	  thisqual=0;
	  for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	    if(predictlevel[actseqtype]==bestpred
	       || predictlevel[actseqtype]==bestpred-1) {
	      tagstr << " / " << ReadGroupLib::getShortNameOfSequencingType(actseqtype)
		     << ": (";
	      if(thisbase=='!'){
		thisbase=predictedbases[actseqtype];
	      }else{
		thisbase=dptools::calcIUPACConsensus(thisbase,
						     predictedbases[actseqtype]);
	      }
	      tagstr << static_cast<uint16>(predictlevel[actseqtype])
		     << ") " << predictedbases[actseqtype]
		     << " " << static_cast<uint16>(predictedquality[actseqtype]);
	    }
	  }
	  tagstr << " || " << thisbase << " " << static_cast<uint16>(thisqual);
	}
	CEBUG("TAGstr1: " << tagstr.str() << '\n');
	addTagToConsensus(actcontigpos,
			  actcontigpos,
			  '=',
			  multitag_t::getIdentifierStr(tagtype).c_str(),
			  tagstr.str().c_str(),
			  true);
      }else{
	// well, shootout: take the best one
	tagstr << "Solved mismatch: ";

	for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	  if(predictlevel[actseqtype]==bestpred) {
	    tagstr << ' ' << ReadGroupLib::getShortNameOfSequencingType(actseqtype)
		   << ": (";
	    thisbase=predictedbases[actseqtype];
	    thisqual=predictedquality[actseqtype];

	    tagstr << static_cast<uint16>(predictlevel[actseqtype])
		   << ") " << thisbase
		   << ' ' << static_cast<uint16>(thisqual);

	    break;
	  }
	}
	CEBUG("TAGstr2: " << tagstr.str() << '\n');
	addTagToConsensus(actcontigpos,
			  actcontigpos,
			  '=',
			  multitag_t::getIdentifierStr(CON_tagentry_idSTMS).c_str(),
			  tagstr.str().c_str(),
			  true);
      }
    }else{
      // no mismatch in type results, good
      // quality = sum(.75*quality of each seqtype) or best qual
      CEBUG("No mismatch. ");

      uint32 compqual=0;
      uint32 bestqual=0;
      bool hassomething=false;
      bool has_n=false;
      bool has_x=false;
      for(uint32 actseqtype=0; actseqtype<ReadGroupLib::getNumSequencingTypes(); actseqtype++){
	if(predictedbases[actseqtype] == 'N') {
	  has_n=true;
	  continue;
	}
	if(predictedbases[actseqtype] == 'X') {
	  has_x=true;
	  continue;
	}
	if(predictedbases[actseqtype] != '!') {

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
	  switch(actseqtype){
	  case ReadGroupLib::SEQTYPE_SANGER :
	  case ReadGroupLib::SEQTYPE_454GS20 :
	  case ReadGroupLib::SEQTYPE_IONTORRENT :
	  case ReadGroupLib::SEQTYPE_PACBIOHQ :
	  case ReadGroupLib::SEQTYPE_PACBIOLQ :
	  case ReadGroupLib::SEQTYPE_TEXT :
	  case ReadGroupLib::SEQTYPE_SOLEXA :
	  case ReadGroupLib::SEQTYPE_ABISOLID: {
	    hassomething=true;
	    thisbase=predictedbases[actseqtype];
	    compqual+=(static_cast<uint32>(predictedquality[actseqtype])*75)/100;
	    if(predictedquality[actseqtype]>bestqual) bestqual=predictedquality[actseqtype];
	    break;
	  }
	  default: {
	    // do nothing???
	    //thisbase='x';
	  }
	  }
	}
      }

      if(!hassomething){
	if(has_x) {
	  thisbase='x';
	  hassomething=true;
	}
	if(has_n) {
	  thisbase='N';
	  hassomething=true;
	}
      }

      if(hassomething){
	BUGIFTHROW(thisbase=='@',"Coverage present but base is '@' ... this is not healthy.\n");
      }else{
	BUGIFTHROW(thisbase!='@',"No coverage present but base is not '@' ... this is not healthy.\n");
	thisbase=missingcoveragechar;
	thisqual=0;
      }

      if(bestqual>compqual) compqual=bestqual;
      if(compqual>90) compqual=90;
      thisqual=static_cast<base_quality_t>(compqual);
    }

    CEBUG("thisbase " << thisbase << "\tthisqual " << (uint16) thisqual << endl);

    thisbase=tolower(thisbase);

    if(thisqual<minqual && thisbase!=missingcoveragechar){
      thisbase='N';
      thisqual=0;
      CEBUG("Minqual not reached, changed to: thisbase " << thisbase << "\tthisqual " << (uint16) thisqual << endl);
    }

    if(thisqual>90) thisqual=90;

    target+=thisbase;
    qual.push_back(thisqual);

    // calc the adjustments
    if(targetadjustments!= nullptr && from==0 && to==CON_counts.size()) {
      if(thisbase=='*') {
	targetadjustments->push_back(-1);
      }else{
	targetadjustments->push_back(unpaddedposcounter);
	unpaddedposcounter++;
      }

      BUGIFTHROW(targetadjustments->size()!=target.size(), "gna1");
      BUGIFTHROW(targetadjustments->size()!=qual.size(), "gna2");

    }

/*
    // dump out .tcs file to ostr if given
    if(ostr!=nullptr) {
      //throw Notify(Notify::INTERNAL, THISFUNC, "must be adapted to multiple read types");

      // name and position
      // don't call paddedPos2UnpaddedPos as this would lead to recursion!
      *ostr << setw(20) << left << getContigName()
	    << setw(9) << right << actcontigpos
	    << setw(9) << CON_adjustments.back()
	    << " | " << thisbase
	    << setw(3) << static_cast<uint16>(thisqual)
	    << " |" << setw(5) << read_pcrIs_in_col.size();

      tcs_totalgroupcount.clear();
      tcs_totalgroupcount.resize(emptygroups.size(),0);
      tcs_totalgroupqual.clear();
      tcs_totalgroupqual.resize(emptygroups.size(),0);
      for(uint32 actseqtype=0; actseqtype<groupsvec.size(); actseqtype++){
	for(uint32 grpi=0; grpi<groupsvec[actseqtype].size(); grpi++){
	  tcs_totalgroupcount[grpi]+=groupsvec[actseqtype][grpi].urdids.size();
	  tcs_totalgroupqual[grpi]+=groupsvec[actseqtype][grpi].groupquality;
	  if(tcs_totalgroupqual[grpi]>90) tcs_totalgroupqual[grpi]=90;
	}
      }
      for(uint32 grpi=0; grpi<tcs_totalgroupcount.size(); grpi++){
	*ostr << setw(5) << tcs_totalgroupcount[grpi];
      }

      *ostr << " |";
      uint32 numgroupswithqual=0;
      for(uint32 grpi=0; grpi<tcs_totalgroupcount.size(); grpi++){
	if(tcs_totalgroupcount[grpi]) {
	  *ostr << setw(3) << tcs_totalgroupqual[grpi];
	  numgroupswithqual++;
	}else{
	  *ostr << " --";
	}
      }

      // TODO: different characters for different cases?
      char bstatus='?';
      if(!dptools::isValidBase(thisbase)
	 && dptools::isValidIUPACBase(thisbase)){
	bstatus='M';
      }else if(thisbase=='*'){
	if(thisqual>=40){
	  bstatus=':';
	}else if(thisqual<30){
	  bstatus='m';
	}else if(numgroupswithqual>2){
	  bstatus='m';
	}else{
	  bstatus=':';
	}
      }else if(thisqual<30 && numgroupswithqual>1){
	bstatus='m';
      }else{
	bstatus=':';
      }
      if(!tcs_aconstagpositions.empty() && tcs_dconstagpositions[actcontigpos]){
	bstatus='$';
      }
      if(bstatus!=':'){
	*ostr << " | !" << bstatus << " |";
      }else{
	*ostr << " |  : |";
      }

      if(!tcs_aconstagpositions.empty() && tcs_aconstagpositions[actcontigpos]){
	vector<consensustag_t>::const_iterator I=CON_consensus_tags.begin();
	bool doneoutput=false;
	for(; I!=CON_consensus_tags.end(); I++) {
	  if((I->from <= actcontigpos) && ((I->to) >= actcontigpos)) {
	    if(doneoutput){
	      *ostr << ' ' << I->getIdentifierStr();
	    }else{
	      *ostr << " \"" << I->getIdentifierStr();
	      doneoutput=true;
	    }
	  }
	}
	if(doneoutput) *ostr << "\"";

      }
      *ostr << "\n";
    }
*/

    mict_restofloop+=diffsuseconds(us_loop);
  }
  mict_totalloop=diffsuseconds(us_start);

  cout << "mict_fin        " << mict_fin << endl;
  cout << "mict_pre        " << mict_pre << endl;
  cout << "mict_shadow     " << mict_shadow << endl;
  cout << "mict_fallout    " << mict_fallout << endl;
  cout << "mict_newin      " << mict_newin << endl;
#ifndef CLOCK_STEPS_CONS
  cout << "mict_helper1    " << mict_helper1 << endl;
#else
  cout << "mict_helper1    " << CON_us_steps_cons[USCLOCONS_H1_TOTAL] << endl;
  cout << "  H1_MGROUPS    " << CON_us_steps_cons[USCLOCONS_H1_MGROUPS] << endl;
  cout << "  H1_PCRI       " << CON_us_steps_cons[USCLOCONS_H1_PCRI] << endl;
  cout << "    H1_CP2URP   " << CON_us_steps_cons[USCLOCONS_H1_CP2URP] << endl;
  cout << "  H1_EGROUP     " << CON_us_steps_cons[USCLOCONS_H1_EGROUP] << endl;
  cout << "  H1_GQUAL      " << CON_us_steps_cons[USCLOCONS_H1_GQUAL] << endl;
  cout << "  H1_CALLH2     " << CON_us_steps_cons[USCLOCONS_H1_CALLH2] << endl;
#endif
  cout << "mict_newin      " <<   mict_newin << endl;
  cout << "mict_restofloop " <<   mict_restofloop << endl;
  cout << "mict_totalloop  " <<   mict_totalloop << endl;

  // seconds on contig with ~4m reads
  //mict_pre                  55
  //mict_shadow           536072
  //mict_fallout        5 855515
  //mict_newin            330100
  //mict_helper1       28 934855
  //mict_restofloop       288773
  //mict_totalloop     35 875607


  CEBUG("target (first 100): " << target.substr(0,100) << endl);

  cout << "done." << endl;

  FUNCEND();

  return;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *
 *************************************************************************/

int8 Contig::rateGoodnessLevelOfConsensus(cccontainer_t::const_iterator ccI, nngroups_t & group, uint32 numpossiblebases, uint8 seqtype)
{
  FUNCSTART("int8 Contig::rateGoodnessLevelOfConsensus(cccontainer_t::const_iterator ccI, nngroups_t & group, uint32 numpossiblebases, uint8 seqtype)");

  // basically, we'll prefer read types that have forward and
  //  complement direction present.

  CEBUG("predicting level for " << group.base << " seqtype " << static_cast<uint16>(seqtype) << '\n' << group << '\n');
  int8 level=-1;

  if(seqtype==ReadGroupLib::SEQTYPE_ABISOLID){
    MIRANOTIFY(Notify::INTERNAL, "Type ABI SOLiD needs more support 3.");
  }


  bool hasfr=(group.forwarddircounter>0) && (group.complementdircounter>0);

#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif

  switch(seqtype) {
  case ReadGroupLib::SEQTYPE_SANGER :{
    if(hasfr &&
       group.groupquality >= 35){
      level=0;
    }else if(hasfr
	     && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=3
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_454GS20 :{
    if(hasfr
       && group.urdids.size() >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 8
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=6
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_IONTORRENT :{
    if(hasfr
       && group.urdids.size() >= 8
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 6
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=4
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_PACBIOHQ :
  case ReadGroupLib::SEQTYPE_PACBIOLQ :{
    // No info. atm use same as 454
    // TODO: change that for pacbio hq and lq
    if(hasfr
       && group.urdids.size() >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 8
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=6
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_TEXT :{
    // TODO: that's just a first shot at atm, would need to refine
    if(hasfr
       && group.urdids.size() >= 8
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && group.urdids.size() >= 6
	      && group.groupquality >= 25){
      level=1;
    } else if(hasfr){
      level=2;
    } else if(group.urdids.size()>=4
	      && numpossiblebases == 1){
      level=3;
    } else {
      level=4;
    }
    break;
  }
  case ReadGroupLib::SEQTYPE_SOLEXA :{
    // Solexa special case: we need to take the merged reads
    //  into account
    uint32 realsize=group.urdids.size();
    if(ccI->getOriginalBBChar()==group.base){
      realsize+=ccI->bbcounts[0];
      if(ccI->bbcounts[0]>=6){
	hasfr=true;
      }
    }
    if(hasfr
       && realsize >= 12
       && group.groupquality >= 35){
      level=0;
    } else if(hasfr
	      && realsize >= 8
	      && group.groupquality >= 30){
      level=1;
    } else if(hasfr
	      && group.groupquality >= 25){
      level=2;
    } else if(group.groupquality >= 20
	      && realsize>=4
	      && numpossiblebases == 1){
      level=3;
    } else if(group.groupquality >= 10){
      level=4;
    } else {
      level=5;
    }
    break;
  }
  default : {
    cout << "seqtype: " << seqtype << endl;
    MIRANOTIFY(Notify::INTERNAL, "Uknown seqtype to rate?");
  }
  }

  FUNCEND();
  return level;
}

//#define CEBUG(bla)




/*************************************************************************
 *
 * makes a temporary char * of a portion of the consensus
 * from and to positions in the consensus
 * from inclusive
 * to exclusive
 * appending a 0 char
 *
 * returns whether a base not defined by *original* backbone or, if backbone,
 *  non-ACGT* was put into consensus
 *  Used to define whether read can be mapped or not.
 *
 *************************************************************************/

bool Contig::makeTmpConsensus(int32 from, int32 to, bool tmpconsfrombackbone)
{
//#define CEBUG(bla)   {cout << bla; cout.flush();}

  FUNCSTART("void Contig::makeTmpConsensus(int32 from, int32 to)");

  definalise();
  CEBUG("\nFrom: " << from << "\tTo: " << to<<"\tCON_counts.size(): "<<CON_counts.size());

  BUGIFTHROW(from>to,"from>to?");
  BUGIFTHROW(from<0, "from < 0 ?");

  if(to>CON_counts.size()) to=CON_counts.size();
  uint32 len_tmpcons=to-from;

  if(CON_2tmpcons.capacity()<2000 || CON_2tmpcons.capacity()<len_tmpcons){
    CON_2tmpcons.reserve(max(len_tmpcons,static_cast<uint32>(2000)));
  }
  CON_2tmpcons.resize(len_tmpcons);

  bool hasNonBBMappable=!tmpconsfrombackbone;
  {
    auto toptr=CON_2tmpcons.begin();
    auto ccI=CON_counts.cbegin();
    BOUNDCHECK(from, 0, CON_counts.size()+1);
    advance(ccI, from);

    for(uint32 i=from; i<to ;++i, ++toptr, ++ccI){
      *toptr=0;
      if(tmpconsfrombackbone && ccI->getOriginalBBChar()!='@'){
	*toptr=ccI->getBBChar();
	hasNonBBMappable |= !dptools::isValidACGTStarBase(ccI->getOriginalBBChar());
	if(ccI->i_backbonecharupdated!='@'
	   && ccI->i_backbonecharorig != ccI->i_backbonecharupdated) hasNonBBMappable=true;
      }else{
	hasNonBBMappable=true;

	ccctype_t maximum= max(ccI->A, max(ccI->C, max(ccI->G, ccI->T)));
	uint8 counts=0;
	//CEBUGF(ccI->A << "\t" << ccI->C << "\t" << ccI->G << "\t" << ccI->T << "\t" << ccI->N << "\t" << ccI->star << "\n");

	// is any ACGT set?
	if(maximum >0 && maximum > ccI->star) {
	  if(ccI->A==maximum){
	    counts++;
	    *toptr='A';
	  }
	  if(ccI->C==maximum){
	    counts++;
	    *toptr='C';
	  }
	  if(ccI->G==maximum){
	    counts++;
	    *toptr='G';
	  }
	  if(ccI->T==maximum){
	    counts++;
	    *toptr='T';
	  }
	  if(counts>1) {
	    *toptr='N';
	  }

	  //// can be somewhat problematic with 454 data
	  //// calls the base until 50/50, then the gap
	  //if(maximum/4 < ccI->star) *toptr='*';

	  // this prefers to call gaps
	  // calls the base until 1/3 base, 2/3 gap, then the gap
	  // this should help, together with the "expected gap" # in
	  //  alignments, to further reduce to a maximum this kind of
	  //  base jiggling in homopolymers
	  //
	  //     ...*AAAAAAAAA...
	  //     ...*AAAAAAAAA...
	  //     ...AAAAAAAAA*...
	  //     ...AAAAAAAAA*...
	  //     ...AAAAAAAAA*...
	  //     ...*AAAAAAAAA...
	  //     ...*AAAAAAAAA...
	  //     ...AAAAAAAAA*...
	  //     ...*AAAAAAAAA...

	  if(maximum/4 < (ccI->star)*2) {
	    switch(*toptr){
	    case 'A': {
	      *toptr='1';
	      break;
	    }
	    case 'C': {
	      *toptr='2';
	      break;
	    }
	    case 'G': {
	      *toptr='3';
	      break;
	    }
	    case 'T': {
	      *toptr='4';
	      break;
	    }
	    default: {
	      *toptr='*';
	    }
	    }
	  }
	} else {
	  if(unlikely(ccI->total_cov==0)){
	    // BaCh 30.11.2012
	    // should normally never happen, certainly not in de-novo
	    // but the two-pass mapping may have this at the end of the contigs after first pass
	    //  (should I decide not to go the chompFront() / chompBack() after 1st pass)
	    //
	    // treat it like a base (well, will be X)
	    *toptr='N';
	  }else if((ccI->star >= ccI->X)
	     && (ccI->star >= ccI->N)){
	    *toptr='*';
	  } else if(ccI->N){
	    *toptr='N';
	  }else{
	    *toptr='X';
	  }
	}
      }
      if(*toptr==0){
	MIRANOTIFY(Notify::INTERNAL,"Ooooops? makeTmpConsensus encountered the unexpected situation of an uncalled base? Please contact the author immediately.");
      }
    }
  }
  CEBUG("Tmp_cons: >>>" <<CON_2tmpcons << "<<<" << endl);

  FUNCEND();

  return hasNonBBMappable;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Contig::updateBackboneConsensus()
{
  makeTmpConsensus(0,CON_counts.size(),false);
  auto ccI=CON_counts.begin();
  auto c2tc=CON_2tmpcons.c_str();

  // TODO: check for coverage, check for fwd,rev
  rcci_t rcci(this);
  {
    vector<int32> allowedstrainids; // empty would be all ids
    vector<uint8> allowedreadtypes;
    rcci.init(allowedstrainids,
	      allowedreadtypes,
	      false,            // don't take rails
	      false,           // nor backbones
	      false);   // nor reads without readpool-reads
  }

  for(; ccI != CON_counts.end(); ++ccI, ++c2tc, rcci.advance()){
    if(rcci.getPCRIsInCol().size()>=6){
      uint32 fwd=0;
      uint32 rev=0;
      for(auto & pcrI : rcci.getPCRIsInCol()){
	if(pcrI.getReadDirection()>0){
	  ++fwd;
	}else{
	  ++rev;
	}
      }
      if(fwd>=2 && rev>=2){
	if(ccI->i_backbonecharorig!='@'){
	  ccI->i_backbonecharupdated=*c2tc;
	}
      }
    }
  }
}
