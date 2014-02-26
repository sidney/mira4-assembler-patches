/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2011 and later by Bastien Chevreux
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 *
 */

#include "skim.H"

#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include "errorhandling/errorhandling.H"


using namespace std;



//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Skim::skimStreamPrepare(ReadPool & rp,
			     uint8  bph,
			     uint8  hss,
			     const char * additionalregexp)
{
  FUNCSTART("uint32 Skim::skimGo( ... )");

  init();

  SKIM3_readpool=&rp;

  if(sizeof(vhash_t)==4){
    if(bph>14) bph=14;
  }
  if(sizeof(vhash_t)==8){
    if(bph>30) bph=30;
  }
  SKIM3_basesperhash=bph;
  SKIM3_hashsavestepping=hss;

  fillTagStatusInfoOfReads();

  prepareSkim(0, rp.size(), SKIM3_vhraparray,false);

  FUNCEND();
  return;
}
//#define CEBUG(bla)






/*************************************************************************
 *
 * seqtype -1 == all sequencing types, else only the given
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void Skim::findAdaptorRightClip(ReadPool & searchpool, vector<int32> & results, int8 seqtype, uint32 minhashes, uint32 numthreads)
{
  FUNCSTART("void Skim::checkForAdaptor(Read & actread, cfh_threaddata_t & cfhd)");

  results.clear();
  results.resize(searchpool.size(),-1);

  SKIM3_farc_searchpool=&searchpool;
  SKIM3_farc_results=&results;
  SKIM3_farc_minhashes=minhashes;
  SKIM3_farc_seqtype=seqtype;

  startMultiThreading(1,numthreads,10000,0,searchpool.size(),
		      boost::bind( &Skim::farcThreadsDataInit, this, _1 ),
		      boost::bind( &Skim::farcThreadLoop, this, _1 ));

  FUNCEND();
}

int32 Skim::findAdaptorRightClip(Read & actread, uint32 minhashes, readid_t & ridofadapfound, int32 threadid)
{
  FUNCSTART("int32 Skim::findAdaptorRightClip(Read & actread, uint32 minhashes. int32 threadid)");

  int32 ret=0;
  if(threadid<0){
    ret=findAdaptorRightClip_internal(actread,minhashes,ridofadapfound,SKIM3_farcdata_fornonmultithread);
  }else{
    BUGIFTHROW(threadid>=SKIM3_farcd_vector.size(),"Oooooops, trying to use thread id " << threadid << ", but prepared only for " << SKIM3_farcd_vector.size() << " threads???");
    ret=findAdaptorRightClip_internal(actread,minhashes, ridofadapfound, SKIM3_farcd_vector[threadid]);
  }
  return ret;
}

#define CEBUG2(bla)
//#define CEBUG2(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush();}

int32 Skim::findAdaptorRightClip_internal(Read & actread, uint32 minhashes, readid_t & ridofadapfound, farc_threaddata_t & farcd)
{
  FUNCSTART("void Skim::checkForAdaptor(Read & actread, cfh_threaddata_t & cfhd)");

  CEBUG("farc_i: " << actread.getName() << endl);

  if(SKIM3_vashortcuts_begin.empty() || SKIM3_vashortcuts_end.empty()) return -1;
  if(!actread.hasValidData()) return -1;
  uint32 slen=actread.getLenClippedSeq();
  if(slen<SKIM3_basesperhash) return -1;

  // don't really need to clear out this re-used vector
  //   - singlereadvhraparray needs to be big enough to be
  //     written into by transformSeqToVariableHash()
  if(farcd.singlereadvhraparray.size() < slen){
    farcd.singlereadvhraparray.resize(slen);
  }

  farcd.readhashmatches.clear();

  vector<vhrap_t>::iterator srvaI=farcd.singlereadvhraparray.begin();

  vector<Read::bposhashstat_t> & bposhashstats=const_cast<vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
  uint32 hashesmade;

  uint32 actreadid=0xffffffff;

  {
    int32 bfpos=0;
    int32 bfposinc=1;

    hashesmade=transformSeqToVariableHash(
      actreadid,
      actread,
      actread.getClippedSeqAsChar(),
      slen,
      SKIM3_basesperhash,
      srvaI,
      false,
      1,
      farcd.unused_tagmaskvector,
      bposhashstats,
      bfpos,
      bfposinc
      );
  }
  farcd.singlereadvhraparray.resize(hashesmade);

  CEBUG("hashesmade: " << hashesmade << endl);


  srvaI=farcd.singlereadvhraparray.begin();

  vector<vhrap_t>::const_iterator lowerbound;
  vector<vhrap_t>::const_iterator upperbound;
  for(; srvaI != farcd.singlereadvhraparray.end(); srvaI++){
    lowerbound=SKIM3_vashortcuts_begin[srvaI->vhash & SKIM3_MAXVHASHMASK];
    upperbound=SKIM3_vashortcuts_end[srvaI->vhash & SKIM3_MAXVHASHMASK];

    // "SKIM3_empty_vector_vhrap_t.end()" is the "nullptr" replacement
    if(SKIM3_completevhraparray_end != lowerbound){
      if(SKIM3_basesperhash>12){
	// with more than 12 bases in a hash, the vhrap array is
	//  subdivided
	pair<vector<vhrap_t>::const_iterator, vector<vhrap_t>::const_iterator>
	  p=equal_range(lowerbound,
			upperbound,
			*srvaI,
			compareVHRAPArrayElem_);
	lowerbound=p.first;
	upperbound=p.second;
      }

      for(;lowerbound!=upperbound; lowerbound++){
	//CEBUG("/// " << actreadid << '\t' << lowerbound->readid << '\n');
	//CEBUG("/// take!\n");
	farcd.readhashmatches.resize(farcd.readhashmatches.size()+1);
	farcd.readhashmatches.back().rid2=lowerbound->readid;
	farcd.readhashmatches.back().hashpos1=srvaI->hashpos;
	farcd.readhashmatches.back().hashpos2=lowerbound->hashpos;
	farcd.readhashmatches.back().eoffset=srvaI->hashpos - lowerbound->hashpos;
	farcd.readhashmatches.back().bhashstats=srvaI->bhashstats;

	CEBUG2("added: " << farcd.readhashmatches.back());
      }
    }
  }

  int32 retvalue=-1;

  if(!farcd.readhashmatches.empty()){
    checkForPotentialAdaptorHits(1, actreadid, actread, farcd.tmpmatchwith, farcd.readhashmatches);

    if(!farcd.tmpmatchwith.empty()){
      CEBUG2("Hits of: " << actread.getName() << endl);
      vector<matchwithsorter_t>::const_iterator ssmwsI=farcd.tmpmatchwith.begin();
      vector<matchwithsorter_t>::const_iterator leftmostokI=farcd.tmpmatchwith.end();
      vector<matchwithsorter_t>::const_iterator largestokI=farcd.tmpmatchwith.end();
      int32 leftmostpos=100000000;
      uint32 largesthash=0;
      for(; ssmwsI != farcd.tmpmatchwith.end(); ++ssmwsI){
	CEBUG2(actread.getName() << " to " << SKIM3_readpool->getRead(ssmwsI->otherid).getName() << " : " << *ssmwsI);
	// have at least an estimated 50% identity, false positive adaptor identifications normally have well below 50%
	if(ssmwsI->numhashes>=minhashes && ssmwsI->percent_in_overlap>=50){
	  int32 numhashes=ssmwsI->numhashes;
	  if(ssmwsI->eoffset < 0) numhashes+=ssmwsI->eoffset;
	  if(numhashes>=minhashes){
	    if(ssmwsI->eoffset < leftmostpos){
	      CEBUG2("newleft!!!!!!\n");
	      leftmostpos=ssmwsI->eoffset;
	      leftmostokI=ssmwsI;
	    }
	    if(numhashes>largesthash){
	      CEBUG2("newlarge!!!!!!\n");
	      largesthash=numhashes;
	      largestokI=ssmwsI;
	    }
	  }
	}
      }
//      if(leftmostokI!=farcd.tmpmatchwith.end()){
//	//cout << "Chosen " << SKIM3_readpool->getRead(leftmostokI->otherid).getName() << ": " << *leftmostokI;
//	ridofadapfound=leftmostokI->otherid;
//	retvalue=leftmostokI->eoffset;
      if(largestokI!=farcd.tmpmatchwith.end()){
	//cout << "Chosen " << SKIM3_readpool->getRead(largestokI->otherid).getName() << ": " << *largestokI;
	ridofadapfound=largestokI->otherid;
	retvalue=largestokI->eoffset;
	if(retvalue<0) retvalue=0;
	if(leftmostokI!=largestokI){
	  //cout << "Leftmost != Biggest " << SKIM3_readpool->getRead(largestokI->otherid).getName() << ": " << *largestokI;
	}
      }
    }
    //selectPotentialHitsForSave2(direction, actreadid,
    //				cfhd);

  }

  return retvalue;
}

//#define CEBUG(bla)





/*************************************************************************
 *
 * a simplified version of checkForPotentialHits()
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG_extra_cFPH

void Skim::checkForPotentialAdaptorHits(const int8 direction, const uint32 actreadid, Read & actread, vector<matchwithsorter_t> & tmpmatchwith, vector<readhashmatch_t> & readhashmatches)
{
  //bool dodebug=false;

  CEBUG("Potential hits of " << actread.getName() << " (" << actreadid << ") (" << static_cast<int16>(direction) << '/' << actread.getLenClippedSeq() << ")\n----------------\n");
  CEBUG(actread << endl);
  CEBUG("----------------\n");

  tmpmatchwith.clear();

  // readhashmatches should not be empty ... normally.
  // but new method to deal with megahubs reduces this vector, keeping only
  //  'approximately normal' frequencies. Which in turn means: some vectors
  //  might be completely emptied
  // so, if it is empty, return immediately
  if(readhashmatches.empty()) return;

  sort(readhashmatches.begin(), readhashmatches.end(), sortreadhashmatch_t_);

  vector<readhashmatch_t>::const_iterator sI=readhashmatches.begin();

  uint32 countid=sI->rid2;
  while(sI != readhashmatches.end()){

    uint32 rid2=sI->rid2;
    uint16 oldhashpos=sI->hashpos1;
    uint16 hp1min=0xffff;
    uint16 hp1max=0;

    uint16 hp2min=0xffff;
    uint16 hp2max=0;

    int32  eoffsetmin=0x7fffffff;
    int32  eoffsetmax=0x80000000;
    int32 oldeoffset=sI->eoffset;

    uint32 numhashes=0;

    vector<readhashmatch_t>::const_iterator sIS=sI;
    for(;sI != readhashmatches.end() && sI->rid2 == countid; sI++){
      CEBUG(*sI);

      // this ensures that the eoffset between two following
      //  entries may not differ by too much (2 bases here for adaptor search)
      // IF they do, then this is treated like a different hit
      //  by breaking the loop
      if(abs(sI->eoffset - oldeoffset) > 2){
	CEBUG("BREAKER!\n");
	break;
      }
      numhashes++;

      hp1min=min(hp1min,sI->hashpos1);
      hp1max=max(hp1max,sI->hashpos1);
      eoffsetmin=min(eoffsetmin,sI->eoffset);
      eoffsetmax=max(eoffsetmax,sI->eoffset);
      oldeoffset=sI->eoffset;

      hp2min=min(hp2min,sI->hashpos2);
      hp2max=max(hp2max,sI->hashpos2);

#ifdef CEBUG_extra_cFPH
      {
	boost::mutex::scoped_lock lock(SKIM3_coutmutex);
	CEBUG(sI->rid2
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getName()
	      << "\t" << SKIM3_readpool->getRead(sI->rid2).getLenClippedSeq()
	      << "\t" << sI->eoffset
	      << "\t" << sI->hashpos1
	      << "\t" << oldhashpos
	      << '\n');
      }
#endif

      oldhashpos=sI->hashpos1;
    }

    int32 maxoverlap;

    // adjust min positions for the hash length
    hp1min-=(SKIM3_basesperhash-1);
    hp2min-=(SKIM3_basesperhash-1);

    int32 eoffsetmean=eoffsetmin+(eoffsetmax-eoffsetmin)/2;

    // calc max overlap
    // currently only for one offset
    if(eoffsetmean<0){
      maxoverlap=min(SKIM3_readpool->getRead(rid2).getLenClippedSeq()+eoffsetmean,actread.getLenClippedSeq());
    }else{
      maxoverlap=min(actread.getLenClippedSeq()-eoffsetmean,SKIM3_readpool->getRead(rid2).getLenClippedSeq());
    }

    // correct the maxoverlap by the modulo of the hash steps as the
    //  border hashes will be found only in 1/(hash stepping) cases
    maxoverlap=maxoverlap-(maxoverlap%SKIM3_hashsavestepping);

    // hashe3soverlap is not the number of hashes in the overlap,
    // but the length of the overlap
    int32 hashesoverlap=hp1max-hp1min+1;

    int32 perc=100*hashesoverlap/maxoverlap;


    int32 minpercentrequired=0;

    // look a bit closer at potential perfect matches
    if(perc == 100){
      if(eoffsetmin != eoffsetmax){
	// this could not be: have at least two different expected offsets
	//  and a 100% coverage. Side effects from intra-read repeats
	// therefore, make sure this does not get through as a 100% match
	perc=99;
      }else if((numhashes-1)*SKIM3_hashsavestepping+SKIM3_basesperhash < maxoverlap){
	// maxoverlap covers the whole potential overlap, but
	//  there are not enough hashes supporting for 100% match
	//  (base mismatch somewhere)
	// reduce the percentage to show it's not a perfect match
	perc=99;
      }
    }else if(eoffsetmin == eoffsetmax){
      if(perc>100) {
	perc=100;
      }else{
	uint32 maxnumhashes=((maxoverlap-1-SKIM3_basesperhash)/SKIM3_hashsavestepping)+1;
	if(perc>=minpercentrequired && numhashes==maxnumhashes){
	  CEBUG("maxnumhashes 100% saver: "  << perc << '\n');
	  perc=100;
	}
      }
    }

    // 13 was a bit too small, still produces rare false positives even if percent_in_overlap >= 60
    // only workable countermeasure: 16 (and that is also more in-line with minimm overlaps one
    //  should use)
    int32 minoverlaprequired=16;

#ifdef CEBUG_extra_cFPH
    {
      boost::mutex::scoped_lock lock(SKIM3_coutmutex);
      CEBUG("eomin: " << eoffsetmin << "\teomax: " << eoffsetmax
	    << "\tmor: " << minoverlaprequired
	    << "\tho: " << hashesoverlap
	    << "\t%: " << perc
	    << "\t%<: " << minpercentrequired << endl);
    }
#endif

    // we take the hit if the overlap percentage is above threshold
    if(hashesoverlap >= minoverlaprequired
       && perc>=minpercentrequired){

//#define CEBUG(bla)   {if(actreadid==273252 && rid2==273250) cout << bla; cout.flush();}
//#define CEBUG(bla)   {cout << bla; cout.flush();}
      matchwithsorter_t tmp;
      tmp.otherid=rid2;
      tmp.eoffset=eoffsetmean;

      if(perc>100) perc=100;
      tmp.percent_in_overlap=perc;
      tmp.numhashes=numhashes;
      tmp.estimscore=0;
      tmp.taken=false;

      tmpmatchwith.push_back(tmp);

      CEBUG("Pushing possible hit with offset: " << tmp.eoffset << endl
	    << rid2
	    << "\t" << actreadid
	    << "\t" << SKIM3_readpool->getRead(rid2).getLenClippedSeq()
	    << "\t" << hp1min
	    << "\t" << hp1max
	    << "\t" << eoffsetmin
	    << "\t" << eoffsetmax
	    << "\t" << maxoverlap
	    << "\t" << hashesoverlap
	    << "\t" << numhashes
	    << "\t" << minoverlaprequired
	    << "\t" << perc << '%'
	    << '\n');
//#define CEBUG(bla)

    }
    if(sI!=readhashmatches.end()) countid=sI->rid2;
  }

}
//#define CEBUG(bla)



void Skim::prepareForMultithreadFarc(uint32 numthreads)
{
  FUNCSTART("void Skim::prepareForMultithreadFarc()");

  SKIM3_numthreads=numthreads;
  farcThreadsDataInit(SKIM3_numthreads);
}

void Skim::farcThreadsDataInit(const uint32 numthreads)
{
  FUNCSTART("void Skim::cfhThreadsDataInit(const uint32 numthreads)");

  SKIM3_farcd_vector.resize(numthreads);
  for(uint32 ti=0; ti<numthreads;++ti){

    SKIM3_farcd_vector[ti].readhashmatches.clear();
    SKIM3_farcd_vector[ti].readhashmatches.reserve(2000);
    SKIM3_farcd_vector[ti].singlereadvhraparray.clear();
    SKIM3_farcd_vector[ti].singlereadvhraparray.reserve(2000);
    SKIM3_farcd_vector[ti].unused_tagmaskvector.clear();
    //SKIM3_farcd_vector[ti].unused_tagmaskvector.reserve(2000);
    SKIM3_farcd_vector[ti].tmpmatchwith.clear();
    SKIM3_farcd_vector[ti].tmpmatchwith.reserve(2000);

  }
  FUNCEND();
}

void Skim::farcThreadLoop(const uint32 threadnr)
{
  FUNCSTART("void Skim::threadloop(const uint32 threadnr)");

  // threads need their own try() catch() block

  try {
    CEBUG("Thread: " << threadnr << " starting.\n");

    BUGIFTHROW(threadnr>=SKIM3_farcd_vector.size(),"threadnr>=SKIM3_farcd_vector.size()???");
    farc_threaddata_t & farcd=SKIM3_farcd_vector[threadnr];

    farcd.readhashmatches.clear();
    farcd.singlereadvhraparray.clear();
    farcd.unused_tagmaskvector.clear();
    farcd.tmpmatchwith.clear();

    readid_t dummy=0; // in this version, we do not give back the read id of the adaptor found, but need a variable to call the internal routine

    // we'll jump out with a break;
    while(true){
      {
	boost::mutex::scoped_lock mylock(SKIM3_mutex);
	CEBUG("Thread " << threadnr << " waiting ...\n");
	while(!SKIM3_threadcontrol[threadnr].flag_datavalid
	      && ! SKIM3_threadcontrol[threadnr].flag_endthread){
	  SKIM3_master2slavesignal.wait(mylock);
	}
      }
      if(SKIM3_threadcontrol[threadnr].flag_datavalid){
	CEBUG("Thread " << threadnr << " working on " << SKIM3_threadcontrol[threadnr].from << " to " << SKIM3_threadcontrol[threadnr].to << "\n");

	for(uint32 readi=SKIM3_threadcontrol[threadnr].from; readi<SKIM3_threadcontrol[threadnr].to; ++readi){
	  if(SKIM3_farc_seqtype < 0
	     || SKIM3_farc_searchpool->getRead(readi).getSequencingType() == SKIM3_farc_seqtype){
	    int32 clip=findAdaptorRightClip_internal(SKIM3_farc_searchpool->getRead(readi),SKIM3_farc_minhashes,dummy, farcd);
	    if(clip>=0){
	      boost::mutex::scoped_lock lock(SKIM3_resultfileoutmutex);
	      (*SKIM3_farc_results)[readi]=clip;
	    }
	  }
	}

	boost::mutex::scoped_lock mylock(SKIM3_mutex);
	SKIM3_threadcontrol[threadnr].flag_datavalid=false;

	SKIM3_slave2mastersignal.notify_one();
      }else if(SKIM3_threadcontrol[threadnr].flag_endthread){
	CEBUG("Thread " << threadnr << "  exiting.\n");
	break;
      }
    }

  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

//#define CEBUG(bla)
