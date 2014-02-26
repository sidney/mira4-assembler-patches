/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 and later by Bastien Chevreux
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


#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#include <boost/filesystem.hpp>

#include "errorhandling/errorhandling.H"

#include "util/machineinfo.H"
#include "util/dptools.H"
#include "util/fileanddisk.H"

#include "mira/hashstats.H"
#include "mira/skim.H"
#include "mira/seqtohash.H"
#include "mira/readgrouplib.H"



using namespace std;


//#define CEBUGFLAG

#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

#ifndef PUBLICQUIET
#define CLOCKSTEPS
#endif

#ifdef CLOCKSTEPS
#define TEBUG(bla)   {cout << bla; cout.flush();}
#else
#define TEBUG(bla)
#endif

//#define CEBUG(bla)   {cout << bla; cout.flush();}

#include "boost/format.hpp"
using boost::format;


size_t HashStatistics::HS_numelementsperbuffer=0;


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void HashStatistics::setHashFrequencyRatios(double freqest_minnormal,
					    double freqest_maxnormal,
					    double freqest_repeat,
					    double freqest_heavyrepeat,
					    double freqest_crazyrepeat,
					    uint32 nastyrepeatratio,
					    uint32 nastyrepeatcoverage)
{
  HS_freqest_minnormal=freqest_minnormal;
  HS_freqest_maxnormal=freqest_maxnormal;
  HS_freqest_repeat=freqest_repeat;
  HS_freqest_heavyrepeat=freqest_heavyrepeat;
  HS_freqest_crazyrepeat=freqest_crazyrepeat;
  HS_nastyrepeatratio=nastyrepeatratio;
  HS_nastyrepeatcoverage=nastyrepeatcoverage;
}




/*************************************************************************
 *
 * all steps until a usable hash statistics file is in memory
 * Note: does not delete the final hash statistics file on disk (only the
 *  temporary files)
 *
 * Returns explicitly:
 *   nothing
 *
 * Returns implicitly:
 *   hashstatfilename: name of final hash statistics file on disk
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void HashStatistics::prepareHashStatistics(const string & directory, ReadPool & rp, bool checkusedinassembly, bool onlyagainstrails, bool alsosavesinglehashes, bool fwdandrev, uint32 fwdrevmin, uint8  basesperhash, uint32 millionhashesperbuffer, bool rarekmerearlykill, string & hashstatfilename)
{
  HS_readpoolptr=&rp;
  HS_hs_basesperhash=basesperhash;

  vector<string> hashfilenames;
  vector<size_t> elementsperfile;

  dateStamp(cout);

  cout << "Writing temporary hstat files:\n";
  hashes2disk(hashfilenames,elementsperfile,
	      rp,checkusedinassembly,fwdandrev,fwdrevmin,
	      basesperhash,
	      millionhashesperbuffer,
	      rarekmerearlykill,
	      directory);

  dateStamp(cout);

  cout << "\nAnalysing hstat files:\n";
  size_t numhashstats=
    createHashStatisticsFile(hashstatfilename,
			     hashfilenames,
			     elementsperfile,
			     fwdrevmin,
			     alsosavesinglehashes,
			     directory);

  cout << "\n";

  dateStamp(cout);

  cout << "clean up temporary stat files..."; cout.flush();
  // clean up temporary stat files
  for(uint32 hfni=0; hfni<hashfilenames.size();hfni++){
    removeFile(hashfilenames[hfni],true);
  }
  // but not this one, needed by mirabait
  // TODO: make configurable?
  //removeFile("hashstat.bin",true);

  dateStamp(cout); cout.flush();

  loadHashStatistics(*HS_readpoolptr,hashstatfilename,HS_hs_basesperhash);

  dateStamp(cout);

  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/*

  No Bloom filter. They work exactly as advertised (keep single hashes out of
  the buffers), but are a major disappointment in terms if speed: influence is
  negligible (1 to max 2% faster) for Solexa data. The additional memory is
  not worth it.

 */

#if __GNUC__ >= 3
#define prefetchwrite(p)     __builtin_prefetch((p), 1, 0)
#else
#define prefetchwrite(p)
#endif

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void HashStatistics::hashes2disk(vector<string> & hashfilenames, vector<size_t> & elementsperfile, ReadPool & rp, bool checkusedinassembly, bool fwdandrev, uint32 fwdrevmin, uint8 basesperhash, uint32 millionhashesperbuffer, bool rarekmerearlykill, const string & directory)
{
  FUNCSTART("void HashStatistics::hashes2disk(uint32 basesperhash)");

  const size_t upperbases=2;

  hashfilenames.clear();

  BUGIFTHROW(basesperhash==0,"basesperhash == 0 ???");
  BUGIFTHROW(upperbases>=basesperhash,"upperbases (" << upperbases << ") >=basesperhash " << static_cast<uint16>(basesperhash) << ") ???");

  size_t numfiles=1<<(upperbases*2);
  size_t rightshift=(basesperhash-upperbases)*2;

  CEBUG("bph: " << static_cast<uint16>(basesperhash) << ".\n");
  CEBUG("Must create " << numfiles << " files.\n");
  CEBUG("Rightshift:" << rightshift << '\n');
  CEBUG("sizeof(vhash_t): " << sizeof(vhash_t) << '\n');

  // define how many elements to reserve
  if(HS_numelementsperbuffer==0){
    // default fallback values
    HS_numelementsperbuffer=millionhashesperbuffer*1048576;    // 16m elements, 16b each, 16 buffers == 4 GiB
    if(sizeof(void *)==4){
      // on 32 bit systems, be careful with memory
      HS_numelementsperbuffer=1048576/2;
    }else{
      auto freemem=MachineInfo::getMemAvail();
      cout << "freemem: " << freemem << endl;
      if(freemem>0){
	uint64 tnumhashes=0;
	for(size_t rpi=0; rpi<HS_readpoolptr->size(); ++rpi){
	  auto & actread = HS_readpoolptr->getRead(rpi);
	  if(actread.hasValidData()
	     && !actread.isBackbone()
	     && !actread.isRail()
	     && actread.getLenClippedSeq()>=basesperhash){
	    tnumhashes+=actread.getLenClippedSeq()-basesperhash+1;
	  }
	}
	cout << "TNH: " << tnumhashes << endl;
	const double fillratio=1.5;
	double xmillionelem=static_cast<double>(tnumhashes)/(fillratio*1024*1024*sizeof(hashstat_t));
	cout << "XME 1: " << xmillionelem << endl;
	if(xmillionelem>16.0) {
	  xmillionelem=16;
	}else if(xmillionelem<0.1){
	  xmillionelem=0.1;
	}
	cout << "XME 2: " << xmillionelem << endl;
	HS_numelementsperbuffer=xmillionelem*1024*1024;
	cout << "NEPB 1: " << HS_numelementsperbuffer << endl;
	const uint64 eightgib=8589934592;
	while(HS_numelementsperbuffer>2*1024*1024){    // equivalent to 512 MB
	  uint64 memneeded=HS_numelementsperbuffer*16*sizeof(hashstat_t);
	  if(freemem-memneeded >= eightgib) break;
	  HS_numelementsperbuffer/=2;
	}
	cout << "NEPB 2: " << HS_numelementsperbuffer << endl;
      }
    }
  }

  vector<FILE *> hashfiles(numfiles);
  vector<vector<hashstat_t> > hashfilebuffer(numfiles);
  for(size_t i=0; i<numfiles; ++i){
    string fname=directory+"/stattmp"+str(format("%x") % i )+".bin";
    hashfilenames.push_back(fname);
    hashfiles[i]=fopen(fname.c_str(), "w");
    hashfilebuffer[i].reserve(HS_numelementsperbuffer);
  }

  elementsperfile.clear();
  elementsperfile.resize(numfiles,0);

  hashstat_t tmpdh;
  tmpdh.count=1;
  tmpdh.hasmultipleseqtype=false;

#ifdef CLOCKSTEPS
  timeval tvfill;
  gettimeofday(&tvfill,nullptr);
#endif

  ProgressIndicator<int32> P(0, rp.size());

  // We will use prefetch in the loops below, therefore make sure we do not prefetch memory
  //  which we do not own by making sure the loops flush the buffer before reaching
  //  the capacity of the buffer
  const size_t capacityflush=hashfilebuffer[0].capacity()-2;

  // both memory write prefetches save ~15 to 20% time (well, 1s for 4m Solexa reads at 100bp)

  for(uint32 actreadid=0; actreadid<rp.size(); actreadid++){
    P.progress(actreadid);

    //if(actreadid>100) return;

    Read & actread= rp.getRead(actreadid);

    // Has been taken out as hash statistics now also used for mirabait
    // TODO: check whether this has big influence on "normal" assembly jobs
    //  !!! it has ... for mapping assemblies !!!

    if(!actread.hasValidData()
       || (checkusedinassembly && !actread.isUsedInAssembly())) continue;

    CEBUG("hname: " << actread.getName() << endl);

    const char * namestr=actread.getName().c_str();
    uint32 slen=actread.getLenClippedSeq();

    if(slen<basesperhash) continue;

    tmpdh.seqtype=actread.getSequencingType();
    tmpdh.hasfwd=true;
    tmpdh.hasrev=false;

    size_t hashfilesindex;

    // force this thing with an old cast to be a uint8 pointer
    // static_cast<> just barfs *sigh*
    const uint8 * seq=(const uint8 *) actread.getClippedSeqAsChar();

    SEQTOHASH_LOOPSTART(vhash_t);
    {
      tmpdh.vhash=acthash;
      tmpdh.lowpos=seqi-(basesperhash-1);
      hashfilesindex=tmpdh.vhash>>rightshift;
      CEBUG("Want to write fwd: " << tmpdh << " to " << hashfilesindex << endl);
      BUGIFTHROW(hashfilesindex>=hashfiles.size(),"hashfilesindex>=hashfiles.size() ???");

      if(hashfilebuffer[hashfilesindex].size()==capacityflush){
#ifdef CLOCKSTEPS
	timeval now;
	gettimeofday(&now,nullptr);
#endif
	elementsperfile[hashfilesindex]+=writeCompressedHFB(hashfilebuffer[hashfilesindex],
							    fwdrevmin,
							    hashfiles[hashfilesindex],
							    false,
							    rarekmerearlykill);
#ifdef CLOCKSTEPS
	timeval after,diff;
	gettimeofday(&after,nullptr);
	timersub(&after,&now,&diff);
	timeradd(&diff,&tvfill,&now);
	tvfill=now;
#endif
      }
      hashfilebuffer[hashfilesindex].push_back(tmpdh);
#ifndef _GLIBCXX_DEBUG
      // _GLIBCXX_DEBUG will barf on the [size()+1], but in normal operation we are allowed to
      //   do that as prefetching on non-existent memory is silently ignored
      prefetchwrite(&(hashfilebuffer[hashfilesindex][hashfilebuffer[hashfilesindex].size()+1]));
#endif
    }
    SEQTOHASH_LOOPEND;


    if(fwdandrev){
      tmpdh.hasfwd=false;
      tmpdh.hasrev=true;

      seq=(const uint8 *) actread.getClippedComplementSeqAsChar();

      SEQTOHASH_LOOPSTART(vhash_t);
      {
	tmpdh.vhash=acthash;
	tmpdh.lowpos=slen-seqi+1;
	hashfilesindex=tmpdh.vhash>>rightshift;
	CEBUG("Want to write fwd: " << tmpdh << " to " << hashfilesindex << endl);
	BUGIFTHROW(hashfilesindex>=hashfiles.size(),"hashfilesindex>=hashfiles.size() ???");

	if(hashfilebuffer[hashfilesindex].size()==capacityflush){
#ifdef CLOCKSTEPS
	  timeval now;
	  gettimeofday(&now,nullptr);
#endif
	  elementsperfile[hashfilesindex]+=writeCompressedHFB(hashfilebuffer[hashfilesindex],
							      fwdrevmin,
							      hashfiles[hashfilesindex],
							      false,
							      rarekmerearlykill);
#ifdef CLOCKSTEPS
	  timeval after,diff;
	  gettimeofday(&after,nullptr);
	  timersub(&after,&now,&diff);
	  timeradd(&diff,&tvfill,&now);
	  tvfill=now;
#endif
	}
	hashfilebuffer[hashfilesindex].push_back(tmpdh);
#ifndef _GLIBCXX_DEBUG
	// _GLIBCXX_DEBUG will barf on the [size()+1], but in normal operation we are allowed to
	//   do that as prefetching on non-existent memory is silently ignored
	prefetchwrite(&(hashfilebuffer[hashfilesindex][hashfilebuffer[hashfilesindex].size()+1]));
#endif
      }
      SEQTOHASH_LOOPEND;

    }
  }

  P.finishAtOnce();
  cout << "done\n";

  TEBUG("\nTiming fill HFB: " << diffsuseconds(tvfill) << endl);

  cout << "Flushing buffers to disk:\n";
  P.reset(0,numfiles);
  for(size_t i=0; i<numfiles; ++i){
    P.progress(i);
    elementsperfile[i]+=writeCompressedHFB(hashfilebuffer[i],
					   fwdrevmin,
					   hashfiles[i],
					   true,
					   rarekmerearlykill);
    fclose(hashfiles[i]);
  }
  P.finishAtOnce();
  cout << "done\n";

  //dateStamp(cout);
  //exit(100);

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

size_t HashStatistics::writeCompressedHFB(vector<hashstat_t> & hfb, uint32 fwdrevmin, FILE * fileptr, bool force, bool rarekmerearlykill)
{
  FUNCSTART("size_t HashStatistics::writeCompressedHFB(vector<hashstat_t> & hfb, FILE * fileptr)");
  size_t retvalue=0;
  if(hfb.size()){
    compressHashStatBufferInPlace(hfb, fwdrevmin, !rarekmerearlykill);
    if(force || hfb.size()>=hfb.capacity()*2/3){
      CEBUG("Write buffer " << &hfb << " " << 100*hfb.size()/hfb.capacity() << endl);
      if(myFWrite(&(hfb[0]),sizeof(hashstat_t),hfb.size(),fileptr) != hfb.size()){
	MIRANOTIFY(Notify::FATAL, "Could not write anymore to hash file. Disk full? Changed permissions?");
      }
      retvalue=hfb.size();
      hfb.clear();
    }else{
      CEBUG("No write buffer " << &hfb << " " << 100*hfb.size()/hfb.capacity() << endl);
    }
  }
  return retvalue;
}

///*************************************************************************
// *
// * sorter to sort from low to high in vhash (first) and direction
// *
// *************************************************************************/
//
//inline bool Skim__sortDiskHashComparator_(const hashstat_t & a,
//					  const hashstat_t & b);
//inline bool Skim__sortDiskHashComparator_(const hashstat_t & a,
//					  const hashstat_t & b)
//{
//  if(a.vhash==b.vhash){
//    return a.dir<b.dir;
//  }
//  return a.vhash < b.vhash;
//}

/*************************************************************************
 *
 * sorter to sort from low to high in vhash, then on fwd flag
 *
 *************************************************************************/

inline bool HashStat__sortDiskNewHashComparator_(const hashstat_t & a,
					     const hashstat_t & b);
inline bool HashStat__sortDiskNewHashComparator_(const hashstat_t & a,
					  const hashstat_t & b)
{
  if(a.vhash==b.vhash){
    return a.hasfwd<b.hasfwd;
  }
  return a.vhash < b.vhash;
}

/*************************************************************************
 *
 * sorter to sort from low to high in count
 *
 *************************************************************************/

inline bool HashStat__sortHashStatComparatorByCount_(const hashstat_t & a,
						 const hashstat_t & b);
inline bool HashStat__sortHashStatComparatorByCount_(const hashstat_t & a,
						 const hashstat_t & b)
{
  return a.count < b.count;
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void HashStatistics::compressHashStatBufferInPlace(vector<hashstat_t> & hsb, uint32 fwdrevmin, bool alsosavesinglehashes)
{
  FUNCSTART("void HashStatistics::compressHashStatBufferInPlace(vector<hashstat_t> & hsb)");

  if(hsb.empty()) return;

  CEBUG("Sorting " << hsb.size() << " elements ..."); cout.flush();

#ifdef CLOCKSTEPS
  timeval tv,tvtotal;
  gettimeofday(&tv,nullptr);
  tvtotal=tv;
#endif

  sort(hsb.begin(), hsb.end(), HashStat__sortDiskNewHashComparator_);
  CEBUG("done.\n");
  TEBUG("\nTiming sort HFB: " << diffsuseconds(tv) << endl);

#ifdef CLOCKSTEPS
  gettimeofday(&tv,nullptr);
#endif

  uint32 hasforward=0;
  uint32 hasreverse=0;
  bool hasmultipleseqtype=false;
  bool hasforwardthresholdok=false;
  bool hasreversethresholdok=false;
  bool hasfrthresholdok=false;

  uint8 thisseqtype=0;
  uint32 thishashcounter=0;
  uint16 thislowpos=0;
  hashstat_t tmphs;
  auto srcI=hsb.cbegin();
  auto dstI=hsb.begin();
  // setting this leads the very first iteration of the main loop
  //  to set correct values
  vhash_t thishash=(srcI->vhash)-1;
  for(; srcI!=hsb.cend(); ++srcI){
    if(srcI->vhash != thishash){
      // save only hashes that appeared more than once
      if(thishashcounter>1
	 || (thishashcounter==1 && alsosavesinglehashes)){
	tmphs.vhash=thishash;
	tmphs.count=thishashcounter;
	tmphs.lowpos=thislowpos;
	tmphs.seqtype=thisseqtype;
	tmphs.hasfwd=hasforward>0;
	tmphs.hasrev=hasreverse>0;
	tmphs.hasfwdthresholdok=hasforwardthresholdok | (hasforward>=fwdrevmin);
	tmphs.hasrevthresholdok=hasreversethresholdok | (hasreverse>=fwdrevmin);
	tmphs.hasfwdrevthresholdok=hasfrthresholdok | ((hasforward>=fwdrevmin) & (hasreverse>=fwdrevmin));
	tmphs.hasmultipleseqtype=hasmultipleseqtype;
	CEBUG("Write mid to " << dstI-hsb.begin() << " from " << srcI-hsb.begin() << ": " << tmphs << '\n');

	*dstI=tmphs;
	++dstI;
      }
      hasforward=0;
      hasreverse=0;
      hasmultipleseqtype=false;
      hasforwardthresholdok=false;
      hasreversethresholdok=false;
      hasfrthresholdok=false;
      thishash=srcI->vhash;
      thislowpos=srcI->lowpos;
      thishashcounter=0;
      thisseqtype=srcI->seqtype;
      CEBUG("New vhash: " << hex << thishash << dec << endl);
    }else{
      CEBUG("Existing vhash: " << *srcI << endl);
    }
    thishashcounter+=srcI->count;
    if(srcI->hasfwd) {
      if(srcI->hasfwdthresholdok){
	hasforward+=srcI->count;
      }else{
	++hasforward;
      }
    }
    if(srcI->hasrev) {
      if(srcI->hasrevthresholdok){
	hasreverse+=srcI->count;
      }else{
	++hasreverse;
      }
    }
    if(srcI->lowpos < thislowpos) thislowpos=srcI->lowpos;
    if(srcI->seqtype != thisseqtype
       || srcI->hasmultipleseqtype) hasmultipleseqtype=true;
    hasforwardthresholdok|=srcI->hasfwdthresholdok;
    hasreversethresholdok|=srcI->hasrevthresholdok;
    hasfrthresholdok|=srcI->hasfwdrevthresholdok;

    CEBUG("thc: " << thishashcounter << "\thf: " << hasforward << "\thr: " << hasreverse << "\ttlp: " << thislowpos
	  << "\thmst: " << hasmultipleseqtype
	  << "\thfto: " << hasforwardthresholdok << "\thrto: " << hasreversethresholdok
	  << "\thfrto: " << hasfrthresholdok << endl);
  }

  // we're out of the loop, write last elements if there were any
  if(thishashcounter>1
     || (thishashcounter==1 && alsosavesinglehashes)){
    tmphs.vhash=thishash;
    tmphs.count=thishashcounter;
    tmphs.lowpos=thislowpos;
    tmphs.hasfwdrevthresholdok=(hasforward>=fwdrevmin) & (hasreverse>=fwdrevmin);
    tmphs.hasmultipleseqtype=hasmultipleseqtype;
    CEBUG("Write end to " << dstI-hsb.begin() << " from " << srcI-hsb.begin() << ": " << tmphs << '\n');
    *dstI=tmphs;
    ++dstI;
  }

  TEBUG("Timing compress HFB: " << diffsuseconds(tv) << endl);
  TEBUG("Timing compressHashStatBufferInPlace: " << diffsuseconds(tvtotal) << endl);

  hsb.resize(dstI-hsb.begin());

#ifndef PUBLICQUIET
  {
    uint64 numsingle=0;
    uint64 nummulti=0;
    auto eI=hsb.cend();
    for(auto iI=hsb.cbegin(); iI!=eI; ++iI){
      if(iI->count==1){
	++numsingle;
      }else{
	++nummulti;
      }
    }
    cout << "\nnumsingle: " << numsingle << "\nnummulti:  " << nummulti << endl;
  }
#endif

  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * sorts every hashfile and writes a hash statistics file
 *
 * TODO: delete single hash files and clear vectors of filenames and sizes
 *
 * returns:
 *  - by value: elements in hash statistics file
 *  - name of file in the call by reference variable
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

size_t HashStatistics::createHashStatisticsFile(string & hashstatfilename, vector<string> & hashfilenames, vector<size_t> & elementsperfile, uint32 fwdrevmin, bool alsosavesinglehashes, const string & directory)
{
  FUNCSTART("size_t HashStatistics::createHashStatisticsFile(string & hashstatfilename, vector<string> & hashfilenames, vector<size_t> & elementsperfile, ReadPool & rp, bool onlyagainstrails, const string & directory)");

  FILE * fout;
  hashstatfilename=directory+"/hashstat.bin";
  fout=fopen(hashstatfilename.c_str(), "w");

  size_t maxelementsperfile=0;

  for(size_t fi=0; fi< elementsperfile.size(); fi++){
    maxelementsperfile=max(maxelementsperfile,elementsperfile[fi]);
  }

  CEBUG("Max elements per file: " << maxelementsperfile << '\n');

  vector<hashstat_t> hashpool;
  hashpool.reserve(maxelementsperfile+10);

  size_t numhashstats=0;

  ProgressIndicator<int32> P(0, static_cast<int32>(elementsperfile.size()));

  for(size_t fi=0; fi< elementsperfile.size(); fi++){
    P.increaseprogress();

    CEBUG("Loading " << hashfilenames[fi] << endl);
    CEBUG("elements in file: " << elementsperfile[fi] << endl);

    if(elementsperfile[fi]==0) continue;
    hashpool.clear();
    hashpool.resize(elementsperfile[fi]);

    FILE * fin;
    fin=fopen(hashfilenames[fi].c_str(), "r");
    if(myFRead(&hashpool[0],sizeof(hashstat_t),elementsperfile[fi],fin) != elementsperfile[fi]) {
      MIRANOTIFY(Notify::FATAL, "Expected to read " << elementsperfile[fi] << " elements in file " << hashfilenames[fi] << " but read less. Was the file deleted? Disk full?");
    }
    fclose(fin);

    //for(size_t i=0; i<hashpool.size(); i++){
    //  CEBUG(hashpool[i] << '\n');
    //}

    for(size_t i=0; i<hashpool.size(); i++){
      CEBUG(hashpool[i] << '\n');
    }

    compressHashStatBufferInPlace(hashpool,fwdrevmin,alsosavesinglehashes);
    numhashstats+=hashpool.size();

    if(myFWrite(&hashpool[0],sizeof(hashstat_t),hashpool.size(),fout) != hashpool.size()){
      MIRANOTIFY(Notify::FATAL, "Expected to write " << hashpool.size() << " elements in file " << hashstatfilename << " but could not. Was the file deleted? Disk full?");
    }
  }

  fclose(fout);

  P.finishAtOnce();

  FUNCEND();
  return numhashstats;
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * Needs: Name of the hashstat filename
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void HashStatistics::loadHashStatistics(ReadPool & rp, string & hashstatfilename, uint8 basesperhash)
{
  FUNCSTART("void HashStatistics::loadHashStatisticsFile(string & hashstatfilename, uint8 basesperhash)");

  HS_readpoolptr=&rp;
  HS_hs_basesperhash=basesperhash;

  HS_hs_hashstats.clear();
  HS_hs_hsshortcuts.clear();

  HS_avg_freq_corrected=0;
  HS_avg_freq_raw=0;
  HS_avg_freq_taken=0;

  BUGIFTHROW(!fileExists(hashstatfilename),"No hash statistics file " << hashstatfilename << " to load data from?");
  auto fs=boost::filesystem::file_size(hashstatfilename);
  BUGIFTHROW(fs==0,"Empty file " << hashstatfilename << " ?");
  if(fs%sizeof(hashstat_t)){
    MIRANOTIFY(Notify::FATAL, "File probably not a hash stat: " << hashstatfilename);
  }

  HS_hs_hashstats.resize(fs/sizeof(hashstat_t));

  CEBUG("Loading hash stats: " << HS_hs_hashstats.size() << endl);

  FILE * fin;
  fin=fopen(hashstatfilename.c_str(), "r");
  size_t numread=myFRead(&HS_hs_hashstats[0],sizeof(hashstat_t),HS_hs_hashstats.size(),fin);
  if(numread != HS_hs_hashstats.size()){
    MIRANOTIFY(Notify::FATAL, "Expected to read " << HS_hs_hashstats.size() << " elements in hashfile " << hashstatfilename << " but read less (" << numread << "). Was the file deleted? Disk full?");
  }
  fclose(fin);

  if(HS_hs_hashstats.begin() == HS_hs_hashstats.end()) return;

  if(HS_logflag_hashcount){
    string logfile=hashstatfilename+".shouldneverbeseen.hashcount.usort";
    ofstream fout;
    fout.open(logfile.c_str(), ios::out);

    for(auto & hs : HS_hs_hashstats){
      if(hs.hasfwdrevthresholdok) fout << hs.count << "\n";
    }
  }

  sort(HS_hs_hashstats.begin(),HS_hs_hashstats.end(),HashStat__sortHashStatComparatorByCount_);

  if(HS_logflag_hashcount){
    string logfile=hashstatfilename+".shouldneverbeseen.hashcount.sort";
    ofstream fout;
    fout.open(logfile.c_str(), ios::out);

    for(auto & hs : HS_hs_hashstats){
      if(hs.hasfwdrevthresholdok) fout << hs.count << "\n";
    }
  }

  // do this before makeHashStatArrayShortcuts() as mHSAS() resorts
  //  HS_hs_hashstats, but not by count!
  calcAvgHashFreq();

  makeHashStatArrayShortcuts();
}


void HashStatistics::calcAvgHashFreq()
{
  HS_avg_freq_raw=calcMidHashStatIndex(HS_hs_hashstats,0);
  HS_avg_freq_corrected=HS_avg_freq_raw;

  cout << "Raw MHI: " << HS_avg_freq_raw << endl;
  cout << "Raw avg. freq. : " << HS_hs_hashstats[HS_avg_freq_raw].count << endl;

  auto hsthreshold=HS_hs_hashstats.size()-HS_hs_hashstats.size()/10;
  cout << "HSS " << HS_hs_hashstats.size() << "\tHSST: " << hsthreshold << endl;

  //// if mh index is in last 10 % of the hashstats, we have a pretty skewed
  ////  distribution. In that case, recalc without last 10%
  //// TODO: check whether 40 or 50% wouldn't be better.
  if(HS_avg_freq_corrected >= hsthreshold){
    HS_avg_freq_corrected=calcMidHashStatIndex(HS_hs_hashstats,10);
    cout << "Corrected MHI: " << HS_avg_freq_corrected << endl;
    cout << "Corrected avg. freq. : " << HS_hs_hashstats[HS_avg_freq_corrected].count << endl;
    cout << "HSS " << HS_hs_hashstats.size() << "\tHSST: " << (HS_hs_hashstats.size()-HS_hs_hashstats.size()/10) << endl;
  }

//// Not sure this is good (also: sometimes segfaults???)
//  size_t mhi=HS_hs_hashstats.size();
//  uint8 dontcarepercent=0;
//
//  // if mh index is in last 10 % of the hashstats, we have a pretty skewed
//  //  distribution. In that case, recalc without last dontcarepercent%
//  for(;dontcarepercent<100 ; ++dontcarepercent){
//    mhi=calcMidHashStatIndex(HS_hs_hashstats,dontcarepercent);
//    cout << "DCP: " << static_cast<uint16>(dontcarepercent) << endl;
//    cout << "Raw MHI: " << mhi << endl;
//    cout << "Raw avg. freq. : " << HS_hs_hashstats[mhi].count << endl;
//    cout << "HSS " << HS_hs_hashstats.size() << "\tHSST: " << (HS_hs_hashstats.size()-HS_hs_hashstats.size()/10) << endl;
//    if(mhi < (HS_hs_hashstats.size()-HS_hs_hashstats.size()/10)){
//      break;
//    }
//  }


  HS_avg_freq_corrected=HS_hs_hashstats[HS_avg_freq_corrected].count;
  HS_avg_freq_raw=HS_hs_hashstats[HS_avg_freq_raw].count;

  HS_avg_freq_taken=HS_avg_freq_corrected;
  if(HS_avg_freq_taken < HS_avg_freq_min){
    HS_avg_freq_taken=HS_avg_freq_min;
    cout << "Forced avg. freq: " << HS_avg_freq_taken << endl;
  }


  FUNCEND();
  return;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}
size_t HashStatistics::calcMidHashStatIndex(const vector<hashstat_t> & hashstats, size_t dontcarepercent)
{
  FUNCSTART("size_t HashStatistics::calcMidHashStatIndex(const vector<hashstat_t> & hashstats, size_t dontcarepercent)");

  if(hashstats.empty()) return 0;

  size_t firsti=0;
  size_t lasti=hashstats.size();
  if(dontcarepercent){
    firsti=hashstats.size()*dontcarepercent/100;
    lasti-=hashstats.size()*dontcarepercent/100;
  }else{
    // 5% default
    firsti=hashstats.size()/20;
    lasti-=hashstats.size()/20;
  }

  size_t sumhashcounts=0;
  uint32 oldhashcount=hashstats[0].count-1;
  size_t oldsumhashcounts=0;
  for(size_t i=firsti; i<lasti; i++){
    if(hashstats[i].count != oldhashcount){
      BUGIFTHROW(oldhashcount>hashstats[i].count,"haststat array not sorted by count???");
      oldhashcount=hashstats[i].count;
      CEBUG("count: " << oldhashcount << "\tsumhash: " << sumhashcounts << "\tdiff: " << sumhashcounts-oldsumhashcounts << endl);
      oldsumhashcounts=sumhashcounts;
    }
    if(hashstats[i].hasfwdrevthresholdok) sumhashcounts+=hashstats[i].count;
  }
  CEBUG("count: " << oldhashcount << "\tsumhash: " << sumhashcounts << endl);

  // Hmmm, pathological case. Maybe all reads were in the same direction.
  //  simply recalc without the "has fwd/rev" clause
  bool dontusefwdrev=false;
  if(sumhashcounts==0){
    dontusefwdrev=true;
    for(size_t i=firsti; i<lasti; i++){
      sumhashcounts+=hashstats[i].count;
    }
    CEBUG("recalc sumhash: " << sumhashcounts << endl);
  }

  size_t midhashstats=sumhashcounts/2;

  CEBUG("midhashstats: " << midhashstats << endl);

  sumhashcounts=0;
  for(size_t i=firsti; i<lasti; i++){
    if(dontusefwdrev || hashstats[i].hasfwdrevthresholdok) sumhashcounts+=hashstats[i].count;
    if(sumhashcounts>midhashstats) {
      return i;
    }
  }

  FUNCEND();

  return 0;
}
//#define CEBUG(bla)





/*************************************************************************
 *
 * Needs:
 *  - the hash statistics vector (sorted by count)
 *
 *************************************************************************/

void HashStatistics::showHashStatisticsInfo()
{
  FUNCSTART("void HashStatistics::showHashStatisticsInfo()");

  cout << "Hash statistics:\n"
       << "=========================================================\n"
       << "Measured avg. raw frequency coverage: " << HS_avg_freq_raw << endl
       << "Corrected avg. raw frequency coverage: " << HS_avg_freq_corrected;

  if(HS_avg_freq_raw!=HS_avg_freq_corrected){
    cout << "\tSKEWED DISTRIBUTION!";
  }
  cout << '\n';

  if(HS_avg_freq_corrected<HS_avg_freq_min){
    cout << "Forced minimum average frequency: " << HS_avg_freq_min << endl;
  }

  cout << "\nFinal average frequency: " << HS_avg_freq_taken << endl;


  cout << "\nDeduced thresholds:\n"
       << "-------------------"
       << "\nMin normal cov: " << HS_freqest_minnormal*HS_avg_freq_taken
       << "\nMax normal cov: " << HS_freqest_maxnormal*HS_avg_freq_taken
       << "\nRepeat cov: " << HS_freqest_repeat*HS_avg_freq_taken
       << "\nHeavy cov: " << HS_freqest_heavyrepeat*HS_avg_freq_taken
       << "\nCrazy cov: " << HS_freqest_crazyrepeat*HS_avg_freq_taken
       << "\nMask cov: " << HS_nastyrepeatratio*HS_avg_freq_taken
       << "\n\nRepeat ratio histogram:\n"
       << "-----------------------"
       << endl;

  vector<size_t> ratiocounts;
  ratiocounts.reserve(8192);
  for(size_t i=0; i<HS_hs_hashstats.size(); i++){
    uint32 rci=static_cast<uint32>((static_cast<double>(HS_hs_hashstats[i].count) / HS_avg_freq_taken) + 0.5);
    if(rci>=ratiocounts.size()){
      ratiocounts.resize(rci+1,0);
    }
    ratiocounts[rci]++;
  }

  for(size_t i=0; i<ratiocounts.size(); i++){
    if(ratiocounts[i]) cout << i << '\t' << ratiocounts[i] << endl;
  }

  cout << "=========================================================\n\n";

  FUNCEND();

  return;
}



/*************************************************************************
 *
 * sorter to sort from low to high, but lower 24bit grouped
 *
 *
 *************************************************************************/

inline bool HashStat__sortHashStatComparatorByLow24bit_(const hashstat_t & a,
						    const hashstat_t & b)
{
  if((a.vhash & SKIM3_MAXVHASHMASK) != (b.vhash & SKIM3_MAXVHASHMASK)) {
    return (a.vhash & SKIM3_MAXVHASHMASK) < (b.vhash & SKIM3_MAXVHASHMASK);
  }
  return a.vhash < b.vhash;
}

/*************************************************************************
 *
 * sorter to sort from low to high, using a definable mask
 *
 *
 *************************************************************************/

vhash_t HashStatistics__vhashmask;

inline bool HashStatistics__sortHashStatComparatorByMask_(const hashstat_t & a,
							  const hashstat_t & b)
{
  if((a.vhash & HashStatistics__vhashmask) != (b.vhash & HashStatistics__vhashmask)) {
    return (a.vhash & HashStatistics__vhashmask) < (b.vhash & HashStatistics__vhashmask);
  }
  return a.vhash < b.vhash;
}

/*************************************************************************
 *
 * ideally, instead of using mincount, it should be mincount per
 *  direction. NHashStatistics will do that ...
 *
 *************************************************************************/

string laberbla;

void HashStatistics::calcKMerForks(uint32 mincount)
{
  if(HS_hs_hashstats.empty()) return;

  for(auto & hse : HS_hs_hashstats) hse.iskmerfork=false;

  if(HS_hs_basesperhash<17) return;

  auto rollbases=HS_hs_basesperhash-1;
  HashStatistics__vhashmask=1;
  /* *grml* undefined behaviour of left shift for 64 shifts in a 64 bit type makes this cludge necessary */
  /* the same for 32 shift in 32 bit types etc.pp */
  if(HS_hs_basesperhash>=sizeof(vhash_t)*4){
    HashStatistics__vhashmask=0;
  }else{
    HashStatistics__vhashmask<<=(rollbases*2);
  }

  // vhashmask is now, e.g. for bph=31, 00010000000....
  --HashStatistics__vhashmask;
  // vhashmask is now, e.g. for bph=31, 000011111....

  // calc the status on ?..............
  sort(HS_hs_hashstats.begin(), HS_hs_hashstats.end(), HashStatistics__sortHashStatComparatorByMask_);
  laberbla="fwd ";
  ckmf_helper(HashStatistics__vhashmask,mincount);

  // calc the status on ..............?
  HashStatistics__vhashmask<<=2;
  // vhashmask is now, e.g. for bph=31, 0011111....00
  sort(HS_hs_hashstats.begin(), HS_hs_hashstats.end(),  HashStatistics__sortHashStatComparatorByMask_);

  laberbla="rev ";
  ckmf_helper(HashStatistics__vhashmask,mincount);

  sort(HS_hs_hashstats.begin(), HS_hs_hashstats.end(), HashStat__sortHashStatComparatorByLow24bit_);
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void HashStatistics::ckmf_helper(vhash_t hashmask, uint32 mincount)
{
  auto hsI=HS_hs_hashstats.begin();
  auto hsJ=hsI+1;

  CEBUG("KMER DUMP:\n");

  for(; hsJ!=HS_hs_hashstats.end(); ++hsI,++hsJ){
    if(hsI->hasfwdrevthresholdok
       && hsJ->hasfwdrevthresholdok
       && (hsI->vhash&hashmask) == (hsJ->vhash&hashmask)
       && hsI->vhash != hsJ->vhash
       && hsI->count >= mincount
       && hsJ->count >= mincount){
      hsI->iskmerfork=true;
      hsJ->iskmerfork=true;
    }
    CEBUG(laberbla << NHashStatistics::hash2string(hsI->vhash,HS_hs_basesperhash) << ' ' << *hsI << '\n');
    CEBUG(laberbla << NHashStatistics::hash2string(hsJ->vhash,HS_hs_basesperhash) << ' ' << *hsJ << '\n');
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * needs:
 *  - hashstats filled with entries (can be unsorted, will be re-sorted
 *    anyway)
 *
 * returns:
 *  - hashstats array sorted by low 24 bit (low to high), then by vhash
 *  - hsshortcuts_begin and ..._end pointing to start and end of each
 *    low 24 bit group of same value
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void HashStatistics::makeHashStatArrayShortcuts()
{
  FUNCSTART("void HashStatistics::makeHashStatArrayShortcuts()");

  CEBUG("makeHashStatArrayShortcuts: basesperhash: " << static_cast<uint16>(basesperhash) << "\n");

  BUGIFTHROW(HS_hs_basesperhash==0, "HS_hs_basesperhash == 0 ???");

  for(size_t i=0; i<HS_hs_hashstats.size(); i++){
    CEBUG(HS_hs_hashstats[i] << '\n');
  }

  sort(HS_hs_hashstats.begin(), HS_hs_hashstats.end(), HashStat__sortHashStatComparatorByLow24bit_);

  HS_hs_hsshortcuts.clear();
  {
    hsvbendit_t tmpb;
    tmpb.b=HS_hs_hashstats.end();
    tmpb.e=HS_hs_hashstats.end();
    HS_hs_hsshortcuts.resize(
      1<<(min(static_cast<uint8>(12),HS_hs_basesperhash)*2),
      tmpb
      );
  }

  CEBUG("HS_hs_hsshortcuts.size(): " << HS_hs_hsshortcuts.size() << endl);

  auto hsI=HS_hs_hashstats.begin();
  if(hsI==HS_hs_hashstats.end()) return;


  vhash_t acthash= (hsI->vhash & SKIM3_MAXVHASHMASK);
  while(hsI != HS_hs_hashstats.end()){
    CEBUG("begin " << hex << acthash << dec << " is: " << *hsI << endl);
    HS_hs_hsshortcuts[acthash].b=hsI;
    for(;(hsI != HS_hs_hashstats.end()) && ((hsI->vhash & SKIM3_MAXVHASHMASK) == acthash); hsI++) {
      CEBUG("INC\n")
    }
    CEBUG("end " << hex << acthash << dec << " is: " << *hsI << endl);
    HS_hs_hsshortcuts[acthash].e=hsI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << HS_hs_hsshortcuts_end[acthash]-HS_hs_hsshortcuts_begin[acthash] << '\n';
    if(hsI != HS_hs_hashstats.end()) acthash= hsI->vhash & SKIM3_MAXVHASHMASK;
  }

  FUNCEND();
}
//#define CEBUG(bla)


/*************************************************************************
 *
 * comparator for lower_bound below
 *
 *************************************************************************/

inline bool HashStat__compareHashStatHashElem_(const hashstat_t & a,
					       const hashstat_t & b);
inline bool HashStat__compareHashStatHashElem_(const hashstat_t & a,
					       const hashstat_t & b)
{
  return a.vhash < b.vhash;
}

//#define CEBUG(bla)   {cout << bla; cout.flush();}



/*************************************************************************
 *
 *
 *
 *
 *
 *************************************************************************/

 //#define CEBUG(bla)   {cout << bla; cout.flush();}
 //void HashStatistics::assignReadBaseStatistics_MultiThread(uint32 numthreads, ReadPool & rp, size_t avghashcov, vector<hashstat_t> & hashstats, const uint8 basesperhash, vector<vector<hashstat_t>::const_iterator > & hsshortcuts_begin, vector<vector<hashstat_t>::const_iterator > & hsshortcuts_end, bool masknastyrepeats, vector<uint32> & rarekmermasking)
void HashStatistics::assignReadBaseStatistics_MultiThread(uint32 numthreads, bool masknastyrepeats, vector<uint32> & rarekmermasking, uint32 mincountkmerforks)
{
  FUNCSTART("void HashStatistics::assignReadBaseStatistics_MultiThread(uint32 numthreads, bool masknastyrepeats, vector<uint32> & rarekmermasking)");

  if(mincountkmerforks>0) calcKMerForks(mincountkmerforks);

  if(rarekmermasking.size() != ReadGroupLib::getNumSequencingTypes()){
    rarekmermasking.clear();
    rarekmermasking.resize(ReadGroupLib::getNumSequencingTypes(),0);
  }

  arbs_threadsharecontrol_t atsc;

  atsc.from=0;
  atsc.to=HS_readpoolptr->size();
  atsc.todo=0;
  atsc.done=0;
  atsc.stepping=1000;

  // TODO: unneeded now as working on HS_* variables, reorganise
  // vvvvvvvvvvvvvvv
  atsc.rpptr=HS_readpoolptr;
  atsc.avghashcov=HS_avg_freq_taken;
  atsc.hashstatsptr=&HS_hs_hashstats;
  atsc.basesperhash=HS_hs_basesperhash;
  atsc.hsscptr=&HS_hs_hsshortcuts;
  // ^^^^^^^^^

  atsc.masknastyrepeats=masknastyrepeats;
  atsc.rarekmermaskingptr=&rarekmermasking;

  CEBUG("minnormalhashcov: " << atsc.avghashcov << endl);


  //uint32 numthreads=8;
  boost::thread_group workerthreads;
  for(uint32 ti=0; ti<numthreads;++ti){
    workerthreads.create_thread(boost::bind(&HashStatistics::priv_arb_thread, this, ti, &atsc));
  }

  ProgressIndicator<int64> pi(0,HS_readpoolptr->size());
  while(atsc.done!=HS_readpoolptr->size()){
    pi.progress(atsc.done);
    sleep(1);
  }
  pi.finishAtOnce(cout);

  // they normally should all have exited at this point, but be nice and play by the rules
  workerthreads.join_all();

}
//#define CEBUG(bla)

void HashStatistics::priv_arb_thread(uint32 threadnum, arbs_threadsharecontrol_t * tscptr)
{
  FUNCSTART("");

  try{
    int32 from;
    int32 to;
    while(true){
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	if(tscptr->todo >= tscptr->to) break;
	from=tscptr->todo;
	tscptr->todo+=tscptr->stepping;
	if(tscptr->todo > tscptr->to) tscptr->todo = tscptr->to;
	to=tscptr->todo;
      }
      priv_arb_DoStuff(
	*(tscptr->rpptr),
	tscptr->avghashcov,
	*(tscptr->hashstatsptr),
	tscptr->basesperhash,
	*(tscptr->hsscptr),
	tscptr->masknastyrepeats,
	*(tscptr->rarekmermaskingptr),
	from,
	to
	);
      {
	boost::mutex::scoped_lock lock(tscptr->accessmutex);
	tscptr->done+=tscptr->stepping;
	if(tscptr->done > tscptr->to) tscptr->done=tscptr->to;
      }
    }
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
}


//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {if(docebug) {cout << bla; cout.flush();}}
void HashStatistics::priv_arb_DoStuff(ReadPool & rp, size_t avghashcov, vector<hashstat_t> & hashstats, const uint8 basesperhash, vector<hsvbendit_t> & hsshortcuts, bool masknastyrepeats, vector<uint32> & rarekmermasking, int32 fromid, int32 toid)
{
  FUNCSTART("HashStatistics::priv_arb_DoStuff(ReadPool & rp, size_t avghashcov, vector<hashstat_t> & hashstats, const uint8 basesperhash, vector<hsvbendit_t> & hsshortcuts, bool masknastyrepeats, vector<uint32> & rarekmermasking, int32 fromid, int32 toid)");

  //bool docebug=true;

  auto minnormalhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_minnormal);
  auto maxnormalhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_maxnormal);
  auto repeathashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_repeat);
  auto heavyrepthashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_heavyrepeat);
  auto crazyrepthashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_freqest_crazyrepeat);
  auto maskhashcov=static_cast<uint32>(static_cast<double>(avghashcov)*HS_nastyrepeatratio);

  if(HS_nastyrepeatcoverage>0 && HS_nastyrepeatcoverage<maskhashcov){
    maskhashcov=HS_nastyrepeatcoverage;
  }

  vector<vhrap_t> singlereadvhraparray;
  singlereadvhraparray.reserve(10000);

  // we will not use a mask, but
  //  we need to supply an empty one anyway
  vector<uint8> tagmaskvector;

  // stores in each read whether the given hash frequency was seen
  vector<uint8> hasfrequency(8);

  vector<uint8> mcmask;
  mcmask.reserve(10000);

  multitag_t tmpmt(Read::REA_defaulttag_MNRr);

  CEBUG("dostuff from " << fromid << " to " << toid << endl);

  for(int32 actreadid=fromid; actreadid<toid; ++actreadid){
    //if(actreadid>100) return;

    Read & actread= rp.getRead(actreadid);

    CEBUG("dsloop " << fromid << " " << actread.getName() << " " << actread.getLenClippedSeq() << endl);

    // get rid of old values
    actread.clearAllBPosHashStats();
    actread.setHasFreqAvg(false);
    actread.setHasFreqRept(false);
    actread.setHasKMerFork(false);
    actread.deleteTag(tmpmt.identifier);

    // whatever happens: this read was looked upon by this routine, so technically we "have" base hashstats
    actread.setHasBaseHashStats(true);


    if(!actread.hasValidData()
      || !actread.isUsedInAssembly()) continue;

//#define CEBUG(bla)   {if(cebugok) cout << bla; cout.flush();}
//    bool cebugok=false;
//    if(actread.getName()=="E0K6C4E01CTNQI") cebugok=true;

    uint32 slen=actread.getLenClippedSeq();

    if(slen<basesperhash) continue;

    mcmask.clear();
    mcmask.resize(actread.getLenSeq(),0);

    hasfrequency.clear();
    hasfrequency.resize(8,0);

    Read::setCoutType(Read::AS_TEXT);
    CEBUG("### Before ...\n" << actread << endl);

    singlereadvhraparray.resize(slen);
    tagmaskvector.resize(slen,0);

    vector<vhrap_t>::iterator srvaI=singlereadvhraparray.begin();

    vector<Read::bposhashstat_t> & bposhashstats=const_cast<vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
    uint32 hashesmade;

    {
      int32 bfpos=actread.calcClippedPos2RawPos(0);
      int32 bfposinc=1;

      hashesmade=Skim::transformSeqToVariableHash(
	actreadid,
	actread,
	actread.getClippedSeqAsChar(),
	slen,
	basesperhash,
	srvaI,
	false,
	1,
	tagmaskvector,
	bposhashstats,
	bfpos,
	bfposinc
	);
    }
    singlereadvhraparray.resize(hashesmade);

    CEBUG("hashesmade: " << hashesmade << endl);
    CEBUG("maskhashcov: " << maskhashcov << endl);

    vector<hashstat_t>::const_iterator lowerbound;

    vector<hashstat_t>::const_iterator hssearchI;
    srvaI=singlereadvhraparray.begin();

    uint32 rarekmercount=rarekmermasking[actread.getSequencingType()];
    //uint32 rarekmercount=30;

    int32 bfpos1,bfpos2;
    hashstat_t hstmp;
    bool foundit=false;
    bool haskmerfork=false;
    for(; srvaI != singlereadvhraparray.end(); srvaI++){
      CEBUG(*srvaI << '\n');

      lowerbound=hsshortcuts[srvaI->vhash & SKIM3_MAXVHASHMASK].b;
      foundit=false;

      // "HS_empty_vector_hashstat_t.end()" is the "nullptr" replacement
      if(hashstats.end() != lowerbound){
	if(basesperhash>12){
	  // with more than 12 bases in a hash, the array is subdivided
	  hstmp.vhash=srvaI->vhash;
	  hssearchI=lower_bound(lowerbound,
				hsshortcuts[srvaI->vhash & SKIM3_MAXVHASHMASK].e,
				hstmp,
				HashStat__compareHashStatHashElem_);
	  if(hssearchI != hashstats.end()
	     && hssearchI->vhash == srvaI->vhash) foundit=true;
	}else{
	  hssearchI=lowerbound;
	  foundit=true;
	}
      }else{
	CEBUG("---------- NO LB HIT??? -------\n");
      }

      if(foundit) {
	CEBUG("VHRAP: " << *srvaI << '\n');
	CEBUG("HashStat: " << *hssearchI << '\n');
	CEBUG("srvaI->hashpos: " << srvaI->hashpos << '\n');

	bfpos1=actread.calcClippedPos2RawPos(srvaI->hashpos-(basesperhash-1));
	bfpos2=bfpos1+basesperhash-1;

	CEBUG("b bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	CEBUG("b bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	bposhashstats[bfpos1].fwd.setValid();
	bposhashstats[bfpos2].rev.setValid();

	if(hssearchI->hasfwdrevthresholdok) {
	  //bhs|=Read::BFLAGS_CONFIRMED_FWDREV;
	  CEBUG("Set ConfFWDREV\n");
	  bposhashstats[bfpos1].fwd.setConfirmedFwdRev();
	  bposhashstats[bfpos2].rev.setConfirmedFwdRev();
	}
	if(hssearchI->iskmerfork) {
	  haskmerfork=true;
	  bposhashstats[bfpos1].fwd.setKMerFork();
	  bposhashstats[bfpos2].rev.setKMerFork();
	}
	if(hssearchI->lowpos<=4){
	  //bhs|=Read::BFLAGS_SEENATLOWPOS;
	  CEBUG("Set SeenAtLowPos\n");
	  bposhashstats[bfpos1].fwd.setSeenAtLowPos();
	  bposhashstats[bfpos2].rev.setSeenAtLowPos();
	}
	if(hssearchI->hasmultipleseqtype){
	  //bhs|=Read::BFLAGS_CONFIRMED_MULTIPLESEQTYPE;
	  CEBUG("Set ConfMultSeqType\n");
	  bposhashstats[bfpos1].fwd.setConfirmedMultipleSeqType();
	  bposhashstats[bfpos2].rev.setConfirmedMultipleSeqType();
	}
	uint8 frequency=2;
	if(hssearchI->count == 1){
	  frequency=1;
	}else if(hssearchI->count <= rarekmercount ){
	  // maybe additional checks ... ?
	  frequency=1;
	}else if(hssearchI->count<minnormalhashcov) {
	  frequency=2;
	}else if(hssearchI->count>=minnormalhashcov
	   && hssearchI->count<=maxnormalhashcov) {
	  frequency=3;
	  //}else if(hssearchI->count > minnormalhashcov*20){
	}else if(hssearchI->count > crazyrepthashcov){
	  frequency=7;
	}else if(hssearchI->count > heavyrepthashcov){
	  frequency=6;
	}else if(hssearchI->count>=repeathashcov){
	  frequency=5;
	}else{
	  frequency=4;
	}
	CEBUG("Set frequency: " << static_cast<uint16>(frequency) << endl);

	if(maskhashcov>0 && hssearchI->count>=maskhashcov){
	  CEBUG("mcmask " << bfpos1 << "\t" << bfpos1+basesperhash << '\n');
	  for(uint32 j=0; j<basesperhash; j++){
	    mcmask[bfpos1+j]=1;
	  }
	}

	CEBUG("a1 bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	CEBUG("a1 bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	bposhashstats[bfpos1].fwd.setFrequency(frequency);
	bposhashstats[bfpos2].rev.setFrequency(frequency);

	CEBUG("a2 bfpos1: " << bfpos1 << '\t' << bposhashstats[bfpos1] << endl);
	CEBUG("a2 bfpos2: " << bfpos2 << '\t' << bposhashstats[bfpos2] << endl);

	hasfrequency[frequency]=1;

	//cout.flush();
	//actread.setBaseFlagsInClippedSequence(bhs,
	//				      srvaI->hashpos-(basesperhash-1),
	//				      basesperhash);
	//actread.setHasBaseFlags(true);
      }
    }

    if(hasfrequency[3]){
      actread.setHasFreqAvg(true);
    }
    if(hasfrequency[5] || hasfrequency[6] || hasfrequency[7]){
      actread.setHasFreqRept(true);
    }
    actread.setHasKMerFork(haskmerfork);

    //Read::setCoutType(Read::AS_TEXT);
    //CEBUG("### After ...\n" << actread << endl);


    // BaCh 07.04.2009 Bad Idea!!!
    // BaCh 12.07.2009 Why? Forgot ... :-(
    //// the fwd/rev of a read now looks like this (e.g.)
    //// (for better viewing dot == 0)
    ////
    //// f   ..........2222222233333....355555....................
    //// r   ................2222222....33333355555...............
    ////
    //// in dubio pro reo and to allow for potential matches,
    //// do this:
    ////
    //// f   ..........2222222233333....355555->..................
    //// r   ..............<-2222222....33333355555...............
    ////
    //// so that this
    ////
    //// f   ..........2222222233333....35555555555...............
    //// r   ..........2222222222222....33333355555...............
    ////
    //// is generated
    ////
    ////
    //
    //{
    //  uint32 bfposi=0;
    //  for(; bfposi<bposhashstats.size() && bposhashstats[bfposi].fwd.getFrequency()==0; bfposi++) {};
    //  uint32 bfpose=bfposi;
    //  for(; bfpose<bposhashstats.size() && bposhashstats[bfpose].rev.getFrequency()==0; bfpose++) {};
    //  if(bfposi<bposhashstats.size() && bfpose<bposhashstats.size()){
    //	for(uint32 i=bfposi; i<bfpose; i++){
    //	  bposhashstats[i].fwd=bposhashstats[bfpose].rev;
    //	}
    //  }
    //
    //  bfposi=bposhashstats.size()-1;
    //  for(; bfposi>0 && bposhashstats[bfposi].rev.getFrequency()==0; bfposi--) {};
    //  bfpose=bfposi;
    //  for(; bfpose>0 && bposhashstats[bfpose].fwd.getFrequency()==0; bfpose--) {};
    //  if(bfposi>0){
    //	for(uint32 i=bfposi; i>bfpose; i--){
    //	  bposhashstats[i].fwd=bposhashstats[bfpose].rev;
    //	}
    //  }
    //}


    // go through multicopy array and set MNRr tags for
    //  consecutive positions in read tagged as multicopy
    if(masknastyrepeats){
      bool inrun=false;
      uint32 runstart=0;
      uint32 pos=0;
      for(; pos<mcmask.size(); pos++){
	CEBUG("pos: " << pos << '\t' << static_cast<uint16>(mcmask[pos]) << '\t' << inrun << '\n');
	if(mcmask[pos]){
	  if(!inrun){
	    runstart=pos;
	    inrun=true;
	  }
	}else{
	  if(inrun){
	    CEBUG("reprun " << actread.getName() << '\t' << runstart << '\t' << pos-1 << endl);
	    tmpmt.from=runstart;
	    tmpmt.to=pos-1;
	    actread.addTagO(tmpmt);
	    inrun=false;
	  }
	}
      }
      if(inrun){
	CEBUG("reprun " << actread.getName() << '\t' << runstart << '\t' << pos-1 << endl);
	tmpmt.from=runstart;
	tmpmt.to=pos-1;
	actread.addTagO(tmpmt);
      }
    }

    Read::setCoutType(Read::AS_TEXT);
    CEBUG("### After ...\n" << actread << endl);


  }

  FUNCEND();
}
//#define CEBUG(bla)




/*************************************************************************
 *
 * comparator for lower_bound below
 *
 *************************************************************************/

inline bool HashStatistics__compareHashStatHashElem_(const hashstat_t & a,
					   const hashstat_t & b);
inline bool HashStatistics__compareHashStatHashElem_(const hashstat_t & a,
					   const hashstat_t & b)
{
  return a.vhash < b.vhash;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)

//uint32 HashStatistics::checkBaitHit(Read & actread)
uint32 HashStatistics::checkBaitHit(Read & actread, std::vector<vhrap_t> & baiting_singlereadvhraparray, std::vector<uint8> & baiting_tagmaskvector)
{
  //, const uint8 basesperhash, vector<hashstat_t> & hashstats, vector<vector<hashstat_t>::const_iterator > & hsshortcuts_begin, vector<vector<hashstat_t>::const_iterator > & hsshortcuts_end)
  FUNCSTART("checkBaitHit()");

  if(!actread.hasValidData()) return 0;
  uint32 slen=actread.getLenClippedSeq();
  if(slen<HS_hs_basesperhash) return 0;

  // don't really need to clear out these re-used vectors
  //   - tagmask just needs to be empty
  //   - singlereadvhraparray needs to be big enough to be
  //     written into by transformSeqToVariableHash(), will be
  //     resized later on num hashes made
  baiting_tagmaskvector.clear();
  if(baiting_singlereadvhraparray.size() < slen){
    baiting_singlereadvhraparray.resize(slen);
  }

  vector<vhrap_t>::iterator srvaI=baiting_singlereadvhraparray.begin();

  vector<Read::bposhashstat_t> & bposhashstats=const_cast<vector<Read::bposhashstat_t> &>(actread.getBPosHashStats());
  uint32 hashesmade;

  {
    int32 bfpos=0;
    int32 bfposinc=1;

    uint32 actreadid=0;

    hashesmade=Skim::transformSeqToVariableHash(
      actreadid,
      actread,
      actread.getClippedSeqAsChar(),
      slen,
      HS_hs_basesperhash,
      srvaI,
      false,
      1,
      baiting_tagmaskvector,
      bposhashstats,
      bfpos,
      bfposinc
      );
  }
  baiting_singlereadvhraparray.resize(hashesmade);

  CEBUG("hashesmade: " << hashesmade << endl);

  vector<hashstat_t>::const_iterator lowerbound;

  vector<hashstat_t>::const_iterator hssearchI;
  srvaI=baiting_singlereadvhraparray.begin();

  hashstat_t hstmp;
  bool foundit;
  uint32 numhits=0;
  for(; srvaI != baiting_singlereadvhraparray.end(); srvaI++){
    CEBUG(*srvaI << '\n');

    lowerbound=HS_hs_hsshortcuts[srvaI->vhash & SKIM3_MAXVHASHMASK].b;

    foundit=false;

    // "SKIM3_empty_vector_hashstat_t.end()" is the "nullptr" replacement
    if(HS_hs_hashstats.end() != lowerbound){
      if(HS_hs_basesperhash>12){
	// with more than 12 bases in a hash, the array is subdivided
	hstmp.vhash=srvaI->vhash;
	hssearchI=lower_bound(lowerbound,
			      HS_hs_hsshortcuts[srvaI->vhash & SKIM3_MAXVHASHMASK].e,
			      hstmp,
			      HashStatistics__compareHashStatHashElem_);
	if(hssearchI != HS_hs_hashstats.end()
	   && hssearchI->vhash == srvaI->vhash) foundit=true;
      }else{
	hssearchI=lowerbound;
	foundit=true;
      }
    }else{
      CEBUG("---------- NO LB HIT??? -------\n");
    }

    if(foundit) {
      ++numhits;
    }
  }

  //cout << "\nskim Needs redo!\n";
  //exit(0);

  FUNCEND();
  return numhits;
}


const hashstat_t * HashStatistics::findVHash(const hashstat_t & searchval)
{
  FUNCSTART("const hashstat_t * HashStatistics::findVHash(const hashstat_t & searchval)");

  const hashstat_t * ret=nullptr;

  // even if executed a couple of million times, this if takes virtually no time at all
  // so keep it
  BUGIFTHROW(unlikely(HS_hs_hsshortcuts.empty()),"no shortcuts made, not ready for searching?");

  if(likely(!HS_hs_hashstats.empty())
     && HS_hs_hashstats.end() != HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].b){
    auto hsI=HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].b;
    // TODO: test with large & diverse data set effect of prefetch
    prefetchrl(&(*hsI));
    if(HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e-hsI > 1){
      // with more than 12 bases in a hash, the array is subdivided
      // TODO: test with large & diverse data set whether this split in lower_bound
      //  vs. simple while loop is OK
      if(HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e-hsI > 4){
	hsI=lower_bound(hsI,
			HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e, // upperbound
			searchval,
			sortHashStatComparator);
      }else{
	while(hsI!=HS_hs_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e && hsI->vhash!=searchval.vhash){
	  ++hsI;
	}
      }
    }
    if(hsI != HS_hs_hashstats.end()
       && hsI->vhash == searchval.vhash) ret=&(*hsI);
  }

  return ret;
}




/*************************************************************************
 *
 * test
 *
 * implicit return:
 *  - dn_vhashindexes with indexes to all valid vhashes in sequence
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
//#define CEBUG(bla)   {if(docebug) {cout << bla; cout.flush();}}
bool HashStatistics::priv_dn_TestSingleSeq(Read & actread, vector<uint8> & dn_allow, vector<vhash_t> & dn_vhashindexes)
{
  FUNCSTART("bool HashStatistics::priv_dn_TestSingleSeq(Read & actread, vector<uint8> & dn_allow, vector<vhash_t> & dn_vhashindexes)");

  //bool docebug=false;

  const uint8 * seq = reinterpret_cast<const uint8 *>(actread.getClippedSeqAsChar());
  uint64 slen=actread.getLenClippedSeq();

  if(slen<HS_hs_basesperhash) return false;

  const char *  namestr=actread.getName().c_str();

  dn_vhashindexes.clear();
  dn_allow.clear();
  dn_allow.resize(slen,1);

  // TODO: option to have it only on MNRr or given HAFx stretches
  auto bhsI=actread.getBPosHashStats().begin();
  bhsI+=actread.getLeftClipoff();
  for(auto ri=0; ri<slen; ++ri){
    if(bhsI->fwd.getFrequency()<2 || !bhsI->fwd.hasConfirmedFwdRev()){
      dn_allow[ri]=0;
    }
  }


  hashstat_t searchval;
  bool takeread=false;

  auto basesperhash=HS_hs_basesperhash;

  SEQTOHASH_LOOPSTART(vhash_t){

    if(HS_hs_hashstats.end() != HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].b){
      auto hsI=HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].b;
      // TODO: test with large & diverse data set effect of prefetch
      prefetchrl(&(*hsI));
      if(HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e-hsI > 1){
	// with more than 12 bases in a hash, the array is subdivided
	// TODO: test with large & diverse data set whether this split in lower_bound
	//  vs. simple while loop is OK
	if(HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e-hsI > 4){
	  searchval.vhash=acthash;
	  hsI=lower_bound(hsI,
			  HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e, // upperbound
			  searchval,
			  sortHashStatComparator);
	}else{
	  while(hsI!=HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e && hsI->vhash!=acthash){
	    ++hsI;
	  }
	}
      }

      if(hsI != HS_hs_hashstats.end()
	 && hsI->vhash == acthash) {
	// hsI on valid valid hash
	auto hsindex=hsI-HS_hs_hashstats.begin();
	CEBUG("hashfound " << seqi << "\t" << hsindex << endl);
	dn_vhashindexes.push_back(hsindex);
	if(dn_allow[seqi] && HS_diginorm_count[hsindex]<10){
	  takeread=true;
	}
      }else{
	CEBUG("no hash? " << seqi << endl);
      }
    }

  }SEQTOHASH_LOOPEND;

  return takeread;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool HashStatistics::digiNormTestRead(Read & actread, bool forcetake)
{
  FUNCSTART("bool HashStatistics::digiNorm(Read & actread)");

  if(unlikely(HS_diginorm_count.empty())){
    HS_diginorm_count.resize(HS_hs_hashstats.size(),0);
  }

  if(!actread.hasTag(Read::REA_defaulttag_MNRr.identifier)) return true;

  bool takeread=priv_dn_TestSingleSeq(actread,HS_diginorm_allow_s1,HS_diginorm_vhashindexes_s1);

  if(forcetake) takeread=true;

  if(takeread){
    CEBUG("dntr take " << actread.getName() << ": " << HS_diginorm_vhashindexes_s1.size() << endl);
    for(auto hsi : HS_diginorm_vhashindexes_s1){
      ++HS_diginorm_count[hsi];
    }
  }

  return takeread;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * test
 *
 *
 *************************************************************************/


//#define CEBUG(bla)   {cout << bla; cout.flush();}

uint32 HashStatistics::estimDigiNormCov(Read & actread)
{
  FUNCSTART("void HashStatistics::estimDigiNormCov(Read & actread)");

  const uint8 * seq = reinterpret_cast<const uint8 *>(actread.getClippedSeqAsChar());
  uint64 slen=actread.getLenClippedSeq();

  if(slen<HS_hs_basesperhash) return 1;

  const char *  namestr=actread.getName().c_str();

  auto basesperhash=HS_hs_basesperhash;


  double dncmin=10000000.0;
  bool hasnewmin=false;
  //double dncmax=0.0;
  double totaladd=0.0;
  uint32 numtotal=0;


  hashstat_t searchval;

  SEQTOHASH_LOOPSTART(vhash_t){

    if(HS_hs_hashstats.end() != HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].b){
      auto hsI=HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].b;
      // TODO: test with large & diverse data set effect of prefetch
      prefetchrl(&(*hsI));
      if(HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e-hsI > 1){
	// with more than 12 bases in a hash, the array is subdivided
	// TODO: test with large & diverse data set whether this split in lower_bound
	//  vs. simple while loop is OK
	if(HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e-hsI > 4){
	  searchval.vhash=acthash;
	  hsI=lower_bound(hsI,
			  HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e, // upperbound
			  searchval,
			  sortHashStatComparator);
	}else{
	  while(hsI!=HS_hs_hsshortcuts[acthash & HS_MAXVHASHMASK].e && hsI->vhash!=acthash){
	    ++hsI;
	  }
	}
      }

      if(hsI != HS_hs_hashstats.end()
	 && hsI->vhash == acthash) {
	// hsI on valid valid hash

	auto hsindex=hsI-HS_hs_hashstats.begin();

	// *sigh* this is very well possible
	// BUGIFTHROW(HS_diginorm_count[hsindex]==0,actread.getName() << " HS_diginorm_count[hsindex]==0 ???");
	// e.g. in first pass a DGNR tag gets set like so: ...TTTTT
	// and somewhere the read gets cut back like so: ...
	// then the tag still may exist.
	// therefore: if count is 0, then do as if no hash existed

	if(HS_diginorm_count[hsindex]!=0){
	  double actdnc=static_cast<double>(hsI->count)/HS_diginorm_count[hsindex];
	  ++numtotal;
	  totaladd+=actdnc;
	  //if(actdnc>dncmax) dncmax=actdnc;
	  if(actdnc<dncmin) {
	    dncmin=actdnc;
	    hasnewmin=true;
	  }
	  CEBUG(actread.getName() << "\t" << hsI->count << "\t" << HS_diginorm_count[hsindex] << "\t" << actdnc << endl);
	}
      }else{
	//BUGIFTHROW(true,"Can't be?");
	// well ... can be, as the haststats does contain only "valid" hashes where we are sure they're present multiple times etc.,
	//  we might encounter a hash here which is not in hashstats.
	// therefore, this is a singlet event

	// but singlet events were, by default, not accounted for in digiNormTestRead(), therefore
	//  we should not return "1" but simply do nothing and continue calculating
	//return 1;
      }
    }

  }SEQTOHASH_LOOPEND;

  CEBUG("dncstats " << actread.getName() << ": " << dncmin);
  //CEBUG("\t" << dncmax(;
  CEBUG("\t" << totaladd/numtotal);
  CEBUG(endl);

  //if(dncmin<=6) return dncmin;

  if(!hasnewmin) {
    CEBUG("WTH??? " << actread.getName() << " has no dncmin ???\n");
    return 1;
  }

  return static_cast<uint32>((totaladd/numtotal)+.5);

  //return static_cast<uint32>(dncmin+0.5);
}
//#define CEBUG(bla)



/***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 *
 * New hashstats, bloomfilter + streaming
 *
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************
 ***********************************************************************************************************/


uint32 NHashStatistics::HSN_hs_magic=0x4D4C6873;  // magic: "MLhs" MiraLibHashStat




/*************************************************************************
 *
 *
 *
 *************************************************************************/

NHashStatistics::~NHashStatistics()
{
  deleteBloomFilter();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::deleteBloomFilter()
{
  if(HSN_bloomfilter!=nullptr){
    delete HSN_bloomfilter;
    HSN_bloomfilter=nullptr;
  }
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::analyseReadPool(ReadPool & rp)
{
  FUNCSTART("void NHashStatistics::analsyeReadPool(ReadPool & rp)");

  ProgressIndicator<int64> pi(0,rp.size());

  for(uint32 ri=0; ri<rp.size(); ++ri){
    rp[ri].getClippedSeqAsChar();
    rp[ri].getClippedComplementSeqAsChar();
  }

  dateStamp(cout);

  for(uint32 step=1; step<=HSN_needsteps; ++step){
    pi.reset(0,rp.size());
    for(uint32 ri=0; ri<rp.size(); ++ri){
      pi.progress(ri);
      learnSequence(rp[ri].getClippedSeqAsChar(),
		    rp[ri].getLenClippedSeq(),
		    rp[ri].getName().c_str(),
		    0,
		    false);
      learnSequence(rp[ri].getClippedComplementSeqAsChar(),
		    rp[ri].getLenClippedSeq(),
		    rp[ri].getName().c_str(),
		    0,
		    true);
    }
    pi.finishAtOnce();
    cout << endl;
    finaliseStep();
    dateStamp(cout);
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::setupNewAnalysis(const uint8 bfbits, const uint32 bfnumkeys, const uint8 basesperhash, uint16 numsteps)
{
  FUNCSTART("void NHashStatistics::setupNewAnalysis(const uint8 bfbits, const uint32 bfnumkeys,  const uint8 basesperhash, uint16 numsteps");
  BUGIFTHROW(HSN_bloomfilter!=nullptr,"HSN_bloomfilter!=nullptr ??");

  HSN_bloomfilter=new BloomFilter(bfbits,bfnumkeys);
  HSN_basesperhash=basesperhash;

  BUGIFTHROW(numsteps==0,"numsteps==0 ???");
  BUGIFTHROW(numsteps>3,"numsteps " << numsteps << " not in 1,2,3.");

  HSN_needsteps=numsteps;

  if(numsteps==1){
    HSN_step=1001;
    cout << "Counting hashes (quick, slightly inaccurate, 1 pass): step 1" << endl;
  }else if(numsteps==2){
    HSN_step=2001;
    cout << "Counting hashes (quick, accurate, 2 pass): step 1" << endl;
  }else{
    HSN_step=3001;
    cout << "Counting hashes (accurate, savemem, 3 pass): step 1" << endl;
  }

  HSN_hs_sortstatus=0;
  HSN_hs_needsconsolidation=false;

}

void NHashStatistics::learnSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics::learnSequenceStep(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  BUGIFTHROW(HSN_bloomfilter==nullptr,"HSN_bloomfilter==nullptr ???");

  switch(HSN_step){
  case 1001 : {
    learnSequenceQuick1(seqvoid,slen,namestr,seqtype,isreverse,false);
    break;
  }
  case 2001 : {
    learnSequenceStep1(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 2002 : {
    learnSequenceQuick1(seqvoid,slen,namestr,seqtype,isreverse,true);
    break;
  }
  case 3001 : {
    learnSequenceStep1(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 3002 : {
    learnSequenceStep2(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 3003 : {
    learnSequenceStep3(seqvoid,slen,namestr,seqtype,isreverse);
    break;
  }
  case 32678 : {
    BUGIFTHROW(true,"HSN_step 32678, nothing more to learn???");
    break;
  }
  default :{
    BUGIFTHROW(true,"unknown HSN_step "  << static_cast<int16>(HSN_step));
  }
  }
}


void NHashStatistics::finaliseStep()
{
  FUNCSTART("void NHashStatistics::finaliseStep1()");

  if(HSN_step==1){
    cout << "Counting hashes: finalised step 1, switching to step 2" << endl;
    HSN_hsv_hashstats.reserve(HSN_bloomfilter->getNumKMersSeenGE2());
    HSN_bloomfilter->reset();
  }else if(HSN_step==2){
    cout << "Counting hashes: finalising step 2 ..."; cout.flush();
    makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
    cout << " done.\nCounting hashes: step 3" << endl;
  }else if(HSN_step==3){
    cout << "Trimming out erroneous single hashes ..."; cout.flush();
    auto oldsize=HSN_hsv_hashstats.size();
    trimHashStatsByFrequency(-1,-1,2);
    cout << " done. Trimmed " << oldsize-HSN_hsv_hashstats.size() << " hashes, " << HSN_hsv_hashstats.size() << " remaining" << endl;
    HSN_step=32767;
  }else if(HSN_step==1001){
    cout << "quick done\n";
    HSN_step=32767;
  }else if(HSN_step==2001){
    cout << "quick 2.1 done\n";
  }else if(HSN_step==2002){
    cout << "quick 2.2 done\n";
    HSN_step=32767;
  }else{
    BUGIFTHROW(true,"HSN_step is " << HSN_step << " ???");
  }

  ++HSN_step;
}


void NHashStatistics::learnSequenceQuick1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse, bool lookuponly)
{
  FUNCSTART("void NHashStatistics::learnSequenceQuick(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  static hscounts_t tmphs;

  BUGIFTHROW(HSN_step!=1001 && HSN_step!=2002,"HSN_step!=1001 && HSN_step!=2002 ???");

  uint32 countinit=1;
  if(lookuponly) ++countinit;

  auto basesperhash=HSN_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);
  int bfres=0;
  SEQTOHASH_LOOPSTART(vhash_t){
    auto thispos=seqi-basesperhash;
    if(unlikely(isreverse)){
      thispos=slen-1-seqi;
    }
    if(unlikely(thispos>1020)) thispos=1020;
    thispos/=4;
    auto umhsI=HSN_hsum_hashstats.find(acthash);
    if(umhsI!=HSN_hsum_hashstats.end()){
      if(thispos < umhsI->second.lowposd4){
	umhsI->second.lowposd4=static_cast<uint8>(thispos);
      }
      if(umhsI->second.seqtype!=seqtype){
	seqtype=0xf;
      }
      if(unlikely(isreverse)){
	if(unlikely(++(umhsI->second.rcount)==0)) --(umhsI->second.rcount);
      }else{
	if(unlikely(++(umhsI->second.fcount)==0)) --(umhsI->second.fcount);
      }
    }else{
      if(lookuponly){
	bfres=HSN_bloomfilter->isNonUnique(acthash);
      }else{
	bfres=HSN_bloomfilter->addVHash(acthash);
      }
      if(bfres==1){
	// make new in unordered map!
	tmphs.lowposd4=thispos;
	tmphs.seqtype=seqtype;
	if(isreverse){
	  tmphs.fcount=0;
	  tmphs.rcount=countinit;
	}else{
	  tmphs.fcount=countinit;
	  tmphs.rcount=0;
	}
	HSN_hsum_hashstats[acthash]=tmphs;
      }
      // hmmm ... should not happen here
      BUGIFTHROW(bfres==2,"bfres==2 ???");
    }
  }SEQTOHASH_LOOPEND;
}


void NHashStatistics::learnSequenceQuick2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
}


void NHashStatistics::learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics::learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  BUGIFTHROW(HSN_step!=1 && HSN_step!=2001,"HSN_step!=1 && HSN_step!=2001 ???");
//  HSN_bloomfilter->addSequenceToBloomfield(seqvoid, slen, HSN_basesperhash, namestr);

  auto basesperhash=HSN_basesperhash;
  for(uint32 xxi=0; xxi<2;++xxi){
    const uint8 * seq=static_cast<const uint8 *>(seqvoid);
    SEQTOHASH_LOOPSTART(vhash_t);
    if(xxi){
      (void) HSN_bloomfilter->addVHash(acthash);
    }else{
      HSN_bloomfilter->prefetchVHash(acthash);
    }
    SEQTOHASH_LOOPEND;
  }
}

void NHashStatistics::learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics::learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)");

  static nhashstat_t tmphs;

  BUGIFTHROW(HSN_step!=2,"HSN_step!=2 ???");

  auto basesperhash=HSN_basesperhash;
  for(uint32 xxi=0; xxi<2;++xxi){
    const uint8 * seq=static_cast<const uint8 *>(seqvoid);
    SEQTOHASH_LOOPSTART(vhash_t);
    if(xxi){
      if(HSN_bloomfilter->addVHash(acthash)==1){
	BUGIFTHROW(HSN_hsv_hashstats.size()==HSN_hsv_hashstats.capacity(),"HSN_hsv_hashstats.size()==hstable.capacity() ???");
	tmphs.vhash=acthash;
	HSN_hsv_hashstats.push_back(tmphs);
      }
    }else{
      HSN_bloomfilter->prefetchVHash(acthash);
    }
    SEQTOHASH_LOOPEND;
  }
}


void NHashStatistics::learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse)
{
  FUNCSTART("void NHashStatistics::learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse))");

  BUGIFTHROW(HSN_step!=3,"HSN_step!=3 ???");

  nhashstat_t tmphs;

  // we expect to have much more kmers occurring more than once than single kmers
  // therefore, do not use a query to the bloomfilter before using find() on the
  //  vector<hashstat_t>, it's just a waste of time

  auto basesperhash=HSN_basesperhash;
  const uint8 * seq=static_cast<const uint8 *>(seqvoid);
  SEQTOHASH_LOOPSTART(vhash_t){
    tmphs.vhash=acthash;
    auto hptr=const_cast<nhashstat_t *>(findVHash(tmphs));
    if(likely(hptr!=nullptr)){
      auto thispos=seqi-basesperhash;
      if(unlikely(isreverse)){
	thispos=slen-1-seqi;
      }
      if(unlikely(thispos>1020)) thispos=1020;
      thispos/=4;
      if(unlikely(hptr->hsc.fcount==0 && hptr->hsc.rcount==0)){
	hptr->hsc.seqtype=seqtype;
	hptr->hsc.lowposd4=static_cast<uint8>(thispos);
      }else if(thispos < hptr->hsc.lowposd4){
	hptr->hsc.lowposd4=static_cast<uint8>(thispos);
      }

      if(unlikely(isreverse)){
	if(unlikely(++(hptr->hsc.rcount)==0)) --(hptr->hsc.rcount);
      }else{
	if(unlikely(++(hptr->hsc.fcount)==0)) --(hptr->hsc.fcount);
      }
    }
  }SEQTOHASH_LOOPEND;

}

#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)

const NHashStatistics::nhashstat_t * NHashStatistics::findVHash(const nhashstat_t & searchval)
{
  const nhashstat_t * ret=nullptr;

  // even if executed a couple of million times, this if takes virtually no time at all
  // so keep it
  if(unlikely(HSN_hsv_hsshortcuts.empty())) {
    makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
  }

  if(likely(!HSN_hsv_hashstats.empty())
     && HSN_hsv_hashstats.end() != HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].b){
    auto hsI=HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].b;
    // TODO: test with large & diverse data set effect of prefetch
    prefetchrl(&(*hsI));
    if(HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e-hsI > 1){
      // with more than 12 bases in a hash, the array is subdivided
      // TODO: test with large & diverse data set whether this split in lower_bound
      //  vs. simple while loop is OK
      if(HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e-hsI > 4){
	hsI=lower_bound(hsI,
			HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e, // upperbound
			searchval,
			sortHashStatComparator);
      }else{
	while(hsI!=HSN_hsv_hsshortcuts[searchval.vhash & HS_MAXVHASHMASK].e && hsI->vhash!=searchval.vhash){
	  ++hsI;
	}
      }
    }
    if(hsI != HSN_hsv_hashstats.end()
       && hsI->vhash == searchval.vhash) ret=&(*hsI);
  }

  return ret;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics::trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  if(HSN_hsum_hashstats.empty()){
    trimHashVStatsByFrequency(minfwd,minrev,mintotal);
  }else{
    trimHashMStatsByFrequency(minfwd,minrev,mintotal);
  }

}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics::trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HSN_hsv_hashstats.begin();
  auto dstI=srcI;

  for(; srcI!=HSN_hsv_hashstats.end(); ++srcI){
    bool ok=true;
    if((minfwd>=0 && srcI->hsc.fcount<minfwd)
       || (minrev>=0 && srcI->hsc.rcount<minrev)
       || (mintotal>=0 && srcI->hsc.fcount+srcI->hsc.rcount < mintotal)){
      ok=false;
      CEBUG("rm\t");
    }else{
      CEBUG("keep\t");
    }
    CEBUG(srcI-HSN_hsv_hashstats.begin() << "\t" << *srcI << endl);
    *dstI=*srcI;
    if(ok)++dstI;
  }
  HSN_hsv_hashstats.resize(dstI-HSN_hsv_hashstats.begin());
  HSN_hsv_hsshortcuts.clear();
  HSN_hs_dist.clear();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void NHashStatistics::trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)
{
  FUNCSTART("void NHashStatistics::trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal)");

  auto srcI=HSN_hsum_hashstats.begin();

  for(; srcI!=HSN_hsum_hashstats.end();){
    bool ok=true;
    if((minfwd>=0 && srcI->second.fcount<minfwd)
       || (minrev>=0 && srcI->second.rcount<minrev)
       || (mintotal>=0 && srcI->second.fcount+srcI->second.rcount < mintotal)){
      ok=false;
      CEBUG("rm\t");
    }else{
      CEBUG("keep\t");
    }
    CEBUG(srcI-HSN_hsum_hashstats.begin() << "\t" << *srcI << endl);
    if(ok) {
      ++srcI;
    }else{
      srcI=HSN_hsum_hashstats.erase(srcI);
    }
  }
}


/*************************************************************************
 *
 * needs:
 *  - hashstats filled with entries (can be unsorted, will be re-sorted
 *    anyway)
 *
 * returns:
 *  - hashstats array sorted by low 24 bit (low to high), then by vhash
 *  - elements .b and .e in hsshortcuts pointing to start and end of each
 *    low 24 bit group of same value
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

void NHashStatistics::makeNHashStatArrayShortcuts(vector<nhashstat_t> & hashstats, const uint8 basesperhash, vector<hsvbendit_t> & hsshortcuts)
{
  FUNCSTART("void HashStatistics::makeNHashStatArrayShortcuts(vector<hashstat_t> & hashstats, const uint8 basesperhash, vector<hsvbendit_t> & hsshortcuts)");

  CEBUG("makeNHashStatArrayShortcuts: basesperhash: " << static_cast<uint16>(basesperhash) << "\n");

  BUGIFTHROW(basesperhash==0, "basesperhash == 0 ???");

  for(size_t hsi=0; hsi<hashstats.size(); ++hsi){
    CEBUG(hashstats[hsi] << '\n');
  }

  sortLow24Bit(hashstats,HSN_hs_sortstatus);

  hsshortcuts.clear();
  auto hsI=hashstats.cbegin();

  {
    hsvbendit_t tmpb;
    tmpb.b=hashstats.end();
    tmpb.e=hashstats.end();
    hsshortcuts.resize(
      1<<(min(static_cast<uint8>(12),basesperhash)*2),
      tmpb
      );
  }

  if(hsI==hashstats.end()) return;

  CEBUG("hsshortcuts.size(): " << hsshortcuts.size() << endl);

  vhash_t acthash= (hsI->vhash & HS_MAXVHASHMASK);
  while(hsI != hashstats.end()){
    CEBUG("begin " << hex << acthash << dec << " is: " << *hsI << endl);
    hsshortcuts[acthash].b=hsI;
    for(;(hsI != hashstats.end()) && ((hsI->vhash & HS_MAXVHASHMASK) == acthash); hsI++) {
      CEBUG("INC\n")
    }
    CEBUG("end " << hex << acthash << dec << " is: " << *hsI << endl);
    hsshortcuts[acthash].e=hsI;
    //cout << "vhash: " << hex << acthash << "\t" << dec << hsshortcuts_end[acthash]-hsshortcuts_begin[acthash] << '\n';
    if(hsI != hashstats.end()) acthash= hsI->vhash & HS_MAXVHASHMASK;
  }

  FUNCEND();
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::dumpHealth(ostream & fout)
{
  FUNCSTART("void NHashStatistics::dumpHealth(ostream & fout)");
  BUGIFTHROW(HSN_bloomfilter==nullptr,"HSN_bloomfilter==nullptr ???");
  fout << *HSN_bloomfilter;
}



/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::moveStatCountMapToVector()
{
  if(HSN_hsum_hashstats.empty()) return;

  nhashstat_t tmphs;

  HSN_hsv_hashstats.clear();
  HSN_hsv_hashstats.reserve(HSN_hsum_hashstats.size());
  for(auto & hsme : HSN_hsum_hashstats){
    tmphs.vhash=hsme.first;
    tmphs.hsc=hsme.second;
    HSN_hsv_hashstats.push_back(tmphs);
  }
  HSN_hsum_hashstats.clear();
  makeNHashStatArrayShortcuts(HSN_hsv_hashstats, HSN_basesperhash, HSN_hsv_hsshortcuts);
}

/*************************************************************************
 *
 *
 *************************************************************************/

uint64 NHashStatistics::calcHashDistrib(vector<uint64> & hsdist)
{
  FUNCSTART("void NHashStatistics::calcHashDistrib(vector<uint64> & hsdist)");

  hsdist.clear();

  if(HSN_hsv_hashstats.empty()) moveStatCountMapToVector();

  uint64 maxc=0;
  uint64 totalc=0;
  for(auto & hsve : HSN_hsv_hashstats){
    totalc+=static_cast<uint64>(hsve.hsc.fcount+hsve.hsc.rcount);
    maxc=max(maxc,static_cast<uint64>(hsve.hsc.fcount+hsve.hsc.rcount));
  }
  hsdist.resize(maxc+1,0);
  for(auto & hse : HSN_hsv_hashstats){
    ++hsdist[hse.hsc.fcount+hse.hsc.rcount];
  }

  return totalc;
}


/*************************************************************************
 *
 *
 *************************************************************************/

vector<uint64> & NHashStatistics::getHashDistrib()
{
  if(HSN_hs_dist.empty()) calcHashDistrib(HSN_hs_dist);
  return HSN_hs_dist;
}


/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::dumpHashDistrib(ostream & ostr)
{
  FUNCSTART("void NHashStatistics::dumpHashDistrib(vector<uint64> & hsdist, ostream & ostr) const");

  auto & hsdist=getHashDistrib();

  uint64 totalc=0;
  uint64 hsdi=0;
  for(auto hsde : hsdist){
    totalc+=hsde*hsdi;
    ++hsdi;
  }

  double dtotalc=static_cast<double>(totalc);
  uint64 cumc=0;
  hsdi=0;
  for(auto hsde : hsdist){
    cumc+=hsde*hsdi;
    double frac=static_cast<double>(cumc)/dtotalc;
    ostr << hsdi
	 << '\t' << hsde
	 << '\t' << hsde*hsdi
	 << '\t' << cumc
	 << '\t' << frac
	 << '\n';
    ++hsdi;
  }
}


/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::saveHashStatistics(const string & filename, bool deleteoldfile)
{
  FUNCSTART("void NHashStatistics::saveHashStatistics(const string & filename, bool deleteoldfile)");

  ofstream fout;
  openFileForAppend(filename,fout,deleteoldfile);
  try{
    if(!fout){
      MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is the disk full? Are permissions set right?");
    }
    saveHashStatistics(fout);
  }
  catch(Notify n){
    cout << "Error for file " << filename << endl;
    n.handleError(THISFUNC);
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::loadHashStatistics(const string & filename)
{
  FUNCSTART("void NHashStatistics::loadHashStatistics(const string & filename)");

  ifstream fin;
  try{
    fin.open(filename.c_str(),ios::in);
    if(!fin){
      MIRANOTIFY(Notify::FATAL,"Could not open " << filename << ", is it present? Are permissions set right?");
    }
    loadHashStatistics(fin);
  }
  catch(Notify n){
    cout << "Error while loading file " << filename << endl;
    n.handleError(THISFUNC);
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::saveHashStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics::saveHashStatistics(ostream & ostr)");

  if(HSN_hsum_hashstats.empty()){
    saveHashVStatistics(ostr);
  }else{
    saveHashMStatistics(ostr);
  }

}

/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::saveHashVStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics::saveHashVStatistics(ostream & ostr)");

  ostr.write(reinterpret_cast<const char *>(&HSN_hs_magic),4);
  ostr.put(2); // version
  ostr.put(static_cast<uint8>(HSN_basesperhash));
  ostr.put(HSN_hs_sortstatus);
  ostr.put(0); // padding
  uint64 hssize=HSN_hsv_hashstats.size();
  ostr.write(reinterpret_cast<const char *>(&hssize),sizeof(hssize));
  if(!HSN_hsv_hashstats.empty()){
    ostr.write(reinterpret_cast<const char *>(&HSN_hsv_hashstats[0]),
	       sizeof(nhashstat_t)*hssize);
  }
  if(ostr.bad()){
    MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics (1). Disk full? Changed permissions?");
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::saveHashMStatistics(ostream & ostr)
{
  FUNCSTART("void NHashStatistics::saveHashMStatistics(ostream & ostr)");

  nhashstat_t tmphs;
  ostr.write(reinterpret_cast<const char *>(&HSN_hs_magic),4);
  ostr.put(2); // version
  ostr.put(static_cast<uint8>(HSN_basesperhash));
  ostr.put(0); // sorting: not sorted
  ostr.put(0); // padding
  uint64 hssize=HSN_hsum_hashstats.size();
  ostr.write(reinterpret_cast<const char *>(&hssize),sizeof(hssize));
  for(auto & hsume : HSN_hsum_hashstats){
    tmphs.vhash=hsume.first;
    tmphs.hsc.lowposd4=hsume.second.lowposd4;
    tmphs.hsc.fcount=hsume.second.fcount;
    tmphs.hsc.seqtype=hsume.second.seqtype;
    tmphs.hsc.rcount=hsume.second.rcount;
    ostr.write(reinterpret_cast<const char *>(&tmphs), sizeof(nhashstat_t));
  }
  if(ostr.bad()){
    MIRANOTIFY(Notify::FATAL, "Could not save anymore the hash statistics. Disk full? Changed permissions?");
  }
}

/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::loadHashStatistics(istream & istr)
{
  FUNCSTART("bool NHashStatistics::loadHashStatistics(istream & istr)");

  auto localmagic=HSN_hs_magic;
  istr.read(reinterpret_cast<char *>(&localmagic),4);
  if(istr.gcount()!=4
     || localmagic!=HSN_hs_magic) {
    MIRANOTIFY(Notify::FATAL,"No magic found?\n");
  }
  uint8 tmpbyte=0;
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(tmpbyte!=2) {
    MIRANOTIFY(Notify::FATAL,"Not version 2?\n");
  }
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(tmpbyte==0 || tmpbyte>32){
    MIRANOTIFY(Notify::FATAL,"Invalid kmer size " << static_cast<uint16>(tmpbyte) << " ???");
  }
  if(HSN_basesperhash!=0 && !HSN_hsv_hashstats.empty() && tmpbyte!= HSN_basesperhash){
    MIRANOTIFY(Notify::FATAL,"Current hashstat kmer size is " << static_cast<uint16>(HSN_basesperhash)
	       << ", but kmer size in data to load is " << static_cast<uint16>(tmpbyte)
	       << " ???\n");
   }else{
    HSN_basesperhash=tmpbyte;
  }
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);
  if(!HSN_hsv_hashstats.empty()){
    HSN_hs_sortstatus=0;
    HSN_hs_needsconsolidation=true;
    MIRANOTIFY(Notify::FATAL,"Appending to existing hashstat not implemented yet\n");
  }else{
    HSN_hs_sortstatus=tmpbyte;
  }
  // padd byte
  istr.read(reinterpret_cast<char *>(&tmpbyte),1);

  uint64 numelem=0;
  istr.read(reinterpret_cast<char *>(&numelem),8);
  if(numelem){
    auto oldsize=HSN_hsv_hashstats.size();
    HSN_hsv_hashstats.resize(HSN_hsv_hashstats.size()+numelem);
    HSN_hsv_hsshortcuts.clear();
    istr.read(reinterpret_cast<char *>(&HSN_hsv_hashstats[oldsize]),numelem*sizeof(nhashstat_t));
    if(istr.gcount()!=numelem*sizeof(nhashstat_t)){
      MIRANOTIFY(Notify::FATAL,"Expected to read " << numelem*sizeof(nhashstat_t) << " bytes, but got " << istr.gcount() << endl);
    }
  }
}


/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::hash2string(vhash_t hash, uint8 basesperhash, std::string & str)
{
  static char acgtc[4]={'A','C','G','T'};

  str.clear();
  str.resize(basesperhash,' ');
  auto srI=str.rbegin();
  for(auto ci=0; ci<basesperhash; ++ci, ++srI){
    *srI=acgtc[hash&3];
    hash>>=2;
  }
}


/*************************************************************************
 *
 *
 *************************************************************************/

void NHashStatistics::dumpHashCount(ostream & ostr)
{
  FUNCSTART("void NHashStatistics::dumpHashCount(ostream & ostr)");

  string tmpstr;
  for(auto & hse : HSN_hsv_hashstats){
    hash2string(hse.vhash,HSN_basesperhash,tmpstr);
    cout << tmpstr
	 << '\t' << hse.hsc.fcount
	 << '\t' << hse.hsc.rcount
	 << '\t' << hse.hsc.fcount + hse.hsc.rcount
	 << '\n';
  }
}
