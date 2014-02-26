/*
 * Written by Bastien Chevreux (BaCh)
 * Copyright (C) 2007 by Bastien Chevreux
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
 */

#ifndef _bas_hashstats_h_
#define _bas_hashstats_h_

#include <unordered_map>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/regex.hpp>

#include "stdinc/defines.H"
#include "stdinc/stlincludes.H"

#include "util/progressindic.H"
#include "mira/readpool.H"
#include "mira/bloomfilter.H"
#include "ads.H"



// now almost same as hashstat? takes up 16 bytes anyway (like diskhash)
struct hashstat_t {
  vhash_t vhash;       // vhash,
  uint32  count;       // how often this hash appeared TODO: somewhat broken for hash statistics with mkfr>1, would need separate rev counter
  uint16  lowpos;      // lowest hash position (lowest pos: 0)
  uint8   seqtype;     // sequencing type

  bool    hasfwd:1;
  bool    hasfwdthresholdok:1;
  bool    hasrev:1;
  bool    hasrevthresholdok:1;

  bool    hasfwdrevthresholdok:1;
  bool    hasmultipleseqtype:1;

  bool    iskmerfork:1;

  //bool operator<(const vhrap_t & other) const {return vhash < other.vhash;};
  friend std::ostream & operator<<(std::ostream &ostr, const hashstat_t & v){
    ostr << "hash: " << std::hex << v.vhash << std::dec
	 << "\tcount: " << v.count
	 << "\tlpos: " << v.lowpos
	 << "\tseqtype: " << static_cast<uint16>(v.seqtype)
	 << "\tfwd: " << v.hasfwd
	 << "\tfwdto: " << v.hasfwdthresholdok
	 << "\trev: " << v.hasrev
	 << "\trevto: " << v.hasrevthresholdok
	 << "\tfrto: " << v.hasfwdrevthresholdok
	 << "\tmseq: " << v.hasmultipleseqtype
	 << "\tkmrf: " << v.iskmerfork;
    return ostr;
  }

  hashstat_t() {
    uint64 * ptr=&vhash;
    *ptr=0;
    ++ptr;
    *ptr=0;
  }
};


struct vhrap_t {
  vhash_t vhash;       // vhash,
  uint32 readid;       // readid
  uint16 hashpos;      // (hash)position  (lowest pos: basesperhash-1)

  Read::bhashstat_t bhashstats; // baseflags for this hash

  //bool operator<(const vhrap_t & other) const {return vhash < other.vhash;};
  friend std::ostream & operator<<(std::ostream &ostr, const vhrap_t & v){
    ostr << "hash: " << std::hex << v.vhash
	 << "\trid: " <<std::dec << v.readid
	 << "\thpos: " << v.hashpos
	 << "\tbhs: " << v.bhashstats;
    return ostr;
  }
};



class HashStatistics
{
private:
  struct hsvbendit_t {
    std::vector<hashstat_t>::const_iterator b;
    std::vector<hashstat_t>::const_iterator e;
  };

  static size_t HS_numelementsperbuffer;

  ReadPool * HS_readpoolptr;

  double HS_freqest_minnormal;
  double HS_freqest_maxnormal;
  double HS_freqest_repeat;
  double HS_freqest_heavyrepeat;
  double HS_freqest_crazyrepeat;
  uint32 HS_nastyrepeatratio;
  uint32 HS_nastyrepeatcoverage;

  uint8 HS_hs_basesperhash;
  std::vector<hashstat_t> HS_hs_hashstats;
  std::vector<hsvbendit_t> HS_hs_hsshortcuts;

  // and the avg frequency
  size_t HS_avg_freq_corrected;
  size_t HS_avg_freq_raw;
  size_t HS_avg_freq_min;

  size_t HS_avg_freq_taken;


  // size of ReadGroupLib::getNumSequencingTypes()
  // kmers having frequency < given value are marked as frequency 1 (or 0 ???) so that
  //  they can be masked later on in data processing
  //
  // assignReadBaseStatistics will initialise it to {0, 0, ...} if not of correct size
  std::vector<uint32> HS_rarekmermasking;


  /*
    Baiting
   */

  std::vector<vhrap_t> HS_baiting_singlereadvhraparray;
  std::vector<uint8>   HS_baiting_tagmaskvector;

  bool HS_logflag_hashcount;

  /*
    Multithreading assign read base statistics
  */

  struct arbs_threadsharecontrol_t {
    boost::mutex accessmutex;
    int32 from;
    int32 to;
    int32 todo;
    int32 done;
    int32 stepping;

    // need to go via this as the boost:bind does not like a "ReadPool &" as parameter
    // and it also cannot have more than 9 parameters in total ... we'd have more with the below
    ReadPool * rpptr;
    size_t avghashcov;

    std::vector<hashstat_t> * hashstatsptr;
    std::vector<hsvbendit_t> * hsscptr;
    std::vector<uint32> * rarekmermaskingptr;

    uint8 basesperhash;
    bool masknastyrepeats;
  };

  //
  // Digital normalisation
  //

  std::vector<size_t>  HS_diginorm_count;

  std::vector<uint8>   HS_diginorm_allow_s1;
  std::vector<uint8>   HS_diginorm_allow_s2;
  std::vector<vhash_t> HS_diginorm_vhashindexes_s1;
  std::vector<vhash_t> HS_diginorm_vhashindexes_s2;

  /*
    Logging, debugging
  */


public:


  //Functions
private:
  void hashes2disk(std::vector<std::string> & hashfilenames,
		   std::vector<size_t> & elementsperfile,
		   ReadPool & rp,
		   bool checkusedinassembly,
		   bool fwdandrev,
		   uint32 fwdrevmin,
		   uint8  basesperhash,
		   uint32 millionhashesperbuffer,
		   bool rarekmerearlykill,
		   const std::string & directory);
  size_t writeCompressedHFB(std::vector<hashstat_t> & hfb,
			    uint32 fwdrevmin,
			    FILE * fileptr,
			    bool force,
			    bool rarekmerearlykill);
  void compressHashStatBufferInPlace(std::vector<hashstat_t> & hsb,
				     uint32 fwdrevmin,
				     bool alsosavesinglehashes);
  size_t createHashStatisticsFile(std::string & hashstatfilename,
				  std::vector<std::string> & hashfilenames,
				  std::vector<size_t> & elementsperfile,
				  uint32 fwdrevmin,
				  bool alsosavesinglehashes,
				  const std::string & directory);

  void calcAvgHashFreq();

  size_t calcMidHashStatIndex(const std::vector<hashstat_t> & hashstats,
			      size_t dontcarepercent);

  void makeHashStatArrayShortcuts();
  inline static bool sortHashStatComparator(const hashstat_t & a, const hashstat_t & b){
    return a.vhash < b.vhash;
  }

  void calcKMerForks(uint32 mincount);
  void ckmf_helper(vhash_t HashStatistics__vhashmask, uint32 mincount);

  void priv_arb_thread(uint32 threadnum, arbs_threadsharecontrol_t * tscptr);
  void priv_arb_DoStuff(
    ReadPool & rp,
    size_t avgcov,
    std::vector<hashstat_t> & hashstats,
    const uint8 basesperhash,
    std::vector<hsvbendit_t> & hsshortcuts,
    bool masknastyrepeats,
    std::vector<uint32> & minkmer,
    int32 fromid,
    int32 toid
    );

  void correctReadBaseStatisticsByRMB(ReadPool & rp,
				      const uint8 basesperhash);

  bool priv_dn_TestSingleSeq(Read & actread, std::vector<uint8> & dn_allow, std::vector<vhash_t> & dn_vhashindexes);


public:
  HashStatistics() {
    setHashFrequencyRatios(.4,1.6,2.0,8.0,20.0,static_cast<uint32>(100),0);
    HS_logflag_hashcount=false;
    HS_avg_freq_corrected=0;
    HS_avg_freq_raw=0;
    HS_avg_freq_min=0;
    HS_avg_freq_taken=0;
    HS_hs_basesperhash=0;
  };
  HashStatistics(HashStatistics const &other);
  ~HashStatistics() {};

  HashStatistics const & operator=(HashStatistics const & other);
  friend std::ostream & operator<<(std::ostream &ostr, HashStatistics const &theskim);

  void discard();

  bool hasStatistics() const { return !HS_hs_hashstats.empty();}

  size_t getAvgHashFreqCorrected() const { return HS_avg_freq_corrected;};
  size_t getAvgHashFreqRaw() const { return HS_avg_freq_raw;};

  void setAvgHashFreqMinimum(size_t m) { HS_avg_freq_min=m;};

  void setHashFrequencyRatios(double freqest_minnormal,
			      double freqest_maxnormal,
			      double freqest_repeat,
			      double freqest_heavyrepeat,
			      double freqest_crazyrepeat,
			      uint32 nastyrepeatratio,
			      uint32 nastyrepeatcoverage);

  void prepareHashStatistics(const std::string & directory,
			     ReadPool & rp, bool checkusedinassembly, bool onlyagainstrails, bool alsosavesinglehashes,
			     bool fwdandrev,
			     uint32 fwdrevmin,
			     uint8  basesperhash,
			     uint32 millionhashesperbuffer,
			     bool rarekmerearlykill,
			     std::string & hashstatfilename);

  void assignReadBaseStatistics_MultiThread(uint32 numthreads,
					    bool masknastyrepeats,
					    std::vector<uint32> & minkmer,
					    uint32 mincountkmerforks
    );

  //void showHashStatisticsInfo(std::string & filename);
  void showHashStatisticsInfo();

  void loadHashStatistics(ReadPool & rp, std::string & hashstatfilename, uint8 basesperhash);

  // use this one for multi-thread baiting on the same HashStatistics object
  uint32 checkBaitHit(Read & actread, std::vector<vhrap_t> & baiting_singlereadvhraparray, std::vector<uint8> & baiting_tagmaskvector);
  // use this one for single thread baiting on the same HashStatistics object
  inline uint32 checkBaitHit(Read & actread) {return checkBaitHit(actread,HS_baiting_singlereadvhraparray,HS_baiting_tagmaskvector);}

  const hashstat_t * findVHash(const hashstat_t & searchval);

  void digiNormReset() { HS_diginorm_count.clear();}
  bool digiNormTestRead(Read & actread, bool force);
  uint32 estimDigiNormCov(Read & actread);
};


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


class NHashStatistics
{
#define HS_MAXVHASHMASK 0xFFFFFFLL
public:
  // new hashstat
  struct hscounts_t {
    uint8   lowposd4;      // lowest hash position/4 (lowest pos: 0, 255 means >=1020)
    uint32  fcount:24;
    bool    b1:1;          // 4 bits, unused atm
    bool    b2:1;
    bool    b3:1;
    bool    b4:1;
    uint8   seqtype:4;     // sequencing type
    uint32  rcount:24;

    hscounts_t() {
      uint64 * ptr=reinterpret_cast<uint64 *>(&lowposd4);
      *ptr=0;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const hscounts_t & v){
      ostr << "fcount: " << v.fcount
	   << "\trcount: " << v.rcount
	   << "\tlowposd4: " << static_cast<uint32>(v.lowposd4)
	   << "\tseqtype: " << static_cast<uint16>(v.seqtype)
	;
      return ostr;
    }
  };


  struct nhashstat_t {
    vhash_t    vhash;       // vhash,
    hscounts_t hsc;

    nhashstat_t() {
      vhash=0;
    }

    friend std::ostream & operator<<(std::ostream &ostr, const nhashstat_t & v){
      ostr << "hash: " << std::hex << v.vhash << std::dec
	   << '\t' << v.hsc
	;
      return ostr;
    }


  };

private:
  // begin and end iterators
  struct hsvbendit_t {
    std::vector<nhashstat_t>::const_iterator b;
    std::vector<nhashstat_t>::const_iterator e;
  };

private:

  static uint32 HSN_hs_magic;

  BloomFilter * HSN_bloomfilter;

  std::vector<nhashstat_t>                  HSN_hsv_hashstats;
  std::vector<hsvbendit_t>                  HSN_hsv_hsshortcuts;

  std::unordered_map<vhash_t,hscounts_t>    HSN_hsum_hashstats;

  std::vector<uint64>                                    HSN_hs_dist;
  uint8 HSN_hs_sortstatus;
  bool HSN_hs_needsconsolidation;

  uint8 HSN_basesperhash;

  std::string HSN_directory; // dir to write hashstat to

  uint16 HSN_needsteps;
  uint16 HSN_step;


private:
  void makeNHashStatArrayShortcuts(std::vector<nhashstat_t> & nhashstats,
				   const uint8 basesperhash,
				   std::vector<hsvbendit_t> & hsshortcuts);
  inline static bool sortHashStatComparatorByLow24bit(const nhashstat_t & a, const nhashstat_t & b){
    if((a.vhash & HS_MAXVHASHMASK) != (b.vhash & HS_MAXVHASHMASK)) {
      return (a.vhash & HS_MAXVHASHMASK) < (b.vhash & HS_MAXVHASHMASK);
    }
    return a.vhash < b.vhash;
  }
  inline static bool sortHashStatComparator(const nhashstat_t & a, const nhashstat_t & b){
    return a.vhash < b.vhash;
  }
  void sortLow24Bit(std::vector<nhashstat_t> & hashstats, uint8 & sortstatus) {
    std::sort(hashstats.begin(), hashstats.end(), sortHashStatComparatorByLow24bit);
    sortstatus=1;
    HSN_hsv_hsshortcuts.clear();
  }
  void sortLexicographically(std::vector<nhashstat_t> & hashstats, uint8 & sortstatus) {
    std::sort(hashstats.begin(), hashstats.end(), sortHashStatComparator);
    sortstatus=2;
    HSN_hsv_hsshortcuts.clear();
  }
  uint64 calcHashDistrib(std::vector<uint64> & hsdist);

  void learnSequenceQuick1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse, bool lookuponly);
  void learnSequenceQuick2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep1(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep2(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void learnSequenceStep3(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);

  void saveHashVStatistics(std::ostream & ostr);
  void saveHashMStatistics(std::ostream & ostr);

  void trimHashVStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);
  void trimHashMStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);

public:
  NHashStatistics() : HSN_bloomfilter(nullptr), HSN_basesperhash(0), HSN_hs_sortstatus(0) {};
  ~NHashStatistics();

  void setupNewAnalysis(const uint8  bfbits,
			const uint32 bfnumkeys,
			const uint8  basesperhash,
			uint16 numsteps);
  void analyseReadPool(ReadPool & rp);
  void deleteBloomFilter();

  //void prefetchSequence(const void * seqvoid, uint64 slen, const char * namestr);
  void learnSequence(const void * seqvoid, uint64 slen, const char * namestr, uint8 seqtype, bool isreverse);
  void finaliseStep();

  const nhashstat_t * findVHash(const nhashstat_t & searchval);

  void trimHashStatsByFrequency(int32 minfwd, int32 minrev, int32 mintotal);

  std::vector<uint64> & getHashDistrib();
  void dumpHashDistrib(std::ostream & ostr);

  void dumpHealth(std::ostream & fout);
  void dumpHashCount(std::ostream & ostr);

  void saveHashStatistics(std::ostream & ostr);
  void loadHashStatistics(std::istream & istr);

  void saveHashStatistics(const std::string & filename, bool deleteoldfile);
  void loadHashStatistics(const std::string & filename);

  static void hash2string(vhash_t hash, uint8 basesperhash, std::string & str);
  static inline std::string hash2string(vhash_t hash, uint8 basesperhash) { std::string tmp; hash2string(hash,basesperhash,tmp); return tmp;}

  void sortLow24Bit(){
    sortLow24Bit(HSN_hsv_hashstats, HSN_hs_sortstatus);
  }
  void sortLexicographically(){
    sortLexicographically(HSN_hsv_hashstats, HSN_hs_sortstatus);
  }

  void moveStatCountMapToVector();


};


#endif
