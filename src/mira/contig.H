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


#ifndef _mira_contig_h_
#define _mira_contig_h_

#include <iostream>
#include <iomanip>
#include <string>

#include <boost/unordered_set.hpp>

#include "stdinc/stlincludes.H"

#include "stdinc/defines.H"

#include "mira/align.H"
#include "mira/hdeque.H"
#include "mira/parameters.H"
#include "mira/pcrcontainer.H"
#include "mira/readpool.H"
#include "mira/sam_collect.H"

#include "util/dptools.H"






// structure used to give back information on GenBank tags in the contig
struct gbfsummary_t {
  PlacedContigReads::const_iterator pcrI;          // pcrI of read with this feature

  uint32 cfrom;              // position of feature in contig
  uint32 cto;

  int32  strainid;
  std::string identifier;  // TODO: switch that to multitag_t::mte_id_t ??
  std::string locustag;
  std::string gene;
  std::string function;
  std::string ecnumber;
  std::string product;
  std::string note;
  int8   translationtable;
  int8   codonstart;
  int8   direction;
  bool   mustbetranslated;
  bool   isgene;             // Fgen, FCDS or one of the main RNAs

  friend std::ostream & operator<<(std::ostream &ostr, const gbfsummary_t & gbfs){
    ostr << "cfrom: " << gbfs.cfrom
	 << "\ncto: " << gbfs.cto
	 << "\npcrI: " << gbfs.pcrI
	 << "\nstrainid: " << gbfs.strainid
	 << "\nidentifier: " << gbfs.identifier
	 << "\nlocustag: " << gbfs.locustag
	 << "\ngene: " << gbfs.gene
	 << "\nfunction: " << gbfs.function
	 << "\necnumber: " << gbfs.ecnumber
	 << "\nproduct: " << gbfs.product
	 << "\nnote: " << gbfs.note
	 << "\ntranslationtable: " << static_cast<uint16>(gbfs.translationtable)
	 << "\ncodonstart: " << static_cast<uint16>(gbfs.codonstart)
	 << "\ndirection: " << static_cast<uint16>(gbfs.direction)
	 << "\nmustbetranslated: " << gbfs.mustbetranslated
	 << "\nisgene: " << gbfs.isgene << std::endl;
    return ostr;
  }

  gbfsummary_t(PlacedContigReads::const_iterator & pcrIp) : pcrI(pcrIp), cfrom(-1), cto(-1), strainid(-1),translationtable(11),codonstart(1),direction(0),mustbetranslated(false),isgene(false) {};
};


// structure used to give back information on SequenceOntology tags in the contig
struct contigSOtag_t {
  multitag_t multitag;                // positions are now contig positions!!!
  int32  readindex;          // index of read with this feature in CON_reads

  friend std::ostream & operator<<(std::ostream &ostr, const contigSOtag_t & csot){
    ostr << "mt: " << csot.multitag
	 << "readindex: " << csot.readindex << std::endl;
    return ostr;
  }
};


struct coverageinfo_t {
  uint64  min;
  uint64  max;
  double  mean;
  double  median;
  double  stddev;

  coverageinfo_t() {clear();}
  void clear() {min=0xffffffffffffffff; max=0; mean=0; median=0; stddev=0;}

  friend std::ostream & operator<<(std::ostream &ostr, const coverageinfo_t & ci){
    ostr << "CINFO: " << ci.min << '\t' << ci.max
	 << '\t' << ci.mean
	 << '\t' << ci.median
	 << '\t' << ci.stddev
	 << std::endl;
    return ostr;
  }
};

struct tagcoverageinfo_t {
  contigSOtag_t  csot;
  coverageinfo_t ccinfo;

  std::string comparator_text;
  double comparator_factor;
};


class Contig
{
public:
  enum{AS_TEXT=0, AS_HTML, AS_CAF, AS_MAF, AS_ACE, AS_FASTA, AS_FASTAQUAL,
	 AS_FASTAPADDED, AS_FASTAPADDEDQUAL, AS_GAP4DA, AS_TCS,
	 AS_DEBUG};

  enum addreadcode{ENOERROR=0,
		   ENOTCALLED,
		   ENOALIGN,
		   EDROPINRELSCORE,
		   ETEMPLATEDIRECTION,
		   ETEMPLATESIZELT,
		   ETEMPLATESIZEGT,
		   ESEGMENTPLACEMENT,
		   ESRMBMISMATCH,
		   EDANGERZONE,
		   ESPECIALSRADDFAIL,
		   EREFIDNOTALLOWED,
		   EMAXCOVERAGEREACHED,
		   EFORCEDGROWTHNOTREACHED,
		   EGROWTHNOTALLOWED,
		   EUNSPECIFIED
  };




// Structure used for initialising a complete contig in initialiseContig(...)
  struct contig_init_read_t {
    int32 id;
    int32 offset_in_contig;
    int32 read_lclip;
    int32 read_rclip;
    int32 direction;
  };



// we can currently merge 2 seqtypes: SOlexa, SOLiD
// SOLiD currently not used
#define NUMMERGESEQTYPES 2

  typedef uint32 ccctype_t;  // consensus_counts_count_type

  struct consensus_counts_t{
    ccctype_t A;             // ACGT are extended counters
    ccctype_t C;             // ccctype_t is enough for a coverage of 16k reads.
    ccctype_t G;
    ccctype_t T;

    ccctype_t N;             // N, X, star and coverage are normal counters.
    ccctype_t X;
    ccctype_t star;

    ccctype_t seqtype_cov[ReadGroupLib::SEQTYPE_END];  /* for each seqtype, count the
					     * coverage at this position. This
					     * number also includes the 100%
					     * mapped reads!
					     */
    ccctype_t total_cov;

    uint16 baselock;     // if > 0 then only one of A, C, T, G may be set, otherwise it's a misassembly

    uint16 snplock;

    // New for Short read mappings, but will also be used for normal
    //  backbone mappings
    //

    uint8 forcemergearea;  /* if >0, short reads falling in this area will always
			      be merged, even if not mapping at 100% */


    // Currently, these MUST be at the end of this struct:
    //  MIRA depends on it by using CON_concounts_zero and
    //  CON_concounts_zero_nobb to initialise quickly or rebuilding
    //  a consensus
    ccctype_t bbcounts[NUMMERGESEQTYPES];    /* coverage of reads supporting 100%
					     backbone 0 is Solexa, 1 is ABI
					     SOLiD */
    base_quality_t bbbestquals[NUMMERGESEQTYPES]; /* best quality of reads
						     supporting 100% backbone
						     0 is Solexa, 1 is ABI
						     SOLiD */
    uint8 bbstrains[NUMMERGESEQTYPES];   /* bitmap. which strains are mapped
					    at this pos max. 8 strains
					    possible. Strains are in 76543210
					    order */
    char   i_backbonecharorig;  /* original backbone consensus character (uppercase)
				   gaps are *, positions not initially covered by backbone
				   '@' */
    char   i_backbonecharupdated;  /* updated backbone consensus character (uppercase)
				      if '@', take i_backbonecharorig, else this */

    inline char getBBChar() const {
      if(i_backbonecharupdated!='@') return i_backbonecharupdated;
      return i_backbonecharorig;
    }
    inline char getOriginalBBChar() const {
      return i_backbonecharorig;
    }

    // that's more for debugging
    friend std::ostream & operator<<(std::ostream &ostr, const consensus_counts_t & cc){
      std::cout << "A: " << cc.A
	   << "\tC: " << cc.C
	   << "\tG: " << cc.G
	   << "\tT: " << cc.T
	   << "\tN: " << cc.N
	   << "\tX: " << cc.X
	   << "\t*: " << cc.star
	   << "\ttcov: " << cc.total_cov
	   << "\tblock: " << cc.baselock
	   << "\tslock: " << cc.snplock
	//<< "\n  arr"
	   << "\tbbco: " << cc.i_backbonecharorig
	   << "\tbbcu: " << cc.i_backbonecharupdated
#if CPP_READ_SEQTYPE_END != 8
#error "This code is made for 8 sequencing types, adapt!"
#endif
	   << "\tcSAN: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_SANGER]
	   << "\tc454: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_454GS20]
	   << "\tcION: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_IONTORRENT]
	   << "\tcPBH: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_PACBIOHQ]
	   << "\tcPBL: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_PACBIOLQ]
	   << "\tcTXT: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_TEXT]
	   << "\tcSXA: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_SOLEXA]
	   << "\tcSID: " << cc.seqtype_cov[ReadGroupLib::SEQTYPE_ABISOLID]
	   << "\tbcSXA: " << cc.bbcounts[0]
	   << "\tbcSID: " << cc.bbcounts[1]
	   << "\tbqSXA: " << static_cast<ccctype_t>(cc.bbbestquals[0])
	   << "\tbqSID: " << static_cast<ccctype_t>(cc.bbbestquals[1])
	   << "\tsSXA: " << static_cast<ccctype_t>(cc.bbstrains[0])
	   << "\tsSID: " << static_cast<ccctype_t>(cc.bbstrains[1]);
      return ostr;
    }
  };

  // Contig exports the container with consensus_counts_t
  //  as type:  Contig::cccontainer_t

  typedef HDeque<consensus_counts_t> cccontainer_t;
  //typedef IndexedDeque<consensus_counts_t> cccontainer_t;




// we need an own tag type for the consensus as additional data must be
//  stored
// lets use simple inheritance for that (first time for MIRA where
//  it is useful :-)
// Consensus tags have a quality identifier for each base (and *)
//  this is useful for example for SRMBs and SNPs
// This thing is a class of its own, especially
//  due to the copy operator for "upgrading" normal tags
// Lazy: no good constructors/destructors
// Contrary to the groups_t struct (farther below), only bases
//  which form a valid group have their quality set >0 here.

  class consensustag_t: public multitag_t {
  public:
    bool additionalinfo_initialised;

    base_quality_t qualACGTStar[5];  // qualities for A, C, G, T, *

    //base_quality_t qualA;
    //base_quality_t qualC;
    //base_quality_t qualG;
    //base_quality_t qualT;
    //base_quality_t qualStar;

/*
    // TODO: find out where used, multitag_t???
    // upgrading copy operator: when copying a tag_t to contigtag_t
    Contig::consensustag_t const & operator=(tag_t const & other){
      from=other.from;
      to=other.to;
      identifier=multitag_t::newIdentifier(other.identifier);
      comment=multitag_t::newComment(other.comment);
      setStrand(other.strand);
      phase=3;

      additionalinfo_initialised=false;

      qualACGTStar[0]=0;
      qualACGTStar[1]=0;
      qualACGTStar[2]=0;
      qualACGTStar[3]=0;
      qualACGTStar[4]=0;

      //qualA=0;
      //qualC=0;
      //qualG=0;
      //qualT=0;
      //qualStar=0;

      return *this;
    }
*/

    Contig::consensustag_t const & operator=(multitag_t const & other){
      //from=other.from;
      //to=other.to;
      //identifier=other.identifier;
      //comment=other.comment;
      //strand=other.strand;
      //phase=other.phase;

      multitag_t::operator=(other);

      additionalinfo_initialised=false;

      qualACGTStar[0]=0;
      qualACGTStar[1]=0;
      qualACGTStar[2]=0;
      qualACGTStar[3]=0;
      qualACGTStar[4]=0;

      //qualA=0;
      //qualC=0;
      //qualG=0;
      //qualT=0;
      //qualStar=0;

      return *this;
    }

    // lazy constructor (need to define as copy-constructor is defined)
    consensustag_t() {
      additionalinfo_initialised=false;

      qualACGTStar[0]=0;
      qualACGTStar[1]=0;
      qualACGTStar[2]=0;
      qualACGTStar[3]=0;
      qualACGTStar[4]=0;
    }

  };

  struct constats_t {
    bool statsvalid;

    bool contains_long_repeats_only;  // copy of CON_contains_long_repeats_only for saving
    int8 islargecontig; // -1=no; 0=undefined; 1=yes

    uint32 conlength;
    uint32 conlength_nogap;

    uint32 AinC;
    uint32 CinC;
    uint32 GinC;
    uint32 TinC;
    uint32 NinC;
    uint32 XinC;
    uint32 starInC;
    uint32 IUPACinC;
    uint32 FunnyInC;
    uint32 NinR;
    uint32 XinR;
    uint32 starInR;

    double gccontent;

    uint32 numSRMc;
    uint32 numWRMc;
    uint32 numSTMU;
    uint32 numSTMS;

    uint32 total_reads;
    uint32 readsperst[ReadGroupLib::SEQTYPE_END];    // EXCLUDING rail reads!!!
    uint64 totalbasesperst[ReadGroupLib::SEQTYPE_END];
    uint32 numreads_noqual;
    uint32 numreads_withqual;

    uint32 numnocoverage;

    uint32 max_coverage;
    uint32 max_covperst[ReadGroupLib::SEQTYPE_END];

    double avg_coverage;
    double avg_covperst[ReadGroupLib::SEQTYPE_END];

    uint32 avg_conqual;
  };


  struct errorstatus_t {
    int32         code;                   // the addread error code
    std::vector<int32> reads_affected;

    inline errorstatus_t() : code(ENOERROR) {};
    inline void reset() { code=ENOERROR; reads_affected.clear(); }
    void dumpStatus(bool longmsg=false, const char * additionalmsg=nullptr);
  };


  struct templateguessinfo_t {
    uint32 tsize_seen;   // outer template size seen
    int8  splace_seen;  // segment placement deduced
    ReadGroupLib::ReadGroupID rgid;

    friend std::ostream & operator<<(std::ostream &ostr, const templateguessinfo_t & tgi){
      ostr << "tgi rgid: " << tgi.rgid.getLibId()
	   << "\t" << static_cast<int16>(tgi.splace_seen)
	   << "\t" << tgi.tsize_seen;
      return ostr;
    }
  };



  // pacbio dark strobe edit
  struct pbdse_t {
    uint32 rid;               // rid in readpool
    int32  rposs;              // between read pos
    int32  rpose;              // and read pos
    int32  changeestim;        // >0: insert Ns; <0 delete Ns

    friend std::ostream & operator<<(std::ostream &ostr, const pbdse_t & pbdse){
      ostr << "rid: " << pbdse.rid << "\trposs: " << pbdse.rposs << "\trpose: " << pbdse.rpose << "\tce: " << pbdse.changeestim << std::endl;
      return ostr;
    }
  };

  // static variables
public:
  // Note: due to possible static initialization fiasco,
  //       the REA_tagentry_* variables are defined/initialised
  //       in "mira/multitag.C" and not in "mira/contig.C"
  static const multitag_t::mte_id_t CON_tagentry_idEmpty;
  static const multitag_t::mte_id_t CON_tagentry_idALUS;
  static const multitag_t::mte_id_t CON_tagentry_idREPT;
  static const multitag_t::mte_id_t CON_tagentry_idSRMc;
  static const multitag_t::mte_id_t CON_tagentry_idWRMc;
  static const multitag_t::mte_id_t CON_tagentry_idSAOc;
  static const multitag_t::mte_id_t CON_tagentry_idSROc;
  static const multitag_t::mte_id_t CON_tagentry_idSIOc;
  static const multitag_t::mte_id_t CON_tagentry_idPSHP;
  static const multitag_t::mte_id_t CON_tagentry_idED_D;
  static const multitag_t::mte_id_t CON_tagentry_idED_C;
  static const multitag_t::mte_id_t CON_tagentry_idED_I;
  static const multitag_t::mte_id_t CON_tagentry_idESDN;
  static const multitag_t::mte_id_t CON_tagentry_idSTMS;
  static const multitag_t::mte_id_t CON_tagentry_idSTMU;
  static const multitag_t::mte_id_t CON_tagentry_idUNSc;
  static const multitag_t::mte_id_t CON_tagentry_idIUPc;
  static const multitag_t::mte_id_t CON_tagentry_idMCVc;
  static const multitag_t::mte_id_t CON_tagentry_idDGPc;

  static const multitag_t::mte_id_t CON_tagentry_idSOFApolyA_signal_sequence;

  static const multitag_t::mte_co_t CON_tagentry_coEmpty;

  // special flag for debugging sessions. If set, the code should
  //  abort via a sigsev so that a debugger can catch it
  // should never ever happen in distributed packes, always commment out
  //  the code which performs the bombing
  static bool CON_abortflag;

private:
  static const consensus_counts_t CON_concounts_zero;
  static const consensus_counts_t CON_concounts_zero_nobb;

  static uint32 CON_id_counter;
  static bool CON_static_ok;

  static uint32 CON_railcounter;

  static uint32 CON_cer_idcounter;     // Coverage Equivalent Read id counter

  //static std::vector< std::string >   CON_danger_zones_ids;
  //static std::vector< std::string >   CON_baselock_ids;
  //static std::vector< std::string >   CON_snplock_ids;

  static std::vector<multitag_t::mte_id_t>   CON_danger_zones_ids;
  static std::vector<multitag_t::mte_id_t>   CON_baselock_ids;
  static std::vector<multitag_t::mte_id_t>   CON_snplock_ids;



  static uint8 CON_outtype;

  static bool  CON_outputrails;

  static boost::unordered_set<std::string> CON_cebugnames;
  static bool CON_mastercebugflag;
  //Variables
private:
  bool CON_cebugflag; // only for special debugging in getConsensus() etc.

  std::vector<MIRAParameters> * CON_miraparams;


  uint32 CON_id;
  std::string CON_nameprefix;

  mutable std::string CON_name;  /* if this string is non-empty, then it's the name of
			       the contig. If it is empty, the name is made out of
			       the name prefix, some characters and the id.
			       see getContigName() for exact name scheme */

  bool CON_finalised;

  ReadPool * CON_readpool;
  PlacedContigReads         CON_reads;
  cccontainer_t             CON_counts;
  std::multiset<int32>           CON_templates_present;
  std::vector<consensustag_t>    CON_consensus_tags;

  // if set and if an element >0, use this for checking free coverage
  //  when trying to add a read
  std::vector<ccctype_t>            CON_targetcoverageperst;


  /*
    This vector is here to allow selective beef up the contig when mapping
    short reads.
    If active, it has the size of the readpool (CON_readpool) and each
    element is either 1 (allowed) or 0 (not allowed), telling which
    reads are allowed or not to be a reference read when adding new
    reads to a contig.
    If dirty, size is 0.
   */
  std::vector<bool> CON_allowedrefids;


  //// fast cons function cached
  //uint32 CON_len_tmpcons;
  //char * CON_tmpcons;

  std::string CON_2tmpcons;

  /* Two flags that tell the contig how to add new reads:
   */

  bool CON_tmpcons_from_backbone;   /* if true, builds tmp consensus from
				       backbone chars if available */
  bool CON_mergenewsrreads;         /* if true, tries to use read merging
				       for short reads (Solexa, SOLiD */
  bool CON_hasforcemergeareas;      /* if true, some areas in the backbone
				       ask for forced merging */

  // Special conditions when aligning short reads (Solexa, SOLiD)
  bool CON_specialsraddconditions;   /* if true, applies special checks
					when adding short reads */
  int32 CON_ssrc_maxtotalerrors; // = to the sum of following two
  int32 CON_ssrc_maxgaps;        /* max number of gaps allowed for adding a
				    short read when using special SR
				    conditions */
  int32 CON_ssrc_maxmismatches;  /* max number of mismatches allowed for
				    adding a short read when using special SR
				    conditions */

  uint32 CON_nummergedreads_perseqtype[NUMMERGESEQTYPES];

  //
  // tracker for (clipped) longest non-rail read which was ever added to contig
  //
  // needed for efficient finding of lower_bound pcrIs in the CON_reads PlacedContigReads
  //
  // Note: should reads in the future be removed via a call to deleteRead() outside
  //  of the whole addRead() functionality, this value might not track the longest
  //  clipped read in the contig, but the longest reads which was ever added.
  //
  int32 CON_longestreadseen;
  int32 CON_longestrailseen;  // same, but only for rails
  int32 CON_longestnonbbreadseen;  // same as longestreadseee, but without backbones

  // counter for num reads per strain. kept up to date by contig
  //
  // BaCh 21.01.2012: with the new readsperreadgroup below, this has become
  //  a bit redundant. Keep for the time being, maybe faster in queries
  std::vector<uint32> CON_readsperstrain;

  // counter for num reads per readgroup. kept up to date by contig
  std::vector<uint32> CON_readsperreadgroup;

  // tracker for marker positions
  // MPs are positions which can be freely defined and where the contig will adjust
  //  the positions while growing the contig
  // Used atm only for con_msr_keependsunmapped (left and right MP), so simplistic version
  std::vector<int32> CON_markerpositions;
  // these two point into CON_markerpositions for initial positions 0 (left) and
  //  and "contigsize" (right)
  int32         CON_mpindex_msrkceu_left;     // a value of -1 equals to "unused"
  int32         CON_mpindex_msrkceu_right;    //

  /*
    The "fixed" consensus is a consensus that the contig gets initialised with
    (e.g. loading from CAF). If set, then fetching consensus sequence and
    qualities will yield this instead of dynamically generated values from the
    CON_allcons* and CON_straincons* variables below.

    Gets trashed when the contig changes)
   */
  std::string                 CON_fixedconsseq;
  std::vector<base_quality_t> CON_fixedconsqual;

  // new
  // when consenus is calculated, store here for every strain
  // it will be up to the caller to trash these values for recalc with other parameters
  // exception: if contig gets changed (added read, removed read etc.pp), these get trashed
  //  automatically by implicitly called definalise()
  // strains with no read in a given contig will still have the vectors filled (char '@', qual 0)

  int32          CON_conscalc_mincov;
  base_quality_t CON_conscalc_minqual;
  char           CON_conscalc_missingchar;

  std::string                  CON_allconsseq;
  std::vector<base_quality_t>  CON_allconsqual;
  std::vector<int32>           CON_alladjustments;

  std::vector<std::string>                  CON_strainconsseq;
  std::vector<std::vector<base_quality_t> > CON_strainconsqual;
  std::vector<std::vector<int32> >          CON_strainadjustments;



  std::vector<int32> CON_last_dangerous_overlaps;


  constats_t CON_stats;


  struct edit454command_t {
    PlacedContigReads::const_iterator pcrI;
    uint32 urdid;                            // just needed for comparator
    uint32 readpos;
    char base;

    inline edit454command_t(PlacedContigReads::const_iterator & pcrIp, uint32 urdidp, uint32 readposp, char basep) :
      pcrI(pcrIp), urdid(urdidp), readpos(readposp), base(basep) {};
  };


  // set by pathfinder if the assembly mode is to keep long repeats
  //  separated and this contig started with a multicopy read
  // used for definig contig name
  bool CON_contains_long_repeats_only;

  // I don't think I'll include a tribool for that :-)
  // -1 = no, 0=unknown, 1=yes
  int8 CON_contains_majority_digitallynormalised_reads;

  enum{USCLO_PRE=0,
       USCLO_DIR,
       USCLO_XCUT,
       USCLO_TEMPL1,
       USCLO_SWALIGN,
       USCLO_PREINSCHK,
       USCLO_INSCON,
       USCLO_INSCONM,
       USCLO_UPDBLOCKS,
       USCLO_DELREAD,
       USCLO_ANRMBZ,
       USCLO_GRACP,
       USCLO_END};

  // track timing
  std::vector<suseconds_t> CON_us_steps;

#define CLOCK_STEPS_IRIC
  enum{USCLOIRIC_INSGLCC=0,
       USCLOIRIC_INSGLARO,
       USCLOIRIC_INSGLACT,
       USCLOIRIC_INSGLTOT,
       USCLOIRIC_PR,
       USCLOIRIC_TEMPL,
       USCLOIRIC_INDEX,
       USCLOIRIC_BIGLCCINS,
       USCLOIRIC_BIGLINTERPOL,
       USCLOIRIC_BIGLUPDTAGS,
       USCLOIRIC_BIGLFPCRI,
       USCLOIRIC_BIGLLIGAPINREADS,
       USCLOIRIC_BIGLSHIFTREADS,
       USCLOIRIC_BIGLTOT,
       USCLOIRIC_INSGR,
       USCLOIRIC_UCV,
       USCLOIRIC_END
  };

#define CLOCK_STEPS_DRFC
  enum{
    USCLODRFC_UBL,
    USCLODRFC_UCV,
    USCLODRFC_ITF,
    USCLODRFC_ITB,
    USCLODRFC_CCEF,
    USCLODRFC_CCEB,
    USCLODRFC_SDT, // shift delete tags
    USCLODRFC_SR,
    USCLODRFC_SMP,
    USCLODRFC_DT,
    USCLODRFC_RR,
    USCLODRFC_DSOC,
    USCLODRFC_TOTAL,
    USCLODRFC_END
  };


//#define CLOCK_STEPS_CONS  // careful with this one, 30% speed impact
  enum{
    USCLOCONS_H1_MGROUPS,
    USCLOCONS_H1_CP2URP,
    USCLOCONS_H1_PCRI,
    USCLOCONS_H1_EGROUP,
    USCLOCONS_H1_GQUAL,
    USCLOCONS_H1_CALLH2,
    USCLOCONS_H1_TOTAL,
    USCLOCONS_END
  };


  // track timing for insert read in contig
  std::vector<suseconds_t> CON_us_steps_iric;
  // track timing for delete read from contig
  std::vector<suseconds_t> CON_us_steps_drfc;
  // track timing for consensus generation
  std::vector<suseconds_t> CON_us_steps_cons;

  // track number of delete calls
  size_t CON_track_numins;
  size_t CON_track_numdels;

public:


  //Functions
private:
  void setCEBUGFlag(int32 id1, int32 id2);
  void foolCompiler();
  void init();
  void zeroVars();
  void definalise();

  //bool outputOrderSortCmp_(const out_order_t & a, const out_order_t & b);

  void addRead_wrapped(
    std::vector<Align> & aligncache,
    const AlignedDualSeqFacts * initialadsf,
    int32 refid,
    int32 newid,
    int32 direction_frnid,
    bool  newid_ismulticopy,
    int32 forcegrow,
    templateguessinfo_t & templateguess,
    errorstatus_t & errstat);

  void updateCountVectors(const int32 from,
			  const int32 len,
			  std::vector<char>::const_iterator updateseq,
			  const uint32 seqtype,
			  bool bla,
			  int32 coveragemultiplier);
  void rebuildConCounts();
  PlacedContigReads::const_iterator insertReadInContig(const AlignedDualSeq & ads,
						       uint32 coffset,
						       int32 direction_frnid,
						       int32 direction_refid,
						       int32 coveragemultiplier);
  bool insertMappedReadInContig(const AlignedDualSeq & ads,
				const uint8 newreadseqtype,
				const uint32 coffset,
				const int32 direction_frnid,
				const int32 coveragemultiplier,
				const bool forcemerge);
  int32 analyseADSForCuttingBackCERMap(const AlignedDualSeq & ads, int32 direction_frnid);

  void finalise();
  void checkContig();
  void updateTagBaseInserted(uint32 contigpos);
  void updateTagBaseDeleted(uint32 contigpos);


  void shiftMarkerPositions(int32 pos, int32 offset);

  // TODO: see whether can be replace by the looking at sorted multitag_t vector?
  static bool consensustag_t_comparator(const consensustag_t & t1,
					const consensustag_t & t2)
    {
      if(t1.from < t2.from) return true;
      if(t1.from > t2.from) return false;
      if(t1.to < t2.to) return true;
      if(t1.to > t2.to) return false;
      return (t1.getIdentifierStr() < t2.getIdentifierStr());
    };


  // rid high to low, within rid startpos high to low
  static bool pbdse_t_comparator(const pbdse_t & p1,
				 const pbdse_t & p2)
    {
      if(p1.rid == p2.rid){
	return (p1.rposs < p2.rposs);
      }
      return p1.rid < p2.rid;
    };

  // sort ascending, first by contigfrom, then to
  //  then sort by identifier: Fsrc, Fgen, FCDS, rest
  static bool gbfsummary_t_comparator(const gbfsummary_t & a,
				      const gbfsummary_t & b)
    {
      if(a.cfrom!=b.cfrom) return (a.cfrom<b.cfrom);
      if(a.cto!=b.cto) return (a.cto<b.cto);
      if(a.identifier!=b.identifier) {
      	if(a.identifier.empty()) return true;
      	if(b.identifier.empty()) return false;
      	if(a.identifier=="Fsrc") return true;
      	if(b.identifier=="Fsrc") return false;
      	if(a.identifier=="Fgen") return true;
      	if(b.identifier=="Fgen") return false;
      	if(a.identifier=="FCDS") return true;
      }

      return false;
    };


  static bool edit454command_t_comparator(const edit454command_t & a,
					  const edit454command_t & b)
    {
      if(a.urdid == b.urdid) {
	// sort within a read by DESCENDING! readpos order
	return a.readpos > b.readpos;
      }
      return a.urdid < b.urdid;
    };

  ///////////////////////////////////////////////////////

  class nngroups_t {
  public:
    char base;
    base_quality_t groupquality;

    bool valid;

    uint32 forwarddircounter;
    uint32 complementdircounter;

    // these vectors are 'linked': nth element of one corresponds to nth of another

    // OUCH: PlacedContigReads need rethinking IDs! (maybe export unsorted ids of
    //  ReadContainer in PlacedContigReads?

    std::vector<int32> urdids;
    std::vector<base_quality_t> quals;
    std::vector<int8>  directions;

    friend std::ostream & operator<<(std::ostream &ostr, const nngroups_t &nngroup){
      ostr << "    base " <<  nngroup.base << '\n';
      ostr << "    valid " <<  nngroup.valid << '\n';
      ostr << "    forwarddircounter " <<  nngroup.forwarddircounter << '\n';
      ostr << "    complementdircounter " <<  nngroup.complementdircounter << '\n';
      ostr << "    groupquality " <<  static_cast<uint32>(nngroup.groupquality) << '\n';
      ostr << "    urdids.size " <<  nngroup.urdids.size() << '\n';
      ostr << "    quals.size " <<  nngroup.quals.size() << '\n';
      ostr << "    directions.size " <<  nngroup.directions.size() << '\n';
      return ostr;
    }
  };

  ////////////////////////////////////////////////////////////////////////

  struct pos_rep_col_t {
    int32 contigpos;
    int32 numvalidgroups;
    std::vector<nngroups_t> groups;
    //std::vector<int32> allids;   not needed anyway (just for lazy debugging output)
    std::string type;

    bool is_dangerous;
    bool tagged;
  };
  void mpr__tagColumn(pos_rep_col_t & prc);
  ///////////////////////////////////////////////////////


  //void calcGroupQual(const groups_t & g);


  void buildMaskShadow(std::vector<int8> & maskshaddow,
		       std::vector<multitag_t::mte_id_t> & masktagtypes,
		       bool onlybackbone);
  void zeroStats();
  void calcStats();
  void priv_calcNumDigitallyNormalisedReads();

  class rcci_t {
    // ReadColContigIterator
    //
    // a vector (read_ids_in_col) keeps the ids of the CON_reads that are
    //  covering a specific position of the contig
    // reserving 1000 position should be enough 99.9% of all cases,
    //  is automatically extended by STL if needed.

    Contig * RCCI_contig;

    // read pcrI per sequencing type per strain
    std::vector<PlacedContigReads::const_iterator> RCCI_pcrIs_in_col;
    // iteraor to the new PCRIs which appeared while calling
    //  the last update().
    // in conjunction with advance, works only if advance(n=...) has
    //  n==1
    std::vector<PlacedContigReads::const_iterator>::iterator RCCI_newPCRIsonlastupdate;

    //std::vector<out_order_t>::const_iterator outorderI;
    PlacedContigReads::const_iterator RCCI_mpcrI;
    uint32 RCCI_actcontigpos;

    std::vector<int32> RCCI_allowedstrainids;
    std::vector<uint8> RCCI_allowedreadtypes;

    bool RCCI_takerails;
    bool RCCI_takebackbones;
    bool RCCI_takereadswithoutreadpool;

  public:
    rcci_t(Contig * cptr) : RCCI_contig(cptr), RCCI_mpcrI(cptr->CON_reads.begin()) {};

    void init(const std::vector<int32> & allowedstrainids,
	      const std::vector<uint8> & allowedsequencingtypes,
	      bool takerails,
	      bool takebackbones,
	      bool takereadswithoutreadpool);
    void update();
    inline void advance(uint32 dist=1) {
      for(auto disti=dist; disti>0; --disti){
	++RCCI_actcontigpos;
	update();
      }
      if(dist!=1) RCCI_newPCRIsonlastupdate=RCCI_pcrIs_in_col.end();
    }

    inline const std::vector<PlacedContigReads::const_iterator> & getPCRIsInCol() const { return RCCI_pcrIs_in_col;}
    inline std::vector<PlacedContigReads::const_iterator>::iterator getItToNewPCRIs() const { return RCCI_newPCRIsonlastupdate;}
    inline uint32 getContigPos() const { return RCCI_actcontigpos;}
  };

  class ercci_t {
    // extendedReadColContigIterator
    //
    // per sequencing type and per strain ids,
    //  a vector keeps the ids of the reads that are
    //  covering a specific position of the contig
    // reserving 1000 position should be enough 99.9% of all cases,
    //  is automatically extended by STL if needed.

    Contig * ERCCI_contig;

    // read pcrI per sequencing type per strain
    std::vector<std::vector<std::vector<PlacedContigReads::const_iterator> > > ERCCI_pcrI_st_st;

    //std::vector<out_order_t>::const_iterator outorderI;
    PlacedContigReads::const_iterator ERCCI_mpcrI;    // main pcrI
    uint32 ERCCI_actcontigpos;

    bool ERCCI_takerails;
    bool ERCCI_takebackbones;

  public:
    ercci_t(Contig * cptr) : ERCCI_contig(cptr), ERCCI_mpcrI(cptr->CON_reads.begin()) {};

    void init(bool takerails,
	      bool takebackbones,
	      uint32 numstrains);
    void update();
    void advance();

    inline const std::vector<std::vector<std::vector<PlacedContigReads::const_iterator>>> & getPCRIstst() const { return ERCCI_pcrI_st_st;}
    inline uint32 getContigPos() const { return ERCCI_actcontigpos;}
  };



  ////////////////////////////////////////////////////////////////////////
  struct nnpos_rep_col_t {
    // the next two are linked, i.e., 1st value of each must
    //  be seen together, 2nd ... etc.
    std::string groupbases;
    std::vector<base_quality_t> groupquals;

    // which urdids are in these groups
    std::vector<uint32> urdids;

    multitag_t::mte_id_t type;       // fille with Read::REA_tagentry_id...
    int32 contigpos;

    bool is_dangerous;
    bool tagged;

    friend std::ostream & operator<<(std::ostream &ostr, const nnpos_rep_col_t & prc){
      ostr << "PRC:\n";
      ostr << "    is_dangerous: " <<  prc.is_dangerous << '\n';
      ostr << "    tagged: " <<  prc.tagged << '\n';
      ostr << "    groupbases: " <<  prc.groupbases << '\n';
      ostr << "    type: " <<  multitag_t::getIdentifierStr(prc.type) << '\n';
      ostr << "    contigpos: " <<  prc.contigpos << '\n';
      return ostr;
    }
  };

public:
  struct repeatmarker_stats_t {
    uint32 numSRMs;
    uint32 numWRMs;
    uint32 numSNPs;

    void init() {numSRMs=0;numWRMs=0;numSNPs=0;}
    repeatmarker_stats_t() {init();}
  };
private:
  ////////////////////////////////////////////////////////////////////////



  void getRailsAsReadsAffected(const int32 refid,
			       std::vector<int32> & reads_affected,
			       const int32 xcut,
			       const int32 ycut);

  bool checkFreeCoverageForAddingRead(
    const uint8 newreadseqtype,
    const int32 xcut,
    const int32 ycut);

  inline uint8 getBBStrainMask(int32 strainid) {
    int32 strainmask=0;
    if(strainid>7 || strainid<0) {
      throw Notify(Notify::FATAL, "uint8 Contig::getBBStrainMask(int32 strainid)", "Mapping >=8 Solexa strains? Sorry, not possible yet.");
    }
    BITSET(strainid, strainmask);
    return static_cast<uint8>(strainmask);
  }
  void getMappedBBStrainIDsFromMask(std::vector<int32> & strains, uint8 mask);


  /****************************************************
   * contig_output.C
   ****************************************************/

  void dumpAsTextOrHTML(std::ostream &ostr,
			const uint8 outtype,
			const std::string & consseq,
			const std::vector<base_quality_t> & consqual,
			const int32 frompos,
			const int32 topos);

  void dumpAsDebug(std::ostream & ostr);
  void dumpAsTCS(std::ostream & ostr);
  void buildFASTAHeader(std::ostringstream & ostr) const;
  void dumpAsFASTA(std::ostream & ostr, bool padded);
  void dumpAsFASTAQUAL(std::ostream & ostr, bool padded);
  void dumpAsCAF(std::ostream & ostr);
  void dumpAsMAF(std::ostream & ostr);
  int32 helper_dumpAsACE_BSLines(std::string & consseq, bool alsodump, std::ostream & ostr);
  void dumpAsACE(std::ostream & ostr);

  void dumpReplay(std::ofstream & eout,
		  const AlignedDualSeqFacts * initialadsf,
		  int32 refid,
		  int32 newid,
		  int32 direction_frnid,
		  bool newid_ismulticopy,
		  int32 forcegrow);

  /****************************************************
   * contig_covanalysis.C
   ****************************************************/

  void findPeaks_helper(ccctype_t avgcov, ccctype_t threshold, std::vector<uint8> & peakindicator);

  /****************************************************
   * contig_featureinfo.C
   ****************************************************/

  void myappend(std::string & a, const std::string & b, const std::string & concatstring) const;
  void calcSOTagCoverage_helper(const contigSOtag_t & csot,
				tagcoverageinfo_t & tci,
				std::vector<uint64> & covvals,
				const tagcoverageinfo_t & ctci);

  /****************************************************
   * contig_analysis.C
   ****************************************************/

  bool analyseRMBZones(PlacedContigReads::const_iterator pcrI);
  void checkForLockedBase(const uint32 readpos,
			  PlacedContigReads::const_iterator pcrI,
			  int32 & baselock,
			  int32 & snplock);
  void checkForLockedBaseComplement(const uint32 readpos,
				    PlacedContigReads::const_iterator pcrI,
				    int32 & baselock,
				    int32 & snplock);

  void calcGroupQual(const nngroups_t & g);
  void nmpr_firstfillin(const ercci_t & ercci,
			const std::vector<int8> & maskshaddow,
			const std::vector<multitag_t::mte_id_t> & masktagtypes,
			std::vector<std::vector<std::vector<nngroups_t> > > & groups_st_st);
  void nmpr_rategroups(std::vector<std::vector<std::vector<nngroups_t> > > & groups_st_st,
		       cccontainer_t::const_iterator ccI);
  void nmpr_secondfillin(const ercci_t & ercci,
			 const std::vector<int8> & maskshaddow,
			 const std::vector<multitag_t::mte_id_t> & masktagtypes,
			 std::vector<std::vector<std::vector<nngroups_t> > > & groups_st_st);
  void nmpr_cautiousMultiSeqTypeTagging(const ercci_t & ercci,
					const rcci_t & rcci,
					const std::vector<std::vector<std::vector<nngroups_t> > > & groups_st_st,
					const nnpos_rep_col_t & emptyprc,
					std::vector<bool> & readsmarkedsrm,
					repeatmarker_stats_t & repstats);
  void nmpr_evaluateOneSeqType(const uint32 actseqtype,
			       const ercci_t & ercci,
  			       const std::vector<std::vector<nngroups_t> > & groups_st,
			       std::vector<nnpos_rep_col_t> & newprc,
			       const nnpos_rep_col_t & emptyprc);
  uint32 nmpr_appendPRCFieldsWithGroupsOfOneStrain(
    const std::vector<nngroups_t> & groups,
    nnpos_rep_col_t & newprc);
  void nmpr_tagColumn(nnpos_rep_col_t & prc,
		      const rcci_t & rcci,
		      std::vector<bool> & readsmarkedsrm,
		      repeatmarker_stats_t & repstats);
  void csbrm_fillin_groups_stst(
    const ercci_t & ercci,
    const std::vector<int8> & maskshadow,
    const std::vector<multitag_t::mte_id_t> & masktagtypes,
    std::vector<std::vector<std::vector<nngroups_t> > > & groups_st_st);
  bool csbrm_checkReadsForHomopolymers(uint32 actcontigpos,
				       uint32 maxrunallowed,
				       const std::vector<int32> & checkthesereads);
  bool csbrm_checkForCleanBlock(cccontainer_t::iterator ccI);
  uint32 getBaseRunLength(const Read & read,
			  const uint32 readpos,
			  const char base,
			  uint32 & from,
			  uint32 & to,
			  uint32 & zeroqualcounts,
			  bool useclipsasends=false);
  void reduceConsensusTags(const std::vector<multitag_t::mte_id_t> & mastertags,
			   const std::vector<multitag_t::mte_id_t> & deletecandidates);


  /****************************************************
   * contig_edit.C
   ****************************************************/
  void priv_rratcp_collectReadsToDelete(ccctype_t minthresh,
					std::vector<uint8> & peakindicator,
					std::vector<ccctype_t> & virtcoverage,
					bool remunpaired,
					bool pairmustbeincontig,
					std::unordered_set<readid_t> & readsremoved);
  bool priv_rratcp_checkDeletable(PlacedContigReads::const_iterator pcrI,
				  ccctype_t minthresh,
				  std::vector<uint8> & peakindicator,
				  std::vector<ccctype_t> & virtcoverage);


  /****************************************************
   * contig_consensus.C
   ****************************************************/

  bool makeTmpConsensus(int32 from, int32 to, bool tmpconsfrombackbone);

  void makeIntelligentConsensus_helper3(
    char & thisbase,
    base_quality_t & thisqual,
    const std::vector<nngroups_t> & groups,
    const std::vector<char> & IUPACbasegroups);

  void makeIntelligentConsensus_helper2_calcSangerQual(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups,
    const std::vector<char> & columnbases,
    const base_quality_t maxqual,
    const uint32 maxcount);

  void makeIntelligentConsensus_helper2_calcSOLEXA(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    cccontainer_t::const_iterator ccI,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups,
    const base_quality_t maxqual,
    const uint32 maxcount,
    int32 strainidtotake);

  void makeIntelligentConsensus_helper2_calc454GS20(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups);

  void makeIntelligentConsensus_helper2_calcText(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups);

  void makeIntelligentConsensus_helper2_calcIonTorrent(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups);

  void makeIntelligentConsensus_helper2_calcPACBIOHQ(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    const std::vector<nngroups_t> & groups,
    std::vector<char> & IUPACbasegroups);

  void makeIntelligentConsensus_helper1(
    char & thisbase,
    base_quality_t & thisqual,
    const uint32 actcontigpos,
    cccontainer_t::const_iterator ccI,
    const int32 mincoverage,
    std::vector<nngroups_t> & groups,
    std::vector<nngroups_t> & maskedshadowgroups,
    std::vector<char> & IUPACbasegroups,
    std::vector<char> & columnbases,
    const std::vector<PlacedContigReads::const_iterator> & read_pcrIs_in_col,
    std::vector<int8> & maskshadow,
    uint8 actreadtype,
    int32 strainidtotake,
    char missingcoveragechar
    );
  void makeIntelligentConsensus(std::string & target,
				std::vector<base_quality_t> & qual,
				std::vector<int32> * targetadjustments,
				std::vector<int32> * targetadjustments_bla,
				int32 from,
				int32 to,
				int32 mincoverage=0,
				base_quality_t minqual=0,
				int32 strainidtotake=-1,
				char missingcoveragechar='@');
//				std::ostream * ostr=nullptr,
//				bool contagsintcs=true);

  int8 rateGoodnessLevelOfConsensus(cccontainer_t::const_iterator ccI,
				    nngroups_t & groups,
				    uint32 numpossiblebases,
				    uint8 seqtype);



/* TODO: reactivate with PlacedContigReads
  int deepThought(const contigread_t &aRead, int offset, char base);
		  //		  SCF_buffer &theBuffer);
  bool analyseRMBZones(const contigread_t & ric, std::vector<int32> & RMBdangeroverlaps);
  bool analyseDangerZones(const contigread_t & ric);
  //bool analyseAllZones(const contigread_t & ric);
*/


/*************************************************************************
 *
 *
 *  PUBLIC
 *
 *
 *************************************************************************/

public:
  Contig(std::vector<MIRAParameters> * params, ReadPool & readpool);
  Contig(Contig const &other);
  ~Contig();

  static void setMasterCEBUGFlag(bool f) {CON_mastercebugflag=f;};

  Contig const & operator=(Contig const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Contig const &con);

  void saveMem();

  const constats_t & getStats() { calcStats(); return CON_stats; };

  void setParams(std::vector<MIRAParameters> * params);

  uint32 getContigID() const {return CON_id;};
  void   setContigID(uint32 id) {CON_id=id;};
  void   setContigNamePrefix(const std::string & pname) {CON_nameprefix=pname;};
  const std::string & getContigName() const;
  void   setContigName(const std::string & name);

  uint32  getNumBackbones() const;

  void trashConsensusCache(bool qualonly) {
    nukeSTLContainer(CON_fixedconsqual);
    if(!qualonly) nukeSTLContainer(CON_fixedconsseq);
  };

  uint32 getNumReadsPerReadGroup(uint32 rgnum) const {
    if(rgnum>=CON_readsperreadgroup.size()) return 0;
    return CON_readsperreadgroup[rgnum];
  }
  uint32 getNumReadsPerStrain(uint32 sid) const {
    if(sid>=CON_readsperstrain.size()) return 0;
    return CON_readsperstrain[sid];
  }

  uint32 getNumReadsInContig() const {return static_cast<uint32>(CON_reads.size());};
  uint32 getContigLength() const {return static_cast<uint32>(CON_counts.size());};

  bool getLongRepeatStatus() const {return CON_contains_long_repeats_only;}
  void setLongRepeatStatus(bool b) {CON_contains_long_repeats_only=b;}
  void calcDigitalNormalisationStatus();

  void getStrainsOfBackbone(std::vector<int32> & sob) const;

  void discard();

  void useBackbone4TmpConsensus(bool b) {CON_tmpcons_from_backbone=b;};
  void mergeNewSRReads(bool b) {CON_mergenewsrreads=b;};
  void setSpecialSRAddConditions(const int32 maxtotalerrors,
				 const int32 maxgaps,
				 const int32 maxmismatches);

  void setContigCoverageTarget(std::vector<uint32>);
  void addRead(
    std::vector<Align> & aligncache,
    const AlignedDualSeqFacts * initialadsf,
    int32 refid,
    int32 newid,
    int32 direction_frnid,
    bool  newid_ismulticopy,
    int32 forcegrow,
    templateguessinfo_t & templateguess,
    errorstatus_t & errstat);
  void addFirstRead(int32 id, int8 direction);
  void coutAddReadTimings();


  void addConsensusSequenceAsReadToContig(int32 strainid);

  PlacedContigReads::const_iterator deleteRead(PlacedContigReads::const_iterator pcrI);
  const std::vector<int32> & getLastDangerzoneOverlapIDS() const {return CON_last_dangerous_overlaps;};

  void initialiseContig(const std::list<contig_init_read_t> & rlist,
			std::vector<multitag_t> & tags,
			const std::string & contigname,
			std::string & fixedseq,
			std::vector<base_quality_t> & fixedqual);

  char getBaseInRead(int32 contigposition,
		     const PlacedContigReads::const_iterator & pcrI) const;
  base_quality_t getQualityInRead(int32 contigposition,
				  const PlacedContigReads::const_iterator & pcrI) const;
  void changeBaseInRead(char base,
			int32 contigposition,
			int32 readindex,
			base_quality_t quality=BQ_UNDEFINED);
  void insertBaseInRead(char base,
			int32 contigposition,
			int32 readindex,
			base_quality_t quality=BQ_UNDEFINED);
  void deleteBaseInRead(int32 contigposition, int32 readindex);

  int32 chompFront(int32 maxchecklen, bool doshiftreads=true);
  int32 chompBack(int32 maxchecklen);
  void trimMapOverhang();

  void adjustReadOffsets(int32 contigposition, int32 offset);

  void addTagToRead(uint32 contigpositionfrom,
		    uint32 contigpositionto,
		    int32 readindex,
		    multitag_t & mt);

  void addTagToConsensus(const uint32 contigpositionfrom,
			 const uint32 contigpositionto,
			 const char strand,
			 const char * identifier,
			 const char * comment,
			 const bool doublecheck,
			 const bool additionalinfo=false,
			 const base_quality_t qualA=0,
			 const base_quality_t qualC=0,
			 const base_quality_t qualG=0,
			 const base_quality_t qualT=0,
			 const base_quality_t qualStar=0);

  bool hasConsensusTag(const multitag_t::mte_id_t identifier) const;

  void sortConsensusTags() {sort(CON_consensus_tags.begin(),CON_consensus_tags.end(),consensustag_t_comparator);};
  const std::vector<consensustag_t> & getConsensusTags() const {return CON_consensus_tags;};
  void clearConsensusTags() {CON_consensus_tags.clear();};
  void deleteTagsInReads(const multitag_t::mte_id_t identifier);

  void updateStatsFromConsensusTags(bool countSRMcs,
				    bool countWRMcs,
				    bool countIUPACs,
				    bool countSTMUs,
				    bool countSTMSs);

  static void setIDCounter(uint32 num) { CON_id_counter=num;}
  const cccontainer_t & getConsensusCounts() const
    { return CON_counts;}
  const PlacedContigReads & getContigReads() const
    { return CON_reads;}

  inline PlacedContigReads::const_iterator getFirstPCRIForReadsCoveringPosition(int32 pos1) const
    {
      auto pcrI(CON_reads.begin());
      if(pos1>=CON_longestreadseen+1) {
	pcrI=CON_reads.getPCRIForReadsStartingAtPos(pos1-CON_longestreadseen-1);
      }
      return pcrI;
    }
  inline PlacedContigReads::const_iterator getFirstPCRIForNonBBReadsCoveringPosition(int32 pos1) const
    {
      auto pcrI(CON_reads.begin());
      if(pos1>=CON_longestnonbbreadseen+1) {
	pcrI=CON_reads.getPCRIForReadsStartingAtPos(pos1-CON_longestnonbbreadseen-1);
      }
      return pcrI;
    }
  inline PlacedContigReads::const_iterator getFirstPCRIForRailsCoveringPosition(int32 pos1) const
    {
      auto pcrI(CON_reads.begin());
      if(pos1>=CON_longestrailseen+1) {
	pcrI=CON_reads.getPCRIForReadsStartingAtPos(pos1-CON_longestrailseen-1);
      }
      return pcrI;
    }

  void getReadORPIDsAtContigPosition(std::vector<int32> & vec,
						int32 pos1,
						int32 pos2=-1) const;
  void getPCRIteratorsAtContigPosition(std::vector<PlacedContigReads::const_iterator> & vec,
				       int32 pos1,
				       int32 pos2=-1) const;
  // TODO: obsolete!
  // This is present in EdIt only now, need to change EdIt
  //
  std::vector<int32> & getCRIndicesAtContigPosition(std::vector<int32> & vec,
                                               int32 pos1,
                                               int32 pos2=-1) const;
  //
  //
  //





  int32 getRealReadPos(const uint32 contigpos, const PlacedContigReads::const_iterator & pcrI) const;

  int32 paddedPos2UnpaddedPos(uint32 padpos);

  void interpolateSRMValuesInCONcounts(cccontainer_t::iterator ccI);
  void initialiseBaseLocks();
  void updateReadTagsFromReadPool();
  void updateBaseLocks(PlacedContigReads::const_iterator pcrI, bool addiftrue);

  void transposeReadSRMTagsToContig();

  void setupAsBackBoneContig();

  size_t addRails(const uint32 raillength,
		const uint32 railoverlap,
		const std::string & straintxt,
		const bool forcestrainset,
		const std::string & railfromstrain,
		const bool simulateonly);

  size_t addSubsequences(const uint32 raillength,
			 const uint32 railoverlap,
			 const std::string & straintxt,
			 const bool forcestrainset,
			 const std::string & railfromstrain,
			 const bool asbackbone,
			 const std::string & bbseq,
			 const std::vector<base_quality_t> & bbqualvec,
			 const std::string & backbonename,
			 const bool initccbbvalues,
			 const bool simulateonly);

  void removeRails();
  void stripToBackbone();
  static void resetCERNumbering() {CON_cer_idcounter=0;};
  void transformCERMappingsToCoverageReads();
  void recalcTemplateIDsAndStrainPresent();

  void exchangeReadIDs(std::vector<uint32> & newids);

  /****************************************************
   * contig_output.C
   ****************************************************/

  static void setCoutType(uint8 type);

  void stats(std::ostream &ostr);
  void dumpAsText(std::ostream &ostr,
		  const int32 frompos = -1,
		  const int32 topos = -1);
  void dumpAsHTML(std::ostream &ostr,
		  const int32 frompos = -1,
		  const int32 topos = -1);

  void dumpConReads();

  void dumpStrainAsFASTAQUAL(std::ostream & fout,
			     std::ostream & qout,
			     uint32 mincoverage,
			     base_quality_t minqual,
			     int32 strainidtotake,
			     bool padded,
			     bool fillholesinstrain);
  void dumpStrainAsFASTQ(std::ostream & fout,
			 uint32 mincoverage,
			 base_quality_t minqual,
			 int32 strainidtotake,
			 bool padded,
			 bool fillholesinstrain);

  void dumpAsSAM(std::ostream & ostr, const SAMCollect & samc, bool alsobackbone);
  void saveAsGAP4DA(const std::string & dirname, std::ostream & fofnstr);

  static void dumpContigReadList_Head(std::ostream &ostr);
  void dumpContigReadList_Body(std::ostream &ostr);
  static void dumpContigStatistics_Head(std::ostream &ostr);
  void dumpContigStatistics_Body(std::ostream &ostr);
  static void dumpReadTagList_Head(std::ostream &ostr);
  void dumpReadTagList_Body(std::ostream &ostr);
  static void dumpConsensusTagList_Head(std::ostream &ostr);
  void dumpConsensusTagList_Body(std::ostream &ostr);
  static void dumpMAF_Head(std::ostream &ostr);
  static void dumpTCS_Head(std::ostream &ostr);
  void dumpTCS_Body(std::ostream &ostr);
  static void dumpWiggle_Head(std::ostream &ostr);
  void dumpWiggle_Body(std::ostream &ostr, std::string & cons, bool aspadded);
  void dumpGCWiggle_Body(std::ostream &ostr, std::string & cons);

  void dumpTagsAsGFF3(std::ostream & ostr, std::string & padcons);
  void dumpTagsAsGFF3(std::ostream & ostr) {std::string dummy; dumpTagsAsGFF3(ostr,dummy);};

  void dumpStatus(std::ostream & ostr);
  void debugDump();

  static void setOutputRails(bool b) {CON_outputrails=b;}


  /****************************************************
   * contig_edit.C
   ****************************************************/
  void reduceReadsAtCoveragePeaks(ccctype_t avgcov, std::vector<uint8> & peakindicator, std::unordered_set<readid_t> & readsremoved);

  bool hasEditableOvercallData() const;
  bool hasSeqTypeData(uint8 seqtype) const;

  uint32 createPacBioDarkStrobeEdits(std::list<pbdse_t> & pbdsev);
  uint32 deleteStarOnlyColumns() { return deleteStarOnlyColumns(0,CON_counts.size(),false,1);}
  uint32 deleteStarOnlyColumns(int32 from,
			     int32 to,
			     bool alsononly=false,
			     uint32 mincov=1);
  uint32 editTrickyOvercalls(const bool onlysetPSHPtags,
				const bool noSRMreadsallowed,
				std::vector<bool> & readsmarkedsrm);
  uint32 edit454checkTrickies(const char basehypo,
			      const uint32 actcontigpos,
			      const std::vector<PlacedContigReads::const_iterator> & pcrIs_in_col,
			      std::vector<edit454command_t> & editcommands,
			      const bool onlysetPSHPtags,
			      const bool noSRMreadsallowed);
  uint32 editSingleDiscrepancyNoHAFTag(std::vector<bool> & readsmarkedsrm, uint8 editmode);
  void editPBSledgeHammer(std::vector<bool> & readsmarkedsrm,
			  uint32 & numcoledits,
			  uint32 & numreadedits);

  void blindContig();
  void upDownCase(base_quality_t threshold);


  /****************************************************
   * contig_featureinfo.C
   ****************************************************/

  void getGBFSummary(std::list<gbfsummary_t> & allGBfeatures,
		     const std::vector<multitag_t::mte_id_t> & allowedfeatures,
		     const std::vector<multitag_t::mte_id_t> & forbiddenfeatures,
		     bool simulateintergenics,
		     bool alsomiraconstags) const;
  void concatAllGBFInfoForLocus(const std::list<gbfsummary_t> & allGBfeatures,
				std::list<gbfsummary_t>::const_iterator gbfsI,
				const std::string & concatstring,
				std::string & gene,
				std::string & function,
				std::string & ecnumber,
				std::string & product,
				std::string & note) const;

  void getSeqOntTags(std::list<contigSOtag_t> & allSOfeatures);
  void calcSOTagCoverage(const std::list<contigSOtag_t> & features, std::list<tagcoverageinfo_t> & result);


  /****************************************************
   * contig_analysis.C
   ****************************************************/

  std::pair<int32,int32> findBestNonMisassembledRange();
  void trimContigToRange(uint32 from, uint32 to);
  bool tctr_helper1(uint32 from, uint32 to, bool simulateonly);
  void tctr_helper2();
  void tctr_helper3();

  void allowedRefIDs_forbidAll();
  void allowedRefIDs_allowAll();
  void allowedRefIDs_allowSomeRailsOnShortReadCoverage(int32 belowcoverage,
						       int32 strainid,
						       uint8 seqtype);
  const std::vector<bool> & getAllowedRefIDs() const {return CON_allowedrefids;};
  void newMarkPossibleRepeats(repeatmarker_stats_t & repstats,
			      std::vector<bool> & readsmarkedsrm);
  void markFeaturesByConsensus(bool markSROs, bool markUNSs, bool markIUPACs);
  bool analyseReadCoverage(std::vector<uint32> & maxcovperread,
			   std::vector<uint8> & multicopies,
			   std::vector<uint32> covperst);    // no reference here!
  void largeGapRepeatMarker(uint32 mingaplen,
			    uint32 mingroupsize,
			    bool needbothstrands);
  void codonSingleBaseRepeatMarker(uint32 mingroupsize,
				   repeatmarker_stats_t & repstats,
				   std::vector<bool> & readsmarkedsrm);
  bool shouldHaveTPartnerInContig(PlacedContigReads::const_iterator pcrI,
				  PlacedContigReads::const_iterator opcrI);


  /****************************************************
   * contig_pairconsistency.C
   ****************************************************/
  std::pair<int32,int32> findBestPairConsistencyRange();
  void priv_cprForRGID(ReadGroupLib::ReadGroupID rgid,
		       std::vector<uint64> & spanner,
		       std::vector<uint64> & spanner12,
		       std::vector<uint64> & nonspanner);

  /****************************************************
   * contig_covanalysis.C
   ****************************************************/

  void collectCoverage(uint32 from, uint32 to, std::vector<uint64> & covvals);
  void collectCoverage(std::vector<uint64> & covvals) {collectCoverage(0,CON_counts.size()-1,covvals); }
  void calcStatsOnContainer(coverageinfo_t & tci, std::vector<uint64> & covvals) const;
  void calcSecondOrderStatsOnContainer(coverageinfo_t & tci, const std::vector<uint64> & covvals) const;
  void findPeaks(ccctype_t avgcov, std::vector<uint8> & peakindicator);


  /****************************************************
   * contig_consensus.C
   ****************************************************/

  void calcConsensi(int32 mincoverage=0,
		    base_quality_t minqual=0,
		    char missingcoveragechar='X');
  void newConsensusGet(std::string & target,
		       std::vector<base_quality_t> & qual,
		       int32 strainidtotake=-1);
  void ensureConsensus(int32 strainidtotake=-1);

  void updateBackboneConsensus();


};


#endif
