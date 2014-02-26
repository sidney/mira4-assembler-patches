// to be included by assembly, pathfinder ...

#ifndef _mira_overlapedges_h_
#define _mira_overlapedges_h_

#include "stdinc/defines.H"
#include <ostream>

struct newedges_t{
  int32 rid1;                   // number of read in seqpool
  int32 linked_with;                   // number of read in seqpool
  uint32 best_weight;                   // score ratio^2*score. Make sure it is >0 !!! (for pathfinder)
  uint32 adsfindex;              // index of adsfacts elements this refers to
  int16 direction;

  bool pf_allowquickoverlap:1;   // if AS_allowquickoverlap for bot rids is true

  bool pf_banned:1;              // temp use by pathfinder: banned overlap

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1; // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >=5

  friend std::ostream & operator<<(std::ostream &ostr, const newedges_t & e){
    ostr << "NE:\t" << e.rid1
	 << '\t' << e.linked_with
	 << '\t' << e.best_weight
	 << '\t' << e.adsfindex
	 << "\tdir " << e.direction
	 << "\tban " << e.pf_banned
	 << "\tsg  " << e.ol_stronggood
	 << "\twg  " << e.ol_weakgood
	 << "\tbaf " << e.ol_belowavgfreq
	 << "\tnrp " << e.ol_norept
	 << "\trep " << e.ol_rept
	 << '\n';
    return ostr;
  }
};


struct skimedges_t{
  int32  rid1;                   // number of read in seqpool
  int32  linked_with;                   // number of read in seqpool
  int32  eoffset;

  uint32 skimweight;                   // score ratio^2*score
  size_t skimindex;         // index of line in original skim file

  int8   rid1dir;  // new for reduceSkimHits2
  int8   rid2dir;
  uint8  scoreratio;

  bool ol_stronggood:1;  // frequency: 2*bph-1 pos at 3, thereof bph-1 contiguous
  bool ol_weakgood:1;   // frequency: bph-1 positions contiguous at 3
  bool ol_belowavgfreq:1; // frequency: bph-1 positions contiguous at <=3
  bool ol_norept:1;      // nothing >3 (but can contain 1 (single hashes == errors)
  bool ol_rept:1;      // bph-1 positions >= 5 freq

  friend std::ostream & operator<<(std::ostream &ostr, const skimedges_t & e){
    ostr << "SE:\t" << e.rid1
	 << '\t' << e.linked_with
	 << '\t' << e.eoffset
	 << '\t' << e.skimweight
	 << '\t' << e.skimindex
	 << "\td1 " << static_cast<int16>(e.rid1dir)
	 << "\td2 " << static_cast<int16>(e.rid2dir)
	 << "\tsr " << static_cast<uint16>(e.scoreratio)
	 << "\tsg  " << e.ol_stronggood
	 << "\twg  " << e.ol_weakgood
	 << "\tbaf " << e.ol_belowavgfreq
	 << "\tnrp " << e.ol_norept
	 << "\trep " << e.ol_rept
	 << '\n';
    return ostr;
  }
};


#endif
