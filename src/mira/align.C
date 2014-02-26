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


#include "mira/align.H"

#include "errorhandling/errorhandling.H"

using namespace std;


#ifdef CEBUGFLAG
#define CEBUG(bla)   {cout << bla; cout.flush();}
#define CEBUGF(bla)  {cout << bla; cout.flush();}
#else
#define CEBUG(bla)
#define CEBUGF(bla)
#endif

#define CLOCK_STEPS1


#if __GNUC__ >= 3
#define prefetchrs(p)     __builtin_prefetch((p), 0, 0)
#define prefetchws(p)     __builtin_prefetch((p), 1, 0)
#define prefetchrl(p)     __builtin_prefetch((p), 0, 3)
#define prefetchwl(p)     __builtin_prefetch((p), 1, 3)
#else
#define prefetchrs(p)
#define prefetchws(p)
#define prefetchrl(p)
#define prefetchwl(p)
#endif


uint64 Align::AL_alloccount=0;



/*************************************************************************
 *
 *
 *
 *************************************************************************/

Align::Align(MIRAParameters * params): Dynamic(params)
{
  FUNCSTART("Align::Align(MIRAParameters * params): Dynamic(params)");

  AL_miraparams=params;

  AL_tmpads=nullptr;
  AL_alseq1=nullptr;
  AL_alseq2=nullptr;
  AL_as12size=0;

  init();
  resetTimings();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::init()
{
  FUNCSTART("void Align::init()");

  AL_no_solutions=0;
  AL_no_diff_solutions=0;
  AL_max_relscore=0;

  AL_mpcache_dyn_score_multiplier=0;
  AL_mpcache_dyn_score_gap=0;
  AL_mpcache_al_max_cutoff=0;
  AL_mpcache_al_min_score=0;
  AL_mpcache_al_min_overlap=0;
  AL_mpcache_al_min_relscore=0;

  AL_mpset_al_min_overlap=0;
  AL_mpset_al_min_relscore=0;

  AL_valid=131;

  FUNCEND();
  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

Align::~Align()
{
  FUNCSTART("Align::~Align()");

  if(AL_tmpads != nullptr) delete AL_tmpads;
  if(AL_alseq1 != nullptr) delete [] AL_alseq1;
  if(AL_alseq2 != nullptr) delete [] AL_alseq2;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::acquireSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, int32 id1, int32 id2, int8 id1dir, int8 id2dir, bool calcwithoffset, int32 expectedoffset)
{
  FUNCSTART("Align::acquireSequences(const char * seq1, uint32 len1, const char * seq2, uint32 len2, int32 id1, int32 id2, int8 id1dir, int8 id2dir, bool calcwithoffset, int32 expectedoffset)");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  Dynamic::setSequences(seq1, len1, seq2, len2, calcwithoffset, expectedoffset);

  AL_id1=id1;
  AL_id2=id2;
  AL_id1dir=id1dir;
  AL_id2dir=id2dir;

  init();

#ifdef CLOCK_STEPS1
  AL_timing_acquires+=diffsuseconds(tv);
#endif

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::discard()
{
  FUNCSTART("Align::discard()");

  if(AL_valid != 131) {
    cerr << "AL_valid not valid?!\n";
    exit(0);
  }

  Dynamic::discard();

  if(AL_alseq1!=nullptr) delete [] AL_alseq1;
  if(AL_alseq2!=nullptr) delete [] AL_alseq2;
  AL_alseq1=nullptr;
  AL_alseq2=nullptr;
  AL_as12size=0;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::prepareAlign(list<AlignedDualSeq> * adslist)
{
  FUNCSTART("void Align::prepareAlign(list<AlignedDualSeq> * adslist)");

  if(!DYN_valid){
    MIRANOTIFY(Notify::INTERNAL, ": Programming error. Tried to align without proper initialisation of Align-object.");
  }

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  AL_adslist=adslist;
  if(AL_tmpads==nullptr) AL_tmpads= new AlignedDualSeq(AL_miraparams);

  AL_align_maxlen=DYN_len_seq1+DYN_len_seq2+1;

  if(AL_as12size < AL_align_maxlen){
    if(AL_alseq1 != nullptr) delete [] AL_alseq1;
    if(AL_alseq2 != nullptr) delete [] AL_alseq2;

    AL_as12size=AL_align_maxlen;
    if(AL_as12size<2000) AL_as12size=2000;

    AL_alseq1= new char[AL_as12size];
    AL_alseq2= new char[AL_as12size];
    AL_alloccount+=2;
  }

  AL_no_solutions=0;
  AL_no_diff_solutions=0;
  AL_max_relscore=0;

  AL_error_hit_band=false;

  // clear the aligned sequences
  // TODO: check: do I need this???
  memset(AL_alseq1,0,AL_align_maxlen);
  memset(AL_alseq2,0,AL_align_maxlen);

  AL_allen=AL_align_maxlen-1;

  AL_seq1ptr=DYN_sequence1+DYN_len_seq1-1;
  AL_seq2ptr=DYN_sequence2+DYN_len_seq2-1;


#ifdef CLOCK_STEPS1
  AL_timing_prepalign+=diffsuseconds(tv);
#endif

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::setRAlignParams()
{
  align_parameters const & AL_params = AL_miraparams->getAlignParams();
  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

  if(AL_mpset_al_min_relscore!=0){
    AL_mpcache_al_min_relscore=AL_mpset_al_min_relscore;
  }else{
    AL_mpcache_al_min_relscore=AL_params.al_min_relscore;
  }
  if(AL_mpset_al_min_overlap!=0){
    AL_mpcache_al_min_overlap=AL_mpset_al_min_overlap;
  }else{
    AL_mpcache_al_min_overlap=AL_params.al_min_overlap;
  }

  AL_mpcache_al_min_score=AL_params.al_min_score;
  AL_mpcache_al_max_cutoff=AL_params.al_max_cutoff;

  AL_mpcache_dyn_score_gap=DYN_params.dyn_score_gap;
  AL_mpcache_dyn_score_multiplier=DYN_params.dyn_score_multiplier;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::simpleAlign(list<AlignedDualSeq> * adslist, bool enforce_clean_ends, bool dontpenalisengaps)
{
  FUNCSTART("void Align::simpleAlign(list<AlignedDualSeq> * adslist, bool enforce_clean_ends, bool dontpenalisengaps)");

  AL_enforce_clean_ends=enforce_clean_ends;
  AL_dont_penalise_ngaps=dontpenalisengaps;

  prepareAlign(adslist);

  AL_new_solution=1;
  AL_cutoff_counter=0;

  setRAlignParams();

  rAlign(DYN_len_seq1, DYN_len_seq2,'d',false);

  //delete AL_tmpads;

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Align::fullAlign(list<AlignedDualSeq> * adslist, bool enforce_clean_ends, bool dontpenalisengaps)
{
  FUNCSTART("void Align::fullAlign(list<AlignedDualSeq> * adslist, bool enforce_clean_ends, bool dontpenalisengaps)");

#ifdef CLOCK_STEPS1
  timeval tv;
  gettimeofday(&tv,nullptr);
#endif

  AL_enforce_clean_ends=enforce_clean_ends;
  AL_dont_penalise_ngaps=dontpenalisengaps;

  prepareAlign(adslist);

  setRAlignParams();

  termAlign();

#ifdef CLOCK_STEPS1
  AL_timing_fullalign+=diffsuseconds(tv);
#endif

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define AL_TATRACE
//#define CEBUG(bla)   {cout << bla; cout.flush();}

//#define CEBUG(bla)   {if(AL_id1==-1 && AL_id2==14200) cout << bla; cout.flush();}


// looks like multiplier of 10000 rather than 100 is more sure

void Align::termAlign()
{
  FUNCSTART("void Align::termAlign()");

  align_parameters const & AL_params = AL_miraparams->getAlignParams();

  // FIXME: bad fix. if both are 1, termAlign will read across memborders
  // Albeit this should never happen, I'll fix termAlign later
  if(AL_params.al_min_score==1 && AL_params.al_min_overlap == 1){
    align_parameters & AL_ncparams = const_cast<align_parameters &>(AL_miraparams->getAlignParams());
    AL_ncparams.al_min_score=2;
    AL_ncparams.al_min_overlap=2;
  }

  // simple tactics now: just take the best value from last column
  //  or row

  int32 maxvalrow=-1;
  int32 maxvalcol=maxvalrow;

  uint32 maxvalrowpos=0;
  uint32 maxvalcolpos=0;
  {
    CEBUG("\n\nLast row (left to right):\n");
    for(uint32 i=1; i<=DYN_len_seq2; ++i) {
      int32 vAL_actpos=DYN_simmatrix[((DYN_len_seq1)*(DYN_len_seq2+1)+i)];
      CEBUG("R-Looking: " << i << "\tval:" << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/min(i,DYN_len_seq1);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos>=maxvalrow){
	maxvalrow=vAL_actpos;
	maxvalrowpos=i;
      }
    }

    CEBUG("Last column:\n");
    for(uint32 i=1; i<=DYN_len_seq1; ++i) {
      int32 vAL_actpos=(DYN_simmatrix[i*(DYN_len_seq2+1)+DYN_len_seq2]);
      CEBUG("C-Looking: " << i << "\tval: " << vAL_actpos << '\n');
      int32 relscore=(vAL_actpos*10000)/min(i,DYN_len_seq2);
      if(relscore >= AL_params.al_min_relscore
	 && vAL_actpos>=AL_params.al_min_score
	 && vAL_actpos >= maxvalcol){
	maxvalcol=vAL_actpos;
	maxvalcolpos=i;
      }
    }
  }

  CEBUG("Maxvalrow: " << maxvalrow << "\tpos" << maxvalrowpos << "\n");
  CEBUG("Maxvalcol: " << maxvalcol << "\tpos" << maxvalcolpos << "\n");

  AL_new_solution=1;
  AL_cutoff_counter=0;
  //if(maxvalrow==maxvalcol
  //   && maxvalrowpos+1==DYN_len_seq2
  //   && maxvalcolpos+1==DYN_len_seq1){
  //  // we start in the corner
  //  rAlign(DYN_len_seq1, DYN_len_seq2,'d');
  //}else
  if(maxvalrow>maxvalcol){
    if(maxvalrowpos>0) {
      CEBUG("maxvalrowpos: " << maxvalrowpos << '\t' << "DYN_len_seq2: " << DYN_len_seq2 << endl);
      for(uint32 i=maxvalrowpos; i<DYN_len_seq2; ++i){
	AL_alseq1[--AL_allen]=' ';
	AL_alseq2[AL_allen]=*AL_seq2ptr--;
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(DYN_len_seq1, maxvalrowpos,'d',false);
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  } else {
    if(maxvalcolpos>0) {
      CEBUG("maxvalcolpos: " << maxvalcolpos << '\t' << "DYN_len_seq1: " << DYN_len_seq1 << endl);
      for(uint32 i=maxvalcolpos; i<DYN_len_seq1; ++i){
	AL_alseq1[--AL_allen]=*AL_seq1ptr--;
	AL_alseq2[AL_allen]=' ';
      }
      CEBUG("DoAlign");
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      rAlign(maxvalcolpos, DYN_len_seq2,'d',false);
#ifdef CLOCK_STEPS1
      AL_timing_raligntot+=diffsuseconds(tv);
#endif
    }
  }
  CEBUG("I'm out! AL_adslist->size(): " << AL_adslist->size() << endl;);

  FUNCEND();
}

//#define CEBUG(bla)





/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define FUNCSTART(bla)  static const char * THISFUNC = bla"  ";
//#define FUNCTRACE(bla) { cout << THISFUNC << bla; cout.flush();}
//#define FUNCEND()

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void Align::rAlign(uint32 i, uint32 j, char lastdir, bool hadn)
{
  FUNCSTART("void Align::rAlign(uint32 i, uint32 j, char lastdir)");

//  align_parameters const & AL_params = AL_miraparams->getAlignParams();
//  dynamic_parameters const & DYN_params = DYN_miraparams->getDynamicParams();

  uint32 mll=DYN_len_seq2+1;

  //  prefetchrl(&DYN_match_matrix[static_cast<uint8>(*AL_seq1ptr)][static_cast<uint8>(*AL_seq2ptr)]);
  //  prefetchrl(&DYN_simmatrix[i*mll+j-1]);
  //  prefetchrl(&DYN_simmatrix[i*mll+j]);
  //  prefetchrl(&DYN_simmatrix[(i-1)*mll+j-1]);
  //  prefetchrl(&DYN_simmatrix[(i-1)*mll+j]);

  CEBUG("Dong!\n");
  CEBUG("i: " << i << "\tj: " << j << "\tAL_allen: " << AL_allen<< endl);
  CEBUG("s[i,j]:" << DYN_simmatrix[i*(DYN_len_seq2+1)+j]<< endl);
  if(i>0 && j>0){
    CEBUG("s[i-1,j-1]:" << DYN_simmatrix[((i-1)*(DYN_len_seq2+1))+j-1] << endl);
  }

  bool hasn=false;
  if(AL_seq1ptr>=DYN_sequence1 && *AL_seq1ptr){
    CEBUG("seq1ptr points on: " << *AL_seq1ptr << endl);
    hasn=(*AL_seq1ptr=='N');
  }else{
    CEBUG("seq1ptr points on: nullptr\n");
  }
  if(AL_seq2ptr>=DYN_sequence2 && *AL_seq2ptr){
    CEBUG("seq2ptr points on: " << *AL_seq2ptr << endl);
    hasn=hasn | (*AL_seq2ptr=='N');
  }else{
    CEBUG("seq2ptr points on: nullptr\n");
  }
  if(unlikely(AL_allen>AL_align_maxlen)) {
    cerr << "allen:: "<< AL_allen;
    MIRANOTIFY(Notify::INTERNAL, ": FOOOOOO!.") ;
  }

  if(unlikely(AL_cutoff_counter==AL_mpcache_al_max_cutoff)) {
    CEBUG("back... (because of cutoff)\n");
    return;
  }
  if(unlikely(AL_error_hit_band)) {
    CEBUG("back... (because of band hit)\n");
    return;
  }

  if(unlikely(i==0 && j==0)) {
    if(AL_seq1ptr==DYN_sequence1-1 && AL_seq2ptr==DYN_sequence2-1){
#ifdef CLOCK_STEPS1
      timeval tv;
      gettimeofday(&tv,nullptr);
#endif
      AL_tmpads->acquireSequences(AL_alseq1+AL_allen, AL_alseq2+AL_allen, AL_id1, AL_id2, AL_id1dir, AL_id2dir, AL_enforce_clean_ends, AL_dont_penalise_ngaps);
#ifdef CLOCK_STEPS1
      AL_timing_ra_adsacquire+=diffsuseconds(tv);
#endif
      CEBUG("Solution\n");
      CEBUG(*AL_tmpads);

      if(AL_tmpads->getScore() >= static_cast<int32>(AL_mpcache_al_min_score / AL_mpcache_dyn_score_multiplier)
	 && AL_tmpads->getOverlapLen() >= static_cast<uint32>(AL_mpcache_al_min_overlap)
	 && AL_tmpads->getScoreRatio() >= static_cast<int32>(AL_mpcache_al_min_relscore)){
#ifdef CLOCK_STEPS1
	gettimeofday(&tv,nullptr);
#endif
	AL_adslist->push_back(*AL_tmpads);
#ifdef CLOCK_STEPS1
	AL_timing_ra_adslist+=diffsuseconds(tv);
#endif
      }

      ++AL_cutoff_counter;
      ++AL_no_solutions;
      if(AL_new_solution){
	AL_new_solution=0;
	++AL_no_diff_solutions;
      }
    }
    CEBUG("back...\n");
    return;
  }

  if(i==0){

    CEBUG("i=0 left ...\n");

    AL_alseq1[--AL_allen]=' ';
    AL_alseq2[AL_allen]=*AL_seq2ptr--;
    rAlign(i,j-1,'l',hasn);
//    AL_alseq1[AL_allen]='^';
//    AL_alseq2[AL_allen++]='^';
    AL_allen++;
    AL_seq2ptr++;
  }else if(j==0){

    CEBUG("j=0 up ...\n");

    AL_alseq1[--AL_allen]=*AL_seq1ptr--;
    AL_alseq2[AL_allen]=' ';
    rAlign(i-1,j,'u',hasn);
//    AL_alseq1[AL_allen]='%';
//    AL_alseq2[AL_allen++]='%';
    AL_allen++;
    AL_seq1ptr++;
  } else if(j-(DYN_leftbandx+i) < 5              // bugfix parentheses; 5 as a gap-bridge while waiting for new code ft_pcbiolow / rle
	    || DYN_rightbandx+i-j < 5){
    AL_error_hit_band=true;
    CEBUG("HIT BAND LIMIT: " << j-(DYN_leftbandx+i) << " " << DYN_rightbandx+i-j << endl);
  }else{
    int32 vgl=DYN_match_matrix[static_cast<uint8>(*AL_seq1ptr)][static_cast<uint8>(*AL_seq2ptr)];

    CEBUG("vgl: " << vgl << endl);

    // precompute possibilities, optimised for memory access
    // at the same time, check whether we hit a band limit
    if(DYN_simmatrix[i*mll+j-1]==DYN_BANDLIMIT) AL_error_hit_band=true;
    bool leftok=(DYN_simmatrix[i*mll+j-1]+AL_mpcache_dyn_score_gap==DYN_simmatrix[i*mll+j]);
    bool diagok=(DYN_simmatrix[(i-1)*mll+j-1]+vgl==DYN_simmatrix[i*mll+j]);
    if(DYN_simmatrix[(i-1)*mll+j]==DYN_BANDLIMIT) AL_error_hit_band=true;
    bool upok=(DYN_simmatrix[(i-1)*mll+j]+AL_mpcache_dyn_score_gap==DYN_simmatrix[i*mll+j]);

    // This is a special rule for strobe sequencing
    // If a stretch of 'N' occurs and there's the choice between
    //  a match and a gap, choose gap
    // this prevents things like:
    //    ..nnnnnnnaaaaaaaaat.....
    //    ..nnnnnnn*****nnnnt.....
    // and make it
    //    ..nnnnnnnaaaaaaaaat.....
    //    ..nn*****nnnnnnnnnt.....
    //
    // i.e., align gaps to N

    // BaCh 21.04.2010: not that good this idea.
    // - new PacBio dark strobe editing strategy doesn't need this
    // - no directly visible improvement
    // -> remove.
    //if(*AL_seq1ptr=='N'
    //   && *AL_seq2ptr=='N'
    //   && diagok
    //   && (leftok || upok)){
    //  diagok=false;
    //}

    // prevent
    //    ..cccctccaccgaatgcctaa
    //    ..cccctcc****a-----------------
    // and make it
    //    ..cccctccaccgaatgcctaa
    //    ..cccctcca****-----------------

    if(hadn
       && diagok
       && (leftok || upok)){
      diagok=false;
    }


    // now, look which possibilities we have
    //  if a diagonal and (up or left) are equal, then take
    //  the same direction as last time
    // this prevents things like:
    //    ..tggaaaaaaaaat.....
    //    ..t*g*********t.....
    // and make it
    //    ..tggaaaaaaaaat.....
    //    ..tg**********t.....


    //if(diagok){
      if(upok){
	if(lastdir=='u') {
	  diagok=false;
	}
      }
      if(leftok){
	if(lastdir=='l') {
	  diagok=false;
	  upok=false;
	}
      }
    //}

    if(diagok){

      CEBUG("diagonal ...\n");

      AL_alseq1[--AL_allen]=*AL_seq1ptr--;
      AL_alseq2[AL_allen]=*AL_seq2ptr--;
      rAlign(i-1,j-1,'d',hasn);
//	AL_alseq1[AL_allen]='/';
//	AL_alseq2[AL_allen++]='/';
      AL_allen++;
      AL_seq1ptr++;
      AL_seq2ptr++;
    }

    if(upok && !AL_error_hit_band){

      CEBUG("up ...\n");

      AL_alseq1[--AL_allen]=*AL_seq1ptr--;
      AL_alseq2[AL_allen]='*';
      rAlign(i-1,j,'u',hasn);
//	AL_alseq1[AL_allen]='$';
//	AL_alseq2[AL_allen++]='$';
      AL_allen++;
      AL_seq1ptr++;
    }

    if(leftok && !AL_error_hit_band){

      CEBUG("left ...\n");

      AL_alseq1[--AL_allen]='*';
      AL_alseq2[AL_allen]=*AL_seq2ptr--;
      rAlign(i,j-1,'l',hasn);
//	AL_alseq1[AL_allen]='#';
//	AL_alseq2[AL_allen++]='#';
      AL_allen++;
      AL_seq2ptr++;
    }
  }

  CEBUG("back...\n");


  FUNCEND();
}



void Align::coutWhatWasGiven()
{
  Dynamic::coutWhatWasGiven();

  cout << "Align\n------\n";

}


void Align::resetTimings()
{
  DYN_timing_seqcopy=0;
  DYN_timing_bswmatrix=0;
  DYN_timing_bswm_setup=0;
  DYN_timing_bswm_p1=0;
  DYN_timing_bswm_p2a=0;
  DYN_timing_bswm_p2b=0;
  DYN_timing_bswm_p3=0;
  DYN_timing_bswm_cleanband=0;
  AL_timing_acquires=0;
  AL_timing_fullalign=0;
  AL_timing_prepalign=0;
  AL_timing_raligntot=0;
  AL_timing_ra_adsacquire=0;
  AL_timing_ra_adslist=0;
}

void Align::dumpTimings()
{
  cout << "Align timing DYN seqcpy : " << DYN_timing_seqcopy << endl;
  cout << "Align timing DYN bsw su : " << DYN_timing_bswm_setup << endl;
  cout << "Align timing DYN bsw p1 : " << DYN_timing_bswm_p1 << endl;
  cout << "Align timing DYN bsw p2a: " << DYN_timing_bswm_p2a << endl;
  cout << "Align timing DYN bsw p2b: " << DYN_timing_bswm_p2b << endl;
  cout << "Align timing DYN bsw p3 : " << DYN_timing_bswm_p3 << endl;
  cout << "Align timing DYN bsw cb : " << DYN_timing_bswm_cleanband << endl;

  cout << "Align timing DYN bsw    : " << DYN_timing_bswmatrix << endl;
  cout << "Align timing AL acqu s  : " << AL_timing_acquires << endl;
  cout << "Align timing AL full    : " << AL_timing_fullalign << endl;
  cout << "Align timing AL prep    : " << AL_timing_prepalign << endl;
  cout << "Align timing AL ralignt : " << AL_timing_raligntot << endl;
  cout << "Align timing AL ralignc : " << AL_timing_raligntot-AL_timing_ra_adsacquire-AL_timing_ra_adslist << endl;
  cout << "Align timing AL ads a   : " << AL_timing_ra_adsacquire << endl;
  cout << "Align timing AL ads s   : " << AL_timing_ra_adslist << endl;
}
