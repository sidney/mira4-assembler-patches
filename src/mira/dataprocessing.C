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


#include "mira/dataprocessing.H"

#include <boost/unordered_map.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


#include "mira/assembly.H"
#include "mira/hashstats.H"


using namespace std;


//#define CEBUG(bla)   {if(CEBUGFLAG) {cout << bla; cout.flush();}}
#define CEBUG(bla)



std::vector<DataProcessing::poolskim_t> DataProcessing::DP_adapskims;
boost::mutex DataProcessing::DP_ps_changemutex;

HashStatistics DataProcessing::DP_phix174hashstatistics;
bool DataProcessing::DP_px174hs_init=false;
boost::mutex DataProcessing::DP_px174hs_changemutex;         // exclusive mutex for write access to hashstatistics, we need that as it's a static variable



bool DataProcessing::priv_staticInitialiser()
{
  DP_freqnames.push_back(Read::REA_tagentry_idHAF0);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF1);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF2);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF3);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF4);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF5);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF6);
  DP_freqnames.push_back(Read::REA_tagentry_idHAF7);
  return true;
}


DataProcessing::DataProcessing(std::vector<MIRAParameters> * params) : DP_miraparams_ptr(params), DP_tmpmtpolyAT(Read::REA_defaulttag_SOFApolyA_sequence)
{
  DP_threadid=-1;
  DP_tmpvu8.reserve(16300); // bit less than 16kb

  // vector with enough capacity so that it does not get reallocated
  // -> multiple threads won't get their data removed under them during the run
  DP_adapskims.reserve(1024);

};



DataProcessing::~DataProcessing()
{
  stopLogging();

//  for(auto & ps : DP_adapskims){
//    if(ps.poolptr!=nullptr) {
//      delete ps.poolptr;
//      delete ps.skimptr;
//    }
//  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics & hsd, vector<uint8> * debrisreasonptr)
{
  FUNCSTART("void DataProcessing::performDigitalNormalisation_Pool(ReadPool & rp, HashStatistics & hsd)");

  uint32 numtaken=0;
  uint32 numnormout=0;

  vector<bool> normdone(rp.size(),false);
  vector<bool> normout(rp.size(),false);
  vector<bool> normthisrg(rp.size(),false);

  // do the normalisation for every readgroup so that we independently get reads from every rg
  for(auto rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);

    hsd.digiNormReset();
    normthisrg.clear();
    normthisrg.resize(rp.size(),false);
    cout << "\nReadgroup " << rgi << ":\n";
    ProgressIndicator<int32>  pi(0,rp.size()*2);

    // step 0 : only reads with freq >= 2, fwd, rev and no N
    // step 1 : remaining reads
    for(uint32 step=0; step<2; ++step){
      for(uint32 rpi=0; rpi<rp.size(); ++rpi){
	pi.increaseprogress();
	if(normdone[rpi]) continue;
	auto & actread=rp[rpi];
	//Read::setCoutType(Read::AS_TEXT);
	//cout << "### bla\n";
	//cout << actread << endl;
	bool lookatread=true;
	if(step==0){
	  auto bhsI=actread.getBPosHashStats().begin();
	  bhsI+=actread.getLeftClipoff();
	  for(auto ri=0; ri<actread.getLenClippedSeq(); ++ri){
	    if(bhsI->fwd.getFrequency()<2 || !bhsI->fwd.hasConfirmedFwdRev()) lookatread=false;
	  }
	  // getClippedSeqAsChar may throw if the MIRA clipping set the left cutoff to the length of the sequence
	  // too lazy to get things otherwise;
	  if(actread.getLenClippedSeq()>0){
	    auto sptr=actread.getClippedSeqAsChar();
	    auto eptr=sptr+actread.getLenClippedSeq();
	    for(; sptr!=eptr; ++sptr){
	      if(toupper(*sptr)=='N') {
		lookatread=false;
		break;
	      }
	    }
	  }
	}else if(step==1){
	  // further tests?
	  lookatread=true;
	}else{
	  BUGIFTHROW(true,"Oooops, step " << step << " not foreseen?");
	}
	if(lookatread){
	  normthisrg[rpi]=true;
	  normdone[rpi]=true;
	  auto taken=hsd.digiNormTestRead(actread,false);
	  if(taken){
	    ++numtaken;
	    CEBUG("Kept " << actread.getName() << endl);
	    auto & tv = const_cast<std::vector<multitag_t> &>(actread.getTags());
	    bool hasolddgnr=false;
	    for(auto & te : tv){
	      if(te.identifier == Read::REA_tagentry_idDGNr
		 && !te.getCommentStr().empty()){
		hasolddgnr=true;
		break;
	      }
	    }
	    if(!hasolddgnr){
	      for(auto & te : tv){
		if(te.identifier == Read::REA_tagentry_idMNRr){
		  bool founddouble=false;
		  for(auto & tf : tv){
		    if(tf.identifier == Read::REA_tagentry_idDGNr
		      && tf.from==te.from
		       && tf.to==te.to){
		      founddouble=true;
		      break;
		    }
		  }
		  if(!founddouble) te.identifier = Read::REA_tagentry_idDGNr;
		}
	      }
	    }
	    actread.deleteTag(Read::REA_tagentry_idMNRr);
	  }else{
	    ++numnormout;
	    CEBUG("NormOut " << actread.getName() << endl);
	    normout[rpi]=true;
	    actread.setRQClipoff(0);
	  }
	}
      }
    }
    pi.finishAtOnce();

    cout << "Calculating replacement coverage";
    uint32 chkall=0;
    for(uint32 rpi=0; rpi<rp.size(); ++rpi){
      if(normthisrg[rpi]){
	auto & actread=rp[rpi];
	auto & tv = const_cast<std::vector<multitag_t> &>(actread.getTags());
	for(auto & te : tv){
	  if(te.identifier == Read::REA_tagentry_idDGNr){
	    CEBUG(actread.getName() << "\tDGNr: " << te.to-te.from+1 << "\t" << actread.getLenClippedSeq() << endl);
	    CEBUG(actread.getName() << "\tLIML: " << actread.getLeftClipoff() << '\t' << te.from << endl;);
	    CEBUG(actread.getName() << "\tLIMR: " << actread.getRightClipoff() << '\t' << te.to << endl);
	    //Read::setCoutType(Read::AS_TEXT);
	    //CEBUG(actread);

	    bool wouldold=false;
	    if(te.to-te.from+1>=actread.getLenClippedSeq()) wouldold=true;
	    if(actread.getLeftClipoff()>=te.from
	       && actread.getRightClipoff()-1 <= te.to){
	      ++chkall;
	      CEBUG("Next read:\n");
	      auto repcov=hsd.estimDigiNormCov(actread);
	      CEBUG("repcov: " << repcov << endl);
	      if(repcov>1){
		auto newtag=te;
		// WARNING: with this we'll probably break the for(auto & te ...) functionality
		// we MUST get out of the loop afterwards with a "break"!
		actread.deleteTag(Read::REA_tagentry_idDGNr);
		string comnum(boost::lexical_cast<string>(repcov));
		newtag.setCommentStr(comnum);
		actread.addTagO(newtag);
		CEBUG("\nDGN repcov " << actread.getName() << ": " << repcov << endl);
		//cout << "DGN repcov " << actread.getName() << ": " << repcov << endl;
		break;
	      }
	    }else{
	      if(wouldold){
		CEBUG("Gaaaaaaaaah, would have taken???\n");
	      }
	      // remove any DGNr coverage number
	      te.setCommentStr("");
	    }
	  }
	}
      }
    }
    CEBUG("Chkall: " << chkall << endl);
  }

  cout << "\nDigital normalisation: removed " << numnormout << " reads.\n";

  if(debrisreasonptr != nullptr){
    auto & db=*debrisreasonptr;
    BUGIFTHROW(db.size()!=rp.size(),"db.size()!=rp.size() ???");
    for(auto rpi=0; rpi<rp.size(); ++rpi){
      if(db[rpi]==0 && normout[rpi]){
	db[rpi]=Assembly::DEBRIS_DIGITAL_NORMALISATION;
      }
    }
  }

}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsureAdapRegexes(ReadGroupLib::ReadGroupID rgid)
{
  if(DP_adapres.size()>rgid.getLibId()
     && DP_adapres[rgid.getLibId()].areinit) return;

  if(DP_adapres.size()<=rgid.getLibId()) DP_adapres.resize(rgid.getLibId()+1);

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
    static const char regexfile[] = {
#include "adaptorsregex.solexa.xxd.H"
      ,0
    };
    addAdapRegexes(rgid,regexfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    static const char regexfile[] = {
#include "adaptorsregex.iontor.xxd.H"
      ,0
    };
    addAdapRegexes(rgid,regexfile);
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::addAdapRegexes(ReadGroupLib::ReadGroupID rgid, const char * regexfile)
{
  FUNCSTART("void priv_constructorAdapRes(uint8 seqtype, const char * regexfile)");

  BUGIFTHROW(rgid.getLibId()>=ReadGroupLib::getNumReadGroups(),"Oooops, readgroupid " << static_cast<uint16>(rgid.getLibId()) << " is unknown?");

  if(DP_adapres.size()<=rgid.getLibId()) DP_adapres.resize(rgid.getLibId()+1);

  CEBUG("prepping regexp for " << rgid.getLibId() << endl);

  masterslavere_t tmpmsre;
  istringstream tmpis(regexfile);
  string line;
  while(true){
    getline(tmpis,line);
    if(tmpis.eof()) break;
    if(line[0]=='>'){
      DP_adapres[rgid.getLibId()].adapres.push_back(tmpmsre);
      line.erase(0,1);         // get away the ">"
      boost::trim(line);
      if(!line.empty()){
	boost::to_upper(line);
	DP_adapres[rgid.getLibId()].adapres.back().masterre=boost::regex(line);
	DP_adapres[rgid.getLibId()].adapres.back().hasmaster=true;
      }
    }else{
      BUGIFTHROW(DP_adapres[rgid.getLibId()].adapres.empty(),"Oooops, no master expression found?");
      boost::to_upper(line);
      DP_adapres[rgid.getLibId()].adapres.back().slaveres.push_back(boost::regex(line));
    }
  }
  DP_adapres[rgid.getLibId()].areinit=true;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsureAdapSkims(ReadGroupLib::ReadGroupID rgid)
{
  if(DP_adapskims.size()>rgid.getLibId()
     && DP_adapskims[rgid.getLibId()].skimptr) return;

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA){
    static const char adapfile[] = {
#include "adaptorsforclip.solexa.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    static const char adapfile[] = {
#include "adaptorsforclip.iontor.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }else if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_454GS20){
    static const char adapfile[] = {
#include "adaptorsforclip.454.xxd.H"
      ,0
    };
    priv_constructorSkimPool(rgid,DP_adapskims,7,adapfile);
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::priv_constructorSkimPool(ReadGroupLib::ReadGroupID rgid, std::vector<poolskim_t> & skimpool, const uint8 basesperhash, const char * adapfile)
{
  FUNCSTART("void DataProcessing::priv_constructorAdapPool(uint8 seqtype, const char * adapfile)");
  BUGIFTHROW(rgid.getLibId()>=ReadGroupLib::getNumReadGroups(),"Oooops, readgroupid " << static_cast<uint16>(rgid.getLibId()) << " is unknown, have " << ReadGroupLib::getNumReadGroups() << " read groups ?");

  boost::mutex::scoped_lock lock(DP_ps_changemutex);
  if(skimpool.size()<=rgid.getLibId()) skimpool.resize(rgid.getLibId()+1);

  if(skimpool[rgid.getLibId()].poolptr==nullptr){
    CEBUG("prepping adappool " << rgid.getLibId() << endl);
    skimpool[rgid.getLibId()].poolptr= new ReadPool(DP_miraparams_ptr);

    istringstream tmpis(adapfile);
    string nline,sline;
    bool addedsomething=false;
    while(true){
      getline(tmpis,nline);
      if(tmpis.eof()) break;
      boost::trim(nline);
      nline.erase(0,1);         // get away the ">"
      if(!nline.empty()){
	getline(tmpis,sline);
	boost::trim(sline);
	boost::to_upper(sline);
	size_t ereadidx=0;
	for( ; ereadidx < skimpool[rgid.getLibId()].poolptr->size(); ++ereadidx){
	  //CEBUG("COMP: " << skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getName()<< "\t" << skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getSeqAsChar() << "\t" << sline << endl);
	  if(skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getLenSeq() == sline.size()
	     && strncmp(skimpool[rgid.getLibId()].poolptr->getRead(ereadidx).getSeqAsChar(),sline.c_str(),sline.size())==0) {
	    //CEBUG("BINGO!\n");
	    break;
	  }
	}
	if(ereadidx==skimpool[rgid.getLibId()].poolptr->size()){
	  CEBUG("New for " << nline<<endl);
	  ereadidx=skimpool[rgid.getLibId()].poolptr->provideEmptyRead();
	  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(ereadidx);
	  actread.disallowAdjustments();
	  actread.setName(nline);
	  if(tmpis.eof()) break;
	  actread.setSequenceFromString(sline);
	  addedsomething=true;
	  actread.setUsedInAssembly(false);
	}else{
	  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(ereadidx);

	  CEBUG("Extend for " << actread.getName() << ": " << nline<<endl);

	  actread.setName(actread.getName()+"_/_"+nline);
	  //actread.setUsedInAssembly(true);
	}
      }
    }

    //for(auto idx=0; idx < skimpool[rgid.getLibId()].poolptr->size(); ++idx){
    //  Read & actread=skimpool[rgid.getLibId()].poolptr->getRead(idx);
    //  if(actread.isUsedInAssembly()){
    //	cout << actread.getName() << endl;
    //  }
    //}
    //exit(0);

    if(addedsomething){
      CEBUG("prepping skim " << rgid.getLibId() << endl);
      if(skimpool[rgid.getLibId()].skimptr!=nullptr) delete skimpool[rgid.getLibId()].skimptr;

      // IMPORTANT: keep assignment of skimpool[rgid.getLibId()].skimptr
      //  as last thing to be done in this branch, i.e., when the Skim is
      //  completely initialised
      // Reason: the checks in priv_ensureAdapSkims() do not use any mutexes, just
      //  the vector size AND the skimptr being nullptr. It could be that
      //  the constructing thread therefore did not finish constructing
      //  the Skim but that another thread would jump forward using it because
      //  everything points to that it's ready
      auto skimptr = new Skim();
      skimptr->skimStreamPrepare(*skimpool[rgid.getLibId()].poolptr,basesperhash,1);
      skimptr->prepareForMultithreadFarc((*DP_miraparams_ptr)[0].getSkimParams().sk_numthreads);
      skimpool[rgid.getLibId()].skimptr = skimptr;

      CEBUG("Done\n");
    }
  }
}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::priv_EnsurePhiX174Statistics()
{
  if(DP_px174hs_init) return;

  static const char adapfile[] = {
#include "seqforfilter_phix174.solexa.xxd.H"
    ,0
  };


  boost::mutex::scoped_lock lock(DP_px174hs_changemutex);
  if(DP_px174hs_init) return;

  ReadPool baitrp(DP_miraparams_ptr);

  istringstream tmpis(adapfile);
  string line;
  bool addedsomething=false;
  while(true){
    getline(tmpis,line);
    if(tmpis.eof()) break;
    line.erase(0,1);         // get away the ">"
    if(!line.empty()){
      size_t ereadidx=baitrp.provideEmptyRead();
      Read & actread=baitrp[ereadidx];
      actread.disallowAdjustments();
      actread.setName(line);
      getline(tmpis,line);
      //CEBUG("For " << actread.getName() << ": " << line<<endl);
      if(tmpis.eof()) break;
      actread.setSequenceFromString(line);
      addedsomething=true;
    }
  }

  if(addedsomething){
    string dummyfn;
    DP_phix174hashstatistics.prepareHashStatistics((*DP_miraparams_ptr)[0].getDirectoryParams().dir_tmp,baitrp,false,false,true,true,1,31,1,false,dummyfn);
    boost::mutex tmpm;
    boost::mutex::scoped_lock lock(tmpm);
    DP_px174hs_init=true;
    CEBUG("Done\n");
  }

}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::startLogging(const std::string filename, bool newfile)
{
  FUNCSTART("void DataProcessing::startLogging(const std::string filename, bool newfile)");
  stopLogging();
  if(!filename.empty()){
    DP_logname=newfile;
    if(newfile){
      DP_logfout.open(filename.c_str(), ios::out|ios::trunc);
    }else{
      DP_logfout.open(filename.c_str(), ios::out|ios::app);
    }
    if(!DP_logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open " << filename << " for logging.");
    }
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::stopLogging()
{
  if(DP_logfout.is_open()){
    DP_logfout.close();
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::performRareKMERMasking_Pool(ReadPool & rpool, const string & logprefix)
{
  auto & miraparams = *DP_miraparams_ptr;

  bool needtomask=false;
  for(auto & mp : miraparams){
    if(mp.getAssemblyParams().as_clipmask_rarekmers>0){
      needtomask=true;
    }
  }

  if(!needtomask) return;

  uint8 basesperhash=miraparams[0].getSkimParams().sk_basesperhash;

  cout << "Rare kmer masking ... ";cout.flush();
  for(uint32 rpi=0; rpi<rpool.size(); ++rpi){
    Read & actread= rpool[rpi];
    if(!actread.hasValidData()
       || !actread.isUsedInAssembly()
       || miraparams[actread.getSequencingType()].getAssemblyParams().as_clipmask_rarekmers==0 ) continue;
    performRareKMERMasking_Read(actread,basesperhash, logprefix);

  }

  cout << "done\n";
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::performRareKMERMasking_Read(Read & actread, uint8 basesperhash, const string & logprefix)
{
  FUNCSTART("void DataProcessing::performRareKMERMasking_Read(Read & actread, uint8 basesperhash, const string & logprefix)");

  if(actread.getLenSeq()==0) return;

  if(!actread.hasBaseHashStats()){
    Read::setCoutType(Read::AS_TEXT);
    cout << actread << endl;
    MIRANOTIFY(Notify::FATAL,"!actread.hasBaseHashStats() ??? ");
  }

  DP_tmpvu8.clear();
  DP_tmpvu8.resize(actread.getLenSeq(),0);

  vector<Read::bposhashstat_t>::const_iterator bhsI=actread.getBPosHashStats().begin();
  vector<Read::bposhashstat_t>::const_iterator bhsE=actread.getBPosHashStats().end();

  vector<uint8>::iterator tfI=DP_tmpvu8.begin();
  vector<uint8>::iterator tfE=DP_tmpvu8.end();

  bool mustmask=false;
  for(uint32 readpos=0; bhsI!= bhsE; ++bhsI, ++tfI, ++readpos){
    if(bhsI->fwd.getFrequency()==1){
      vector<uint8>::iterator ttfI=tfI;
      uint32 tmpreadpos=readpos;
      for(uint32 i=0; i<basesperhash && ttfI!=tfE; i++, ttfI++, ++tmpreadpos){
	if(tmpreadpos>=actread.getLeftClipoff() && tmpreadpos<=actread.getRightClipoff()){
	  mustmask=true;
	  *ttfI=1;
	}
      }
    }
  }

  if(mustmask){
    // "gnawing" (basesperhash/2-2 bases from each masked stretch if that stretch
    //  is >= basesperhash long, but not at read ends)
    {
      auto numgnaw=basesperhash/2;
      if(numgnaw>2) numgnaw-=2;

      CEBUG("NEW rgna: " << actread.getName() << " " << actread.getLenSeq() << " " << numgnaw << endl);

      auto sI=DP_tmpvu8.begin();

      // jump over first masked stretch if it starts at read start
      while(sI!=DP_tmpvu8.end() && *sI) {++sI;}

      CEBUG("after first jump: " << sI-DP_tmpvu8.begin());

      for(; sI!=DP_tmpvu8.end();){
	//jump over unmarked stretch
	for(;sI!=DP_tmpvu8.end() && *sI==0; ++sI){};
	CEBUG("After unmarked jump: " << sI-DP_tmpvu8.begin() << endl);

	// find end of marked stretch
	auto eI=sI;
	for(;eI!=DP_tmpvu8.end() && *eI; ++eI){};
	CEBUG("After marked jump: " << eI-DP_tmpvu8.begin() << endl);
	CEBUG("size: " << eI-sI << endl);

	// gnaw if needed and / or set iterator to begin of next search
	if(eI-sI >= basesperhash && eI!=DP_tmpvu8.end()){
	  CEBUG("GNAW IT\n");
	  for(uint8 i=0; i<numgnaw; ++i, ++sI) {*sI=0;}
	  auto jI=eI;
	  for(uint8 i=0; i<numgnaw; ++i) {--eI; *eI=0;}
	  sI=jI;
	}else{
	  sI=eI;
	}
      }
    }

    // and log changes
    bool mustoutputprefix=true;
    for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
      if(DP_tmpvu8[readpos]){
	auto start=readpos;
	while(readpos<DP_tmpvu8.size() && DP_tmpvu8[readpos]) ++readpos;
	if(readpos>start) --readpos;
	if(mustoutputprefix){
	  DP_logfout << logprefix
		     << " rare kmer "
		     << actread.getName();
	  mustoutputprefix=false;
	}
	DP_logfout << "\t["
		   << start
		   << ".."
		   << readpos
		   << ']';
      }
    }
    if(!mustoutputprefix) DP_logfout << '\n';

    // masking with 'x'

    for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
      if(DP_tmpvu8[readpos]) actread.changeBaseInSequence('x',0,readpos);
    }

    auto gaplen=basesperhash;
    maskClips_Read(actread, logprefix,gaplen,gaplen,gaplen);

    // TODO: first mask with X so that maskClips_Read() can do it's work, then use 'n' is not
    //   really elegant. Alternative would be yet another parameter to maskClips_Read()
    for(uint32 readpos=0; readpos<DP_tmpvu8.size(); ++readpos){
      if(DP_tmpvu8[readpos]) actread.changeBaseInSequence('n',0,readpos);
    }
  }
}
//#define CEBUG(bla)






/*************************************************************************
 *
 * expects reads to have baseflags set
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}

struct tmpbhentry_t{
  uint32 from;
  uint32 to;
  uint8 freq;
};

void DataProcessing::buntifyReadsByHashFreq_Pool(ReadPool & rp, uint8 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq()");

  cout << "Buntifying reads";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  for(uint32 actid=0; actid<rp.size(); actid++){
    Read & actread=rp[actid];

    buntifyReadsByHashFreq_Read(rp[actid],basesperhash);
  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)


void DataProcessing::buntifyReadsByHashFreq_Read(Read & actread, uint8 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq_Read(Read & actread, uint8 basesperhash)");

  //Read::setCoutType(Read::AS_TEXT);
  //cout << actread;

  // remove old hash frequence tags
  for(uint32 i=0; i<DP_freqnames.size(); i++){
    actread.deleteTag(DP_freqnames[i]);
  }

  if(actread.hasValidData()
     && actread.hasBaseHashStats()){

    static multitag_t tmpmt("","","MIRA");

    DP_tmpvu8.clear();
    DP_tmpvu8.resize(actread.getLenSeq(),0);

    vector<Read::bposhashstat_t>::const_iterator bhsI=actread.getBPosHashStats().begin();
    vector<Read::bposhashstat_t>::const_iterator bhsE=actread.getBPosHashStats().end();
    vector<uint8>::iterator tfI=DP_tmpvu8.begin();
    vector<uint8>::iterator tfE=DP_tmpvu8.end();

    priv_buntifyHelper(2, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(3, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(4, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(5, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(6, basesperhash, bhsI, bhsE, tfI, tfE);
    priv_buntifyHelper(7, basesperhash, bhsI, bhsE, tfI, tfE);

    vector<tmpbhentry_t> telist;
    telist.reserve(20);

    //{
    //	cout << "bfr: " << actread.getName() << endl;
    //	for(uint32 i=0;i<DP_tmpvu8.size(); i++){
    //	  cout << "i: " << i << '\t' << static_cast<uint16>(DP_tmpvu8[i]) << endl;
    //	}
    //}

    uint32 from=0;
    uint32 to=0;
    for(; from<actread.getLenSeq(); from=to+1){
      to=from;
      uint8 actfreq=DP_tmpvu8[to];
      for(; to<actread.getLenSeq() && DP_tmpvu8[to]==actfreq; to++) {} ;
      to--;
      if(actfreq>0){
	telist.resize(telist.size()+1);
	telist.back().from=from;
	telist.back().to=to;
	telist.back().freq=actfreq;
      }
    }

    // for first or last entry, do not put tags for frequencies
    //  >=2 if their length is < basesperhash
    // BaCh 03.06.2011: hmmm, why not. OK, makes CAF/MAF bigger, but else?
    for(uint32 ti=0; ti<telist.size(); ti++){
      bool settag=true;
//	if(telist[ti].freq>=2 &&
//	   (ti==0 || ti==telist.size()-1)){
//	  if(telist[ti].to - telist[ti].from < basesperhash-1){
//	    settag=false;
//	  }
//	}
      if(settag) {
	tmpmt.identifier=DP_freqnames[telist[ti].freq];
	tmpmt.from=telist[ti].from;
	tmpmt.to=telist[ti].to;
	actread.addTagO(tmpmt);
      }
    }
  }

}
//#define CEBUG(bla)


void DataProcessing::priv_buntifyHelper(uint8 allowedfreq, uint8 basesperhash, vector<Read::bposhashstat_t>::const_iterator bhsI, vector<Read::bposhashstat_t>::const_iterator bhsE, vector<uint8>::iterator tfI, vector<uint8>::iterator tfE)
{
  for(; bhsI!= bhsE; bhsI++, tfI++){
    uint8 actfreq=bhsI->fwd.getFrequency();
    if(allowedfreq==actfreq){
      if(actfreq>0){
	vector<uint8>::iterator ttfI=tfI;
	for(uint32 i=0; i<basesperhash && ttfI!=tfE; i++, ttfI++){
	  *ttfI=actfreq;
	}
      }
    }
  }
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/


void DataProcessing::addKMerForkTags_Pool(ReadPool & rp, uint8 basesperhash)
{
  FUNCSTART("void DataProcessing::addKMerForkTags_Pool(ReadPool & rp, uint8 basesperhash)");

  cout << "Adding fork tags";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  static multitag_t tmpmt("","","MIRA");
  tmpmt.identifier=Read::REA_tagentry_idKMRF;

  for(uint32 actid=0; actid<rp.size(); actid++){
    Read & actread=rp[actid];

    //Read::setCoutType(Read::AS_TEXT);
    //cout << actread;

    // remove old KMRF tags
    actread.deleteTag(tmpmt.identifier);

    if(actread.hasValidData()
       && actread.hasBaseHashStats()){

      DP_tmpvu8.clear();
      DP_tmpvu8.resize(actread.getLenSeq(),0);

      auto bhsI=actread.getBPosHashStats().cbegin();
      auto bhsE=actread.getBPosHashStats().cend();
      auto tfI=DP_tmpvu8.begin();
      auto tfE=DP_tmpvu8.end();

      for(; bhsI!= bhsE; bhsI++, tfI++){
	if(bhsI->fwd.isKMerFork()){
	  auto ttfI=tfI;
	  for(uint32 i=0; i<basesperhash && ttfI!=tfE; ++i, ++ttfI){
	    *ttfI=1;
	  }
	}
	if(bhsI->rev.isKMerFork()){
	  auto ttfI=tfI;
	  for(uint32 i=0; i<basesperhash; ++i, --ttfI){
	    *ttfI=1;
	    if(ttfI!=DP_tmpvu8.begin()) break;
	  }
	}
      }

      uint32 from=0;
      uint32 to=0;
      for(; from<actread.getLenSeq(); from=to+1){
	to=from;
	if(DP_tmpvu8[to]){
	  for(; to<actread.getLenSeq() && DP_tmpvu8[to]; ++to) {} ;
	  to--;
	  tmpmt.from=from;
	  tmpmt.to=to;
	  actread.addTagO(tmpmt);
	}
      }
    }

    //if(actid==273250 || actid==273252){
    //  Read::setCoutType(Read::AS_TEXT);
    //  cout << actread;
    //}

  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/


void DataProcessing::performKMERRepeatTagging_Pool(ReadPool & rp, uint8 basesperhash)
{
  FUNCSTART("void DataProcessing::buntifyReadsByHashFreq()");

  cout << "Adding RMB tags by fork";
  if(rp.size()>500000) cout << " (this may take a while)";
  cout << " ... "; cout.flush();

  for(uint32 actid=0; actid<rp.size(); actid++){
    if(1){
      performKMERRepeatTagging_Read(rp[actid],basesperhash);
    }
  }

  cout << "done." << endl;

  FUNCEND();

}
//#define CEBUG(bla)

void DataProcessing::performKMERRepeatTagging_Read(Read & actread, uint8 basesperhash)
{
  FUNCSTART("void DataProcessing::addBla_Read(Read & actread, uint8 basesperhash)");
  if(actread.hasValidData()
     && actread.hasBaseHashStats()
     && actread.getLenSeq() >= 2*basesperhash){

    static multitag_t tmpmt("","addBla","MIRA");
    tmpmt.identifier=Read::REA_tagentry_idCRMr;

    DP_tmpvu8.clear();
    DP_tmpvu8.resize(actread.getLenSeq(),0);

    auto bhsE=actread.getBPosHashStats().cend();
    auto bhsIf=actread.getBPosHashStats().cbegin();
    uint32 tagpos=basesperhash-1;
    auto bhsIr=bhsIf+2*tagpos;
    auto tfI=DP_tmpvu8.begin()+tagpos;


    for(; bhsIr!= bhsE; ++bhsIf, ++bhsIr, ++tfI, ++tagpos){
      if(bhsIf->fwd.isKMerFork()
	 && bhsIr->rev.isKMerFork()){
	*tfI=1;
      }
    }

    uint32 runcount=0;
    for(tfI=DP_tmpvu8.begin(); tfI!=DP_tmpvu8.end(); ++tfI){
      if(*tfI){
	++runcount;
      }else{
	if(runcount==1){
	  tmpmt.from=tfI-1-DP_tmpvu8.begin();
	  tmpmt.to=tmpmt.from;
	  actread.addTagO(tmpmt);
	}
	runcount=0;
      }
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::clipBadSolexaEnds_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::clipBadSolexaEnds_Pool(ReadPool & rp, const string & logprefix)");

  for(uint32 i=0;i<rp.size();++i){
    Read & actread=rp[i];
    if(actread.hasValidData()
       && actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)
       && !(actread.isBackbone()
	    || actread.isRail())){
      clipBadSolexaEnds_Read(actread,logprefix);
    }
  }
}


/*************************************************************************
 *
 * TODO: not really good, rethink that for eukaryotes
 *
 *************************************************************************/

void DataProcessing::clipBadSolexaEnds_Read(Read & actread, const string & logprefix)
{
  FUNCSTART("void DataProcessing::clipBadSolexaEnds_Read(Read & actread, const string & logprefix)");

  // invalidate all Solexa reads that have a stretch of 20 A (or T)
  //  or if has stretch >= 12 and non-A (or T) bases < 20%
  // N in-between do not reset the counter
  // invalidate by setting left seq vec to length of read

  int32 runindex=actread.getLeftClipoff();
  char actbase=' ';
  uint32 bcount=0;

  uint32 arun=0;
  uint32 maxarun=0;
  uint32 nona=0;

  uint32 trun=0;
  uint32 maxtrun=0;
  uint32 nont=0;

  for(; runindex<actread.getRightClipoff(); runindex++) {
    actbase=static_cast<char>(toupper(actread.getBaseInSequence(runindex)));
    if(actbase!='N'){
      bcount++;
      if(actbase=='A'){
	arun++;
	if(arun>maxarun) maxarun=arun;
	nont++;
	trun=0;
      }else if(actbase=='T'){
	trun++;
	if(trun>maxtrun) maxtrun=trun;
	nona++;
	arun=0;
      }else{
	nona++;
	nont++;
	arun=0;
	trun=0;
      }
    }
  }
  if(maxarun>=20){
    actread.setLSClipoff(actread.getLenSeq());
    DP_logfout << logprefix << " bad solexa end: A hard "
	       << actread.getName()
	       << '\n';
  }else if(maxarun>=12){
    uint32 ratio= static_cast<uint32>((static_cast<double>(100.0)/bcount)*nona);
    if(ratio<20) {
      actread.setLSClipoff(actread.getLenSeq());
      DP_logfout << logprefix << " bad solexa end: A soft "
		 << actread.getName()
		 << '\n';
    }
  }

  if(maxtrun>=20){
    actread.setLSClipoff(actread.getLenSeq());
    DP_logfout << logprefix << " bad solexa end: T (hard) "
	       << actread.getName()
	       << '\n';
  }else if(maxtrun>=12){
    uint32 ratio= static_cast<uint32>((static_cast<double>(100.0)/bcount)*nont);
    if(ratio<20) {
      actread.setLSClipoff(actread.getLenSeq());
      DP_logfout << logprefix << " bad solexa end: T (soft) "
		 << actread.getName()
		 << '\n';
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * clip all lowercase at the end of reads
 *
 *************************************************************************/

void DataProcessing::lowerCaseClipping_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClipping_Pool(ReadPool & rp, const string & logprefix)");

  uint64 totallen=0;
  uint64 lowercaselen=0;
  for(uint32 i=0;i<rp.size();i++){
    Read & actread=rp.getRead(i);
    if(actread.hasValidData()
       && ((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_front
	    || (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_back)
       && !(actread.isBackbone()
	    || actread.isRail())){
      totallen+=actread.getLenClippedSeq();
      int32 runindex=actread.getLeftClipoff();
      for(; runindex<actread.getRightClipoff(); ++runindex){
	if(islower(actread.getBaseInSequence(runindex))) lowercaselen++;
      }
    }
  }

  if(totallen==lowercaselen) {
    cout << "Lowercase clip: all sequences to be clipped are lowercase?! Failsafe: no clipping performed.\n";
    return;
  }

  for(uint32 i=0;i<rp.size();i++){
    Read & actread=rp.getRead(i);
    if(actread.hasValidData()
       && !(actread.isBackbone()
	    || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_front){
	lowerCaseClippingFront_Read(actread,logprefix);
      }
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_lowercase_back){
	lowerCaseClippingBack_Read(actread,logprefix);
      }
    }
  }

  FUNCEND();
}

/*************************************************************************
 *
 * clip all lowercase at the end of reads
 *
 *************************************************************************/

void DataProcessing::lowerCaseClippingFront_Read(Read & actread, const string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClippingFront_Read(Read & actread, const string & logprefix)");

  // TODO: for streaming, implement check counter by rgid for all sequence lowercase

  int32 runindex=actread.getLeftClipoff();
  for(; runindex<actread.getRightClipoff(); ++runindex){
    char ab=actread.getBaseInSequence(runindex);
    if(!islower(ab)
       && ab != 'N'
       && ab != 'X') break;
  }
  // TODO: 01.01.2013 check this, changed != to >
  if(runindex>actread.getLeftClipoff()) {
    actread.setLSClipoff(runindex);
    DP_logfout << logprefix << " changed left (lowercase) "
	       << actread.getName() << " to " << actread.getLeftClipoff() << '\n';
  }
}


void DataProcessing::lowerCaseClippingBack_Read(Read & actread, const string & logprefix)
{
  FUNCSTART("void DataProcessing::lowerCaseClippingBack_Read(Read & actread, const string & logprefix)");

  // TODO: for streaming, implement check counter by rgid for all sequence lowercase

  int32 runindex=actread.getRightClipoff()-1;
  for(; runindex>=actread.getLeftClipoff() && islower(actread.getBaseInSequence(runindex)); --runindex) ;
  // TODO: implement jumping over N,X like for front
  // TODO: 01.01.2013 really check this, changed != to <
  if(runindex<actread.getRightClipoff()-1) {
    actread.setRSClipoff(runindex+1);
    DP_logfout << logprefix << " changed right (lowercase) "
	       << actread.getName() <<  " to " << actread.getRightClipoff() << '\n';
  }
  //cout << actread;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::qualClips_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::qualClips_Pool(ReadPool & rp, const string & logprefix)");

  cout << "Starting qual clips: ";

  for(uint32 i=0;i<rp.size();i++){
    Read & r=rp.getRead(i);
    if(r.hasValidData()
       && !(r.isBackbone()
	    || r.isRail())){
      if((*DP_miraparams_ptr)[r.getSequencingType()].getAssemblyParams().as_clip_quality) {
	qualClips_Read(r,logprefix);
      }
    }
  }
}


void DataProcessing::qualClips_Read(Read & actread, const string & logprefix)
{
  FUNCSTART("void DataProcessing::qualClips_Read(Read & actread, const string & logprefix)");

  int32 oldlq=actread.getLQClipoff();
  int32 oldrq=actread.getRQClipoff();
  actread.performQualityClip(
    (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_minqual,
    (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_winlen);
  int changed=2;
  if(oldlq>actread.getLQClipoff()){
    // old LQ was more conservative, change back
    actread.setLQClipoff(oldlq);
    --changed;
  }
  if(oldrq<actread.getRQClipoff()){
    // old RQ was more conservative, change back
    actread.setRQClipoff(oldrq);
    --changed;
  }
  if(changed){
    DP_logfout << logprefix
	       << " changed qual. "
	       << actread.getName()
	       << "\tfrom: " << oldlq << ' ' << oldrq << "\tto: "
	       << actread.getLQClipoff()
	       << ' '
	       << actread.getRQClipoff() << '\n';
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::maskClips_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::maskClips_Pool(ReadPool & rp, const string & logprefix)");

  cout << "Starting qual clips: ";

  for(uint32 i=0;i<rp.size();i++){
    Read & r=rp.getRead(i);
    if(r.hasValidData()
       && !(r.isBackbone()
	    || r.isRail())){
      if((*DP_miraparams_ptr)[r.getSequencingType()].getAssemblyParams().as_clip_maskedbases) {
	maskClips_Read(r,logprefix);
      }
    }
  }
}


void DataProcessing::maskClips_Read(Read & actread, const string & logprefix, int32 gapsize, int32 maxfrontgap, int32 maxendgap)
{
  FUNCSTART("void DataProcessing::maskClips_Read(Read & actread, const string & logprefix, int32 gapsize, int32 maxfrontgap, int32 maxendgap)");

  if(gapsize<0) gapsize=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_gapsize;
  if(maxfrontgap<0) maxfrontgap=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_maxfrontgap;
  if(maxendgap<0) maxendgap=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_maskedbase_maxendgap;

  auto oldlm=actread.getLMClipoff();
  auto oldrm=actread.getRMClipoff();
  auto oldlc=actread.getLeftClipoff();
  auto oldrc=actread.getRightClipoff();
  actread.setClipoffsToMaskedChars(
    gapsize,
    maxfrontgap,
    maxendgap,
    false);
  actread.setClipoffsToMaskedChars(
    1,
    1,
    1,
    true);
  int changed=2;
  if(oldlm>=actread.getLMClipoff()){
    // old LM was more conservative, change back
    if(oldlm>actread.getLMClipoff()) actread.setLMClipoff(oldlm);
    --changed;
  }
  if(oldrm<=actread.getRMClipoff()){
    // old RM was more conservative, change back
    if(oldrm<actread.getRMClipoff()) actread.setRMClipoff(oldrm);
    --changed;
  }
  if(changed){
    DP_logfout << logprefix
	       << " changed mask. "
	       << actread.getName()
	       << "\tfrom: " << oldlm << ' ' << oldrm << " (" << oldlc << " " << oldrc << ")\tto: "
	       << actread.getLMClipoff()
	       << ' '
	       << actread.getRMClipoff() << " (" << actread.getLeftClipoff() << " " << actread.getRightClipoff() << ")\n";
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

void DataProcessing::minimumQualityThreshold_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::minimumQualityThreshold_Pool(ReadPool & rp, const string & logprefix)");

  cout << "Starting minimum quality threshold clip ... "; cout.flush();

  uint32 numkilled=0;

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if(!minimumQualityThreshold_Read(actread,logprefix)){
	++numkilled;
      }
    }
  }
  cout << "done. Killed " << numkilled << " reads.\n";

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool DataProcessing::minimumQualityThreshold_Read(Read & actread, const string & logprefix)
{
  FUNCSTART("bool DataProcessing::minimumQualityThreshold_Read(Read & actread, const string & logprefix)");

  base_quality_t minqual=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_minthreshold;
  uint32 minnum=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_quality_numminthreshold;
  const vector<base_quality_t> & quals=actread.getQualities();
  vector<base_quality_t>::const_iterator qI=quals.begin();
  bool mustkill=true;
  for(; qI != quals.end(); ++qI){
    if(*qI>=minqual
       && --minnum==0) {
      mustkill=false;
      break;
    }
  }

  if(mustkill){
    actread.setLQClipoff(actread.getLenSeq());
    actread.setRQClipoff(actread.getLenSeq());
    DP_logfout << logprefix << ' ' << actread.getName() << ": min qual threshold not met, killed\n";
  }

  FUNCEND();
  return !mustkill;
}






/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::minimumLeftClip_Pool(ReadPool & rp, const string & logprefix, bool qual, bool seqvec, bool mask)
{
  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_ensureminimumleftclipoff){
	minimumLeftClip_Read(actread,logprefix,qual,seqvec,mask);
      }
    }
  }
}

void DataProcessing::minimumLeftClip_Read(Read & actread, const string & logprefix, bool qual, bool seqvec, bool mask)
{
  auto oldlc=actread.getLeftClipoff();
  if(oldlc < (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minslrequired){
    if(qual) actread.setLQClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);
    if(seqvec) actread.setLSClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);
    if(mask) actread.setLMClipoff((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minqlsetto);

    DP_logfout << logprefix
	       << " changed minleft. "
	       << actread.getName()
	       << "\tLeft: "
	       << oldlc
	       << "\t -> "
	       << actread.getLeftClipoff()
	       << '\n';
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::minimumRightClip_Pool(ReadPool & rp, const string & logprefix, bool qual, bool seqvec, bool mask)
{
  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_ensureminimumrightclipoff){
	minimumRightClip_Read(actread,logprefix,qual,seqvec,mask);
      }
    }
  }
}

void DataProcessing::minimumRightClip_Read(Read & actread, const string & logprefix, bool qual, bool seqvec, bool mask)
{
  auto oldrc=actread.getRightClipoff();
  if(actread.getLenSeq()-oldrc < (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minsrrequired){
    int32 newr=actread.getLenSeq()-(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_minsrrequired;
    if(qual) actread.setRQClipoff(newr);
    if(seqvec) actread.setRSClipoff(newr);
    if(mask) actread.setRMClipoff(newr);

    DP_logfout << logprefix
	       << " changed minRight. "
	       << actread.getName()
	       << "\tRightt: "
	       << oldrc
	       << "\t -> "
	       << actread.getRightClipoff()
	       << '\n';
  }
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void DataProcessing::badSequenceSearch_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::badSequenceSearch_Pool(ReadPool & rp, const string & logprefix)");

  cout << "Performing search for bad sequence quality ... "; cout.flush();

  for(uint32 ri=0;ri<rp.size();++ri){
    Read & actread=rp[ri];
    if(actread.hasValidData()
       && actread.hasQuality()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality){
      }
    }
  }
}

void DataProcessing::badSequenceSearch_Read(Read & actread, const string & logprefix)
{
  uint32 winlen=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality_winlen;
  base_quality_t minqual=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_badstretchquality_minqual;

  const vector<base_quality_t> & bquals=actread.getQualities();
  int32 runi=actread.getLeftClipoff();
  int32 endi=actread.getRightClipoff();

  uint32 qualsbelow=0;
  bool foundbad=false;
  for(; runi < endi; runi++){
    if(bquals[runi] < minqual){
      ++qualsbelow;
      if(qualsbelow >= winlen){
	foundbad=true;
	break;
      }
    }else{
      qualsbelow=0;
    }
  }

  if(foundbad) {
    int32 newrclip=runi-qualsbelow+1;
    int32 shortened=actread.getRightClipoff()-newrclip;
    //cout << actread.getName() << " has bad stretch, shortening by " << actread.getRightClipoff()-newrclip << '\n';
    if(newrclip < actread.getLQClipoff()) newrclip=actread.getLQClipoff();
    actread.setRQClipoff(newrclip);
    DP_logfout << logprefix << " bad seq. "
	       << actread.getName()
	       << "\tShortened by " << shortened
	       << "\tNew right: "
	       << actread.getRQClipoff()
	       << '\n';
  }

  FUNCEND();
  return;
}





/*************************************************************************
 *
 * clip poly-A in forward and poly-T in reverse direction
 * or: clip only after the poly-stretches, and tag the stretches with Fpas
 *
 * If poly stretches are kept, they are also "cleaned", i.e., bases are
 *  forced to be A (or T) in the detected stretch
 *
 *************************************************************************/

void DataProcessing::clipPolyATAtEnds_Pool(ReadPool & rp, const string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyATAtEnds_Pool(ReadPool & rp, const string & logprefix)");

  cout << "Clipping or tagging poly A/T stretches at ends of reads ... ";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_polyat){
	clipPolyATAtEnds_Read(actread,logprefix);
      }
    }
  }
}

void DataProcessing::clipPolyATAtEnds_Read(Read & actread, const string & logprefix)
{
  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  auto & as_params=(*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  uint32 mincount=as_params.as_clip_polyat_len;
  uint32 maxbad=as_params.as_clip_polyat_maxerrors;
  int32 grace=static_cast<int32>(as_params.as_clip_polyat_maxgap);
  bool keepstretch=as_params.as_clip_polyat_keeppolystretch;

  // search poly-a in forwarddirection
  {
    int32 lpolystart=-1;
    int32 rpolyend=-1;
    if(searchPolyBaseFrom5Prime(actread,'a',lpolystart,rpolyend,mincount,maxbad,grace)){
      if(keepstretch){
	actread.setRMClipoff(rpolyend+1);
	CEBUG("setting rm " << rpolyend+1 << endl);
	for(auto i=lpolystart; i<=rpolyend; ++i){
	  if(toupper(actread.getBaseInSequence(i))!='A'){
	    actread.changeBaseInSequence('a',0,i);
	  }
	}
      }else{
	actread.setRMClipoff(lpolystart);
	CEBUG("setting rm " << lpolystart << endl);
      }
      CEBUG("taggingl " << lpolystart << " " << rpolyend << endl);

      DP_logfout << logprefix << " poly-A fwd. "
		 << actread.getName()
		 << "\tMask right: "
		 << actread.getRMClipoff()
		 << '\n';

      DP_tmpmtpolyAT.from=lpolystart;
      DP_tmpmtpolyAT.to=rpolyend;
      actread.addTagO(DP_tmpmtpolyAT);
    }
  }

  // search poly-t in reverse direction
  {
    int32 lpolystart=-1;
    int32 rpolyend=-1;

    if(searchPolyBaseFrom3Prime(actread,'t',lpolystart,rpolyend,mincount,maxbad,grace)){
      if(keepstretch){
	actread.setLMClipoff(lpolystart);
	CEBUG("setting lm " << lpolystart << endl);
	for(auto i=lpolystart; i<=rpolyend; ++i){
	  if(toupper(actread.getBaseInSequence(i))!='T'){
	    actread.changeBaseInSequence('t',0,i);
	  }
	}
      }else{
	actread.setLMClipoff(rpolyend+1);
	CEBUG("setting lm " << rpolyend+1 << endl);
      }
      CEBUG("taggingl " << lpolystart << " " << rpolyend << endl);

      DP_logfout << logprefix << " poly-T rev. "
		 << actread.getName()
		 << "\tMask left: "
		 << actread.getLMClipoff()
		 << '\n';

      DP_tmpmtpolyAT.from=lpolystart;
      DP_tmpmtpolyAT.to=rpolyend;
      actread.addTagO(DP_tmpmtpolyAT);
    }
  }
}




/*************************************************************************
 *
 * Search poly-base (mincount length and maximum maxbad other bases) from left
 *  side of read (with 'grace' length grace if not encountered),
 *
 *  return:
 *    - true if found and return left and right coordinates in lpolystart and
 *      rpolyend
 *    - false if not found (lpolystart and rpolyend undefined)
 *
 *************************************************************************/

bool DataProcessing::searchPolyBaseFrom5Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)
{
  FUNCSTART("bool DataProcessing::searchPolyBaseFrom5Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)");


  BUGIFTHROW(!dptools::isValidACGTBase(polybase),"Ummm ... " << polybase << " is not ACGT?");
  BUGIFTHROW(grace<0,"grace (" << grace << ") < 0 ?");
  BUGIFTHROW(maxbad>=mincount,"maxbad (" << maxbad << ") >= mincount (" << mincount << ") ?");

  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  lpolystart=-1;
  rpolyend=-1;

  int32 runindex=actread.getLeftClipoff();
  int32 lastgoodrunindex=runindex;
  char actbase=' ';
  bool found=false;

  for(; grace >=0 && runindex<actread.getRightClipoff(); ++runindex, --grace) {
    actbase=actread.getBaseInSequence(runindex);
    CEBUG("###1 : " << grace << " " << runindex << "\t" << actbase << endl);
    if(dptools::areBasesContained(polybase,actbase)){
      lpolystart=runindex;
      lastgoodrunindex=runindex;
      uint32 acount=0;
      uint32 othercount=0;
      char cbase;
      for(; lastgoodrunindex<actread.getRightClipoff(); lastgoodrunindex++){
	cbase=actread.getBaseInSequence(lastgoodrunindex);
	if(dptools::areBasesContained(polybase,cbase)){
	  acount++;
	}else if(tolower(cbase)!='n'){
	  othercount++;
	  if(othercount>maxbad) break;
	}
      }
      if(acount>=mincount) {
	found=true;
	// get off non-poly characters as far as possible
	if(lastgoodrunindex==actread.getRightClipoff()) lastgoodrunindex--;
	while(lastgoodrunindex>runindex && !dptools::areBasesContained(polybase,actread.getBaseInSequence(lastgoodrunindex))) lastgoodrunindex--;
	rpolyend=lastgoodrunindex;
	break;
      }
      lpolystart=-1;
    }else{
      lpolystart=-1;
    }
  }

  if(rpolyend >=0 && lpolystart != -1) {
    FUNCEND();
    return true;
  }

  FUNCEND();
  return false;
}



/*************************************************************************
 *
 * Search poly-base (mincount length and maximum maxbad other bases) from left
 *  side of read (with 'grace' length grace if not encountered),
 *
 *  return:
 *    - true if found and return left and right coordinates in lpolystart and
 *      rpolyend
 *    - false if not found (lpolystart and rpolyend undefined)
 *
 * TODO: how dumb to have an own function instead of a templated function
 *       shared with 5p and working on iterators
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush();}
bool DataProcessing::searchPolyBaseFrom3Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)
{
  FUNCSTART("bool DataProcessing::searchPolyBaseFrom3Prime(Read & actread, const char polybase, int32 & lpolystart, int32 & rpolyend, const uint32 mincount, const uint32 maxbad, int32 grace)");


  BUGIFTHROW(!dptools::isValidACGTBase(polybase),"Ummm ... " << polybase << " is not ACGT?");
  BUGIFTHROW(grace<0,"grace (" << grace << ") < 0 ?");
  BUGIFTHROW(maxbad>=mincount,"maxbad (" << maxbad << ") >= mincount (" << mincount << ") ?");

  CEBUG(actread.getName() << endl);
  CEBUG(actread << endl);

  lpolystart=-1;
  rpolyend=-1;


  int32 runindex=actread.getRightClipoff()-1;
  int32 lastgoodrunindex=runindex;
  char actbase=' ';
  bool found=false;

  CEBUG("Reverse " << actread.getName() << '\n');

  for(; grace >=0 && runindex>=actread.getLeftClipoff(); --runindex, --grace) {
    actbase=static_cast<char>(tolower(actread.getBaseInSequence(runindex)));
    CEBUG("###1 : " << grace << " " << runindex << "\t" << actbase << endl);
    if(dptools::hasNucleicAcidInIUPAC(polybase,actbase)){
      rpolyend=runindex;
      lastgoodrunindex=runindex;
      uint32 tcount=0;
      uint32 othercount=0;
      uint32 runcount=0;
      char cbase;
      char dbase;
      for(; lastgoodrunindex>=actread.getLeftClipoff(); --lastgoodrunindex, ++runcount){
	cbase=actread.getBaseInSequence(lastgoodrunindex);
	CEBUG("###2 : " << runcount << " " << lastgoodrunindex << "\t" << cbase << " " << tcount << " " << othercount << endl);
	if(dptools::areBasesContained(polybase,cbase)){
	  tcount++;
	  if(othercount>0 && runcount>=mincount){
	    dbase=actread.getBaseInSequence(lastgoodrunindex+mincount);
	    if(dptools::areBasesContained(polybase,dbase)){
	      --othercount;
	    }
	  }
	}else if(tolower(cbase)!='n'){
	  othercount++;
	  if(othercount>maxbad) break;
	}
      }
      if(tcount>=mincount) {
	CEBUG("Found tcount\n");
	found=true;
	// get off non-t characters as far as possible
	if(lastgoodrunindex<actread.getLeftClipoff()) lastgoodrunindex++;
	while(lastgoodrunindex<rpolyend && !dptools::areBasesContained(polybase,actread.getBaseInSequence(lastgoodrunindex))) lastgoodrunindex++;
	lpolystart=lastgoodrunindex;
	break;
      }
      rpolyend=-1;
    }else{
      rpolyend=-1;
    }
  }

  CEBUG("LPOLYSTART: " << lpolystart << "\tRPOLYEND: " << rpolyend << endl);

  if(lpolystart >=0 && rpolyend != -1) {
    FUNCEND();
    return true;
  }

  FUNCEND();
  return false;
}
//#define CEBUG(bla)








/*************************************************************************
 *
 * clip poly-base at right end of read
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::clipPolyBaseAtEnd_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyBaseAtEnd_Pool(ReadPool & rpool, const std::string & logprefix)");

  cout << "Clipping dubious poly-base stretches at end of reads ... ";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      if((*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams().as_clip_3ppolybase_len){
	clipPolyBaseAtEnd_Read(actread,logprefix);
      }
    }
  }
}


void DataProcessing::clipPolyBaseAtEnd_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::clipPolyBaseAtEnd_Read(Read & actread, const std::string & logprefix)");

  assembly_parameters const & as_params= (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  CEBUG(actread.getName() << endl);

  uint32 mincount=as_params.as_clip_3ppolybase_len;
  if(mincount==0){
    MIRANOTIFY(Notify::FATAL, "-AS:c3ppmsl may not be 0");
  }
  if(actread.getLenClippedSeq() < mincount) return;

  uint32 maxbad=as_params.as_clip_3ppolybase_maxerrors;
  int32 grace=static_cast<int32>(as_params.as_clip_3ppolybase_maxgap);

  // first guess which base might be a polybase
  //
  // count occurences of bases in last (mincount+grace or mincount?) positions of read
  // the largest count >=30% of real bases (no 'n') wins

  Read::setCoutType(Read::AS_FASTA);
  CEBUG(actread << endl);
  Read::setCoutType(Read::AS_TEXTSHORT);
  CEBUG(actread << endl);

  DP_tmpau32_128['a']=0;
  DP_tmpau32_128['c']=0;
  DP_tmpau32_128['g']=0;
  DP_tmpau32_128['n']=0;
  DP_tmpau32_128['t']=0;

  int32 runindex=actread.getRightClipoff()-1;
  for(uint32 ri=0; ri<mincount && runindex>=actread.getLeftClipoff(); --runindex, ++ri) {
    char actbase=static_cast<char>(tolower(actread.getBaseInSequence(runindex)));
    ++DP_tmpau32_128[actbase];
  }

  CEBUG("CV: " << DP_tmpau32_128['a'] << " " << DP_tmpau32_128['c'] << " " << DP_tmpau32_128['g'] << " " << DP_tmpau32_128['t'] << endl);

  uint32 realbases=DP_tmpau32_128['a']+DP_tmpau32_128['c']+DP_tmpau32_128['g']+DP_tmpau32_128['t'];
  uint32 maxreal=max(DP_tmpau32_128['a'],max(DP_tmpau32_128['c'],max(DP_tmpau32_128['g'],DP_tmpau32_128['t'])));

  CEBUG("RB: " << realbases << "\tMR: " << maxreal << endl);

  char tentativepolybase='?';
  if(realbases>0 && 100*maxreal/realbases >= 30){
    CEBUG("MRThresh\n");
    for(uint32 testi=0; testi<4; ++testi){
      if(DP_tmpau32_128["acgt"[testi]]==maxreal){
	CEBUG("MRThreshHit\n");
	tentativepolybase="acgt"[testi];
	break;
      }
    }
  }

  // so, if a tentative polybase was found, try to find some clips and clip if found

  if(tentativepolybase!='?') {
    int32 lpolystart=-1;
    int32 rpolyend=-1;
    CEBUG("looking...\n");
    if(searchPolyBaseFrom3Prime(actread,tentativepolybase,lpolystart,rpolyend,mincount,maxbad,grace)){
      actread.setRMClipoff(lpolystart);
      CEBUG("setting rm " << lpolystart << endl);

      DP_logfout << logprefix << " poly-base " << tentativepolybase << " at end "
		 << actread.getName()
		 << "\tMask right: "
		 << actread.getRMClipoff()
		 << '\n';
    }
  }
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush();}
void DataProcessing::adaptorRightClip_Pool(ReadPool & rp, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::adaptorRightClip_Pool(ReadPool & rpool, const std::string & logprefix)");

  cout << "Searching for sequencing adaptors.\n";
  cout.flush();

  for(uint32 actid=0; actid < rp.size(); ++actid){
    Read & actread = rp[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      adaptorRightClip_Read(actread,logprefix);
    }
  }
}

void DataProcessing::adaptorRightClip_Read(Read & actread, const std::string & logprefix)
{
  FUNCSTART("void DataProcessing::adaptorRightClip_Read(Read & actread, const std::string & logprefix)");

  priv_EnsureAdapRegexes(actread.getReadGroupID());
  BUGIFTHROW(actread.getReadGroupID().getLibId()>=DP_adapres.size(),"Huh? no re lib " << actread.getReadGroupID().getLibId());

  priv_EnsureAdapSkims(actread.getReadGroupID());

  //assembly_parameters const & as_params= (*DP_miraparams_ptr)[actread.getSequencingType()].getAssemblyParams();
  auto oldrsclip=actread.getRSClipoff();
  auto newclip=-1;
  if(DP_adapskims[actread.getReadGroupID().getLibId()].skimptr==nullptr){
    cout << "Bah? nullptr for " << actread.getReadGroupID().getLibId() << "?" << endl;
    return;
  }
  readid_t ridadapfound=-1;
  if(DP_adapskims[actread.getReadGroupID().getLibId()].skimptr){
    newclip=DP_adapskims[actread.getReadGroupID().getLibId()].skimptr->findAdaptorRightClip(actread,9,ridadapfound,DP_threadid);
  }
  if(newclip>=0){
    ++DP_stats.cadapright;
    actread.setRSClipoff(newclip);
    DP_logfout << logprefix << " "
	       << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
	       << " adaptor: ";
    if(ridadapfound>=0) {
      DP_logfout << DP_adapskims[actread.getReadGroupID().getLibId()].poolptr->getRead(ridadapfound).getName() << " in ";
    }
    DP_logfout << actread.getName()
	       << " changed right clip from " << oldrsclip << " to " << newclip << "\n";
  }else{
    string seq(actread.getSeqAsChar());
    boost::to_upper(seq);

    boost::match_results<std::string::const_iterator> what;
    boost::match_flag_type flags = boost::match_default;
    std::string::const_iterator start, end;

    for(auto & msre : DP_adapres[actread.getReadGroupID().getLibId()].adapres){
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
	      ++DP_stats.cadaprightpartial;
	      actread.setRSClipoff(what.position());
	      DP_logfout << logprefix << " "
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


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void DataProcessing::seqMatchPhiX174_Read(Read & actread, const std::string & logprefix, bool filter)
{
  FUNCSTART("void DataProcessing::seqMatchPhiX174_Read(Read & actread, const std::string & logprefix)");

  BUGIFTHROW(!DP_px174hs_init,"Phi X 174 search structure not initialised.");

  auto numbaithits=DP_phix174hashstatistics.checkBaitHit(actread,DP_baiting_singlereadvhraparray,DP_baiting_tagmaskvector);
  if(numbaithits>10){
    ++DP_stats.cphix174;
    DP_logfout << logprefix << " "
	       << ReadGroupLib::getNameOfSequencingType(actread.getSequencingType())
	       << " phix174 in "
	       << actread.getName()
	       << " ... killed read\n";
    if(filter) actread.setRSClipoff(0);
  }
}


void DataProcessing::stdTreatmentPool_SingleThread(DataProcessing & dp, ReadPool & rpool, vector<uint8> * debrisreasonptr, std::string & logprefix, bool progress, int32 fromid, int32 toid)
{
  FUNCSTART("void DataProcessing::stdTreatmentPool_SingleThread(DataProcessing & dp, ReadPool & rpool, std::string & logprefix, bool progress, int32 fromid, int32 toid)");

  if(fromid<0) fromid=0;
  if(toid<0) toid=rpool.size();
  BUGIFTHROW(fromid>toid,"fromid>toid ?");
  BUGIFTHROW(toid>rpool.size(),"toid>rpool.size()?");

  unique_ptr<ProgressIndicator<int64> > pi;
  if(progress) pi=std::unique_ptr<ProgressIndicator<int64>>(new ProgressIndicator<int64>(fromid,toid));
  for(uint32 actid=fromid; actid < toid; ++actid){
    if(progress) pi->increaseprogress();
    Read & actread = rpool[actid];
    if(actread.hasValidData()
       && !(actread.isBackbone() || actread.isRail())){
      priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_SHORTONLOAD);

      auto & asp = (*(rpool.getMIRAParams()))[actread.getSequencingType()].getAssemblyParams();

      if(asp.as_search_phix174){
	dp.seqMatchPhiX174_Read(actread,logprefix,asp.as_filter_phix174);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_PHIX174);
      }
      if(asp.as_clip_badsolexaends && actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	dp.clipBadSolexaEnds_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_BADSOLEXAEND);
      }
      if(asp.as_clip_knownadaptorsright){
	dp.adaptorRightClip_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_KNOWNADAPTORRIGHT);
      }
      if(asp.as_clip_quality_minthreshold){
	dp.minimumQualityThreshold_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_QUALMINTHRESHOLD);
      }
      if(asp.as_clip_lowercase_front){
	dp.lowerCaseClippingFront_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_LOWERCASEFRONT);
      }
      if(asp.as_clip_lowercase_back){
	dp.lowerCaseClippingBack_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_LOWERCASEBACK);
      }
      if(asp.as_clip_quality){
	dp.qualClips_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_QUALCLIPS);
      }
      if(asp.as_clip_maskedbases){
	dp.maskClips_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MASKEDBASES);
      }
      bool mlc=asp.as_clip_ensureminimumleftclipoff;
      if(asp.as_clip_badstretchquality){
	if(mlc){
	  dp.maskClips_Read(actread,logprefix);
	  mlc=false;
	  priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MASKEDBASES);
	}
	dp.badSequenceSearch_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_BADSEQUENCESERACH);
      }
      if(asp.as_clip_3ppolybase){
	dp.clipPolyBaseAtEnd_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_POLYBASEATEND);
      }
      if(asp.as_clip_polyat){
	dp.clipPolyATAtEnds_Read(actread,logprefix);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_POLYAT);
      }
      if(mlc){
	dp.minimumLeftClip_Read(actread,logprefix,true,false,false);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MINLEFTCLIP);
      }
      if(asp.as_clip_ensureminimumrightclipoff){
	dp.minimumRightClip_Read(actread,logprefix,false,true,false);
	priv_stp_helperDebris(rpool,actread,actid,debrisreasonptr,Assembly::DEBRIS_CLIP_MINRIGHTCLIP);
      }

    }
  }
  if(progress) pi->finishAtOnce();
}

void DataProcessing::priv_stp_helperDebris(ReadPool & rpool, Read & actread, int32 rid, vector<uint8> * debrisreasonptr, uint8 reason)
{
  FUNCSTART("void DataProcessing::priv_stp_helperDebris(ReadPool & rpool, Read & actread, int32 rid, vector<uint8> * debrisreasonptr, uint8 reason)");

  if(debrisreasonptr == nullptr) return;
  if(debrisreasonptr->empty()) return;
  BUGIFTHROW(rid>=debrisreasonptr->size(),"rid " << rid << " >= debrisreasonptr->size() " << debrisreasonptr->size());
  if(actread.getLenClippedSeq() < (*(rpool.getMIRAParams()))[actread.getSequencingType()].getAssemblyParams().as_minimum_readlength){
    if((*debrisreasonptr)[rid]==0) (*debrisreasonptr)[rid]=reason;
  }
}


void DataProcessing::stdTreatmentPool_MultiThread(DataProcessing & dpcollector, vector<unique_ptr<DataProcessing>> & dpv, ReadPool & rpool, vector<uint8> * debrisreasonptr, std::string & logprefix, bool progress, int32 fromid, int32 toid)
{
  FUNCSTART("void DataProcessing::stdTreatmentPool_MultiThread(vector<DataProcessing> & dpv, ReadPool & rpool, std::string & logprefix, bool progress)");

  if(fromid<0) fromid=0;
  if(toid<0) toid=rpool.size();
  BUGIFTHROW(dpv.empty(),"dpv.empty() ?");
  BUGIFTHROW(fromid>toid,"fromid>toid ?");
  BUGIFTHROW(toid>rpool.size(),"toid>rpool.size()?");

  dpv[0]->priv_EnsurePhiX174Statistics();

  threadsharecontrol_t tsc;

  tsc.from=fromid;
  tsc.to=toid;
  tsc.todo=fromid;
  tsc.done=fromid;
  tsc.stepping=1000;

  uint32 numthreads=dpv.size();
  boost::thread_group workerthreads;
  for(uint32 ti=0; ti<numthreads;++ti){
    dpv[ti]->setThreadID(ti);
    workerthreads.create_thread(boost::bind(&DataProcessing::priv_stdTreatmentThread, ti, &tsc, &(*dpv[ti]), &rpool, debrisreasonptr, &logprefix));
  }

  ProgressIndicator<int64> pi(fromid,toid);
  while(tsc.done!=toid){
    pi.progress(tsc.done);
    sleep(1);
  }
  pi.finishAtOnce(cout);

  // they normally should all have exited at this point, but be nice and play by the rules
  workerthreads.join_all();

  // collect all stats
  for(auto & dpvp : dpv){
    dpcollector.DP_stats.cphix174+=dpvp->DP_stats.cphix174;
    dpcollector.DP_stats.cadapright+=dpvp->DP_stats.cadapright;
    dpcollector.DP_stats.cadaprightpartial+=dpvp->DP_stats.cadaprightpartial;
  }
}


void DataProcessing::priv_stdTreatmentThread(uint32 threadnum, threadsharecontrol_t * tscptr, DataProcessing * dpptr, ReadPool * rpoolptr, vector<uint8> * debrisreasonptr, std::string * logprefixptr)
{
  FUNCSTART("void DataProcessing::priv_stdTreatmentThread(uint32 threadnum, threadsharecontrol_t * tscptr, DataProcessing * dpptr, ReadPool * rpoolptr, std::string * logprefixptr)");

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
      stdTreatmentPool_SingleThread(*dpptr,*rpoolptr,debrisreasonptr,*logprefixptr,false,from,to);
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
