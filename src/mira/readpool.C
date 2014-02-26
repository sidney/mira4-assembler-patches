/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 1997-2000 by the German Cancer Research Center (Deutsches
 *   Krebsforschungszentrum, DKFZ Heidelberg) and Bastien Chevreux
 *   and Thomas Pfisterer
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


#include "mira/readpool.H"

#include "errorhandling/errorhandling.H"

#include "io/generalio.H"
#include "io/fasta.H"
#include "io/fastq-lh.H"
#include "io/phd.H"
#include "io/ncbiinfoxml.H"

#include "util/fileanddisk.H"
#include "util/progressindic.H"

#include "caf/caf.H"

#include "mira/gbf_parse.H"
#include "mira/gff_parse.H"
#include "mira/maf_parse.H"

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

KSEQ_INIT(gzFile, gzread)

using namespace std;

#define CEBUG(bla)

string ReadPool::RP_missingfastaqual_resolvemsg("name your files with a .fna postfix");

const base_quality_t ReadPool::RP_sxa2phredmap[256]= {
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 1, 1, 2, 2, 3,
  3, 4, 4, 5, 5, 6, 7, 8,
  9, 10, 10, 11, 12, 13, 14, 15,
  16, 17, 18, 19, 20, 21, 22, 23,
  24, 25, 26, 27, 28, 29, 30, 31,
  32, 33, 34, 35, 36, 37, 38, 39,
  40, 41, 42, 43, 44, 45, 46, 47,
  48, 49, 50, 51, 52, 53, 54, 55,
  56, 57, 58, 59, 60, 61, 62, 63,
  64, 65, 66, 67, 68, 69, 70, 71,
  72, 73, 74, 75, 76, 77, 78, 79,
  80, 81, 82, 83, 84, 85, 86, 87,
  88, 89, 90, 91, 92, 93, 94, 95,
  96, 97, 98, 99, 100, 101, 102, 103,
  104, 105, 106, 107, 108, 109, 110, 111,
  112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 124, 125, 126, 127,
  128, 129, 130, 131, 132, 133, 134, 135,
  136, 137, 138, 139, 140, 141, 142, 143,
  144, 145, 146, 147, 148, 149, 150, 151,
  152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 163, 164, 165, 166, 167,
  168, 169, 170, 171, 172, 173, 174, 175,
  176, 177, 178, 179, 180, 181, 182, 183,
  184, 185, 186, 187, 188, 189, 190, 191
};



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::ReadContainer const & ReadPool::ReadContainer::operator=(ReadPool::ReadContainer const & other)
{
  // 1 to 1 copy, including eventually free elements so that
  // operator[] has the same result on source and copy object
  if(this != &other){
    RC_thepool.clear();
    for(auto rptr : other.RC_poolrptr){
      RC_thepool.push_back(*rptr);
      RC_poolrptr.push_back(&(RC_thepool.back()));
    }
    RC_releasedidx=other.RC_releasedidx;
  }
  //cout << "Copied ReadContainer size " << size() << " and with " << getNumActiveReads() << " active reads\n";
  return *this;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::ReadPool(vector<MIRAParameters> * params)
{
  FUNCSTART("ReadPool::ReadPool()");

  REP_thepool3.clear();
  //REP_filenames=nullptr;
  REP_miraparams=params;

  REP_allownameindex=false;
  REP_nameindex.clear();

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::discard()
{
  FUNCSTART("ReadPool::discard()");

  //REP_miraparams=nullptr;
  REP_thepool3.clear();
  REP_nameindex.clear();

  //delete [] REP_filenames;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
ReadPool::~ReadPool()
{
  FUNCSTART("ReadPool::~ReadPool()");

  discard();

  FUNCEND();
}


//void ReadPool::setMIRAParameters(MIRAParameters * params)
//{
//  FUNCSTART("void ReadPool::setMIRAParameters(MIRAParameters * params);");
//
//  REP_miraparams=params;
//
//  FUNCEND();
//}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

size_t ReadPool::estimateMemoryUsage() const
{
  FUNCSTART("size_t ReadPool::estimateMemoryUsage()");

  size_t ret=sizeof(ReadPool);
  // TODO: estimateMemoryUsage for deque
  //ret+=estimateMemoryUsageOfContainer(REP_thepool2,false);

  FUNCEND();
  return ret;
}


/*************************************************************************
 *
 * makes template IDs for the reads in the pool
 *
 * returns
 *  true if there are usable templates (less templates than valid reads)
 *  false otherwise
 *
 *************************************************************************/

bool ReadPool::makeTemplateIDs(bool verbose)
{
  FUNCSTART("void ReadPool::makeTemplateIDs()");

// TODO: bla auf backbone und rails eingehen

  if(verbose){
    rpDateStamp();
    cout << endl;
  }

  vector<int32> tidcounter;
  vector<int32> tid_firstpartner;
  tidcounter.resize(size(), 0);
  tid_firstpartner.resize(size(), -1);

  // will we need to check template ends?
  // rationale: if no info about template ends was available in the
  //  ancillary data, the no_te_check will ensure that we'll then still use
  //  the template information
  bool notecheck=true;
  for(size_t i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).getTemplateSegment()!=0) {
      notecheck=false;
    }
  }


  typedef boost::unordered_map<std::string, int32> strintmap;
  strintmap tnmap;
  strintmap::iterator tnI;

  int32 acttid=0;
  int32 validreads=0;

  uint32 outlines_fatal=0;
  uint32 outlines_warn=0;
  uint32 ol_emptytn=0;
  uint32 ol_gt2reads=0;
  uint32 ol_unknownsegment=0;
  uint32 ol_unknownsegmentfixed=0;
  uint32 ol_samesegment=0;
  for(uint32 rid=0; rid<size();++rid){
      //cout << "acttid: " << acttid << "\t";
    if(getRead(rid).hasValidData()==false) continue;
    ++validreads;
    if(getRead(rid).getTemplate().empty()) {
      getRead(rid).setTemplateID(acttid);
      tidcounter[acttid]++;
      ++acttid;
      if(getRead(rid).getReadGroupID().hasTemplateInfo()){
	// Ooooops ... empty template but readgroupinfo?
	if(ol_emptytn<200){
	  ++outlines_fatal;
	  cout << "Read " << getRead(rid).getName() << " has template info given, but the internal template name is empty? This is fishy!" << endl;
	  if(++ol_emptytn==200){
	    cout << "More than 200 cases like the above, will not report more.\n";
	  }
	}
      }
    }else{
      tnI=tnmap.find(getRead(rid).getTemplate());
      if(tnI!=tnmap.end()){
	getRead(rid).setTemplateID(tnI->second);
	int32 firstpartner=tid_firstpartner[tnI->second];
	if(++tidcounter[tnI->second]==2){
	  getRead(rid).setTemplatePartnerID(firstpartner);
	  getRead(firstpartner).setTemplatePartnerID(rid);

	  // look out for unknown segments ("0"), try to fix if possible, else dump error message
	  uint8 numunknownseg=(getRead(rid).getTemplateSegment()==0)+(getRead(firstpartner).getTemplateSegment()==0);
	  if(numunknownseg){
	    bool couldfix=false;
	    if(numunknownseg==1){
	      Read * tofix=&getRead(rid);
	      Read * other=&getRead(firstpartner);
	      uint8 segmenttoset=0;
	      if(tofix->getTemplateSegment()>0) swap(tofix,other);
	      if(other->getTemplateSegment()==1){
		segmenttoset=255;
	      }else if(other->getTemplateSegment()==255){
		segmenttoset=1;
	      }
	      if(segmenttoset){
		couldfix=true;
		tofix->setTemplateSegment(segmenttoset);
		ol_unknownsegmentfixed++;
		if(ol_unknownsegmentfixed<200){
		  cout << "Read " << tofix->getName() << ": fixed unrecognised template segment.\n";
		}
		if(++ol_unknownsegmentfixed==200){
		  cout << "More than 200 cases like the above, will not report more.\n";
		}
	      }
	    }
	    if(!couldfix){
	      if(ol_unknownsegment<200){
		++outlines_fatal;
		if(getRead(rid).getTemplateSegment()==0){
		  cout << "DNA template " << tnI->first << ", read " << getRead(rid).getName() << ": template segment not recognised.\n";
		}
		if(getRead(firstpartner).getTemplateSegment()==0){
		  cout << "DNA template " << tnI->first << ", read " << getRead(firstpartner).getName() << ": template segment not recognised.\n";
		}
		if(++ol_unknownsegment==200){
		  cout << "More than 200 cases like the above, will not report more.\n";
		}
	      }
	    }
	  }

	  // the folowing is almost impossible when reading FASTA/FASTQ and other simple files,
	  //  but it may happen in CAF/MAF (and maybe SAM?)
	  if(getRead(rid).getTemplateSegment()==getRead(firstpartner).getTemplateSegment()){
	    if(ol_samesegment<200){
	      ++outlines_fatal;
	      cout << "Reads " << getRead(rid).getName() << " and " << getRead(firstpartner).getName() << " have the same template segments.\n";
	      if(++ol_samesegment==200){
		cout << "More than 200 cases like the above, will not report more.\n";
	      }
	    }
	  }

	}else{
	  // Ooooops ... not good: more than two reads for this template
	  getRead(rid).setTemplatePartnerID(-1);
	  getRead(firstpartner).setTemplatePartnerID(-1);
	  if(ol_gt2reads<200){
	    ++outlines_fatal;
	    cout << tidcounter[tnI->second] << " ";
	    cout << "DNA template " << tnI->first << " has more than two reads, template info not used. Read at fault: " << getRead(rid).getName() << endl;
	    if(++ol_gt2reads==200){
	      cout << "More than 200 cases like the above, will not report more.\n";
	    }
	  }
	}
      }else{
	getRead(rid).setTemplateID(acttid);

	tid_firstpartner[acttid]=rid;

	tidcounter[acttid]++;
	tnmap[getRead(rid).getTemplate()]=acttid;
	++acttid;
      }
    }
  }

  if(outlines_fatal>0 && (*REP_miraparams)[0].getNagAndWarnParams().nw_check_templateproblems){
    static string emsg="Problems found with the template data, see output log for more info. This points some serious problem either with read naming (like unrecognised read naming scheme) or broken template information, please fix your input files!\nOr switch off this warning with -NW:ctp=no (but you'll do this at own risk!)";
    if((*REP_miraparams)[0].getNagAndWarnParams().nw_check_templateproblems==NWSTOP) {
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << emsg << endl;
    }
  }

  if(verbose){
    cout << "Generated " << acttid << " unique DNA template ids for " << validreads << " valid reads.\n";
  }

  FUNCEND();
  return acttid!=validreads;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadPool::checkTemplateIDs(const string & errmsg)
{
  FUNCSTART("void ReadPool::checkTemplateIDs(string & errmsg)");
  for(uint32 ri=0;ri< REP_thepool3.size(); ++ri){
    if(REP_thepool3.getRead(ri).getTemplatePartnerID()>=0
       && REP_thepool3.getRead(ri).getTemplateID() != REP_thepool3.getRead(REP_thepool3.getRead(ri).getTemplatePartnerID()).getTemplateID()){
      cout << "Ouch, template problem for read " << ri << " " << REP_thepool3.getRead(ri).getName() << ", dumping readpool for debug\n";
      dumpAs(cout,Read::AS_TEXT,true);
      BUGIFTHROW(true, errmsg);
    }
  }
}

/*************************************************************************
 *
 * return true if no duplicate was found
 *
 *************************************************************************/

bool ReadPool::checkForDuplicateReadNames()
{
  FUNCSTART("bool ReadPool::checkForDuplicateReadNames(bool verbose)");

  typedef boost::unordered_set<std::string> strset;
  strset rnset;
  strset::iterator rnI;

  bool allok=true;
  uint32 outlines=0;
  for(uint32 i=0; i<size();i++){
    if(getRead(i).hasValidData()==false) continue;

    rnI=rnset.find(getRead(i).getName());
    if(rnI!=rnset.end()){
      allok=false;
      if(outlines<2000){
	cout << "Read " << *rnI << " is present more than once in the data set. Did you load a file twice in the manifest? Is a read present more than once in your file(s)?\n";
	if(++outlines==2000){
	  cout << "More than 2000 cases like the above, will not report more. Fix your input!\n";
	}
      }
    }else{
      rnset.insert(getRead(i).getName());
    }
  }

  if(outlines>0 && (*REP_miraparams)[0].getNagAndWarnParams().nw_check_duplicatereadnames){
    static string emsg="Some read names were found more than once (see log above). This usually hints to a serious problem with your input and should really, really be fixed. You can choose to ignore this error with '-NW:cdrn=no', but this will almost certainly lead to problems with result files (ACE and CAF for sure, maybe also SAM) and probably to other unexpected effects.";
    if((*REP_miraparams)[0].getNagAndWarnParams().nw_check_duplicatereadnames==NWSTOP){
      MIRANOTIFY(Notify::FATAL,emsg);
    }else{
      cout << "WARNING!\n" << emsg << endl;
    }
  }

  FUNCEND();
  return allok;
}





/*************************************************************************
 *
 * FASTA: optfilename2 is name of quality file
 * FASTQ: if optfilename2 not empty: has FASTQ qual offset as ASCII string
 *
 *
 *
 *************************************************************************/

size_t ReadPool::loadData_rgid(const string & filetype, const string & filename1, const string & optfilename2, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadData_rgid(const string & filetype, const string & filename1, const string & optfilename2, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))");

  if(filetype=="fastq"){
    base_quality_t fqqualoffset=0;
    if(!optfilename2.empty()){
      fqqualoffset=static_cast<base_quality_t>(atoi(optfilename2.c_str()));
    }
    loadDataFromFASTQ_rgid(filename1, fqqualoffset, rgid, countonly, callback);
  }else if(filetype=="fasta"){
    loadDataFromFASTA_rgid(filename1, rgid, true, optfilename2, countonly, callback);
  }else if(filetype=="fna" || filetype == "fastanoqual" || filetype=="fa"){ // FASTA, but no qual. fna to be used by users, fastanoqual internally
    loadDataFromFASTA_rgid(filename1, rgid, false, "", countonly, callback);
  }else if(filetype=="gbf"){
    loadDataFromGBF_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="gff3"){
    loadDataFromGFF3_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="fofnexp"){
    loadDataFromFOFNEXP_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="exp"){
    loadDataFromEXP_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="caf"){
    loadDataFromCAF_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="maf"){
    loadDataFromMAF_rgid(filename1, rgid, countonly, callback);
  }else if(filetype=="scf"){
    // scf is actually just a data directory
    cout << "Setting to " << filename1 << endl;
    const_cast<ReadGroupLib::ReadGroupID &>(rgid).setDataDir(filename1);
  }else if(filetype.empty()){
    MIRANOTIFY(Notify::FATAL,"While trying to load data for readgroup, an empty string was given as file type (should have been fasta, fastq, etc.pp).\nReadgroup:" << rgid << endl);
  }else{
    MIRANOTIFY(Notify::FATAL,"While trying to load data for readgroup, type " << filetype << " not known.");
  }

  FUNCEND();

  //TODO: should we return the number of leaded reads?
  return 0;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

size_t ReadPool::loadDataFromFOFNEXP_rgid(const string & fofn, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadDataFromFOFNEXP(const string & fofn, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))");

  string justfilename, justpath;
  splitFullPathAndFileName(fofn,justpath,justfilename);

  cout << "Loading file of filenames: " << fofn << endl;


  // Load the file of filenames
  ifstream fin;
  fin.open(fofn.c_str(), ios::in|ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << fofn);
  }

  if(fin.rdstate()){
    MIRANOTIFY(Notify::FATAL, "Failed to open file: " << fofn);
  }

  streamoff len_fofn=fin.tellg();
  if(len_fofn==1){
    MIRANOTIFY(Notify::FATAL, "Zero length file: " << fofn);
  }
  fin.seekg(0, ios::beg);

  vector<string> names;
  {
    string filename, dummy;
    while(GeneralIO::readKeyValue(fin, filename, dummy)){
      names.push_back(filename);
    }
  }
  fin.close();

  bool stopprocessing=false;
  {
    //stringhash_t M;
    //pair<stringhash_t::const_iterator, stringhash_t::const_iterator> p;

    typedef boost::unordered_set<std::string> strset;
    strset namemap;
    strset::iterator nI;

    for(uint32 i=0; i< names.size(); i++){
      nI=namemap.find(names[i]);
      if(nI!=namemap.end()){
	cout << "WARNING: file " << names[i] << " is present more than once in your file of filenames." << endl;
	stopprocessing=true;
      }else{
	namemap.insert(names[i]);
      }
    }
  }
  cout << "done." << endl;

  if(stopprocessing){
    MIRANOTIFY(Notify::FATAL, "Some reads lead to unrecoverable error: duplicate names. Aborting, see log above for further information.");
  }

  if(countonly) return names.size();

  rpDateStamp();

  size_t no_files_ok=0;
  {
    ProgressIndicator<int64> P(0, names.size()-1);

    cout << "Loading EXP files: " << endl;

    //string completename;
    for(uint32 namei=0; namei<names.size(); ++namei){
      P.progress(namei);
      if(names[namei].empty()) continue;
      try{
	Read & actread = getRead(provideEmptyRead());
	actread.setReadGroupID(rgid);

	actread.loadDataFromEXP(names[namei],justpath);
	actread.transferSVTagsToClip(20,60);

	no_files_ok++;

	if(callback!=nullptr) {
	  (*callback)(*this);
	}
      }
      catch(Flow){
      }
      catch(Notify n){
	stopprocessing=true;
	n.handleError(THISFUNC);
      }
    }
    P.finishAtOnce();

    cout << "\nDone.\n";

    cout << "There haven been " << names.size() << " files given, " ;
    cout << no_files_ok << " of which have loaded ok.\n";
  }

  rpDateStamp();
  cout << endl;

  FUNCEND();

  return names.size();;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

size_t ReadPool::loadDataFromEXP_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadEXP(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))");

  uint32 no_files_ok=0;
  try{
    Read & actread = getRead(provideEmptyRead());
    actread.setReadGroupID(rgid);

    actread.loadDataFromEXP(filename,"");
    actread.transferSVTagsToClip(20,60);

    no_files_ok++;

    if(callback!=nullptr) {
      (*callback)(*this);
    }
  }
  catch(Flow){
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }

  FUNCEND();

  return no_files_ok;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/
size_t ReadPool::loadDataFromFASTQ_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadDataFromFASTQ_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))");
  return loadDataFromFASTQ_rgid(filename,
				-1,
				rgid,
				countonly,
				callback);
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

size_t ReadPool::loadDataFromFASTQ_rgid(const string & filename, base_quality_t fastqoffset, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadDataFromFASTQ_rgid(const string & filename, int8 fastqoffset, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))");

  BUGIFTHROW(callback!=nullptr && fastqoffset<33,"Cannot handle automatic guessing of FASTQ type with callback after every read, sorry. callback: " << callback << "fastqoffset: " << static_cast<uint16>(fastqoffset));

  size_t numseqsloaded=0;

  bool fatalloaderror=false;
  bool qualerror=false;

  streamsize fsize;
  {
    ifstream fin;
    fin.open(filename.c_str(), ios::in|ios::ate);
    if(!fin){
      MIRANOTIFY(Notify::FATAL,"Could not open FASTQ file '" << filename << "'. Is it present? Is it readable? Did you want to load your data in another format?");
    }
    fsize=fin.tellg();
  }

  gzFile fp;
  kseq_t *seq;

  rpDateStamp();

  base_quality_t minqualseen=255;

  fp = gzopen(filename.c_str(), "r");
  if(fp==Z_NULL) {
    MIRANOTIFY(Notify::FATAL,"Could not open FASTQ file '" << filename << "' though it was possible just moments ago? Was it deleted?");
  }

  seq = kseq_init(fp);

  vector<base_quality_t> bq;
  bq.reserve(1000);

  size_t oldrpsize=size();
  base_quality_t seenminqual=255;
  base_quality_t seenmaxqual=0;
  uint32 longestread=0;


  cout << "Loading data from FASTQ file: " << filename << endl;
  cout << "(sorry, no progress indicator for that, possible only with zlib >=1.34)" << endl;

  //ProgressIndicator<int64> P(0, fsize,10000);

  int l;
  while ((l = kseq_read(seq)) >= 0) {
    //P.progress(gzoffset(fp));
    ++numseqsloaded;
    if(countonly) continue;

    Read & actread = getRead(provideEmptyRead());
    actread.setReadGroupID(rgid);

    bool maybesolexa=(rgid.getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA);

    if(rgid.getSequencingType() == ReadGroupLib::SEQTYPE_TEXT
       && seq->comment.l>0){
      uint32 numcolons=0;
      char * sptr=seq->comment.s;
      for(; *sptr; ++sptr){
	if(*sptr==':') ++numcolons;
      }
      if(numcolons==3){
	sptr=seq->comment.s;
	if(sptr[1]==':'
	   && (sptr[2]=='Y' || sptr[2]=='N')
	   && sptr[3]==':') {
	  maybesolexa=true;
	}
      }
    }

    if(maybesolexa){
      actread.disallowAdjustments();
    }

    if(maybesolexa){
      string tmpname(seq->name.s);
      string::size_type bpos = tmpname.rfind("/");
      //cout << "tmpname: " << tmpname << endl;
      if (bpos == string::npos && seq->comment.l>0) {
	//cout << "No / for " << tmpname << " ... need to make one:" << seq->comment.s << endl;
	char * colonptr=seq->comment.s;
	for(; *colonptr!=0; ++colonptr){
	  if(*colonptr==':') break;
	}
	if(*colonptr){
	  tmpname+='/';
	  colonptr=seq->comment.s;
	  while(*colonptr!=':') {
	    tmpname+=*colonptr;
	    ++colonptr;
	  }
	}
      }
      actread.setName(tmpname);
    }else{
      actread.setName(seq->name.s);
    }

    if(actread.getName().empty()){
      cout << "Ouch, there's a read without a name? This is illegal. The sequence\n  "
	   << seq->seq.s
	   << "\nmust have a name!\n";
      fatalloaderror=true;
    }

    if(seq->seq.l==0){
      actread.setValidData(false);
    }else{
      actread.setSequenceFromString(seq->seq.s);
      longestread=min(longestread,actread.getLenSeq());
      bq.clear();
      if(seq->qual.l){
	if(seq->qual.l != seq->seq.l){
	  cout << actread.getName()
	       << ": different number of quality values than bases?\n";
	  qualerror=true;
	}else{
	  const uint8 * qi = reinterpret_cast<const uint8 *>(seq->qual.s);
	  bool qualok=true;
	  for(;*qi; qi++) {
	    if(*qi<33 || *qi>164){
	      cout << "Read " << actread.getName() << ": invalid quality " << static_cast<uint16>(*qi) << '\n';
	      qualok=false;
	    }
	    bq.push_back(*qi);
	    if(*qi<seenminqual) seenminqual=*qi;
	    if(*qi>seenmaxqual) seenmaxqual=*qi;
	  }
	  if(qualok) {
	    actread.setQualities(bq);
	  }else{
	    qualerror=true;
	  }
	}
      }else{
	if(fastqoffset<33){
	  // ooops, trying to guess automatically ... not good if there's no sequence
	  // most probable nowadays: Sanger style FASTQ
	  bq.resize(seq->seq.l,rgid.getDefaultQual()+33);
	}else{
	  bq.resize(seq->seq.l,rgid.getDefaultQual()+fastqoffset);
	}
	actread.setQualities(bq);
	actread.setQualityFlag(rgid.getDefaultQual()>0);
      }
    }

    if(callback!=nullptr) {
      // if we have a callback, the translation of ascii values to base qualities
      //  must be done now ... no guessing what the file might be.
      for(auto & qv : const_cast<vector<base_quality_t> &>(actread.getQualities())){
	qv-=fastqoffset;
      }

      (*callback)(*this);
    }
  }
  //P.finishAtOnce();
  cout << "\n";

  if(l!=-1){
    cout << "Whoooops, something seems fishy with the last sequence loaded, the FASTQ parser returned " << l << " instead of the expected -1.\n";
    if(numseqsloaded>0 && !countonly && size()>0){
      cout << "\nLast read which seemed OK: " << getRead(size()-1).getName() << endl;
      MIRANOTIFY(Notify::FATAL,"FASTQ seems broken, see log above.\n");
    }
  }

  kseq_destroy(seq);
  gzclose(fp);

  if(qualerror){
    MIRANOTIFY(Notify::FATAL,"Unrecoverable error while loading data from FASTQ (see output above) Fix your input please.");
  }

  cout << "\nDone.\n";

  cout << "Loaded " << numseqsloaded << " reads, " ;
  //cout << num_reads_qual_ok << " of which have quality accounted for.\n";

  rpDateStamp();

  if(fatalloaderror) {
    MIRANOTIFY(Notify::FATAL, "Fatal error encountered during load of data (see log), aborting.\n") ;
  }

  if(callback==nullptr) {
    cout << "Looking at FASTQ type ... ";

    base_quality_t qualcorrector=0;
    bool needoldsxamapping=false;

    if(fastqoffset<33){
      // longestread criterion catches everything non-Solexa (including 454, or contigs) as Illumina
      //  switched to Sanger style before attaining 200bp
      if(seenminqual<59 || seenmaxqual <= 73 || seenminqual==seenmaxqual || longestread>=200){
	cout << "guessing FASTQ-33 (Sanger)\n";
	qualcorrector=33;
      }else if(seenminqual>=59 && seenminqual<64){
	cout << "guessing FASTQ-59 (old Solexa)\n";
	needoldsxamapping=false;
      }else{
	cout << "guessing FASTQ-64 (Illumina)\n";
	qualcorrector=64;
      }
    }else{
      cout << "told it to be FASTQ-" << static_cast<uint16>(fastqoffset) << '\n';
      qualcorrector=fastqoffset;
    }

    cout << "Running quality values adaptation ... "; cout.flush();
    if(needoldsxamapping){
      for(size_t rpi=oldrpsize; rpi<size(); ++rpi){
	for(auto & qv : const_cast<vector<base_quality_t> &>(getRead(rpi).getQualities())){
	  qv=RP_sxa2phredmap[qv];
	}
      }
    }else{
      for(size_t rpi=oldrpsize; rpi<size(); ++rpi){
	for(auto & qv : const_cast<vector<base_quality_t> &>(getRead(rpi).getQualities())){
	  qv-=qualcorrector;
	}
      }
    }
    cout << "done.\n";
  }

  FUNCEND();

  return numseqsloaded;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

size_t ReadPool::loadDataFromFASTA_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, const bool wantsqualfiletoexist, const string & qualfilename, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("size_t ReadPool::loadDataFromFASTA_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, const bool wantsqualfiletoexist, const string & qualfilename, void (*callback)(ReadPool &))");

  size_t numseqsloaded=0;

  ifstream fin;

  bool callcallbackforeachread=true;
  bool hasqualfile=false;
  if(qualfilename.size()!=0){
    fin.open(qualfilename.c_str(), ios::in|ios::ate);
    if(!fin){
      cout << "Could not find FASTA quality file " << qualfilename.c_str();
      if(wantsqualfiletoexist){
	cout << ", aborting. If you want to work without qualities, " << RP_missingfastaqual_resolvemsg;
	MIRANOTIFY(Notify::FATAL, "File not found: " << qualfilename);
      }else{
	cout << ", using default values for these reads.\n";
      }
    }else{
      hasqualfile=true;
      callcallbackforeachread=false;
    }
    if(!fin.tellg() && wantsqualfiletoexist){
      MIRANOTIFY(Notify::FATAL, "FASTA quality file " << qualfilename << " has zero length? Seems fishy.");
    }
    fin.close();
    fin.clear();
  }else{
    if(wantsqualfiletoexist){
      MIRANOTIFY(Notify::FATAL, "FASTA quality file expected to exist, but no quality filename given???");
    }
  }

  fin.open(filename.c_str(), ios::in|ios::ate);
  if(!fin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }
  if(!fin.tellg()){
    MIRANOTIFY(Notify::FATAL, "Zero length file: "  << filename);
  }

  bool fatalloaderror=false;

  // these two are for the automatic conversion of Solexa base scores
  //  to quality values
  bool sxa_foundnegativevalue=false;

  streamsize fsize=fin.tellg();

  size_t formerpoolsize=size();

  // qualset looks whether the loaded fastas had quality values
  vector<bool> qualset(size(),true);

  rpDateStamp();

  numseqsloaded=0;
  {
    cout << "Loading data from FASTA file:\n";
    fin.clear();
    fin.seekg(0, ios::beg);
    ProgressIndicator<streamsize> P(0, fsize-1,1000);
    FASTA thefasta;
    while(1){
      thefasta.loadNextSeq(fin);
      if(P.delaytrigger()) P.progress(fin.tellg());
      if(thefasta.testIfEmpty()) {
	// no more sequences.
	break;
      }
      ++numseqsloaded;
      if(countonly) continue;
      Read & actread = getRead(provideEmptyRead());
      actread.setReadGroupID(rgid);
      actread.setName(thefasta.getSeqName());

      if(rgid.getSequencingType() == ReadGroupLib::SEQTYPE_SOLEXA){
	actread.disallowAdjustments();
      }

      if(thefasta.getSequence().empty()){
	actread.setValidData(false);
	cout << "\nWarning: " << actread.getName() << " has no bases?! This usually points at some error in the processing of data before it arrives to MIRA.\n";
      }else{
	actread.setSequenceFromString(thefasta.getSequence());
      }


      if(REP_miraparams!=nullptr) {
	// standard quality at first, may be overwritten later in load stage
	base_quality_t bdq=rgid.getDefaultQual();

	if(actread.hasValidData()){
	  actread.setQualities(bdq);
	  actread.setQualityFlag(false);
	}
      }

      if(callcallbackforeachread && callback!=nullptr) {
	(*callback)(*this);
      }
    }
    P.finishAtOnce();
    cout << "\n";
  }

  fin.close();

  rpDateStamp();

  if(countonly) {
    return numseqsloaded;
  }

  //stringhash_t M;
  //pair<stringhash_t::const_iterator, stringhash_t::const_iterator> p;

  typedef boost::unordered_map<std::string, int32> strintmap;
  strintmap rnmap;
  strintmap::iterator rnI;

  // make a quick hash lookup of read names
  {
    bool haserror=false;
    for(uint32 i=0; i<size();i++){
      if(getRead(i).hasValidData()==false) continue;
      if(getRead(i).getName().size()==0) continue;
      rnI=rnmap.find(getRead(i).getName());
      if(rnI!=rnmap.end()){
	//cout << "uh oh ... double?";
	haserror=true;
	rnI->second+=1;
      }else{
	rnmap[getRead(i).getName()]=1;;
      }
    }

    if(haserror){
      for(rnI=rnmap.begin(); rnI!=rnmap.end(); ++rnI){
	if(rnI->second > 1){
	  cout << "Error: read name " << rnI->first << " present " << rnI->second << " times in readpool!\n";
	}
      }
      MIRANOTIFY(Notify::FATAL, "Read names not unique (either in this file or together with files loaded earlier): " << filename);

    }

    // now re-fill the rnmap with the read-ids
    rnmap.clear();
    cout << "rnm size: " << rnmap.size() << endl;
    for(uint32 i=0; i<size();i++){
      if(getRead(i).hasValidData()==false) continue;
      if(getRead(i).getName().size()==0) continue;
      rnmap[getRead(i).getName()]=i;
    }

    //cout << "\n---\n";
    //rnI=rnmap.begin();
    //for(; rnI != rnmap.end(); ++rnI){
    //  cout << "rnif" << rnI->first << "\trnis: " << rnI->second << endl;
    //}
    //cout << "---\n";
  }

  //dumpAsMAF(cout);

  size_t oldqualsetsize=qualset.size();

  qualset.resize(size(),false);
  int32 num_reads_qual_ok=0;
  if(hasqualfile && size()) {
    FASTA thefasta;

    fin.clear();
    fin.open(qualfilename.c_str(), ios::in|ios::ate);

    if(!fin){
      cout << "Could not find FASTA quality file " << qualfilename.c_str();
      if(wantsqualfiletoexist){
	cout << ", aborting. If you want to work without qualities, " << RP_missingfastaqual_resolvemsg;
	fatalloaderror=true;
      }else{
	cout << ", using default values for these reads.\n";
      }
    }
    if(!fin){
      cout << "Could not find FASTA quality file " << qualfilename.c_str() << " although I found it before loading the FASTA data???\nStrange ... please check what could have happened.\n";
      fatalloaderror=true;
    } else {
      cout << "Loading quality data from FASTA quality file " << qualfilename << ":\n";

      fsize=fin.tellg();
      if(fsize==0) fsize=1;
      fin.seekg(0, ios::beg);

      ProgressIndicator<streamsize> P(0, fsize-1,1000);

      while(1){
	thefasta.loadNextINTSeq(fin,255);
	if(P.delaytrigger()) P.progress(fin.tellg());
	if(thefasta.testIfEmpty()) break;
	try{
	  rnI=rnmap.find(thefasta.getQualName());
	  if(rnI!=rnmap.end()){
	    //cout << "Wanna set rnif" << rnI->first << "\trnis: " << rnI->second << endl;
	    if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA || rgid.getSequencingType()==ReadGroupLib::SEQTYPE_TEXT){
	      if(!sxa_foundnegativevalue) {
		const vector<int32> & svalues=thefasta.getINTValues();
		vector<int32>::const_iterator sI=svalues.begin();
		if(!sxa_foundnegativevalue){
		  for(;sI!=svalues.end(); sI++){
		    if(*sI<0) {
		      sxa_foundnegativevalue=true;
		    }
		  }
		  if(sxa_foundnegativevalue){
		    // first time a negative value has been encountered. I.e.,
		    //  reads loaded previously in this file were also old solexa
		    //  score. Need to correct those!
		    cout << "Whooops ... FASTA quality values <0? That's the very old Solexa scoring scheme.\n";
		    for(size_t rpi=oldqualsetsize; rpi<qualset.size(); ++rpi){
		      if(qualset[rpi]){
			vector<base_quality_t> & qv=const_cast<vector<base_quality_t> &>(getRead(rpi).getQualities());
			vector<base_quality_t>::iterator qvI=qv.begin();
			for(; qvI != qv.end(); ++qvI){
			  *qvI=RP_sxa2phredmap[*qvI+64];
			}
		      }
		    }
		  }
		}
	      }
	      if(sxa_foundnegativevalue){
		getRead(rnI->second).setQualities(0);
		vector<base_quality_t> & qv=const_cast<vector<base_quality_t> &>(getRead(rnI->second).getQualities());
		vector<base_quality_t>::iterator qvI=qv.begin();
		const vector<int32> & iv=thefasta.getINTValues();
		vector<int32>::const_iterator ivI=iv.begin();
		for(; qvI != qv.end(); ++qvI, ++ivI){
		  if(ivI==iv.end()){
		    MIRANOTIFY(Notify::FATAL,"Less values in quality for read " << getRead(rnI->second).getName() << " than there were bases?");
		  }
		  *qvI=RP_sxa2phredmap[*ivI+64];
		}
	      }else{
		getRead(rnI->second).setQualities(thefasta.getQualities());
	      }
	    }else{
	      getRead(rnI->second).setQualities(thefasta.getQualities());
	    }
	    qualset[rnI->second]=true;
	    num_reads_qual_ok++;
	  }else{
	    // TODO: make a WARNINGS file
	    cout << "Warning: " << thefasta.getQualName() << " has quality values, but was not present in sequence file?!\n";
	  }
	}
	catch(Notify n){
	  if(n.gravity==Notify::FATAL) fatalloaderror=true;
	  n.setGravity(Notify::WARNING);
	  n.handleError(THISFUNC);
	}
      }
      P.finishAtOnce();

      cout << "\n";

      //if(seqtype==ReadGroupLib::SEQTYPE_SOLEXA && !sxa_foundnegativevalue){
      //	cerr << "No negative values found in FASTA quality file of Solexa data??\n";
      //}
    }
  }else{
    cout << "No FASTA quality file given, using default qualities for all reads just loaded." << endl;
  }

  rpDateStamp();

  if(REP_miraparams!=nullptr) {
    for(size_t i=formerpoolsize; i<qualset.size(); i++){
      if(qualset[i]==false
	 && getRead(i).hasValidData()
	 && !(getRead(i).isBackbone() || getRead(i).isRail())) {
	if(hasqualfile) cout << getRead(i).getName() << " has no valid qualities, using default.\n";
      }
    }
  }


  cout << "\nDone.\n";

  cout << "Loaded " << numseqsloaded << " reads with "
       << num_reads_qual_ok << " reads having quality accounted for.\n";

  if(fatalloaderror) {
    MIRANOTIFY(Notify::FATAL, "Fatal error encountered during load of data (see log), aborting.\n") ;
  }

  if(!callcallbackforeachread && callback!=nullptr) {
    (*callback)(*this);
  }

  FUNCEND();
  return numseqsloaded;
}


/*************************************************************************
 *
 * GBF is one of the formats which may contain several sequences
 *  the object already loads them completely into memory
 *
 *************************************************************************/

void ReadPool::loadDataFromGBF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("void ReadPool::loadDataFromGBF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly=false, void (*callback)(ReadPool &))");

  GBF thegbf;
  thegbf.load(filename);
  thegbf.transferGeneInfoToCDSInfo();

  vector<multitag_t> tag2mt;
  for(uint32 si=0; si<thegbf.getNumSequences(); si++){
    Read & actread=getRead(provideEmptyRead());
    actread.setReadGroupID(rgid);

    //cout << "Read GBF " << i << endl;
    //cout << "\tName   : " << thegbf.getSequenceName(i) << endl;
    //cout << "\tLenseq : " << thegbf.getSequence(i).size() << endl;
    //cout << "\tNumtags: " << thegbf.getTags(i).size() << endl;

    actread.setName(thegbf.getSequenceName(si));
    actread.setSequenceFromString(thegbf.getSequence(si));
    actread.setTags(thegbf.getTags(si));
    actread.setQualities(rgid.getDefaultQual());

    if(callback!=nullptr) {
      (*callback)(*this);
    }
  }

  FUNCEND();
}

/*************************************************************************
 *
 * GFF3 is one of the formats which may contain either no sequence, one sequence or several sequences
 * If sequences are present, the object already loads them completely into memory
 *
 *************************************************************************/

void ReadPool::loadDataFromGFF3_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("void ReadPool::loadDataFromGFF3_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &));");

  GFFParse thegff(this);
  thegff.loadFile(filename);

  for(uint32 i=0; i<thegff.getNumSequences(); i++){
    if(thegff.getSequence(i).size()==0 && !thegff.getTags(i).empty()){
      // special case: empty sequence but tags means: the GFF had descriptions
      //  but no sequence for it. One of the already loaded reads however has a name
      //  matching, so simply add the tags to that read
      auto rid=getReadIndex(thegff.getSequenceName(i));
      BUGIFTHROW(rid<0,"rid<0 ? should not happen here");
      Read & actread=getRead(rid);
      for(auto & ttag : thegff.getTags(i)){
	actread.addTagO(ttag);
      }
    }else{
      Read & actread=getRead(provideEmptyRead());
      actread.setReadGroupID(rgid);

      string minft_strainname;
      string minft_seqtypename;
      string minft_machinetype;
      int8 minft_tplacementcode;
      bool minft_isbb;
      bool minft_israil;
      bool minft_isCER;

      // ugly, but we know what we do and do not want to copy around the tags
      // unnecessarily
      vector<multitag_t> & gff3tags=const_cast<vector<multitag_t> & >(thegff.getTags(i));

      if(Read::extractMINFTagInfo(gff3tags,
				  thegff.getSequenceName(i),
				  minft_strainname,
				  minft_seqtypename,
				  minft_machinetype,
				  minft_tplacementcode,
				  minft_isbb,
				  minft_israil,
				  minft_isCER)){

	string dummy_empty;

	uint8 st=ReadGroupLib::stringToSeqType(minft_seqtypename);
	if(st==ReadGroupLib::SEQTYPE_END && !minft_seqtypename.empty()){
	  MIRANOTIFY(Notify::FATAL,"File " << filename << " has invalid sequencing type '" << minft_seqtypename << "' in MIRA MIT2 tag.\n");
	}

	ReadGroupLib::ReadGroupID rgidt=ReadGroupLib::searchExactRGMatch(
	  dummy_empty,
	  st,
	  -1,
	  -1,
	  minft_tplacementcode,
	  minft_strainname,
	  minft_isbb,
	  minft_israil,
	  minft_isCER,
	  dummy_empty,
	  dummy_empty,
	  dummy_empty);

	if(rgidt.isDefaultNonValidReadGroupID()){
	  rgidt=ReadGroupLib::newReadGroup();
	  rgidt.setGroupName(dummy_empty);
	  rgidt.setSequencingType(st);
	  rgidt.setStrainName(minft_strainname);
	  rgidt.setBackbone(minft_isbb);
	  rgidt.setRail(minft_israil);
	  rgidt.setCoverageEquivalentRead(minft_isCER);
	  rgidt.setMachineType(minft_machinetype);
	}
      }

      actread.setName(thegff.getSequenceName(i));
      actread.setSequenceFromString(thegff.getSequence(i));
      actread.setTags(thegff.getTags(i));
      actread.setQualities(rgid.getDefaultQual());
    }

    if(callback!=nullptr) {
      (*callback)(*this);
    }
  }

  allowNameIndex(false);

  FUNCEND();
}




/*************************************************************************
 *
 * MAF is one of the formats which may contain contigs and/or reads
 * When loaded via the readpool, all reads are loaded, but no contigs
 *  created.
 * E.g.: fle contains 1 contig (2 reads) and one 1 read without contig
 *  -> 3 reads are added to readpool
 *
 *************************************************************************/

void ReadPool::loadDataFromMAF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("void ReadPool::loadDataFromMAF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly=false, void (*callback)(ReadPool &))");

  MAFParse mafp(this, nullptr, REP_miraparams);
  vector<uint32> dummy;
  mafp.load(filename,
	    rgid.getSequencingType(),
	    1,
	    dummy,
	    false,
	    nullptr,
	    callback
	  );

  FUNCEND();
}


/*************************************************************************
 *
 * CAF is one of the formats which may contain contigs and/or reads
 * When loaded via the readpool, all reads are loaded, but no contigs
 *  created.
 * E.g.: fle contains 1 contig (2 reads) and one 1 read without contig
 *  -> 3 reads are added to readpool
 *
 *************************************************************************/

void ReadPool::loadDataFromCAF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly, void (*callback)(ReadPool &))
{
  FUNCSTART("void ReadPool::loadDataFromCAF_rgid(const string & filename, const ReadGroupLib::ReadGroupID rgid, bool countonly=false, void (*callback)(ReadPool &))");

  CAF cafp(this, nullptr, REP_miraparams);
  vector<uint32> dummy;
  cafp.load(filename,
	    rgid.getSequencingType(),
	    1,
	    dummy,
	    false,
	    nullptr,
	    callback
	  );

  FUNCEND();
}




/*************************************************************************
 *
 * loads names of reads from external file and deletes them from pool
 * if invertselection true, then delets those not in the file
 *
 * BEWARE SIDE EFFECT: deletes all reads from pool that have invalid data
 *
 * the data must be in a key (a value in one line)
 * key can be: read name, or filename of exp read or file name of caf read
 *   (may not contain spaces, sorry)
 *
 * line with # as first nonwhitespace character are comments and read over
 *
 *************************************************************************/

void ReadPool::deleteReadsByName(const string & nfile, bool invertselection)
{
  FUNCSTART("void ReadPool::InvalidateReadsByName(const string & nfile)");

  BUGIFTHROW(true,"this needs to be adapted to indirection");

//  //stringhash_t M;
//  typedef boost::unordered_map<std::string, uint32> strintmap;
//  strintmap rnmap;
//  strintmap::iterator rnI;
//
//  for(uint32 i=0; i<size();i++){
//    if(!REP_thepool[i].getName().empty()) {
//      rnmap[REP_thepool[i].getName()]=i;
//    }
//  }
//
//  ifstream fin;
//  fin.open(nfile.c_str(), ios::in|ios::ate);
//  if(!fin){
//    MIRANOTIFY(Notify::FATAL, "File not found: " << nfile);
//  }
//  fin.seekg(0, ios::beg);
//
//  if(invertselection){
//    for(uint32 i=0; i<size();i++){
//      REP_thepool[i].setValidData(false);
//    }
//  }
//
//  string readname, dummy;
//  while(GeneralIO::readKeyValue(fin, readname,dummy)){
//    rnI=rnmap.find(readname);
//    if(rnI!=rnmap.end()) {
//      if(invertselection){
//	REP_thepool[rnI->second].setValidData(true);
//      }else{
//	REP_thepool[rnI->second].setValidData(false);
//      }
//    }
//  }
//  fin.close();
//
//  if(size()){
//    for(int32 i=size()-1; i >= 0; i--){
//      if(!REP_thepool[i].hasValidData()){
//	vector<Read>::iterator I=REP_thepool.begin();
//	advance(I,i);
//	REP_thepool.erase(I);
//      }
//    }
//  }

  FUNCEND();
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::mergeXMLTraceInfo(const string & xmlfile)
{
  FUNCSTART("void ReadPool::mergeXMLTraceInfo(const string & filename)");

  rpDateStamp();

  cout << "Merging data from XML trace info file " << xmlfile << " ...";
  cout.flush();

  //make this called!

  string id454="454";
  uint32 numfound=0;

  NCBIInfoXML nix;

  list<NCBIInfoXML::ncbitraceelements_t> traces;

  try{
    nix.readXMLFile(xmlfile, traces);
  }
  catch(Notify n) {
    cout << "\n\n\nMIRA tried to load a XML TRACEINFO file containing ancillary data, but failed.\n"
      "Loading ancillary data when using FASTA files as input is\n"
      "really,\n"
      "        really,\n"
      "                REALLY encouraged, and therefore MIRA sets this as default.\n"
      "\nHowever, if you are really sure that you do not want to load ancillary data\n"
      "in TRACEINFO files, you can switch it off.\n"
      "Either use '<technology>_SETTINGS -LR:mxti=no' (e.g. SANGER_SETTING -LR:mxti=no),\n"
      "or use the '-notraceinfo' quickswitch to kill loading traceinfo files for all\n"
      "types of sequencing technologies. (place it after -fasta and -job quickswitches)\n\n\n";
    n.handleError(THISFUNC);
  }


  cout << "Num reads: " << traces.size() << endl;

  typedef boost::unordered_map<std::string, int32> strintmap;
  strintmap rnmap;
  strintmap::iterator rnI;

  cout << "Building hash table ... "; cout.flush();

  for(uint32 i=0; i<size();i++){
    if(!getRead(i).getName().empty()) {
      rnmap[getRead(i).getName()]=i;
    }
  }
  cout << "done." << endl;

  string acttracename;
  list<uint32>::const_iterator E;
  list<string>::const_iterator ECD;

  list<NCBIInfoXML::ncbitraceelements_t>::const_iterator T=traces.begin();

  for(;T!=traces.end(); T++) {
    rnI=rnmap.end();

    E=T->elements.begin();
    ECD=T->elements_cdata.begin();
    bool found=false;
    for(;!found && E!=T->elements.end(); E++, ECD++) {
      if((*E == NCBIInfoXML::NCBIXML_TRACE_NAME
	  || *E == NCBIInfoXML::NCBIXML_TI)
	 && !ECD->empty()){
	if(*E == NCBIInfoXML::NCBIXML_TRACE_NAME){
	  acttracename=*ECD;
	} else if(*E == NCBIInfoXML::NCBIXML_TI){
	  acttracename="gnl|ti|"+*ECD;
	}
	rnI=rnmap.find(acttracename);

	if(rnI!=rnmap.end()){
	  numfound++;
	  found=true;
	}
      }
    }

    if(found){
      int32 idoffound=rnI->second;
      // cout << "Found " << REP_thepool[idoffound].getName() << endl;
      // Read::setCoutType(Read::AS_TEXTCLIPS);
      // cout << REP_thepool[idoffound];

      int32 insertsize=-1;
      int32 insertstdev=-1;
      //int32 inssizemin=-1;
      //int32 inssizemax=-1;

      uint8 seqtype=ReadGroupLib::SEQTYPE_SANGER;

      E=T->elements.begin();
      ECD=T->elements_cdata.begin();
      for(;E!=T->elements.end(); E++, ECD++) {
	switch(*E) {
	case NCBIInfoXML::NCBIXML_TRACE_NAME : {
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_QUALITY_LEFT  : {
	  getRead(idoffound).setLQClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_QUALITY_RIGHT  : {
	  getRead(idoffound).setRQClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_VECTOR_LEFT  : {
	  getRead(idoffound).setLSClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_CLIP_VECTOR_RIGHT  : {
	  getRead(idoffound).setRSClipoff(atoi(ECD->c_str()));
	  break;
	}
	case NCBIInfoXML::NCBIXML_TEMPLATE_ID  : {
	  getRead(idoffound).setTemplate(ECD->c_str());
	  //cout<< *ECD << endl;
	  break;
	}
	case NCBIInfoXML::NCBIXML_TRACE_END  : {
	  if(strlen(ECD->c_str())>0){
	    switch(toupper(ECD->c_str()[0])){
	    case 'F': {
	      getRead(idoffound).setTemplateSegment(1);
	      break;
	    }
	    case 'R': {
	      getRead(idoffound).setTemplateSegment(255);
	      break;
	    }
	    case 'U' : // fall through
	    case 'N' : {
	      getRead(idoffound).setTemplateSegment(0);
	      break;
	    }
	    default : {
	      MIRANOTIFY(Notify::FATAL, "Illegal trace_end, it's not one of F/R/N(U) (or empty): " << ECD->c_str(););
	    }
	    }
	  }
	  break;
	}
	default : {
	  // Ooooops?
	}
	}
      }
      // cout << "After:\n";
      // cout << getRead(idoffound);

    }
  }

  makeTemplateIDs();

  cout << "Done merging XML data, matched " << numfound << " reads." << endl;

  rpDateStamp();

  FUNCEND();
  return;
}


/*************************************************************************
 *
 * ugly and slow, but works and is fast enough
 *
 *
 *************************************************************************/
//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void ReadPool::mergeSSAHA2SMALTVecScreenData(const string & ssahafile, bool issmalt, const string & logname, const string & logprefix)
{
  FUNCSTART("void ReadPool::mergeSSAHA2VecScreenData(const string & ssahafile, bool issmalt, const string & logname, const string & logprefix)");

  cout << "Merging vector screen data from ";
  if(issmalt){
    cout << "SMALT";
  }else{
    cout << "SSAHA2";
  }
  cout << " results file " << ssahafile << ":\n";

  CEBUG("Building hash table ... "); cout.flush();

  typedef boost::unordered_map<std::string, int32> strmap;
  strmap rnmap;
  strmap::iterator rnI;

  for(uint32 rpi=0; rpi<size();++rpi){
    if(!getRead(rpi).getName().empty()) {
      rnmap[getRead(rpi).getName()]=rpi;
    }
  }
  CEBUG("done." << endl);

  ofstream logfout;
  if(!logname.empty()){
    logfout.open(logname.c_str(), ios::out|ios::app);
    if(!logfout){
      MIRANOTIFY(Notify::FATAL, "Could not open log for appending: " << logname << "\nPossible causes: Disk full? Changed permissions? Directory deleted?");
    }
  }

  ifstream ssahafin;
  ssahafin.open(ssahafile.c_str(), ios::in|ios::ate);
  if(!ssahafin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << ssahafile);
  }
  streampos sfinsize=ssahafin.tellg();
  ssahafin.seekg(0, ios::beg);

  ProgressIndicator<streamsize> P (0, sfinsize,1000);

  uint32 sd_score;
  string sd_readname;
  string sd_vecname;
  uint32 sd_rfrom;
  uint32 sd_rto;
  uint32 sd_vfrom;
  uint32 sd_vto;
  string sd_dir;
  uint32 sd_totalmatchsize;
  float  sd_percentmatch;
  uint32 sd_rlen;

  string token;
  string alstring;
  if(issmalt){
    alstring="alignment:";
  }else{
    alstring="ALIGNMENT:";
  }

  bool haserrors=false;

  while(!ssahafin.eof()){
    ssahafin >> token;
    if(ssahafin.eof()) break;
    if(P.delaytrigger()) P.progress(ssahafin.tellg());
    if(token.compare(0,alstring.size(),alstring) != 0) {
      getline(ssahafin,token);
      continue;
    }
    ssahafin >> sd_score >> sd_readname;

    if(ssahafin.eof()) break;

    // *sigh* allow for empty names
    sd_vecname.clear();
    {
      bool loopit=true;
      char tmp;

      ssahafin.get(tmp);
      loopit=true;
      do{
	ssahafin.get(tmp);
	if(ssahafin.eof()) break;
	if(tmp==' ' || tmp=='\t'){
	  loopit=false;
	}else{
	  sd_vecname.push_back(tmp);
	}
      }while(loopit);
    }

    if(ssahafin.eof()) break;

    ssahafin >> sd_rfrom
	     >> sd_rto
	     >> sd_vfrom
	     >> sd_vto
	     >> sd_dir
	     >> sd_totalmatchsize
	     >> sd_percentmatch
	     >> sd_rlen;

    if(ssahafin.eof()) break;

    CEBUG(sd_readname << '\t' << sd_rfrom << '\t' << sd_rto << '\n');

    bool foundname=false;
    rnI=rnmap.find(sd_readname);
    if(rnI==rnmap.end()) {
      CEBUG("Not found: " << sd_readname << endl);
      continue;
    }
    uint32 foundreadid=rnI->second;
    if(!getRead(foundreadid).hasValidData()) continue;

    Read actread(getRead(foundreadid));
    assembly_parameters const & as_params= (*REP_miraparams)[actread.getSequencingType()].getAssemblyParams();

    if(actread.getLenSeq() != sd_rlen){
      if(actread.isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA) && actread.getLenSeq() != sd_rlen+1) {
	cout << "\nError! The length of read " << actread.getName()
	     << " (" << actread.getLenSeq()
	     << ") does not match the length given in the SSAHA2/SMALT file ("
	     << sd_rlen << ")\nSSAHA2 line:"
	     << ' ' << token
	     << ' ' << sd_score
	     << ' ' << sd_readname
	     << ' ' << sd_vecname
	     << ' ' << sd_rfrom
	     << ' ' << sd_rto
	     << ' ' << sd_vfrom
	     << ' ' << sd_vto
	     << ' ' << sd_dir
	     << ' ' << sd_totalmatchsize
	     << ' ' << sd_percentmatch
	     << ' ' << sd_rlen << endl;
	haserrors=true;
      }
    }

    CEBUG("SSAHA2/SMALT line:"
	  << ' ' << token
	  << ' ' << sd_score
	  << " r: " << sd_readname
	  << " v: " << sd_vecname
	  << " # " << sd_rfrom
	  << ' ' << sd_rto
	  << ' ' << sd_vfrom
	  << ' ' << sd_vto
	  << ' ' << sd_dir
	  << ' ' << sd_totalmatchsize
	  << ' ' << sd_percentmatch
	  << ' ' << sd_rlen << endl);

    //Read::setCoutType(Read::AS_FASTA);
    //CEBUG(actread);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    // in SSAHA2 output, from rfrom may be > rto for reverse matches
    // swap in these cases
    if(sd_rfrom > sd_rto) swap(sd_rfrom,sd_rto);

    for(uint32 i=sd_rfrom-1; i<sd_rto; i++){
      bool domask=false;
      if(as_params.as_clip_ssahamerge_strictfrontclip >0
	 || as_params.as_clip_ssahamerge_strictendclip >0){
	if(as_params.as_clip_ssahamerge_strictfrontclip >0
	   && static_cast<int32>(i)<as_params.as_clip_ssahamerge_strictfrontclip) domask=true;
	if(as_params.as_clip_ssahamerge_strictendclip>0
	   && i>=actread.getLenSeq()-as_params.as_clip_ssahamerge_strictendclip) domask=true;
      }else{
	domask=true;
      }
      if(domask) actread.changeBaseInSequence('X',0,i);
    }
    //Read::setCoutType(Read::AS_FASTA);
    //CEBUG(actread);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    actread.setClipoffsToMaskedChars(as_params.as_clip_ssahamerge_gapsize,
				     as_params.as_clip_ssahamerge_maxfrontgap,
				     as_params.as_clip_ssahamerge_maxendgap,
				     false);
    //Read::setCoutType(Read::AS_CLIPPEDFASTA);
    //CEBUG(actread);

    if(actread.getLMClipoff() > getRead(foundreadid).getLSClipoff()){
      getRead(foundreadid).setLSClipoff(actread.getLMClipoff());
      CEBUG("clippyl\n");
      if(!logname.empty()){
	logfout << logprefix << " SSAHA2/SMALT clip left "
		<< actread.getName()
		<< " to: "
		<< getRead(foundreadid).getLSClipoff() << '\n';
      }
    }else{
      if(!logname.empty()){
	logfout << logprefix << "unchanged SSAHA2/SMALT clip left "
		<< actread.getName()
		<< " stays: "
		<< getRead(foundreadid).getLSClipoff() << '\n';
      }
    }
    if(actread.getRMClipoff() < getRead(foundreadid).getRSClipoff()){
      getRead(foundreadid).setRSClipoff(actread.getRMClipoff());
      CEBUG("clippyr\n");
      if(!logname.empty()){
	logfout << logprefix << " SSAHA2/SMALT clip right "
		<< actread.getName()
		<< " to: "
		<< getRead(foundreadid).getRSClipoff() << '\n';
      }
    }else{
      if(!logname.empty()){
	logfout << logprefix << "unchanged SSAHA2/SMALT clip right "
		<< actread.getName()
		<< " stays: "
		<< getRead(foundreadid).getRSClipoff() << '\n';
      }
    }

    //Read::setCoutType(Read::AS_TEXTSHORT);
    //CEBUG(getRead(foundreadid));
  }
  P.finishAtOnce();

  ssahafin.close();

  if(!logname.empty()){
    logfout.close();
  }

  cout << "\nDone merging SSAHA2 vector screen data." << endl;

  if(haserrors){
    MIRANOTIFY(Notify::FATAL,"There were errors in the SSAHA2 data, most probably the sequences used to screen are different from\nthe ones loaded now (see log above). Sorry, MIRA has to abort, please check your data.");
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
Read & ReadPool::getRead(const string & name)
{
  return getRead(getReadIndex(name));
}





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/
void ReadPool::dumpAs(ostream & ostr, uint8 astype, bool alsoinvalids) const
{
  FUNCSTART("void ReadPool::dumpAs(ostream & ostr, uint8 astype, bool alsoinvalids) const)");

  if(astype==Read::AS_MAF){
    for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
      // use dumpReadGroupAsMAF() instead saveReadGroupAsMAF!
      ReadGroupLib::dumpReadGroupAsMAF(rgi,ostr);
    }
  }
  Read::setCoutType(astype);

  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    try{
      if(getRead(i).hasValidData() || alsoinvalids) ostr << getRead(i);
    }
    catch(Notify n){
      cout << "Ouch, ouch, ouch, ouch ... ouch. Error while output of read " << i << ", " << getRead(i).getName() << endl;
      n.handleError(THISFUNC);
    }
  }

  FUNCEND();
}


/*************************************************************************
 *
 * like above, but using save()
 * This is for miraconvert
 *
 *************************************************************************/
void ReadPool::saveAsMAF(ostream & ostr, bool alsoinvalids) const
{
  FUNCSTART("void ReadPool::dumpAs(ostream & ostr, uint8 astype, bool alsoinvalids) const)");

  for(uint32 rgi=1; rgi<ReadGroupLib::getNumReadGroups(); ++rgi){
    ReadGroupLib::saveReadGroupAsMAF(rgi,ostr);
  }

  Read::setCoutType(Read::AS_MAF);
  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    try{
      if(getRead(i).hasValidData() || alsoinvalids) ostr << getRead(i);
    }
    catch(Notify n){
      cout << "Ouch, ouch, ouch, ouch ... ouch. Error while output of read " << i << ", " << getRead(i).getName() << endl;
      n.handleError(THISFUNC);
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
void ReadPool::dumpPoolInfo(ostream & ostr) const
{
  FUNCSTART("void ReadPool::dumpPoolInfo(ostream & ostr)");

  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).hasValidData()) {
      ostr << i << '\t' << getRead(i).getName() << '\n';
    }else{
      ostr << i << "\tinvalid\n";
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
void ReadPool::dumpAsEXPs(string & dirname) const
{
  FUNCSTART("void ReadPool::dumpAsEXPs(string & dirname) const");

  if(ensureDirectory(dirname,true)){
    MIRANOTIFY(Notify::FATAL, "Could not make sure that directory '" << dirname << "' exists, aborting MIRA.");
  }

  ofstream fofnout((dirname+"/fofn").c_str(), ios::out | ios::trunc);

  string dummyAP="";
  for(uint32 i=0; i<REP_thepool3.size(); ++i){
    if(getRead(i).hasValidData()) {
      ofstream expout((dirname+"/"+getRead(i).getName()+".exp").c_str(), ios::out | ios::trunc);
      (const_cast<Read &>(getRead(i))).dumpAsGAP4DA(expout, dummyAP);

      expout.close();

      fofnout << getRead(i).getName() << ".exp" << endl;
    }
  }

  fofnout.close();

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void ReadPool::refreshNameIndex()
{
  FUNCSTART("void ReadPool::refreshNameIndex()");
  BUGIFTHROW(!REP_allownameindex,"Pool not allowed to make a name index?");

  REP_nameindex.clear();

  for(uint32 ri=0; ri<size(); ++ri){
    boost::unordered_map<string, uint32>::iterator mI=REP_nameindex.find(getRead(ri).getName());
    if(mI!=REP_nameindex.end()){
      MIRANOTIFY(Notify::FATAL,"Read '" << mI->first << "' present at least twice in read pool? Not good, you may not have more than one read with the same name.");
    }
    REP_nameindex[getRead(ri).getName()]=ri;
  }

  FUNCEND();
}

//  // returns (first) read with name i.p.
//  Read & getRead(const string & name);
//

/*************************************************************************
 *
 * returns id of (first) read with name or -1
 *
 *
 *************************************************************************/

int32 ReadPool::getReadIndex(const string & name)
{
  if(REP_nameindex.empty()) refreshNameIndex();

  boost::unordered_map<string, uint32>::iterator mI=REP_nameindex.find(name);
  if(mI!=REP_nameindex.end()){
    return mI->second;
  }
  return -1;
}



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void ReadPool::adjustIllegalQualities(base_quality_t bq)
{
  for(size_t i=0; i<REP_thepool3.size(); ++i){
    vector<base_quality_t> & bqv =const_cast<vector<base_quality_t>&>(getRead(i).getQualities());
    bool mustadjust=true;
    for(auto tbq: bqv){
      if(tbq<=100) {
	mustadjust=false;
	break;
      }
    }
    if(mustadjust){
      auto s=bqv.size();
      bqv.clear();
      bqv.resize(s,bq);
    }
  }

}
