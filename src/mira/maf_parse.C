/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2009 and later by Bastien Chevreux
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
 *
 */


#include "mira/maf_parse.H"

#include <boost/algorithm/string.hpp>


#include "errorhandling/errorhandling.H"
#include "util/progressindic.H"

#include "io/annotationmappings.H"


using namespace std;


#define CEBUG(bla)



MAFParse::MAFParse(ReadPool * rpool, list<Contig> * clist, vector<MIRAParameters> * mp)
{
  FUNCSTART("MAFParse::MAFParse(ReadPool * rpool, list<Contig> * clist, vector<MIRAParameters> * mp)");
  BUGIFTHROW(rpool==nullptr,"rpool==nullptr???");

  MAF_readpool   = rpool;
  MAF_contiglist = clist;
  MAF_miraparams = mp;
}


MAFParse::~MAFParse() {
  //discard();
}



//void MAFParse::discard() {
//}


size_t MAFParse::countReadsBeforeLoad(const string & fileName, size_t & maxlinelen)
{
  FUNCSTART("void MAFParse::countReadsBeforeLoad()");

  ifstream mafin;

  size_t numseqsloaded=0;
  maxlinelen=0;

  mafin.open(fileName.c_str(), ios::in|ios::ate);
  if(!mafin) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " not found for loading.");
  }
  if(!mafin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " is empty.");
  }

  ProgressIndicator<streamsize> P(0, mafin.tellg(),5000);

  mafin.seekg(0, ios::beg);

  string actline;
  string acttoken;

  actline.reserve(10000);

  while(!mafin.eof()){
    if(mafin.eof()) break;
    getline(mafin,actline);

    if(actline.size()>=2
       && actline[0]=='R'
       && actline[1]=='D') {
      numseqsloaded++;
    }
    if(P.delaytrigger()) P.progress(mafin.tellg());
  }
  P.finishAtOnce();

  mafin.close();

  maxlinelen=actline.capacity();

  FUNCEND();

  return numseqsloaded;
}



void MAFParse::cleanupHeaderData()
{
  MAF_vmajor=-1;
  MAF_vminor=-1;
}

void MAFParse::cleanupContigData()
{
  MAF_contig_assembledfrom.clear();
  MAF_contig_sequence.clear();
  MAF_contig_qualities.clear();
  MAF_contig_taglist.clear();

  MAF_contig_name.clear();

  MAF_contig_numreads=0;
  MAF_contig_len=0;
}

void MAFParse::cleanupReadData()
{
  MAF_read_sequence.clear();
  MAF_read_qualities.clear();
  MAF_read_align_origin.clear();
  MAF_read_taglist.clear();

  MAF_read_name.clear();
  MAF_read_scf_file.clear();
  MAF_read_template.clear();
  MAF_read_base_caller.clear();
  MAF_read_sequencing_vector.clear();
  MAF_read_strain.clear();
  MAF_read_machinetype.clear();

  MAF_read_len=-1;
  MAF_read_insert_size_min=-1;
  MAF_read_insert_size_max=-1;

  MAF_read_ql=-1;
  MAF_read_qr=-1;
  MAF_read_cl=-1;
  MAF_read_cr=-1;
  MAF_read_sl=-1;
  MAF_read_sr=-1;

  MAF_read_strand_given='N';
  MAF_read_tsegment_given=0;
  MAF_read_seqtype=ReadGroupLib::SEQTYPE_SANGER;

  MAF_read_isbackbone=false;
  MAF_read_israil=false;
  MAF_read_isCER=false;

  MAF_read_seenATline=false;

  MAF_readpoolid=-1;
}


/*
  seqtype = default seqtype of sequences if not encoded in the CAF
  loadaction:
    //  0 = count only
    //  1 = load
  lrperseqtype = longest read per seqtype

  returns:
    1) number of sequences loaded
    2) when loading: size of longest read per seqtype in lrperseqtype
 */

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
size_t MAFParse::load(const string & fileName, const uint8 seqtype, const uint8 loadaction, vector<uint32> & lrperseqtype, bool recalcconsensus, void (*ccallback)(list<Contig> &, ReadPool &), void (*rcallback)(ReadPool &), bool isVerbose)
{
  FUNCSTART("void MAFParse::load()");

  BUGIFTHROW(loadaction>1,"loadaction>1??");

  static const string cpsHVersion("@Version");
  static const string cpsHProgram("@Program");
  static const string cpsHReadGroup("@ReadGroup");
  static const string cpsHEndReadGroup("@EndReadGroup");
  static const string cpsHRG("@RG");

  static const string cpsRS("RS");
  static const string cpsRG("RG");
  static const string cpsRQ("RQ");
  static const string cpsRD("RD");
  static const string cpsLR("LR");
  static const string cpsSV("SV");
  static const string cpsTN("TN");
  static const string cpsDI("DI");
  static const string cpsTF("TF");
  static const string cpsTT("TT");
  static const string cpsTS("TS");
  static const string cpsSF("SF");
  static const string cpsBC("BC");
  static const string cpsSL("SL");
  static const string cpsSR("SR");
  static const string cpsQL("QL");
  static const string cpsQR("QR");
  static const string cpsCL("CL");
  static const string cpsCR("CR");
  static const string cpsAO("AO");
  static const string cpsRT("RT");
  static const string cpsST("ST");
  static const string cpsSN("SN");
  static const string cpsMT("MT");
  static const string cpsIB("IB");
  static const string cpsIC("IC");
  static const string cpsIR("IR");
  static const string cpsER("ER");
  static const string cpsCS("CS");
  static const string cpsCQ("CQ");
  static const string cpsCO("CO");
  static const string cpsNR("NR");
  static const string cpsLC("LC");
  static const string cpsCT("CT");
  static const string cpsSLSL("//");
  static const string cpsAT("AT");
  static const string cpsEC("EC");
  static const string cpsBSBS("\\\\");

  MAF_ccallbackfunc=ccallback;
  MAF_rcallbackfunc=rcallback;
  MAF_recalcconsensus=recalcconsensus;

  size_t numseqsloaded=0;
  size_t maxlinelen=0;
  if(loadaction==0){
    cout << "Counting reads:\n";
    numseqsloaded=countReadsBeforeLoad(fileName,maxlinelen);
    return numseqsloaded;
  }

  cout << "Loading MAF " << fileName << " :\n";

  MAF_lrperseqtype.clear();
  MAF_lrperseqtype.resize(ReadGroupLib::SEQTYPE_END,0);

  ifstream mafin;

  mafin.open(fileName.c_str(), ios::in|ios::ate);
  if(!mafin) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " not found for loading.");
  }
  if(!mafin.tellg() ) {
    MIRANOTIFY(Notify::FATAL, "MAF file " << fileName << " is empty.");
  }

  ProgressIndicator<streamsize> P(0, mafin.tellg(),5000);

  mafin.seekg(0, ios::beg);

  string acttoken;
  string actline;

  size_t linenumber=0;

  if(maxlinelen>10000) {
    actline.reserve(maxlinelen);
  }else{
    actline.reserve(10000);
  }

  MAF_isinread=false;
  MAF_isincontig=false;
  MAF_isinreadgroup=false;

  cleanupHeaderData();
  cleanupReadData();
  cleanupContigData();

  try {
    while(true){
      ++linenumber;
      mafin >> acttoken;
      if(mafin.eof()) break;

      CEBUG("l: " << linenumber << "\tt: ###" << acttoken << "###" << endl);

      if(acttoken.empty()) continue;


/* here for read*/

      if(acttoken==cpsRD){
	// read name
	parseLineRD(mafin,acttoken,actline);
      }else if(acttoken==cpsRG){
	// Read Group
	parseLineRG(mafin,acttoken,actline);
      }else if(acttoken==cpsRS){
	// Read Sequence
	parseLineRS(mafin,acttoken,actline);
      }else if(acttoken==cpsRQ){
	// Read Qualities
	parseLineRQ(mafin,acttoken,actline);
      }else if(acttoken==cpsLR){
	// length read
	parseLineLR(mafin,acttoken,actline);
      }else if(acttoken==cpsSV){
	// sequencing vector
	parseLineSV(mafin,acttoken,actline);
      }else if(acttoken==cpsTN){
	// template name
	parseLineTN(mafin,acttoken,actline);
      }else if(acttoken==cpsDI){
	// Direction (strand)
	parseLineDI(mafin,acttoken,actline);
      }else if(acttoken==cpsTF){
	// template insize from
	parseLineTF(mafin,acttoken,actline);
      }else if(acttoken==cpsTT){
	// template insize to
	parseLineTT(mafin,acttoken,actline);
      }else if(acttoken==cpsTS){
	// template segment
	parseLineTS(mafin,acttoken,actline);
      }else if(acttoken==cpsSF){
	// Sequencing File
	parseLineSF(mafin,acttoken,actline);
      }else if(acttoken==cpsBC){
	// base caller
	parseLineBC(mafin,acttoken,actline);
      }else if(acttoken==cpsSL){
	//
	parseLineSL(mafin,acttoken,actline);
      }else if(acttoken==cpsSR){
	//
	parseLineSR(mafin,acttoken,actline);
      }else if(acttoken==cpsQL){
	//
	parseLineQL(mafin,acttoken,actline);
      }else if(acttoken==cpsQR){
	//
	parseLineQR(mafin,acttoken,actline);
      }else if(acttoken==cpsCL){
	//
	parseLineCL(mafin,acttoken,actline);
      }else if(acttoken==cpsCR){
	//
	parseLineCR(mafin,acttoken,actline);
      }else if(acttoken==cpsAO){
	//
	parseLineAO(mafin,acttoken,actline);
      }else if(acttoken==cpsRT){
	//
	parseLineRT(mafin,acttoken,actline);
      }else if(acttoken==cpsST){
	//
	parseLineST(mafin,acttoken,actline);
      }else if(acttoken==cpsSN){
	//
	parseLineSN(mafin,acttoken,actline);
      }else if(acttoken==cpsMT){
	//
	parseLineMT(mafin,acttoken,actline);
      }else if(acttoken==cpsIB){
	//
	parseLineIB(mafin,acttoken,actline);
      }else if(acttoken==cpsIC){
	//
	parseLineIC(mafin,acttoken,actline);
      }else if(acttoken==cpsIR){
	//
	parseLineIR(mafin,acttoken,actline);
      }else if(acttoken==cpsER){
	//
	parseLineER(mafin,acttoken,actline);

/* here for contig*/
      }else if(acttoken==cpsCS){
	// Consensus Sequence
	parseLineCS(mafin,acttoken,actline);
      }else if(acttoken==cpsCQ){
	// Consensus Qualities
	parseLineCQ(mafin,acttoken,actline);
      }else if(acttoken==cpsCO){
	// COntig name
	parseLineCO(mafin,acttoken,actline);
      }else if(acttoken==cpsNR){
	// Num Reads
	parseLineNR(mafin,acttoken,actline);
      }else if(acttoken==cpsLC){
	// Length Contig
	parseLineLC(mafin,acttoken,actline);
      }else if(acttoken==cpsCT){
	// Contig Tag
	parseLineCT(mafin,acttoken,actline);
      }else if(acttoken==cpsSLSL){
	// start of contig reads

      }else if(acttoken==cpsBSBS){
	// end of contig reads

      }else if(acttoken==cpsAT){
	// Assembled From
	parseLineAT(mafin,acttoken,actline);
      }else if(acttoken==cpsEC){
	// end contig
	parseLineEC(mafin,acttoken,actline);


/* here for header*/

      }else if(acttoken==cpsHVersion){
	// file version
	parseLineHeaderVersion(mafin,acttoken,actline);
      }else if(acttoken==cpsHProgram){
	// simply read the rest of the line and do nothing
	getline(mafin,acttoken);
      }else if(acttoken==cpsHReadGroup){
	// file version
	parseLineHeaderReadGroup(mafin,acttoken,linenumber);
      }else if(acttoken==cpsHEndReadGroup
	       || acttoken==cpsHRG){
	cout << "File " << fileName << ": around line " << linenumber
	     << "\n" << acttoken << " occurred without being in an open read group?\n";
	MIRANOTIFY(Notify::FATAL, "Error while reading MAF file.");
      }else{
	cout << "File " << fileName << ": around line " << linenumber
	     << "\ndid not recognize token " << acttoken << '\n';
	MIRANOTIFY(Notify::FATAL, "Error while reading MAF file.");
      }

      if(P.delaytrigger()) P.progress(mafin.tellg());
    }

    P.finishAtOnce();
    cout << endl;

    if(MAF_isinread){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open read: " << MAF_read_name << " .... file truncated?");
    }
    if(MAF_isincontig){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open contig: " << MAF_contig_name << " ... file truncated?");
    }
    if(MAF_isinreadgroup){
      MIRANOTIFY(Notify::FATAL, "MAF file ends without closing an open readgroup. File truncated?");
    }
  }
  catch(Notify n){
    cout << "\nError around line " << linenumber << " of file " << fileName << endl;
    cout << "Last contig name read: " << MAF_contig_name << endl;
    cout << "Last read name read: " << MAF_read_name << endl;
    n.handleError(THISFUNC);
  }

  lrperseqtype=MAF_lrperseqtype;

  FUNCEND();

  return numseqsloaded;
}
//#define CEBUG(bla)


void MAFParse::checkParseIsInReadGroup(string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInReadGroup(string & acttoken)");
  if(!MAF_isinreadgroup) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in readgroup (@ReadGroup line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsNotInReadGroup(string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsNotInReadGroup(string & acttoken)");
  if(MAF_isinreadgroup) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while in readgroup (@EndReadGroup line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsInRead(string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInRead(string & acttoken)");
  if(!MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in read (RD line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsNotInRead(string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsNotInRead(string & acttoken)");
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while bein in read (ER line missing?)");
  }
  FUNCEND();
}

void MAFParse::checkParseIsInContig(string & acttoken)
{
  FUNCSTART("void MAFParse::checkParseIsInContig(string & acttoken)");
  if(!MAF_isincontig) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while not in contig (CO line missing?)");
  }
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered " << acttoken << " line while being in read (RD line not closed by ER?)");
  }
  FUNCEND();
}




void MAFParse::parseLineRD(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineRD(ifstream & mafin, string & acttoken, string & actline)");
  if(MAF_isinread) {
    MIRANOTIFY(Notify::FATAL,"Encountered new " << acttoken << " line when the previous read " << MAF_read_name << " was not closed with 'ER'");
  }
  checkParseIsNotInReadGroup(acttoken);

  cleanupReadData();
  mafin >> MAF_read_name;
  MAF_isinread=true;

  FUNCEND();
}

void MAFParse::parseLineRG(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineRG(ifstream & mafin, string & acttoken, string & actline)");
  checkParseIsInRead(acttoken);

  mafin >> MAF_tmp_str;
  int32 dummy=atoi(MAF_tmp_str.c_str());
  BUGIFTHROW(dummy<0 || dummy >65535,"Line RG: id must be >=0 and <= 65535, but " << dummy << " was given.");
  BUGIFTHROW(dummy>=MAF_readgroup_externalidmapper.size()+1,"Line RG: id of " << dummy << " was given, but not readgroup with this id was defined (@RG ID)");
  BUGIFTHROW(MAF_readgroup_externalidmapper[dummy].isDefaultNonValidReadGroupID(),"Line RG: id of " << dummy << " was given, but not readgroup with this id was defined (@RG ID)");
  MAF_read_rgid=MAF_readgroup_externalidmapper[dummy];

  FUNCEND();
}

void MAFParse::parseLineRS(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineRS(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInRead(acttoken);

  if(!MAF_read_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered RS line when there already was one for read " << MAF_read_name);
  }

  mafin >> actline;

  MAF_read_sequence.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_read_sequence.push_back(*seq);
  }

  FUNCEND();
}

void MAFParse::parseLineRQ(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineRQ(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInRead(acttoken);

  if(!MAF_read_qualities.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered RQ line when there already was one for read " << MAF_read_name);
  }

  mafin >> actline;

  MAF_read_qualities.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_read_qualities.push_back(*seq-33);
  }

  FUNCEND();
}

void MAFParse::parseLineLR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_len;
  if(MAF_read_len>actline.capacity()) actline.reserve(MAF_read_len+10);
}

void MAFParse::parseLineSV(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_sequencing_vector;
}

void MAFParse::parseLineTN(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_template;
}

void MAFParse::parseLineDI(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_strand_given;
}

void MAFParse::parseLineTF(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_insert_size_min;
}

void MAFParse::parseLineTT(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_insert_size_max;
}

void MAFParse::parseLineTS(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  uint32 dummy;     // sigh ... detour to read a number and NOT a 'char';
  mafin >> dummy;
  MAF_read_tsegment_given=static_cast<uint8>(dummy);
}

void MAFParse::parseLineSF(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_scf_file;
}

void MAFParse::parseLineBC(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_base_caller;
}

void MAFParse::parseLineSL(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_sl;
  MAF_read_sl--;
}

void MAFParse::parseLineSR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_sr;
}

void MAFParse::parseLineQL(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_ql;
  MAF_read_ql--;
}

void MAFParse::parseLineQR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_qr;
}

void MAFParse::parseLineCL(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_cl;
  MAF_read_cl--;
}

void MAFParse::parseLineCR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_cr;
}

void MAFParse::parseLineAO(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineAO(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInRead(acttoken);

  if(MAF_read_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ": sequence (SQ line) must be defined before AO line");
  }
  if(MAF_read_align_origin.empty()){
    MAF_read_align_origin.resize(MAF_read_sequence.size(),-1);
  }

  int32 seqfrom, seqto, origfrom, origto;
  mafin >> seqfrom;
  mafin >> seqto;
  mafin >> origfrom;
  mafin >> origto;

  //cout << "xxx " << seqfrom << " " << seqto << " " << origfrom << " " << origto << "\n";

  if(seqfrom<1
     || seqto<1
     || origfrom<1
     || origto<1){
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  values may not be <1");
  }

  int32 seqinc=1;
  int32 originc=1;
  if(seqto<seqfrom) seqinc=-1;
  if(origto<origfrom) originc=-1;


  if (abs(seqto - seqfrom) != abs(origto - origfrom)) {
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  distance between seqfrom/to (" << seqfrom << " " << seqto << ") is unequal to originalfrom/to (" << origfrom << " " << origto << ")");
  }
  if (max(seqfrom, seqto) > static_cast<int32>(MAF_read_align_origin.size())) {
    MIRANOTIFY(Notify::FATAL,"While reading AO line for read " << MAF_read_name << ":  seqfrom/to (" << seqfrom << " " << seqto << ") is larger than size of read (" << MAF_read_align_origin.size() << ")");
  }

  int32 seqi=seqfrom;
  int32 origi=origfrom;
  for(int32 loopi=0; loopi < abs(seqto-seqfrom)+1; loopi++){
    if(seqi-1 <0 || seqi-1 >= MAF_read_align_origin.size()) {
      MIRANOTIFY(Notify::FATAL,"While reading AO line for read: " << MAF_read_name << " with AO values " << seqfrom << " " << seqto << " " <<origfrom << " " << origto << "\nThis makes for an illegal alignment of the sequence to the original, wrong values in this line?\n");
    }

    MAF_read_align_origin[seqi-1] = origi - 1;
    seqi += seqinc;
    origi += originc;
  }

  FUNCEND();
}

void MAFParse::parseLineRT(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);

  MAF_read_taglist.resize(MAF_read_taglist.size()+1);
  if(MAF_vmajor<2){
    parseTagData(mafin,acttoken,MAF_read_taglist.back());
  }else{
    parseTagDataV2(mafin,acttoken,MAF_read_taglist.back());
  }

  //cout << "Stored Rtag: " << MAF_read_taglist.back();
  //MAF_read_taglist.back().dumpDebug(cout);
}

void MAFParse::parseTagData(ifstream & mafin, string & acttoken, multitag_t & targettag)
{
  FUNCSTART("void MAFParse::parseTagData(ifstream & mafin, string & acttoken, multitag_t & tag)");

  multitag_t tmptag;

  mafin >> MAF_tmp_str;

  if(!AnnotationMappings::isValidGFF3SOEntry(MAF_tmp_str)){
    string soident(AnnotationMappings::translateGAP4feat2SOfeat(MAF_tmp_str));
    if(soident.empty()){
      soident=AnnotationMappings::translateXGAP4feat2SOfeat(MAF_tmp_str);
      if(soident.empty()){
	tmptag.setIdentifierStr(MAF_tmp_str);
      }
    }
    if(!soident.empty()){
      tmptag.setIdentifierStr(soident);
    }
  }else{
    tmptag.setIdentifierStr(MAF_tmp_str);
  }

  mafin >> tmptag.from;
  mafin >> tmptag.to;

  if(tmptag.from<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tmptag line " << acttoken << ": (" << tmptag.from << " " << tmptag.to << ") -> " << tmptag.from << " is <1, not allowed.");
  }
  if(tmptag.to<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tmptag line " << acttoken << ": (" << tmptag.from << " " << tmptag.to << ") -> " << tmptag.to << " is <1, not allowed.");
  }

  tmptag.from-=1;
  tmptag.to-=1;


  if(tmptag.from<=tmptag.to){
    tmptag.setStrand('+');
  }else{
    tmptag.setStrand('-');
    swap(tmptag.from, tmptag.to);
  }

  // comment may be present or not
  char nextchar;
  mafin.get(nextchar);
  if(nextchar=='\n') {
    targettag=tmptag;
    return;
  }
  if(nextchar=='\r') {
    // also eat \n
    mafin.get(nextchar);
    targettag=tmptag;
    return;
  }
  getline(mafin,MAF_tmp_str);
  tmptag.setCommentStr(MAF_tmp_str);

//  cout << "BEFORE: ";
//  tmptag.dumpDebug(cout);
  Read::upgradeOldTagToMultitagWithGFF3(tmptag,targettag);
//  cout << "AFTER: ";
//  targettag.dumpDebug(cout);

  if(targettag.getSourceStr().empty()){
    if(AnnotationMappings::translateSOfeat2GAP4feat(targettag.getIdentifierStr()).empty()
       && AnnotationMappings::translateSOfeat2XGAP4feat(targettag.getIdentifierStr()).empty()){
      targettag.source = multitag_t::MT_tagsrcentry_idMIRA;
    }else{
      targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
    }
    targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
  }else if(targettag.source == multitag_t::MT_tagsrcentry_idMIRA){
    if(!AnnotationMappings::translateSOfeat2SOID(targettag.getIdentifierStr()).empty()){
      targettag.source = multitag_t::MT_tagsrcentry_idGFF3;
    }
  }

  FUNCEND();
}



void MAFParse::parseTagDataV2(ifstream & mafin, string & acttoken, multitag_t & targettag)
{
  // Raw speed to parse tab delimited line
  // both boost::split and boost::tokenizer are 40% to 50% slower.
  //
  // Cheating big time: the line string gets rewritten, replacing tabs by 0 and
  //  writing substring pointer (char *) to a fixed array
  //
  // When throwing with MAF_tmp_str in the message, we need to re-tabify beforehand!

  FUNCSTART("void MAFParse::parseTagData(ifstream & mafin, string & acttoken, multitag_t & tag)");

  targettag.commentisgff3=true;

  {
    char nextchar;
    mafin.get(nextchar); // eat away the following \t
  }
  getline(mafin,MAF_tmp_str);

  static char * sarr[8];

  //////////// Overwrite \t with 0, populate the array pf char * to subparts of line
  char * sptr=const_cast<char *>(MAF_tmp_str.c_str());
  char * rptr=sptr;

  uint32 numtabs=0;
  while(*rptr){
    if(*rptr=='\t'){
      *rptr=0;
      if(numtabs>8){
	for(auto & x : MAF_tmp_str) if(x==0) x='\t';
	MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found more");
      }
      sarr[numtabs]=sptr;
      ++numtabs;
      sptr=rptr+1;
    }
    ++rptr;
  }
  if(rptr!=MAF_tmp_str.c_str()){
    sarr[numtabs]=sptr;
    ++numtabs;
  }

  if(numtabs<4 || numtabs>7){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found " << numtabs);
  }

  for(;numtabs<8;++numtabs) sarr[numtabs]=nullptr;
  ///// done

  static string acts;
  acts=sarr[0];

  if(AnnotationMappings::isMIRAEntry(acts)){
    targettag.setIdentifierStr(acts);
  }else{
    const string * soident=&AnnotationMappings::translateGAP4feat2SOfeat(acts);
    if(soident->empty()){
      soident=&AnnotationMappings::translateXGAP4feat2SOfeat(acts);
    }
    if(soident->empty()){
      targettag.setIdentifierStr(acts);
    }else{
      targettag.setIdentifierStr(*soident);
    }
  }

  targettag.from=atoi(sarr[1]);
  targettag.to=atoi(sarr[2]);

  if(targettag.from<1){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.from << " is <1, not allowed.");
  }
  if(targettag.to<1){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.to << " is <1, not allowed.");
  }

  targettag.from-=1;
  targettag.to-=1;

  if(*sarr[3]==0){
    for(auto & x : MAF_tmp_str) if(x==0) x='\t';
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": the entry for the strand is empty, not allowed.");
  }

  targettag.setStrand(*sarr[3]);

  if(sarr[4]!=nullptr){
    targettag.setSourceStr(sarr[4]);

    if(sarr[5]!=nullptr){
      // hmmmmmm .... something was very wrong with the previous implementation
      switch(*sarr[5]){
      case '.':
      case '3':{
	targettag.phase=3;
	break;
      }
      case '0':{
	targettag.phase=0;
	break;
      }
      case '1':{
	targettag.phase=1;
	break;
      }
      case '2':{
	targettag.phase=2;
	break;
      }
      default : {
	for(auto & x : MAF_tmp_str) if(x==0) x='\t';
	MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << "\t" << MAF_tmp_str << "\nthe entry for the strand (" << *sarr[5] << ") is not 0, 1, 2, 3 or .");
      }
      }

      if(sarr[6]!=nullptr){
	targettag.setCommentStr(sarr[6]);
      }
    }
  }

  FUNCEND();
}


void MAFParse::parseLineST(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineST(ifstream & mafin, string & acttoken, string & actline)");
  checkParseIsInRead(acttoken);
  mafin >> actline;

  MAF_read_seqtype=ReadGroupLib::stringToSeqType(actline);

  if(MAF_read_seqtype==ReadGroupLib::SEQTYPE_END){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in tag " << acttoken << ": unknown sequencing type '" << actline << "'?");
  }
  FUNCEND();
}

void MAFParse::parseLineSN(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_strain;
}

void MAFParse::parseLineMT(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_machinetype;
}

void MAFParse::parseLineIB(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_isbackbone;
}

void MAFParse::parseLineIC(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_isCER;
}

void MAFParse::parseLineIR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  mafin >> MAF_read_israil;
}

void MAFParse::parseLineER(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInRead(acttoken);
  addReadToReadPool();
  MAF_isinread=false;
  if(!MAF_isincontig && MAF_rcallbackfunc!=nullptr) {
    (*MAF_rcallbackfunc)(*MAF_readpool);
  }
  MAF_read_rgid.resetLibId();
}


void MAFParse::parseLineCO(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineCO(ifstream & mafin, string & acttoken, string & actline)");
  if(MAF_isincontig){
    MIRANOTIFY(Notify::FATAL, "Seen new CO line while previous CO was not closed by EC");
  }
  checkParseIsNotInReadGroup(acttoken);

  cleanupContigData();
  mafin >> MAF_contig_name;
  MAF_isincontig=true;
  FUNCEND();
}

void MAFParse::parseLineCS(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineSQ(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInContig(acttoken);

  if(!MAF_contig_sequence.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered CS line when there already was one for contig " << MAF_contig_name);
  }

  mafin >> actline;

  MAF_contig_sequence.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_contig_sequence.push_back(*seq);
  }

  FUNCEND();
}

void MAFParse::parseLineCQ(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineCQ(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInContig(acttoken);

  if(!MAF_contig_qualities.empty()){
    MIRANOTIFY(Notify::FATAL,"Encountered CQ line when there already was one for contig " << MAF_contig_name);
  }

  mafin >> actline;

  MAF_contig_qualities.reserve(actline.size());
  const char * seq=actline.c_str();
  for(; *seq; ++seq){
    MAF_contig_qualities.push_back(*seq-33);
  }

  FUNCEND();
}

void MAFParse::parseLineNR(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInContig(acttoken);
  mafin >> MAF_contig_numreads;
}

void MAFParse::parseLineLC(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInContig(acttoken);
  mafin >> MAF_contig_len;

  if(MAF_contig_len>actline.capacity()) actline.reserve(MAF_contig_len+10);
}

void MAFParse::parseLineCT(ifstream & mafin, string & acttoken, string & actline)
{
  checkParseIsInContig(acttoken);

  MAF_contig_taglist.resize(MAF_contig_taglist.size()+1);
  if(MAF_vmajor<2){
    parseTagData(mafin,acttoken,MAF_contig_taglist.back());
    // for consensus tags, change strand to '='
    // BaCh 21.01.2012: Why?
    //MAF_contig_taglist.back().setStrand('=');
  }else{
    parseTagDataV2(mafin,acttoken,MAF_contig_taglist.back());
  }

  //cout << "Stored Ctag: ";
  //MAF_contig_taglist.back().dumpDebug(cout);
}

void MAFParse::parseLineAT(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineAT(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInContig(acttoken);
  checkParseIsNotInRead(acttoken);

  if(MAF_readpoolid<0){
    MIRANOTIFY(Notify::FATAL, "Seen AT line but no read in contig defined before? (RD/ER block in a CO block)");
  }
  if(MAF_read_seenATline){
    MIRANOTIFY(Notify::FATAL, "Seen AT line, but either no read before or multiple AT lines.");
  }
  MAF_read_seenATline=true;

  int32 cfrom,cto,rfrom,rto;
  Contig::contig_init_read_t tmpcr;
  int8 direction=1;

  mafin >> cfrom;
  mafin >> cto;
  mafin >> rfrom;
  mafin >> rto;

  if(cfrom > cto){
    direction=-1;
    tmpcr.offset_in_contig=cto-1;
  }else{
    tmpcr.offset_in_contig=cfrom-1;
  }

  if(rfrom>rto){
    tmpcr.read_rclip=rfrom;
    tmpcr.read_lclip=rto-1;
    tmpcr.direction= -direction;
  }else{
    tmpcr.read_rclip=rto;
    tmpcr.read_lclip=rfrom-1;
    tmpcr.direction= direction;
  }

  tmpcr.id=MAF_readpoolid;
  MAF_contig_assembledfrom.push_back(tmpcr);

  FUNCEND();
}

void MAFParse::parseLineEC(ifstream & mafin, string & acttoken, string & actline)
{
  FUNCSTART("void MAFParse::parseLineEC(ifstream & mafin, string & acttoken, string & actline)");

  checkParseIsInContig(acttoken);
  MAF_isincontig=false;

  if(MAF_contiglist!=nullptr){
    {
      Contig dummy(MAF_miraparams, *MAF_readpool);
      MAF_contiglist->push_back(dummy);
    }

    try{
      if(MAF_recalcconsensus){
	string dummy1;
	vector<base_quality_t> dummy2;
	MAF_contiglist->back().initialiseContig(MAF_contig_assembledfrom,
						MAF_contig_taglist,
						MAF_contig_name,
						dummy1,dummy2);
      }else{
	string dummy1;
	dummy1.reserve(MAF_contig_sequence.size()+2);
	for(auto dnaI=MAF_contig_sequence.cbegin(); dnaI != MAF_contig_sequence.end(); dnaI++) dummy1+=*dnaI;
	MAF_contiglist->back().initialiseContig(MAF_contig_assembledfrom,
						MAF_contig_taglist,
						MAF_contig_name,
						dummy1,
						MAF_contig_qualities);
      }
    }
    catch(Notify n){
      cout << "Error for contig " << MAF_contig_name << endl;
      n.handleError(THISFUNC);
    }

    try{
      if(MAF_ccallbackfunc!=nullptr) {
	(*MAF_ccallbackfunc)(*MAF_contiglist, *MAF_readpool);
      }
    }
    catch(Notify n){
      cout << "Error while calling callback!\n";
      cout << "Error for contig " << MAF_contig_name << endl;
      n.handleError(THISFUNC);
    }
  }

  FUNCEND();
}



void MAFParse::checkReadData()
{
  FUNCSTART("void MAFParse::checkReadData()");

  BUGIFTHROW(MAF_read_seqtype>=ReadGroupLib::SEQTYPE_END, "Illegal seqtype in checkReadData()???");

  if(MAF_read_len>=0
     && MAF_read_sequence.size() != MAF_read_len){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": size of sequence (" << MAF_read_sequence.size() << ") is not equal to size given in LR line (" << MAF_read_len << ")");
  }
  if(!MAF_read_qualities.empty()){
    if(MAF_read_sequence.size() != MAF_read_qualities.size()){
      MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": size of sequence (" << MAF_read_sequence.size() << ") is not equal to size of qualities (" << MAF_read_qualities.size() << ")");
    }
  }else{
    // sequence but no qualities ... then fake some according to the sequencing type
    // give them the standard qual for this sequencing type
    MAF_read_qualities.resize(MAF_read_sequence.size(),MAF_read_rgid.getDefaultQual());
  }

  if(!MAF_read_align_origin.empty() && MAF_read_align_origin.size() != MAF_read_sequence.size()){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << ": the align to origin (AO) data led to a larger or smaller array that the length of the sequence?");
  }

  if(MAF_read_ql < 0) MAF_read_ql = 0;
  if(MAF_read_sl < 0) MAF_read_sl = 0;
  if(MAF_read_cl < 0) MAF_read_cl = 0;

  if(MAF_read_qr < 0) MAF_read_qr = MAF_read_sequence.size();
  if(MAF_read_sr < 0) MAF_read_sr = MAF_read_sequence.size();
  if(MAF_read_cr < 0) MAF_read_cr = MAF_read_sequence.size();

  FUNCEND();
}


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void MAFParse::addReadToReadPool()
{
  FUNCSTART("void MAFParse::addReadToReadPool()");

  checkReadData();

  MAF_readpoolid=MAF_readpool->size();
  Read & newread = MAF_readpool->getRead(MAF_readpool->provideEmptyRead());

  ReadGroupLib::ReadGroupID rgid;

  if(MAF_vmajor>=2 && MAF_read_rgid.isDefaultNonValidReadGroupID()){
    MIRANOTIFY(Notify::FATAL,"Read " << MAF_read_name << " has no RG line to define the read group. This is needed for MAF version 2");
  }
  if(MAF_read_rgid.isDefaultNonValidReadGroupID()){
    string rgname;
    rgid=ReadGroupLib::searchExactRGMatch(
      rgname,
      MAF_read_seqtype,
      MAF_read_insert_size_min,
      MAF_read_insert_size_max,
      ReadGroupLib::SPLACE_UNKNOWN,
      MAF_read_strain,
      MAF_read_isbackbone,
      MAF_read_israil,
      MAF_read_isCER,
      MAF_read_sequencing_vector,
      MAF_read_machinetype,
      MAF_read_base_caller);

    if(rgid.isDefaultNonValidReadGroupID()){
      rgid=ReadGroupLib::newReadGroup();
      rgid.setGroupName(rgname);
      rgid.setSequencingType(MAF_read_seqtype);
      rgid.setInsizeFrom(MAF_read_insert_size_min);
      rgid.setInsizeTo(MAF_read_insert_size_max);
      rgid.setSegmentPlacement("?");
      rgid.setStrainName(MAF_read_strain);
      rgid.setBackbone(MAF_read_isbackbone);
      rgid.setRail(MAF_read_israil);
      rgid.setCoverageEquivalentRead(MAF_read_isCER);
      rgid.setSeqVecName(MAF_read_sequencing_vector);
      rgid.setMachineType(MAF_read_machinetype);
      rgid.setBaseCaller(MAF_read_base_caller);
    }
  }else{
    rgid=MAF_read_rgid;
  }

  //cout << "Setting rgid to:\n" << rgid << endl;

  if(rgid.getSequencingType()==ReadGroupLib::SEQTYPE_SOLEXA
     || rgid.getSequencingType()==ReadGroupLib::SEQTYPE_IONTORRENT){
    newread.disallowAdjustments();
  }

  if(newread.usesAdjustments() && MAF_read_align_origin.empty() && !MAF_read_sequence.empty()){
    // no AO line given, i.e., no insertions/deletions
    // create vector which represents that
    MAF_read_align_origin.resize(MAF_read_sequence.size());
    uint32 num=0;
    for(auto & x : MAF_read_align_origin){
      x=num++;
    }
  }

  newread.initialiseRead(false,
			 false,
			 true,     // always padded
			 rgid,
			 MAF_read_sequence,
			 MAF_read_qualities,
			 MAF_read_align_origin,
			 MAF_read_taglist,
			 MAF_read_name,
			 MAF_read_scf_file,
			 MAF_read_ql,
			 MAF_read_qr,
			 MAF_read_sl,
			 MAF_read_sr,
			 MAF_read_cl,
			 MAF_read_cr);


  if(MAF_read_tsegment_given!=0){
    newread.setTemplateSegment(MAF_read_tsegment_given);
  }else if(MAF_read_strand_given!='N'){
    if(MAF_read_strand_given=='F'){
      newread.setTemplateSegment(1);
    }else{
      newread.setTemplateSegment(255);
    }
  }
  if (!MAF_read_template.empty()) {
    newread.setTemplate(MAF_read_template);
  }

  CEBUG("ER read is:\n" << newread << endl);
}
#define CEBUG(bla)



void MAFParse::parseLineHeaderVersion(ifstream & mafin, string & acttoken, string & actline)
{
  mafin >> MAF_tmp_str;
  MAF_vmajor=atoi(MAF_tmp_str.c_str());
  mafin >> MAF_tmp_str;
  MAF_vminor=atoi(MAF_tmp_str.c_str());
}

void MAFParse::parseLineHeaderReadGroup(ifstream & mafin, string & acttoken, size_t & linenumber)
{
  FUNCSTART("void MAFParse::parseLineHeaderReadGroup(ifstream & mafin, string & acttoken, size_t & linenumber)");

  checkParseIsNotInReadGroup(acttoken);

  MAF_readgroup_rgid=ReadGroupLib::newReadGroup();
  parseReadGroup(mafin,MAF_readgroup_rgid,MAF_readgroup_externalidmapper,linenumber);
  MAF_readgroup_rgid.fillInSensibleDefaults();
  MAF_readgroup_rgid.resetLibId();

  FUNCEND();
}


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void MAFParse::parseReadGroup(ifstream & mafin, ReadGroupLib::ReadGroupID & rgid, vector<ReadGroupLib::ReadGroupID> & externalidmapper, size_t & linenumber)
{
  FUNCSTART("void MAFParse::parseReadGroup(ifstream & mafin, ReadGroupLib::ReadGroupID & rgid, size_t & linenumber)");

  string mafline;
  vector<string> mafsplit;

  while(true){
    ++linenumber;
    getline(mafin,mafline);
    CEBUG("MLRG: " << mafline << endl);
    if(mafin.eof()) break;
    if(mafline.empty()) continue;
    boost::split(mafsplit, mafline, boost::is_any_of("\t"));
    if(mafsplit.empty()) continue;
    if(mafsplit.size()==1){
      if(mafsplit[0]=="@EndReadGroup") break;
      cout << "\nOuch, erroneous line: " << mafline << endl;
      MIRANOTIFY(Notify::FATAL,"Did not find a tab character in line " << linenumber << " and keyword is not @EndReadGroup? Something is wrong.");
    }

    string & rgtoken=mafsplit[1];
    CEBUG("read rgtoken #" << rgtoken << "#\n");

    if(rgtoken=="isbackbone"){
      rgid.setBackbone(true);
      continue;
    }else if(rgtoken=="israil"){
      rgid.setRail(true);
      continue;
    }else if(rgtoken=="iscoverageequivalent"){
      rgid.setCoverageEquivalentRead(true);
      continue;
    }

    if(mafsplit.size()<3){
      MIRANOTIFY(Notify::FATAL,"Line " << mafline << "\nexpected at least 3 elements, found " << mafsplit.size() << endl);
    }

    string & rgval1=mafsplit[2];
    CEBUG("rgval1 #" << rgval1 << "#\n");

    if(rgtoken=="name"){
      rgid.setGroupName(rgval1);
    }else if(rgtoken=="segmentnaming"
	     || rgtoken=="templatenaming"){
      //TODO:  implement templatenaming
    }else if(rgtoken=="ID"){
      int32 dummy=atoi(rgval1.c_str());
      if(dummy<0 || dummy >65535){
	MIRANOTIFY(Notify::FATAL,"Line @RG ID: id must be >=0 and <= 65535, but " << dummy << " was given.");
      }
      if(dummy>=externalidmapper.size()){
	externalidmapper.resize(dummy+1);
      }
      externalidmapper[dummy]=rgid;
    }else if(rgtoken=="technology"){
      rgid.setSequencingType(rgval1);
    }else if(rgtoken=="strainname"){
      rgid.setStrainName(rgval1);
    }else if(rgtoken=="segmentplacement"
	     || rgtoken=="templateplacement"){
      if(!rgid.setSegmentPlacement(rgval1)){
	MIRANOTIFY(Notify::FATAL,"Line @RG segmentplacement: did not recognise '" << rgval1 << "' as valid placement code.");
      }
    }else if(rgtoken=="templatesize"){
      int32 dummy=atoi(rgval1.c_str());
      rgid.setInsizeFrom(dummy);
      dummy=atoi(mafsplit[3].c_str());
      rgid.setInsizeTo(dummy);
    }else if(rgtoken=="machinetype"){
      rgid.setMachineType(rgval1);
    }else if(rgtoken=="basecaller"){
      rgid.setBaseCaller(rgval1);
    }else if(rgtoken=="dye"){
      rgid.setDye(rgval1);
    }else if(rgtoken=="primer"){
      rgid.setPrimer(rgval1);
    }else if(rgtoken=="clonevecname"){
      rgid.setCloneVecName(rgval1);
    }else if(rgtoken=="seqvecname"){
      rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="adaptorleft"){
//    rgid.set(rgval1);
    }else if(rgtoken=="adaptorright"){
//    rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="adaptorsplit"){
//    rgid.setSeqVecName(rgval1);
    }else if(rgtoken=="datadir"){
      rgid.setDataDir(rgval1);
    }else if(rgtoken=="datafile"){
      rgid.setDataFile(rgval1);
    }else{
      MIRANOTIFY(Notify::FATAL,"For line @RG: did not recognize token " << rgtoken);
    }
  }

  FUNCEND();
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////        Obsolete        ///////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

void MAFParse::parseTagDataV2(ifstream & mafin, string & acttoken, multitag_t & targettag)
{
  FUNCSTART("void MAFParse::parseTagData(ifstream & mafin, string & acttoken, multitag_t & tag)");

  static const auto anyof=boost::is_any_of("\t");

  targettag.commentisgff3=true;

  MAF_tmp_strv.clear();
  char nextchar;
  mafin.get(nextchar); // eat away the following \t
  getline(mafin,MAF_tmp_str);

  boost::split(MAF_tmp_strv, MAF_tmp_str, anyof);

  if(MAF_tmp_strv.size() < 4 || MAF_tmp_strv.size() > 7){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found " << MAF_tmp_strv.size());
  }


  if(AnnotationMappings::isMIRAEntry(MAF_tmp_strv[0])){
    targettag.setIdentifierStr(MAF_tmp_strv[0]);
  }else{
    const string * soident=&AnnotationMappings::translateGAP4feat2SOfeat(MAF_tmp_strv[0]);
    if(soident->empty()){
      soident=&AnnotationMappings::translateXGAP4feat2SOfeat(MAF_tmp_strv[0]);
    }
    if(soident->empty()){
      targettag.setIdentifierStr(MAF_tmp_strv[0]);
    }else{
      targettag.setIdentifierStr(*soident);
    }
  }

  targettag.from=atoi(MAF_tmp_strv[1].c_str());
  targettag.to=atoi(MAF_tmp_strv[2].c_str());

  if(targettag.from<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.from << " is <1, not allowed.");
  }
  if(targettag.to<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.to << " is <1, not allowed.");
  }

  targettag.from-=1;
  targettag.to-=1;

  if(MAF_tmp_strv[3].empty()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": the entry for the strand is empty, not allowed.");
  }

  targettag.setStrand(MAF_tmp_strv[3].front());

  if(MAF_tmp_strv.size() >=5){
    targettag.setSourceStr(MAF_tmp_strv[4]);

    if(MAF_tmp_strv.size() >=6){
      if(!MAF_tmp_strv[5].empty()){
	switch(MAF_tmp_strv[5].front()){
	case '.':{
	  targettag.phase=3;
	  break;
	}
	case '1':{
	  targettag.phase=2;
	  break;
	}
	case '2':{
	  targettag.phase=1;
	  break;
	}
	case '3':{
	  targettag.phase=0;
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << "\t" << MAF_tmp_str << "\nthe entry for the strand (" << MAF_tmp_strv[5] << ") is not 1, 2, 3 or .");
	}
	}
      }

      if(MAF_tmp_strv.size() ==7){
	targettag.setCommentStr(MAF_tmp_strv[6]);
      }
    }
  }

  FUNCEND();
}



#include<boost/tokenizer.hpp>

void MAFParse::parseTagDataV2(ifstream & mafin, string & acttoken, multitag_t & targettag)
{
  FUNCSTART("void MAFParse::parseTagData(ifstream & mafin, string & acttoken, multitag_t & tag)");

  targettag.commentisgff3=true;

  {
    char nextchar;
    mafin.get(nextchar); // eat away the following \t
  }
  getline(mafin,MAF_tmp_str);

  static const boost::char_separator<char> sep("\t", "", boost::keep_empty_tokens);
  boost::tokenizer<boost::char_separator<char>> tok(MAF_tmp_str,sep);

  auto tokI=tok.begin();

  if(tokI==tok.end()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found none");
  }


  if(AnnotationMappings::isMIRAEntry(*tokI)){
    targettag.setIdentifierStr(*tokI);
  }else{
    const string * soident=&AnnotationMappings::translateGAP4feat2SOfeat(*tokI);
    if(soident->empty()){
      soident=&AnnotationMappings::translateXGAP4feat2SOfeat(*tokI);
    }
    if(soident->empty()){
      targettag.setIdentifierStr(*tokI);
    }else{
      targettag.setIdentifierStr(*soident);
    }
  }

  if(++tokI==tok.end()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found one");
  }
  targettag.from=atoi(tokI->c_str());
  if(++tokI==tok.end()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found two");
  }
  targettag.to=atoi(tokI->c_str());

  if(targettag.from<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.from << " is <1, not allowed.");
  }
  if(targettag.to<1){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": (" << targettag.from << " " << targettag.to << ") -> " << targettag.to << " is <1, not allowed.");
  }

  targettag.from-=1;
  targettag.to-=1;

  if(++tokI==tok.end()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " for " << acttoken << "\t" << MAF_tmp_str << "\nExpected at between 4 and 7 tab delimited values, but found three");
  }
  if(tokI->empty()){
    MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << ": the entry for the strand is empty, not allowed.");
  }

  targettag.setStrand(tokI->front());

  if(++tokI!=tok.end()){
    targettag.setSourceStr(*tokI);

    if(++tokI!=tok.end()){
      if(!tokI->empty()){
	switch(tokI->front()){
	case '.':{
	  targettag.phase=3;
	  break;
	}
	case '1':{
	  targettag.phase=2;
	  break;
	}
	case '2':{
	  targettag.phase=1;
	  break;
	}
	case '3':{
	  targettag.phase=0;
	  break;
	}
	default : {
	  MIRANOTIFY(Notify::FATAL, "Error in " << MAF_read_name << " in targettag line " << acttoken << "\t" << MAF_tmp_str << "\nthe entry for the strand (" << *tokI << ") is not 1, 2, 3 or .");
	}
	}
      }

      if(++tokI!=tok.end()){
	targettag.setCommentStr(*tokI);
      }
      // hmmm, cannot really check whether ++tokI is end now, tokenizer throws :-(
    }
  }

  FUNCEND();
}

*/
