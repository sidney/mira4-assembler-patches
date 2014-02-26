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



#include "modules/mod_convert.H"

#include <boost/lexical_cast.hpp>

#include "io/generalio.H"
#include "caf/caf.H"
#include "mira/maf_parse.H"

#include "modules/misc.H"
#include "version.H"



using namespace std;


vector<MIRAParameters> ConvPro::CP_Pv;

string ConvPro::CP_fromtype;
list<string> ConvPro::CP_totype;
list<ofstream *> ConvPro::CP_ofs;

string ConvPro::CP_infile;
string ConvPro::CP_outbasename;

string ConvPro::CP_renamesequences;
string ConvPro::CP_renamenamescheme;

bool ConvPro::CP_splitcontigs2singlefiles=false;
bool ConvPro::CP_deletestaronlycolumns=false;
bool ConvPro::CP_blinddata=false;
bool ConvPro::CP_fillholesinstraingenomes=false;
bool ConvPro::CP_makecontigs=false;
bool ConvPro::CP_extractreadsinsteadcontigs=false;
bool ConvPro::CP_hardtrim=false;

string ConvPro::CP_namefile;
bool ConvPro::CP_sortbyname=false;
bool ConvPro::CP_keepnamesfromfile=true;

bool ConvPro::CP_mustdeletetargetfiles=true;

bool ConvPro::CP_specialtestcode=false;
bool ConvPro::CP_filter2readgroup=false;

base_quality_t ConvPro::CP_minqual=0;
bool ConvPro::CP_needsquality=true;
base_quality_t ConvPro::CP_defaultqual=30;

char ConvPro::CP_recalcconopt=' ';
char ConvPro::CP_recalcfeatureopt=' ';

uint32 ConvPro::CP_minbasecoverage=0;

uint32 ConvPro::CP_mincontiglength=0;
bool ConvPro::CP_minlengthisclipped=false;
uint32 ConvPro::CP_mincontigcoverage=1;
uint32 ConvPro::CP_minnumreads=0;

list<Contig> ConvPro::CP_clist;   // needed for CAF conversion (and GBF)
AssemblyInfo ConvPro::CP_assemblyinfo;


uint64 ConvPro::CP_readrenamecounter=1;

GFFSave ConvPro::CP_gffsave;
SAMCollect ConvPro::CP_samcollect;



ConvPro::~ConvPro()
{
  closeOpenStreams(CP_ofs);
}

void ConvPro::usage()
{
  cout << "miraconvert\t(MIRALIB version " << MIRAVERSION << ")\n"
    "Author:  Bastien Chevreux\t(bach@chevreux.org)\n"
    "Purpose: convert assembly and sequencing file types.\n\n";
  cout << "Usage:\n"
    "miraconvert [-f <fromtype>] [-t <totype> [-t <totype> ...]]\n"
    "\t[-aChimMsuZ]\n"
    "\t[-AcflnNoqrtvxXyz {...}]\n"
    "\t{infile} {outfile} [<totype> <totype> ...]\n\n";
  cout << "Options:\n";
  cout <<
    "\t-f <fromtype>\tload this type of project files, where fromtype is:\n"
    "\t   caf\t\t a complete assembly or single sequences from CAF\n"
    "\t   maf\t\t a complete assembly or single sequences from CAF\n"
    "\t   fasta\t sequences from a FASTA file\n"
    "\t   fastq\t sequences from a FASTQ file\n"
    "\t   gbf\t\t sequences from a GBF file\n"
    "\t   phd\t\t sequences from a PHD file\n"
    "\t   fofnexp\t sequences in EXP files from file of filenames\n";
  cout << "\t-t <totype>\twrite the sequences/assembly to this type (multiple\n"
    "\t\t\tmentions of -t are allowed):\n"
    "\t   ace\t\t sequences or complete assembly to ACE\n"
    "\t   caf\t\t sequences or complete assembly to CAF\n"
    "\t   maf\t\t sequences or complete assembly to MAF\n"
    "\t   sam\t\t complete assembly to SAM\n"
    "\t   samnbb\t like above, but leaving out reference (backbones) in mapping assemblies\n"
    "\t   gbf\t\t sequences or consensus to GBF\n"
    "\t   gff3\t\t consensus to GFF3\n"
    "\t   wig\t\t assembly coverage info to wiggle file\n"
    "\t   gcwig\t assembly gc content info to wiggle file\n"
    "\t   fasta\t sequences or consensus to FASTA file (qualities to\n"
    "\t\t\t  .qual)\n"
    "\t   fastq\t sequences or consensus to FASTQ file\n"
    "\t   exp\t\t sequences or complete assembly to EXP files in\n"
    "\t\t\t  directories. Complete assemblies are suited for gap4\n"
    "\t\t\t  import as directed assembly.\n"
    "\t\t\t  Note: using caf2gap to import into gap4 is recommended though\n"
    "\t   text\t\t complete assembly to text alignment (only when -f is\n"
    "\t\t\t  caf, maf or gbf)\n"
    "\t   html\t\t complete assembly to HTML (only when -f is caf, maf or\n"
    "\t\t\t  gbf)\n"
    "\t   tcs\t\t complete assembly to tcs\n"
    "\t   hsnp\t\t surrounding of SNP tags (SROc, SAOc, SIOc) to HTML\n"
    "\t\t\t (only when -f is caf, maf or gbf)\n"
    "\t   asnp\t\t analysis of SNP tags\n"
    "\t\t\t (only when -f is caf, maf or gbf)\n"
    "\t   cstats\t contig statistics file like from MIRA\n"
    "\t\t\t (only when source contains contigs)\n"
    "\t   crlist\t contig read list file like from MIRA\n"
    "\t\t\t (only when source contains contigs)\n"
    "\t   maskedfasta\t reads where sequencing vector is masked out\n"
    "\t\t\t (with X) to FASTA file (qualities to .qual)\n"
    "\t   scaf\t\t sequences or complete assembly to single sequences CAF\n";

  cout << "\t-a\t\tAppend to target files instead of rewriting\n";

  cout <<
    "\t-A <string>\tString with MIRA parameters to be parsed\n"
    "\t\t\t Useful when setting parameters affecting consensus\n"
    "\t\t\t calling like -CO:mrpg etc.\n"
    "\t\t\t E.g.: -a \"454_SETTINGS -CO:mrpg=3\"\n";

  cout <<
    "\t-b\t\tBlind data\n"
    "\t\t\t Replaces all bases in reads/contigs with a 'c'\n";

  cout << "\t-C\t\tPerform hard clip to reads\n"
    "\t\t\t When reading formats which define clipping points, will\n"
    "\t\t\t  save only the unclipped part into the result file.\n"
    "\t\t\t Applies only to files/formats which do not contain\n"
    "\t\t\t  contigs.\n";

  cout <<
    "\t-d\t\tDelete gap only columns\n"
    "\t\t\t When output is contigs: delete columns that are\n"
    "\t\t\t  entirely gaps (like after having deleted reads during\n"
    "\t\t\t  editing in gap4 or similar)\n"
    "\t\t\t When output is reads: delete gaps in reads\n";

  cout <<
    "\t-F\t\tFilter to read groups\n"
    "\t\t\t Special use case, do not use yet.\n";


  cout <<
    "\t-m\t\tMake contigs (only for -t = caf or maf)\n"
    "\t\t\t Encase single reads as contig singlets into the CAF/MAF\n"
    "\t\t\t file.\n";
  cout <<
    "\t-n <filename>\twhen given, selects only reads or contigs given by\n"
    "\t\t\t name in that file.\n";
  cout <<
    "\t-i\t\twhen -n is used, inverts the selection\n";
  cout <<
    "\t-o\t\tfastq quality Offset (only for -f = 'fastq')\n"
    "\t\t\t Offset of quality values in FASTQ file. Default of 0\n"
    "\t\t\t tries to automatically recognise.\n";

  cout <<
    "\t-Q <quality>\tSet default quality for bases in file types without quality values\n"
    "\t\t\t Furthermore, do not stop if expected quality files are missing (e.g. '.fasta')\n";

  cout <<
    "\t-R <name>\tRename contigs/singlets/reads with given name string\n"
    "\t\t\t to which a counter is appended.\n"
    "\t\t\t Known bug: will create duplicate names if input\n"
    "\t\t\t  contains contigs/singlets as well as free reads, i.e.\n"
    "\t\t\t  reads not in contigs nor singlets.\n";

  cout <<
    "\t-S <name>\t(name)Scheme for renaming reads, important for paired-ends\n"
    "\t\t\t Only 'solexa' is currently supported.\n";


  // TODO: re-adapt these switches to >2.9.8
//  cout << "\t-s <filename>\twhen loading assemblies from files that do not contain\n";
//  cout << "\t\t\t strain information (e.g. CAF), load the strain\n";
//  cout << "\t\t\t information from this file. (2 columns, tab delimited:\n";
//  cout << "\t\t\t readname left, strain name right)\n";
//

//

  cout <<
    "\n\t--------------------------------------------------------\n"
    "\tThe following switches work only when input (CAF or MAF)\n"
    "\tcontains contigs. Beware: CAF and MAf can also contain\n"
    "\tjust reads.\n"
    "\t--------------------------------------------------------\n\n";

  // TODO: check if ok for >2.9.8
  cout <<
    "\t-M\t\tDo not extract contigs (or their consensus), but the\n"
    "\t\t\t  sequence of the reads they are composed of.\n";
  cout <<
    "\t-N <filename>\tlike -n, but sorts output according to order given\n"
    "\t\t\t in file.\n";
  cout <<
    "\t-r [cCqf]\tRecalculate consensus and / or consensus quality values\n"
    "\t\t\t and / or SNP feature tags.\n"
    "\t\t\t 'c' recalc cons & cons qualities (with IUPAC)\n"
    "\t\t\t 'C' recalc cons & cons qualities (forcing non-IUPAC)\n"
    "\t\t\t 'q' recalc consensus qualities only\n"
    "\t\t\t 'f' recalc SNP features\n"
    "\t\t\t Note: only the last of cCq is relevant, f works as a\n"
    "\t\t\t  switch and can be combined with cQq (e.g. \"-r C -r f\")\n"
    "\t\t\t Note: if the CAF/MAF contains multiple strains,\n"
    "\t\t\t recalculation of cons & cons qualities is forced, you\n"
    "\t\t\t  can just influence whether IUPACs are used or not.\n";
  cout <<
    "\t-s\t\tsplit output into multiple files instead of creating a\n"
    "\t\t\t single file\n";
  cout <<
    "\t-u\t\t'fillUp strain genomes'\n"
    "\t\t\t Fill holes in the genome of one strain (N or @)\n"
    "\t\t\t with sequence from a consensus of other strains\n"
    "\t\t\t Takes effect only with -r and -t gbf or fasta/q\n"
    "\t\t\t in FASTA/Q: bases filled up are in lower case\n"
    "\t\t\t in GBF: bases filled up are in upper case\n";

  cout <<
    "\t-q <integer>\tDefines minimum quality a consensus base of a strain\n"
    "\t\t\t must have, consensus bases below this will be 'N'\n"
    "\t\t\t Default: 0\n"
    "\t\t\t Only used with -r, and -f is caf/maf and -t is (fasta\n"
    "\t\t\t  or gbf)\n";
//  cout <<
//    "\t-v <integer>\tDefines minimum coverage a consensus base of a strain\n"
//    "\t\t\t must have, bases with coverage below this will be 'N'\n"
//    "\t\t\t Default: 0\n"
//    "\t\t\t Only used with -r, and -t is (fasta\n"
//    "\t\t\t  or gbf)\n";

  cout <<
    "\t-v\t\tPrint version number and exit\n";

  cout <<
    "\t-x <integer>\tMinimum contig or read length\n"
    "\t\t\t When loading, discard all contigs / reads with a\n"
    "\t\t\t length less than this value. Default: 0 (=switched off)\n"
    "\t\t\t Note: not applied to reads in contigs!\n";
  cout <<
    "\t-X <integer>\tSimilar to -x but applies only to reads and\n"
    "\t\t\t then to the clipped length.\n";

  cout <<
    "\t-y <integer>\tMinimum average contig coverage\n"
    "\t\t\t When loading, discard all contigs with an\n"
    "\t\t\t average coverage less than this value.\n"
    "\t\t\t Default: 1\n";

  cout <<
    "\t-z <integer>\tMinimum number of reads in contig\n"
    "\t\t\t When loading, discard all contigs with a\n"
    "\t\t\t number of reads less than this value.\n"
    "\t\t\t Default: 0 (=switched off)\n";


  cout <<
    "\t-l <integer>\twhen output as text or HTML: number of bases shown in\n"
    "\t\t\t one alignment line. Default: 60.\n"
    "\t-c <character>\twhen output as text or HTML: character used to pad\n"
    "\t\t\t endgaps. Default: ' ' (blank)\n";

  cout << "\nExamples:\n"
    "\tmiraconvert source.maf dest.sam\n"
    "\tmiraconvert source.caf dest.fasta wig ace\n"
    "\tmiraconvert -x 2000 -y 10 source.caf dest.caf\n";
}


bool ConvPro::checkForFromType(const string & ftype)
{
  static set<string> ftypes={
    "caf",
    "maf",
    "phd",
    "gbf",
    "gbff",
    "gbfk",
    "gff3",
    "fasta",
    "fastq",
    "fofnexp"
  };
  return ftypes.find(ftype)!=ftypes.end();
}

bool ConvPro::checkForToType(const string & ttype)
{
  static set<string> ttypes={
    "fasta",
    "fastq",
    "maskedfasta",
    "caf",
    "maf",
    "sam",
    "samnbb",
    "ace",
    "scaf",
    "exp",
    "gbf",
    "gff3",
    "tcs",
    "text",
    "txt",
    "html",
    "wiggle",
    "wig",
    "gcwiggle",
    "gcwig",
    "asnp",
    "fcov",
    "hsnp",
    "cstats",
    "crlist",
    "null"
  };
  return ttypes.find(ttype)!=ttypes.end();
}

void ConvPro::checkTypes(const string & fromtype,list<string> & totype)
{
  if(!checkForFromType(fromtype)){
    usage();
    cout << endl;
    cerr << "Unknown or illegal file type '" << fromtype << "' defined as <fromtype>\n";
    exit(1);
  }
  if(CP_totype.empty()){
    CP_totype.push_back(fromtype);
  }
  for(list<string>::iterator ttI= CP_totype.begin(); ttI!=CP_totype.end(); ++ttI){
    if(!checkForToType(*ttI)){
      usage();
      cout << endl;
      cerr << "ConvPro::checkTypes: Unknown or illegal file type '" << *ttI << "' defined as <totype>\n";
      exit(1);
    }
  }
}

// comfort function: parse fromtype and totype from name of infile, outfile
void ConvPro::guessFromAndToType(const string & fnamefrom, string & fromtype, string * fromstem, const string & fnameto, list<string> & totypes, string * tostem)
{
  uint8 ziptype=0;
  string ft;
  string dummyfromstem;
  string dummypathto;
  guessFileAndZipType(fnamefrom,dummypathto,dummyfromstem,ft,ziptype);
  if(fromtype.empty()){
    fromtype.swap(ft);
  }
  if(fromstem != nullptr) {
    fromstem->swap(dummypathto);
    if(!dummypathto.empty() && dummypathto!="/"){
      *fromstem+='/';
    }
    *fromstem+=dummyfromstem;
  }

  string dummytostem;
  ft.clear();
  guessFileAndZipType(fnameto,dummypathto,dummytostem,ft,ziptype);
  if(!ft.empty() && checkForToType(ft)){
    totypes.push_back(ft);
    if(dummytostem.empty()) dummytostem=dummyfromstem;
    if(tostem != nullptr) {
      *tostem=dummypathto;
      if(!dummypathto.empty() && dummypathto!="/"){
	*tostem+='/';
      }
      *tostem+=dummytostem;
    }
  }
  //cout << "dfs: " << dummyfromstem << " dpt: " << dummypathto << " dts: " << dummytostem << endl;
}


void ConvPro::filterToReadGroup(ReadPool & rp)
{
  FUNCSTART("void ConvPro::filterToReadGroup(ReadPool & rp)");

  vector<ofstream> ofspassed(ReadGroupLib::getNumReadGroups());
  vector<ofstream> ofswidow(ReadGroupLib::getNumReadGroups());
  vector<ofstream> ofsdebris(ReadGroupLib::getNumReadGroups());

  deque<string> seenrgnames;
  for(uint32 rglid=1; rglid<ReadGroupLib::getNumReadGroups(); ++rglid){
    auto rgid=ReadGroupLib::getReadGroupID(rglid);
    string rgname(rgid.getGroupName());
    if(rgname.empty()){
      rgname=string("readgroup")+boost::lexical_cast<string>(rglid)+"_passed.fastq";
    }
    for(auto & sn : seenrgnames){
      if(sn==rgname){
	cout << "\nReadgroup '"<<sn<<"' seen more than once in MAF file, this is illegal.\n";
	exit(0);
      }
    }
    seenrgnames.push_back(rgname);
    string pname(rgname+"_passed.fastq");
    openOFStream(ofspassed[rglid],pname,ios::out);
    pname=rgname+"_widow.fastq";
    openOFStream(ofswidow[rglid],pname,ios::out);
    pname=rgname+"_debris.fastq";
    openOFStream(ofsdebris[rglid],pname,ios::out);
  }

  Read::setCoutType(Read::AS_FASTQ);

  rp.makeTemplateIDs(true);
  {
    vector<uint32> sortdummy;
    rp.sortPoolToMIRAStandard(sortdummy);
  }
  vector<bool> hasmultisegmentemplates(ReadGroupLib::getNumReadGroups(),false);
  for(uint32 rpi=0; rpi< rp.size(); ++rpi){
    if(rp[rpi].getTemplatePartnerID()!=-1){
      hasmultisegmentemplates[rp[rpi].getReadGroupID().getLibId()]=true;
    }
  }
  hasmultisegmentemplates[0]=false;

  for(uint32 rpi=0; rpi< rp.size(); ){
    bool allok=true;
    auto lasti=rpi;
    for(; lasti<rp.size() && rp[rpi].getTemplateID()==rp[lasti].getTemplateID(); ++lasti){
      if(CP_mincontiglength>0 && rp[lasti].getLenSeq()<CP_mincontiglength) allok=false;
    }
    auto numseq=lasti-rpi;
    if(hasmultisegmentemplates[rp[rpi].getReadGroupID().getLibId()] && numseq<2) allok=false;
    for(;rpi<lasti; ++rpi){
      if(allok){
	ofspassed[rp[rpi].getReadGroupID().getLibId()] << rp[rpi];
      }else if(numseq==1){
	ofsdebris[rp[rpi].getReadGroupID().getLibId()] << rp[rpi];
      }else{
	if(CP_mincontiglength>0 && rp[rpi].getLenSeq()>=CP_mincontiglength){
	  ofswidow[rp[rpi].getReadGroupID().getLibId()] << rp[rpi];
	}else{
	  ofsdebris[rp[rpi].getReadGroupID().getLibId()] << rp[rpi];
	}
      }
    }
  }

}


#define CEBUG(bla)   //{cout << bla; cout.flush();}
void ConvPro::specialTestCode(list<Contig> & clist, ReadPool & rp)
{
  FUNCSTART("void ConvPro::specialTestCode(list<Contig> & clist, ReadPool & rp)");
  try {
    for(auto cI=clist.begin(); cI!=clist.end(); ++cI){
      cI->findBestNonMisassembledRange();
      cout << "CHecking: " << cI->getContigName() << endl;
      if(cI->getContigLength()>=1000){
	auto p=cI->findBestPairConsistencyRange();
	if(p.first>=0) cout << "Found best range: " << cI->getContigName() << "\t" << p.first << "\t" << p.second << "\t" << p.second-p.first << endl;
      }


      if(1){
	auto range=cI->findBestNonMisassembledRange();
	if(range.first>=0){
	  cout<<"Found misassembly by repeat marker. Best range: " << range.first << ".." << range.second << '\t' << range.second-range.first << endl;
	  cI->trimContigToRange(range.first,range.second);

	}
      }

      if(1){
      	coverageinfo_t cinfo;
      	vector<uint64> covvals;
      	cI->collectCoverage(covvals);
      	cI->calcStatsOnContainer(cinfo,covvals);
      	cout << "1st covnum: " << cinfo << endl;

      	// TODO: perhaps make this dependend of ratio mean vs stddev ?
      	cI->calcSecondOrderStatsOnContainer(cinfo,covvals);
      	cout << "2nd covnum: " << cinfo << endl;

      	vector<uint8> peakindicator;
      	cI->findPeaks(cinfo.median,peakindicator);
      	unordered_set<readid_t> readsremoved;
      	cI->reduceReadsAtCoveragePeaks(cinfo.median,peakindicator,readsremoved);
      	cout << "Removed " << readsremoved.size() << endl;
      }


//      vector<bool> readsmarkedsrm;
//      if(1){
//	cout << "\nMarking tricky 454 / Solexa overcalls in temporary contig.\n";
//	cout << "Marked " << I->editTrickyOvercalls(true,false,readsmarkedsrm) << " reads.\n";
//	bool newreptmarked=Assembly::markRepeats(*I, readsmarkedsrm);
//	cout << "Edited " << I->editTrickyOvercalls(false,false,readsmarkedsrm) << " reads.\n";
//	I->deleteStarOnlyColumns(0, I->getContigLength()-1);
//	I->deleteTagsInReads(Read::REA_defaulttag_PSHP.identifier);
//      }
//      I->markFeaturesByConsensus(true,true,true);
//      I->editSingleDiscrepancyNoHAFTag(readsmarkedsrm);
//      assout::saveReadTagList(*I,
//			      "blabla",
//			      true);


//      uint32 numcoledits=0;
//      uint32 numreadedits=0;
//      for(uint32 i=0; i<5; ++i){
//	I->editPBSledgeHammer(readsmarkedsrm,numcoledits,numreadedits);
//	cout <<"### i " << i << " " << numcoledits << " " << numreadedits << endl;
//	I->deleteStarOnlyColumns(0,I->getContigLength());
//	//I->upDownCase(5);
//      }

    }
  }
//  try{
//    for(auto & cle : clist){
//      auto & cr=cle.getContigReads();
//      for(auto & cre : cr){
//	Read & r=const_cast<Read &>(cre);
//	if(r.isBackbone()){
//	  cout << "Sorting " << r.getName() << endl;
//	  r.sortTagsForGFF3();
//	}
//      }
//    }
//  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
  FUNCEND();
}
#undef CEBUG



void ConvPro::putReadsInContigsAndSave(vector<MIRAParameters> & Pv, ReadPool & rp)
{
  for(uint32 i=0; i<rp.size(); i++) {
    if(!rp[i].hasQuality()){
      base_quality_t bq=30;
      if(rp[i].getReadGroupID().getDefaultQual()<=100) bq=rp[i].getReadGroupID().getDefaultQual();
      rp[i].setQualities(bq);
    }

    Contig con(&Pv, rp);
    CP_clist.push_back(con);
    CP_clist.back().addFirstRead(i,1);
    CP_clist.back().setContigName(rp[i].getName()+"_contig");
    saveContigList_helper(CP_clist, rp);
    CP_clist.clear();
  }
}

void ConvPro::discardShortReads(vector<MIRAParameters> & Pv, ReadPool & rp, uint32 minlength, bool fromclipped)
{
  for(uint32 i=0; i<rp.size(); i++) {
    uint32 len;
    if(fromclipped){
      len=rp[i].getLenClippedSeq();
    }else{
      len=rp[i].getLenSeq();
    }
    if(len<minlength) rp[i].discard();
  }
}



//string ConvPro::createFileNameFromBasePostfixContigAndRead(const string & basename, string & postfix, Contig * actcon, Read * actread)
string ConvPro::createFileNameFromBasePostfixContigAndRead(const string & basename, char * postfix, Contig * actcon, Read * actread)
{
  string filename=basename;
  if(actcon != nullptr){
    if(!filename.empty()) filename+='_';
    filename+=actcon->getContigName();
  }else if(actread != nullptr){
    if(!filename.empty()) filename+='_';
    filename+=actread->getName();
  }
  filename+=postfix;
  return filename;
}


bool ConvPro::contig__nameordercomp(const Contig & a, const Contig & b)
{
  return General::getNameOrder(a.getContigName()) < General::getNameOrder(b.getContigName());
}


void ConvPro::sortContigsByName(list<Contig> & clist)
{
  clist.sort(contig__nameordercomp);
}


void ConvPro::sortPoolByName(ReadPool & rp, string & filename)
{
  FUNCSTART("void ConvPro::sortPoolByName(ReadPool & rp, string & filename)");

  cout << "Sorting pool ..."; cout.flush();

  rp.allowNameIndex(true);

  ifstream fin;
  openIFStream(fin,filename, ios::in);

  vector<uint32> newsortorder;

  string elemname, dummy;
  uint32 numread=0;
  while(GeneralIO::readKeyValue(fin, elemname,dummy)){
    int32 newpos=rp.getReadIndex(elemname);
    if(newpos>=0) newsortorder.push_back(static_cast<uint32>(newpos));
  }
  fin.close();

  rp.sortPoolToGivenOrder(newsortorder);

  rp.allowNameIndex(false);

  cout << "done.\n";

  FUNCEND();
}

void ConvPro::saveContigList_helper(list<Contig> & clist, ReadPool & rp)
{
  FUNCSTART("void ConvPro::saveContigList_helper(list<Contig> & clist, ReadPool & rp)");

  if(CP_specialtestcode) specialTestCode(clist,rp);

  //{
  //  cout << "CLISTSIZE: " << clist.size() << endl;
  //  list<Contig>::iterator cI=clist.begin();
  //  for(; cI != clist.end(); cI++){
  //    cout << "cname: " << cI->getContigName() << endl;
  //  }
  //}

  BUGIFTHROW(!CP_ofs.empty() && CP_ofs.size() != CP_totype.size(), "Ooops? !CP_ofs.empty() && CP_ofs.size() != CP_totype.size() ???");

  list<ofstream *>::iterator ofsI= CP_ofs.begin();
  list<string>::iterator ttI= CP_totype.begin();
  for(; ttI!=CP_totype.end(); ++ttI, ++ofsI){
    if(*ttI=="null"){
      // do nothing
    }else if(*ttI=="scaf"){
      //clear_conandrp=false;
    }else if(*ttI=="hsnp"){
      MIRAParameters::generateProjectOutNames(CP_Pv,CP_outbasename);
      string fn;
      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_snpenvironment,
	  ".html",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_snpenvironment,
	  ".html");
      }
      assout::saveSNPSurroundingAsHTML(clist,fn,CP_mustdeletetargetfiles);
    }else if(*ttI=="cstats"){
      MIRAParameters::generateProjectOutNames(CP_Pv,CP_outbasename);
      string fn;
      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_contigstats,
	  ".txt",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_contigstats,
	  ".txt");
      }
      assout::saveStatistics(clist,fn,CP_mustdeletetargetfiles);
      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_contigtags,
	  ".txt",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_contigtags,
	  ".txt");
      }
      assout::saveConsensusTagList(clist,fn,CP_mustdeletetargetfiles);
    }else if(*ttI=="crlist"){
      MIRAParameters::generateProjectOutNames(CP_Pv,CP_outbasename);
      string fn;
      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_crlist,
	  ".txt",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_crlist,
	  ".txt");
      }
      assout::saveContigReadList(clist,fn,CP_mustdeletetargetfiles);
    }else if(*ttI=="asnp"){
      MIRAParameters::generateProjectOutNames(CP_Pv,CP_outbasename);
      string fn,fa,fs,fc;

      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_snpanalysis,
	  ".txt",
	  &clist.front());
	fa=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featureanalysis,
	  ".txt",
	  &clist.front());
	fs=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featuresummary,
	  ".txt",
	  &clist.front());
	fc=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featuresequences,
	  ".txt",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_snpanalysis,
	  ".txt");
	fa=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featureanalysis,
	  ".txt");
	fs=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featuresummary,
	  ".txt");
	fc=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featuresequences,
	  ".txt");
      }

      assout::saveSNPList(clist,fn,CP_mustdeletetargetfiles);
      assout::saveFeatureAnalysis(clist,rp,
				  fa,fs,fc,
				  CP_mustdeletetargetfiles);

    }else if(*ttI=="fcov"){
      MIRAParameters::generateProjectOutNames(CP_Pv,CP_outbasename);
      string fn,fa,fs,fc;

      if(CP_splitcontigs2singlefiles){
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featurecoverage,
//	  "coveragei",
	  ".txt",
	  &clist.front());
      }else{
	fn=createFileNameFromBasePostfixContigAndRead(
	  CP_Pv[0].getAssemblyParams().as_outfile_stats_featurecoverage,
//	  "coveragei",
	  ".txt");
      }

      assout::saveCoverageInfo(clist,fn,CP_mustdeletetargetfiles);
    }else if(*ttI=="fasta"){
      //CALLGRIND_START_INSTRUMENTATION;
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "");
      }
      assout::saveStrainsAsFASTAQ(clist,
				  rp,
				  bn,
				  false,
				  CP_minbasecoverage,
				  CP_minqual,
				  CP_mustdeletetargetfiles,
				  CP_fillholesinstraingenomes);
    }else if(*ttI=="fastaqual"){
      // fastaqual is "do-nothing" as "fasta" also write fastaqual here!
    }else if(*ttI=="fastq"){
      //CALLGRIND_START_INSTRUMENTATION;
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "");
      }
      assout::saveStrainsAsFASTAQ(clist,
				  rp,
				  bn,
				  true,
				  CP_minbasecoverage,
				  CP_minqual,
				  CP_mustdeletetargetfiles,
				  CP_fillholesinstraingenomes);
    } else if(*ttI=="caf"){
      Contig::setCoutType(Contig::AS_CAF);
      list<Contig>::iterator I=clist.begin();
      for(;I!=clist.end(); I++){
	if(CP_recalcfeatureopt=='f') I->markFeaturesByConsensus(true,true,true);
	if(CP_recalcfeatureopt=='r') {
	  vector<bool> dummy;
	  Contig::repeatmarker_stats_t contigrepstats;
	  Assembly::markRepeats(*I,dummy,contigrepstats);
	}
	bool mustclose=false;
	if(!(*ofsI)->is_open()){
	  string bn;
	  if(CP_splitcontigs2singlefiles){
	    bn=createFileNameFromBasePostfixContigAndRead(
	      CP_outbasename,
	      ".caf",
	      &clist.front());
	  }else{
	    bn=createFileNameFromBasePostfixContigAndRead(
	      CP_outbasename,
	      ".caf",
	      nullptr);
	  }
	  openOFStream(*(*ofsI),bn,ios::out);
	  mustclose=true;
	}
	*(*ofsI) << *I;
	if(mustclose){
	  (*ofsI)->close();
	}
      }
    } else if(*ttI=="sam"){
      BUGIFTHROW(!(*ofsI)->is_open(),"Ooops, SAM stream not open?");
      for(auto & cle : clist){
	cle.dumpAsSAM(*(*ofsI),CP_samcollect,true);
      }
    } else if(*ttI=="samnbb"){
      BUGIFTHROW(!(*ofsI)->is_open(),"Ooops, SAM stream not open?");
      for(auto & cle : clist){
	cle.dumpAsSAM(*(*ofsI),CP_samcollect,false);
      }
    } else if(*ttI=="maf"){
      Contig::setCoutType(Contig::AS_MAF);
      list<Contig>::iterator I=clist.begin();
      for(;I!=clist.end(); I++){
	if(CP_recalcfeatureopt=='f') I->markFeaturesByConsensus(true,true,true);
	if(CP_recalcfeatureopt=='r') {
	  vector<bool> dummy;
	  Contig::repeatmarker_stats_t contigrepstats;
	  Assembly::markRepeats(*I,dummy,contigrepstats);
	}
	if(!(*ofsI)->is_open()){
	  string bn;
	  if(CP_splitcontigs2singlefiles){
	    bn=createFileNameFromBasePostfixContigAndRead(
	      CP_outbasename,
	      ".maf",
	      &clist.front());
	  }else{
	    bn=createFileNameFromBasePostfixContigAndRead(
	      CP_outbasename,
	      ".maf",
	      nullptr);
	  }
	  openOFStream(*(*ofsI),bn,ios::out);
	  Contig::dumpMAF_Head(*(*ofsI));
	}
	*(*ofsI) << *I;
	if(CP_splitcontigs2singlefiles){
	  (*ofsI)->close();
	}
      }
    } else if(*ttI=="html"){
      //cerr << "HTML output currently deactivated in development version!\n";
      //exit(1);
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".html",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".html");
      }
      assout::dumpContigListAsHTML(clist, bn, CP_mustdeletetargetfiles, CP_outbasename);
    } else if(*ttI=="text"
	      || *ttI=="txt"){
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".txt",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".txt");
      }
      assout::saveAsTXT(clist, bn, CP_mustdeletetargetfiles);
    } else if(*ttI=="exp"){
      // outbasename is in this case a directory name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "");
      }
      assout::saveAsGAP4DA(clist,bn,CP_mustdeletetargetfiles);
    } else if(*ttI=="gbf"){
      // outbasename is in this case the basename name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "");
      }
      assout::saveStrainsAsGBF(clist,
			       rp,
			       bn,
			       CP_minqual,
			       CP_fillholesinstraingenomes,
			       CP_mustdeletetargetfiles);
    } else if(*ttI=="gff3"){
      // outbasename is in this case the basename name
      list<Contig>::iterator cI=clist.begin();
      for(; cI != clist.end(); ++cI){
	string bn;
	if(CP_splitcontigs2singlefiles){
	  bn=createFileNameFromBasePostfixContigAndRead(
	    CP_outbasename,
	    "",
	    &clist.front());
	  if(CP_gffsave.is_open()){
	    CP_gffsave.close();
	  }
	  CP_gffsave.open(bn);
	}else{
	  bn=createFileNameFromBasePostfixContigAndRead(
	    CP_outbasename,
	    "");
	  if(!CP_gffsave.is_open()){
	    CP_gffsave.open(bn);
	  }
	}
	CP_gffsave.acquireContig(*cI,rp);
      }
    } else if(*ttI=="ace"){
      // outbasename is in this case the basename name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".ace",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".ace");
      }
      assout::saveAsACE(clist,bn,CP_mustdeletetargetfiles);
    } else if(*ttI=="tcs"){
      // outbasename is in this case the basename name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".tcs",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".tcs");
      }
      assout::saveAsTCS(clist,bn,CP_mustdeletetargetfiles);
    } else if(*ttI=="wiggle" || *ttI=="wig"){
      // outbasename is in this case the basename name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".wig",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  ".wig");
      }
      assout::saveAsWiggle(clist,bn,CP_mustdeletetargetfiles,false);
    } else if(*ttI=="gcwiggle" || *ttI=="gcwig"){
      // outbasename is in this case the basename name
      string bn;
      if(CP_splitcontigs2singlefiles){
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "_gccontent.wig",
	  &clist.front());
      }else{
	bn=createFileNameFromBasePostfixContigAndRead(
	  CP_outbasename,
	  "_gccontent.wig");
      }
      assout::saveAsWiggle(clist,bn,CP_mustdeletetargetfiles,true);
    } else {
      cerr << "\n\n-t " << *ttI << " is not a valid 'to' type when converting contigs (sorry). But maybe something went wrong, please contact the author.\n";
      exit(1);
    }
  }

  if(!CP_splitcontigs2singlefiles){
    CP_mustdeletetargetfiles=false;
  }

  FUNCEND();
}

void ConvPro::saveContigList(list<Contig> & clist, ReadPool & rp)
{
  FUNCSTART("void ConvPro::saveContigList(list<Contig> & clist, ReadPool & rp)");
  bool dosomeoutput=false;

  list<Contig>::iterator cI=clist.begin();
  for(; cI != clist.end(); cI++){
    bool conout=true;

    if(CP_mincontiglength>0
       && cI->getContigLength() < CP_mincontiglength){
      conout=false;
    } else {
      Contig::constats_t constats(cI->getStats());

      //cI->stats(cout);

      if(CP_mincontigcoverage>0
	 && constats.avg_coverage < CP_mincontigcoverage){
	conout=false;
      } else if(CP_minnumreads>0
		&& constats.total_reads < CP_minnumreads){
	conout=false;
      }
    }

    if(conout){
      if(General::hasNames()){
	string cname(cI->getContigName());
	if(!General::checkNamePresence(cname)){
	  conout=false;
	}
	if(!CP_keepnamesfromfile) conout=!conout;
      }
    }

    // delete contigs which should not be output
    // TODO:
    //  would generally be better to have that in some loading callback (and
    //  would work for contigs well enough), but readpool mechanisms would
    //  not at the moment, too primitive
    if(!conout){
      cI=clist.erase(cI);
      if(cI != clist.begin()) --cI;
    }

    dosomeoutput|=conout;
  }

  if(dosomeoutput){
    for(cI=clist.begin(); cI != clist.end(); cI++){
      if(CP_deletestaronlycolumns) {
	cI->deleteStarOnlyColumns(0,cI->getContigLength());
      }
      if(CP_blinddata){
	cI->blindContig();
      }
    }

    Assembly::refreshContigAndReadpoolValuesAfterLoading(rp,clist);

    string dummyempty;
    // TODO: ! make autoconfigure: on several strains, this is needed!
    //  else, let user define via switch
    for(cI=clist.begin(); cI != clist.end(); cI++){
      if(CP_recalcconopt=='c'
	 || CP_recalcconopt=='C'){
	cI->trashConsensusCache(false);
      }
      if(CP_recalcconopt=='q'){
	cI->trashConsensusCache(true);
      }

      CP_assemblyinfo.storeContigStats(cI->getStats(), dummyempty);
    }

    try{
      saveContigList_helper(CP_clist, rp);
    }
    catch(Notify n){
      n.handleError(THISFUNC);
    }
  }
}

void ConvPro::saveReadPool(ReadPool & rp, list<ofstream *> & ofs)
{
  FUNCSTART("void ConvPro::saveReadPool(ReadPool & rp, list<ofstream *> & ofs)");

  rp.adjustIllegalQualities(30);

  if(CP_deletestaronlycolumns) {
    for(uint32 i=0; i<rp.size(); i++) {
      rp.getRead(i).removeGapsFromRead();
    }
  }
  if(CP_blinddata) {
    for(uint32 i=0; i<rp.size(); i++) {
      rp.getRead(i).blindSeqData('c');
    }
  }
  if(General::hasNames()){
    for(uint32 i=0; i<rp.size(); i++) {
      string rname(rp[i].getName());
      if((!General::checkNamePresence(rname) && CP_keepnamesfromfile)
	 || (General::checkNamePresence(rname) && !CP_keepnamesfromfile)){
	rp[i].discard();
      }
    }
  }
  if(CP_hardtrim){
    for(uint32 i=0; i<rp.size(); ++i) {
      rp[i].performHardTrim();
    }
  }
  if(CP_mincontiglength>0 && !CP_filter2readgroup){
    discardShortReads(CP_Pv,rp,CP_mincontiglength,CP_minlengthisclipped);
  }

  if(!CP_renamesequences.empty()){
    if(!CP_renamenamescheme.empty()){
      for(size_t i=1; i<ReadGroupLib::getNumReadGroups(); ++i){
	auto rgid=ReadGroupLib::getReadGroupID(i);
	if(!rgid.setSegmentNaming(CP_renamenamescheme)){
	  BUGIFTHROW(true,"Naming scheme '" << CP_renamenamescheme << "' is unknown.");
	}
      }
      for(uint32 i=0; i<rp.size(); ++i){
	rp[i].setTemplate("");
	rp[i].setTemplateSegment(0);
	rp[i].calcTemplateInfo();
      }
      rp.makeTemplateIDs();
    }

    vector<uint8> readrenamed(rp.size(),0);
    string tmpname;
    for(uint32 rpi=0; rpi<rp.size(); ++rpi){
      if(rp[rpi].getTemplatePartnerID()>=0 && readrenamed[rp[rpi].getTemplatePartnerID()]){
	tmpname=rp[rp[rpi].getTemplatePartnerID()].getTemplate();
      }else{
	tmpname=CP_renamesequences+"_"+boost::lexical_cast<std::string>(CP_readrenamecounter++);
      }
      auto segment=rp[rpi].getTemplateSegment();
      if(rp[rpi].getReadNamingScheme()==ReadGroupLib::SCHEME_SOLEXA){
	if(rp[rpi].getTemplateSegment()==1){
	  tmpname+="/1";
	}else if(rp[rpi].getTemplateSegment()==255){
	  tmpname+="/2";
	}else if(rp[rpi].getTemplateSegment()!=0){
	  BUGIFTHROW(true,"Read " << rp[rpi].getName() << " has illegal segment " << static_cast<uint16>(rp[rpi].getTemplateSegment()) << " ???");
	}
      }
      rp[rpi].setTemplate("");
      rp[rpi].setTemplateSegment(0);
      rp[rpi].setName(tmpname);
      readrenamed[rpi]=1;
    }
  }

  if(CP_makecontigs) {
    putReadsInContigsAndSave(CP_Pv, rp);
  }else if(CP_filter2readgroup){
    filterToReadGroup(rp);
  }else{
    list<string>::iterator ttI= CP_totype.begin();
    list<ofstream *>::iterator ofsI= ofs.begin();
    for(; ttI!=CP_totype.end(); ++ttI, ++ofsI){
      bool needofsIclose=false;
      // BaCh 16.01.2012: changed != below to == ... was probably a simple programming error (I hope)
      if(*ttI=="gff3"){
	if(!CP_splitcontigs2singlefiles){
	  BUGIFTHROW(!(*(*ofsI)).is_open(), *ttI << " file stream not open???");
	}
	openOFStream((*(*ofsI)),rp[0].getName()+'.'+*ttI,ios::out);
	needofsIclose=true;
      }
      if(*ttI=="fasta"){
	// double indirection because iterator needs one and it is a list of ofstream pointers ...
	rp.dumpAs(*(*ofsI),Read::AS_FASTA,false);
      } else if(*ttI=="fastaqual"){
	rp.dumpAs(*(*ofsI),Read::AS_FASTAQUAL,false);
      } else if(*ttI=="maskedfasta"){
	rp.dumpAs(*(*ofsI),Read::AS_MASKEDMASKFASTA,false);
      } else if(*ttI=="maskedfastaqual"){
	rp.dumpAs(*(*ofsI),Read::AS_MASKEDMASKFASTAQUAL,false);
      } else if(*ttI=="fastq"){
	rp.dumpAs(*(*ofsI),Read::AS_FASTQ,false);
      } else if(*ttI=="caf" || *ttI=="scaf" ){
	rp.dumpAs(*(*ofsI),Read::AS_CAF,false);
      } else if(*ttI=="maf"){
	//rp.dumpAs(*(*ofsI),Read::AS_MAF,false);
	rp.saveAsMAF(*(*ofsI),false);
      } else if(*ttI=="gff3"){
	for(uint32 rpi=0; rpi<rp.size(); ++rpi){
	  bool mustclose=false;
	  if(!CP_gffsave.is_open()){
	    mustclose=true;
	    CP_gffsave.open(rp[rpi].getName());
	  }
	  CP_gffsave.acquireRead(rp[rpi]);
	  if(mustclose) CP_gffsave.close();
	}
      } else {
	cout.flush();
	cerr << "\n\n-t " << *ttI << " is not a valid type for saving a readpool (internal)!\n";
	//usage();
	exit(1);
      }
      if(needofsIclose){
	(*(*ofsI)).close();
      }
    }
  }

  FUNCEND();
}

uint32 ConvPro::openOFSlist(Contig * optcontig, list<ofstream *> & ofs)
{
  FUNCSTART("uint32 ConvPro::openOFSlist(Contig * optcontig, list<ofstream *> & ofs)");
  BUGIFTHROW(CP_totype.empty(), " CP_totype.empty() ???");

  uint32 mustclose=0;
  ofstream * ofstmp;

  for(list<string>::iterator ttI= CP_totype.begin(); ttI!=CP_totype.end(); ++ttI){

    cout << "opening " << *ttI << endl;

    ofstmp=new ofstream;
    ofs.push_back(ofstmp);
    if(*ttI=="fasta"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="fastaqual"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta.qual",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="maskedfasta"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="maskedfastaqual"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta.qual",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="fastq"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fastq",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="caf" || *ttI=="scaf" ){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".caf",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="maf"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".maf",optcontig),ios::out);
      Contig::dumpMAF_Head(*(ofs.back()));
      ++mustclose;
    } else if(*ttI=="sam" || *ttI=="samnbb"){
      openOFStream(*ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".sam",optcontig),ios::out);
      ++mustclose;
    } else if(*ttI=="gff3"){
      if(CP_gffsave.is_open()) CP_gffsave.close();
      CP_gffsave.open(createFileNameFromBasePostfixContigAndRead(CP_outbasename,".gff3",optcontig));
    }
  }
  return mustclose;
}

void ConvPro::closeOFSList(uint32 howmany, list<ofstream *> & ofs)
{
  FUNCSTART("uint32 ConvPro::closeOFSList(uint32 howmany)");
  BUGIFTHROW(howmany>ofs.size(),"howmany>ofs.size() ???");
  for(uint32 i=0; i<howmany; ++i){
    delete ofs.back();
    ofs.pop_back();
  }
  FUNCEND();
}


void ConvPro::cafmafload_callback(list<Contig> & clist, ReadPool & rp)
{
  FUNCSTART("void ConvPro::cafmafload_callback(list<Contig> & clist, ReadPool & rp)");
  BUGIFTHROW(clist.empty() && rp.size()==0,"clist.empty() && rp.size()==0");
  {
    list<Contig>::iterator cI=clist.begin();
    for(; cI != clist.end(); ++cI){
      if(!CP_renamesequences.empty()){
	cI->setContigName("");
	cI->setContigNamePrefix(CP_renamesequences);
      }
    }
  }
  if(!clist.empty() && !CP_extractreadsinsteadcontigs){
    saveContigList(clist,rp);
  }else{
    list<ofstream *> ofs;
    uint32 mustclose=0;
    if(CP_splitcontigs2singlefiles && !clist.empty()){
      mustclose=openOFSlist(&clist.front(),ofs);
      saveReadPool(rp,ofs);
    }else{
      saveReadPool(rp,CP_ofs);
    }
    closeOFSList(mustclose,ofs);
  }

  Read::trashReadNameContainer();
  clist.clear();
  rp.discard();
}

void ConvPro::readpoolload_callback(ReadPool & rp)
{
  saveReadPool(rp,CP_ofs);

  Read::trashReadNameContainer();
  rp.discard();
}


void ConvPro::closeOpenStreams(list<ofstream *> & ofsl)
{
  list<ofstream *>::iterator ofsI= ofsl.begin();
  for(; ofsI!=ofsl.end(); ++ofsI){
    delete *ofsI;
  }
}

void ConvPro::openOFStream(std::ofstream & ofs, std::string fname, std::ios_base::openmode mode)
{
  FUNCSTART("void ConvPro::openOFStream(std::ofstream & ofs, std::string fname, std::ios_base::openmode mode)");
  ofs.open(fname,mode);
  if(ofs.fail()){
    MIRANOTIFY(Notify::FATAL,"File " << fname << " could not be opened for writing. Possible causes: non-existing or write protected directory; disk quota exceeded; others.");
  }
}

void ConvPro::openIFStream(std::ifstream & ifs, std::string fname, std::ios_base::openmode mode)
{
  FUNCSTART("void ConvPro::openOFStream(std::ofstream & ofs, std::string fname, std::ios_base::openmode mode)");
  ifs.open(fname,mode);
  if(ifs.fail()){
    MIRANOTIFY(Notify::FATAL,"File " << fname << " could not be opened for reading. Possible causes: non-existing directory; not allowed to read directory or file; others.");
  }
}


int ConvPro::mainConvPro(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int mainConvPro(int argc, char ** argv)");

  int c;
  extern char *optarg;
  extern int optind;

  //base_quality_t fqqualoffset=0;
  string fqqualoffset("33");

  string strainfile="";

  int32 linelen=60;
  char  endgap_fillchar=' ';

  string miraparams;

  //"CZihumMsl:r:c:f:t:s:q:n:N:v:x:X:y:z:o:a:"
  const char pstring[]=
    "abCdFhimMsuvZ"
    "A:c:f:l:n:N:o:q:Q:r:R:S:t:x:X:y:z:";

  while (1){
    c = getopt(argc, argv, pstring);
    if(c == -1) break;

    switch (c) {
    case 'v': {
      cout << MIRAVERSION << endl;
      exit(0);
    }
    case '?': {
      cerr << "Use '-h' to get a short online help." << endl;
      exit(100);
    }
    case 'h': {
      usage();
      exit(0);
    }

    case 'a': {
      CP_mustdeletetargetfiles=false;
      break;
    }
    case 'A': {
      miraparams=optarg;
      break;
    }
    case 'b': {
      CP_blinddata=true;
      break;
    }
    case 'c': {
      string egfc=optarg;
      if(egfc.size()!=1){
	usage();
	cout << endl;
	cerr << "ERROR: -c must be a single character\n";
	exit(1);
      }
      endgap_fillchar=egfc[0];
      break;
    }
    case 'C': {
      CP_hardtrim=true;
      break;
    }
    case 'd': {
      CP_deletestaronlycolumns=true;
      break;
    }
    case 'f': {
      CP_fromtype=optarg;
      break;
    }
    case 'F': {
      CP_filter2readgroup=true;
      break;
    }
    case 'i': {
      CP_keepnamesfromfile=false;
      break;
    }
    case 'l': {
      linelen=atoi(optarg);
      if(linelen <= 0) {
	usage();
	cout << endl;
	cerr << "ERROR: -l must be >=0\n";
	exit(1);
      }
      break;
    }
    case 'm': {
      CP_makecontigs=true;
      CP_extractreadsinsteadcontigs=false;
      break;
    }
    case 'M': {
      CP_extractreadsinsteadcontigs=true;
      CP_makecontigs=false;
      break;
    }
    case 'n': {
      CP_namefile=optarg;
      break;
    }
    case 'N': {
      CP_namefile=optarg;
      CP_sortbyname=true;
      break;
    }
    case 'o': {
      fqqualoffset=optarg;
      break;
    }
    case 'q': {
      CP_minqual=atoi(optarg);
      if(CP_minqual >100) {
	usage();
	cout << endl;
	cerr << "ERROR: -q must be <= 100\n";
	exit(1);
      }
      break;
    }
    case 'Q': {
      auto tmpv=atoi(optarg);
      if(tmpv<0 || tmpv > 100) {
	usage();
	cout << endl;
	cerr << "ERROR: -Q must be  0 <= value <= 100\n";
	exit(1);
      }
      CP_defaultqual=tmpv;
      CP_needsquality=false;
      break;
    }
    case 'r': {
      string rrr=optarg;
      for(size_t si=0; si<rrr.size(); si++){
	switch(rrr[si]){
	case 'c' :
	case 'C' :
	case 'q' : {
	  CP_recalcconopt=rrr[si];
	  break;
	}
	case 'f' :
	case 'r' : {
	  CP_recalcfeatureopt=rrr[si];
	  break;
	}
	default : {
	  cerr << "ERROR: -r must be one of c, C, q, f, r\n";
	  usage();
	  exit(1);
	}
	}
      }
      break;
    }
    case 'R': {
      CP_renamesequences=optarg;
      break;
    }
    case 's': {
      CP_splitcontigs2singlefiles=true;
      break;
    }
    case 'S': {
      CP_renamenamescheme=optarg;
      break;
    }
//    case 'v': {
//      CP_minbasecoverage=atoi(optarg);
//      break;
//    }
    case 't': {
      CP_totype.push_back(optarg);
      break;
    }
    case 'u': {
      CP_fillholesinstraingenomes=true;
      break;
    }
    case 'x': {
      CP_mincontiglength=atoi(optarg);
      break;
    }
    case 'X': {
      CP_mincontiglength=atoi(optarg);
      CP_minlengthisclipped=true;
      break;
    }
    case 'y': {
      CP_mincontigcoverage=atoi(optarg);
      break;
    }
    case 'z': {
      CP_minnumreads=atoi(optarg);
      break;
    }
    case 'Z': {
      CP_specialtestcode=true;
      break;
    }
    default : {
      cout << "Oooops? Known but unhandled option " << c << " ?\n";
      exit(100);
    }
    }
  }

  if(argc-optind < 1) {
    usage();
    cout << endl;
    cerr << argv[0] << ": " << "Missing infile and out-basename as arguments!\n";
    exit(1);
  }

  if(argc-optind < 2) {
    usage();
    cout << endl;
    cerr << argv[0] << ": " << "Missing either infile or out-basename as arguments!\n";
    exit(1);
  }

  CP_infile=argv[optind++];
  CP_outbasename=argv[optind++];

  if(CP_infile=="--help"){
    usage();
    exit(0);
  }

  guessFromAndToType(CP_infile, CP_fromtype, nullptr,
		     CP_outbasename, CP_totype, &CP_outbasename);

  // anything additional on the command line is treated as an additional totype
  for(; optind < argc; ++optind){
    CP_totype.push_back(argv[optind]);
  }

  if(CP_fromtype=="gbk" || CP_fromtype=="gbff"){
    CP_fromtype="gbf";
  }

  // uniquify totypes
  CP_totype.sort();
  CP_totype.erase(unique(CP_totype.begin(),CP_totype.end()),CP_totype.end());

  checkTypes(CP_fromtype,CP_totype);

  MIRAParameters::setupStdMIRAParameters(CP_Pv);
  if(!miraparams.empty()){
    cout << "Parsing special MIRA parameters: " << miraparams << endl;
    // switch off the checking of technology presence in the MIRAParameters object,
    //  this makes almost no sense. See bug report for -A "SOLEXA_SETTINGS -CO:fnicpst=yes"
    CP_Pv[0].getNonConstSpecialParams().sp_parse_checktechnologypresence=false;
    MIRAParameters::parse(miraparams.c_str(),CP_Pv,false);
    cout << "Ok.\n";
  }

  CP_Pv[0].getNonConstContigParams().con_output_text_cpl=linelen;
  CP_Pv[0].getNonConstContigParams().con_output_html_cpl=linelen;
  CP_Pv[0].getNonConstContigParams().con_output_text_gapfill=endgap_fillchar;
  CP_Pv[0].getNonConstContigParams().con_output_html_gapfill=endgap_fillchar;
  CP_Pv[0].getNonConstNagAndWarnParams().nw_check_templateproblems=NWNONE;

  ReadPool thepool(&CP_Pv);
  thepool.setMissingFASTAQualFileResolveMsg("use -Q");

  CP_assemblyinfo.setLargeContigSize(CP_mincontiglength);
  CP_assemblyinfo.setLargeTotalCov(CP_mincontigcoverage);

  if(CP_recalcconopt=='C'){
    for(uint32 i=0; i< CP_Pv.size(); i++){
      CP_Pv[i].setContigForceNonIUPAC(true,true);
    }
  }

  if(!CP_namefile.empty()){
    General::makeSelectionStringSet(CP_namefile);
  }

  cout << "Loading from " << CP_fromtype << ", saving to:";
  ofstream * ofstmp;
  for(list<string>::iterator ttI= CP_totype.begin(); ttI!=CP_totype.end(); ++ttI){
    ofstmp=new ofstream;
    CP_ofs.push_back(ofstmp);
    cout << ' ' << *ttI;
    if(*ttI=="fasta"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta"),ios::out);
      }
      CP_totype.push_back("fastaqual");
    } else if(*ttI=="fastaqual"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta.qual"),ios::out);
      }
    } else if(*ttI=="maskedfasta"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta"),ios::out);
      }
      CP_totype.push_back("maskedfastaqual");
    } else if(*ttI=="maskedfastaqual"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fasta.qual"),ios::out);
      }
    } else if(*ttI=="fastq"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".fastq"),ios::out);
      }
    } else if(*ttI=="caf" || *ttI=="scaf" ){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".caf"),ios::out);
      }
    } else if(*ttI=="maf"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".maf"),ios::out);
	Contig::dumpMAF_Head(*(CP_ofs.back()));
      }
    } else if(*ttI=="sam" || *ttI=="samnbb"){
      if(!CP_splitcontigs2singlefiles){
	openOFStream(*CP_ofs.back(),createFileNameFromBasePostfixContigAndRead(CP_outbasename,".sam"),ios::out);
      }
    } else if(*ttI=="gff3"){
      if(!CP_splitcontigs2singlefiles){
	CP_gffsave.open(createFileNameFromBasePostfixContigAndRead(CP_outbasename,"").c_str());
      }
    } else if(*ttI=="hsnp"){
    } else if(*ttI=="asnp"){
    } else if(*ttI=="fcov"){
    } else if(*ttI=="cstats"){
    } else if(*ttI=="crlist"){
    } else if(*ttI=="html"){
    } else if(*ttI=="text"){
    } else if(*ttI=="txt"){
    } else if(*ttI=="exp"){
    } else if(*ttI=="gbf"){
    } else if(*ttI=="gff3"){
    } else if(*ttI=="ace"){
    } else if(*ttI=="tcs"){
    } else if(*ttI=="wiggle" || *ttI=="wig"){
    } else if(*ttI=="gcwiggle" || *ttI=="gcwig"){
    } else if(*ttI=="null"){
    } else {
      BUGIFTHROW(true,"should never arrive here!");
    }
  }
  cout << '\n';

  auto cpofsI=CP_ofs.begin();
  for(list<string>::iterator ttI= CP_totype.begin(); ttI!=CP_totype.end(); ++ttI, ++cpofsI){
    if(*ttI=="sam" || *ttI=="samnbb"){
      if(CP_fromtype=="maf"){
	cout << "Collecting basic SAM info from MAF file" << endl;
	CP_samcollect.processMAF(CP_infile);
	CP_samcollect.createSAMHeader();
	*(*cpofsI) << CP_samcollect.SAMC_headerstring;
	ReadGroupLib::discard();
      }else{
	cout.flush();
	cerr << "\n\ncan only convert MAF to SAM for the time being, sorry\n";
	exit(1);
      }
    }
  }

  try{
    if(CP_fromtype=="caf" || CP_fromtype=="maf") {
      void (*usecallback)(list<Contig> &, ReadPool &) = cafmafload_callback;
      if(CP_sortbyname){
	usecallback=nullptr;
      }

      if(CP_fromtype=="caf") {
	CAF tcaf(&thepool, &CP_clist, &CP_Pv);
	vector<uint32> dummy;
	tcaf.load(CP_infile.c_str(),
		  ReadGroupLib::SEQTYPE_SANGER,
		  1,
		  dummy,
		  false,
		  usecallback
	  );
      }else if(CP_fromtype=="maf") {
	MAFParse mafp(&thepool, &CP_clist, &CP_Pv);
	vector<uint32> dummy;
	mafp.load(CP_infile.c_str(),
		  ReadGroupLib::SEQTYPE_SANGER,
		  1,
		  dummy,
		  false,
		  usecallback,
		  nullptr
	  );
      }
      if(!CP_clist.empty() && usecallback==nullptr){
	sortContigsByName(CP_clist);
	cafmafload_callback(CP_clist,thepool);
      }else if(thepool.size()!=0){
	if(CP_sortbyname) sortPoolByName(thepool,CP_namefile);
	cafmafload_callback(CP_clist,thepool);
      }
    }else{
      if(CP_fromtype=="fasta"
	 || CP_fromtype=="fastq"
	 || CP_fromtype=="fofnexp"
	 || CP_fromtype=="gbf"
	 || CP_fromtype=="gff3"){
	cout << "Loading data from " << CP_fromtype << " ...";

	string fn2;
	if(CP_fromtype=="fasta"){
	  fn2=CP_infile+".qual";
	}else if(CP_fromtype=="fastq"){
	  fn2=fqqualoffset;
	}

	ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
	rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);
	rgid.setDefaultQual(CP_defaultqual);

	string loadtype(CP_fromtype);
	if(loadtype=="fasta" && !CP_needsquality){
	  loadtype="fastanoqual";
	}

	bool canusecallback=true;
	if(CP_sortbyname){
	  canusecallback=false;
	}else{
	  int32 dummy=-1;
	  if(!fqqualoffset.empty()){
	    dummy=atoi(fqqualoffset.c_str());
	  }
	  if(0==dummy){
	    cout << "\nFor guessing the FASTQ offset, more RAM will be used as the full file needs to load into memory.";
	    canusecallback=false;
	  }
	}

	if(canusecallback){
	  thepool.loadData_rgid(loadtype, CP_infile, fn2, rgid, false, readpoolload_callback);
	}else{
	  thepool.loadData_rgid(loadtype, CP_infile, fn2, rgid, false, nullptr);
	  if(CP_sortbyname){
	    sortPoolByName(thepool,CP_namefile);
	  }
	  readpoolload_callback(thepool);
	}
      } else {
	cerr << "\n\n-f " << CP_fromtype << " is not a valid from type! (simple pool)\n";
	//usage();
	exit(1);
      }
      cout << " done.\n";
    }
  }
  catch(Notify n){
    // Need to close by hand as handleError() will perform a hard exit
    closeOpenStreams(CP_ofs);
    n.handleError(THISFUNC);
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }
  catch(...){
    cout.flush();
    cerr.flush();
    cerr << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    abort();
  }

  cout << "\nData conversion process finished, no obvious errors encountered.\n";

  FUNCEND();
  return 0;
}
