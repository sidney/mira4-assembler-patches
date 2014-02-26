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

#include <getopt.h>

#include "modules/mod_bait.H"
#include "modules/mod_convert.H"

#include "caf/caf.H"
#include "mira/maf_parse.H"
#include "version.H"


using namespace std;


vector<MIRAParameters> MiraBait::MB_Pv;

string MiraBait::MB_fromtype;
list<string> MiraBait::MB_totype;
list<ofstream *> MiraBait::MB_ofs;

string MiraBait::MB_infile;
string MiraBait::MB_baitfile;
string MiraBait::MB_outbasename;

bool   MiraBait::MB_loadhashstat=false;
bool   MiraBait::MB_deletestaronlycolumns=false;
bool   MiraBait::MB_mustdeletetargetfiles=true;
bool   MiraBait::MB_inversehit=false;
bool   MiraBait::MB_fwdandrev=true;
uint32 MiraBait::MB_numbaithits=1;

list<Contig> MiraBait::MB_clist;   // needed for CAF conversion (and GBF)

HashStatistics MiraBait::MB_hashstatistics;;

uint64 MiraBait::MB_numreadsread=0;
uint64 MiraBait::MB_numreadswritten=0;

MiraBait::~MiraBait()
{
  ConvPro::closeOpenStreams(MB_ofs);
}

void MiraBait::usage()
{
  cout << "mirabait\t(MIRALIB version " << MIRALIBVERSION << ")\n";
  cout << "Author: Bastien Chevreux\t(bach@chevreux.org)\n\n";

  cout << "... baiting ...\n";
  cout << "Usage:\n";
  //cout << "\tconvert_project [-f <fromtype>] [-t <totype>] [-s strainfile] [-q] infile outfile\n\n";
  cout << "mirabait [-f <fromtype>] [-t <totype> [-t <totype> ...]] [-ikLor] baitfile infile <basename_for_outfile(s)>\n\n";
  cout << "Options:\n";
  cout << "\t-f <fromtype>\tload this type of project files, where fromtype is:\n"
    "\t   caf\t\t sequences from CAF\n"
    "\t   maf\t\t sequences from MAF\n"
    "\t   phd\t\t sequences from a PHD\n"
    "\t   gbf\t\t sequences from a GBF\n"
    "\t   fasta\t sequences from a FASTA\n"
    "\t   fastq\t sequences from a FASTQ\n";
  cout << "\t-t <totype>\twrite the sequences to this type (multiple mentions\n"
    "\t\t\tof -t are allowed):\n"
    "\t   fasta\t sequences to FASTA\n"
    "\t   fastq\t sequences to FASTQ\n"
    "\t   caf\t\t sequences to CAF\n"
    "\t   maf\t\t sequences to MAF\n"
    "\t   txt\t\t sequence names to text file\n";

  cout << "\n"
    "\t-L\t\tDo not compute hash statistics, load from baitfile\n"
    "\t\t\t (which the must be a mirabait hash statistics file).\n";

  cout << "\n"
    "\t-k\t\tk-mer, length of bait in bases (<32, default=31)\n"
    "\t-n\t\tMin. number of k-mer baits needed (default=1)\n"
    "\t-i\t\tInverse hit: writes only sequences that do not hit bait\n"
    "\t-r\t\tNo checking of reverse complement direction\n";

  cout << "\n"
    "\t-o\t\tfastq quality Offset (only for -f = 'fastq')\n"
    "\t\t\t Offset of quality values in FASTQ file. Default: 33\n"
    "\t\t\t A value of 0 tries to automatically recognise.\n";



//  cout << "\t-a <string>\tString with MIRA parameters to be parsed\n"
//    "\t\t\t Useful when setting parameters affecting consensus\n"
//    "\t\t\t calling like -CO:mrpg etc.\n"
//    "\t\t\t E.g.: -a \"454_SETTINGS -CO:mrpg=3\"\n";

  cout << "\nExamples:\n"
    "\t...\n"
    "\t...\n";
}


void MiraBait::checkTypes(string & fromtype,list<string> & totype)
{
  if(fromtype.empty()){
    fromtype="fastq";
  }
  if(fromtype=="gbk" || fromtype=="gbff"){
    fromtype="gbf";
  }
  if(!(fromtype=="caf"
       || fromtype=="maf"
       || fromtype=="phd"
       || fromtype=="gbf"
       || fromtype=="exp"
       || fromtype=="fasta"
       || fromtype=="fastq"
       )){
    usage();
    cout << endl;
    cerr << "Unknown or illegal file type '" << fromtype << "' defined as <fromtype>\n";
    exit(1);
  }
  if(MB_totype.empty()){
    if(fromtype=="caf"
       || fromtype=="maf"
       || fromtype=="fasta"
       || fromtype=="fastq"
      ){
      MB_totype.push_back(fromtype);
    }else{
      MB_totype.push_back("fastq");
    }
  }
  for(list<string>::iterator ttI= MB_totype.begin(); ttI!=MB_totype.end(); ++ttI){
    if(*ttI=="scaf") *ttI="caf";
    if(!(*ttI=="fasta"
	 || *ttI=="fastq"
	 || *ttI=="caf"
	 || *ttI=="maf"
	 || *ttI=="txt"
	 )){
      usage();
      cout << endl;
      cerr << "MiraBait::checkTypes(): Unknown or illegal file type '" << *ttI << "' defined as <totype>\n";
      exit(1);
    }
  }
}

// Note: clears the readpool after saving!
void MiraBait::saveReadPool(ReadPool & rp)
{
  MB_numreadsread+=rp.size();
  MB_numreadswritten+=rp.size();

  // first, bait all reads. Those who bite, discard.
  for(uint32 i=0; i<rp.size(); ++i){
    if((MB_hashstatistics.checkBaitHit(rp[i]) >= MB_numbaithits) ^ !MB_inversehit){
      rp[i].discard();
      --MB_numreadswritten;
    }
  }

  // then save the read pool
  list<string>::iterator ttI= MB_totype.begin();
  list<ofstream *>::iterator ofsI= MB_ofs.begin();
  for(; ttI!=MB_totype.end(); ++ttI, ++ofsI){
    if(*ttI=="fasta"){
      // double indirection because iterator needs one and it is a list of ofstream pointers ...
      rp.dumpAs(*(*ofsI),Read::AS_FASTA,false);
    } else if(*ttI=="fastaqual"){
      rp.dumpAs(*(*ofsI),Read::AS_FASTAQUAL,false);
    } else if(*ttI=="fastq"){
      rp.dumpAs(*(*ofsI),Read::AS_FASTQ,false);
    } else if(*ttI=="caf" || *ttI=="scaf" ){
      rp.dumpAs(*(*ofsI),Read::AS_CAF,false);
    } else if(*ttI=="maf"){
      rp.dumpAs(*(*ofsI),Read::AS_MAF,false);
    } else if(*ttI=="txt"){
      rp.dumpAs(*(*ofsI),Read::AS_READNAME,false);
    } else {
      cout.flush();
      cerr << "\n\n-t " << *ttI << " is not a valid type when the source file does not contain a full assembly!\n";
      //usage();
      exit(1);
    }
  }
}


void MiraBait::cafmafload_callback(list<Contig> & clist, ReadPool & rp)
{
  // TODO: check if needed
  Assembly::refreshContigAndReadpoolValuesAfterLoading(rp,clist);

  saveReadPool(rp);

  Read::trashReadNameContainer();
  clist.clear();
  rp.discard();
}

void MiraBait::readpoolload_callback(ReadPool & rp)
{
  // TODO: check if needed (slows loading by ~30 to 50%
//  rp.makeTemplateIDs(false);
//  rp.makeStrainIDs(false);

  saveReadPool(rp);

  Read::trashReadNameContainer();
  rp.discard();
}



int MiraBait::mainMiraBait(int argc, char ** argv)
{
  //CALLGRIND_STOP_INSTRUMENTATION;

  FUNCSTART("int mainMiraBait(int argc, char ** argv)");

  int c;
  extern char *optarg;
  extern int optind;


  string fqqualoffset="33";

  string path;
  string convertprog;
  splitFullPathAndFileName(argv[0],path,convertprog);

  string miraparams;

  uint8 basesperhash=31;

  while (true){
    static struct option long_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"version", no_argument,         0, 'v'},
	{"loadhsf", no_argument,    0, 'L'},
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "hdiLrvf:t:o:a:k:n:",
		     long_options, &option_index);
    if(c == -1) break;

    switch (c) {
    case 'a': {
      miraparams=optarg;
      break;
    }
    case 'f': {
      MB_fromtype=optarg;
      break;
    }
    case 'L': {
      MB_loadhashstat=true;
      break;
    }
    case 'n': {
      MB_numbaithits=atoi(optarg);
      break;
    }
    case 'k': {
      uint64 bla=atoi(optarg);
      if(bla>31) bla=31;
      basesperhash=bla;
      break;
    }
    case 't': {
      MB_totype.push_back(optarg);
      break;
    }
    case 'o': {
      fqqualoffset=optarg;
      break;
    }
    case 'd': {
      MB_deletestaronlycolumns=true;
      break;
    }
    case 'i': {
      MB_inversehit=true;
      break;
    }
    case 'r': {
      MB_fwdandrev=false;
      break;
    }
    case 'h':
    case '?': {
      usage();
      exit(0);
    }
    case 'v':
      cout << MIRAVERSION << endl;
      exit(0);
    default : {}
    }
  }

  if(argc-optind < 1) {
    cerr << argv[0] << ": " << "Missing baitfile, infile and out-basename as arguments!\n";
    usage();
    exit(1);
  }

  if(argc-optind < 3) {
    cerr << argv[0] << ": " << "Missing one of baitfile, infile or out-basename as argument!\n";
    usage();
    exit(1);
  }

  if(argc-optind > 3) {
    cerr << argv[0] << ": " << "Whoops, found more than baitfile, infile and out-basename as arguments left on the command line!\n";
    cerr << "Unparsed command line: ";
    for(;optind<argc;optind++) cerr <<argv[optind] << " ";
    cerr << endl;
    usage();
    exit(1);
  }

  MB_baitfile=argv[optind++];
  MB_infile=argv[optind++];
  MB_outbasename=argv[optind];

  if(MB_baitfile=="--help"){
    usage();
    exit(0);
  }

  ConvPro::guessFromAndToType(MB_infile,MB_fromtype,nullptr,MB_outbasename,MB_totype,nullptr);

  if(MB_fromtype=="gbk" || MB_fromtype=="gbff"){
    MB_fromtype="gbf";
  }

  checkTypes(MB_fromtype,MB_totype);

  MIRAParameters::setupStdMIRAParameters(MB_Pv);
  if(!miraparams.empty()){
    cout << "Parsing special MIRA parameters: " << miraparams << endl;
    MIRAParameters::parse(miraparams.c_str(),MB_Pv,false);
    cout << "Ok.\n";
  }

  uint32 baitpoolsize=0;
  {
    ReadPool baitrp(&MB_Pv);
    if(MB_loadhashstat){
      cout << "Loading from existing hashstat file."
	"\nHope you used the same k there as you are now, else results will be wrong!\n";
      MB_hashstatistics.loadHashStatistics(baitrp,
					   MB_baitfile,
					   basesperhash);
    }else{
      cout << "Loading baits ...";
      //baitrp.loadDataFromFASTA(MB_baitfile,1, dummy, false,"",false);
      ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
      rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);
      baitrp.loadData_rgid("fastanoqual",MB_baitfile,MB_baitfile,rgid,false,nullptr);

      cout << "baitrp.size(): " << baitrp.size() << endl;

      string dummyfn;
      MB_hashstatistics.prepareHashStatistics(".",baitrp,false,false,true,MB_fwdandrev,1,basesperhash,
					      MB_Pv[0].getHashStatisticsParams().hs_million_hashes_per_buffer,
					      false,
					      dummyfn);
      baitpoolsize=baitrp.size();
    }
  }

  ReadPool loadrp(&MB_Pv);

  cout << "Loading from " << MB_fromtype << ", saving to:";
  ofstream * ofstmp;
  for(list<string>::iterator ttI= MB_totype.begin(); ttI!=MB_totype.end(); ++ttI){
    cout << ' ' << *ttI;
    ofstmp=new ofstream;
    MB_ofs.push_back(ofstmp);
    if(*ttI=="fasta"){
      MB_ofs.back()->open((MB_outbasename + ".fasta").c_str(), ios::out);
    } else if(*ttI=="fastq"){
      MB_ofs.back()->open((MB_outbasename + ".fastq").c_str(), ios::out);
    } else if(*ttI=="caf" || *ttI=="scaf" ){
      MB_ofs.back()->open((MB_outbasename + ".caf").c_str(), ios::out);
    } else if(*ttI=="maf"){
      MB_ofs.back()->open((MB_outbasename + ".maf").c_str(), ios::out);
    } else if(*ttI=="txt"){
      MB_ofs.back()->open((MB_outbasename + ".txt").c_str(), ios::out);
    } else {
      cout.flush();
      cerr << "\n\n-t " << *ttI << " is not a valid type\n";
      //usage();
      exit(1);
    }
  }
  cout << '\n';

  try{
    if(MB_fromtype=="caf") {
      CAF tcaf(&loadrp, &MB_clist, &MB_Pv);
      vector<uint32> dummy;
      tcaf.load(MB_infile.c_str(),
		ReadGroupLib::SEQTYPE_SANGER,
		1,
		dummy,
		false,
		cafmafload_callback,
		nullptr
	);
    }else if(MB_fromtype=="maf") {
      MAFParse mafp(&loadrp, &MB_clist, &MB_Pv);
      vector<uint32> dummy;
      mafp.load(MB_infile.c_str(),
		ReadGroupLib::SEQTYPE_SANGER,
		1,
		dummy,
		false,
		cafmafload_callback,
		nullptr
	);
    }else{

      uint32 dummy=0;
      if(MB_fromtype=="fasta"
	 || MB_fromtype=="fastq"
	 || MB_fromtype=="gbf"
	 || MB_fromtype=="gff3"){
	cout << "Loading data from " << MB_fromtype << " ...";

	string fn2;
	if(MB_fromtype=="fasta"){
	  fn2=MB_infile+".qual";
	}else if(MB_fromtype=="fastq"){
	  fn2=fqqualoffset;
	}
	string loadtype(MB_fromtype);
	if(loadtype=="fasta"){
	  loadtype="fastanoqual";
	}

	ReadGroupLib::ReadGroupID rgid=ReadGroupLib::newReadGroup();
	rgid.setSequencingType(ReadGroupLib::SEQTYPE_TEXT);
	loadrp.loadData_rgid(loadtype, MB_infile, fn2, rgid, false, readpoolload_callback);
      } else {
	cerr << "\n\n-f " << MB_fromtype << " is not a valid from type!\n";
	//usage();
	exit(1);
      }
      cout << " done.\n";
    }
  }
  catch(Notify n){
    // Need to close by hand as handleError() will perform a hard exit
    ConvPro::closeOpenStreams(MB_ofs);
    n.handleError("main");
  }
  catch(Flow f){
    cerr << "Unexpected exception: Flow()\n";
  }
  catch(...){
    cerr << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    abort();
  }

  cout << "\nBaiting process finished.\n\n";
  if(baitpoolsize>0){
    cout << "Number of bait sequences:   " << baitpoolsize << endl;
  }
  cout << "Number of sequences baited: " << MB_numreadswritten << " / " << MB_numreadsread << " (" << fixed << setprecision(2) << 100.0f*MB_numreadswritten/MB_numreadsread << "%)\n";

  FUNCEND();
  return 0;
}

