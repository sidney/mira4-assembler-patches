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

#include <iostream>
#include <string>

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

#include "util/fileanddisk.H"
#include "mira/assembly.H"
#include "mira/parameters.H"
#include "mira/manifest.H"

#include "modules/misc.H"


#include "version.H"

using namespace std;




bool   opt_mira_wrapped=false;          // 'secret' option
bool   opt_mira_resumeassembly=false;
string opt_mira_cwd;

uint32 opt_mira_numthreads=0;
bool   opt_mira_checkmanifestonly=false;
bool   opt_mira_checkmanifestandfilesonly=false;


void miraParseCmdLine(int argc, char ** argv)
{
  // that loop is straight from the GNU getopt_long example
  // http://www.gnu.org/s/hello/manual/libc/Getopt-Long-Option-Example.html
  while (1){
    static struct option long_options[] =
      {
	{"help",  no_argument,           0, 'h'},
	{"resume", no_argument,          0, 'r'},
	{"version", no_argument,         0, 'v'},
	{"mcheck", no_argument,          0, 'm'},
	{"mdcheck", no_argument,         0, 'M'},
	{"cwd", required_argument,       0, 'c'},
	{"threads", required_argument,       0, 't'},

	{"nowrap", no_argument,       0, 'W'},
	{"job", optional_argument,       0, ' '},        // catch old command line
	{"project", optional_argument,   0, ' '},        // catch old command line
	{0, 0, 0, 0}
      };
    /* getopt_long stores the option index here. */
    int option_index = 0;

    int c = getopt_long (argc, argv, "WhmMrvc:t: ",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 'W': {
      opt_mira_wrapped=true;
      break;
    }
    case 'M':
      opt_mira_checkmanifestandfilesonly=true;
      // fall through
    case 'm':
      opt_mira_wrapped=true;            // we need this: checks of manifest should not be run wrapped
      opt_mira_checkmanifestonly=true;
      break;
    case 'c':
      if(optarg){
	opt_mira_cwd=optarg;
      }else{
	cout << "Missing directory name for option -c / --cwd=" << endl;
	exit(100);
      }
      break;
    case 'h':
      cout << "mira\t\tMIRALIB version " << MIRAVERSION << "\n"
	"Author:\t\tBastien Chevreux (bach@chevreux.org)\n"
	"Purpose:\tassemble sequencing data.\n\n";

      dumpStdMsg();
      cout << "Usage:\n"
	"mira [options] manifest_file [manifest_file ...]\n";
      cout << "\nOptions:\n";
      cout <<
	"  -c / --cwd=\t\tdirectory\tChange working directory\n"
	"  -h / --help\t\t\t\tPrint short help and exit\n"
	"  -m / --mcheck\t\t\t\tOnly check the manifest file, then exit.\n"
	"  -M / --mdcheck\t\t\tLike -m, but also check existence of\n"
	"\t\t\t\t\t data files.\n"
	"  -r / --resume\t\t\t\tResume/restart an interrupted assembly\n"
	"  -t / --threads=\tinteger\t\tForce number of threads (overrides\n"
	"\t\t\t\t\t equivalent -GE:not manifest entry)\n"
	"  -v / --version\t\t\tPrint version and exit\n"
	;
      exit(0);
    case 'r':
      opt_mira_resumeassembly=true;
      break;
    case 't':
      opt_mira_numthreads=atoi(optarg);
      break;
    case 'v':
      cout << MIRAVERSION << endl;
      exit(0);
    case ' ':
      cout << "It looks like you are using the old command line format of MIRA 3.4.x and earlier."
	"\nPlease look up in the manual on how to use manifest files for MIRA 3.9.x and later.\n";
      exit(0);
    default:
      abort ();
    }
  }

  if (optind == argc) {
    cout << "You did not specify a manifest file to load?\n";
    exit(100);
  }
}

void miraMain_wrapped(int argc, char ** argv)
{
  FUNCSTART("void mira(int argc, char ** argv)");

  cout << "This is MIRA " MIRAVERSION ".\n\n";

  cout << "Please cite: Chevreux, B., Wetter, T. and Suhai, S. (1999), Genome Sequence\nAssembly Using Trace Signals and Additional Sequence Information.\nComputer Science and Biology: Proceedings of the German Conference on\nBioinformatics (GCB) 99, pp. 45-56.\n\n";

  if(!opt_mira_cwd.empty() && chdir(opt_mira_cwd.c_str())){
    cout << "Changing working directory to '" << opt_mira_cwd << "' failed, system message is: " << strerror(errno) << endl;
    exit(100);
  }

  dumpStdMsg();

  Manifest manifest;
  {
    bool tmprflag=opt_mira_resumeassembly;
    if(opt_mira_checkmanifestonly) tmprflag|=!opt_mira_checkmanifestandfilesonly;
    for(; optind < argc; ++optind) {
      manifest.loadManifestFile(argv[optind],tmprflag);
    }
  }

  // if we are just to check the manifest (and maybe data) and we came this far: exit
  if(opt_mira_checkmanifestonly) {
    cout << "\nManifest looks OK";
    if(opt_mira_checkmanifestandfilesonly){
      cout << " and the data files referenced within were found";
    }
    cout << ".\n";
    exit(0);
  }

  cout << manifest;

  vector<MIRAParameters> Pv;
  MIRAParameters::setupStdMIRAParameters(Pv);

  MIRAParameters::generateProjectNames(Pv,manifest.getProjectName());

  string mparams(manifest.getFullMIRAParameterString());
  //cout << "Seen parameters in manifest: " << mparams << endl;
  MIRAParameters::parse(mparams, Pv);

  // some users make the error to use "mira" instead of miraSearchESTSNPs
  // this code takes care of it:
  // if start_step >0, then the miraSearchESTSNPs pipeline is used,
  //  else the normal mira
  if(Pv[0].getSpecialParams().sp_est_startstep){
    MIRANOTIFY(Notify::FATAL,"Oooooops? You called the 'mira' executable but have parameters set for the"
	       "EST-SNP-Search pipeline set. For this, you have to use the 'miraSearchESTSNPs'"
	       "executable (sorry).");
  }

  cout << "\nParameters parsed without error, perfect.\n";

  if(opt_mira_numthreads>0){
    cout << "Overriding number of threads via '-t' with " << opt_mira_numthreads << endl;
    Pv[0].getNonConstAssemblyParams().as_numthreads=opt_mira_numthreads;
    Pv[0].getNonConstSkimParams().sk_numthreads=opt_mira_numthreads;
  }

  cout << '\n';

  MIRAParameters::postParsingChanges(Pv);
  MIRAParameters::dumpAllParams(Pv, cout);



  //exit(0);

  {
    Assembly as(manifest, Pv, opt_mira_resumeassembly);

    as.loadSequenceData();

    //doAbort();

    if(Pv[0].getAssemblyParams().as_filecheck_only==false){
      try {
	as.assemble();
	cout << "\n\n";
	as.setEverythingWentFine(true);
      }
      catch(const std::bad_alloc & e){
	cout << "Ouch, out of memory detected.\n";
	as.dumpMemInfo();
	throw;
      }
      catch(...){
	throw;
      }
    }
  }
  //Read::dumpStringContainerStats(cout);

  cout << "\n\nEnd of assembly process, thank you for using MIRA." << endl;

  return;
}

void miraMain_launchWrapped(int argc, char ** argv)
{
  // find out where the called binary is installed, then recreate a full command line
  string calledbin;
  findLocationOfSelfBinary(calledbin);

  string newcmdline(calledbin);
  newcmdline+=" -W";
  for(int32 i=1; i<argc; i++){
    newcmdline+=' ';
    newcmdline+=argv[i];
  }
  //cout << "Launching: " << newcmdline << endl;

  int ret=0;
  ret=system(newcmdline.c_str());
  if(ret<0){
    cout << "Ooooops? Somehow failed to launch: " << newcmdline << endl;
    cout << "Using rescue assembly mode, but extraction of large contigs will not be done automatically.\n" << endl;
    miraMain_wrapped(argc,argv);
  }else if(ret>0){
    // do nothing, just exit as the wrapped call already printed out all there is to say in case of an error
    cout << "Failure, wrapped MIRA process aborted.\n";
    exit(100);
  }else{
    bool couldextractlargecontigs=false;

    Manifest manifest;

    for(; optind < argc; ++optind) {
      manifest.loadManifestFile(argv[optind],opt_mira_resumeassembly);
    }
    vector<MIRAParameters> Pv;
    MIRAParameters::setupStdMIRAParameters(Pv);
    MIRAParameters::generateProjectNames(Pv,manifest.getProjectName());

    MIRAParameters::parse(manifest.getFullMIRAParameterString(), Pv, false);

    auto & as_fixparams=Pv[0].getAssemblyParams();
    auto & di_fixparams=Pv[0].getDirectoryParams();

    string mcprog(boost::filesystem::path(calledbin).parent_path().native()+"/miraconvert");
    if(system((mcprog+" -v >/dev/null 2>&1").c_str())){
      // try the 'd' version (just my personal setup on some machines)
      mcprog+='d';
      if(system((mcprog+" -v >/dev/null 2>&1").c_str())){
	mcprog="miraconvert";
	if(system((mcprog+" -v >/dev/null 2>&1").c_str())){
	  mcprog.clear();
	}
      }
    }

    string lcfile(di_fixparams.dir_info+"/"+as_fixparams.as_outfile_stats_largecontigs+".txt");
    if(fileExists(lcfile)){
      //cout << "found large contigs file " << lcfile << ", extracting now ..." << endl;
      string convertstr;
      if(as_fixparams.as_output_html) convertstr+=" -t html";
      if(as_fixparams.as_output_txt) convertstr+=" -t txt";
      if(as_fixparams.as_output_caf) convertstr+=" -t caf";
      if(as_fixparams.as_output_maf) convertstr+=" -t maf";
      if(as_fixparams.as_output_tcs) convertstr+=" -t tcs";
      if(as_fixparams.as_output_ace) convertstr+=" -t ace";
      if(as_fixparams.as_output_gff3) convertstr+=" -t gff3";
      if(as_fixparams.as_output_wiggle) convertstr+=" -t wig";
      if(as_fixparams.as_output_fasta) convertstr+=" -t fasta";
      if(as_fixparams.as_output_gap4da) convertstr+=" -t exp";

      string relativelcfile("../"+manifest.getProjectName()+"_d_info/"+as_fixparams.as_outfile_stats_largecontigs+".txt");
      convertstr+=" -n "+relativelcfile;

      convertstr+=" -A \"--job="+manifest.getJob()+","+manifest.getTechnologies()+" ";
      // just to be really sure ...
      for(auto pc : manifest.getParameters()){
	if(pc=='"'){
	  convertstr+="\\\"";
	}else if(pc<32 || pc==0x7f){
	  // do nothing
	}else{
	  convertstr+=pc;
	}
      }
      convertstr+="\" ";

      convertstr+=as_fixparams.as_outfile_MAF+".maf";
      convertstr+=" "+manifest.getProjectName()+"_LargeContigs_out";

      {
	string elcshname=di_fixparams.dir_results+"/extractLargeContigs.sh";
	ofstream fout(elcshname);
	fout << "#!/bin/sh"
	  "\n"
	  "\necho"
	  "\necho This is how MIRA extracted 'large contigs' from the total assembly for you:"
	  "\necho  it made a list of contigs \\(see info file " << relativelcfile << "\\)"
	  "\necho  which reached a certain length \\(usually 500, see -MI:lcs=...\\) and had at least 1/3 of"
	  "\necho  the average coverage \\(per sequencing technology\\)."
	  "\necho"
	  "\necho Note that you can redefine what large contigs are for you by simply using any"
	  "\necho combination of -n, -x, -y and -z parameters of 'miraconvert' instead of only the"
	  "\necho '-n' parameter as used in this example."
	  "\necho"
	  "\necho You can follow the progress of the conversion in the file \"ec.log\""
	  "\necho"
	  "\n\nmiraconvert " << convertstr << " >ec.log 2>&1"
	  "\n"
	  "\nif [ $? -eq 0 ];then"
	  "\n   rm ec.log"
	  "\n   echo Finished, all done."
	  "\nelse"
	  "\n   tail -50 ec.log"
	  "\n   echo"
	  "\n   echo Ooops, something went wrong. Please consult the file 'ec-log' in this directory."
	  "\nfi"
	  ;
	chmod(elcshname.c_str(),S_IRWXU|S_IRGRP|S_IROTH);
      }

      if(mcprog.empty()){
	cout << "Could not find executable 'miraconvert' for extracting large contigs?\n";
      }else{
	if(chdir(di_fixparams.dir_results.c_str())){
	  cout << "Changing working directory to '" << di_fixparams.dir_results << "' failed, system message is: " << strerror(errno) << endl;
	}else{
	  cout << "Starting extraction of large contigs now (this may take a while, please be patient)." << endl;
	  if(system((mcprog+" "+convertstr+" >ec.log 2>&1").c_str())){
	    cout << "\nExtraction failed, please consult " << di_fixparams.dir_results << "/ec.log"
	      "\nfor more information on why.\n";
	  }else{
	    cout << "\nDone.\n";
	    couldextractlargecontigs=true;
	    removeFile("ec.log",false);
	  }
	}
      }
      if(!couldextractlargecontigs){
	cout << "Large contigs could not be extracted, sorry.\n"
	  "\nDON'T PANIC! In the directory\n    '" << di_fixparams.dir_results << "'"
	  "\nyou will find a script called"
	  "   'extractLargeContigs.sh'\n"
	  "which shows you how to extract large contigs from your assembly. Either"
	  "read it to understand how it's done or ... just run it :-)\n";
      }
      chdir("../..");
    }

    string rcfile(di_fixparams.dir_results+"/_tmprecreate");
    if(fileExists(rcfile)){
      string convertstr;
      convertstr+=as_fixparams.as_outfile_MAF+".maf";
      convertstr+=" "+as_fixparams.as_outfile_CAF+".caf";
      if(mcprog.empty()){
	cout << "Could not find executable 'miraconvert' for re-creating CAF results file?\n";
      }else{
	if(chdir(di_fixparams.dir_results.c_str())){
	  cout << "Changing working directory to '" << di_fixparams.dir_results << "' failed, system message is: " << strerror(errno) << endl;
	}else{
	  cout << "Starting recreating CAF results file (please be patient)." << endl;
	  if(system((mcprog+" "+convertstr+" >rc.log 2>&1").c_str())){
	    cout << "\nOooops, that failed, please consult " << di_fixparams.dir_results << "/rc.log"
	      "\nfor more information on why.\n";
	  }else{
	    cout << "\nDone.\n";
	    removeFile("rc.log",false);
	    removeFile(rcfile,false);
	  }
	}
      }
      chdir("../..");
    }

    cout << "\n\nWe're done here, thank you for using MIRA." << endl;
  }
}

int miraMain(int argc, char ** argv)
{
  FUNCSTART("void miraMain(int argc, char ** argv)");

  miraParseCmdLine(argc,argv);
  if(opt_mira_wrapped){
    miraMain_wrapped(argc,argv);
  }else{
    miraMain_launchWrapped(argc,argv);
  }

  return 0;
}
