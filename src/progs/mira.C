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


#include <sys/resource.h>     // defines rlimit for getMoreStack()

#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"
#include "util/fileanddisk.H"
#include "progs/quirks.H"

#include "modules/mod_mira.H"
#include "modules/mod_bait.H"
#include "modules/mod_convert.H"
#include "modules/mod_memestim.H"
#include "modules/mod_tagsnp.H"

#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

using namespace std;


// This strange double declaration/definition is no error
// If the declaration is not there, then for some strange
//  reason the gcc linker (at least 4.8) will not resolve
//  the extern declaration in modules/misc.C
// This occurs only for "char[]" types, not for normal types
//  like int etc.
//
// Very strange

extern const char compileinfo[];
const char compileinfo[] = {
#include "compileinfo.itxt.xxd.H"
};



void doAbort()
{
#ifndef PUBLICQUIET
  Read::dumpStringContainerStats(cout);
#endif

  cout << "\n\nFor general help, you will probably get a quicker response on the\n"
    "    MIRA talk mailing list\n"
    "than if you mailed the author directly.\n"
    "\nTo report bugs or ask for features, please use the SourceForge ticketing\nsystem at:\n"
    "\thttp://sourceforge.net/p/mira-assembler/tickets/\n"
    "This ensures that requests do not get lost.\n";
  abort();
}



// usual linux stack size of 8Mb will lead to segfaults in very long
//  alignments (>15-18k) in align.C
// therefore, get some more stack, atm 64 Mb

void getMoreStack()
{
  struct rlimit rl;
  const rlim_t wantstacksize=64L*1024L*1024L;

  auto result = getrlimit(RLIMIT_STACK, &rl);
  if (result == 0){
    //cout << "Have cur " << rl.rlim_cur << endl;
    //cout << "Have max " << rl.rlim_max << endl;
    if(rl.rlim_cur<wantstacksize){
      rl.rlim_cur=wantstacksize;
      result = setrlimit(RLIMIT_STACK, &rl);
    }
  }else{
    cout << "could not query stack size?\n";
  }
}

int main(int argc, char ** argv)
{
#ifdef MIRAMEMORC
  MemORC::setChecking(true);
#endif

  fixQuirks();
  getMoreStack();

#ifdef TIMERESTRICTED
  {
    struct stat st;
    int rc=stat(dir_params.dir_log.c_str(),&st);
    struct tm *mytm;
    mytm=localtime(&st.st_mtime);

    if(mytm->tm_year > TR_MAXYEAR
       || (mytm->tm_year == TR_MAXYEAR
	   && mytm->tm_mon >= TR_OUT_MAXMONTH)){
      cerr << "\n\nThis version of MIRA is definitively old, please get a newer version of the assembler.\n";
      exit(0);
    }
  }
#endif

  //cpu_set_t mask;
  //cout << "####" << sizeof(cpu_set_t) << endl;
  //sched_getaffinity(0,sizeof(cpu_set_t),&mask);
  //
  //cout << "Affinity: " << (uint8 *) &mask[0] << endl;


  string path;
  string miraprog;

  splitFullPathAndFileName(argv[0],path,miraprog);
  boost::to_lower(miraprog);

  try{
    if(miraprog=="mira"
      || miraprog=="mira4"){
      miraMain(argc,argv);
//    } else if(miraprog=="mirasearchestsnps"){
//      //miraEST(argc,argv);
//    } else if(miraprog=="miraclip"){
//      //miraClip(argc,argv);
//    } else if(miraprog=="mirapre"){
//      //miraPre(argc,argv);
//    } else if(miraprog=="tagsnp"){
//      tagsnp t;
//      t.mainTagSNP(argc, argv);
    } else if(miraprog=="miramem"
	      ||miraprog=="mira4mem"){
      miraMemEstimate(argc, argv);
    } else if(miraprog=="dbgreplay"){
      dbgReplayMain(argc, argv);
    }else if(miraprog=="mirabait"
	     || miraprog=="mira4bait"){
      MiraBait m;
      m.mainMiraBait(argc, argv);
    }else if(miraprog=="miraconvert"
	     || miraprog=="mira4convert"
	     || miraprog=="convert_project"
	     || miraprog=="convert_projectd"){
      if(miraprog.front()=='c'){
	cout << "convert_project is a deprecated name for miraconvert, please switch to new name.\n";
      }
      ConvPro cp;
      cp.mainConvPro(argc, argv);
    } else {
      cout << miraprog << " is a non-recognised program name of MIRA.\n"
	"The programs SHOULD be named either\n"
	"\"mira\", \"miraconvert\", \"miramem\", \"mirabait\""
	"\nAssuming 'mira'\n" << endl;

      miraMain(argc,argv);
    }
  }
  catch(Notify n){
    n.handleError("main");
  }
  catch(Flow f){
    cout << "INTERNAL ERROR: Unexpected exception: Flow()\n";
    doAbort();
  }
  catch(const std::bad_alloc & e){
    cout << "Out of memory detected, exception message is: ";
    cout << e.what() << endl;

    if(sizeof(size_t) == sizeof(int32)){
      cout << "\nYou are running a 32 bit executable. Please note that the maximum"
	"\ntheoretical memory a 32 bit programm can use (be it in Linux, Windows or"
	"\nother) is 4 GiB, in practice less: between 2.7 and 3.3 GiB. This is valid"
	"\neven if your machine has hundreds of GiB."
	"\nShould your machine have more that 4 GiB, use a 64 bit OS and a 64 bit"
	"\nversion of MIRA.";
    }

    cout << "\n\nIf you have questions on why this happened, please send the last 1000"
      "\nlines of the output log (or better: the complete file) to the author"
      "\ntogether with a short summary of your assembly project.\n\n";

    doAbort();
  }
  catch(const ios_base::failure & e){
    cout << "Failure in IO stream detected, exception message is: "
	 << e.what() << endl
	 << "\nWe perhaps ran out of disk space or hit a disk quota?\n";
    doAbort();
  }
  catch (exception& e)
  {
    cout << "A 'standard' exception occured (that's NOT normal):\n" << e.what() << "\n\nIf the cause is not immediatly obvious, please contact: bach@chevreux.org\n\n";
    doAbort();
  }
  catch(...){
    cout << "Unknown exception caught, aborting the process.\n\nPlease contact: bach@chevreux.org\n\n";
    doAbort();
  }

#ifndef PUBLICQUIET
  Read::dumpStringContainerStats(cout);
#endif

  return 0;
}
