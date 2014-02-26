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


#include "assembly.H"

// BOOST
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>


#include "util/dptools.H"
#include "util/fileanddisk.H"
#include "util/machineinfo.H"
#include "util/misc.H"
#include "util/progressindic.H"

#include "errorhandling/errorhandling.H"
#include "mira/dynamic.H"
#include "mira/ads.H"
#include "mira/structs.H"
#include "mira/contig.H"

#include "caf/caf.H"

#ifdef MIRAMEMORC
#include "memorc/memorc.H"
#endif

#if 0
#include <valgrind/memcheck.h>
#define VALGRIND_LEAKCHECK
#endif


using namespace std;

#define CEBUG(bla)

// cs1 for normal clocking ('user compatible' as is does not disturb)
//  cs2 for extensive clocking output, more for analysis of MIRA behaviour
//  during development

#ifndef PUBLICQUIET
#define CLOCK_STEPS2
#endif
#define CLOCK_STEPS2


//#define TRACKMEMUSAGE 1
#define TRACKMEMUSAGE 0




void Assembly::test()
{
  if(AS_readpool[2].getLSClipoff()>0){
    Read::setCoutType(Read::AS_TEXT);
    cout << AS_readpool[2];
    exit(10);
  }

//  uint32 s=5000000;
//  AS_readpool.reserve(s);
//  for(uint32 i=0;i<s;i++) AS_readpool.addNewEmptyRead();
//  AS_permanent_overlap_bans.resize(s);
//  AS_istroublemaker.resize(s,0);
//  loadAlignmentsFromFile(-1,"miratmp","", ".ads_pass");
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

Assembly::Assembly(Manifest & manifest, vector<MIRAParameters> & params, bool resumeassembly): AS_readpool(&params), AS_dataprocessing(&params)
{
  FUNCSTART("Assembly::Assembly(MIRAParameters * params)");

  AS_manifest=manifest;
  AS_miraparams=params;
  AS_resumeasembly=resumeassembly;

  init();

#ifdef TIMERESTRICTED
  cout << "Compiled on " << __DATE__ << ". Will run until "<< TR_OUT_MAXDAY << "." << TR_OUT_MAXMONTH << "." << TR_OUT_MAXYEAR << " (dd.mm.yyyy)\n";
#endif

  setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);

  // For resuming assemblies, do purge the results directory!
  if(resumeassembly && ensureDirectory(AS_miraparams[0].getDirectoryParams().dir_results, true)){
    MIRANOTIFY(Notify::FATAL,"Could not delete and recreate results directory?");
  }
  // purge the the remaining directories (if we're not resuming)
  ensureStandardDirectories(!resumeassembly);

  // after ensureStandardDirectories() as usually located in info directory
  AS_warnings.setOutputPath(AS_miraparams[0].getDirectoryParams().dir_info+"/"+AS_miraparams[0].getAssemblyParams().as_outfile_stats_warnings);

  setContigBuiltCallback();

  AS_assemblyinfo.setLargeContigSize(500);

  AS_systemmemory=MachineInfo::getMemTotal();

  //makeTmpDir();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::init()
{
  AS_steps.resize(ASNUMOFSTEPS);

  for(int8 i=ASNUMOFSTEPS-1; i>=0;i--){
    AS_steps[i]=0;
  }

  AS_tmptag_CRMr=Read::REA_defaulttag_CRMr;

  // initialise some other variables
  zeroVars();
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::zeroVars()
{
  AS_num_reads_valid=0;
  AS_num_reads_too_small=0;
  AS_numADSFacts_fromalignments=0;
  AS_numADSFacts_fromshreds=0;
  AS_seqtypespresent.clear();
  AS_hasbackbones=false;
  AS_needsskimfornastyrepeats=false;

  AS_hashstat_avghashfreq=0;

  AS_doneskimchimera=false;
  AS_resumeisok=false;

#ifdef TIMERESTRICTED
  AS_timesup=false;
#endif

  AS_shouldrun_nfs_check=true;

  AS_coveragetotal=0;

  AS_everythingwentfine=false;

  //TODO: rest
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

Assembly::~Assembly()
{
  FUNCSTART("Assembly::~Assembly()");

  dumpMemInfo();

  AS_warnings.dumpWarnings();

  cout << "Dynamic s allocs: " << Dynamic::DYN_alloccounts << endl;
  cout << "Dynamic m allocs: " << Dynamic::DYN_alloccountm << endl;
  cout << "Align allocs: " << Align::AL_alloccount << endl;

  discard();

  if(AS_everythingwentfine && AS_miraparams[0].getAssemblyParams().as_output_removetmpdir){
    removeDirectory(AS_miraparams[0].getDirectoryParams().dir_tmp,false,false);
  }

  // TODO: scandir on result and remove all if no results?

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO not complete

void Assembly::discard()
{
  FUNCSTART("Assembly::discard()");

  AS_readpool.discard();
  AS_contigs.clear();
  AS_bbcontigs.clear();
  //AS_ok_for_assembly.clear();

  nukeSTLContainer(AS_adsfacts);
  nukeSTLContainer(AS_confirmed_edges);

  nukeSTLContainer(AS_used_ids);
  nukeSTLContainer(AS_multicopies);
  nukeSTLContainer(AS_hasmcoverlaps);
  nukeSTLContainer(AS_maxcoveragereached);
  nukeSTLContainer(AS_steps);
  nukeSTLContainer(AS_istroublemaker);
  nukeSTLContainer(AS_allrmbsok);
  nukeSTLContainer(AS_probablermbsnotok);
  nukeSTLContainer(AS_weakrmbsnotok);

  AS_permanent_overlap_bans.nuke();

  nukeSTLContainer(AS_readhitmiss);
  nukeSTLContainer(AS_readhmcovered);
  nukeSTLContainer(AS_count_rhm);
  nukeSTLContainer(AS_clipleft);
  nukeSTLContainer(AS_clipright);

  init();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

// TODO not complete

void Assembly::dmi_dumpALine(ostream & ostr, const char * desc, size_t numelem, size_t bytes_size, size_t freecapacity, size_t lostbyalign)
{

  ostr << setw(30) << desc
       << setw(10) << numelem;

  {
    ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(bytes_size), ostrstr);
    ostr << setw(12) << ostrstr.str();
  }
  {
    ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(freecapacity), ostrstr);
    ostr << setw(12) << ostrstr.str();
  }
  {
    ostringstream ostrstr;
    byteToHumanReadableSize(static_cast<double>(lostbyalign), ostrstr);
    ostr << setw(12) << ostrstr.str();
  }

  ostr << endl;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::dumpMemInfo()
{
  FUNCSTART("Assembly::dumpMemInfo()");

  cout << "\n\n========================== Memory self assessment ==============================\n";

  size_t bytes_size = 0;
  size_t tmp_bytes_size=0;

  size_t numelem=0;
  size_t tmp_numelem=0;

  size_t freecapacity=0;
  size_t tmp_freecapacity=0;

  size_t lostbyalign=0;
  size_t tmp_lostbyalign=0;

  // we currently do not use these
  (void) numelem;
  (void) freecapacity;
  (void) lostbyalign;


  if(sizeof(size_t) == sizeof(int32)){
    cout << "Running in 32 bit mode.\n\n";
  }else{
    cout << "Running in 64 bit mode.\n\n";
  }

  dumpFile("/proc/meminfo",cout);
  cout << '\n';
  dumpFile("/proc/self/status",cout);

  cout << "\nInformation on current assembly object:\n\n";

  cout << "AS_readpool: " << AS_readpool.size() << " reads.\n";
  cout << "AS_contigs: " << AS_contigs.size() << " contigs.\n";
  cout << "AS_bbcontigs: " << AS_bbcontigs.size() << " contigs.\n";

  bytes_size+=AS_readpool.estimateMemoryUsage();
  cout << "Mem used for reads: " << bytes_size  << " (";
  byteToHumanReadableSize(static_cast<double>(bytes_size), cout);
  cout << ")\n\nMemory used in assembly structures:\n"
       << setw(52) << "Eff. Size" << setw(12) << "Free cap." << setw(12) << "LostByAlign" << endl;


  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_writtenskimhitsperid,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_writtenskimhitsperid: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skim_edges,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skim_edges: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_adsfacts,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_adsfacts: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_confirmed_edges,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_confirmed_edges: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_permanent_overlap_bans,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_permanent_overlap_bans: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readhitmiss,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readhitmiss: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readhmcovered,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readhmcovered: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_count_rhm,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_count_rhm: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_clipleft,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_clipleft: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_clipright,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_clipright: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_used_ids,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_used_ids: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_multicopies,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_multicopies: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_hasmcoverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_hasmcoverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_maxcoveragereached,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_maxcoveragereached: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_coverageperseqtype,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_coverageperseqtype: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_istroublemaker,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_istroublemaker: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_isdebris,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_isdebris: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_needalloverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_needalloverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readsforrepeatresolve,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readsforrepeatresolve: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_allrmbsok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_allrmbsok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_probablermbsnotok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_probablermbsnotok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_weakrmbsnotok,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_weakrmbsnotok: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_readmaytakeskim,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_readmaytakeskim: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimstaken,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimstaken: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numskimoverlaps,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_numskimoverlaps: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numleftextendskims,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_numleftextendskims: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_numrightextendskims,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_rightextendskims: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimleftextendratio,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimleftextendratio: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_skimrightextendratio,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_skimrightextendratio: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;

  tmp_bytes_size=estimateMemoryUsageOfContainer(AS_usedtmpfiles,true,tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  dmi_dumpALine(cout,"AS_usedtmpfiles: ",tmp_numelem,tmp_bytes_size,tmp_freecapacity,tmp_lostbyalign);
  bytes_size+=tmp_bytes_size;


  cout << "Total: " << bytes_size << " (";
  byteToHumanReadableSize(static_cast<double>(bytes_size), cout);
  cout << ")";
  cout << "\n\n================================================================================\n";

  FUNCEND();
}

//void Assembly::dmiAddBytesWithCapacity(size_t & bcap, size_t & bsize)
//{
//}

/*************************************************************************
 *
 * go through list of previous files and delete old with same base name
 *  but different file name
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }

uint32 Assembly::cleanupOldFile(const string & basename, const string & filename)
{
  FUNCSTART("uint32 Assembly::cleanupOldFile(const string & basename, const string & filename)");

  CEBUG("\nCOF: ###" << basename << "### and ###"<<filename<<"###\n");

  boost::system::error_code ec;

  uint32 numdeleted=0;
  list<usedtmpfiles_t>::iterator ulfI=AS_usedtmpfiles.begin();
  while(ulfI != AS_usedtmpfiles.end()){
    //cout << "cOF: " << ulfI->basename << '\t' << ulfI->filename << endl;
    if(ulfI->basename == basename
       && ulfI->filename != filename) {
      numdeleted++;

      // the following two checks seem pointless as normally the list would
      //  arrive at an end sometime
      // however, I've seen one case where this loop became endless,
      //  hence this foolguard
      if(numdeleted>100) {
	cerr << "\n\nOUCH! something strange ... tried more than 100 deletes of " << basename << " ... list size is " << AS_usedtmpfiles.size() << '\n';
      }
      if(numdeleted>1200) {
	cerr << "\n\nOUCH! something weird ... tried more than 1200 deletes of " << basename << " ... list size is " << AS_usedtmpfiles.size() << '\n';
	cerr << "We'll stop that here.\n";
	return numdeleted;
      }

      removeFile(ulfI->filename,true);

      if(boost::filesystem::exists(ulfI->filename,ec)){
	cerr << "WARNING: Could not delete old file " + ulfI->filename
	     << "\nThis can have a number of different reasons, none of them"
	     << "\nwarranting an abort, but this is strange anyway.\n\n";
	ulfI++;
      }else{
	ulfI=AS_usedtmpfiles.erase(ulfI);
      }
    }else{
      ulfI++;
    }
  }

  FUNCEND();
  return numdeleted;
}
//#define CEBUG(bla)


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

string Assembly::buildFileName(int32 version, const string & prefix, const string & postfix, const string & basename, const string & suffix, const string & dirname, bool removeold)
{
  FUNCSTART("string Assembly::buildFileName(int32 version, const string & prefix, const string & postfix, const string & basename, const string & suffix, const string & dirname, bool removeold)");

  ostringstream ostr;

  if(version>=0){
    ostr << AS_miraparams[0].getDirectoryParams().dir_tmp << "/";
  } else if(!dirname.empty()){
    ostr << dirname << "/";
  }

  ostr << prefix << basename << postfix ;

  if(version>=0){
    ostr << "." << version;
  }
  ostr << suffix;

  string filename=ostr.str();
  string newbasename(basename+suffix);

  if(removeold && AS_miraparams[0].getAssemblyParams().as_output_removerollovertmps) {
    cleanupOldFile(newbasename,filename);
  }

  bool mustadd=true;
  list<usedtmpfiles_t>::const_iterator ulfI=AS_usedtmpfiles.begin();
  for(;ulfI != AS_usedtmpfiles.end(); ulfI++){
    if(ulfI->basename == newbasename
       && ulfI->filename == filename) {
      mustadd=false;
      break;
    }
  }

  if(mustadd){
    usedtmpfiles_t ulf;

    AS_usedtmpfiles.push_back(ulf);
    AS_usedtmpfiles.back().basename=newbasename;
    AS_usedtmpfiles.back().filename=filename;
  }

  //{
  //  cout << "IVEGOT: " << AS_usedtmpfiles.size() << endl;
  //  //cout << "IVEGOT: " << filename << "IVEGOT" << endl;
  //  list<usedtmpfiles_t>::const_iterator ulfI=AS_usedtmpfiles.begin();
  //  for(;ulfI != AS_usedtmpfiles.end(); ulfI++){
  //    cout << "ulf: " << ulfI->basename << '\t' << ulfI->filename << endl;
  //  }
  //}

  FUNCEND();
  return filename;
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::ensureStandardDirectories(bool purge){
  FUNCSTART("void Assembly::ensureStandardDirectories(bool purge)");

  auto & dparams = AS_miraparams[0].getNonConstDirectoryParams();

  string existingtmpdir;
  // special purge logic to handle eventually existing symlinks to -DI:trt
  if(!dparams.dir_tmp_symlink.empty()
     && boost::filesystem::exists(dparams.dir_tmp_symlink)
     && boost::filesystem::is_symlink(dparams.dir_tmp_symlink)){
    existingtmpdir=boost::filesystem::read_symlink(dparams.dir_tmp_symlink).native();
  }else if(boost::filesystem::exists(dparams.dir_tmp)
	   && boost::filesystem::is_symlink(dparams.dir_tmp)){
    existingtmpdir=boost::filesystem::read_symlink(dparams.dir_tmp).native();
  }
  if(purge
     && !existingtmpdir.empty()){
    boost::filesystem::remove_all(existingtmpdir);
    existingtmpdir.clear();
  }

  // make sure the main directories exist or are created
  if(ensureDirectory(dparams.dir_top, purge, true, false)
     || ensureDirectory(dparams.dir_results, purge, true, false)
     || ensureDirectory(dparams.dir_info, purge, true, false)
     || ensureDirectory(dparams.dir_checkpoint, purge, true, false)){

    MIRANOTIFY(Notify::FATAL, "Could not make sure that a needed directory exists (see log above for more info), aborting MIRA.");
  }

  // make sure the tmp directory exists or is created
  // either as directory or as symlinked directory (-DI:trt=...)
  if(dparams.dir_tmp_symlink.empty()){
    if(ensureDirectory(dparams.dir_tmp, purge, true, false)){
      MIRANOTIFY(Notify::FATAL, "Could not make sure that the MIRA tmp directory exists, aborting.");
    }
  }else{
    if(!existingtmpdir.empty() && !boost::filesystem::is_directory(existingtmpdir)){
      // existing ist not a directory??? Should be very, very ... very rare
      // but if yes, let's purge it
      boost::filesystem::remove_all(existingtmpdir);
      existingtmpdir.clear();
    }
    if(existingtmpdir.empty()){
      existingtmpdir=dparams.dir_tmp+"_XXXXXX";
      auto * ptr=mkdtemp(const_cast<char *>(existingtmpdir.c_str()));
      if(ptr==nullptr){
	perror(static_cast<string>("Could not create directory for temporary MIRA data "+existingtmpdir).c_str());
      }
      // contrary to what "man mkdtemp" may lead to believe, mkdtemp() does not return NULL (failure)
      //  if given a path to a non-existing directory ("/bla/bla_XXXXXX")
      // need to check ourselves
      if(!boost::filesystem::exists(existingtmpdir)){
	MIRANOTIFY(Notify::FATAL, "Could not create MIRA tmp directory \"" << existingtmpdir << "\": is some part of the path not existing or access protected?");
      }
    }

    dparams.dir_tmp=existingtmpdir;

    if(!boost::filesystem::exists(dparams.dir_tmp_symlink)){
      boost::filesystem::create_symlink(dparams.dir_tmp,dparams.dir_tmp_symlink);
    }

    cout << "Symlink " << dparams.dir_tmp_symlink << " now pointing to " << dparams.dir_tmp << endl;;
  }

  checkForNFSMountOnTmpDir();

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::checkForNFSMountOnTmpDir()
{
  FUNCSTART("void Assembly::checkForNFSMountOnTmpDir()");

  if(!AS_shouldrun_nfs_check) return;

  auto res=checkForNFSMountOnDirectory(AS_miraparams[0].getDirectoryParams().dir_tmp, true);

  cout << '\n';
  if(res==0){
    cout << "Tmp directory is not on a NFS mount, good.\n\n";
  }else if(res==1){
    cout << "\nMake sure " << AS_miraparams[0].getDirectoryParams().dir_tmp
	 << " is *NOT* on a NFS mount or else MIRA will run *very* slowly.\n";
  }else{
    cout << "\n\n\n\n\nWARNING WARNING WARNING!\n\n"
      "It looks like the directory MIRA uses for temporary files\n    " << AS_miraparams[0].getDirectoryParams().dir_tmp <<
      "\nis on a NFS (Network File System) mount. This will slow down MIRA *considerably*\n"
      "... by about a factor of 10!\n\n"
      "If you don't want that, you have three possibilities:\n\n"
      "1) RECOMMENDED! Use -DI:trt to redirect the tmp directory somewhere else on a\n"
      "   local disk or even SSD.\n"
      "2) ALSO POSSIBLE: put the whole project somewhere else on your file system.\n"
      "3) ABSOLUTELY NOT RECOMMENDED AT ALL: use \"-NW:cnfs=warn\" to tell MIRA not\n"
      "   to stop when it finds the tmp directory on NFS.\n\n"
      "If you do not know what NFS is and which directory to use in \"-DI:trt\", ask\n"
      "your local system administrator to guide you.\n\n";

    if(AS_miraparams[0].getNagAndWarnParams().nw_check_nfs==NWSTOP){
      MIRANOTIFY(Notify::FATAL,"Tmp directory is on a NFS mount ... but we don't want that.");
    }
  }

  AS_shouldrun_nfs_check=false;

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::preassembleTasks(bool usereadextension, bool clipvectorleftovers)
{
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  {
    string tmpfname;
    tmpfname=buildFileName(0,"","",
			  as_fixparams.as_tmpf_clippings,
			   ".txt","",0);

    // doing it twice catches a few outliers missed the first time
    string logprefix="proposed cutback 1a: ";
    uint64 numclipped=performNewProposedCutbackClips(tmpfname,logprefix);
    if(numclipped>0){
      logprefix="proposed cutback 1b: ";
      performNewProposedCutbackClips(tmpfname,logprefix);
    }else{
      cout << "No bases clipped in first pec round, skipping second round.\n";
    }
    dumpSomeStatistics();
  }

  //performSnapshot(0);

  //performHashEditing();

#if TRACKMEMUSAGE
  cout << "\ndmi pre 00\n";
  dumpMemInfo();
#endif

/*
  // find out whether there are SRMr or CRMr tags that need attention
  bool initialrepeatmarkerspresent=false;
  for(uint32 i=0; i< AS_readpool.size() && !initialrepeatmarkerspresent; i++){
    if(AS_readpool[i].isUsedInAssembly()){
      if(AS_readpool[i].hasTag(Read::REA_tagentry_idSRMr)){
	initialrepeatmarkerspresent=true;
      }else if(AS_readpool[i].hasTag(Read::REA_tagentry_idCRMr)){
	initialrepeatmarkerspresent=true;
      }
    }
  }

  if(initialrepeatmarkerspresent){
    performHashAnalysis(false, 0,"","");

    //performSnapshot(0);

#if TRACKMEMUSAGE
    cout << "\ndmi pre 10\n";
    dumpMemInfo();
#endif

    cout << "Repeat markers found in data loaded, performing an initial step of marker propagation.\nInitial marker propagation 1/2\n";
    findPossibleOverlaps(0, "", "_initialRMpropagation");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 20\n";
    dumpMemInfo();
#endif

    AS_steps[ASADSLISTOK]=0;
    makeAlignments(Assembly::ma_needSRMrOrTwoCRMr, true, false, 0, "", "", "initialRMpropagation1");
#if TRACKMEMUSAGE
    cout << "\ndmi pre 30\n";
    dumpMemInfo();
#endif

    AS_steps[ASADSLISTOK]=0;
    cout << "Initial marker propagation 2/2\n";
    makeAlignments(Assembly::ma_needSRMrOrTwoCRMr, true, false, 0, "", "", "initialRMpropagation2");
    AS_steps[ASADSLISTOK]=0;
  }
*/

#if TRACKMEMUSAGE
  cout << "\ndmi pre 50\n";
  dumpMemInfo();
#endif

  //if(AS_454dosimpleedit) editSimple454Overcalls(0);

  if(clipvectorleftovers
	 || (usereadextension && as_fixparams.as_readextension_firstpassnum == 0)){
    cout << "Pre-assembly alignment search for read extension and / or vector clipping:\n";

    findPossibleOverlaps(0, "", "_preassembly");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 61\n";
    dumpMemInfo();
#endif

    // do not use the 100% trans rule for read extension and clipping!
    // needed vectors would not get filled
    makeAlignments(Assembly::ma_takeall, false, false, 0, "", "_preassembly1");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 62a\n";
    dumpMemInfo();
#endif

    loadAlignmentsFromFile(0, "", "_preassembly1");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 62b\n";
    dumpMemInfo();
#endif

    if(usereadextension) {
      cout << "Pre-assembly read extension:\n";
      extendADS(0, "", "_preassembly1");
      AS_needsskimfornastyrepeats=true;
#if TRACKMEMUSAGE
      cout << "\ndmi pre 62c\n";
      dumpMemInfo();
#endif
    }
    if(clipvectorleftovers) {
      cout << "Pre-assembly vector clipping\n";
      performSeqVectorClippings();
      AS_needsskimfornastyrepeats=true;
#if TRACKMEMUSAGE
      cout << "\ndmi pre 62d\n";
      dumpMemInfo();
#endif
    }

    dumpSomeStatistics();

    // we need to throw away all permbans that might have appeared in
    //  this pre-assembly sequence massage

    AS_permanent_overlap_bans.nuke();
    AS_permanent_overlap_bans.resize(AS_readpool.size());

#if TRACKMEMUSAGE
    cout << "\ndmi pre 63\n";
    dumpMemInfo();
#endif

    performHashAnalysis(false, false, 0, "", "_preassembly2");

#if TRACKMEMUSAGE
    cout << "\ndmi pre 64\n";
    dumpMemInfo();
#endif

    {
      string tmpfname;
      tmpfname=buildFileName(0,"","",
			    as_fixparams.as_tmpf_clippings,
			     ".txt","",false);
      string logprefix="proposed cutback preassembly: ";

      performNewProposedCutbackClips(tmpfname,logprefix);
      dumpSomeStatistics();
#if TRACKMEMUSAGE
      cout << "\ndmi pre 64b\n";
      dumpMemInfo();
#endif
    }

    //nukeSTLContainer(AS_adsfacts);
    //nukeSTLContainer(AS_confirmed_edges);
    AS_adsfacts.clear();
    AS_confirmed_edges.clear();

#if TRACKMEMUSAGE
    cout << "\ndmi pre 65a\n";
    dumpMemInfo();
#endif
//    findPossibleOverlaps(0, "", "_preassembly2");
//#if TRACKMEMUSAGE
//    cout << "\ndmi pre 65b\n";
//    dumpMemInfo();
//#endif
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

bool Assembly::checkTerminationRequest()
{
  struct stat st;
  string fname=AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate";
  int rc=stat(fname.c_str(),&st);
  if(rc==0) {
    string command="mv "
      +AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate"
      +" "
      +AS_miraparams[0].getDirectoryParams().dir_checkpoint+"/terminate_acknowledged";
    int dummy=system(command.c_str());
    return true;
  }
  return false;
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::assemble()
{
  FUNCSTART("Assembly::assemble()");

#ifdef MIRAMEMORC
  cout.flush();
  MemORC::statistics();
  MemORC::checkAllMemBlocks();
#endif

  basicDataChecks();

  uint32 startpass=1;
  AS_guessedtemplatevalues=false;

  ensureStandardDirectories(false);

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  //directory_parameters const & dir_params= AS_miraparams->getDirectoryParams();
  edit_parameters const & ed_params= AS_miraparams[0].getEditParams();

  AS_basesperhash_base=AS_miraparams[0].getSkimParams().sk_basesperhash;

  AS_donequickdenovocoveragecheck=false;

  // fill quick lookup which sequencing types are present as well as whether backbones are there
  AS_seqtypespresent.clear();
  AS_seqtypespresent.resize(ReadGroupLib::SEQTYPE_END,false);
  for(uint32 rgi=0; rgi< ReadGroupLib::getNumReadGroups(); ++rgi){
    auto rgid=ReadGroupLib::getReadGroupID(rgi);
    if(rgid.isBackbone()){
      AS_hasbackbones=true;
    }
    if(!rgid.isRail() && !rgid.isBackbone()){
      AS_seqtypespresent[rgid.getSequencingType()]=true;
    }
  }

  if (as_fixparams.as_numrmbbreakloops <1) {
    const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams()).as_numrmbbreakloops=1;
    cout << "Number of RMB break loops <1, setting to 1\n";
  }
  if (AS_hasbackbones
      && as_fixparams.as_startbackboneusage_inpass > static_cast<int32>(as_fixparams.as_numpasses)) {
    const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams()).as_startbackboneusage_inpass=as_fixparams.as_numpasses;
    cout << "Start of backbone usage > number of passes, correcting start to " << as_fixparams.as_numpasses << endl;
  }

  if(AS_systemmemory==0 && as_fixparams.as_automemmanagement){
    cout << "Can't find info about system or process memory, switching off automatic"
      "\nmemory management.\n";
    AS_miraparams[0].getNonConstAssemblyParams().as_automemmanagement=false;
  }

  AS_assemblyinfo.setLargeContigSize(AS_miraparams[0].getSpecialParams().mi_as_largecontigsize);
  AS_assemblyinfo.setLargeContigSizeForStats(AS_miraparams[0].getSpecialParams().mi_as_largecontigsize4stats);

  if((AS_seqtypespresent[ReadGroupLib::SEQTYPE_454GS20]
      || AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]
      || AS_seqtypespresent[ReadGroupLib::SEQTYPE_IONTORRENT])
     && (as_fixparams.as_output_gap4da
	 || as_fixparams.as_output_tmp_gap4da
	 || as_fixparams.as_output_exttmp_gap4da)){
    cout << "454, Solexa or IonTorrent data present, switching off GAP4DA type output results (you *DO NOT* want millions of files in a directory, really.)\n";
    MIRAParameters::parseQuickmode("-OUT:org=no:otg=no:oetg=no",
				   "", AS_miraparams);
  }

  // look for template ids found
  AS_maxtemplateid=-1;
  for(uint32 rpi=0; rpi<AS_readpool.size(); ++rpi){
    if(AS_readpool[rpi].getTemplateID()>AS_maxtemplateid) AS_maxtemplateid=AS_readpool[rpi].getTemplateID();
  }
  cout << "PRED MAXTID " << AS_maxtemplateid << endl;

  // allocate or reserve memory that is quite static from the size,
  //  reducing memory fragmentations at least a bit
  {
#if TRACKMEMUSAGE
    cout << "\ndmi as_init 00\n";
    dumpMemInfo();
#endif
    if(AS_maxtemplateid>=0){
      AS_templateguesses.reserve(AS_maxtemplateid+1);
    }

    AS_multicopies.reserve(AS_readpool.size());  // do NOT init multicopies!
    AS_hasmcoverlaps.reserve(AS_readpool.size()); // does not need init

    AS_permanent_overlap_bans.nuke();
    AS_permanent_overlap_bans.resize(AS_readpool.size());

    AS_used_ids.resize(AS_readpool.size(),0);
    AS_clipleft.resize(AS_readpool.size(),0);
    AS_clipright.resize(AS_readpool.size(),0);
    AS_istroublemaker.resize(AS_readpool.size(),0);
    AS_debrisreason.resize(AS_readpool.size(),0);
    AS_isdebris.resize(AS_readpool.size(),0);
    AS_needalloverlaps.resize(AS_readpool.size(),false);
    AS_maxcoveragereached.resize(AS_readpool.size(),0);

    //AS_allowquickoverlap.resize(AS_readpool.size(),false);

    // these 3 only temp filled, but block anyway to reduce fragmentation
    AS_allrmbsok.reserve(AS_readpool.size());
    AS_probablermbsnotok.reserve(AS_readpool.size());
    AS_weakrmbsnotok.reserve(AS_readpool.size());
#if TRACKMEMUSAGE
    cout << "\ndmi as_init 10\n";
    dumpMemInfo();
#endif
  }

  AS_needsskimfornastyrepeats=true;

  bool clipvectorleftovers=false;
  bool usereadextension=false;
  for(uint32 i=0; i<ReadGroupLib::SEQTYPE_END; i++){
    if(AS_seqtypespresent[i]){
      if(AS_miraparams[i].getAssemblyParams().as_clip_possible_vectors){
	clipvectorleftovers=true;
      }
      if(AS_miraparams[i].getAssemblyParams().as_use_read_extension){
	usereadextension=true;
      }
    }
  }

  AS_resumeisok=false;
  if(AS_resumeasembly){
    AS_resumeisok=true;
    try{
      loadSnapshotData(startpass);
    }
    catch(Notify n){
      cout << "Error while loading snapshot metadata, resuming assembly is not possible, sorry\n";
      n.handleError(THISFUNC);
    }
    // set a couple of values
    if(startpass>1) AS_doneskimchimera=true;
  }else{
    preassembleTasks(usereadextension,clipvectorleftovers);
    performSnapshot(1);
  }

  EDITParameters eparams;
#ifdef MIRA_HAS_EDIT
  //  eparams.setDoEval();
  eparams.setVerbose(0);
  eparams.setShowProgress(true);
#else
#endif

  uint32 actpass=startpass;
  if(as_fixparams.as_numpasses==0){
    // 0 passes? Then the user does not want an assembly
    // Just calculate some statistics and the hash analysis

    // currently not needed I think, just to be sure
    AS_resumeisok=false;

    cout << "You have selected to have 0 passes on the assembly, therefore"
      "\njust running a couple of hash statistics and info for read repeats."
      "\n(also performing rare kmer clips if wished)\n";
    dumpSomeStatistics();
    performHashAnalysis(true, true, 0, "", "_pass");

    dumpSomeStatistics();
    performSnapshot(1);
  }else{
    bool rerunSKIM=true;
    for(; actpass<=as_fixparams.as_numpasses; actpass++){

      //if(actpass==2) {
      //	Contig::setMasterCEBUGFlag(true);
      //}else{
      //	Contig::setMasterCEBUGFlag(false);
      //}

      //dumpMemInfo();
      dumpFile("/proc/self/status",cout);

      cout << "\n\nPass: " << actpass << " / " << as_fixparams.as_numpasses << endl;

      AS_warnings.dumpWarnings();

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK1 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK1 debugging end\n";

      if(actpass==20) {
	cout << "\n==MEMTRACK exiting\n";
	dumpMemInfo();
	exit(0);
      }
#endif
#ifdef MIRAMEMORC
      cout.flush();
      MemORC::statistics();
      MemORC::checkAllMemBlocks();
#endif

      AS_contigs.clear();

#ifdef MIRA_HAS_EDIT
      if(actpass<as_fixparams.as_numpasses) {
	// eventually strict editing in first passes
	//  eparams.setStrictEvaluation(false);
	eparams.setStrictEvaluation(true);
	eparams.setConfirmationThreshold(static_cast<float>(0.8));
      } else {
	// eventually lazy editing in last passes
	eparams.setStrictEvaluation(ed_params.ed_strict_editing_mode);
	eparams.setConfirmationThreshold(static_cast<float>(ed_params.ed_confirmation_threshold/100.0));
      }
#endif

      // logic to increase the -SK:bph
      // used the AS_hashstat_avghashfreq computed in earlier skim passes
      if(actpass>=2
	 && AS_miraparams[0].getSkimParams().sk_bph_increasestep>0
	 && AS_hashstat_avghashfreq>AS_miraparams[0].getSkimParams().sk_bph_coveragethreshold){
	// we'll increase maximum of 3 per pass
	uint32 maxbph=AS_basesperhash_base+(actpass-1)*3;
	uint32 newbph=AS_basesperhash_base+(actpass-1)*static_cast<double>(AS_hashstat_avghashfreq)/static_cast<double>(AS_miraparams[0].getSkimParams().sk_bph_coveragethreshold);
	if(newbph>maxbph) newbph=maxbph;
	if(newbph>31) newbph=31;
	AS_miraparams[0].getNonConstSkimParams().sk_basesperhash=newbph;
	cout << "Automatic -SK:bph increased to " << newbph << endl;
      }


      if(as_fixparams.as_skimeachpass
	 || rerunSKIM                  // forced by new SRMs in previous pass
	 || usereadextension
	 || clipvectorleftovers
	 || AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]){

	dumpSomeStatistics();

#if TRACKMEMUSAGE
	cout << "\ndmi 20\n";
	dumpMemInfo();
#endif
	// as something might have changed in the reads, redo
	//  the hash analysis
	performHashAnalysis(actpass==1, actpass==1, actpass, "", "_pass");


#if TRACKMEMUSAGE
	cout << "\ndmi 30\n";
	dumpMemInfo();
#endif

	findPossibleOverlaps(actpass, "", "_pass");
      }

      rerunSKIM=false;

      // TODO: setting that here is messy
      AS_steps[ASVECTORSCLIPPED]=0;

#if TRACKMEMUSAGE
      cout << "\ndmi 40\n";
      dumpMemInfo();
#endif

      {
	string signalfile(buildFileName(actpass, "", "_pass",
					AS_miraparams[0].getAssemblyParams().as_tmpf_signal_mainalignments,
					".ok"));
	if(!AS_resumeasembly || !AS_resumeisok || !fileExists(signalfile)){
	  AS_resumeisok=false;
	  makeAlignments(Assembly::ma_takeall, false, true, actpass, "", "_pass");
	  saveResumeDataMA(actpass, "", "_pass");
	  ofstream fout(signalfile.c_str());  // create checkpoint signal file for main alignments
	}else{
	  cout << "Resume assembly: alignments already present, good.\n";
	  loadResumeDataMA(actpass, "", "_pass");
	}
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 50\n";
      dumpMemInfo();
#endif

      loadAlignmentsFromFile(actpass, "", "_pass");

#if TRACKMEMUSAGE
      cout << "\ndmi 60\n";
      dumpMemInfo();
#endif

      // count all SRMr tags in reads
      // idea: after building contigs, run the repeat resolver with
      //  all reads that have new tags

      vector<uint32> xrmrcount(AS_readpool.size(),0);
      for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	xrmrcount[rnr]=AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idSRMr);
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 70\n";
      dumpMemInfo();
#endif

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK2 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK2 debugging end\n";
#endif

      if(as_fixparams.as_dateoutput) dateStamp(cout);

      bool foundrepeats=buildFirstContigs(actpass,
					  eparams,
					  (actpass==as_fixparams.as_numpasses));

#ifdef VALGRIND_LEAKCHECK
      cout << "\n==MEMTRACK3 debugging start\n";
      dumpMemInfo();
      VALGRIND_DO_LEAK_CHECK
	cout << "\n==MEMTRACK3 debugging end\n";
#endif

#if TRACKMEMUSAGE
      cout << "\ndmi 80\n";
      dumpMemInfo();
#endif

      // try to guess a good value for minimum coverage
      //  for large contigs
      AS_assemblyinfo.calcCurrentInfo();
      {
	auto avc=AS_assemblyinfo.ASI_avgcoverage[1];
	if(avc==0){
	  avc=AS_assemblyinfo.ASI_avgcoverage[0];
	}
	double multval=0.5d;
	if(avc<40) multval=1.0d/3;
	AS_assemblyinfo.setLargeTotalCov(
	  static_cast<uint32>(.5+avc*multval)
	  );
	for(uint8 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	  avc=AS_assemblyinfo.ASI_avgcoverage_perst[1][st];
	  if(avc==0) avc=AS_assemblyinfo.ASI_avgcoverage_perst[0][st];
	  AS_assemblyinfo.setLargeContigCovPerST(
	    static_cast<uint32>(.5+avc*multval),
	    st
	    );
	}
      }
      saveAssemblyInfo();

      // save large contigs info only for genome de-novo,
      //  i.e., no backbones present and we're using genomic pathfinder (not EST)
      // ah, and hunt for smile shape
      if(AS_bbcontigs.empty() && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
	saveLargeContigsInfo();
	warnAtSmileCoverage();
      }

#if TRACKMEMUSAGE
      cout << "\ndmi 90\n";
      dumpMemInfo();
#endif

      // some things are worth doing only if it's not the last pass
      if(actpass!=as_fixparams.as_numpasses){

	if(foundrepeats){
	  cout << "Repeats found during contig building, adding additional alignment iteration\nfor quick repeat resolving.\n";
	  AS_steps[ASADSLISTOK]=0;
	  //makeAlignments(Assembly::ma_needSRMrOrTwoCRMr, true, actpass, "", "", "repeat_resolve");

	  // Step 1 of repeat resolving
	  // define which reads should be taken into the alignments phase
	  //  (== those who got additional SRMr tags)
	  AS_readsforrepeatresolve.clear();
	  AS_readsforrepeatresolve.resize(AS_readpool.size(),false);
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    if(AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idSRMr) != xrmrcount[rnr]) AS_readsforrepeatresolve[rnr]=true;
	  }

	  // prepare for second round of repeat resolve
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    xrmrcount[rnr]=AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idCRMr);
	  }

#if TRACKMEMUSAGE
	  cout << "\ndmi a0\n";
	  dumpMemInfo();
#endif
	  makeAlignments(Assembly::ma_needRRFlag, true, false, actpass, "", "", "repeat_resolve");
	  AS_steps[ASADSLISTOK]=0;

	  // Step 2 of repeat resolving
	  // define which reads should be taken into the alignments phase
	  //  (== those who got additional CRMr tags in alignment phase above)

	  AS_readsforrepeatresolve.clear();
	  AS_readsforrepeatresolve.resize(AS_readpool.size(),false);
	  for(uint32 rnr=0; rnr<AS_readpool.size(); rnr++){
	    if(AS_readpool.getRead(rnr).countTags(Read::REA_tagentry_idCRMr) != xrmrcount[rnr]) AS_readsforrepeatresolve[rnr]=true;
	  }

#if TRACKMEMUSAGE
	  cout << "\ndmi b0\n";
	  dumpMemInfo();
#endif
	  makeAlignments(Assembly::ma_needRRFlagAndBothCRMr, true, false, actpass, "", "", "repeat_resolve");

#if TRACKMEMUSAGE
	  cout << "\ndmi c0\n";
	  dumpMemInfo();
#endif
	  nukeSTLContainer(AS_readsforrepeatresolve);

#if TRACKMEMUSAGE
	  cout << "\ndmi d0\n";
	  dumpMemInfo();
#endif
	  rerunSKIM=true;
	}

	if(usereadextension
	   && actpass >= as_fixparams.as_readextension_firstpassnum
	   && actpass <= as_fixparams.as_readextension_lastpassnum) {
	  extendADS(actpass, "", "_pass");
	  AS_needsskimfornastyrepeats=true;
	  rerunSKIM=true;
	}
	if(clipvectorleftovers) {
	  performSeqVectorClippings();
	  AS_needsskimfornastyrepeats=true;
	  rerunSKIM=true;
	}
      }

      if(checkTerminationRequest()){
	cout << "Seen termination request by user\n";
	if(actpass + 2 < as_fixparams.as_numpasses){
	  assembly_parameters & ap=const_cast<assembly_parameters &>(AS_miraparams[0].getAssemblyParams());
	  cout << "Changing number of passes from " << ap.as_numpasses << " to ";
	  ap.as_numpasses=actpass+2;
	  cout << ap.as_numpasses << endl;
	}else{
	  cout << "No further action necessary, will terminate anyway in 2 passes max.\n";
	}
      }

#if TRACKMEMUSAGE
      cout << "\ndmi e0\n";
      dumpMemInfo();
#endif

      performSnapshot(actpass+1);

    }

    if(AS_hasbackbones && AS_guessedtemplatevalues){
      priv_hackMergeTwoResultMAFs();
    }

  }

  AS_warnings.dumpWarnings();


#ifdef MIRAMEMORC
  cout.flush();
  MemORC::statistics();
  MemORC::checkAllMemBlocks();
#endif

  FUNCEND();
  return;
}

/*************************************************************************
 *
 * In mapping assemblies, the results for guessing template numbers
 *  (size, segment placement) are available only after the MAF has been written
 * Therefore, need to rewrite header of MAF (and CAF completely)
 * Easiest way out: MAF header from chkpoint, body of 'usual' result MAF
 *
 * delete CAF, have outer caller recreate it from new MAF
 *
 *************************************************************************/

void Assembly::priv_hackMergeTwoResultMAFs()
{
  string headermaf(buildDefaultCheckpointFileName(AS_miraparams[0].getAssemblyParams().as_infile_chkptMAF));
  if(!fileExists(headermaf)) return;

  string bodymaf(getMAFFilename());
  if(!fileExists(bodymaf)) return;

  string newmaf(bodymaf+"tmp");

  std::ifstream fin(headermaf, ios::in);
  if(!fin.is_open()){
    cout << "MAFmerge: Could not open headermaf " << headermaf << endl;
    return;
  }

  std::ofstream  fout(newmaf,ios::out);
  if(!fout.is_open()){
    cout << "MAFmerge: Could not open newmaf " << newmaf << endl;
    return;
  }

  string actline;
  actline.reserve(1000);
  while(!fin.eof()){
    getline(fin,actline);
    if(!actline.empty()){
      if(actline[0]=='@'
	 || actline[0]=='#'){
	fout << actline << '\n';
      }else{
	break;
      }
    }
  }
  fin.close();
  fin.open(bodymaf, ios::in);
  if(!fin.is_open()){
    cout << "MAFmerge: Could not open bodymaf " << bodymaf << endl;
    return;
  }
  while(!fin.eof()){
    getline(fin,actline);
    if(actline.empty()) continue;
    if(actline[0]!='@'
       && actline[0]!='#') break;
  }
  fout << actline << '\n';
  fout << fin.rdbuf();

  if(fout.bad()){
    cout << "MAFmerge: Could not finish copying to " << newmaf << endl;
  }
  fin.close();
  fout.close();

  fileRename(newmaf,bodymaf);
  string caffile(getCAFFilename());
  if(fileExists(caffile)){
    removeFile(caffile,true);
    fout.open(AS_miraparams[0].getDirectoryParams().dir_results+"/_tmprecreate");
    fout << caffile << endl;
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::performSnapshot(uint32 actpass)
{
  FUNCSTART("void Assembly::performSnapshot(uint32 actpass)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  cout << "Performing snapshot " << actpass << endl;
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  try{
    boost::filesystem::rename(AS_miraparams[0].getDirectoryParams().dir_checkpoint,AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp);
  }
  catch(boost::filesystem::filesystem_error fse){
    MIRANOTIFY(Notify::FATAL,"Could not rename checkpoint directory in preparation of saving new checkpoint?\nError message is " << fse.what());
  }
  catch(...){
    MIRANOTIFY(Notify::FATAL,"Could not rename checkpoint directory in preparation of saving new checkpoint?\nUnspecified exception?");
  }
  if(ensureDirectory(AS_miraparams[0].getDirectoryParams().dir_checkpoint, false)){
    // Oooops? Cannot create new? Well, roll back the rename and exit
    try{
      boost::filesystem::rename(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp,AS_miraparams[0].getDirectoryParams().dir_checkpoint);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Could not rollback initial checkpoint rename?\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Could not rollback initial checkpoint rename?\nUnspecified exception?");
    }
    MIRANOTIFY(Notify::FATAL,"Could not create new snapshot directory? Disk full? changed permissions?");
  }
  try{
    ssdReadPool(buildDefaultCheckpointFileName(as_fixparams.as_infile_chkptMAF));
    ssdPassInfo(buildDefaultCheckpointFileName("passInfo.txt"),actpass);
    ssdMaxCovReached(buildDefaultCheckpointFileName("maxCovReached.txt"));
    ssdBannedOverlaps(buildDefaultCheckpointFileName("bannedOverlaps.txt"));
  }
  catch(...){
    // Oooops? Error while writing new checkpoint data? Rollback ...
    try{
      boost::filesystem::remove_all(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not delete failed checkpoint directory. Contact the author.\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not delete failed checkpoint directory. Contact the author.\nnUnspecified exception?");
    }
    try{
      boost::filesystem::rename(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp,AS_miraparams[0].getDirectoryParams().dir_checkpoint);
    }
    catch(boost::filesystem::filesystem_error fse){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not rollback checkpoint rename after failed snapshot?\nError message is " << fse.what());
    }
    catch(...){
      MIRANOTIFY(Notify::FATAL,"Now, this is embarassing: could not rollback checkpoint rename after failed snapshot?\nUnspecified exception?");
    }
    MIRANOTIFY(Notify::FATAL,"Could not correctly write new snapshot data. Disk full? Changed permissions?");
  }
  try{
    // now move the files which are not written anew
    // then remove checkpoint tmp¨
    if(fileExists(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp+"/static_hashstat.bin")){
      fileRename(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp+"/static_hashstat.bin",buildDefaultCheckpointFileName("static_hashstat.bin"));
    }
    boost::filesystem::remove_all(AS_miraparams[0].getDirectoryParams().dir_checkpoint_tmp);
  }
  catch(...){
  }

  if(as_fixparams.as_dateoutput) dateStamp(cout);

  FUNCEND();
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadSnapshotData(uint32 & actpass)
{
  FUNCSTART("void Assembly::loadSnapshotData(uint32 & actpass)");

  // readpool was loaded before, don't load here
  actpass=lsdPassInfo(buildDefaultCheckpointFileName("passInfo.txt"));
  lsdMaxCovReached(buildDefaultCheckpointFileName("maxCovReached.txt"));
  lsdBannedOverlaps(buildDefaultCheckpointFileName("bannedOverlaps.txt"));

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::ssdReadPool(const string & filename)
{
  FUNCSTART("void Assembly::ssdReadPool(const string & filename)");
  ofstream fout(filename,ios::out|ios::trunc);
  Contig::dumpMAF_Head(fout);
  AS_readpool.dumpAs(fout,Read::AS_MAF,true);
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot readpool?");
  }
}

void Assembly::ssdPassInfo(const string & filename, uint32 actpass)
{
  FUNCSTART("void Assembly::ssdPassInfo(const string & filename, uint32 actpass)");
  ofstream fout(filename.c_str(),ios::out|ios::trunc);
  fout << actpass << endl;
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot actpass?");
  }
}

void Assembly::ssdMaxCovReached(const string & filename)
{
  FUNCSTART("void Assembly::ssdMaxCovReached(const string & filename)");
  ofstream fout(filename.c_str(),ios::out|ios::trunc);
  for(auto mc : AS_maxcoveragereached) fout << mc << endl;
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot maxcov?");
  }
}

void Assembly::ssdBannedOverlaps(const string & filename)
{
  FUNCSTART("void Assembly::ssdBannedOverlaps(const string & filename)");
  // 1st line: size of AS_permanent_overlap_bans
  // nth line: id1 id2 id2 id2 ...
  ofstream fout(filename.c_str(),ios::out|ios::trunc);
  fout << AS_permanent_overlap_bans.size() << endl;
  for(size_t rid=0; rid < AS_permanent_overlap_bans.size(); ++rid){
    if(!AS_permanent_overlap_bans.bop[rid].empty()){
      fout << rid;
      for(auto id2 : AS_permanent_overlap_bans.bop[rid]){
	fout << '\t' << id2;
      }
      fout << '\n';
    }
  }
  fout.close();
  if(fout.fail()){
    MIRANOTIFY(Notify::FATAL,"Could not write snapshot banned overlaps?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

uint32 Assembly::lsdPassInfo(const string & filename)
{
  FUNCSTART("uint32 Assembly::lsdPassInfo(const string & filename)");

  ifstream fin(filename.c_str(), ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }
  string tmp;
  fin >> tmp;
  auto n=atoll(tmp.c_str());
  if(n<0){
    MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
  }
  return static_cast<uint32>(n);
}

void Assembly::lsdMaxCovReached(const string & filename)
{
  FUNCSTART("uint32 Assembly::lsdMaxCovReached(const string & filename)");

  ifstream fin(filename.c_str(), ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }
  AS_maxcoveragereached.clear();
  string tmpline;
  while(getline(fin,tmpline)){
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    auto n=atoll(tmpline.c_str());
    AS_maxcoveragereached.push_back(static_cast<uint32>(n));
  }
}

void Assembly::lsdBannedOverlaps(const string & filename)
{
  FUNCSTART("void Assembly::lsdBannedOverlaps(const string & filename)");

  ifstream fin(filename.c_str(), ios::in);
  if(!fin.good()) {
    MIRANOTIFY(Notify::FATAL,"Did not find " << filename);
  }

  size_t bopsize=0;
  string tmpline;
  vector<string> substrs;
  if(getline(fin,tmpline)){
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty first line in " << filename << " ??");
    }
    boost::split(substrs, tmpline, boost::is_any_of(" \t"));
    if(substrs.size()>1){
      MIRANOTIFY(Notify::FATAL,"first line in " << filename << " should have one element only");
    }
    auto n=atoll(substrs[0].c_str());
    if(n<0){
      cout << "Faulty first line, has negative values:\n" << tmpline << endl;
      MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
    }
    bopsize=static_cast<size_t>(n);
  }else{
    MIRANOTIFY(Notify::FATAL,"Error reading first line of " << filename);
  }
  if(bopsize==0){
    MIRANOTIFY(Notify::FATAL,"Ooooooops, first line should have positive number: " << filename);
  }
  if(bopsize!=AS_readpool.size()){
    MIRANOTIFY(Notify::FATAL,"Ooooooops, : size of read pool in " << filename  << " (" << bopsize << ") is not equal to size of current readpool (" << AS_readpool.size() << ") ???");
  }
  AS_permanent_overlap_bans.nuke();
  AS_permanent_overlap_bans.resize(bopsize);
  vector<uint32> tmpvec;
  while(getline(fin,tmpline)){
    substrs.clear();
    boost::trim(tmpline);
    if(tmpline.empty()){
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    boost::split(substrs, tmpline, boost::is_any_of(" \t"));
    if(substrs.size()<2){
      cout << "Faulty line, need more elements:\n" << tmpline << endl;
      MIRANOTIFY(Notify::FATAL,"empty line in " << filename << " is not expected");
    }
    auto sI=substrs.begin();
    size_t rid=0;
    {
      auto n=atoll(sI->c_str());
      if(n<0){
	cout << "Faulty line, has negative first value:\n" << tmpline << endl;
	MIRANOTIFY(Notify::FATAL,"Ooooooops, first value should value >= 0: " << filename);
      }
      rid=static_cast<size_t>(n);
    }
    ++sI;
    for(; sI!=substrs.end(); ++sI){
      auto n=atoll(sI->c_str());
      if(n<0){
	cout << "Faulty line, has negative values:\n" << tmpline << endl;
	MIRANOTIFY(Notify::FATAL,"negative value in " << filename << " is not expected");
      }
      tmpvec.push_back(static_cast<uint32>(n));
    }
    AS_permanent_overlap_bans.bop[rid].swap(tmpvec);
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveResumeDataFPO(int32 version, const string & prefix, const string & postfix)
{
  FUNCSTART("void Assembly::saveResumeDataFPO()");

  string wellconnectedfilename(buildFileName(version,
					 prefix,
					 postfix,
					 AS_miraparams[0].getAssemblyParams().as_tmpf_wellconnected,
					 ".bin"));
  if(!saveVector(AS_wellconnected,wellconnectedfilename)){
    MIRANOTIFY(Notify::FATAL, "Error while writing file " << wellconnectedfilename << ". Disk full? Changed permissions?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadResumeDataFPO(int32 version, const string & prefix, const string & postfix)
{
  FUNCSTART("void Assembly::loadResumeDataFPO()");

  string wellconnectedfilename(buildFileName(version,
					 prefix,
					 postfix,
					 AS_miraparams[0].getAssemblyParams().as_tmpf_wellconnected,
					 ".bin"));
  if(!loadVector(AS_wellconnected,wellconnectedfilename,0)){
    MIRANOTIFY(Notify::FATAL, "Error while reading file " << wellconnectedfilename << ". Is the file present and correct? Are permissions right?");
  }
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveResumeDataMA(int32 version, const string & prefix, const string & postfix)
{
  FUNCSTART("void Assembly::saveResumeDataMA()");

  //  need to save current overlap bans to tmp as makeAlignments could have created new ones
  ssdBannedOverlaps(buildFileName(version, "", "_pass",
				  AS_miraparams[0].getAssemblyParams().as_tmpf_banned_overlaps,
				  ".txt"));

  string filename(buildFileName(version, "", "_pass",
				AS_miraparams[0].getAssemblyParams().as_tmpf_istroublemaker,
				".bin"));
  bool allok=false;
  if(saveVector(AS_istroublemaker,filename)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_needalloverlaps,
			   ".bin");
    if(saveVector(AS_needalloverlaps,filename)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_multicopies,
			     ".bin");
      if(saveVector(AS_multicopies,filename)){
	filename=buildFileName(version, "", "_pass",
			       AS_miraparams[0].getAssemblyParams().as_tmpf_hasmcoverlap,
			       ".bin");
	if(saveVector(AS_hasmcoverlaps,filename)){
	  allok=true;
	}
      }
    }
  }
  if(!allok){
    MIRANOTIFY(Notify::FATAL, "Error while reading file " << filename << ". Is the file present and correct? Are permissions right?");
  }
}

/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::loadResumeDataMA(int32 version, const string & prefix, const string & postfix)
{
  FUNCSTART("void Assembly::loadResumeDataMA()");
  //  need to load banned overlaps
  lsdBannedOverlaps(buildFileName(version, "", "_pass",
				  AS_miraparams[0].getAssemblyParams().as_tmpf_banned_overlaps,
				  ".txt"));

  string filename(buildFileName(version, "", "_pass",
				AS_miraparams[0].getAssemblyParams().as_tmpf_istroublemaker,
				".bin"));
  bool allok=false;
  if(loadVector(AS_istroublemaker,filename,0)){
    filename=buildFileName(version, "", "_pass",
			   AS_miraparams[0].getAssemblyParams().as_tmpf_needalloverlaps,
			   ".bin");
    if(loadVector(AS_needalloverlaps,filename,0)){
      filename=buildFileName(version, "", "_pass",
			     AS_miraparams[0].getAssemblyParams().as_tmpf_multicopies,
			     ".bin");
      if(loadVector(AS_multicopies,filename,0)){
	filename=buildFileName(version, "", "_pass",
			       AS_miraparams[0].getAssemblyParams().as_tmpf_hasmcoverlap,
			       ".bin");
	if(loadVector(AS_hasmcoverlaps,filename,0)){
	  allok=true;
	}
      }
    }
  }
  if(allok){
    // checks?
  }else{
    MIRANOTIFY(Notify::FATAL, "Error while reading file " << filename << ". Is the file present and correct? Are permissions right?");
  }
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::saveResults()
{
  FUNCSTART("void Assembly::saveResults()");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  saveDebrisList();
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  saveStatistics();
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  saveReadTagList();
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  saveConsensusTagList();
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  saveContigReadList();
  if(as_fixparams.as_dateoutput) dateStamp(cout);
  if(as_fixparams.as_output_caf){
    saveAsCAF();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_maf){
    saveAsMAF();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_wiggle){
    saveAsWiggle();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_gap4da){
    saveAsGAP4DA();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_fasta) {
    saveAsFASTA();
    if(ReadGroupLib::getNumOfStrains()>1) saveStrainsAsFASTAQUAL();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_tcs) {
    saveAsTCS();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  saveSNPList();
  //saveFeatureAnalysis();
  if(as_fixparams.as_output_txt){
    saveAsTXT();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_ace){
    saveAsACE();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }
  if(as_fixparams.as_output_html){
    saveAsHTML();
    if(as_fixparams.as_dateoutput) dateStamp(cout);
  }

  FUNCEND();
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::makeNewReadPoolFromContigs()
{
  FUNCSTART("void Assembly::makeNewReadPoolFromContigs()");


  int32 ccounter=0;
  for(auto cI=AS_contigs.begin(); cI!=AS_contigs.end(); ++cI, ++ccounter){

    auto & cr = cI->getContigReads();

    for(auto crI = cr.begin(); crI!=cr.end(); ++crI){
      if(crI.getORPID() >= 0){
	if(crI->checkRead()){
	  cout << "Precheck failed: " << endl;
	  cout << *crI;
	  MIRANOTIFY(Notify::FATAL, crI->checkRead()) ;
	}

	AS_readpool[crI.getORPID()]=*crI;

	if(AS_readpool[crI.getORPID()].checkRead()){
	  cout << "Postcheck1 failed: " << endl;
	  cout << AS_readpool[crI.getORPID()];
	  MIRANOTIFY(Notify::FATAL, AS_readpool[crI.getORPID()].checkRead()) ;
	}

	const_cast<Read &>(*crI).discard();

	if(AS_readpool[crI.getORPID()].checkRead()){
	  cout << "Postcheck2 failed: " << endl;
	  cout << AS_readpool[crI.getORPID()];
	  MIRANOTIFY(Notify::FATAL, AS_readpool[crI.getORPID()].checkRead()) ;
	}
      }
    }
  }
  AS_contigs.clear();

  FUNCEND();
}



/*************************************************************************
 *
 * returns whether new strong repeat markers (SRMs) were found for
 *  any contig built in any stage
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }

bool Assembly::buildFirstContigs(const int32 passnr, const EDITParameters & eparams, const bool lastpass)
{
  FUNCSTART("void Assembly::buildFirstContigs()");

  CEBUG("BFC: " << passnr << "\t" << lastpass << endl);

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  directory_parameters const & dir_params= AS_miraparams[0].getDirectoryParams();
  edit_parameters const & ed_params= AS_miraparams[0].getEditParams();

  AS_deleteoldresultfiles=true;

  AS_bfcstats.clear();
  AS_bfcstats.resize(2); // 0 for non-rep contigs, 1 for rep contigs

  AS_contigs.clear();
  Contig::setIDCounter(1);
  Contig::resetCERNumbering();

  AS_assemblyinfo.zeroInfo();

  AS_templateguesses.clear();
  if(AS_maxtemplateid>=0){
    AS_templateguesses.resize(AS_maxtemplateid+1);
    cout << "TGS: " << AS_templateguesses.size() << endl;
  }

  vector<Align> aligncache;
  setupAlignCache(aligncache);

  // AS_hasmcoverlaps will be initialised by Pathfinder
  AS_hasmcoverlaps.clear();

  // Keep track of coverage information for contigs so to be able to make
  //  average contig coverage predictions
  vector<vector<uint32> > covperstpercon(ReadGroupLib::SEQTYPE_END);
  vector<uint32> covtotalpercon;

  // initially, overlaps of every read can be reduced
  if(AS_needalloverlaps.size()==0){
    AS_needalloverlaps.resize(AS_readpool.size(),false);
  }

  // filter out all reads that have no overlap at this stage
  // for this, first set all reads to "used" and "singlet", then
  //   re-allow those that have at least one SW ovelap

  AS_used_ids.clear();
  AS_used_ids.resize(AS_readpool.size(),1);
  AS_isdebris.clear();
  AS_isdebris.resize(AS_readpool.size(),DEBRIS_NOOVERLAP);

  // overlay AS_debrisreason
  {
    auto dstI=AS_isdebris.begin();
    auto srcI=AS_debrisreason.begin();
    for(; srcI!=AS_debrisreason.end(); ++dstI, ++srcI){
      if(*srcI) *dstI=*srcI;
    }
  }

  // take back all reads with overlaps
  {
    vector<newedges_t>::const_iterator ceI=AS_confirmed_edges.begin();
    for(;ceI!=AS_confirmed_edges.end();ceI++){
      AS_isdebris[ceI->rid1]=0;
      AS_used_ids[ceI->rid1]=0;
    }
  }

  // Now, go through all the singlets and take them back into
  //  assembly if they have special MIRA tags attached
  // Also remove backbones & rails from the debris list
  //  and set them to unused
  {
    for(uint32 i=0; i< AS_readpool.size(); i++){
      if(AS_readpool[i].hasTag(Read::REA_tagentry_idSRMr)
	 || AS_readpool[i].hasTag(Read::REA_tagentry_idCRMr)
	 || AS_readpool[i].hasTag(Read::REA_tagentry_idWRMr)
	 || AS_readpool[i].hasTag(Read::REA_tagentry_idSAOr)
	 || AS_readpool[i].hasTag(Read::REA_tagentry_idSROr)
	 || AS_readpool[i].hasTag(Read::REA_tagentry_idSIOr)) {
	AS_isdebris[i]=0;
	AS_used_ids[i]=0;
      }
      if(AS_readpool[i].isBackbone()
	 || AS_readpool[i].isRail()) {
	AS_used_ids[i]=0;
	AS_isdebris[i]=0;
      }
    }
  }


  // rollback digitally normalised reads to "used, debris"
  for(uint32 rpi=0; rpi< AS_readpool.size(); ++rpi){
    if(AS_debrisreason[rpi]==DEBRIS_DIGITAL_NORMALISATION){
      AS_isdebris[rpi]=DEBRIS_DIGITAL_NORMALISATION;
      AS_used_ids[rpi]=1;
    }
  }


  // get total length of backbones (if any)
  uint32 totalbblen=0;
  for(list<Contig>::const_iterator cI=AS_bbcontigs.begin(); cI!=AS_bbcontigs.end(); cI++){
    totalbblen+=cI->getContigLength();
  }

  ofstream fout;
  if(as_fixparams.as_tmpf_unused_ids.size()!=0){
    fout.open((dir_params.dir_tmp+"/"+as_fixparams.as_tmpf_unused_ids).c_str(), ios::out);
    fout.close();
  }

  // outside precomputed lowerbound of oedges, for PathFinder
  // is used lateron in constructStepByStep() of the Pathfinder
  //
  // however, every once in a while (every 10k to 100k reads used from
  //  the pool), the overlap edges vector will be compressed (unecessary
  //  edges thrown out) and therefore this vector will be emptied (to be reconstructed
  //  by pathfinder)
  //
  // On mapping 8m Solexas to 70 contigs (4MB), time goes down from 234 minutes to
  //  221 minutes (Core i940)
  //
  // effect of compression not as big as hoped for, but every little bit helps

  vector<vector<newedges_t>::iterator> tmp_lowerbound_oedges;

  // TODO: currently filled by pathfinder, maybe good to have that already in skim?
  AS_hasreptoverlap.clear();
  AS_hasnoreptoverlap.clear();

  // PathFinder object can be created ouside loop and re-used

  PPathfinder qaf(&AS_miraparams,
		  &AS_readpool,
		  &AS_confirmed_edges,
		  &AS_adsfacts,
		  &aligncache,
		  &AS_used_ids,
		  &AS_multicopies,
		  &AS_hasmcoverlaps,
		  &AS_hasreptoverlap,
		  &AS_hasnoreptoverlap,
		  &AS_istroublemaker,
		  &AS_wellconnected,
		  &tmp_lowerbound_oedges,
		  &AS_templateguesses);

  // this vector will hold the read IDs added by pathfinder to contig
  // + backbone + rail IDs
  // used after call to *_constructStepByStep() to set AS_used_ids[] elements
  vector<int32> tmp_ids_in_contig;

  uint32 overlapcompressstepping=AS_readpool.size()/5;
  if(overlapcompressstepping<10000) overlapcompressstepping=10000;
  //if(overlapcompressstepping<100) overlapcompressstepping=100;
  if(overlapcompressstepping>100000) overlapcompressstepping=100000;
  if(overlapcompressstepping>AS_readpool.size()) overlapcompressstepping=AS_readpool.size();
  uint32 nextoverlapcompress=overlapcompressstepping;
  cout << "overlapcompressstepping: " << overlapcompressstepping
       << "\nnextoverlapcompress: " << nextoverlapcompress << endl;

  bool foundSRMs=false;
  uint32 numsingletssincecleanup=0;

  bool maykillintermediatesinmglets=true;
  bool shouldmovesmallclusterstodebris=false;
  for(uint32 st=0;st < AS_miraparams.size(); st++){
    if(AS_miraparams[st].getAssemblyParams().as_savesimplesingletsinproject) maykillintermediatesinmglets=false;
    if(AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig>1) shouldmovesmallclusterstodebris=true;
  }

  if(shouldmovesmallclusterstodebris) bfc_moveSmallClustersToDebris();


#ifdef CLOCK_STEPS2
  timeval tv;
  timeval tvloop;
  timeval tvtotal;
  gettimeofday(&tvtotal,nullptr);
#endif

  //uint32 countingunused=AS_used_ids.size();

  uint32 trackingunused=AS_used_ids.size();
  for(uint32 i=0; i<AS_used_ids.size(); i++){
    if(AS_used_ids[i]) --trackingunused;
  }

  uint32 numcontigs=1;
  // bug: if someone specifically sets as_maxcontigsperpass to 2^32-1, then
  //  this loop never runs.
  for(;trackingunused>0; ++numcontigs){
    CEBUG("bfc 1\n");
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    cout << '\n';

    bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);


//    // REMOVEME: paranoia check
//    for(uint32 ri=0; ri<AS_used_ids.size(); ++ri){
//      BUGIFTHROW(AS_isdebris[ri] && AS_used_ids[ri]==0, "AS_isdebris[ri] && AS_used_ids[ri]==0  for ri " << ri << "\t" << AS_readpool[ri].getName());
//    }



    // jump out of for loop if max number of contigs was reached
    if(as_fixparams.as_maxcontigsperpass>0 && numcontigs==as_fixparams.as_maxcontigsperpass+1) break;

    CEBUG("bfc 2\n");
    if(trackingunused>0){

#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
      tvloop=tv;
#endif

      //// compress the overlap edges when needed
      //if(AS_readpool.size()-unused > nextoverlapcompress){
      //	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      //	cout << "Compressing overlap edges ..."; cout.flush();
      //	vector<newedges_t>::const_iterator srcI=AS_confirmed_edges.begin();
      //	vector<newedges_t>::iterator dstI=AS_confirmed_edges.begin();
      //	for(; srcI != AS_confirmed_edges.end(); ++srcI){
      //	  if(AS_used_ids[srcI->rid1] == 0
      //	     || AS_used_ids[srcI->linked_with] == 0){
      //	    *dstI=*srcI;
      //	    ++dstI;
      //	  }
      //	}
      //	AS_confirmed_edges.resize(dstI-AS_confirmed_edges.begin());
      //	nextoverlapcompress=AS_readpool.size()-unused+overlapcompressstepping;
      //	cout << "done.\n";
      //	if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
      //}

      Contig::setIDCounter(numcontigs);
      // TODO: change wrt multiple MIRAparams
      Contig buildcon(&AS_miraparams, AS_readpool);

      //vector<int8> tmpused=AS_used_ids;

#ifdef CLOCK_STEPS2
      cout << "Timing BFC prelim1: " << diffsuseconds(tv) << endl;
#endif

      CEBUG("bfc 3\n");
      list<Contig>::const_iterator bbContigI=AS_bbcontigs.begin();
      if(numcontigs <= AS_bbcontigs.size()){
	advance(bbContigI,numcontigs-1);
      } else {
	bbContigI=AS_bbcontigs.end();
      }

      CEBUG("bfc 4\n");
      // set newreptmarked to true just to get into the for pass
      //  value is changed within pass
      bool newreptmarked=true;
      bool wasovercalledited=false;
      bool wasmajorovercalledited=false;

      bool mastermayeditovercalls=ed_params.ed_mira_automatic_contic_editing; // may be changed in one of the iter
      bool mayeditovercalls=mastermayeditovercalls; // recomputed after each contig build, just default init here

      // build a contig, repetitively until
      //  ... maximum number of iterations has passed
      //  ... or no now repeats were marked
      //  ... or no 454 edits happened in last iteration

      // Note: if we´re in last pass, iterate at least once if
      //  454 edits were made! (disregarding as_fixparams.as_numrmbbreakloops,
      //  maxiter gets adapted in loop then

      uint32 maxiter=as_fixparams.as_numrmbbreakloops;

      // Note: assemblies with a lot of passes (>=4) will
      //  get only one loop in the first pass. The reason: first pass already
      //  discovers a lot of repeats that are dealt with later on
      //  in makeAlignments(). It's faster to let makeAlignments() deal
      //  with "discovering" and marking reads than to loop here.

      if(as_fixparams.as_mark_repeats){
	if(as_fixparams.as_numpasses >= 4 && passnr==1){
	  maxiter=1;
	}
	// same thing: passes >= 6, pass 2 leads to maxiter = 2
	if(as_fixparams.as_numpasses >= 6 && passnr==2
	  && maxiter > 2){
	  maxiter=2;
	}
      }

      list<Contig::pbdse_t> pbdsev;

      bool markrepeatsduringstore=true;
      bool contignotok=false;
      bool continueiter=true;

#ifdef CLOCK_STEPS2
      gettimeofday(&tv,nullptr);
#endif
      buildcon.discard();
#ifdef CLOCK_STEPS2
      cout << "Timing BFC discard con: " << diffsuseconds(tv) << endl;
#endif

      string basename_forextsave;
      {
	ostringstream ostr;
	ostr << dir_params.dir_tmp << "/miratmp.pass_" << passnr << "_cb" << numcontigs << "_";

	basename_forextsave=ostr.str();
      }

      bool contigtrimmed=false;
      CEBUG("bfc 5\n");
      //for(uint32 iter=0; iter < maxiter && continueiter; ++iter){
      uint32 iter=0;
      do{
	string basename_forextsave_iter;
	{
	  ostringstream ostr;
	  ostr << basename_forextsave << "i" << iter << "_";

	  basename_forextsave_iter=ostr.str();
	}

	CEBUG("bfc 6/"<<iter << '\n');
	if(iter>0) bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);

	if(AS_hasbackbones
	   && passnr>=as_fixparams.as_startbackboneusage_inpass
	   && bbContigI != AS_bbcontigs.end()) {

#ifdef CLOCK_STEPS2
	  gettimeofday(&tv,nullptr);
#endif

	  // if we're in an iter>0, the contig is not empty
	  // as, we're re-initialising with fresh bb contig,
	  //  we need to do some housekeeping
	  if(buildcon.getContigLength()>0){
	    // release all reads
	    //Contig::setCoutType(Contig::AS_TEXT);
	    //cout << "Throwing away\n" << buildcon;
	    auto & cr=buildcon.getContigReads();
	    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	      if(pcrI.getORPID()>=0){
		BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"bbrebuild ! AS_used_ids[...] ???");
		AS_used_ids[pcrI.getORPID()]=0;
		++trackingunused;
	      }
	    }
	  }

	  // new contig: initialise what's needed
	  buildcon=*bbContigI;
	  buildcon.setContigID(numcontigs);

	  {
	    //Contig::setCoutType(Contig::AS_TEXT);
	    //cout << "registering \n" << buildcon;

	    // track the reads
	    auto & cr=buildcon.getContigReads();
	    for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	      if(pcrI.getORPID()>=0){
		BUGIFTHROW(AS_used_ids[pcrI.getORPID()],"register AS_used_ids[" << pcrI.getORPID() << "]==" << static_cast<uint16>(AS_used_ids[pcrI.getORPID()]) << " ???");
		AS_used_ids[pcrI.getORPID()]=1;
		//--trackingunused;
	      }
	    }
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC copy bbcon: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  // re-initialising the baselocks is necessary as reads might
	  //  have got new SRMc/WRMc tags in previous iterations
	  //  these are not known in the initial backbone contig, so
	  //  they must be made known
	  // TODO: 11.10.2012 not sure whether still good
	  buildcon.initialiseBaseLocks();

	  // tell contig to use backbone characters when possible for
	  //  tmp consensus
	  // TODO: make configurable?
	  buildcon.useBackbone4TmpConsensus(true);

	  // and by default try to merge short reads (Solexa, SOLiD)
	  buildcon.mergeNewSRReads(true);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC bbsetup remain: " << diffsuseconds(tv) << endl;
#endif

	} else {
	  // do nothing
	}

	CEBUG("bfc 7/"<<iter << '\n');

	buildcon.setContigNamePrefix(
	  AS_miraparams[0].getContigParams().con_nameprefix);

	if(iter==0){
	  cout << "Building new contig " << numcontigs;
	  if(buildcon.getNumBackbones()) cout << " from backbone " << buildcon.getContigName();
	  cout << endl;
	}else{
	  cout << "Rebuilding contig " << numcontigs << endl;
	}

	//Contig::setCoutType(Contig::AS_DEBUG);
	//cout << con;

	if(!AS_coverageperseqtype.empty()){
	  cout << "Setting contig coverage targets to: ";
	  for(uint8 ii=0; ii<AS_coverageperseqtype.size(); ii++){
	    cout << '\t' << AS_coverageperseqtype[ii];
	  }
	  cout << endl;
	  buildcon.setContigCoverageTarget(AS_coverageperseqtype);
	}

	if(as_fixparams.as_dateoutput) dateStamp(cout);

	cout << "Unused reads: " << trackingunused << endl;
	cout.flush();

	CEBUG("bfc 8/"<<iter << '\n');

	bfc_callPathfinder(passnr,iter,trackingunused,shouldmovesmallclusterstodebris,
			   buildcon,qaf);

	CEBUG("bfc 9/"<<iter << '\n');

	if(buildcon.getNumReadsInContig()==0){
	  cout << "WARNING WARNING WARNING! no reads in contig?!?!?!" << endl;
	  trackingunused=0;
	  break;
	}

	bool contigmeetsrequirements=bfc_checkIfContigMeetsRequirements(buildcon);
	if(!contigmeetsrequirements){
	  newreptmarked=false;

	  // this contig won't be taken, but contig numbering would
	  //  go up in next loop.
	  // Countermeasure: decrease contig number by one to have the
	  //  loop increase it afterwards, effectivly re-using same contig number
	  --numcontigs;

	  uint32 numdebs=0;
	  for(auto crI=buildcon.getContigReads().begin();crI!=buildcon.getContigReads().end();++crI){
	    if(crI.getORPID() >= 0){
	      BUGIFTHROW(AS_isdebris[crI.getORPID()]>0,"Ooooops, read is already debris? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
	      BUGIFTHROW(!AS_used_ids[crI.getORPID()]>0,"Debris from an unused read in contig? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
	      ++numdebs;
	      AS_isdebris[crI.getORPID()]=DEBRIS_TINYCONTIG;
	    }
	  }
	  cout << "Discarding" << endl;
	  buildcon.discard();

	  cout << "\nContig does not meet requirement of minimum reads per contig."
	    "\nMoved " << numdebs << " reads to debris." << endl;

	  // break out of the iter loop
	  break;
	}

	cout << "\n\nFinished building." << endl;

	if(buildcon.getContigLength()>5 && buildcon.getNumReadsInContig()>1){
	  if(as_fixparams.as_put_asswithmira_tags){
	    buildcon.addTagToConsensus(0,
				       4,
				       '=',
				       "MIRA",
				       "Assembled with MIRA",
				       false);
	  }
	}

	try {
	  if(buildcon.getNumReadsInContig()>1){
	    if(as_fixparams.as_dateoutput) dateStamp(cout);
	    if(buildcon.getContigLength()>100000){
	      cout << "Calculating statistics (this may take a while)." << endl;
	    }

#ifdef CLOCK_STEPS2
	    gettimeofday(&tv,nullptr);
#endif
	    // see whether we can define new multicopies
	    if(as_fixparams.as_automatic_repeat_detection){
	      buildcon.analyseReadCoverage(AS_maxcoveragereached,
				      AS_multicopies,
				      AS_coverageperseqtype);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC analysereadcov: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif

	    buildcon.setCoutType(Contig::AS_TEXT);
	    buildcon.stats(cout);

#ifdef CLOCK_STEPS2
	    cout << "Timing BFC cout constats: " << diffsuseconds(tv) << endl;
#endif
	    if(as_fixparams.as_dateoutput) dateStamp(cout);
	  }
	}
	catch (Notify n) {
	  n.setGravity(Notify::FATAL);
	  n.handleError(n.tif);
	}
	catch (...) {
	  cerr << "Darn, error with that contig. See darn.fasta.\n";
	  Read::setCoutType(Read::AS_CLIPPEDFASTA);
	  for(auto & cre : buildcon.getContigReads()) {
	    cout << cre;
	  }
	  abort();
	}

	CEBUG("bfc 10/"<<iter << '\n');

	// saving pre-edit
	saveExtTmpContig(buildcon,(basename_forextsave_iter+"pre"));

	bool wasedited=false;
	newreptmarked=false;

#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif


	bool goforit=true;
/*
  PacBio trials with strobed reads ... MIRA throws on real, non-strobed data
  furthermore, not sure whether this is really doable (accuracy ~85% *sigh*)

	pbdsev.clear();
	if(buildcon.hasSeqTypeData(ReadGroupLib::SEQTYPE_PACBIOHQ) && !(lastpass && iter==maxiter-1)){
	  uint32 maxcorrect=buildcon.createPacBioDarkStrobeEdits(pbdsev);

	  CEBUG("PacBio maxcorrect: " << maxcorrect << endl);

	  if(maxcorrect>5) {
	    goforit=false;
	    wasmajorovercalledited=true;
	  }

	  // another iteration is needed to get a valid contig
	  //if(maxiter<2) maxiter=2;
	}
*/

	mayeditovercalls= mastermayeditovercalls & buildcon.hasEditableOvercallData();

	if(goforit){
	  // handling of tricky overcalls (normally 454 & Ion)
	  vector<bool> readsmarkedsrm;

	  // mark areas with tricky overcalls
	  // in those areas, no repeat marker may be set
	  // do this always
	  if(mayeditovercalls && buildcon.getNumReadsInContig() >1){
	    cout << "P " << passnr << ", marked " << buildcon.editTrickyOvercalls(true,false,readsmarkedsrm) << " reads.\n";
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC edit tricky1: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  markrepeatsduringstore=true;
	  if(as_fixparams.as_mark_repeats
	     && !as_fixparams.as_mark_repeats_onlyinresult
	     && buildcon.getNumReadsInContig() >1){
	    Contig::repeatmarker_stats_t contigrepstats;
	    newreptmarked=markRepeats(buildcon, readsmarkedsrm, contigrepstats);
	    AS_bfcstats[buildcon.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
	    foundSRMs|=newreptmarked;
	    markrepeatsduringstore=false;
	  }

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC mark reps: " << diffsuseconds(tv) << endl;
	  gettimeofday(&tv,nullptr);
#endif

	  CEBUG("bfc 11/"<<iter << '\n');

	  // edit only when no misassembled repeats found
	  //  this is to prevent the editor to try something foolish on
	  //  misassembled things
	  // same applies to overcall editing
	  //

	  wasovercalledited=false;

//	  if(buildcon.getNumReadsInContig() <= 8) mayeditovercalls=false;
//	  cout << "\nXXXXXXXXXXXXXX " << newreptmarked
//	       << " " << mayeditovercalls
//	       << " " << buildcon.getNumReadsInContig() << endl;
	  if(!newreptmarked && mayeditovercalls && buildcon.getNumReadsInContig() >1){
	    uint32 nummarks=buildcon.editTrickyOvercalls(false,false,readsmarkedsrm);
	    cout << "Edited " << nummarks << " reads.\n";
	    AS_bfcstats[buildcon.getLongRepeatStatus()].numeditovercall+=nummarks;
	    wasovercalledited=nummarks>0;

	    // major edit is edits in >=5% of the reads
	    wasmajorovercalledited=(nummarks>0 && nummarks >= buildcon.getNumReadsInContig()/20);

	    if(wasovercalledited){
	      AS_needsskimfornastyrepeats=true;
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      cout << "Deleting superfluous gap columns (1) ... ";
	      cout.flush();
	      auto tmpnumdel=buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done, deleted " << tmpnumdel << " columns.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC edit tricky2: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif
	  }

	  // BaCh 29.04.2011
	  // However, editSingleDiscrepancyNoHAFTag() is conservative enough to allow
	  //  that kind of edits all the time
	  //

	  if(ed_params.ed_mira_automatic_contic_editing
	     && ed_params.ed_kmer_singlets
	     && buildcon.getNumReadsInContig() >4){
	    // editmode: increasingly less conservative (but still conservative enough)
	    uint32 editmode=passnr-1;
	    if(editmode>3) editmode=3;
	    // careful, editmode 3 not suited for EST/RNASeq!
	    if(!AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
	       && editmode>2) {
	      editmode=2;
	    }
	    uint32 numedits=buildcon.editSingleDiscrepancyNoHAFTag(readsmarkedsrm,editmode);
	    cout << "\nEdited " << numedits << " positions.\n";
	    if(numedits>0){
	      AS_bfcstats[buildcon.getLongRepeatStatus()].numedithashfreq+=numedits;
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      AS_needsskimfornastyrepeats=true;
	      cout << "Deleting superfluous gap columns (2) ... ";
	      cout.flush();
	      auto tmpnumdel=buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done, deleted " << tmpnumdel << " columns.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    }
#ifdef CLOCK_STEPS2
	    cout << "Timing BFC edit single discrepancy, no HAF: " << diffsuseconds(tv) << endl;
	    gettimeofday(&tv,nullptr);
#endif
	  }


/*
  Trials with original, raw PB reads (non-CCS, non-CLR), currently on hold

	  // PacBio LQ edits
	  if(buildcon.getNumReadsInContig() >4 && buildcon.hasSeqTypeData(ReadGroupLib::SEQTYPE_PACBIOLQ)){
	    uint32 numcoledits=0;
	    uint32 numreadedits=0;
	    buildcon.editPBSledgeHammer(readsmarkedsrm,numcoledits,numreadedits);
	    //buildcon.deleteStarOnlyColumns(0,I->getContigLength());
	    if(numreadedits){
	      cout << "PacBio low quality SledgeHammer edits: " << numcoledits << " c-edits, " << numreadedits << " r-edits.\n";
	    }
	  }
*/

	  // if we're in last pass and 454 / pacbio dark strobe edits were made, loop at least
	  //  once (switching off the 454 / pacbio edits for next pass)
	  if(lastpass
	     && wasovercalledited
	     && iter==maxiter-1){
	    if(maxiter>=1) {
	      mastermayeditovercalls=false;
	    }
	    maxiter++;
	  }

	  CEBUG("bfc 12/"<<iter << '\n');

#ifdef CLOCK_STEPS2
	  gettimeofday(&tv,nullptr);
#endif
	  // get rid of all PSHP tags
	  buildcon.deleteTagsInReads(Read::REA_defaulttag_PSHP.identifier);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC delPSHP: " << diffsuseconds(tv) << endl;
#endif

	  CEBUG("bfc 13/"<<iter << '\n');

#ifdef MIRA_HAS_EDIT
	  if(!newreptmarked && ed_params.ed_automatic_contic_editing!=0){
	    cout << "Editing temporary contig: ";
	    if (buildcon.getNumReadsInContig() >1) {
	      cout << endl;
#ifdef CLOCK_STEPS2
	      gettimeofday(&tv,nullptr);
#endif
	      wasedited=true;

	      editContigBack(buildcon, const_cast<EDITParameters &>(eparams));

	      //ScfBuffer::statistics();
	      //ScfBuffer::show();
	      //ScfBuffer::discard();
	      //ScfBuffer::statistics();
	      //ScfBuffer::show();

	      cout << "done.\n";
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      cout << "Deleting superfluous gap columns ... ";
	      cout.flush();
	      buildcon.deleteStarOnlyColumns(0, buildcon.getContigLength()-1);
	      cout << "done.\n";

#ifdef CLOCK_STEPS2
	      cout << "Timing BFC editconback and more: " << diffsuseconds(tv) << endl;
#endif
	      Contig::setCoutType(Contig::AS_TEXT);
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	      buildcon.stats(cout);
	      if(as_fixparams.as_dateoutput) dateStamp(cout);
	    } else {
	      cout << "(only 1 read in contig, no editing)\n";
	    }
	  }
#endif


	  CEBUG("bfc 14/"<<iter << '\n');

	  // saving again if rept-marked or edited
	  if(newreptmarked || wasedited || wasovercalledited) {
	    if(wasovercalledited) {
	      saveExtTmpContig(buildcon,(basename_forextsave_iter+"post454"));
	    } else {
	      saveExtTmpContig(buildcon,(basename_forextsave_iter+"post"));
	    }
	  } else {
	    if ( as_fixparams.as_output_exttmp_fasta
		 || as_fixparams.as_output_exttmp_ace
		 || as_fixparams.as_output_exttmp_gap4da
		 || as_fixparams.as_output_exttmp_caf) {
	      cout << "No edit and no new repeat found, not saving extra temporary contig again.\n";
	    }
	  }
	}


	CEBUG("bfc 15/"<<iter << '\n');

	CEBUG("bfc 16/"<<iter << '\n');

	// Transfer all the reads fron the new contig into readpool
	//  if no misassembly was detected (or 454 editing happened)
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	if(!newreptmarked || wasovercalledited || !pbdsev.empty()) {
	  transferContigReadsToReadpool(buildcon, pbdsev, passnr);

#ifdef CLOCK_STEPS2
	  cout << "Timing BFC rp transfer: " << diffsuseconds(tv) << endl;
#endif
	} else {
	  cout << " Transfering read tags to readpool." << endl;
	  transferContigReadTagsToReadpool(buildcon, bbContigI);
#ifdef CLOCK_STEPS2
	  cout << "Timing BFC crtag2rp transfer: " << diffsuseconds(tv) << endl;
#endif
	}
	cout << "Done." << endl;

	CEBUG("bfc 17/"<<iter << '\n');

	contigtrimmed=false;
	if(passnr>1){
	  // reason: on highly repetitive data, reads trimmed away have SRMr tags,
	  // leading to long assembly times in later contigs. Not trimming away the
	  // reads in first pass makes the 1st pass finish quicker and arrive to the
	  // Smith-Waterman repeat disambiguation stage (SRMr->CRMr tags)
	  contignotok=bfc_trimDenovoIfNecessary(buildcon,foundSRMs,basename_forextsave_iter,trackingunused);
	  if(contignotok) ++AS_bfcstats[buildcon.getLongRepeatStatus()].numdisassemblies;
	  contigtrimmed=true;
	}

	++iter;

	continueiter=(iter < maxiter) && (newreptmarked | wasmajorovercalledited | contignotok);

	CEBUG("I have newreptmarked " << newreptmarked << " wasmajorovercalledited " << wasmajorovercalledited << " contignotok " << contignotok << "\tcontinueiter: " << continueiter << endl);

	bool nukeexistingcontig=false;
	if(newreptmarked){
	  cout << "Identified misassembled reads in contig.\n";
	}
	if(wasmajorovercalledited){
	  cout << "Had many overcall edits in reads.\n";
	  nukeexistingcontig=true;
	}
	if(!contigmeetsrequirements){
	  BUGIFTHROW(true,"We should never be here: contigmeetsrequirements");
	}

	if(continueiter && buildcon.getNumBackbones()>0){
	  cout << "Has backbones, no additional iteration allowed.\n";
	  continueiter=false;
	}

	if(lastpass) nukeexistingcontig=false;

	if(continueiter){
	  cout << "Need to loop contig building\n";
	  if(nukeexistingcontig){
	    cout << "Iteration will nuke contig built so far.\n";
	    //auto & cr=buildcon.getContigReads();
	    //for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	    //  if(pcrI.getORPID()>=0){
	    //	BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"nec ! AS_used_ids[...] ???");
	    //	AS_used_ids[pcrI.getORPID()]=0;
	    //	trackingunused+=1;
	    //  }
	    //}

	    for(auto & rid : qaf.getRIDsKnownInContig()){
	      if(rid>=0
		 && !AS_readpool[rid].isBackbone()
		 && !AS_readpool[rid].isRail()){
		if(AS_used_ids[rid]){
		  trackingunused+=1;
		  AS_used_ids[rid]=0;
		}
	      }
	    }

	    buildcon.discard();
	  }else{
	    cout << "Iteration will keep contig built so far.\n";
	  }
	}

      }while(continueiter);

      bfc_sanityCheckASUSEDIDS(trackingunused,numcontigs);

      // no contig? Then it was discarded, restart building one completely anew
      if(buildcon.getNumReadsInContig() == 0) continue;

      // this here to handle cases repeats should be marked only in results
      // TODO: check whether not to remove this parameter / option at all
      if(pbdsev.empty()
	 && lastpass
	 && as_fixparams.as_mark_repeats
	 && as_fixparams.as_mark_repeats_onlyinresult
	 && buildcon.getNumReadsInContig() >1){
	CEBUG("bfc 18\n");
#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	vector<bool> dummy;

	Contig::repeatmarker_stats_t contigrepstats;
	newreptmarked=markRepeats(buildcon, dummy,contigrepstats);
	AS_bfcstats[buildcon.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
	foundSRMs|=newreptmarked;
	markrepeatsduringstore=false;
#ifdef CLOCK_STEPS2
	cout << "Timing BFC markrep during store: " << diffsuseconds(tv) << endl;
#endif
      }

      if(buildcon.getNumReadsInContig() >0){
	if(newreptmarked){
	  cout << "\nAccepting probably misassembled contig, ";
	  if(contigtrimmed) {
	    cout << " but kept only best, non-problematic part.\n";
	  }else{
	    cout << " keeping as is.\n";
	  }
	}

	CEBUG("bfc 19\n");

#ifdef CLOCK_STEPS2
	gettimeofday(&tv,nullptr);
#endif
	bfc_storeContig(buildcon,numcontigs,markrepeatsduringstore,passnr,lastpass);
#ifdef CLOCK_STEPS2
	cout << "Timing BFC store con: " << diffsuseconds(tv) << endl;
#endif

	{
	  if(buildcon.getContigLength()>=5000){
	    cout << "Contig coverage analysis ";
	    if(buildcon.getContigLength()>=100000){
	      cout << "(this may take a while) ";
	    }
	    cout << "... "; cout.flush();
	    const Contig::constats_t & cs=buildcon.getStats();
	    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
	      covperstpercon[st].push_back(static_cast<uint32>(cs.avg_covperst[st]+.5));
	      //cout << "cccccccccccc: " << covperstpercon[st].back() << endl;
	    }
	    covtotalpercon.push_back(static_cast<uint32>(cs.avg_coverage+.5));
	  }
	}

	// did the caller to assemble() ask for a callback for each contig built?
	// if we're on the last past, then call it
	if(lastpass && AS_contigbuilt_callbackfunc!=nullptr){
	  (*AS_contigbuilt_callbackfunc)(buildcon, AS_readpool);
	}

	if(as_fixparams.as_spoilerdetection) {
	  if(!as_fixparams.as_spdetect_lastpassonly
	     || passnr==as_fixparams.as_numpasses-1) {
	    huntSpoilSports(buildcon);
	  }
	}
      }

      // if we were in mapping mode, on the last contig and not also building
      //  new contigs, move all remaining reads to debris

      if( !as_fixparams.as_backbone_alsobuildnewcontigs
	  && bbContigI != AS_bbcontigs.end()){
	++bbContigI;  // don't bother to decrease after using this ... will be re-init in next loop anyway
	if(bbContigI == AS_bbcontigs.end()){
	  uint32 adddebris=0;
	  for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
	    if(AS_used_ids[uid]==0){
	      ++adddebris;
	      --trackingunused;
	      AS_isdebris[uid]=DEBRIS_NOTMAPPED;
	      AS_used_ids[uid]=1;
	    }
	  }
	  cout << "Last backbone mapped, not building new ones.\nMoved " << adddebris << " remaining reads to debris.\n";
	}
      }


    }  // if(trackingunused>0) ...


    //unused=0;
#ifdef CLOCK_STEPS2
    cout << "Timing BFC loop total: " << diffsuseconds(tvloop) << endl;
#endif
  }   //for(;trackingunused>0; ++numcontigs) ...


  // force a sanity check
  bfc_sanityCheckASUSEDIDS(trackingunused,0);

  // unused reads? we jumped out of contig creation
  // define unused reads as debris

  if(trackingunused>0){
    for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
      if(AS_used_ids[uid]==0){
	AS_isdebris[uid]=DEBRIS_ABORTEDCONTIGCREATION;
	AS_used_ids[uid]=1;
      }
    }
  }

  if(lastpass) {
    saveDebrisList();
  }else{
    saveDebrisList(passnr, "", "_pass");
  }

  cout << "\nBuildstats - RM positions        :\t"
       << AS_bfcstats[0].numnewsrm << '\t' << AS_bfcstats[1].numnewsrm
       << "\nBuildstats - overcall edits      :\t"
       << AS_bfcstats[0].numeditovercall << '\t' << AS_bfcstats[1].numeditovercall
       << "\nBuildstats - hash edits          :\t"
       << AS_bfcstats[0].numedithashfreq << '\t' << AS_bfcstats[1].numedithashfreq
       << "\nBuildstats - contig disassemblies:\t"
       << AS_bfcstats[0].numdisassemblies << '\t' << AS_bfcstats[1].numdisassemblies
       << endl << endl;

  analyseTemplateGuesses();
  if(lastpass && passnr==1){
    analyseTemplateGuesses();
  }

  // reads that are debris or singlets apparently need every chance
  //  they can get to align, therefore subsequent passes should not
  //  reduce the overlaps
  // Change: for Solexa and Ion reads ... sorry, we'll just miss out
  {
    for(uint32 i=0; i<AS_needalloverlaps.size(); i++){
      if(AS_isdebris[i] && AS_readpool[i].getSequencingType() != ReadGroupLib::SEQTYPE_SOLEXA
	 && AS_isdebris[i] && AS_readpool[i].getSequencingType() != ReadGroupLib::SEQTYPE_IONTORRENT) AS_needalloverlaps[i]=true;
    }

    list<Contig>::const_iterator clI=AS_contigs.begin();
    vector<int32> contigids;
    for(; clI!=AS_contigs.end(); clI++){
      if(clI->getNumReadsInContig()==1){
	clI->getReadORPIDsAtContigPosition(contigids,0,0);
	AS_needalloverlaps[contigids[0]]=true;
      }
    }
  }

  // calculate median of average contig coverage
  // TODO: this is calculated only once because of problems
  //  with decreasing average/median in subsequent passes
  //  see if it can be improved
  if(as_fixparams.as_uniform_read_distribution
     && passnr+1>=as_fixparams.as_urd_startinpass
     && !covperstpercon[0].empty() && AS_coverageperseqtype.empty()){
    AS_coverageperseqtype.clear();

    cout << "Setting coverage analysis values for uniform read distribution:\n";
    for(uint32 st=0; st<ReadGroupLib::SEQTYPE_END; st++){
      sort(covperstpercon[st].begin(),covperstpercon[st].end());
      AS_coverageperseqtype.push_back(covperstpercon[st][covperstpercon[st].size()/2]);
      cout << '\t' << ReadGroupLib::getNameOfSequencingType(st) << " coverage:\t" << AS_coverageperseqtype.back() << '\n';
    }
    sort(covtotalpercon.begin(),covtotalpercon.end());
    AS_coveragetotal=covtotalpercon[covtotalpercon.size()/2];
  }

  AS_steps[ASCONTIGSOK]=1;

  AS_used_ids.clear();

  //  saveAsCAF();

  FUNCEND();

  return foundSRMs;
}
#define CEBUG(bla)



/*************************************************************************
 *
 * numexpected: number of reads expected to be unused
 * numcontigs is a simple "don't do that all the time" switch, but every 100th contig
 *
 *************************************************************************/

void Assembly::bfc_sanityCheckASUSEDIDS(uint32 numexpected, uint32 numcontigs)
{
  // REMOVEME from production code once stable

  FUNCSTART("void Assembly::bfc_sanityCheckASUSEDIDS(uint32 numexpected)");

  if(AS_miraparams[0].getSpecialParams().mi_extra_flag1 || numcontigs%100==0){
#ifdef CLOCK_STEPS2
    timeval tv;
    gettimeofday(&tv,nullptr);
#endif

    uint32 count=0;
    for(auto a : AS_used_ids) if(!a) ++count;

#ifdef CLOCK_STEPS2
    cout << "Timing BFC unused: " << diffsuseconds(tv) << endl;
    cout << "CUnused: " << count << endl;
    cout << "TUnused: " << numexpected << endl;
    cout << "AS_used_ids.size(): " << AS_used_ids.size() << endl;
#endif

    BUGIFTHROW(count!=numexpected,"count " << count << " != expected " << numexpected);
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::bfc_callPathfinder(const int32 passnr, const uint32 iter, uint32 & trackingunused, bool shouldmovesmallclusterstodebris,Contig & buildcon, PPathfinder & qaf)
{
  FUNCSTART("void Assembly::bfc_callPathfinder(const int32 passnr, const uint32 iter, uint32 trackingunused, bool shouldmovesmallclusterstodebris,Contig & buildcon, PPathfinder & qaf)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

#ifdef CLOCK_STEPS2
  timeval tv;
#endif

  // in which mode are we? de-novo or mapping?
  bool assemblymode_mapping=false;
  if(AS_hasbackbones
     && passnr >= AS_miraparams[0].getAssemblyParams().as_startbackboneusage_inpass
     && ! AS_miraparams[0].getAssemblyParams().as_backbone_alsobuildnewcontigs){
    assemblymode_mapping=true;
  }

  cout << iter << "\tKnown 1: " << qaf.getRIDsKnownInContig().size() << endl;
  if(iter==0
     || buildcon.getContigReads().size()==0
     || buildcon.getContigLength()==0){
    CEBUG("iter 0, PF init new contig\n");
    qaf.prepareForNewContig(buildcon);
  }else if(assemblymode_mapping){
    CEBUG("iter n, mapping, PF resync contig\n");
    qaf.resyncContig();
  }
  cout << "Known 2: " << qaf.getRIDsKnownInContig().size() << endl;

  CEBUG("assemblymode_mapping: " << assemblymode_mapping << '\n');

  bool wantbootstrap=false;
  //if(assemblymode_mapping){
  //  for(uint32 st=0; st<AS_seqtypespresent.size(); ++st){
  //    if(AS_seqtypespresent[st] && AS_miraparams[st].getAssemblyParams().as_backbone_bootstrapnewbackbone) wantbootstrap=true;
  //  }
  //}
  if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA] && AS_miraparams[ReadGroupLib::SEQTYPE_SOLEXA].getAssemblyParams().as_backbone_bootstrapnewbackbone) wantbootstrap=true;

  if(assemblymode_mapping && AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]){
    if(wantbootstrap){
      CEBUG("mapping & solexa bootstrap\n");
      bfc_cp_mapWithSolexa(buildcon,qaf);

      cout << "Looking at what to throw away ... "; cout.flush();
      priv_removePotentiallyWrongBaseInserts(buildcon);

      cout << "stripping ... ";cout.flush();
      buildcon.stripToBackbone();
      for(auto & rid : qaf.getRIDsKnownInContig()){
	if(rid>=0
	   && !AS_readpool[rid].isBackbone()
	   && !AS_readpool[rid].isRail()){
	  AS_used_ids[rid]=0;
	}
      }
      // chomp is needed here:
      //  removing reads from the contig may leave overhangs at the ends which are not
      //  covered by backbone. The alignment routines in Contig::addRead_wrapped() will not
      //  cope well with that as the calculation of the expected offset is then wrong
      //  (correct for indels *in* the reference, but not made for pseudo-indels at the
      //  ends of the contig)
      cout << "done, chomping ... ";cout.flush();
      buildcon.chompFront(-1);
      buildcon.chompBack(-1);

      cout << "done, synching ... ";cout.flush();
      qaf.resyncContig();
      cout << "done\n";
    }

    bfc_cp_mapWithSolexa(buildcon,qaf);

    // TODO: hack until there's a routine that only clears tags set by the
    //  makeIntelligentConsensus() functions.
    buildcon.clearConsensusTags();

    CEBUG("TU before " << trackingunused << endl);
    CEBUG("RIDs known " << qaf.getRIDsKnownInContig().size() << endl);
    trackingunused-=qaf.getRIDsKnownInContig().size();
    CEBUG("TU after " << trackingunused << endl);

  }else if(assemblymode_mapping){
    CEBUG("bfccp2" << endl);
    qaf.map();
    buildcon.coutAddReadTimings();
    cout << "Known 3: " << qaf.getRIDsKnownInContig().size() << endl;
    trackingunused-=qaf.getRIDsKnownInContig().size();
  }else{
#ifdef CLOCK_STEPS2
    gettimeofday(&tv,nullptr);
#endif

    CEBUG("use general pathfinder: " << AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms << '\n');

    bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

    auto hacktracking=buildcon.getNumReadsInContig();
    if(buildcon.getNumReadsInContig()==0){
      CEBUG("bfccp3" << endl);
    }else{
      CEBUG("no bfccp3" << endl);
	  // so if the contig has some reads already? We're in an iter loop
	  // bump trackingunused by number of reads because after calling PPathfinder,
	  //  they'll be in the number of reads the pathfinder knows of.
    }

    qaf.denovo();

#ifdef CLOCK_STEPS2
    cout << "Timing BFC paf construct: " << diffsuseconds(tv) << endl;
#endif
    buildcon.coutAddReadTimings();

    cout << iter << "\tKnown 3: " << qaf.getRIDsKnownInContig().size() << endl;
    trackingunused+=hacktracking;
    trackingunused-=qaf.getRIDsKnownInContig().size();

    bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

    // Trigger moving small clusters to debris if the pathfinder had to refill its
    //  start cache or when the start cache containes only singlets
    if(shouldmovesmallclusterstodebris
       && (qaf.startCacheRanDry() || qaf.startCacheHasSinglets())){
      cout << "Triggering additional cluster check:";
      if(shouldmovesmallclusterstodebris) cout << " shouldmovesmallclusterstodebris";
      if(qaf.startCacheRanDry()) cout << " startCacheRanDry";
      if(qaf.startCacheHasSinglets()) cout << " startCacheHasSinglets";
      cout << '\n';
      trackingunused-=bfc_moveSmallClustersToDebris();
    }
  }

  bfc_sanityCheckASUSEDIDS(trackingunused, buildcon.getContigID());

  // REMOVEME: paranoia check
  for(auto crI=buildcon.getContigReads().begin();crI!=buildcon.getContigReads().end();++crI){
    if(crI.getORPID() >= 0){
      BUGIFTHROW(!AS_used_ids[crI.getORPID()],"Ooooops, read is in contig but not used? " << crI.getORPID() << " " << AS_readpool[crI.getORPID()].getName() << endl);
    }
  }
}
#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::bfc_cp_mapWithSolexa(Contig & buildcon, PPathfinder & qaf)
{
  FUNCSTART("void Assembly::bfc_cp_mapWithSolexa(Contig & buildcon, PPathfinder & qaf)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();

  bool hasotherst=false;
  //AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]){
  for(uint8 st=0; st<AS_seqtypespresent.size(); ++st){
    if(st!=ReadGroupLib::SEQTYPE_SOLEXA && AS_seqtypespresent[st]) hasotherst=true;
  }

  buildcon.setSpecialSRAddConditions(-1,-1,-1);
  qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_SOLEXA);
  for(uint32 coel=28; coel>=4; coel-=4){
    cout << "Gogo: coel " << coel << endl;
    qaf.setWantsCleanOverlapEnds(coel);
    qaf.setMinTotalNonMatches(1);
    qaf.map();
    buildcon.updateBackboneConsensus();
  }
  qaf.setWantsCleanOverlapEnds(0);
  qaf.setMinTotalNonMatches(0);

  CEBUG("bfccp1" << endl);

  cout << "Gogo: 100% mapping\n";
  buildcon.setSpecialSRAddConditions(0,0,0);
  qaf.map();
  buildcon.coutAddReadTimings();

  if(hasotherst){
    cout << "Gogo: add others clean ends\n";
    qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_END);
    qaf.setWantsCleanOverlapEnds(16);
    qaf.setMinTotalNonMatches(0);
    // SRAddcondition is still 0,0,0 and cleanOverlap of 16 was also
    //  already done for Solexas, so no new Solexa should get added here
    qaf.map();
  }

  qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_SOLEXA);
  qaf.setWantsCleanOverlapEnds(0);
  qaf.setMinTotalNonMatches(0);
  if(AS_miraparams[0].getAlignParams().al_solexahack_maxerrors>0 && qaf.getReadAddAttempts()>0) {
    cout << "Gogo: mapping 1 mismatch\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    buildcon.setSpecialSRAddConditions(1,0,1);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();

    cout << "Gogo: mapping 1 gap\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    buildcon.setSpecialSRAddConditions(1,-1,0);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();
  }

  if(AS_miraparams[0].getAlignParams().al_solexahack_maxerrors>1 && qaf.getReadAddAttempts()>0) {
    cout << "Gogo: mapping 2 mismatches\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    buildcon.setSpecialSRAddConditions(2,0,2);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();

    cout << "Gogo: mapping 1 gap, 1 mismatch\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    buildcon.setSpecialSRAddConditions(2,1,1);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();

    cout << "Gogo: mapping 2 errors (==remaining 2 gaps)\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();
  }

  for(uint32 numerr=3; numerr<=AS_miraparams[0].getAlignParams().al_solexahack_maxerrors; numerr++){
    if(qaf.getReadAddAttempts()==0) break;
    cout << "Gogo: mapping all " << numerr << " errors\n";
    if(as_fixparams.as_dateoutput) dateStamp(cout);
    buildcon.setSpecialSRAddConditions(numerr,-1,-1);
    qaf.map();
    buildcon.coutAddReadTimings();
    buildcon.updateBackboneConsensus();
  }

  if(hasotherst){
    cout << "Gogo: mapping whatever left\n";
    qaf.setAllowedSeqTypeForMapping(ReadGroupLib::SEQTYPE_END);
    qaf.map();
    buildcon.updateBackboneConsensus();
  }

}
#define CEBUG(bla)



/*************************************************************************
 *
 * Move clusters smaller than wished minimum number of reads per contig
 *  to debris
 * However, not rails and not backbones!
 *
 * returns number of reads pushed to debris
 *
 *************************************************************************/

 //#define CEBUG(bla)   {cout << bla; cout.flush(); }
uint32 Assembly::bfc_moveSmallClustersToDebris()
{

  cout << "Moving small clusters to debris:\n";

  uint32 totaldebris=0;

  vector<uint32> numreadsperst(ReadGroupLib::getNumSequencingTypes(),0);
  vector<int32> clusteridperread;
  vector<list<int32> > readinclusterlist;

  clusterUnassembledReads(clusteridperread,readinclusterlist, AS_used_ids);

  uint32 clustercount=0;

  for(size_t ricli=0; ricli<readinclusterlist.size(); ricli++){
    if(!readinclusterlist[ricli].empty()){
      uint32 totalreadsincluster=0;
      numreadsperst.clear();
      numreadsperst.resize(ReadGroupLib::getNumSequencingTypes(),0);
      {
	list<int32>::const_iterator rI=readinclusterlist[ricli].begin();
	for(; rI != readinclusterlist[ricli].end(); ++rI) {
	  numreadsperst[AS_readpool[*rI].getSequencingType()]++;
	  totalreadsincluster++;
	}
      }
      bool takecluster=false;
      for(size_t st=0; st<numreadsperst.size(); ++st){
	if(numreadsperst[st]>0
	   && totalreadsincluster>=AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig){
	  takecluster=true;
	}
      }
      // the above also kill mapping reads
      // therefore, first check whether an unsued rail is part of that cluster
      //  if yes, then don't kill cluster
      if(!takecluster){
	list<int32>::const_iterator rI=readinclusterlist[ricli].begin();
	for(; rI != readinclusterlist[ricli].end(); ++rI) {
	  if(AS_readpool[*rI].isRail()
	     || AS_readpool[*rI].isBackbone()){

	    // should have this ... need to rework how as_usedids is filled (only after contig is made)
	    //&& !AS_used_ids[*rI]){

	    takecluster=false;
	    break;
	  }
	}
      }
      if(!takecluster){
	CEBUG("Killing cluster: " << ricli);
	list<int32>::const_iterator rI=readinclusterlist[ricli].begin();
	for(; rI != readinclusterlist[ricli].end(); ++rI) {
	  CEBUG(" " << *rI << AS_readpool[*rI].getName());
	  AS_isdebris[*rI]=DEBRIS_TINYCLUSTER;
	  AS_used_ids[*rI]=1;
	  ++totaldebris;
	}
	CEBUG('\n');
      }
    }
  }

  // clusterUnassembledReads() did not return orphans as cluster list
  // therefore, look for unused read ids with no cluster number and also
  //  put them into debris

  for(size_t uid=0; uid<AS_used_ids.size(); ++uid){
    if(AS_used_ids[uid]==0 && clusteridperread[uid]==-1
       && !AS_readpool[uid].isRail()
       && !AS_readpool[uid].isBackbone()){
      CEBUG("Killing orphan: " << uid << endl);
      AS_isdebris[uid]=DEBRIS_TINYCLUSTERORPHAN;
      AS_used_ids[uid]=1;
      ++totaldebris;
    }
  }

  cout << "\nDone. " << totaldebris << " reads moved to debris.\n";

  return totaldebris;
}
//#define CEBUG(bla)



/*************************************************************************
 *
 * Checks if Contig meets specified requirements (atm: num of reads)
 *
 * If not, mark the reads in the contig as debris and empty the contig
 *
 *************************************************************************/

bool Assembly::bfc_checkIfContigMeetsRequirements(Contig & con)
{
  FUNCSTART("bool Assembly::bfc_checkIfContigMeetsRequirements(Contig & con)");

  bool contigok=false;

  vector<uint32> numreadsperst(ReadGroupLib::getNumSequencingTypes(),0);
  uint32 totalreadsincon=0;

  auto crI=con.getContigReads().begin();
  for(; crI!=con.getContigReads().end(); ++crI){
    if(crI.getORPID() >= 0){
      if(crI->isBackbone()){
	contigok=true;
	break;
      }
      ++numreadsperst[crI->getSequencingType()];
      ++totalreadsincon;
    }
  }

  if(!contigok){
    for(size_t st=0; st<numreadsperst.size(); ++st){
      if(numreadsperst[st]>0
	 && totalreadsincon>=AS_miraparams[st].getAssemblyParams().as_minimum_readspercontig){
	contigok=true;
      }
    }
  }

  return contigok;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::bfc_markRepReads(Contig & con)
{
  multitag_t tagtstR;
  tagtstR.setIdentifierStr("tstR");
  tagtstR.source=multitag_t::MT_tagsrcentry_idMIRA;

  auto & conreads=con.getContigReads();
  auto crI = conreads.begin();

  for(;crI != conreads.end(); crI++){
    if(crI.getORPID() >= 0
       && AS_multicopies[crI.getORPID()]) {
      cout << "xxxxxxxxxxxx mark " << crI.getORPID() << endl;
      Read & nonconstread = const_cast<Read &>(*crI);
      int32 rc=nonconstread.getRightClipoff()-1;
      if(rc<nonconstread.getLeftClipoff()) rc=nonconstread.getLeftClipoff();
      tagtstR.from=nonconstread.getLeftClipoff();
      tagtstR.to=rc;
      nonconstread.addTagO(tagtstR);
    }
  }
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::priv_tmpcheckroutine(Contig & buildcon)
{
  uint32 index=1;
  cout << "ptcr AS_used_ids[" << index << "]=" << static_cast<uint16>(AS_used_ids[index]) << endl;
  auto & cr=buildcon.getContigReads();
  bool foundidx=false;
  for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
    if(pcrI.getORPID()==index) {
      foundidx=true;
      break;
    }
  }
  cout << "Found idx: " << foundidx << endl;
}

/*************************************************************************
 *
 * return true if trimmed due to misassembly
 *
 *************************************************************************/

bool Assembly::bfc_trimDenovoIfNecessary(Contig & buildcon, bool foundSRMs, const string & basename_forextsave, uint32 & trackingunused)
{
  FUNCSTART("bool Assembly::bfc_trimDenovoIfNecessary(Contig & buildcon, bool foundSRMs, uint32 & trackingunused)");

  bool trimmedmisassembly=false;

  // if we have not been mapping:
  //  1) get misassembled parts out by looking at SRMc tags and trimming back to best range
  //  2) for genome assemblies with pairs: pair analysis and break contig at misassembled sites
  //  3) for genome assemblies: coverage analysis and remove reads in overcovered areas
  if(buildcon.getNumBackbones()==0){

    //priv_tmpcheckroutine(buildcon);
    if(foundSRMs){
      auto range=buildcon.findBestNonMisassembledRange();
      if(range.first>=0){
	cout<<"Found misassembly by repeat marker. Best range: " << range.first << ".." << range.second << '\t' << range.second-range.first << endl;
	trimmedmisassembly=true;
	auto & cr=buildcon.getContigReads();
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"srm ! AS_used_ids[...] ???");
	    AS_used_ids[pcrI.getORPID()]=0;
	  }
	}

	saveExtTmpContig(buildcon,(basename_forextsave+"_pretrimsrm"));

	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=cr.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	buildcon.trimContigToRange(range.first,range.second);

	// rewrite the AS_used_ids with the shortened
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    AS_used_ids[pcrI.getORPID()]=1;
	    --trackingunused;
	  }
	}
	cout << "New trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_posttrimsrm"));
      }
    }

    //priv_tmpcheckroutine(buildcon);

    // pair consistency analysis and contig breaking for genome assemblies, denovo
    if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
       && AS_miraparams[0].getSpecialParams().mi_extra_flag3
       && (!AS_hasbackbones || AS_miraparams[0].getAssemblyParams().as_backbone_alsobuildnewcontigs)){
      auto range=buildcon.findBestPairConsistencyRange();
      if(range.first>0 || range.second < buildcon.getContigLength()){
	cout<<"Found misassembly by pair consistency. Best range: " << range.first << ".." << range.second << '\t' << range.second-range.first << endl;
	trimmedmisassembly=true;
	auto & cr=buildcon.getContigReads();
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    BUGIFTHROW(!AS_used_ids[pcrI.getORPID()],"pair consistency ! AS_used_ids[...] ???");
	    AS_used_ids[pcrI.getORPID()]=0;
	  }
	}

	saveExtTmpContig(buildcon,(basename_forextsave+"_pretrimpair"));

	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=cr.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	buildcon.trimContigToRange(range.first,range.second);

	// rewrite the AS_used_ids with the shortened
	for(auto pcrI=cr.begin(); pcrI!=cr.end(); ++pcrI){
	  if(pcrI.getORPID()>=0){
	    AS_used_ids[pcrI.getORPID()]=1;
	    --trackingunused;
	  }
	}
	cout << "New trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_posttrimpair"));
      }
    }

    // priv_tmpcheckroutine(buildcon);

    // coverage analysis and coverage reduction for genome assemblies, denovo
    // note:
    // has a bug for data which was digitally normalised ...
    //  ... and besides, does absolutely not make sense to do this on
    //  digitally normalised data. Therefore, not done if diginorm active.
    if(AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms
       && AS_miraparams[0].getSpecialParams().mi_extra_flag2
       && !AS_miraparams[0].getHashStatisticsParams().hs_apply_digitalnormalisation
       && (!AS_hasbackbones || AS_miraparams[0].getAssemblyParams().as_backbone_alsobuildnewcontigs)){
      saveExtTmpContig(buildcon,(basename_forextsave+"_prered"));

      Contig::ccctype_t avgcovused=AS_coveragetotal;  // may be 0
      coverageinfo_t cinfo;
      vector<uint64> covvals;
      buildcon.collectCoverage(covvals);
      buildcon.calcStatsOnContainer(cinfo,covvals);
      cout << "1st covnum: " << cinfo << endl;

      // TODO: perhaps make this dependend of ratio mean vs stddev ?
      buildcon.calcSecondOrderStatsOnContainer(cinfo,covvals);
      cout << "2nd covnum: " << cinfo << endl;
      if(cinfo.median>2*avgcovused) avgcovused=cinfo.median;
      cout << "Using: " << avgcovused << endl;

      vector<uint8> peakindicator;
      buildcon.findPeaks(avgcovused,peakindicator);
      unordered_set<readid_t> readsremoved;
      buildcon.reduceReadsAtCoveragePeaks(avgcovused,peakindicator,readsremoved);
      cout << "Coverageremove: " << readsremoved.size() << endl;

      // if reads were removed, get the tracking corrected
      if(!readsremoved.empty()){
	for(auto & rid : readsremoved){
	  BUGIFTHROW(!AS_used_ids[rid],"rere ! AS_used_ids[rid] ???");
	  AS_used_ids[rid]=0;
	}
	cout << "Old trackingunused: " << trackingunused<< endl;
	trackingunused+=readsremoved.size();
	cout << "Intermediate trackingunused: " << trackingunused<< endl;

	saveExtTmpContig(buildcon,(basename_forextsave+"_postred"));
      }
    }
    //priv_tmpcheckroutine(buildcon);
  }

  return trimmedmisassembly;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void Assembly::bfc_storeContig(Contig & con, uint32 & numcontigs, const bool mustmarkrepeats, const int32 passnr, const bool lastpass)
{
  FUNCSTART("void Assembly::bfc_storeContig(Contig & con, uint32 & numcontigs, const bool mustmarkrepeats, const int32 passnr, const bool lastpass)");

  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  auto & conreads=con.getContigReads();

  // look whether we store singlets in contig or not

  bool storecontig=true;

  if(con.getNumReadsInContig()==1){
    const assembly_parameters & as_rt_params = AS_miraparams[conreads.begin()->getSequencingType()].getAssemblyParams();

    storecontig=as_rt_params.as_savesimplesingletsinproject;
    if(conreads.begin()->hasTag(Read::REA_tagentry_idSRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idCRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idWRMr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSROr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSAOr)
       || conreads.begin()->hasTag(Read::REA_tagentry_idSIOr)) {
      if(conreads.begin().getORPID() >= 0) AS_needalloverlaps[conreads.begin().getORPID()]=true;
      storecontig|=as_fixparams.as_savetaggedsingletsinproject;
     }
  }

  if(as_fixparams.as_backbone_trimoverhangingreads){
    con.trimMapOverhang();
  }

  // store contig names only for de-novo genome assemblies
  // used for listing large contigs at end of assembly
  string cnameforstore;
  if(!as_fixparams.as_assemblyjob_mapping
     && AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms){
    cnameforstore=con.getContigName();
  }

  // TODO: U13 hat 0 werte im ASCII?
  //Contig::setCoutType(Contig::AS_DEBUG);
  //cout << "Debug in bfc_storeContig()\n" << con;

  if(storecontig){
    if(as_fixparams.as_mark_repeats&&mustmarkrepeats){
      vector<bool> dummy;
      Contig::repeatmarker_stats_t contigrepstats;
      markRepeats(con, dummy,contigrepstats);
      AS_bfcstats[con.getLongRepeatStatus()].numnewsrm+=contigrepstats.numSRMs;
    }

    //bfc_markRepReads(con);

    cout << "Storing contig ... "; cout.flush();
    cout << as_fixparams.as_mark_repeats << mustmarkrepeats;
    if(AS_hasbackbones){
      con.removeRails();

      // TODO: ask contig whether it has mappings
      //if(as_fixparams.as_loadSOLEXA || as_fixparams.as_loadSOLID){
      if(AS_seqtypespresent[ReadGroupLib::SEQTYPE_SOLEXA]
	 || AS_seqtypespresent[ReadGroupLib::SEQTYPE_ABISOLID]) {
	cout << "Transforming CER mappings." << endl;
	con.transformCERMappingsToCoverageReads();
	cout << "done transforming CER mappings." << endl;
	//assout::saveAsMAF(con, getMAFFilename()+".bla", AS_deleteoldresultfiles);
      }
    }

    con.markFeaturesByConsensus(true, true, true);
    // transfer important tags to readpool
    transferContigReadTagsToReadpool(con,AS_bbcontigs.end());

    con.updateStatsFromConsensusTags(true,true,true,true,true);
    // store the contig information
    AS_assemblyinfo.storeContigStats(con.getStats(),cnameforstore);


    if(lastpass) {
      assout::saveStatistics(con,
			     getStatisticsFilename(),
			     AS_deleteoldresultfiles);
      assout::saveReadTagList(con,
			      getReadTagListFilename(),
			      AS_deleteoldresultfiles);
      assout::saveConsensusTagList(con,
				   getConsensusTagListFilename(),
				   AS_deleteoldresultfiles);
      assout::saveContigReadList(con,
				 getContigReadListFilename(),
				 AS_deleteoldresultfiles);
      if(as_fixparams.as_output_gff3){
	// store sequence for later
//	  AS_gff3defer_names.push_back(con.getContigName());
//TODO: weiterhier
//	  assout::saveTagsAsGFF3(con, getGFF3Filename(), AS_deleteoldresultfiles);
      }
      if(as_fixparams.as_output_caf){
	cout << "Saving CAF ... "; cout.flush();
	assout::saveAsCAF(con, getCAFFilename(), AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_maf){
	cout << "Saving MAF ... "; cout.flush();
	assout::saveAsMAF(con, getMAFFilename(), AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_gap4da){
	cout << "Saving gap4 direct assembly ... "; cout.flush();
	assout::saveAsGAP4DA(con,getGAP4DAFilename(),AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_fasta) {
	if(ReadGroupLib::getNumOfStrains()>1){
	  cout << "Saving strains as FASTA ... "; cout.flush();
	  assout::saveStrainsAsFASTAQ(con, AS_readpool,
				      buildDefaultResultsFileName(
					-1,"","", "",
					AS_miraparams[0].getAssemblyParams().as_outfile_FASTA,
					""),
				      false,
				      0,0,
				      AS_deleteoldresultfiles,false);
	}else{
	  cout << "Saving FASTA ... "; cout.flush();
	  assout::saveAsFASTA(con,
			      getFASTAFilename(),
			      getFASTAPaddedFilename(),
			      AS_deleteoldresultfiles);
	}
	cout << "done.\n";

      }
      if(as_fixparams.as_output_tcs) {
	cout << "Saving TCS ... "; cout.flush();
	assout::saveAsTCS(con, getTCSFilename(),AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_wiggle) {
	cout << "Saving Wiggle ... "; cout.flush();
	assout::saveAsWiggle(con, getWiggleFilename(),AS_deleteoldresultfiles, false);
	cout << "done.\n";
      }
      // TODO: enable these functions for incremental write
      //saveSNPAnalysis();
      //saveFeatureAnalysis();
      if(as_fixparams.as_output_txt){
	cout << "Saving text ... "; cout.flush();
	assout::saveAsTXT(con,getTXTFilename(),AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_ace){
	cout << "Saving ACE ... "; cout.flush();
	assout::saveAsACE(con,getACEFilename(),AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_html) {
	cout << "Saving HTML ... "; cout.flush();
	assout::dumpContigAsHTML(con,
				 getHTMLFilename(),
				 AS_deleteoldresultfiles,
				 AS_miraparams[0].getAssemblyParams().as_projectname_out);
	cout << "done.\n";
      }
    }else{
      assout::saveStatistics(con,
			     getStatisticsFilename(passnr, "", "_pass"),
			     AS_deleteoldresultfiles);
      assout::saveReadTagList(con,
			      getReadTagListFilename(passnr),
			      AS_deleteoldresultfiles);
      assout::saveConsensusTagList(con,getConsensusTagListFilename(passnr),
				   AS_deleteoldresultfiles);
      assout::saveContigReadList(con,
				 getContigReadListFilename(passnr, "", "_pass"),
				 AS_deleteoldresultfiles);

      if(as_fixparams.as_output_tmp_caf) {
	cout << "Saving temp CAF ... "; cout.flush();
	assout::saveAsCAF(con,
			  getCAFFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_maf) {
	cout << "Saving temp MAF ... "; cout.flush();
	assout::saveAsMAF(con,
			  getMAFFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_gap4da) {
	cout << "Saving temp gap4 direct assembly ... "; cout.flush();
	assout::saveAsGAP4DA(con,
			     getGAP4DAFilename(passnr, "", "_pass"),
			     AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_fasta){
	cout << "Saving temp FASTA ... "; cout.flush();
	assout::saveAsFASTA(con,
			    getFASTAFilename(passnr, "", "_pass"),
			    getFASTAPaddedFilename(passnr, "", "_pass"),
			    AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_txt){
	cout << "Saving temp text ... "; cout.flush();
	assout::saveAsTXT(con,
			  getTXTFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_ace) {
	cout << "Saving temp ACE ... "; cout.flush();
	assout::saveAsACE(con,
			  getACEFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      if(as_fixparams.as_output_tmp_tcs) {
	cout << "Saving temp TCS ... "; cout.flush();
	assout::saveAsTCS(con,
			  getTCSFilename(passnr, "", "_pass"),
			  AS_deleteoldresultfiles);
	cout << "done.\n";
      }
      //if(as_fixparams.as_output_tmp_html) saveAsHTML(passnr, "", "_pass");
      if(as_fixparams.as_output_tmp_html) {
	cout << "Saving temp HTML ... "; cout.flush();
	assout::dumpContigAsHTML(con,
				 getHTMLFilename(passnr, "", "_pass"),
				 AS_deleteoldresultfiles,
				 AS_miraparams[0].getAssemblyParams().as_projectname_out);
	cout << "done.\n";
      }
    }
    cout << "done." << endl;
  }else{
    // store the contig information
    AS_assemblyinfo.storeContigStats(con.getStats(),cnameforstore);

    if(conreads.begin().getORPID() >= 0) AS_isdebris[conreads.begin().getORPID()]=DEBRIS_UNSAVEDSINGLET;
    numcontigs--;
  }

  AS_deleteoldresultfiles=false;

  return;
}




/*************************************************************************
 *
 *
 *
 *************************************************************************/


bool Assembly::markRepeats(Contig & con, vector<bool> & readsmarkedsrm, Contig::repeatmarker_stats_t & repstats)
{
  FUNCSTART("void Assembly::markRepeats(Contig & con)");

  repstats.init();

  cout << "Marking possibly misassembled repeats: ";
  cout.flush();

  con.newMarkPossibleRepeats(repstats, readsmarkedsrm);
  cout << "done step 1, starting step 2:";
  con.codonSingleBaseRepeatMarker(6,repstats, readsmarkedsrm);
  if(repstats.numSRMs>0 || repstats.numWRMs>0 || repstats.numSNPs>0){
    cout << "\nFound\n";
    cout << " - " << repstats.numSRMs << " new Strong RMB (SRMc)\n";
    cout << " - " << repstats.numWRMs << " new Weak RMB (WRMc)\n";
    cout << " - " << repstats.numSNPs << " SNP\npositions tagged.";
    cout.flush();
  }else{
    cout << "done. Found none." << endl;
  }

  FUNCEND();

  return repstats.numSRMs>0;
}





/*************************************************************************
 *
 *
 *
 *************************************************************************/

////#define CEBUG(bla)   {cout << bla; cout.flush(); }
//void Assembly::transferContigReadsToReadpool(const Contig & buildcon, list<Contig::pbdse_t> & pbdsev, int32 passnr)
//{
//  FUNCSTART("void Assembly::transferContigReadsToReadpool(const Contig & buildcon, vector<Contig::pbdse_t> & pbdsev)");
//
//  cout << "Transfering reads to readpool." << endl;
//
//  BUGIFTHROW(true,"need redo for PlacedContigReads");
//
//  const vector<Contig::contigread_t> & cr = buildcon.getContigReads();
//
//  // split one list pbdsev into sublists for each contigread if needed
//
//  CEBUG("tCRTR 1" << endl);
//  vector<list<Contig::pbdse_t> > vpbdsev;
//  if(!pbdsev.empty()){
//    vpbdsev.resize(AS_readpool.size());
//    for(; pbdsev.begin() != pbdsev.end(); ){
//      //cout << "ls: " << pbdsev.size() << "\tSplicing: " << pbdsev.front();
//      size_t trid=pbdsev.front().rid;
//      vpbdsev[trid].splice(vpbdsev[trid].begin(),pbdsev,pbdsev.begin());
//      //cout << "ts: " << vpbdsev[trid].size() << endl;
//    }
//
//    //for(uint32 vi=0;vi<vpbdsev.size(); ++vi){
//    //  cout << "vi: " << vi << endl;
//    //  list<Contig::pbdse_t>::const_iterator pI=vpbdsev[vi].begin();
//    //  for(; pI!=vpbdsev[vi].end(); ++pI){
//    //	cout << "\t" << *pI;
//    //  }
//    //}
//  }
//
//  CEBUG("tCRTR 2" << endl);
//  // create a temporary read with enough capacity to hold the
//  //  largest of the reads to transfer (to prevent re-allocation)
//  //
//  // used for PacBio where elastic N stretches need to
//  //  be corrected
//  // also used to remove the gaps from the reads
//  vector<Contig::contigread_t>::const_iterator crI = cr.begin();
//  Read tmpr;
//  {
//    uint32 reservelen=0;
//    for(;crI!=cr.end();crI++){
//      if(crI->orpid>=0){
//	reservelen=max(reservelen,crI->read.getLenSeq());
//      }
//    }
//    tmpr.reserve(reservelen);
//  }
//
//  CEBUG("tCRTR 3" << endl);
//  // now copy all contig reads to read pool, one by one, using the tmp read
//  //  and performing edits on this tmp read, leaving contig reads untouched
//  crI = cr.begin();
//  for(;crI!=cr.end();crI++){
//    if(crI->orpid>=0){
//      tmpr=crI->read;
//
//      // carefull, maybe there is nothing in the vpbdsev vector!
//      //  then one would get segfaults
//      if(tmpr.isSequencingType(ReadGroupLib::SEQTYPE_PACBIOHQ) && !vpbdsev.empty()){
//	list<Contig::pbdse_t>::const_iterator pRI=vpbdsev[crI->orpid].begin();
//	for(; pRI != vpbdsev[crI->orpid].end(); ++pRI){
//	  CEBUG("Apply: " << *pRI);
//	  BUGIFTHROW(pRI->rid >= AS_readpool.size(), "pRI->rid (" << pRI->rid << ") >= AS_readpool.size() ?");
//	  //Read::setCoutType(Read::AS_TEXT);
//	  //cout << CON_reads[pRI->cri].read;
//	  tmpr.correctNStretch(pRI->rposs,
//			       pRI->rpose,
//			       pRI->changeestim);
//	}
//
//	if(passnr==4){
//	  tmpr.transformGapsToNs();
//	}
//      }
//
//      tmpr.removeGapsFromRead();
//      AS_readpool[crI->orpid]=tmpr;
//    }
//  }
//
//  FUNCEND();
//}
////#define CEBUG(bla)


//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::transferContigReadsToReadpool(const Contig & buildcon, list<Contig::pbdse_t> & pbdsev, int32 passnr)
{
  FUNCSTART("void Assembly::transferContigReadsToReadpool(const Contig & buildcon, vector<Contig::pbdse_t> & pbdsev)");

  cout << "Transfering reads to readpool." << endl;

  // create a temporary read with enough capacity to hold the
  //  largest of the reads to transfer (to prevent re-allocation)
  //
  // used for PacBio where elastic N stretches need to
  //  be corrected
  // also used to remove the gaps from the reads
  Read tmpr;
  {
    uint32 reservelen=0;
    for(auto pcrI=buildcon.getContigReads().begin(); pcrI!=buildcon.getContigReads().end(); ++pcrI){
      if(pcrI.getORPID() >= 0){
	reservelen=max(reservelen,pcrI->getLenSeq());
      }
    }
    tmpr.reserve(reservelen);
  }

  CEBUG("tCRTR 3" << endl);
  // now copy all contig reads to read pool, one by one, using the tmp read
  //  and performing edits on this tmp read, leaving contig reads untouched
  for(auto pcrI=buildcon.getContigReads().begin(); pcrI!=buildcon.getContigReads().end(); ++pcrI){
    if(pcrI.getORPID() >= 0){
      tmpr=*pcrI;
      tmpr.removeGapsFromRead();
      AS_readpool[pcrI.getORPID()]=tmpr;
    }
  }

  FUNCEND();
}
//#define CEBUG(bla)




/*************************************************************************
 *
 *
 *
 *************************************************************************/

//#define CEBUG(bla)   {cout << bla; cout.flush(); }
void Assembly::transferContigReadTagsToReadpool(const Contig & con, const list<Contig>::const_iterator bbContigI)
{
  FUNCSTART("void Assembly::transferContigReadTagsToReadpool(const Contig & con, const list<Contig>::const_iterator bbContigI)");

  cout << "Transfering tags to readpool.\n";

  auto & cr = con.getContigReads();
  auto pcrI = cr.begin();
  uint32 tagnumber=0;

  multitag_t tmpmt;

  try{
    // Transfer the contigread RMB tags into the readpool only!
    // Go through all the contigreads, if they have SRMr or WRMr tags,
    //  check if they are at an edited place.
    // If not, transfer the tag to the readpool (if not already present
    //  there)
    for(;pcrI!=cr.end(); ++pcrI){
      if(pcrI.getORPID()==-1) continue;
      uint32 numtags=pcrI->getNumOfTags();

      //CEBUGF(pcrI->getName() << " has " << numtags << " tags.\n");

      for(tagnumber=0; tagnumber < numtags; tagnumber++) {
	const multitag_t & acttag=pcrI->getTag(tagnumber);
	if(acttag.identifier==Read::REA_tagentry_idSRMr
	   ||acttag.identifier==Read::REA_tagentry_idWRMr
	   ||acttag.identifier==Read::REA_tagentry_idSROr
	   ||acttag.identifier==Read::REA_tagentry_idSIOr
	   ||acttag.identifier==Read::REA_tagentry_idSAOr
	  ) {

	  tmpmt=acttag;

	  //CEBUGF("Tag " << tagnumber << " at " << acttag.from << " is " << acttag.identifier << endl);

	  bool foundedit=false;
	  for(uint32 i=0; i<numtags; i++) {
	    if(pcrI->getTag(i).from==acttag.from
	       && (pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_C.identifier
		   || pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_I.identifier
		   || pcrI->getTag(i).identifier==Read::REA_defaulttag_ED_D.identifier
		 )
	      ) {
	      //CEBUGF("Found " << pcrI->getTag(i).identifier << " at that position, skipping!");
	      foundedit=true;
	      break;
	    }
	  }
	  if(foundedit) continue;


	  // the next if is true for tags at gap bases;
	  //if(adj_tagpos < 0 ) {
	  //  cout << "Tag at gap base!!!!!\n";
	  //  continue;
	  //}

	  if(bbContigI!=AS_bbcontigs.end()){
	    if(pcrI->isBackbone() || pcrI->isRail()) {
	      //cout << "IsBackbone or read\n";

	      int32 adj_tagposl=pcrI->getAdjustmentPosOfReadPos(acttag.from);
	      int32 adj_tagposu=pcrI->getAdjustmentPosOfReadPos(acttag.to);
	      // don't do that for gaps in backbones or reads
	      if(adj_tagposl < 0 && adj_tagposu < 0) continue;
	      if(adj_tagposl < 0) {
		adj_tagposl=adj_tagposu-(acttag.to-acttag.from);
		if(adj_tagposl<0) adj_tagposl=0;
	      }else if(adj_tagposu < 0) {
		adj_tagposu=adj_tagposl+(acttag.to-acttag.from);
		if(adj_tagposu > static_cast<int32>(pcrI->getLenSeq()-1)) adj_tagposu=pcrI->getLenSeq()-1;
	      }

	      // This is an ugly, ugly hack: we set the P|WRMB tag in
	      //  the read of the backbonecontig.
	      // This should be a no no, but it's needed to work
	      //  with backbones

	      int32 bbcrtagposl=pcrI->getReadPosOfAdjustmentPos(adj_tagposl);
	      int32 bbcrtagposu=pcrI->getReadPosOfAdjustmentPos(adj_tagposu);
	      if(bbcrtagposl >= 0 && bbcrtagposu>=0) {
		try {
		  tmpmt.from=bbcrtagposl;
		  tmpmt.to=bbcrtagposu;
		  const_cast<Read &>(*pcrI).addTagO(tmpmt);
		}
		catch (Notify n) {
		  cout << "Tried to transfer tags to bbackbone contig read from:\n";
		  Read::setCoutType(Read::AS_TEXT);
		  cout << *pcrI;
		  cout << "Exiting.\n";
		  n.handleError(THISFUNC);
		}
	      }
	    }
	  }

	  Read & rpr=AS_readpool.getRead(pcrI.getORPID());

	  // TODO: FIXME
	  // this is wrong for reads without adjustments.
	  //  -> no way to detect deletions in those reads!
	  int32 adj_tagposl=pcrI->getLowerNonGapAdjustmentPosOfReadPos(acttag.from);
	  int32 adj_tagposu=pcrI->getUpperNonGapAdjustmentPosOfReadPos(acttag.to);
	  int32 rprtagposl=rpr.getReadPosOfAdjustmentPos(adj_tagposl);
	  int32 rprtagposu=rpr.getReadPosOfAdjustmentPos(adj_tagposu);

	  try {
	    //cout << "Transfering tag for " << rpr.getName()
	    //	 << "\t" << rprtagposl << " " << rprtagposu
	    //	 << "\tadj " << adj_tagposl << " " << adj_tagposu << endl;
	    //CEBUGF("Transfering tag for " << rpr.getName() << "\t" << rprtagposl << " " << rprtagposu << endl);
	    //Read::setCoutType(Read::AS_TEXT);
	    //CEBUGF("Before:\n" << rpr << endl);

	    // add tag only when both positions are >=0 (i.e. not starting/stopping
	    //  on an inserted base, else addTag() would understandably barf
	    if(rprtagposl>=0 && rprtagposu>=0){
	      tmpmt.from=rprtagposl;
	      tmpmt.to=rprtagposu;
	      rpr.addTagO(tmpmt);
	    }
	    //CEBUGF("After:\n" << rpr << endl);
	  }
	  catch (Notify n) {
	    // care about errors only if we have adjustments in the read
	    //  if not, the whole thing is more or less wrong anyway
	    if(rpr.usesAdjustments()){
	      cout << "Tried to transfer tags to readpool read from:\n";
	      Read::setCoutType(Read::AS_TEXT);
	      cout << *pcrI;
	      cout << "Exiting." << endl;
	      n.handleError(THISFUNC);
	    }else{
	      cout << "Dbg: wrong tag transfer for " << rpr.getName() << '\n';
	    }
	  }
	}
      }
    }
  }
  catch(Notify n) {
    cout << "General error while transfering tag " << tagnumber << " to readpool read from:\n";
    Read::setCoutType(Read::AS_TEXT);
    cout << *pcrI;
    cout << "Exiting.\n";
    n.handleError(THISFUNC);
  }

  FUNCEND();
}

//#define CEBUG(bla)



/*************************************************************************
 *
 * New version to search for spoil sports
 * Works with one contig at a time
 *
 * Look through all reads. If end of read with (500) bases of end of contig
 *  and has permbans and no freq >4:
 *     throw out (later version perhaps cut back)
 *
 *
 *************************************************************************/
// idea: keep info whether read has been totally embedded (>1000 bases left, right)
// if true, then less likely to be an error

void Assembly::huntSpoilSports(Contig & chkcon)
{
  cout << "Hunting join spoiler" << endl;

  auto & cr=chkcon.getContigReads();

  const uint32 endrange=30;

  auto pcrI=cr.begin();
  for(; pcrI != cr.end() ; ++pcrI){
    const Read & actread=*pcrI;
    bool atfront=false;
    bool atback=false;
    if(pcrI.getReadStartOffset()<=endrange) atfront=true;
    if(pcrI.getReadStartOffset()+actread.getLenClippedSeq() > chkcon.getContigLength()-endrange) atback=true;
    if(atfront || atback){
      cout << "HSS At end: " << atfront << ' ' << atback << '\t' << actread.getName() << endl;

      if(pcrI.getORPID()>=0){
	// currently: always remove

	// if it has >2 permbans, remove
	//cout << "HSS Permbans: " << AS_permanent_overlap_bans[pcrI.getORPID()].size() << endl;
	//if(AS_permanent_overlap_bans[pcrI.getORPID()].size()>3) {
	//  if(!actread.hasTag(Read::REA_tagentry_idSRMr)
	//     && !actread.hasTag(Read::REA_tagentry_idCRMr)){
	//    //if(true){
	//    cout << "HSS remove " << actread.getName() << endl;
	//
	//    const vector<Read::bposhashstat_t> & bposhashstats=actread.getBPosHashStats();
	//    if(!bposhashstats.empty()){
	//      bool clipfront=false;
	//      bool clipback=false;
	//      if(atfront){
	//	if(pcrI->direction>0){
	//	  clipfront=true;
	//	}else{
	//	  clipback=true;
	//	}
	//      }
	//      if(atback){
	//	if(pcrI->direction>0){
	//	  clipback=true;
	//	}else{
	//	  clipfront=true;
	//	}
	//      }
	//      cout << "HSS clip front back " << clipfront << ' ' << clipback << endl;
	//      if(clipfront){
	//	uint32 bposfrom=actread.getLeftClipoff();
	//	uint32 maxcheck=max(actread.getLenSeq()/4,static_cast<uint32>(50));
	//	uint32 maxto=min(maxcheck,actread.getLenSeq()-1);
	//	bool foundinv=false;
	//	for(; bposfrom<maxto; bposfrom++){
	//	  if(!bposhashstats[bposfrom].fwd.isValid()){
	//	    foundinv=true;
	//	    break;
	//	  }else{
	//	    if(bposhashstats[bposfrom].fwd.getFrequency()>=4) break;
	//	  }
	//	}
	//	if(foundinv){
	//	  if(bposfrom-actread.getLeftClipoff()<50){
	//	    bposfrom=actread.getLeftClipoff()+50;
	//	    if(bposfrom>=actread.getLenSeq()) bposfrom=actread.getLenSeq()-1;
	//	  }
	//	  cout << "HSS moving left " <<  AS_readpool[pcrI.getORPID()].getLQClipoff();
	//	  AS_readpool[pcrI.getORPID()].setLQClipoff(bposfrom);
	//	  cout << " to " << AS_readpool[pcrI.getORPID()].getLQClipoff()
	//	       << "\t" << actread.getName() << endl;
	//	}
	//      }
	//      if(clipback){
	//	uint32 bposto=actread.getRightClipoff();
	//	uint32 maxcheck=max(actread.getLenSeq()/4,static_cast<uint32>(50));
	//	bool foundinv=false;
	//	for(; bposto>0 && maxcheck>0; --bposto, --maxcheck){
	//	  if(!bposhashstats[bposto].rev.isValid()){
	//	    foundinv=true;
	//	    break;
	//	  }
	//	  if(bposhashstats[bposto].rev.getFrequency()>=4) break;
	//	}
	//	if(foundinv){
	//	  if(actread.getRightClipoff()-bposto<50){
	//	    if(actread.getRightClipoff()<50){
	//	      bposto=0;
	//	    }else{
	//	      bposto=actread.getRightClipoff()-50;
	//	    }
	//	  }
	//	  cout << "HSS moving right " <<  AS_readpool[pcrI.getORPID()].getRQClipoff();
	//	  AS_readpool[pcrI.getORPID()].setRQClipoff(bposto);
	//	  cout << " to " << AS_readpool[pcrI.getORPID()].getRQClipoff()
	//	       << "\t" << actread.getName() << endl;
	//	}
	//      }
	//    }
	//  }
	//}
	AS_istroublemaker[pcrI.getORPID()]=true;
      }
    }
  }
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::findPossibleOverlaps(int32 version, const string prefix, const string postfix, const string tmpfname)
{
  FUNCSTART("void Assembly::findPossibleOverlaps(int32 version, const string prefix, const string postfix, const string tmpfname)");

  string signalfilename(buildFileName(version,
				  prefix,
				  postfix,
				  AS_miraparams[0].getAssemblyParams().as_tmpf_signal_findpossibleoverlaps,
				  ".ok"));

  if(!AS_resumeasembly || !AS_resumeisok || !fileExists(signalfilename)){
    cout << "AS_resumeasembly " << AS_resumeasembly << endl;
    cout << "AS_resumeisok " << AS_resumeisok << endl;
    cout << "fileExists(" << signalfilename << ") " << fileExists(signalfilename) << endl;

    AS_resumeisok=false;
    fpo_main(version,prefix,postfix,tmpfname);

    // initialiase well connected with overlap criterion levels
    // log overlap criterion levels if wanted

    AS_wellconnected.clear();
    AS_wellconnected.resize(AS_readpool.size(),false);

    {
      ofstream fout;
      //AS_logflag_oclevel=true;
      if(AS_logflag_oclevel){
	string filename=buildFileName(version, "", "",
				      "elog.oclevel_pass",
				      ".lst");
	fout.open(filename.c_str(), ios::out);
      }

      for(uint32 i=0; i<AS_readpool.size();++i){
	if(AS_logflag_oclevel){
	  fout << AS_readpool[i].getName()
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvl[0][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvr[0][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvl[1][i])
	       << '\t' << static_cast<uint16>(AS_overlapcritlevelvr[1][i])
	       << '\n';
	}
	// for Solexa, critlevel 0 is too harsch (because of special rule to calc
	//  Solexa critlevels in Skim)
	for(uint32 ocvi=0; ocvi<2;++ocvi){
	  if(AS_readpool[i].isSequencingType(ReadGroupLib::SEQTYPE_SOLEXA)){
	    if(AS_overlapcritlevelvl[ocvi][i] <= 2 && AS_overlapcritlevelvr[ocvi][i] <= 2 ){
	      AS_wellconnected[i]=true;
	    }
	  }else if(AS_overlapcritlevelvl[ocvi][i] == 0 && AS_overlapcritlevelvr[ocvi][i] == 0 ){
	    AS_wellconnected[i]=true;
	  }
	}
      }
    }
    // save some needed data
    saveResumeDataFPO(version,prefix,postfix);

    reduceSkimHits4(version, prefix, postfix, tmpfname);

    ofstream fout(signalfilename.c_str());  // create checkpoint signal file for findPossibleOverlaps
  }else{
    cout << "Resume assembly: skim and skim reduction already present, good.\n";
    // set AS_pos?match_filename*
    fpo_buildFileNames(version,prefix,postfix,tmpfname);
    AS_posfmatch_full_filename=AS_posfmatch_filename+".reduced";
    AS_poscmatch_full_filename=AS_poscmatch_filename+".reduced";

    loadResumeDataFPO(version,prefix,postfix);
  }
}

void Assembly::fpo_buildFileNames(int32 version, const string prefix, const string postfix, const string tmpfname)
{
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  if(tmpfname.size()){
    AS_posfmatch_filename=buildFileName(version,
					prefix,
					postfix,
					tmpfname+"f",
					".bin");
    AS_poscmatch_filename=buildFileName(version,
					prefix,
					postfix,
					tmpfname+"c",
					".bin");
    // the following two only to also automatically remove
    //  older versions of the reduced skim files
    buildFileName(version,
		  prefix,
		  postfix,
		  tmpfname+"f",
		  ".bin.reduced");
    buildFileName(version,
		  prefix,
		  postfix,
		  tmpfname+"c",
		  ".bin.reduced");
  }else{
    AS_posfmatch_filename=buildFileName(version,
					prefix,
					postfix,
					as_fixparams.as_tmpf_posmatch+"f",
					".bin");
    AS_poscmatch_filename=buildFileName(version,
					prefix,
					postfix,
					as_fixparams.as_tmpf_posmatch+"c",
					".bin");
    // the following two only to also automatically remove
    //  older versions of the reduced skim files
    buildFileName(version,
		  prefix,
		  postfix,
		  as_fixparams.as_tmpf_posmatch+"f",
		  ".bin.reduced");
    buildFileName(version,
		  prefix,
		  postfix,
		  as_fixparams.as_tmpf_posmatch+"c",
		  ".bin.reduced");
  }
}

void Assembly::fpo_main(int32 version, const string prefix, const string postfix, const string tmpfname)
{
  FUNCSTART("void Assembly::fpo_main(int32 version, const string prefix, const string postfix, const string tmpfname)");

  //directory_parameters const & dir_params= AS_miraparams->getDirectoryParams();
  assembly_parameters const & as_fixparams= AS_miraparams[0].getAssemblyParams();
  skim_parameters const & skim_params= AS_miraparams[0].getSkimParams();

  if(as_fixparams.as_dateoutput) dateStamp(cout);
  cout << '\n';

//  if(AS_needsskimfornastyrepeats && skim_params.sk_masknastyrepeats){
//    AS_needsskimfornastyrepeats=false;
//    markNastyReapeatsWithSkim(version, prefix, postfix, tmpfname);
//
//    // test strategy: reads with MRMr must get all overlaps
//    // BaCh 14.03.2009: BAD STRATEGY
//    //  for high coverage (100x +), too many reads may get that set,
//    //  leading to long run times and high memory usage
//    //for(uint32 i=0; i<AS_readpool.size(); i++) {
//    //  if(AS_readpool[i].hasTag(Read::REA_tagMRMr)) AS_needalloverlaps[i]=1;
//    //}
//  }

  cout << "\nSearching for possible overlaps";


#if TRACKMEMUSAGE
    cout << "\ndmi fpo 00\n";
    dumpMemInfo();
#endif

  // save memory for this step, those structures will have to be recomputed anyway
  //nukeSTLContainer(AS_adsfacts);
  //nukeSTLContainer(AS_confirmed_edges);
  AS_adsfacts.clear();
  AS_confirmed_edges.clear();

  nukeSTLContainer(AS_readhitmiss);
  nukeSTLContainer(AS_readhmcovered);
  nukeSTLContainer(AS_count_rhm);

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 05\n";
    dumpMemInfo();
#endif

  {
    vector<uint32> overlapcounter(AS_readpool.size(),0);

    //vector<uint32> rawhashhitcounter;
    //rawhashhitcounter.resize(AS_readpool.size(),0);

    bool onlyagainstrails=false;


    // very first call will be with version=0 ... pre_assembly
    //  don't set to true there as this analysis is then needed
    //  for multicopy analysis
    if(AS_hasbackbones
       && version >= as_fixparams.as_startbackboneusage_inpass
       && ! as_fixparams.as_backbone_alsobuildnewcontigs){
      onlyagainstrails=true;
    }

    if(onlyagainstrails) cout << " (only against backbone, the progress bar will be skewed)";
    cout << ":\n";

    //string rawhitsfilename;
    string megahubtmpfname;
    {
      fpo_buildFileNames(version,prefix,postfix,tmpfname);
      AS_posfmatch_full_filename=AS_posfmatch_filename;
      AS_poscmatch_full_filename=AS_poscmatch_filename;

      if(tmpfname.size()){
	megahubtmpfname=buildFileName(version,
				      prefix,
				      postfix,
				      tmpfname+"_megahubs",
				      ".lst");
      }else{
	megahubtmpfname=buildFileName(version,
				      prefix,
				      postfix,
				      as_fixparams.as_tmpf_posmatch+"_megahubs",
				      ".lst");
      }

      //cout << "Only against rails? " << onlyagainstrails << endl;

      if(!onlyagainstrails) {
	if(skim_params.sk_basesperhash <=12
	   || (skim_params.sk_basesperhash <=14
	       && skim_params.sk_hashsavestepping <3)){
	  cout << "\n\nWARNING!!!!!!\nYou are not performing a 'mapping only' assembly and the parameters"
	    "\n -SK:bph=" << static_cast<uint16>(skim_params.sk_basesperhash) << " and -SK:hss="
	       << static_cast<uint16>(skim_params.sk_hashsavestepping)
	       << "\nare quite low. If SKIM takes ages, stop this assembly and restart while"
	    "\nincreasing these parameters.\n\n";
	}
      }

      vector<int32> overlaplenrequired;
      vector<int32> prrequired;
      for(uint32 i=0;i<ReadGroupLib::getNumSequencingTypes(); i++){
	overlaplenrequired.push_back(AS_miraparams[i].getAlignParams().al_min_overlap);
	prrequired.push_back(AS_miraparams[i].getSkimParams().sk_percentrequired);
      }

      vector<int32> chuntleftcut;
      vector<int32> chuntrightcut;

      // trigger chimera search in pre-assembly skim;
      if(!AS_doneskimchimera &&
	 (as_fixparams.as_clip_skimchimeradetection || as_fixparams.as_clip_skimjunkdetection)){
	chuntleftcut.resize(1);
	chuntrightcut.resize(1);
	AS_doneskimchimera=true;
      }

      AS_chimeracutflag.clear();

      // force taking strong good overlaps only for denovo genome
      bool forcetakestronggood=!as_fixparams.as_assemblyjob_mapping & AS_miraparams[0].getPathfinderParams().paf_use_genomic_algorithms;

      Skim s2;
      s2.setExtendedLog(AS_miraparams[0].getSpecialParams().mi_extended_log);
      uint32 nummegahubs=s2.skimGo(AS_readpool,
				   AS_posfmatch_filename,
				   AS_poscmatch_filename,
				   megahubtmpfname,
				   AS_permanent_overlap_bans,
				   overlapcounter,
				   AS_writtenskimhitsperid,
				   chuntleftcut,
				   chuntrightcut,
				   AS_overlapcritlevelvl,
				   AS_overlapcritlevelvr,
				   skim_params.sk_numthreads,
				   skim_params.sk_maxhashesinmem,
				   onlyagainstrails,
				   skim_params.sk_alsoskimrevcomp,
				   static_cast<uint8>(skim_params.sk_basesperhash),
				   static_cast<uint8>(skim_params.sk_hashsavestepping),
				   prrequired,
				   overlaplenrequired,
				   skim_params.sk_maxhitsperread,
				   forcetakestronggood
	);

      cout << "Total megahubs: " << nummegahubs << endl;

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 10\n";
    dumpMemInfo();
#endif

      if(nummegahubs){
	cout << "\n\nMIRA has detected megahubs in your data."
	  "This may not be a problem, but most probably is, especially for eukaryotes.\n\n";
	if(100.0/AS_num_reads_valid*nummegahubs > skim_params.sk_maxmegahubratio){
	  cout << "\n\nYou have more than " << skim_params.sk_maxmegahubratio << "% of your reads found to be megahubs."
	    "\n\nYou should check the following:\n\n"
	    "\t1) for Sanger sequences: are all the sequencing vectors masked / clipped?\n"
	    "\t2) for 454 sequences: are all the adaptors masked / clipped?\n\n";
	  if(AS_miraparams[0].getHashStatisticsParams().hs_repeatlevel_in_infofile){
	    cout << "You will find in the info directory a file called\n"
	      "    '*_info_readrepeats.lst',\n"
	      "consult the MIRA manual on how to extract repeat information from there.\n\n";
	  }else{
	    cout << "To learn more on the types of repeats you have and how MIRA\n"
	      " can help you to find them, please consult the manual on the\n"
	      " usage of -HS:rliif and the tmp files they create.\n";
	  }
	  cout << "*ONLY* when you are sure that no (or only a very negligible number) of sequencing"
	    "\nvector / adaptor sequence is remaining, try this:\n\n"
	    "\t3) for organisms with complex repeats (eukaryots & some bacteria):\n";
	  if(!AS_miraparams[0].getHashStatisticsParams().hs_masknastyrepeats) cout << "\t\t- use -HS:mnr=yes\n";
	  cout << "\t\t- reduce the -HS:nrr parameter (divide by 2)\n"
	    // huh?
	    //"4) for EST projects, -SK:nrr will not really work, use -SK:nrr (start at\n"
	    //"   10 and increase in steps of of 10)\n"
	    "\n"
	    "*ONLY* if the above fails, try increasing the -SK:mmhr parameter\n"
	    "Note that the number of present megahubs will increase computation time in\n"
	    "an exponential way, so be careful when changing -SK:mmhr.\n";
	}

	if(100.0/AS_readpool.size()*nummegahubs >= skim_params.sk_maxmegahubratio){
	  cout << "\n\nYou have " << 100.0/AS_readpool.size()*nummegahubs
	       << "% of your reads as megahubs.\n"
	       << "You have set a maximum allowed ratio of: " << skim_params.sk_maxmegahubratio
	       << "\n\nEnding the assembly because the maximum ratio has been reached/surpassed.\n";
	  exit(10);
	}
      }

      if(chuntleftcut.size()){
	string cltmpfname;
	cltmpfname=buildFileName(0,"","",
				as_fixparams.as_tmpf_clippings,
				 ".txt",
				 "",
				 false);
	string logprefix="skim detect: ";
	AS_chimeracutflag.resize(1);
	cutBackPossibleChimeras(cltmpfname, logprefix, chuntleftcut,chuntrightcut,AS_chimeracutflag);
      }

    }

    // in mapping assemblies, correct the matches not being 100%
    if(AS_miraparams[0].getSkimParams().sk_swcheckonbackbones
       && AS_hasbackbones){
      recalcNonPerfectSkimMappingsBySW(version);
    }


    //// log the raw hash hits
    //{
    //  ofstream ofs;
    //  ofs.open(rawhitsfilename.c_str(), ios::out| ios::trunc);
    //  for(size_t rhhci=0; rhhci < rawhashhitcounter.size(); rhhci++){
    //	ofs << rhhci << '\t' << rawhashhitcounter[rhhci] << '\n';
    //  }
    //  ofs.close();
    //}

    // Do this only once
    // TODO: check if ok to do more, i.e. if skim can be adapted to
    //       still count banned overlaps.
    if(AS_multicopies.size()==0){
      string filename;
      if(tmpfname.size()){
	filename=buildFileName(version, prefix, postfix, tmpfname+"_multicopystat", ".txt");
      }else{
	filename=buildFileName(version, prefix, postfix,
			       AS_miraparams[0].getAssemblyParams().as_tmpf_posmatch+"_multicopystat",
			       ".txt");
      }
      ofstream fout;
      fout.open(filename.c_str(), ios::out);
      fout.close();
      //flagMulticopyReads(overlapcounter, filename);
    }

    //if(AS_miraparams[0].getAssemblyParams().as_dateoutput) dateStamp(cout);
    //cout << '\n';
    //exit(0);
  }


// TODO: adapt to data in files instead of the old posXmatch multimaps
/*
#if 0
  {
    vector<int32> cluster;
    cluster.resize(AS_readpool.size(),-1);
    uint32 clustercount=0;

    possible_overlaps_t::const_iterator I=AS_posfmatch.begin();
    while(I!=AS_posfmatch.end()){
      int32 cnum1=cluster[I->first];
      int32 cnum2=cluster[I->second.otherid];
      if(cnum1==-1 && cnum2==-1) {
	cluster[I->first]=clustercount;
	cluster[I->second.otherid]=clustercount;
	clustercount++;
      } else if(cnum1==-1) {
	cluster[I->first]=cluster[I->second.otherid];
      } else if(cnum2==-1) {
	cluster[I->second.otherid]=cluster[I->first];
      } else {
	if (cnum1!=cnum2) {
	  // uh oh ... we have to merge both these clusters
	  // simply change all cnum1 into cnum2 in cluster vector
	  for(uint32 j=0; j<AS_readpool.size(); j++) {
	    if(cluster[j]==cnum1) cluster[j]=cnum2;
	  }
	}
      }

      I++;
    }

    I=AS_poscmatch.begin();
    while(I!=AS_poscmatch.end()){
      int32 cnum1=cluster[I->first];
      int32 cnum2=cluster[I->second.otherid];
      if(cnum1==-1 && cnum2==-1) {
	cluster[I->first]=clustercount;
	cluster[I->second.otherid]=clustercount;
	clustercount++;
      } else if(cnum1==-1) {
	cluster[I->first]=cluster[I->second.otherid];
      } else if(cnum2==-1) {
	cluster[I->second.otherid]=cluster[I->first];
      } else {
	if (cnum1!=cnum2) {
	  // uh oh ... we have to merge both these clusters
	  // simply change all cnum1 into cnum2 in cluster vector
	  for(uint32 j=0; j<AS_readpool.size(); j++) {
	    if(cluster[j]==cnum1) cluster[j]=cnum2;
	  }
	}
      }

      I++;
    }

    string filename;
    if(tmpfname.size()){
      filename=buildFileName(version, prefix, postfix, tmpfname+"_pcluster", ".lst");
    }else{
      filename=buildFileName(version, prefix, postfix,
			     AS_miraparams->getAssemblyParams().as_tmpf_posmatch+"_pcluster",
			     ".lst");
    }

    ofstream fout;
    fout.open(filename.c_str(), ios::out);
    uint32 outputcount=0;
    for(uint32 i=0; i<clustercount; i++) {
      bool found=false;
      for(uint32 j=0; j<AS_readpool.size(); j++) {
	if(cluster[j]==static_cast<int32>(i)) {
	  found=true;
	  fout << outputcount << " " << AS_readpool.getRead(j).getName() << endl;
	}
      }
      if(found) outputcount++;
    }
    for(uint32 j=0; j<AS_readpool.size(); j++) {
      if(cluster[j]==-1) {
	fout << "-1 " << AS_readpool.getRead(j).getName() << endl;
      }
    }
  }
#endif
*/

#if TRACKMEMUSAGE
    cout << "\ndmi fpo 20\n";
    dumpMemInfo();
#endif

  AS_steps[ASADSLISTOK]=0;

  FUNCEND();
}


/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

/*
void Assembly::flagMulticopyReads(const vector<uint32> & overlapcounter, const string & tmpfilename)
{
  AS_multicopies.clear();
  AS_multicopies.resize(AS_readpool.size(),0);

  ofstream fout;
  fout.open(tmpfilename.c_str(), ios::out);

  //cout << "Searching for multicopy reads:\n";

  for(uint8 actseqtype=0; actseqtype<ReadGroupLib::SEQTYPE_END; actseqtype++){
    vector<uint32> sortedoverlapcounter=overlapcounter;

    // set overlapcounter of backbones and railreads and read that are
    //  not of currently analysed sequencing type to 0 so as not to
    //  count them
    uint32 numreadsinactseqtype=0;
    for(uint32 i=0; i<overlapcounter.size(); i++){
      if(AS_readpool[i].isRail()
	 || AS_readpool[i].isBackbone()
	 || AS_readpool[i].getSequencingType()!=actseqtype) {
	sortedoverlapcounter[i]=0;
      }else{
	numreadsinactseqtype++;
      }
    }

    if(numreadsinactseqtype==0) continue;

    //cout << numreadsinactseqtype << " in sequencing type " << static_cast<uint16>(actseqtype) << endl;
    sort(sortedoverlapcounter.begin(), sortedoverlapcounter.end());

    // 5% quantil
    uint32 quantilnumber=5*numreadsinactseqtype/100;
    uint32 ifrom=0;
    // well, start the 5% quantil only at reads which have at least
    //  1 overlap
    while(ifrom<sortedoverlapcounter.size()
	  && sortedoverlapcounter[ifrom]==0) ifrom++;

    ifrom+=quantilnumber;
    uint32 ito=static_cast<uint32>(sortedoverlapcounter.size())-quantilnumber;

    if(ito<=ifrom || ifrom>=sortedoverlapcounter.size()){
      ifrom=0;
      ito=static_cast<uint32>(sortedoverlapcounter.size());
    }

    uint32 nonsinglets=0;
    uint32 totaloverlaps=0;
    for(uint32 i=ifrom; i<ito; i++){
      totaloverlaps+=sortedoverlapcounter[i];
      nonsinglets++;
    }


    if(nonsinglets==0) nonsinglets=1;
    uint32 avgoverlaps=static_cast<uint32>(.5+static_cast<double>(totaloverlaps)/static_cast<double>(nonsinglets));
    // strictly speaking, this is not median. But close enough.
    uint32 medianoverlaps=sortedoverlapcounter[ifrom+((ito-ifrom)/2)];
    //uint32 multicopythreshold=avgoverlaps*2;
    uint32 multicopythreshold=medianoverlaps*2;

    fout << "Hitstatistics (" << ReadGroupLib::getShortNameOfSequencingType(actseqtype) << "): nonsinglets:" << nonsinglets << "\ttotaloverlaps: " << totaloverlaps << "\tavgoverlaps: " << avgoverlaps << "\tmedianoverlaps: " << medianoverlaps << endl;

    // now set multicopy flags for affected reads of this sequencing type
    for(uint32 i=0; i<overlapcounter.size(); i++){
      if(AS_readpool[i].isRail()
	 || AS_readpool[i].isBackbone()
	 || AS_readpool[i].getSequencingType()!=actseqtype) continue;
      if(overlapcounter[i]>multicopythreshold) {
	AS_multicopies[i]=static_cast<uint8>(actseqtype+1);
	if(overlapcounter[i]>multicopythreshold*10) {
	  AS_multicopies[i]=static_cast<uint8>(actseqtype+1+100);
	}
      }
    }
  }

  for(uint32 i=0; i<overlapcounter.size(); i++){
    if(AS_readpool[i].isRail()
       || AS_readpool[i].isBackbone()) continue;
    fout << i << "\t" << AS_readpool[i].getName() << "\t" << overlapcounter[i];
    if(AS_multicopies[i]>100) {
      fout << "\tst: " << static_cast<uint16>(AS_multicopies[i]-101);
      fout << "\tmulticopy / insane (forgotten clone vector?)";
    }else if(AS_multicopies[i]>0) {
      fout << "\tst: " << static_cast<uint16>(AS_multicopies[i]-1);
      fout << "\tmulticopy";
    }
    fout << '\n';
  }
  fout.close();
}
*/





/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

const Read & Assembly::getRead(uint32 index)
{
  FUNCSTART("const Read & Assembly::getRead(uint32 index)");

  if(index>=AS_readpool.size()){
    MIRANOTIFY(Notify::INTERNAL,"index: " << index << " greater than AS_readpool.size():" << AS_readpool.size() << "  (out of bounds)");
  }

  FUNCEND();

  return AS_readpool.getRead(static_cast<int32>(index));
}




/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void Assembly::refreshContigAndReadpoolValuesAfterLoading(ReadPool & rp, list<Contig> & contigs)
{
  FUNCSTART("void Assembly::refreshContigAndReadpoolValuesAfterLoading(ReadPool & rp, list<Contig> & contigs)");

  try{
    rp.makeTemplateIDs();

    // Not needed anymore for readgroup approach: RG keeps track of this
    //rp.makeStrainIDs();

    for(auto & cle : contigs) cle.recalcTemplateIDsAndStrainPresent();
  }
  catch(Notify n){
    n.handleError(THISFUNC);
  }
}






/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////        Obsolete         ///////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////


/*************************************************************************
 *
 * dead code, not used anymore?
 *
 *
 *************************************************************************/

//void Assembly::banReadPairGroups(const vector<int32> & g1, const vector<int32> & g2)
//{
//  FUNCSTART("void Assembly::banReadPairGroups(const vector<int32> & g1, const vector<int32> & g2)");
//
//  int32 id1, id2;
//  newedges_t tmpegde;
//  for(uint32 i1=0; i1 < g1.size(); i1++) {
//    id1=g1[i1];
//    for(uint32 i2=0; i2 < g2.size(); i2++) {
//      id2=g2[i2];
//
//      // skip it if there's alreaddy a permban on it
//      if(AS_permanent_overlap_bans.checkIfBanned(id1,id2)) continue;
//
//      if(AS_readpool.getRead(id1).getTemplatePartnerID()==id2) {
//#ifndef PUBLICQUIET
//	cout << "Do not ban template partners "  << id1 << ": " << AS_readpool.getRead(id1).getName() << "\t" << id2 << ": " << AS_readpool.getRead(id2).getName() << endl;
//#endif
//	  continue;
//      }
//
//#ifndef PUBLICQUIET
//      cout << "Banning: " << id1 << ": " << AS_readpool.getRead(id1).getName() << "\t" << id2 << ": " << AS_readpool.getRead(id2).getName() << endl;
//#endif
//
//      // put both ids in permanent overlap banlist
//      AS_permanent_overlap_bans.insertBan(id1,id2);
//
//      // now remove overlap edges between these reads
//      // first in one direction ...
//      //overlap_edges_t::iterator Irun=AS_edges_forward.lower_bound(id1);
//      tmpegde.rid1=id1;
//      vector<newedges_t>::iterator Irun=lower_bound(AS_confirmed_edges.begin(),
//						    AS_confirmed_edges.begin(),
//						    tmpegde,
//						    Assembly__compareNewEdges_t_);
//      while(Irun != AS_confirmed_edges.end()
//	    && Irun->rid1 == id1) {
//	if(Irun->linked_with==id2) {
//	  // erase() doesn't give back an iterator as advertised?
//	  AS_confirmed_edges.erase(Irun);
//	  //Irun=AS_edges_forward.lower_bound(id1);
//	  Irun=lower_bound(AS_confirmed_edges.begin(),
//			   AS_confirmed_edges.begin(),
//			   tmpegde,
//			   Assembly__compareNewEdges_t_);
//	  // Irun points to element after (or end) now, no need to increment
//	  continue;
//	}
//	Irun++;
//      }
//
//      // .. then in the other one
//      //Irun=AS_edges_forward.lower_bound(id2);
//      tmpegde.rid1=id2;
//      Irun=lower_bound(AS_confirmed_edges.begin(),
//		       AS_confirmed_edges.begin(),
//		       tmpegde,
//		       Assembly__compareNewEdges_t_);
//      while(Irun != AS_confirmed_edges.end()
//	    && Irun->rid1 == id2) {
//	if(Irun->linked_with==id1) {
//	  AS_confirmed_edges.erase(Irun);
//	  //Irun=AS_edges_forward.lower_bound(id2);
//	  Irun=lower_bound(AS_confirmed_edges.begin(),
//			   AS_confirmed_edges.begin(),
//			   tmpegde,
//			   Assembly__compareNewEdges_t_);
//	  // Irun points to element after (or end) now, no need to increment
//	  continue;
//	}
//	Irun++;
//      }
//    }
//  }
//
//  FUNCEND();
//}



/////////////////////////////////////////////////////////////////////////
////////////////////        Dead since a while         //////////////////
/////////////////////////////////////////////////////////////////////////
