/*
 * Written by Bastien Chevreux (BaCh)
 *
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

#include "util/fileanddisk.H"
#include "util/misc.H"

#include <fstream>

#ifdef __APPLE__
#include <mach-o/dyld.h>   // for _NSGetExecutablePath
#endif

#include <glob.h>                  // globWalkPath() uses this


#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

#include "errorhandling/errorhandling.H"


using namespace std;
namespace fs = boost::filesystem;


#define CEBUG(bla)



/*************************************************************************
 *
 *
 *
 *************************************************************************/

size_t myFRead(void *ptr, size_t sizememb, size_t nmemb, FILE *stream)
{
  static size_t maxreadchunk=134217728; // 128 MiB

  uint8 * uptr=static_cast<uint8 *>(ptr);

  size_t totalbytestoread=sizememb*nmemb;
  size_t totalnumbytesread=0;
  while(totalbytestoread){
    auto bytestoread=totalbytestoread;
    if(bytestoread>maxreadchunk){
      bytestoread=maxreadchunk;
    }

    auto numbytesread=fread(uptr,1,bytestoread,stream);
    totalnumbytesread+=numbytesread;
    if(numbytesread!=bytestoread) break;
    totalbytestoread-=bytestoread;
    uptr+=bytestoread;
  }

  return (totalnumbytesread/sizememb);
}

size_t myFWrite(const void *ptr, size_t sizememb, size_t nmemb, FILE *stream)
{
  static size_t maxwritechunk=134217728; // 128 MiB

  const uint8 * uptr=static_cast<const uint8 *>(ptr);

  size_t totalbytestowrite=sizememb*nmemb;
  size_t totalnumbyteswritten=0;
  while(totalbytestowrite){
    auto bytestowrite=totalbytestowrite;
    if(bytestowrite>maxwritechunk){
      bytestowrite=maxwritechunk;
    }

    auto numbyteswritten=fwrite(uptr,1,bytestowrite,stream);
    totalnumbyteswritten+=numbyteswritten;
    if(numbyteswritten!=bytestowrite) break;
    totalbytestowrite-=bytestowrite;
    uptr+=bytestowrite;
  }

  return (totalnumbyteswritten/sizememb);
}


int myFSeek(FILE *stream, size_t offset, int whence)
{
  return fseeko(stream, offset,whence);
}

size_t myFTell(FILE *stream)
{
  return ftello(stream);
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void splitFullPathAndFileName(const string & what, string & path, string & filename)
{
  string::size_type bpos = what.rfind("/");

  if (bpos != string::npos) {
    path=what.substr(0,bpos);
  }
  if(bpos!=what.size()) {
    filename=what.substr(bpos+1,what.size());
  }

  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

void findLocationOfSelfBinary(string & location)
{
  location.clear();

  const char * self=nullptr;

#ifdef __APPLE__
//#error "need to implement findLocationOfSelfBinary _NSGetExecutablePath() (man 3 dyld)"
  char path[16384];
  uint32 size = static_cast<uint32>(sizeof(path));
  if (_NSGetExecutablePath(path, &size) == 0){
    self=path;
  }else{
    cout << "Uhhh, could not determine binary path because too small buffer??? Need to program workaround with new/free.\n";
    exit(100);
  }
#else
  const static string selfpse("/proc/self/exe");         // Linux
  const static string selfpcf("/proc/curproc/file");     // FreeBSD / DragonFlyBSD if they have /proc
  const static string selfpce("/proc/curproc/exe");      // NetBSD

  if(boost::filesystem::exists(selfpse)){
    self=selfpse.c_str();
  }else if(boost::filesystem::exists(selfpcf)){
    self=selfpcf.c_str();
  }else if(boost::filesystem::exists(selfpce)){
    self=selfpce.c_str();
  }

#endif

  if(self!=nullptr){
    boost::filesystem::path bpath(self);
    while(boost::filesystem::is_symlink(bpath)) {
      bpath=boost::filesystem::read_symlink(bpath);
    }
    location=bpath.native();
  }else{
    cout << "Uhhh, could not determine binary path on this OS, needs programming!\n";
    exit(100);
  }

  return;
}


/*************************************************************************
 *
 *
 *
 *************************************************************************/

uint64 countLinesInFile(const string & filename)
{
  FUNCSTART("uint32 countLinesInFile(const string & filename)");

  size_t totallines=0;

  ifstream finfin;
  finfin.open(filename.c_str(), ios::in);
  if(!finfin){
    MIRANOTIFY(Notify::FATAL, "File not found: " << filename);
  }

  string dummy;
  while(!finfin.eof()){
    getline(finfin,dummy);
    if(!finfin.eof()) ++totallines;
  }

  finfin.close();

  FUNCEND();
  return totallines;
}



/*************************************************************************
 *
 *
 *
 *
 *************************************************************************/

void dumpFile(const char * fname, ostream & ostr)
{
  ifstream fin;
  fin.open(fname, ios::in);
  if(fin){
    ostr << "Dump from " << fname << "\n"
	 << "--------------------------------------------------------------------------------\n";
    string dummy;
    while(!fin.eof()){
      getline(fin,dummy);
      if(!fin.eof()) ostr << dummy << '\n';
    }
    fin.close();
    ostr << "--------------------------------------------------------------------------------\n";
  }else{
    ostr << "Could not read file " << fname
	 << "\n--------------------------------------------------------------------------------\n";
  }
}




/*************************************************************************
 *
 * returns 0 if directory is there, -1 if not
 *
 *
 *************************************************************************/

int ensureDirectory(const string & dirname, bool purge, bool verbose, bool waituntilok)
{
  struct stat st;
  int rc=stat(dirname.c_str(),&st);
  if(rc || purge) {
    return purgeCreateDir(dirname, verbose, waituntilok);
  }

  return 0;
}



/*************************************************************************
 *
 * returns 0 if everything went fine (old dir evtually deleted and created
 *   anew)
 * -1 if something unexpected happened
 *
 *************************************************************************/

int purgeCreateDir(const string & dirname, bool verbose, bool waituntilok)
{
  auto ret=removeDirectory(dirname,verbose,waituntilok);
  if(ret==0) ret=createDirectory(dirname,verbose,waituntilok);
  return ret;
}



/*************************************************************************
 *
 * returns 0 if everything went fine (old dir evtually deleted and created
 *   anew)
 * -1 if something unexpected happened
 *
 *************************************************************************/

int removeDirectory(const string & dirname, bool verbose, bool waituntilok)
{
  boost::system::error_code ec;

  bool isok=false;
  uint32 waitseconds=1;
  uint32 waititer=0;
  if(boost::filesystem::is_directory(dirname,ec)){
    if(verbose) {
      cout << "Deleting old directory " << dirname << " ... "; cout.flush();
    }
    while(!isok){
      isok=true;
      // remove_all should not throw according to docs when error_code is used,
      //   but does at least for BOOST 1.50 (on OSX? Not on Linux) if directory is write protected
      // so play it safe for all boost::filesystem functions
      try{
	boost::filesystem::remove_all(dirname,ec);
      }
      catch(boost::filesystem::filesystem_error fse){
	isok=false;
	cout << "\n\nCould not delete directory " << dirname << ", error message is: " << fse.what() << endl;
      }
      catch(...){
	isok=false;
	cout << "\n\nCould not delete directory " << dirname << ", unspecified reason." << endl;
      }
      if(isok){
	try{
	  if(boost::filesystem::is_directory(dirname,ec)){
	    isok=false;
	    cout << "\n\nDirectory " << dirname << " should not exist by now, but it still does?" << endl;
	  }
	}
	catch(...){
	}
      }
      if(!isok){
	if(++waititer>=3) waitseconds=60;
	if(waituntilok){
	  dateStamp(cout);
	  cout << "\nTo not loose eventual results which took a long time to compute, MIRA will wait until"
	    "\neither you have resolved the problem manually or you killed MIRA."
	    "\nNext directory check in " << waitseconds << " seconds.\n";
	  sleep(waitseconds);
	}else{
	  return -1;
	}
      }
    }
    if(verbose) cout << "done.\n";
  }
  return 0;
}



/*************************************************************************
 *
 * returns 0 if everything went fine
 * -1 if something unexpected happened and not waituntilok
 *
 *************************************************************************/

int createDirectory(const string & dirname, bool verbose, bool waituntilok)
{
  boost::system::error_code ec;

  bool isok=false;
  uint32 waitseconds=1;
  uint32 waititer=0;

  if(verbose) {
    cout << "Creating directory " << dirname << " ... ";
  }

  isok=false;
  waitseconds=1;
  waititer=0;

  while(!isok){
    isok=true;
    try{
      boost::filesystem::create_directory(dirname,ec);
      if(ec){
	isok=false;
	cout << "\n\nCould not create directory " << dirname << ": " << ec.message() << endl;
      }
    }
    catch(boost::filesystem::filesystem_error fse){
      isok=false;
      cout << "\n\nCould not create directory " << dirname << ", error message is: " << fse.what() << endl;
    }
    catch(...){
      isok=false;
      cout << "\n\nCould not create directory " << dirname << ", unspecified reason." << endl;
    }
    if(isok && !boost::filesystem::is_directory(dirname,ec)){
      isok=false;
      cout << "\n\nError: directory " << dirname << " should exist by now, but it does not." << endl;
    }
    if(!isok){
      if(++waititer>=3) waitseconds=60;
      if(waituntilok){
	dateStamp(cout);
	cout << "\nTo not loose eventual results which took a long time to compute, MIRA will wait until"
	  "\neither you have resolved the problem manually or you killed MIRA."
	  "\nNext directory check in " << waitseconds << " seconds.\n";
	sleep(waitseconds);
      }else{
	return -1;
      }
    }
  }
  if(verbose) cout << "done.\n";

  return 0;
}




/*************************************************************************
 *
 * returns true if file exists
 *
 *************************************************************************/

bool fileExists(const std::string& filename)
{
    struct stat buf;
    return (stat(filename.c_str(), &buf) != -1);
}



/*************************************************************************
 *
 * cannot easily use boost::filesystem::file_copy()
 * in case BOOST was not compiled with c++0x flags we will get linker errors
 * there IS a workaround with a funny define, but I do not like that:
 *  the following from
 * http://www.robertnitsch.de/notes/cpp/cpp11_boost_filesystem_undefined_reference_copy_file
 * as workaround if BOOST was not compiled with -std=c++0x
 * I'm not sure I like this
 * #define BOOST_NO_SCOPED_ENUMS
 *
 *************************************************************************/

void fileCopy(const std::string& from, const std::string & to)
{
  std::ifstream  src(from, ios::in | ios::binary);
  if(!src.is_open()){
    cout << "fileCopy: Could not open src " << from << endl;
    exit(999);
  }

  std::ofstream  dst(to,ios::out | ios::binary);
  if(!dst.is_open()){
    cout << "fileCopy: Could not open dst " << to << endl;
    exit(999);
  }

  dst << src.rdbuf();
  if(dst.bad()){
    cout << "fileCopy: Could not finish copying to " << to << endl;
    exit(999);
  }
}


/*************************************************************************
 *
 * boost::filesystem::rename will fail if files are on different filesystems
 * let's try to be clever
 *
 *************************************************************************/

void fileRename(const std::string& from, const std::string & to)
{
  try{
    boost::filesystem::rename(from,to);
  }
  catch(...){
    fileCopy(from,to);
    removeFile(from,true);
  }
}



/*************************************************************************
 *
 * returns like boost::filesystem::remove():
 *     false if filename did not exist in the first place, otherwise true.
 *
 *************************************************************************/


bool removeFile(const std::string& filename, bool waituntilok)
{
  boost::system::error_code ec;
  try{
    if(!boost::filesystem::exists(filename,ec)) return false;
  }
  catch(...){
  }

  bool retval=false;
  bool isok=false;
  uint32 waitseconds=1;
  uint32 waititer=0;
  while(!isok){
    isok=true;
    // remove should not throw according to docs when error_code is used, but play it safe:
    //  remove_all should neither, but at least it sometimes does for BOOST 1.50 (see above in purgeCreateDir())
    try{
      retval=boost::filesystem::remove(filename,ec);
    }
    catch(boost::filesystem::filesystem_error fse){
      isok=false;
      cout << "\n\nCould not delete file " << filename << ", error message is: " << fse.what() << endl;
    }
    catch(...){
      isok=false;
      cout << "\n\nCould not delete file " << filename << ", unspecified reason." << endl;
    }
    if(isok){
      try{
	if(boost::filesystem::exists(filename,ec)){
	  isok=false;
	  cout << "\n\nFile " << filename << " should not exist by now, but it still does?" << endl;
	}
      }
      catch(...){
      }
    }
    if(!isok){
      if(++waititer>=3) waitseconds=60;
      if(waituntilok){
	dateStamp(cout);
	cout << "\nTo not loose eventual results which took a long time to compute, MIRA will wait until"
	  "\neither you have resolved the problem manually or you killed MIRA."
	  "\nNext file check in " << waitseconds << " seconds.\n";
	sleep(waitseconds);
      }else{
	return retval;
      }
    }
  }
  return retval;
}



/*************************************************************************
 *
 * returns true if file already existed and we can append
 * opens the ofstream
 *
 *************************************************************************/

bool openFileForAppend(const string & filename, ofstream & fout, bool deleteanyway)
{
  struct stat st;
  if(deleteanyway || stat(filename.c_str(),&st)) {
    CEBUG("Opening " << filename << " and truncating.\n");
    fout.open(filename.c_str(), ios::out | ios::trunc);
    return false;
  }
  CEBUG("Opening " << filename << " and appending.\n");
  fout.open(filename.c_str(), ios::out | ios::app);
  return true;
}



/*************************************************************************
 *
 * Input: filename
 * Output: pathto, stem, filetype, ziptype
 *
 * example:  name            pathto    stem     filetype     ziptype
 *           "foo.bar.gz"              "foo"     "bar"        2
 *           "fo.bar.baz"              "fo.bar"  "baz"        0
 *           "ble/bla.bz2"    ble      "bla"     ""           3
 *           "/ble/bla.z"     /ble     "bla"     ""           1
 *
 *************************************************************************/

void guessFileAndZipType(const std::string & filename, std::string & pathto, std::string & stem, std::string & filetype, uint8 & ziptype)
{
  ziptype=0;
  filetype.clear();
  fs::path mypath(filename);

  // careful: after this, mypath may be different (missing .gz / .bz2)
  if(mypath.has_extension()){
    while(mypath.has_extension()
	  && (boost::iequals(mypath.extension().native(),".gz")
	      || boost::iequals(mypath.extension().native(),".bz2")
	      || boost::iequals(mypath.extension().native(),".z"))){
      if(ziptype==0){
	if(boost::iequals(mypath.extension().native(),".z")){
	  ziptype=1;
	}else if(boost::iequals(mypath.extension().native(),".gz")){
	  ziptype=2;
	}else if(boost::iequals(mypath.extension().native(),".bz2")){
	  ziptype=3;
	}
      }
      mypath=mypath.stem();
    }
    // careful with names just ending with "."
    if(mypath.extension().native().size()>1){
      filetype=mypath.extension().native().substr(1,mypath.extension().native().size());
    }
  }
  pathto=mypath.parent_path().native();
  stem=mypath.stem().native();
  // bug in BOOST: filename entries like "tmp/" give "." as stem while an empty entry is expected
  if(stem==".") stem.clear();

  //cout << "gfazt ## fn: " << filename << " pt: " << pathto << " st: " << stem << " ft: " << filetype << " zt: " << static_cast<uint16>(ziptype) << endl;
}


/*************************************************************************
 *
 * Like guessFileAndZipType(), but sets filetype to "" if type is not
 *  one of the file types MIRA reads
 *
 *************************************************************************/


void getCanonicalFileAndZipType(const std::string & filename, std::string & pathto, std::string & stem, std::string & filetype, uint8 & ziptype)
{
  guessFileAndZipType(filename,pathto,stem,filetype,ziptype);
  if(!filetype.empty()){
    if(filetype == "fastq"
       || filetype == "fasta"
       || filetype == "fa"
       || filetype == "fna"
       || filetype == "caf"
       || filetype == "maf"
       || filetype == "exp"
       || filetype == "gbf"
       || filetype == "gbff"
       || filetype == "gbk"
       || filetype == "gb"
       || filetype == "gff"
       || filetype == "gff3"
       || filetype == "xml"
       || filetype == "ssaha2"
       || filetype == "smalt"
       || filetype == "fofnexp"){
      // TODO: how to do fastanoqual??? probably: .fna
      if(filetype=="gbff" || filetype=="gbk" || filetype=="gb"){
	filetype="gbf";
      }else if(filetype=="gff"){
	filetype="gff3";
      }
    }else{
      filetype.clear();
    }
  }
}



/*************************************************************************
 *
 * for a given path which may contain directories and/or glob-characters
 *  with special meaning like ? or * (e.g. "something/some?re/file.*"),
 *  returns a list of fnft_t (file name/ file type) for every single file
 *  which matches
 *
 * expansion: to allow MIRA to load data which does not have a file-type
 *  postfix (needed e.g. for Galaxy), also parse special double-colon
 *  format borrowed from EMBOSS:   format::path
 * E.g.: "fastq.gz::something/some?re/*.dat"
 * This allows complete misdirections ("fastq::directory/somefile.fasta"),
 *  but then this really is a problem of the user.
 *
 * return: true if error encountered
 *
 * appends(!) list of filenames & deduced types to "fnftl"
 *
 * would be better to implement this as function which uses a callback
 *  for every found file, but what the heck: Sanger is on the way out,
 *  use cases with millions of files will be very rare
 *
 *************************************************************************/

#define CEBUG(bla)   {cout << bla; cout.flush(); }
bool globWalkPath(string dn, list<fnft_t> & fnftl)
{
  FUNCSTART("bool globWalkPath(const string & dn, list<fnft_t>& fnftl)");

  bool hassomeerror=false;

  // parse the double-colon part if present
  string overridetype;

  {
    string tmptype;
    auto dnI=dn.cbegin();
    for(; dnI!=dn.cend(); ++dnI){
      if(*dnI==':' && dnI+1!=dn.cend() && *(dnI+1)==':'){
	auto tlen=dnI-dn.cbegin();
	// to be able to use guessFileAndZipType(), we need to prepend a dot to
	//  the string to be parsed for file types
	tmptype=string(".")+dn.substr(0,tlen);
	dn=dn.substr(tlen+2,dn.size()-(tlen+2));
	//cout << "### " << dn << " ###" << tmptype << "###\n";
	break;
      }
    }
    if(!tmptype.empty()){
      uint8 ziptype;
      string dummypathto;
      string dummystem;
      getCanonicalFileAndZipType(tmptype,dummypathto,dummystem,overridetype,ziptype);
      if(overridetype.empty()){
	tmptype=tmptype.substr(1,999);
	MIRANOTIFY(Notify::FATAL,"You tried to explicitely define a file type via " << tmptype << "::" << dn << ", however the type '" << tmptype << "' is not known to MIRA.");
      }
      //cout << "###" << overridetype << "###\n";
    }
  }

  // unfortunately, boost::filesystem has no glob and boost no glob-like regex
  glob_t globpaths;
  globpaths.gl_offs = 0;
  int globret=glob(dn.c_str(), 0, nullptr, &globpaths);
  //paths.gl_pathv[0] = "ls";
  if(globret!=0){
    cout << "Data '" << dn << "' was not found (neither as file nor as directory) or led to a read error.\n";
    hassomeerror=true;
  }

  if(!hassomeerror){
    for(auto thisgpath=globpaths.gl_pathv; *thisgpath != nullptr; ++thisgpath){
      fs::path path(*thisgpath);
      if(!fs::exists(path)){
	cout << "Problem with this data (file or path): " << *thisgpath << " was not found?" << endl;
	hassomeerror=true;
	continue;
      }
      list<fs::path> filepaths;
      if(fs::is_directory(path)){
	for (fs::directory_iterator dir_itr(path),end_iter; dir_itr !=  end_iter; ++dir_itr) {
	  fs::path subpath(dir_itr->path());
	  if(fs::is_regular_file(subpath) || fs::is_symlink(subpath)){
	    filepaths.push_back(subpath);
	  }
	}
      }else if(fs::is_regular_file(path) || fs::is_symlink(path)){
	filepaths.push_back(path);
      }

      uint8 ziptype=0;
      string filetype;
      string dummypathto;
      string dummystem;
      fnft_t tmpfnft;
      for(auto mypath : filepaths){
	if(!fs::is_regular_file(mypath)) {
	  cout << "mypath is not recognised as regular file???\n";
	  hassomeerror=true;
	  continue;
	}
	CEBUG("Pushing back filename: " << mypath << endl);
	tmpfnft.fn=mypath.native();

	if(overridetype.empty()){
	  getCanonicalFileAndZipType(tmpfnft.fn,dummypathto,dummystem,filetype,ziptype);
	  tmpfnft.ft=filetype;
	}else{
	  tmpfnft.ft=overridetype;
	}
	fnftl.push_back(tmpfnft);

	if(tmpfnft.ft.empty()){
	  cout << "Could not determine filetype from data '" << dn << "' (" << mypath << ")\n";
	  hassomeerror=true;
	}else if(ziptype>0 && filetype!="fastq"){
	  cout << "zipped data only possible for FASTQ at the moment, " << filetype << " currently not supported, sorry\n";
	  hassomeerror=true;
	}
      }
    }
  }

  globfree(&globpaths);

  FUNCEND();
  return hassomeerror;
}
#define CEBUG(bla)



/*************************************************************************
 *
 * returns: 0 == not on nfs, 1 == could not check, 2 == NFS
 *
 * appends(!) list of filenames & deduced types
 *
 *************************************************************************/

int checkForNFSMountOnDirectory(const std::string & dir, bool verbose)
{
  int retval=0;
  string runresult;

#ifdef __APPLE__
  string cmd("df -T nfs " + dir);
  // on MacOSX, the above command is expected to have an empty output if the
  //  filesystem is not on NFS.
  // It actually works, however in those cases the return value of "df" is also
  //  "1" (something failed) and not "0" (all OK).
  //
  // Sooooo ... to check runability, we make things a bit more complicated: check
  //  two runs of "df", once with nfs and once with nonfs. If both results are
  //  empty, something went wrong

  string cmd2("df -T nonfs " + dir);
  string runresult2;

  getSTDOUTFromCommand(cmd,runresult);
  getSTDOUTFromCommand(cmd2,runresult2);

  bool checkrun=!runresult.empty() || !runresult2.empty();

#else
  string cmd("stat -f -L -c %T " + dir);
  bool checkrun=checkRunabilityOfCommand(cmd);
#endif

  if(checkrun){
#ifdef __APPLE__
    if(!runresult.empty()) retval=2;
#else
    if(getSTDOUTFromCommand(cmd,runresult)){
      boost::to_lower(runresult);
      if(boost::find_first(runresult,"nfs")) retval=2;
    }
#endif
  }else{
    retval=1;
    if(verbose){
      cout << "Could not perform NFS check for directory "
	   <<  dir;
#ifdef __APPLE__
      cout << "\n\nFor a check to run smoothly, please make sure the Unix 'df' command is available"
	"\nand understands the following call: " << cmd << "\n\n";
#else
      cout << "\n\nFor a check to run smoothly, please make sure the Unix 'stat' command is available"
	"\nand understands the following call: " << cmd << "\n\n";
#endif
    }
  }

  return retval;
}
