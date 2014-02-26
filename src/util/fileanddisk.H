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


#ifndef _util_fileanddisk_h
#define _util_fileanddisk_h


// for fopen, fclose, fwrite and stat() in template functions below
#include <cstdio>
#include <sys/stat.h>

#include <iostream>

#include <string>
#include <list>
#include <vector>

#include <sys/time.h>

#include "stdinc/defines.H"

#include "stdinc/stlincludes.H"




/*************************************************************************
 *
 * own wrapper functions for buggy OSX system function fread, fwrite, ftell etc.
 *  (problems with data >= 2GiB)
 *
 *************************************************************************/

size_t myFRead(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t myFWrite(const void *ptr, size_t size, size_t nmemb, FILE *stream);
int myFSeek(FILE *stream, size_t offset, int whence);
size_t myFTell(FILE *stream);



/*************************************************************************
 *
 *
 *
 *************************************************************************/

void splitFullPathAndFileName(const std::string & what,
			      std::string & path,
			      std::string & filename
  );
void findLocationOfSelfBinary(std::string & result);

uint64 countLinesInFile(const std::string & filename);
void dumpFile(const char * fname, std::ostream & ostr);


int ensureDirectory(const std::string & dirname, bool purge, bool verbose = true, bool waituntilok=true);
int purgeCreateDir(const std::string & dirname, bool verbose, bool waituntilok);
int removeDirectory(const std::string & dirname, bool verbose, bool waituntilok);
int createDirectory(const std::string & dirname, bool verbose, bool waituntilok);
bool dirExists(const std::string& filename);
bool fileExists(const std::string& filename);
bool removeFile(const std::string& filename, bool waituntilok);
void fileCopy(const std::string& from, const std::string & to);
void fileRename(const std::string& from, const std::string & to);

bool openFileForAppend(const std::string & filename, std::ofstream & fout, bool deleteanyway);

void guessFileAndZipType(const std::string & filename,
			 std::string & pathto,
			 std::string & stem,
			 std::string & filetype,
			 uint8 & ziptype);
void getCanonicalFileAndZipType(const std::string & filename,
				std::string & pathto,
				std::string & stem,
				std::string & filetype,
				uint8 & ziptype);



struct fnft_t {
  std::string fn;      // filename
  std::string ft;      // filetype
};
bool globWalkPath(std::string dn, std::list<fnft_t> & fnftl);

int checkForNFSMountOnDirectory(const std::string & dir, bool verbose);


/*************************************************************************
 *
 * returns false if anything went wrong
 *
 *************************************************************************/

template<typename T>
bool saveVector(const std::vector<T> & vec, const std::string filename)
{
  bool ret=false;
  FILE * fout=fopen(filename.c_str(), "w");
  if(fout){
    if(vec.empty()
       || myFWrite(&vec[0], sizeof(T), vec.size(), fout) == vec.size()) {
      ret=true;
    }
    if(fclose(fout)!=0) ret=false;;
  }
  return ret;
}

/*************************************************************************
 *
 * returns false if anything went wrong
 *
 * howtohandle:
 *  - 0 means resize the vector as needed
 *  - 1: expect the vector to have exactly the needed size
 *
 *************************************************************************/

template<typename T>
bool loadVector(std::vector<T> & vec, const std::string filename, uint8 howtohandle)
{
  FUNCSTART("bool loadVector(std::vector<T> & vec, const std::string filename, uint8 howtohandle)");
  bool ret=false;

  FILE * fin=fopen(filename.c_str(), "r");
  if(fin){
    struct stat filestatus;
    stat(filename.c_str(), &filestatus);
    //std::cout << filestatus.st_size << " bytes\n";

    if(filestatus.st_size%sizeof(T) == 0){
      size_t numelements=filestatus.st_size/sizeof(T);
      if(howtohandle==0){
	vec.resize(numelements);
      } else {
	if(numelements!=vec.size()){
	  std::cout << "Could load " << vec.size() << " elements, but file has " << numelements << " ???";
	  fclose(fin);
	  return false;
	}
      }
      if(vec.empty()){
	ret=true;
      }else{
	if(myFRead(&vec[0], sizeof(T), vec.size(), fin) == vec.size()) ret=true;
      }
      if(fclose(fin)!=0) ret=false;;
    }else{
      std::cout << "File " << filename << " has " << filestatus.st_size << " bytes, that does not fit into a the vector type we want to put data into. Element size: " << sizeof(T) << " means " << filestatus.st_size%sizeof(T) << " leftover bytes.\n";
      ret=false;
    }
  }

  return ret;
}

#endif
