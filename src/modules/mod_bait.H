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


#include <fstream>

#include <string>
#include <vector>
#include <list>

#include "mira/parameters.H"
#include "mira/contig.H"
#include "mira/readpool.H"
#include "mira/hashstats.H"


class MiraBait
{
private:


  static std::vector<MIRAParameters> MB_Pv;

  static std::string MB_fromtype;
  static std::list<std::string> MB_totype;

  static std::list<std::ofstream *> MB_ofs;


  static std::string MB_baitfile;
  static std::string MB_infile;
  static std::string MB_outbasename;

  static bool MB_loadhashstat;
  static bool MB_deletestaronlycolumns;
  static bool MB_inversehit;
  static bool MB_fwdandrev;
  static uint32 MB_numbaithits;

  static bool MB_mustdeletetargetfiles;

  static std::list<Contig> MB_clist;   // needed for CAF conversion (and GBF)

  static HashStatistics MB_hashstatistics;

  static uint64 MB_numreadsread;
  static uint64 MB_numreadswritten;

private:
  static void usage();
  static void checkTypes(std::string & fromtype,std::list<std::string> & totype);
  static void putReadsInContigsAndSave(std::vector<MIRAParameters> & Pv, ReadPool & rp);
  static void specialTestCode(std::list<Contig> & clist, ReadPool & rp);

  static void saveReadPool(ReadPool & rp);
  static void cafmafload_callback(std::list<Contig> & clist, ReadPool & rp);
  static void readpoolload_callback(ReadPool & rp);

public:
  ~MiraBait();

  int mainMiraBait(int argc, char ** argv);

};

