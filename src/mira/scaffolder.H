/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2003 and later by Bastien Chevreux
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

#ifndef _scaffolder_h_
#define _scaffolder_h_

#include <iostream>

#include "stdinc/defines.H"
#include "errorhandling/errorhandling.H"

#include "mira/contig.H"



class Scaffolder
{
  //struct
private:
  struct linkinfo_t {
    int32 readid;       // read id in readpool
    int32 distance;     // distance from border of contig

    linkinfo_t() {
      readid=-1;
      distance=0;
    }
  };

  struct scaffinfo_t {
    std::string contigname;
    std::vector<linkinfo_t> linksleft;
    std::vector<linkinfo_t> linksright;
    bool leftwelllinked;
    bool rightwelllinked;
    bool leftlinktononrep;
    bool rightlinktononrep;

    scaffinfo_t(){
      leftwelllinked=false;
      rightwelllinked=false;
      leftlinktononrep=false;
      rightlinktononrep=false;
    }
  };

  //Variables
private:
  std::list<scaffinfo_t> SCA_scaffinfolist;

public:


  //Functions
private:
  void foolCompiler();
  void init();
  void zeroVars();

public:
  Scaffolder();
  Scaffolder(Scaffolder const &other);
  ~Scaffolder();

  Scaffolder const & operator=(Scaffolder const & other);
  friend std::ostream & operator<<(std::ostream &ostr, Scaffolder const &sca);

  void discard();

  void storeInfoFreshContig(Contig & con);
};


#endif
