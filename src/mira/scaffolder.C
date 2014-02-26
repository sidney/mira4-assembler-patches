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


#include "mira/scaffolder.H"


using namespace std;


// Plain vanilla constructor
Scaffolder::Scaffolder()
{
  FUNCSTART("Scaffolder::Scaffolder()");

  zeroVars();
  init();

  ERROR("Not implemented yet.");
  FUNCEND();
}

void Scaffolder::zeroVars()
{
  FUNCSTART("void Scaffolder::zeroVars()");
  FUNCEND();
}

void Scaffolder::init()
{
  FUNCSTART("void Scaffolder::init()");
  FUNCEND();
}



Scaffolder::~Scaffolder()
{
  FUNCSTART("Scaffolder::~Scaffolder()");

  discard();

  ERROR("Not implemented yet.");
  FUNCEND();
}


void Scaffolder::discard()
{
  FUNCSTART("Scaffolder::discard()");

  zeroVars();

  FUNCEND();
}


//// Copy constructor
////  no discard needed as this object will be freshly created when
////  called through this constructor
//Scaffolder::Scaffolder(Scaffolder const &other)
//{
//  FUNCSTART("Scaffolder::Scaffolder(Scaffolder const &other)");
//
//  ??_valid=0;
//
//  *this=other;                               // call the copy operator
//
//  FUNCEND();
//}
//
//// Copy operator, needed by copy-constructor
//Scaffolder const & Scaffolder::operator=(Scaffolder const & other)
//{
//  FUNCSTART("Scaffolder const & Scaffolder::operator=(Scaffolder const & other)");
//  ERROR("Not implemented yet.");
//  FUNCEND();
//  return *this;
//}

//ostream & operator<<(ostream &ostr, Scaffolder const &???)
//{
//  FUNCSTART("friend ostream & Scaffolder::operator<<(ostream &ostr, const  &???)");
//  ERROR("Not implemented yet.");
//
//  FUNCEND();
//  return ostr;
//}


void Scaffolder::storeInfoFreshContig(Contig & con)
{
  FUNCSTART("void Scaffolder::storeInfoFreshContig(contig & con)");

  {
    scaffinfo_t si;
    SCA_scaffinfolist.push_back(si);
  }
  scaffinfo_t & si=SCA_scaffinfolist.back();
  si.contigname=con.getContigName();

  FUNCEND();
}
