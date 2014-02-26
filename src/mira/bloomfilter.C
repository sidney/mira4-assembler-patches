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

#include "bloomfilter.H"
#include "util/dptools.H"

using namespace std;


// Plain vanilla constructor
BloomFilter::BloomFilter(const uint8 bits, const uint32 numkeys)
{
  FUNCSTART("BloomFilter::BloomFilter()");

  cout << "initialising Bloom filter with " << static_cast<uint16>(bits) << " bits (" << (0x1ull<<bits) << " elements), " << (0x1ull<<bits)/4 << " bytes\n";

  BUGIFTHROW(bits>64,"bits >64 ?");
  BUGIFTHROW(numkeys==0||numkeys>20,"numkeys == " << numkeys << " ???");
  BF_numkeys=numkeys;
  BF_bfaddressmask=(0x1uLL<<bits)-1;
  BF_bloomfield.resize((0x1ull<<bits)/4);
  reset();

  FUNCEND();
}


BloomFilter::~BloomFilter()
{
  FUNCSTART("BloomFilter::~BloomFilter()");

  discard();

  FUNCEND();
}


void BloomFilter::discard()
{
  FUNCSTART("BloomFilter::discard()");

  BF_bloomfield.clear();
  BF_level1count=0;
  BF_level2count=0;
  BF_numuniqkmers=0;
  BF_numkmerseenge2=0;
  BF_numkmerseenge3=0;

  FUNCEND();
}

void BloomFilter::reset()
{
  discard();
  BF_bloomfield.resize(BF_bloomfield.capacity(),0);
}


ostream & operator<<(ostream &ostr, BloomFilter const &bf)
{
  FUNCSTART("friend ostream & BloomFilter::operator<<(ostream &ostr, const  &bf)");

  ostr << "BloomFilter:"
       << "\nmemory: " << bf.BF_bloomfield.size()
       << "\nkeys per kmer: " << bf.BF_numkeys
       << "\nl1 occupancy: " << bf.BF_level1count
       << "\nl2 occupancy: " << bf.BF_level2count
       << "\nprobable num kmers    : " << bf.BF_numuniqkmers
       << "\n thereof num kmers>=2 : " << bf.BF_numkmerseenge2
       << "\ntotal kmers counted>=3: " << bf.BF_numkmerseenge3
       << endl;

  FUNCEND();
  return ostr;
}


