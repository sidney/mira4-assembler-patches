/*
 * Written by Bastien Chevreux (BaCh)
 *
 * Copyright (C) 2011 and later by Bastien Chevreux
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


#ifndef _util_stlimprove_h
#define _util_stlimprove_h

#include "boost/concept/assert.hpp"

template <typename T>
void zerovector(std::vector<T> & vec)
{
  BOOST_CONCEPT_ASSERT((boost::Integer<T>));
  if(vec.size()){
    memset(&vec[0],0,sizeof(T)*vec.size());
  }
}

#endif
