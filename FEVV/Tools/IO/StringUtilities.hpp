// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif


namespace FEVV {
namespace StrUtils {

/**
 * Split a string according to the provided delimiters.
 */
inline
std::vector< std::string >
split(const std::string &str,
      const std::string &delims,
      bool keep_empty_tokens = false)
{
  std::vector< std::string > tokens, non_empty_tokens;

  boost::split(
      tokens,
      str,
      boost::is_any_of(delims)); // boost string split => gives a warning: C4996
                                 // std::copy::unchecked_iterators::_Deprecate

  if(keep_empty_tokens)
    return tokens;
  // TODO : these two codes do the same thing. The second one is faster but
  // leave compilation warning on Windows
  for(std::vector< std::string >::iterator it = tokens.begin();
      it != tokens.end();
      ++it)
    if(!it->empty())
      non_empty_tokens.push_back(*it);
  /*
          std::copy_if(	tokens.begin(),
                                          tokens.end(),
                                          std::back_inserter(non_empty_tokens),
                                          [](std::string s) {return !s.empty();
     }); //vector copy to avoid empty tokens
  */
  return non_empty_tokens;
}

/**
 * Returns true if the two strings are identical.
 */
inline
bool
is_equal(const std::string &str1, const std::string &str2)
{
  return str1.compare(str2) == 0;
}

/**
 * Convert a string into another type.
 */
template< typename ConvertType >
void
convert(const std::string &str, ConvertType &elem)
{
  std::stringstream ss(str);
  ss >> elem;
}

/**
 * Convert a string into another type.
 */
template< typename ConvertType >
ConvertType
convert(const std::string &str)
{
  ConvertType elem;
  convert(str, elem);
  return elem;
}

/**
 * Convert a data into a string.
 */
template< typename ScalarType >
void
convert(const ScalarType &s, std::string &st)
{
  std::ostringstream conv;
  conv << s;
  st = conv.str();
}

/**
 * Convert a data into a string.
 */
template< typename ScalarType >
std::string
convert(const ScalarType &s)
{
  std::string st;
  convert< ScalarType >(s, st);
  return st;
}

} // namespace StrUtils
} // namespace FEVV

