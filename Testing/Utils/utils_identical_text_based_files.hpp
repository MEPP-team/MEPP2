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
#include <vector>
#include <fstream>
#include <iostream>

/*
 * \brief  Compare two text based files in search possible differences.
 * \param  filenameA    Filename of first file.
 * \param  filenameB    Filename of second file.
 * \param  skip         skip lines beginning with the provided strings.
 * \result true when files are identical false otherwise
 */
inline
bool
identical_text_based_files(std::string filename_a,
                           std::string filename_b,
                           const std::vector< std::string > &skip =
                              std::vector< std::string >())
{
  std::ifstream file_a(filename_a);
  if(!file_a)
  {
    std::cout << "Unable to read first file." << filename_a << std::endl;
    return false;
  }

  std::ifstream file_b(filename_b);
  if(!file_b)
  {
    std::cout << "Unable to read second file " << filename_b << std::endl;
    return false;
  }

  bool result = true;
  while((!file_a.eof()) && (!file_b.eof()))
  {
    std::string line_a, line_b;

    // read one line from each file
    getline(file_a, line_a);
    getline(file_b, line_b);

    // skip lines if requested
    bool skip_line = false;
    for(auto &skip_str : skip)
    {
      if((line_a.substr(0, skip_str.size()) == skip_str) &&
         (line_b.substr(0, skip_str.size()) == skip_str))
      {
        skip_line = true;
        break;
      }
    }
    if(skip_line)
    {
      std::cout << "skip line_a '" << line_a << "'" << std::endl;
      std::cout << "skip line_b '" << line_b << "'" << std::endl;
      continue; // stop processing the current line
    }

    // remove DOS end of line extra character if any
    if((!line_a.empty()) && (line_a.back() == '\r'))
      line_a.pop_back();
    if((!line_b.empty()) && (line_b.back() == '\r'))
      line_b.pop_back();

    // compare the 2 lines
    if(line_a != line_b)
    {
      result = false;
      std::cout << "Comparison with reference file failed" << std::endl;
      std::cout << "<line_a '" << line_a << "'" << std::endl;
      std::cout << "---" << std::endl;
      std::cout << ">line_b '" << line_b << "'" << std::endl;
      break;
    }
  }

  file_a.close();
  file_b.close();

  return result;
}
