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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <vector>
#include <algorithm>
#include <iterator>
#include <cctype> // for tolower()

#include <boost/filesystem.hpp>


namespace FEVV {
namespace FileUtils {


/**
 * Returns the file name extension or an empty string if the file name has
 * no extension.
 */
inline
std::string
get_file_extension(const std::string &file_name)
{
  return static_cast< std::string >(boost::filesystem::extension(file_name));
}

/**
 * Returns the file name extension or an empty string if the file name has
 * no extension.
 */
inline
std::string
get_file_extension(const char *file_name)
{
  return get_file_extension(std::string(file_name));
}

/**
 * Returns true if the file name has an extension,
 * aka ends up with '.something'.
 */
inline
bool
has_extension(const std::string &file_name)
{
  return !boost::filesystem::extension(file_name).empty();
}

/**
 * Returns true if the file name has an extension,
 * aka ends up with '.something'.
 */
inline
bool
has_extension(const char *file_name)
{
  return has_extension(std::string(file_name));
}

/**
 * Returns true if the file name ends up with the provided
 * extension.
 */
inline bool
has_extension(const std::string &file_name,
              const std::string &ext_name,
              bool case_sensitive = false)
{
  std::string file_ext = get_file_extension(file_name);

  if(! case_sensitive)
  {
    std::transform(
        file_ext.begin(), file_ext.end(), file_ext.begin(), ::tolower);
  }

  return ext_name.compare(file_ext) == 0;
}

/**
 * Returns true if the file name ends up with the provided
 * extension.
 */
inline bool
has_extension(const char *file_name,
              const char *ext_name,
              bool case_sensitive = false)
{
  return has_extension(
      std::string(file_name), std::string(ext_name), case_sensitive);
}

/**
 * Returns true if the file name ends up with one of the provided
 * extensions.
 */
inline bool
has_extension(const std::string &file_name,
              const std::vector< std::string > &ext_names,
              bool case_sensitive = false)
{
  return std::any_of(
      ext_names.cbegin(), ext_names.cend(), [&](const std::string &ext) {
        return has_extension(file_name, ext, case_sensitive);
      });
}

/**
 * Returns the file name with no path and no extension.
 * Example: for '/path/foo.bar', returns 'foo'.
 */
inline
std::string
get_file_name(const std::string &file_name)
{
  boost::filesystem::path p(file_name);
  boost::filesystem::path file_full_name = p.filename();
  return static_cast< std::string >(file_full_name.stem().string());
}

/**
 * Returns the file name with no path and no extension.
 * Example: for '/path/foo.bar', returns 'foo'.
 */
inline
std::string
get_file_name(const char *file_name)
{
  return get_file_name(std::string(file_name));
}

/**
 * Returns the file name with no path.
 * Example: for '/path/foo.bar', returns 'foo.bar'.
 */
inline
std::string
get_file_full_name(const std::string &file_name)
{
  boost::filesystem::path p(file_name);
  return static_cast< std::string >(p.filename().string());
}

/**
 * Returns the number of lines in the file.
 */
inline
unsigned int
count_file_lines(const std::string &file_name)
{
  std::ifstream myfile(file_name);
  myfile.unsetf(std::ios_base::skipws);

  unsigned int line_count = static_cast< unsigned int >(
      std::count(std::istreambuf_iterator< char >(myfile),
                 std::istreambuf_iterator< char >(),
                 '\n'));
  return line_count + 1;
}

/**
 * Returns the number of lines in the file.
 */
inline
unsigned int
count_file_lines(const char *file_name)
{
  return count_file_lines(std::string(file_name));
}

/**
 * Returns the file path.
 * Example: for '/dir1/dir2/foo.bar', returns '/dir1/dir2'.
 */
inline
std::string
get_parent_directory(const std::string &file_name)
{
  return static_cast< std::string >(
      boost::filesystem::path(file_name).branch_path().string());
}

/**
 * Get a line from the file and fix DOS end of line extra character if any.
 * To be used in the same way as std::getline(stream, string).
 *
 * \param  input  the file-stream to read from
 * \param  str    the string that will be filled with the line being read
 *
 * \return  the file-stream (mimic the return type of std::getline)
 */
inline
std::ifstream &
safe_getline(std::ifstream &input, std::string &str)
{
  std::getline(input, str);

  // remove DOS end of line extra character
  if((!str.empty()) && (str.back() == '\r'))
    str.pop_back();

  return input;
}

/**
 * Get a line from the input stream skipping any commented line.
 * To be used in the same way as std::getline(stream, string).
 *
 * \param  input  the input stream to read from
 * \param  line   the string that will be filled with the line being read
 *
 * \return  true if the read succeeded, else false (aka when EOF)
 */
inline
bool
getline_skip_comment(std::istream &input, std::string &line)
{
  while(std::getline(input, line))
  {
    // skip empty line and comment line
    if(! (line.empty() || line[0] == '#'))
      break;
  }

  return static_cast< bool >(input);
}

/**
 * Get a line from the input stream skipping any commented line.
 * To be used in the same way as std::getline(stream, string).
 *
 * \param  input    the input stream to read from
 * \param  line     the string that will be filled with the line being read
 * \param  line_ss  the stringstream based on the read line
 *
 * \return  true if the read succeeded, else false (aka when EOF)
 */
inline
bool
getline_skip_comment(std::istream &input,
                     std::string &line,
                     std::istringstream &line_ss)
{
  getline_skip_comment(input, line);
  line_ss.clear();
  line_ss.str(line);

  return static_cast< bool >(input);
}


/**
 * Copy a file. Display a warning if destination file exists
 * and do not overwrite it.
 *
 * \param  from  the name of the file to be copied
 * \param  to    the name of the target file
 */
inline
void
copy_file(const std::string &from, const std::string &to)
{
  try
  {
    boost::filesystem::copy_file(from, to);
  }
  catch(const boost::filesystem::filesystem_error &e)
  {
    std::cout << "Copy_file WARNING: " << e.what() << std::endl;
  }
}


} // namespace FileUtils
} // namespace FEVV

