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

#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#include "FEVV/Wrappings/Wrappings_pcl_point_cloud.h"

#include "FEVV/Wrappings/properties.h"

#include "FEVV/Filters/Generic/generic_reader.hpp" // specialization of
#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()

#include <pcl/io/ascii_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

#include <string>
#include <stdexcept> // for std::invalid_argument


namespace FEVV {

// parse PLY header to detect normals and colors
inline
void
parse_ply_header(const std::string &filename,
                 bool &has_normal,
                 bool &has_color)
{
  // open input file
  std::ifstream in(filename);
  if(! in)
  {
    throw std::invalid_argument("read_mesh() error: can not open file " +
                                filename);
  }

  // parse file
  std::string line, word;
  has_normal = false;
  has_color = false;

  do
  {
    std::getline(in, line);
    std::istringstream line_ss(line);

    // look for "property ... nx" and "property ... r"
    line_ss >> word;
    if(word == "property")
    {
      line_ss >> word;
      line_ss >> word;

      if(word == "nx")
      {
        has_normal = true;
      }
      else if(word == "red")
      {
        has_color = true;
      }
    }
  }
  while(word != "end_header");
}

// parse PCD header to detect normals and colors
inline
void
parse_pcd_header(const std::string &filename,
                 bool &has_normal,
                 bool &has_color)
{
  // open input file
  std::ifstream in(filename);
  if(! in)
  {
    throw std::invalid_argument("read_mesh() error: can not open file " +
                                filename);
  }

  // parse file
  std::string line, word;
  has_normal = false;
  has_color = false;

  do
  {
    std::getline(in, line);
    std::istringstream line_ss(line);

    // look for "FIELDS x y z normal_x normal_y normal_z r g b"
    line_ss >> word;
    if(word == "FIELDS")
    {
      while(line_ss >> word)
      {
        if(word == "normal_x")
        {
          has_normal = true;
        }

        if(word == "rgb")
        {
          has_color = true;
        }
      }
      break;
    }
  }
  while(in);
}

} // namespace FEVV


namespace FEVV {
namespace Filters {


/**
 * \brief  Load mesh from file.
 *
 * \param  filename    name of the input mesh file
 * \param  g           mesh object to store the mesh being read
 * \param  pmaps       property maps bag to store the properties
 *                     of the mesh being read
 *
 * \sa     the variant that use the geometry traits provided by the user.
 */
template<>
inline
void
read_mesh< FEVV::PCLPointCloud, FEVV::Geometry_traits< FEVV::PCLPointCloud > >(
    const std::string &filename,
    FEVV::PCLPointCloud &g,
    PMapsContainer &pmaps,
    bool only_pts)
{
#if 1
  bool success = false;
  bool has_normal = false;
  bool has_color = false;

  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    pcl::ASCIIReader reader;
    // note: the reader wants the file to have all fields of
    //       the point datastructure used with reader.setInputFields<...>(),
    //       else it doesn't read anything ; as most often .xyz files
    //       contain only point coordinates X, Y, Z we must
    //       use 'pcl::PointXYZ' with reader.setInputFields<...>()
    //       rather than 'FEVV::PCLEnrichedPoint'
    reader.setInputFields< pcl::PointXYZ >();
    reader.setExtension(".xyz");

    success = (reader.read(filename, g) >= 0);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".pcd"))
  {
    // parse header to detect normals and colors
    parse_pcd_header(filename, has_normal, has_color);

    // load point cloud
    success = (pcl::io::loadPCDFile(filename, g) >= 0);
  }
  else if(FEVV::FileUtils::has_extension(filename, ".ply"))
  {
    // parse header to detect normals and colors
    parse_ply_header(filename, has_normal, has_color);

    // load point cloud
    success = (pcl::io::loadPLYFile(filename, g) >= 0);
  }

  if(success)
  {
    // store Normal and Color property maps into bag
    if(has_normal)
    {
      auto v_nm = make_property_map(FEVV::vertex_normal, g);
      put_property_map(FEVV::vertex_normal, g, pmaps, v_nm);
    }

    if(has_color)
    {
      auto v_cm = make_property_map(FEVV::vertex_color, g);
      put_property_map(FEVV::vertex_color, g, pmaps, v_cm);
    }
  }
  else
  {
    throw std::invalid_argument("read_mesh() error: can not read file " +
                                filename);
  }
#else
  // DBG
  g.points.push_back(FEVV::PCLPoint(-1, -1, 2));
  g.points.push_back(FEVV::PCLPoint(-1,  1, 2));
  g.points.push_back(FEVV::PCLPoint( 1,  1, 2));
  g.points.push_back(FEVV::PCLPoint( 1, -1, 2));
#endif  
}


} // namespace Filters
} // namespace FEVV
