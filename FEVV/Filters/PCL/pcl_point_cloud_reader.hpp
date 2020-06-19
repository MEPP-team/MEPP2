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

#include "FEVV/Tools/IO/FileUtilities.hpp" // for FileUtils::has_extension()


#if defined _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4267) // for VS-2015
#endif

#include <pcl/io/ascii_io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

#if defined _MSC_VER
#pragma warning(pop)
#endif


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

// parse XYZ 1st line to detect normals
inline
void
parse_xyz_1stline(const std::string &filename,
                  bool &has_normal)
{
  // open input file
  std::ifstream in(filename);
  if(! in)
  {
    throw std::invalid_argument("read_mesh() error: can not open file " +
                                filename);
  }

  std::string line, word;
  has_normal = false;

  // load 1st line
  std::getline(in, line);
  std::istringstream line_ss(line);

  // count number of values on line
  unsigned int val_nbr = 0;
  while(line_ss >> word)
  {
    val_nbr++;
  }

  // 3 values -> geometry only
  // 6 values -> geometry + normals
  // other    -> unsupported
  if(val_nbr == 3)
    has_normal = false;
  else if(val_nbr == 6)
    has_normal = true;
  else
  {
    throw std::invalid_argument("read_mesh() error: unsupported XYZ format in file " +
                                filename);
  }
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
inline
void
read_mesh(
    const std::string &filename,
    FEVV::PCLPointCloud &g,
    PMapsContainer &pmaps,
    bool only_pts = false)
{
#if 1
  bool success = false;
  bool has_normal = false;
  bool has_color = false;

  if(FEVV::FileUtils::has_extension(filename, ".xyz"))
  {
    // parse header to detect normals
    parse_xyz_1stline(filename, has_normal);

    // configure PCL reader
    pcl::ASCIIReader reader;
    reader.setExtension(".xyz");
    if(has_normal)
    {
      // note:
      // - 'reader.setInputFields< pcl::PointNormal >()' doesn't
      //   work because pcl::PointNormal type has a 'curvature' field, 
      //   and all field must be present in the file else the reading
      //   fails ; no other predefined PCL Point type matches our needs
      // - the solution seems to use the function
      //   'setInputFields(const std::vector<pcl::PCLPointField>& fields)',
      //   but how to build the fields vector with the right values ?
      // - answer: call 'pcl::getFields< pcl::PointNormal >(fields)'
      //   to get the fields vector for pcl::PointNormal, then remove
      //   the last unwanted field (curvature)
      // - pcl::getFields< pcl::PointNormal >(fields);
      //   -> fields = std::vector of length 7, capacity 8 =
      //      {
      //        {name = "x", offset = 0, datatype = 7 '\a', count = 1},
      //        {name = "y", offse t = 4, datatype = 7 '\a', count = 1},
      //        {name = "z", offset = 8, datatype = 7 '\a', count = 1},
      //        {name = "normal_x", offs et = 16, datatype = 7 '\a', count = 1},
      //        {name = "normal_y", offset = 20, datatype = 7 '\a', count = 1},
      //        {name = "normal_z", offset = 24, datatype = 7 '\a', count = 1},
      //        {name = "curvature", offset = 32, datatype = 7 '\a', count = 1}
      //      }
      //

      // retrieve all pcl::PointNormal fields
      std::vector<pcl::PCLPointField> fields;
      pcl::getFields< pcl::PointNormal >(fields);
      // remove the last unwanted field (curvature)
      fields.pop_back();
      reader.setInputFields(fields);
    }
    else
      reader.setInputFields< pcl::PointXYZ >();

    // read XYZ file
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
