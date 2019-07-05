// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14/Point_set_processing_3/include/CGAL/IO/read_off_points.h $
// $Id: read_off_points.h 2f9408f %aI Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Pierre Alliez and Laurent Saboret

//TODO-elo-fix-licence

#ifndef CGAL_READ_OFF_POINTS_H
#define CGAL_READ_OFF_POINTS_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/IO/io.h>
#include <CGAL/property_map.h>
#include <CGAL/value_type_traits.h>
#include <CGAL/Origin.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Kernel_traits.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <iostream>
#include <sstream>
#include <string>

namespace CGAL {


/**
   \ingroup PkgPointSetProcessing3IO
   Reads points (positions + normals, if available) from a .off ASCII stream.
   The function expects for each point a line with the x y z position,
   optionally followed by the nx ny nz normal.
   Faces are ignored.

   \tparam OutputIteratorValueType type of objects that can be put in `OutputIterator`.
   It is default to `value_type_traits<OutputIterator>::%type` and can be omitted when the default is fine.
   \tparam OutputIterator iterator over output points.

   \param stream input stream.
   \param output output iterator over points.
   \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

   \cgalNamedParamsBegin
     \cgalParamBegin{point_map} a model of `WritablePropertyMap` with value type `geom_traits::Point_3`.
     If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
     \cgalParamBegin{normal_map} a model of `ReadWritePropertyMap` with value type
     `geom_traits::Vector_3`. If this parameter is omitted, normals in the input stream are
     ignored.\cgalParamEnd
     \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
   \cgalNamedParamsEnd

   \return true on success.
*/
template <typename OutputIteratorValueType,
          typename OutputIterator,
#ifdef DOXYGEN_RUNNING
          typename NamedParameters
#else
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS
#endif
>
bool
read_off_points(
    std::istream& stream,
    OutputIterator output,
    bool &normals_found, // MEPP2 PATCH, normal indicator for caller
#ifdef DOXYGEN_RUNNING
    const NamedParameters& np)
#else
    const CGAL_BGL_NP_CLASS& np)
#endif
{
  using boost::choose_param;

  typedef Point_set_processing_3::Fake_point_range<OutputIteratorValueType> PointRange;
  
  // basic geometric types
  typedef typename Point_set_processing_3::GetPointMap<PointRange, CGAL_BGL_NP_CLASS>::type PointMap;
  typedef typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::type NormalMap;
  typedef typename Point_set_processing_3::GetK<PointRange, CGAL_BGL_NP_CLASS>::Kernel Kernel;

  bool has_normals = !(boost::is_same<NormalMap,
                       typename Point_set_processing_3::GetNormalMap<PointRange, CGAL_BGL_NP_CLASS>::NoMap>::value);

  PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
  NormalMap normal_map = choose_param(get_param(np, internal_np::normal_map), NormalMap());
  
  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  // typedef typename value_type_traits<OutputIterator>::type Enriched_point;
  typedef OutputIteratorValueType Enriched_point;

  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;

  if(!stream)
  {
    std::cerr << "Error: cannot open file" << std::endl;
    return false;
  }

  // scan points
  long pointsCount = 0, facesCount = 0, edgesCount = 0; // number of items in file
  int pointsRead = 0; // current number of points read
  int lineNumber = 0; // current line number
  std::string line;
  std::istringstream iss;
  normals_found = false; // MEPP2 PATCH

  while(getline(stream,line))
  {
    iss.clear();
    iss.str(line);

    // Ignore empty lines and comments
    if (line.empty () || line[0] == '#')
      continue;

    lineNumber++;

    // Reads file signature on first line
    if (lineNumber == 1)
    {
      std::string signature;
      if ( !(iss >> signature)
        || (signature != "OFF" && signature != "NOFF") )
      {
        // if wrong file format
        std::cerr << "Incorrect file format line " << lineNumber << " of file" << std::endl;
        return false;
      }
    }

    // Reads number of points on 2nd line
    else if (lineNumber == 2)
    {
      if ( !(iss >> pointsCount >> facesCount >> edgesCount) )
      {
        std::cerr << "Error line " << lineNumber << " of file" << std::endl;
        return false;
      }
    }

    // Reads 3D points on next lines
    else if (pointsRead < pointsCount)
    {
      // Reads position + normal...
      double x,y,z;
      double nx,ny,nz;
      if (iss >> iformat(x) >> iformat(y) >> iformat(z))
      {
        Point point(x,y,z);
        Vector normal = CGAL::NULL_VECTOR;
        // ... + normal...
        if (iss >> iformat(nx))
        {
          // In case we could read one number, we expect that there are two more
          if(iss  >> iformat(ny) >> iformat(nz)){
            normal = Vector(nx,ny,nz);
            normals_found = true; // MEPP2 PATCH
          } else {
            std::cerr << "Error line " << lineNumber << " of file" << std::endl;
            return false;
          }
        }
#if 0 // MEPP2 PATCH
        Enriched_point pwn;
        put(point_map,  pwn, point);  // point_map[&pwn] = point
        if (has_normals)
          put(normal_map, pwn, normal); // normal_map[&pwn] = normal
        *output++ = pwn;
#else // MEPP2 PATCH
        *output++ = point;

        if (has_normals)
          put(normal_map, pointsRead, normal); // normal_map[index] = normal
#endif // MEPP2 PATCH

        pointsRead++;
      }
      // ...or skip comment line
    }
    // Skip remaining lines
  }

  return true;
}

// variant with default output iterator value type
template <typename OutputIterator,
          typename CGAL_BGL_NP_TEMPLATE_PARAMETERS>
bool
read_off_points(
  std::istream& stream, ///< input stream.
  OutputIterator output,
  bool &normals_found, // MEPP2 PATCH, normal indicator for caller
  const CGAL_BGL_NP_CLASS& np)
{
  return read_off_points<typename value_type_traits<OutputIterator>::type>
    (stream, output, normals_found, np);
}


} //namespace CGAL

#endif // CGAL_READ_OFF_POINTS_H
