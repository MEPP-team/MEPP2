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

#define PCL_NO_PRECOMPILE
#include <pcl/pcl_macros.h>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

#include <stdexcept> // for std::runtime_error
#include <ostream>


namespace FEVV {

using PCLKernelType = float;
using PCLColorType = uint8_t;
using PCLEnrichedPoint = pcl::PointXYZRGBNormal;
using PCLPointCloud = pcl::PointCloud< PCLEnrichedPoint >;


// Point, Normal, Color types with constructors are needed to write
// 'put(pm, vd, Point(1, 2, 3))'

struct PCLPoint
{
  // constructors
  PCLPoint() { }

  PCLPoint(PCLKernelType _x, PCLKernelType _y, PCLKernelType _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }

  // convenient operator[] to mimic other data structures
  PCLKernelType operator[](int i)	const
  {
    if(i == 0)
      return x;
    else if(i == 1)
      return y;
    else if(i == 2)
      return z;
    else
    {
      throw std::runtime_error("FEVV::PCLPoint: invalid index in operator[]");
      return 0; // dummy return to avoid a warning
    }
  }


  bool operator==(const FEVV::PCLPoint& rhs)
  {
    return (x == rhs.x  &&  y == rhs.y  &&  z == rhs.z);
  }

  friend std::ostream& operator<<(std::ostream& stream, const PCLPoint& p)
  {
    stream << p.x << " " << p.y << " " << p.z;
    return stream;
  }

  PCLKernelType x;
  PCLKernelType y;
  PCLKernelType z;
};

struct PCLNormal
{
  // constructors
  PCLNormal() { }

  PCLNormal(PCLKernelType _x, PCLKernelType _y, PCLKernelType _z)
  {
    x = _x;
    y = _y;
    z = _z;
  }

  // normal must behave as a vector
  PCLKernelType operator[](int i)	const
  {
    if(i == 0)
      return x;
    else if(i == 1)
      return y;
    else if(i == 2)
      return z;
    else
    {
      throw std::runtime_error("FEVV::PCLNormal: invalid index in operator[]");
      return 0; // dummy return to avoid a warning
    }
  }

  bool operator==(const FEVV::PCLNormal& rhs)
  {
    return (x == rhs.x  &&  y == rhs.y  &&  z == rhs.z);
  }

  friend std::ostream& operator<<(std::ostream& stream, const PCLNormal& n)
  {
    stream << n.x << " " << n.y << " " << n.z;
    return stream;
  }

  PCLKernelType x;
  PCLKernelType y;
  PCLKernelType z;
};

struct PCLColor
{
  // constructors
  PCLColor() { }

  PCLColor(PCLColorType _r, PCLColorType _g, PCLColorType _b)
  {
    r = _r;
    g = _g;
    b = _b;
  }

  // color must behave as a vector
  PCLColorType operator[](int i)	const
  {
    if(i == 0)
      return r;
    else if(i == 1)
      return g;
    else if(i == 2)
      return b;
    else
    {
      throw std::runtime_error("FEVV::PCLColor: invalid index in operator[]");
      return 0; // dummy return to avoid a warning
    }
  }

  bool operator==(const FEVV::PCLColor& rhs)
  {
    return (r == rhs.r  &&  g == rhs.g  &&  b == rhs.b);
  }

  friend std::ostream& operator<<(std::ostream& stream, const PCLColor& c)
  {
    stream << (int) c.r << " " << (int) c.g << " " << (int) c.b;
    return stream;
  }

  PCLColorType r;
  PCLColorType g;
  PCLColorType b;
};

} // namespace FEVV

