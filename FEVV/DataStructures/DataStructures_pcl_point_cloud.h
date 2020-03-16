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

#include <Eigen/Core>
#include <stdexcept> // for std::runtime_error
#include <ostream>


namespace FEVV {

using PCLKernelType    = float;
using PCLColorType     = uint8_t;
using PCLEnrichedPoint = pcl::PointXYZRGBNormal;
using PCLPointCloud    = pcl::PointCloud< PCLEnrichedPoint >;

// Point, Normal, Color types with constructors are needed to write
// code like 'put(pm, vd, Point(1, 2, 3))'
using PCLPoint  = Eigen::Vector3f;
using PCLVector = Eigen::Vector3f;
using PCLColor  = Eigen::Matrix<PCLColorType, 3, 1>; // R,G,B

} // namespace FEVV

