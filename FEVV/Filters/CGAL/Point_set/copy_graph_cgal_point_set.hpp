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

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
#include "FEVV/Filters/Generic/PointCloud/copy_point_cloud.hpp"
#include "FEVV/Filters/Generic/copy_graph.hpp"


namespace FEVV {
namespace Filters {


/**
 * \brief  Specialization of CopyGraphParameters for CGAL Point Set copy from.
 * 
 */
template< typename PointCloudT >
struct CopyGraphParameters< FEVV::CGALPointSet, PointCloudT > :
    CopyPCParameters< FEVV::CGALPointSet, PointCloudT >
{ };


/**
 * \brief  Overloading of copy_graph() for CGAL Point Set copy from.
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename PointCloudT,
          typename Parameters,
            // = CopyGraphParameters< FEVV::CGALPointSet, PointCloudT > >
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const FEVV::CGALPointSet   &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           PointCloudT                &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params,
           const GeometryTraitsS      &gt_s,
           const GeometryTraitsT      &gt_t)
{
  FEVV::Filters::copy_point_cloud(
      pc_s, pmaps_s, pc_t, pmaps_t, params, gt_s, gt_t);
}


/**
 * \brief  Specialization of CopyGraphParameters for CGAL Point Set copy to.
 * 
 */
template< typename PointCloudT >
struct CopyGraphParameters< PointCloudT, FEVV::CGALPointSet > :
    CopyPCParameters< PointCloudT, FEVV::CGALPointSet >
{ };


/**
 * \brief  Overloading of copy_graph() for CGAL Point Set copy to.
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename PointCloudT,
          typename Parameters,
            // = CopyGraphParameters< PointCloudT, FEVV::CGALPointSet > >
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const PointCloudT          &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::CGALPointSet         &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params,
           const GeometryTraitsS      &gt_s,
           const GeometryTraitsT      &gt_t)
{
  FEVV::Filters::copy_point_cloud(
      pc_s, pmaps_s, pc_t, pmaps_t, params, gt_s, gt_t);
}


/**
 * \brief  Specialization of CopyGraphParameters for copy from CGAL Point Set
 *         to CGAL Point Set.
 * 
 */
template< >
struct CopyGraphParameters< FEVV::CGALPointSet , FEVV::CGALPointSet > :
    CopyPCParameters< FEVV::CGALPointSet, FEVV::CGALPointSet >
{ };


/**
 * \brief  Overloading of copy_graph() for copy from CGAL Point Set
 *         to CGAL Point Set.
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename Parameters,
            // = CopyGraphParameters< FEVV::CGALPointSet, FEVV::CGALPointSet > >
          typename GeometryTraitsS,
          typename GeometryTraitsT >
void
copy_graph(const FEVV::CGALPointSet   &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::CGALPointSet         &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params,
           const GeometryTraitsS      &gt_s,
           const GeometryTraitsT      &gt_t)
{
  FEVV::Filters::copy_point_cloud(
      pc_s, pmaps_s, pc_t, pmaps_t, params, gt_s, gt_t);
}


} // namespace Filters
} // namespace FEVV
