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
#include "FEVV/Filters/Generic/PointCloud/copy_point_cloud.hpp"
#include "FEVV/Filters/Generic/copy_graph.hpp"


namespace FEVV {
namespace Filters {


/**
 * \brief  Specialization of CopyGraphParameters for PCL Point Cloud copy from.
 * 
 */
template< typename PointCloudT >
struct CopyGraphParameters< FEVV::PCLPointCloud, PointCloudT > :
    CopyPCParameters< FEVV::PCLPointCloud, PointCloudT >
{ };


/**
 * \brief  Overloading of copy_graph() for PCL Point Cloud (copy from).
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename PointCloudT,
          typename Parameters =
              CopyGraphParameters< FEVV::PCLPointCloud, PointCloudT > >
void
copy_graph(const FEVV::PCLPointCloud  &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           PointCloudT                &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params =
               CopyGraphParameters< FEVV::PCLPointCloud, PointCloudT >())
{
  FEVV::Filters::copy_point_cloud(pc_s, pmaps_s, pc_t, pmaps_t, params);
}


/**
 * \brief  Specialization of CopyGraphParameters for PCL Point Cloud copy to.
 * 
 */
template< typename PointCloudT >
struct CopyGraphParameters< PointCloudT, FEVV::PCLPointCloud > :
    CopyPCParameters< PointCloudT, FEVV::PCLPointCloud >
{ };


/**
 * \brief  Overloading of copy_graph() for PCL Point Cloud (copy to).
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename PointCloudT,
          typename Parameters =
              CopyGraphParameters< PointCloudT, FEVV::PCLPointCloud > >
void
copy_graph(const PointCloudT          &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::PCLPointCloud        &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params =
               CopyGraphParameters< PointCloudT, FEVV::PCLPointCloud >())
{
  FEVV::Filters::copy_point_cloud(pc_s, pmaps_s, pc_t, pmaps_t, params);
}


/**
 * \brief  Specialization of CopyGraphParameters for copy from PCL Point Cloud
 *         to PCL Point Cloud.
 * 
 */
template< >
struct CopyGraphParameters< FEVV::PCLPointCloud , FEVV::PCLPointCloud > :
    CopyPCParameters< FEVV::PCLPointCloud, FEVV::PCLPointCloud >
{ };


/**
 * \brief  Overloading of copy_graph() for PCL Point Cloud (copy from & to).
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename Parameters =
              CopyGraphParameters< FEVV::PCLPointCloud, FEVV::PCLPointCloud > >
void
copy_graph(const FEVV::PCLPointCloud  &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::PCLPointCloud        &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params =
               CopyGraphParameters< FEVV::PCLPointCloud, FEVV::PCLPointCloud >())
{
  FEVV::Filters::copy_point_cloud(pc_s, pmaps_s, pc_t, pmaps_t, params);
}


} // namespace Filters
} // namespace FEVV


// we must explicitely define the cases "copy from PCL to CGALPS" and 
// "copy from CGALPS to PCL" to avoid an "ambiguous call error" because the
// compiler can not choose for example between "copy from PCL" and
// "copy to CGALPS" ; this leads to PCL/CGAL mixed code, which is bad, but we
// have no other choice
#ifdef FEVV_USE_CGAL

#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
#include "FEVV/Wrappings/Graph_traits_cgal_point_set.h" // for vertex_decriptor
#include "FEVV/Wrappings/Geometry_traits_cgal_point_set.h" // for CGALPS Kernel
#include "FEVV/Wrappings/Graph_properties_cgal_point_set.h" // for get(vertex_point_t,...)
#include "FEVV/Wrappings/properties_cgal_point_set.h" // for CGALPS prop maps

namespace FEVV {
namespace Filters {


/**
 * \brief  Specialization of CopyGraphParameters for copy from PCL Point Cloud
 *         to CGAL Point Set.
 * 
 */
template< >
struct CopyGraphParameters< FEVV::PCLPointCloud , FEVV::CGALPointSet > :
    CopyPCParameters< FEVV::PCLPointCloud, FEVV::CGALPointSet >
{ };


/**
 * \brief  Overloading of copy_graph() for copy from PCL Point Cloud
 *         to CGAL Point Set.
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename Parameters =
              CopyGraphParameters< FEVV::PCLPointCloud, FEVV::CGALPointSet > >
void
copy_graph(const FEVV::PCLPointCloud  &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::CGALPointSet         &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params =
               CopyGraphParameters< FEVV::PCLPointCloud, FEVV::CGALPointSet >())
{
  FEVV::Filters::copy_point_cloud(pc_s, pmaps_s, pc_t, pmaps_t, params);
}


/**
 * \brief  Specialization of CopyGraphParameters for copy from CGAL Point Set
 *         to PCL Point Cloud.
 * 
 */
template< >
struct CopyGraphParameters< FEVV::CGALPointSet, FEVV::PCLPointCloud > :
    CopyPCParameters< FEVV::CGALPointSet, FEVV::PCLPointCloud >
{ };


/**
 * \brief  Overloading of copy_graph() for copy from CGAL Point Set
 *         to PCL Point Cloud.
 * 
 */
// note: overloading because function partial specialization is not allowed.
template< typename Parameters =
              CopyGraphParameters< FEVV::CGALPointSet, FEVV::PCLPointCloud > >
void
copy_graph(const FEVV::CGALPointSet   &pc_s,
           const FEVV::PMapsContainer &pmaps_s,
           FEVV::PCLPointCloud        &pc_t,
           FEVV::PMapsContainer       &pmaps_t,
           const Parameters           &params =
               CopyGraphParameters< FEVV::CGALPointSet, FEVV::PCLPointCloud >())
{
  FEVV::Filters::copy_point_cloud(pc_s, pmaps_s, pc_t, pmaps_t, params);
}


} // namespace Filters
} // namespace FEVV

#endif // FEVV_USE_CGAL
