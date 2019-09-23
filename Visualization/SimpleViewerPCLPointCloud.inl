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

#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#include "FEVV/Wrappings/Graph_properties_pcl_point_cloud.h"


namespace FEVV {

//
// SimpleViewer<>::internal_createMesh(...) specialization for FEVV::PCLPointCloud
// with  HalfedgeGraph = FEVV::PCLPointCloud
// and   PointMap = FEVV::PCLPointCloudPointMap
//
template<>
inline void
SimpleViewer::internal_createMesh< FEVV::PCLPointCloud,
                                   FEVV::PCLPointCloudPointMap >(
    osg::Geode *&geode,
    FEVV::PCLPointCloud *_g,
    PMapsContainer *_pmaps,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesL,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometriesP,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_edges,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_vertices,
    std::vector< osg::ref_ptr< osg::Geometry > > &geometries_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &vertexArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArraysF,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &normalsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec3Array > > &tangentsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_edges,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_vertices,
    std::vector< osg::ref_ptr< osg::Vec4Array > > &colorsArrays_normals,
    std::vector< osg::ref_ptr< osg::Vec2Array > > &texcoordsArrays,
    FEVV::PCLPointCloudPointMap *_pm,
    std::string _mesh_file)
{
  internal_createMesh_pointcloud(geode,
                                 _g,
                                 _pmaps,
                                 geometries,
                                 geometriesL,
                                 geometriesP,
                                 geometries_edges,
                                 geometries_vertices,
                                 geometries_normals,
                                 vertexArrays,
                                 vertexArrays_edges,
                                 vertexArrays_vertices,
                                 vertexArrays_normals,
                                 normalsArrays,
                                 normalsArraysF,
                                 normalsArrays_edges,
                                 normalsArrays_vertices,
                                 tangentsArrays,
                                 colorsArrays,
                                 colorsArrays_edges,
                                 colorsArrays_vertices,
                                 colorsArrays_normals,
                                 texcoordsArrays,
                                 _pm,
                                 _mesh_file);
}

} // namespace FEVV