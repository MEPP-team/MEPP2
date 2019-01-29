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
#include "FEVV/Wrappings/Geometry_traits_openmesh.h" // Always in first position for MVS 4996 warnings

#include <OpenMesh/Core/IO/MeshIO.hh>
// Boost properties adaptor
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

//#define USE_GENERIC_READER_WRITER // when activated you need to link to vtk
//libs if VTK is used
#ifdef USE_GENERIC_READER_WRITER
#include "FEVV/Wrappings/Graph_traits_extension_openmesh.h"
#include "FEVV/Wrappings/properties_openmesh.h"
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#endif

#include "FEVV/Tools/IO/FileUtilities.hpp"

#include "FEVV/Filters/Generic/Manifold/calculate_vertices_one_ring_barycenter.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertices_one_ring_angles_based_centroid.hpp"
#include "FEVV/Filters/Generic/reposition_vertices.hpp"

#include <fstream>

#define USE_OPENMESH_DOUBLE_TRAITS
#ifdef USE_OPENMESH_DOUBLE_TRAITS
struct MyTraitsDouble : public OpenMesh::DefaultTraits
{
  typedef OpenMesh::Vec3d Point; // use double-values points
  typedef OpenMesh::Vec3d Normal;
};
#endif

using namespace FEVV;
using namespace FEVV::Filters;

//------------------------------------------------------------------------------
void
test_smoothing_open_mesh(char *filename)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<
#ifdef USE_OPENMESH_DOUBLE_TRAITS
      MyTraitsDouble
#endif
      >
      Mesh;
  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Point Point;
  typedef boost::vector_property_map<
      Point,
      typename boost::property_map< Mesh, boost::vertex_index_t >::const_type >
      BarycenterMap;
  /**********************************************************************************************************/
  Mesh m;
#ifdef USE_GENERIC_READER_WRITER
  FEVV::PMapsContainer pmaps;
  FEVV::Filters::read_mesh(filename, m, pmaps);
#else
  std::ifstream in(filename);
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "reading failed";
    return;
  }
#endif
  /**********************************************************************************************************/
  auto pos_pm = get(boost::vertex_point, m);
  /**********************************************************************************************************/
  BarycenterMap barycenters_pm(get(boost::vertex_index, m));
  const auto pos_pm_const = get(boost::vertex_point, m);
  FEVV::Filters::calculate_vertices_one_ring_barycenter(
      m, pos_pm_const, barycenters_pm, 0.2f);
  FEVV::Filters::reposition_vertices(m, pos_pm, barycenters_pm);

  FEVV::Filters::calculate_vertices_one_ring_angles_based_centroid(
      m, pos_pm, barycenters_pm, 0.2f);
  FEVV::Filters::reposition_vertices(m, pos_pm, barycenters_pm);
  /**********************************************************************************************************/
#ifdef USE_GENERIC_READER_WRITER
  FEVV::Filters::write_mesh("smoothed_om_" +
                                FEVV::FileUtils::get_file_name(filename) +
                                FEVV::FileUtils::get_file_extension(filename),
                            m,
                            pmaps);
#else
  if(!OpenMesh::IO::write_mesh(
         m,
         "smoothed_om_" + FEVV::FileUtils::get_file_name(filename) +
             FEVV::FileUtils::get_file_extension(filename)))
  {
    std::cout << "writing failed";
    return;
  }
#endif
  std::cout << "Done." << std::endl;
}

//------------------------------------------------------------------------------
int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    exit(EXIT_FAILURE);
  }

  test_smoothing_open_mesh(argv[1]);
  return 0;
}
