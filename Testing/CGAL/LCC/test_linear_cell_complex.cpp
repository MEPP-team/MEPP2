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
#include <fstream>

#include <CGAL/Cartesian.h>

// Linear_cell_complex
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"

#include <CGAL/Linear_cell_complex_incremental_builder.h>

#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/scaling.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertices_one_ring_barycenter.hpp"
#include "FEVV/Filters/Generic/reposition_vertices.hpp"
#include "FEVV/Operators/Generic/Manifold/collapse_edge.hpp"
#include "FEVV/Filters/Generic/print_points.hpp"
#include "FEVV/Filters/Generic/print_face_normals.hpp"

using namespace FEVV;
using namespace FEVV::Filters;

typedef CGAL::Cartesian< double > Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_traits< 3, Kernel > MyTraits;

struct Myitem
{
  template< class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Tag_true Darts_with_id;
    typedef CGAL::Cell_attribute_with_point_and_id< Refs > Vertex_attribute;
    typedef CGAL::Cell_attribute_with_id< Refs > Face_attribute;
    typedef CGAL::cpp11::tuple< Vertex_attribute, void, Face_attribute >
        Attributes;
  };
};

typedef CGAL::
    Linear_cell_complex_for_combinatorial_map< 2, 3, MyTraits, Myitem >
        LCC;

//------------------------------------------------------------------------------
void
test_linear_cell_complex(char *filename)
{
  // 1) Load a mesh
  std::ifstream in(filename);
  LCC lcc;
  try
  {
    CGAL::read_off(in, lcc);
    lcc.display_characteristics(std::cout);
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[LCC] Exception caught while reading input file " << filename
              << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(false, "[LCC] Exception caught while reading input file.");
  }

  auto pos_pm = get(boost::vertex_point, lcc);

  // 2) Add face normal property
  typedef boost::property_map< LCC, boost::face_index_t >::const_type
      Face_index_map;
  boost::vector_property_map< Vector, Face_index_map > normals(
      get(boost::face_index, lcc));

  // 3) Compute normals
  calculate_face_normals(lcc // Graph
                         ,
                         pos_pm // map from vertex_descriptor to point
                         ,
                         normals // map from face_descriptor to Vector_3
  );

  // 4) Draw normals
  print_face_normals(lcc, normals);

  // 5) Scale the mesh
  calculate_scaling(lcc, pos_pm, 2.0, 3.0, 4.0);

  // 6) Smoothing
  typedef FEVV::Geometry_traits< LCC > Geometry;
  typedef Geometry::Point Point;
  typedef boost::vector_property_map<
      Point,
      typename boost::property_map< LCC, boost::vertex_index_t >::const_type >
      BarycenterMap;

  BarycenterMap barycenters_pm(get(boost::vertex_index, lcc));
  FEVV::Filters::calculate_vertices_one_ring_barycenter(
      lcc, pos_pm, barycenters_pm);
  FEVV::Filters::reposition_vertices(lcc, pos_pm, barycenters_pm);

  // 7) Display the points
  print_points(lcc, pos_pm);

  // 8) Collapse some edges
  /*std::cout << "Start collapsing edges ..." << std::endl;
  collapse_some_edges(lcc);
  std::ofstream collapse_some_edges_os("collapse_some_edges.off");
  collapse_some_edges_os << M;*/


  std::cout << "Done." << std::endl;
}

//------------------------------------------------------------------------------
int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: test-linear-cell-complex filename;"
                 " filename being an off file."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  test_linear_cell_complex(argv[1]);
  return EXIT_SUCCESS;
}
