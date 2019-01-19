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

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

// Boost properties adaptor
#define CGAL_USE_OM_POINTS
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

#include "FEVV/Wrappings/Geometry_traits_openmesh.h"

using namespace FEVV;

void
usage(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off file." << std::endl;
    std::cout << "Note: start with providing Testing/Data/flat.off as filename"
              << std::endl;
    exit(EXIT_FAILURE);
  }
}

int
main(int narg, char **argv)
{
  usage(narg, argv);

  typedef OpenMesh::PolyMesh_ArrayKernelT< /* MyTraits*/ > Mesh;

  // Reading a mesh (that we assume has some associated geometry)
  Mesh m;
  if(!OpenMesh::IO::read_mesh(m, argv[1]))
    return 1; // OpenMesh failed reading file

  /// [Snippet Accessing Geometry]
  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Point Point;
  const Geometry g(m);

  typedef boost::graph_traits< Mesh >::vertex_iterator vertex_iterator;

  // Under the hood the following auto type corresponds to something like
  //   boost::property_map< Mesh, boost::vertex_point_t >::type
  auto pm = get(boost::vertex_point, m);
  vertex_iterator vi;
  for(vi = vertices(m).first; vi != vertices(m).second; ++vi)
  {
    auto p = pm[*vi]; // Meaning get( pm, *vi );
    p = Point(g.get_x(p), g.get_y(p), g.get_z(p));
    // ... do your stuff with the geometry
  }

  /// [Snippet Accessing Geometry]

  std::cout << "Done." << std::endl;
  return 0;
}
