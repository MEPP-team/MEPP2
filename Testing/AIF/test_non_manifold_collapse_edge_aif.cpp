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
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshWriter.hpp"

//#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/Graph_traits_aif.h"
//#include "FEVV/Wrappings/Graph_properties_aif.h"

#include <fstream>
#include <iostream>

#include "Testing/Utils/utils_retrieve_edge.h"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "FEVV/Operators/AIF/collapse_edge.hpp"
// DBG #include "FEVV/Filters/Generic/print_points.hpp"

using namespace FEVV;
using namespace FEVV::Operators;

//------------------------------------------------------------------------------

void
test_collapse_edge_aif(std::string filename,
                       int source_index,
                       int target_index,
                       const std::string &output_file_name)
{
  typedef FEVV::DataStructures::AIF::AIFMesh Mesh;
  typedef FEVV::DataStructures::AIF::AIFMesh::ptr_mesh ptr_mesh;

  // Load a mesh
  ptr_mesh pm;

  FEVV::DataStructures::AIF::AIFMeshReader in;
  if(!(pm = in.read(filename)))
  {
    std::cout << "failed";
    return;
  }

  typedef boost::graph_traits< Mesh > GraphTraits;
  typedef GraphTraits::edge_descriptor edge_descriptor;

  // DBG print_points(M, get(CGAL::vertex_point, *pm));

  edge_descriptor e = retrieve_edge(*pm, source_index, target_index);
  if(e == GraphTraits::null_edge())
  {
    std::cout << "Failed to retrieve edge from " << source_index << " to "
              << target_index << "." << std::endl;
    std::cout << "Exiting";
    exit(EXIT_FAILURE);
  }

  std::cout << "Collapsing edge " << source_index << " to " << target_index
            << "." << std::endl;

  collapse_edge_keep_target(*pm, e);

  // DBG print_points(*pm, get(CGAL::vertex_point, M));

  FEVV::DataStructures::AIF::AIFMeshWriter out;
  try
  {
    out.write(pm, output_file_name);
  }
  catch(...)
  {
    std::cout << "writing failed" << std::endl;
    return;
  }
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 3 || narg > 5)
  {
    std::cout << "Usage: " << argv[0]
              << " filename int [filenameresult [filenameexpectedresult]]; int "
                 "is either 0, 1, 2 or 3."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  std::string to_treat = argv[1];

  int halfedge_case = std::stoi(std::string(argv[2]));
  std::string output_file_name = std::string(argv[3]);

  if(halfedge_case == 0)
  {
    test_collapse_edge_aif(to_treat, 1, 2, output_file_name); // border halfedge
  }
  else if(halfedge_case == 1)
  {
    test_collapse_edge_aif(to_treat, 2, 1, output_file_name); // border halfedge
  }
  else if(halfedge_case == 2)
  {
    test_collapse_edge_aif(to_treat, 3, 4, output_file_name);
  }
  else if(halfedge_case == 3)
  {
    test_collapse_edge_aif(to_treat, 4, 3, output_file_name);
  }
  else if(halfedge_case == 4)
  {
    test_collapse_edge_aif(to_treat, 5, 7, output_file_name);
  }
  else if(halfedge_case == 5)
  {
    test_collapse_edge_aif(to_treat, 7, 5, output_file_name);
  }
  else if(halfedge_case == 6)
  { // used for particular tests... (not for unit testing, but more for
    // intensive local testing)

    // for tetra.off
    test_collapse_edge_aif(to_treat, 0, 1, output_file_name);
    test_collapse_edge_aif(to_treat, 1, 2, output_file_name);
    test_collapse_edge_aif(to_treat, 0, 2, output_file_name);
    test_collapse_edge_aif(to_treat, 2, 3, output_file_name);

    // for flip14.vtk:
    // testCollapseEdgeAIF(toTreat, 0, 1, outputFileName);
    // testCollapseEdgeAIF(toTreat, 1, 2, outputFileName);
    // testCollapseEdgeAIF(toTreat, 2, 3, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3, 4, outputFileName);
    // testCollapseEdgeAIF(toTreat, 4, 3, outputFileName);

    // for bunny.vtk:
    // testCollapseEdgeAIF(toTreat, 3058, 4830, outputFileName);
    // testCollapseEdgeAIF(toTreat, 274, 3220, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3007, 4231, outputFileName);
    // testCollapseEdgeAIF(toTreat, 4892, 2359, outputFileName);
    // testCollapseEdgeAIF(toTreat, 530, 1282, outputFileName);
    // testCollapseEdgeAIF(toTreat, 609, 3423, outputFileName);
    // testCollapseEdgeAIF(toTreat, 4041, 33, outputFileName);
    // testCollapseEdgeAIF(toTreat, 4710, 4707, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3867, 2009, outputFileName);
    // testCollapseEdgeAIF(toTreat, 2009, 2021, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3897, 2826, outputFileName);
    // testCollapseEdgeAIF(toTreat, 2826, 2638, outputFileName);
    // testCollapseEdgeAIF(toTreat, 2638, 3119, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3119, 3897, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3557, 3681, outputFileName);
    // testCollapseEdgeAIF(toTreat, 3681, 603, outputFileName);
  }
  else
  {
    return 1; // test failed
  }

  if(narg == 5)
  {
    if(!are_meshes_equal(
           output_file_name,
           argv[4],
           true)) // this test currently fails for case 0 [border edge, because
                  // this is not handled for the time being]
      return 1;   // test failed
  }

  return 0; // test succeeded
}
