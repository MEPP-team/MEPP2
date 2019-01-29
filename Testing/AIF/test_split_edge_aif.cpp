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
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Geometry_traits_aif.h"

#include "FEVV/Wrappings/properties_aif.h"

using AIFKernel = FEVV::AIF_mesh_kernel_generator;
using MeshT = FEVV::DataStructures::AIF::AIFMesh;

//---------------- specific code above --------------------

//---------------- generic code below ---------------------

#include "Testing/Utils/utils_retrieve_halfedge.h"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "FEVV/Operators/Generic/Manifold/split_edge.hpp"

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

//---------------- generic code above ---------------------

//---------------- other code below -----------------------
#include <fstream>
#include <iostream>

void
test_split_edge_aif(const std::string &input_file_path,
                    int source_index,
                    int target_index,
                    const std::string &output_file_path)
{
  typedef boost::graph_traits< MeshT > GraphTraits;
  typedef GraphTraits::halfedge_descriptor halfedge_descriptor;

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // some check and filter application
  halfedge_descriptor h = retrieve_halfedge(m, source_index, target_index);
  if(h == GraphTraits::null_halfedge())
  {
    std::cout << "Failed to retrieve edge from " << source_index << " to "
              << target_index << "." << std::endl;
    std::cout << "Exiting";
    exit(EXIT_FAILURE);
  }

  std::cout << "Splitting edge " << source_index << " to " << target_index
            << "." << std::endl;
  FEVV::Operators::split_edge(m, get(boost::vertex_point, m), h);

  // write mesh to file
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);
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
    test_split_edge_aif(to_treat, 1, 2, output_file_name); // border halfedge
  }
  else if(halfedge_case == 1)
  {
    test_split_edge_aif(to_treat, 2, 1, output_file_name); // border halfedge
  }
  else if(halfedge_case == 2)
  {
    test_split_edge_aif(to_treat, 3, 4, output_file_name);
  }
  else if(halfedge_case == 3)
  {
    test_split_edge_aif(to_treat, 4, 3, output_file_name);
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
