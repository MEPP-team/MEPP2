// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
// Graph traits adaptors
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <fstream>
#include <string> // std::stoi

#include "Testing/Utils/utils_retrieve_halfedge.h"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "FEVV/Operators/Generic/Manifold/collapse_edge_euler.hpp"
// DBG #include "FEVV/Filters/print_points.h"

using namespace FEVV;
using namespace FEVV::Operators;

/*
 * \brief  Tests the \link collapse_edge_keep_target test \endlink
 * \param  sourceIndex The index (position) of the source vertex of the edge
 *                     that is seeked for deletion.
 * \param  targetIndex The index (position) of the destination vertex of the
 * edge that is seeked for deletion.
 */
void
test_collapse_edge_polyhedron(std::ifstream &in,
                              int source_index,
                              int target_index,
                              const std::string &output_file_name)
{
  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Polyhedron_3< Kernel, CGAL::Polyhedron_items_with_id_3 >
      Polyhedron;

  // Load a mesh
  Polyhedron p;
  in >> p;

  typedef boost::graph_traits< Polyhedron > GraphTraits;
  typedef typename GraphTraits::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor h = retrieve_halfedge(p, source_index, target_index);
  if(h == GraphTraits::null_halfedge())
  {
    std::cout << "Failed to retrieve edge from " << source_index << " to "
              << target_index << "." << std::endl;
    std::cout << "Exiting";
    exit(EXIT_FAILURE);
  }

  std::cout << "Collapsing edge " << source_index << " to " << target_index
            << "." << std::endl;
  // DBG FEVV::Filters::print_points(p, get(CGAL::vertex_point, p));
  FEVV::Operators::collapse_edge_keep_target_euler(p, h);
  // DBG FEVV::Filters::print_points(p, get(CGAL::vertex_point, p));

  std::ofstream os(output_file_name);
  os << p;
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 3 || narg > 4)
  {
    //    0 --------- 1 --------- 2
    //    |\         / \         /|\
    //    | \       /   \       / | \
    //    |  \     /     \     /  |  \
    //    |   \   /       \   /   |   \
    //    |    \ /         \ /    |    \
    //    |     3 --------- 4     |     5
    //    |    / \         / \    |    /
    //    |   /   \       /   \   |   /
    //    |  /     \     /     \  |  /
    //    | /       \   /       \ | /
    //    |/         \ /         \|/
    //    6 --------- 7 --------- 8

    std::cout << "Usage: a.out filename_a deletion_case filename_b";
    std::cout << " - filename_a filename of file to be treated (off format)"
              << std::endl;
    std::cout << " - deletion case: integer specifying the edge to delete"
              << std::endl;
    std::cout << "      0 for unspecified edge" << std::endl;
    std::cout << "      1 for edge 3 --> 4" << std::endl;
    std::cout << "      2 for edge 4 --> 3" << std::endl;
    std::cout << "      3 for edge 6 --> 3" << std::endl;
    std::cout << "      4 for edge 6 --> 7" << std::endl;
    std::cout << " - filename_b optional reference off file of valid result."
              << std::endl;

    exit(EXIT_FAILURE);
  }


  std::ifstream to_treat(argv[1]);
  if(!to_treat)
  {
    std::cout << "Unable to read file " << argv[1] << std::endl;
    exit(EXIT_FAILURE);
  }

  int deletion_case = std::stoi(std::string(argv[2]));
  std::string output_file_name = std::string(argv[0]) + argv[2] + ".off";

  if(deletion_case == 0)
  {
    test_collapse_edge_polyhedron(
        to_treat, 0, 1, output_file_name); // Whatever comes first
  }
  else if(deletion_case == 1)
  {
    test_collapse_edge_polyhedron(to_treat, 3, 4, output_file_name);
  }
  else if(deletion_case == 2)
  {
    test_collapse_edge_polyhedron(to_treat, 4, 3, output_file_name);
  }
  else if(deletion_case == 3)
  {
    test_collapse_edge_polyhedron(to_treat, 6, 3, output_file_name);
  }
  else if(deletion_case == 4)
  {
    test_collapse_edge_polyhedron(to_treat, 6, 7, output_file_name);
  }
  else
  {
    std::cout << "Unknown deletation case " << argv[2] << std::endl;
  }

  to_treat.close();

  if(narg == 4)
  {
    if(!are_meshes_equal(output_file_name, argv[3], false))
      return 1; // test failed
  }

  return 0;
}
