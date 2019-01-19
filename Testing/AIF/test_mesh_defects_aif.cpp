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
#include "FEVV/Wrappings/Graph_traits_aif.h"

#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshWriter.hpp"

#include "FEVV/Operators/AIF/similarity.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"

using namespace FEVV::DataStructures::AIF;

int
main(int narg, char **argv)
{
  typedef AIFMeshReader reader_type;
  typedef AIFMeshWriter writer_type;
  typedef boost::graph_traits< reader_type::output_type > GraphTraits;

  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  typedef typename GraphTraits::edge_iterator edge_iterator;
  typedef typename GraphTraits::edge_descriptor edge_descriptor;
  typedef typename GraphTraits::face_iterator face_iterator;
  typedef typename GraphTraits::face_descriptor face_descriptor;
  /////////////////////////////////////////////////////////////////////////////
  if(narg < 2 || narg > 4)
  {
    std::cerr << "Cannot proceed arguments. Please use " << argv[0]
              << " meshfile [removeIsolatedElement "
                 "[resolveVerticesWithSimilarIncidentEdges]]"
              << std::endl;
    return -1;
  }
  std::string input_file_path(argv[1]);
  std::string remove_isolated_elements(((narg > 2) ? argv[2] : "n"));
  std::string resolve_vertices_with_similar_incident_edges(
      ((narg > 3) ? argv[3] : "n"));
  reader_type my_reader;
  // reader_type::output_type  m;
  reader_type::ptr_output ptr_mesh;
  if(remove_isolated_elements != "y" && remove_isolated_elements != "n" &&
     remove_isolated_elements != "Y" && remove_isolated_elements != "N")
  {
    std::cerr << "Cannot understand input for removing of not isolated "
                 "element. Please use either \"y\" or \"n\" "
              << std::endl;
    return -1;
  }
  /////////////////////////////////////////////////////////////////////////////
  try
  {
    ptr_mesh = my_reader.read(input_file_path);
  }
  catch(const std::invalid_argument &e)
  {
    std::cerr << "Invalid Argument Exception catch while reading "
              << input_file_path << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(const std::ios_base::failure &e)
  {
    std::cerr << "IOS Failure Exception catch while reading " << input_file_path
              << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[AIF] Exception caught while reading input file "
              << input_file_path << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(false, "[AIF] Exception caught while reading input file.");
  }
  /////////////////////////////////////////////////////////////////////////////
  const AIFMesh::Vector red(255, 0, 0);
  const AIFMesh::Vector green(0, 255, 0);
  const AIFMesh::Vector blue(0, 0, 255);
  const AIFMesh::Vector white(255, 255, 255);
  if(!ptr_mesh->isPropertyMap< AIFVertex::ptr >("v:color"))
  { // create a property map to store vertex colors if not already created
    ptr_mesh->AddPropertyMap< AIFVertex::ptr, AIFMesh::Vector >("v:color");
    if(!ptr_mesh->isPropertyMap< AIFVertex::ptr >("v:color"))
      throw std::runtime_error("Failed to create vertex-color property map.");
  }
  if(!ptr_mesh->isPropertyMap< AIFEdge::ptr >("e:color"))
  { // create a property map to store edge colors if not already created
    ptr_mesh->AddPropertyMap< AIFEdge::ptr, AIFMesh::Vector >("e:color");
    if(!ptr_mesh->isPropertyMap< AIFEdge::ptr >("e:color"))
      throw std::runtime_error("Failed to create edge-color property map.");
  }
  if(!ptr_mesh->isPropertyMap< AIFFace::ptr >("f:color"))
  { // create a property map to store face colors if not already created
    ptr_mesh->AddPropertyMap< AIFFace::ptr, AIFMesh::Vector >("f:color");
    if(!ptr_mesh->isPropertyMap< AIFFace::ptr >("f:color"))
      throw std::runtime_error("Failed to create face-color property map.");
  }
  /////////////////////////////////////////////////////////////////////////////
  // 1 - Count number of isolated/complex/normal elements

  // VERTICES
  auto iterator_pair_v = vertices(*ptr_mesh);
  vertex_iterator vi = iterator_pair_v.first;
  vertex_iterator vi_end = iterator_pair_v.second;
  int nb_isolated_vertices = 0, nb_cut_vertices = 0,
      nb_vertices_with_similar_incident_edges = 0;
  std::vector< vertex_descriptor > v_to_remeove;
  for(; vi != vi_end; ++vi)
  {
    if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_vertex(*vi))
    {
      ++nb_isolated_vertices;

      if(remove_isolated_elements == "y" || remove_isolated_elements == "Y")
        v_to_remeove.push_back(*vi);
      else
      {
        ptr_mesh->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
            "v:color", (*vi)->GetIndex(), red);
      }
    }
    else if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_cut_vertex(*vi))
    {
      ++nb_cut_vertices;
      ptr_mesh->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
          "v:color", (*vi)->GetIndex(), blue);
    }
    else if(FEVV::Operators::contains_similar_incident_edges(
                *vi, *ptr_mesh))
    {

      ++nb_vertices_with_similar_incident_edges;

      if(resolve_vertices_with_similar_incident_edges == "y" ||
         resolve_vertices_with_similar_incident_edges == "Y")
        FEVV::Operators::resolve_similar_incident_edges(*vi, *ptr_mesh);
      else
        ptr_mesh->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
            "v:color", (*vi)->GetIndex(), blue);
    }
    else
    {
      ptr_mesh->SetProperty< AIFVertex::ptr, AIFMesh::Vector >(
          "v:color", (*vi)->GetIndex(), green);
    }
  }

  // EDGES
  auto iterator_pair_e = edges(*ptr_mesh);
  edge_iterator ei = iterator_pair_e.first;
  edge_iterator ei_end = iterator_pair_e.second;
  int nb_isolated_edges = 0, nb_dangling_edges = 0, nb_complex_edges = 0, 
      nb_border_edges = 0;
  std::vector< edge_descriptor > e_to_remeove;
  for(; ei != ei_end; ++ei)
  {
    if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_edge(*ei))
    {
      ++nb_isolated_edges;
      if(remove_isolated_elements == "y" || remove_isolated_elements == "Y")
        e_to_remeove.push_back(*ei);
      else
      {
        ptr_mesh->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
            "e:color", (*ei)->GetIndex(), red);
      }
    }
    else if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_dangling_edge(
                *ei))
    {
      ++nb_dangling_edges;
      ptr_mesh->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
          "e:color", (*ei)->GetIndex(), red);
    }
    else if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_complex_edge(*ei))
    {
      ++nb_complex_edges;
      ptr_mesh->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
          "e:color", (*ei)->GetIndex(), blue);
    }
    else if (FEVV::DataStructures::AIF::AIFTopologyHelpers::is_surface_border_edge(*ei))
    {
      ++nb_border_edges;
      ptr_mesh->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
          "e:color", (*ei)->GetIndex(), white);
    }
    else
    {
      ptr_mesh->SetProperty< AIFEdge::ptr, AIFMesh::Vector >(
          "e:color", (*ei)->GetIndex(), green);
    }
  }

  // FACES
  auto iterator_pair_f = faces(*ptr_mesh);
  face_iterator fi = iterator_pair_f.first;
  face_iterator fi_end = iterator_pair_f.second;
  int nb_isolated_faces = 0, nb_degenerated_faces = 0;
  std::vector< face_descriptor > f_to_remeove;
  for(; fi != fi_end; ++fi)
  {
    if(FEVV::DataStructures::AIF::AIFTopologyHelpers::is_isolated_face(*fi))
    {
      ++nb_isolated_faces;
      if(remove_isolated_elements == "y" || remove_isolated_elements == "Y")
        f_to_remeove.push_back(*fi);
      else
      {
        ptr_mesh->SetProperty< AIFFace::ptr, AIFMesh::Vector >(
            "f:color", (*fi)->GetIndex(), red);
      }
    }
    else
    {
      ptr_mesh->SetProperty< AIFFace::ptr, AIFMesh::Vector >(
          "f:color", (*fi)->GetIndex(), green);
    }

    if (FEVV::DataStructures::AIF::AIFTopologyHelpers::is_degenerated_face(*fi))
      ++nb_degenerated_faces;
  }
  /////////////////////////////////////////////////////////////////////////////
  // 2 - Remove isolated elements
  std::string suffix = "";
  if(remove_isolated_elements == "y" || remove_isolated_elements == "Y")
  {
    if(nb_isolated_vertices != 0 || nb_isolated_edges != 0 ||
       nb_isolated_faces != 0)
    {
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_faces(
          f_to_remeove.begin(), f_to_remeove.end(), ptr_mesh);
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_edges(
          e_to_remeove.begin(), e_to_remeove.end(), ptr_mesh);
      FEVV::DataStructures::AIF::AIFTopologyHelpers::remove_vertices(
          v_to_remeove.begin(), v_to_remeove.end(), ptr_mesh);
      suffix = "_without_isolated_elmt";
      if(resolve_vertices_with_similar_incident_edges == "y" ||
         resolve_vertices_with_similar_incident_edges == "Y")
        suffix += "_and_similar_incident_edges";
    }
    else if(resolve_vertices_with_similar_incident_edges == "y" ||
            resolve_vertices_with_similar_incident_edges == "Y")
      suffix = "_without_similar_incident_edges";
  }
  /////////////////////////////////////////////////////////////////////////////
  // 3 - Save corrected mesh
  writer_type my_writer;
  try
  {
    if(FEVV::FileUtils::get_file_extension(input_file_path) == ".ply")
      my_writer.write(ptr_mesh,
                      FEVV::FileUtils::get_file_name(input_file_path) + suffix +
                          ".off");
    else
      my_writer.write(ptr_mesh,
                      FEVV::FileUtils::get_file_name(input_file_path) + suffix +
                          FEVV::FileUtils::get_file_extension(input_file_path));
  }
  catch(...)
  {
    std::cout << "writing failed";
    return -1;
  }
  /////////////////////////////////////////////////////////////////////////////
  std::cout << "The mesh file named "
            << FEVV::FileUtils::get_file_name(input_file_path) +
                   FEVV::FileUtils::get_file_extension(input_file_path);
  if(nb_isolated_vertices > 0 || nb_cut_vertices > 0 || nb_isolated_edges > 0 ||
     nb_dangling_edges > 0 || nb_complex_edges > 0)
    std::cout << " is not 2-manifold " << std::endl;
  else
    std::cout << " is 2-manifold" << std::endl;

  std::string prefix =
      (remove_isolated_elements == "y" || remove_isolated_elements == "Y")
          ? "Number of removed isolated"
          : "Number of isolated";
  std::cout << prefix << " vertices: " << nb_isolated_vertices << std::endl;
  std::cout << prefix << " edges: " << nb_isolated_edges << std::endl;
  std::cout << prefix << " faces: " << nb_isolated_faces << std::endl;
  /////////////////////////////////////////////////////////////////////////////
  prefix = (resolve_vertices_with_similar_incident_edges == "y" ||
            resolve_vertices_with_similar_incident_edges == "Y")
               ? "Number of resolved vertices with similar incident edges: "
               : "Number of vertices with similar incident edges: ";
  std::cout << prefix << nb_vertices_with_similar_incident_edges << std::endl;
  /////////////////////////////////////////////////////////////////////////////
  std::cout << "Number of cut vertices: " << nb_cut_vertices << std::endl;
  std::cout << "Number of dangling edges: " << nb_dangling_edges << std::endl;
  std::cout << "Number of complex edges: " << nb_complex_edges << std::endl;
  std::cout << "Number of surface border edges: " << nb_border_edges << std::endl;
  std::cout << "Number of degenerated faces: " << nb_degenerated_faces << std::endl;
  /////////////////////////////////////////////////////////////////////////////
  return 0;
}
