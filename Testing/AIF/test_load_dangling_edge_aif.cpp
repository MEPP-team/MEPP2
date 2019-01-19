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
//
// Mesh loading test using either generic reader
// or native readers.

#include <iostream>
#include <fstream>
//---------------------------------------------------------

// uncomment one at a time
//#define  WITH_GENERIC_READER
//#define  WITH_NATIVE_READER

// uncomment one at a time
//#define  USE_POLYHEDRON
//#define  USE_SURFACEMESH
//#define  USE_LCC
//#define  USE_OPENMESH
//#define  USE_AIF

//---------------------------------------------------------
#if defined(FEVV_USE_CGAL) &&                                                  \
    (defined(USE_POLYHEDRON) || defined(USE_SURFACEMESH) || defined(USE_LCC))
#include "FEVV/DataStructures/DataStructures_cgal.h"
#ifdef WITH_GENERIC_READER
#include "FEVV/Wrappings/Geometry_traits_cgal_linear_cell_complex.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Graph_traits_extension_cgal_polyhedron_3.h"
#endif
#endif

#if defined(FEVV_USE_OPENMESH) && defined(USE_OPENMESH)
#include "FEVV/DataStructures/DataStructures_openmesh.h"
#ifdef WITH_GENERIC_READER
#include "FEVV/Wrappings/Geometry_traits_openmesh.h"
#include "FEVV/Wrappings/Graph_traits_extension_openmesh.h"
#endif
#endif

#if defined(FEVV_USE_AIF) && defined(USE_AIF)
#include "FEVV/DataStructures/DataStructures_aif.h"
// \todo Inquire on the asymetry of the inclusion of the Graph_traits*.h
//       between SurfaceMesh and OpenMesh (that only include the GraphTraits
//       when in the WITH_GENERIC_READER section) and AIF.
#include "FEVV/Wrappings/Graph_traits_aif.h"
#include "FEVV/Wrappings/Graph_traits_extension_aif.h"
#ifdef  WITH_GENERIC_READER
#include "FEVV/Wrappings/Geometry_traits_aif.h"
#include "FEVV/Wrappings/Graph_properties_aif.h"
#endif
#endif

#ifdef WITH_GENERIC_READER
#include "FEVV/Filters/Generic/generic_reader.hpp"

std::string reader_type("generic reader");
#else
#ifdef WITH_NATIVE_READER
std::string reader_type("native reader");
#else
std::string reader_type("unkonwn reader");
#endif
#endif

#ifdef USE_POLYHEDRON
using MeshT = FEVV::MeshPolyhedron;
std::string mesh_type("Polyhedron");
#endif

#ifdef USE_SURFACEMESH
using MeshT = FEVV::MeshSurface;
std::string mesh_type("Surface_mesh");
#endif

#ifdef USE_LCC
using MeshT = FEVV::MeshLCC;
std::string mesh_type("LCC");
#endif

#ifdef USE_OPENMESH
using MeshT = FEVV::MeshOpenMesh;
std::string mesh_type("Openmesh");
#endif

#ifdef USE_AIF
using MeshT = FEVV::MeshAIF;
std::string mesh_type("AIF");
#endif

//---------------------------------------------------------

#ifdef WITH_GENERIC_READER

// with generic reader
void
load_mesh(std::string filename)
{
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(filename, m, pmaps_bag);

  std::cout << "the mesh has " << num_vertices(m) << " vertice(s), "
            << num_edges(m) << " edge(s), "
#if !defined(USE_AIF)
            << num_halfedges(m) << " halfedge(s), "
#endif
            << num_faces(m) << " face(s)." << std::endl;

  std::cout << "Done." << std::endl;
}

#else

// with native reader
void
load_mesh(std::string filename)
{
#ifdef USE_AIF
  typedef FEVV::MeshAIFPtr ptr_mesh;
  // Load a mesh
  ptr_mesh pm;

  FEVV::DataStructures::AIF::AIFMeshReader in;
  if(!(pm = in.read(filename)))
  {
    std::cout << "failed";
    return;
  }
  MeshT m(*pm);
#else
#ifdef USE_OPENMESH

  MeshT m;
  if(!OpenMesh::IO::read_mesh(m, filename))
  {
    std::cout << "failed";
    return;
  }

#else
  std::ifstream in(filename);
  if(!in)
  {
    std::cout << "Unable to read file " << filename << std::endl;
    exit(EXIT_FAILURE);
  }
  MeshT m;
#ifndef USE_LCC
  // Polyhedron, Surface_mesh
  in >> m;
#else
  // USE_LCC
  CGAL::load_off(m, in);
#endif // USE_LCC

#endif // USE_OPENMESH

#endif // USE_AIF

  std::cout << "the mesh has " << num_vertices(m) << " vertice(s), "
            << num_edges(m) << " edge(s), "
#if !defined(USE_AIF)
            << num_halfedges(m) << " halfedge(s), "
#endif
            << num_faces(m) << " face(s)." << std::endl;

  std::cout << "Done." << std::endl;
}

#endif


int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0]
              << " filename; filename being an off or vtk file." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout << "Loading mesh using " << reader_type << " and " << mesh_type
            << std::endl;

  load_mesh(argv[1]);

  return 0;
}
