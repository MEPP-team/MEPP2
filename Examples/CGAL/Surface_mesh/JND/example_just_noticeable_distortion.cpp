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

#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"

#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"

#include "FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/just_noticeable_distortion.hpp"
#include "FEVV/Filters/Generic/minmax_map.h"
#include "FEVV/Filters/Generic/color_mesh.h"
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"
#include "FEVV/Filters/Generic/Manifold/calculate_vertex_normals.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <time.h>

//---------------------------------------------------------
//                       Main
//---------------------------------------------------------

// Main: load a mesh, apply the filter, write the mesh
int
main(int narg, char **argv)
{
  typedef boost::graph_traits< FEVV::MeshSurface > GraphTraits;
  typedef typename GraphTraits::vertex_iterator vertex_iterator;
  typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
  // input and output files
  // std::string inputFilePath =
  // "../../../../Testing/Data/CubeTriangleFaces.obj";
  std::string input_file_path = "sphere.obj";
  if(narg == 2 || narg == 3)
  {
    std::string reader = std::string(argv[1]);
    input_file_path = reader;
  }
  std::string output_file_path = "just_noticeable_distortion_filter.output.obj";

  int screen_width = 1920;
  int screen_height = 1080;
  double screen_size = 55.;
  double user_dist = 50.;
  int scene_width = 1080;
  double scene_fov = M_PI * 0.3333;
  int number_of_lights = 128;

  if(narg == 3)
  {
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(argv[2], pt);

    screen_width = pt.get< int >("screen_width", screen_width);
    screen_height = pt.get< int >("screen_height", screen_height);
    screen_size = pt.get< double >("screen_size", screen_size);
    user_dist = pt.get< double >("user_dist", user_dist);
    scene_width = pt.get< int >("scene_width", scene_width);
    scene_fov = M_PI * pt.get< double >("scene_fov", scene_fov);
    number_of_lights = pt.get< int >("number_of_lights", number_of_lights);
  }
  std::cout << "Confing used :" << std::endl;
  std::cout << "\tScreen : " << screen_width << " * " << screen_height << " - "
            << screen_size << std::endl;
  std::cout << "\tUser distance : " << user_dist << std::endl;
  std::cout << "\tScene : " << scene_width << " - " << scene_fov << std::endl;
  std::cout << "\tNumber of lights : " << number_of_lights << std::endl;


  ScreenParam screen(screen_width, screen_height, screen_size);
  UserParam user(user_dist);
  SceneParam scene(scene_width, scene_fov);
  // read mesh from file
  FEVV::MeshSurface m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // retrieve point property map (aka geometry)
  auto pm = get(boost::vertex_point, m);
  // Note: the property maps must be extracted from the
  //       property maps bag, and explicitely passed as
  //       parameters to the filter, in order to make
  //       clear what property is used by the filter

  // retrieve or create vertex-color property map
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                  FEVV::MeshSurface >::pmap_type;
  VertexColorMap v_cm;
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "use existing vertex-color map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-color map" << std::endl;
    v_cm = make_property_map(FEVV::vertex_color, m);
    // store property map in property maps bag
    put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);
  }


  // retrieve or create vertex-normal property map
  using FaceNormalMap =
      typename FEVV::PMap_traits< FEVV::face_normal_t,
                                  FEVV::MeshSurface >::pmap_type;
  FaceNormalMap f_nm;
  if(has_map(pmaps_bag, FEVV::face_normal))
  {
    std::cout << "use existing face-normal map" << std::endl;
    f_nm = get_property_map(FEVV::face_normal, m, pmaps_bag);
  }
  else
  {
    std::cout << "create face-normal map" << std::endl;
    f_nm = make_property_map(FEVV::face_normal, m);
    // store property map in property maps bag
    put_property_map(FEVV::face_normal, m, pmaps_bag, f_nm);
    FEVV::Filters::calculate_face_normals(m, pm, f_nm);
  }


  // retrieve or create vertex-normal property map
  using VertexNormalMap =
      typename FEVV::PMap_traits< FEVV::vertex_normal_t,
                                  FEVV::MeshSurface >::pmap_type;
  VertexNormalMap v_nm;
  if(has_map(pmaps_bag, FEVV::vertex_normal))
  {
    std::cout << "use existing vertex-normal map" << std::endl;
    v_nm = get_property_map(FEVV::vertex_normal, m, pmaps_bag);
  }
  else
  {
    std::cout << "create vertex-normal map" << std::endl;
    v_nm = make_property_map(FEVV::vertex_normal, m);

    put_property_map(FEVV::vertex_normal, m, pmaps_bag, v_nm);
    FEVV::Filters::calculate_vertex_normals(m, pm, f_nm, v_nm);
  }

  using JndTypeMap =
      typename FEVV::Vertex_pmap_traits< FEVV::MeshSurface, double >::pmap_type;
  JndTypeMap jnd_m;

  jnd_m = FEVV::make_vertex_property_map< FEVV::MeshSurface, double >(m);


  // apply filter
  clock_t t_start = clock();
  FEVV::Filters::just_noticeable_distortion_filter(m,
                                                   pm,
                                                   v_nm,
                                                   f_nm,
                                                   jnd_m,
                                                   screen,
                                                   user,
                                                   scene,
                                                   number_of_lights,
                                                   true,
                                                   false);

  std::cout << "Time used to compute JND : "
            << (double)(clock() - t_start) / CLOCKS_PER_SEC << std::endl;


  double max_jnd, min_jnd;

  FEVV::Filters::ColorMeshLUT lut_courbure_clust =
      FEVV::Filters::make_LUT(false);

  FEVV::Filters::compute_min_max_vertices(m, jnd_m, min_jnd, max_jnd);

  FEVV::Filters::color_vertices_from_map(
      m, jnd_m, v_cm, min_jnd, max_jnd, lut_courbure_clust);

  // write mesh to file
  int counter = 0;
  BOOST_FOREACH(vertex_descriptor vd, m.vertices())
  {
    double jnd = get(jnd_m, vd);
    std::cout << "vertex #" << ++counter;
    std::cout << " - Jnd= " << jnd << std::endl;
  }
  FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);


  return 0;
}
