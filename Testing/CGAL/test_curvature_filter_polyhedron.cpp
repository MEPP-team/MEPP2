#include <fstream>
#include <string>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include "FEVV/Wrappings/Geometry_traits_cgal_polyhedron_3.h" // Geometry adaptor
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include "FEVV/Filters/Generic/Manifold/Curvature/curvature.hpp" // for calculate_curvature()
#include "FEVV/Filters/Generic/calculate_face_normals.hpp"

// For datastructure default property map:
#include "FEVV/Wrappings/properties_polyhedron_3.h"

typedef CGAL::Cartesian< double > Kernel;
// typedef Kernel::Vector_3 Vector;
typedef CGAL::Polyhedron_3< Kernel, CGAL::Polyhedron_items_with_id_3 >
    Polyhedron;

//------------------------------------------------------------------------------

// code from Visualization/PluginFilters/curvature/CurvaturePlugin.h
void
curvature(Polyhedron *mesh)
{
  typedef Polyhedron HalfedgeGraph;

  std::cout << "Asking to Curvature mesh ! " << std::endl;

  auto pm = get(boost::vertex_point, *mesh);

  // ---

  using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;

  // Face normal map
  // this one is a standard property map
  auto f_nm = make_property_map(FEVV::face_normal, *mesh);

  FEVV::Filters::calculate_face_normals(*mesh, pm, f_nm);

  // ---

  // Vertex curvature map
  // this one is a NON-standard property map!
  auto v_cm =
      FEVV::make_vertex_property_map< HalfedgeGraph,
                                      FEVV::Filters::v_Curv< HalfedgeGraph > >(
          *mesh);

  /*! \brief Minimum and maximum values of the minimum and maximum curvature
   * fields (usefull for color rendering)*/
  double min_nrm_min_curvature, max_nrm_min_curvature, min_nrm_max_curvature,
      max_nrm_max_curvature;

  bool value_is_geod = false;
  double value_radius = 0.001;
  FEVV::Filters::calculate_curvature( // B) call the filter corresponding to
                                      // your operation
      *mesh,
      v_cm,
      pm,
      f_nm,
      value_is_geod,
      value_radius,
      min_nrm_min_curvature,
      max_nrm_min_curvature,
      min_nrm_max_curvature,
      max_nrm_max_curvature);

  // ---

  std::cout << "Curvature mesh, isGeod: " << value_is_geod
            << " - radius: " << value_radius << "." << std::endl;
}

//------------------------------------------------------------------------------

int
main(int narg, char **argv)
{
  if(narg < 2)
  {
    std::cout << "Usage: " << argv[0] << " mesh_filename" << std::endl;
    exit(EXIT_FAILURE);
  }
  std::string filename(argv[1]);

  std::ifstream in(filename);
  if(!in)
  {
    std::cout << "Unable to read file " << filename << std::endl;
    exit(EXIT_FAILURE);
  }

  Polyhedron p;
  in >> p;


  curvature(&p);


  std::cout << "Done." << std::endl;

  return 0;
}
