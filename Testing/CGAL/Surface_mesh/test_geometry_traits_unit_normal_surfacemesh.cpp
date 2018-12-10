#include <iostream>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_surface_mesh.h"

int
main(int narg, char **argv)
{
  if(narg > 1)
  {
    std::cout << "Usage: " << argv[0] << " (no arguments)." << std::endl;
    exit(EXIT_FAILURE);
  }

  typedef CGAL::Cartesian< double > Kernel;
  typedef CGAL::Surface_mesh< Kernel::Point_3 > Mesh;

  typedef FEVV::Geometry_traits< Mesh > Geometry;
  typedef Geometry::Point Point;
  typedef Geometry::Vector Vector;

  Mesh m;
  Geometry g(m);
  Point p1(0.0, 0.0, 0.0);
  Point p2(1.0, 0.0, 0.0);
  Point p3(0.0, 0.1, 0.0);
  Vector n = g.unit_normal(p1, p2, p3);

  if(n == Vector(0.0, 0.0, 1.0))
  {
    std::cout << "OK." << std::endl;
    return 0;
  }
  std::cout << "Result is " << n << " when (0, 0, 1) was expected."
            << std::endl;
  return 1;
}
