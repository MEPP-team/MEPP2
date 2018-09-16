#ifndef FEVV_WRAPPINGS_H
#define FEVV_WRAPPINGS_H

#include <CGAL/boost/graph/properties.h>

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/Wrappings_cgal_polyhedron_3.h"
#include "FEVV/Wrappings/Wrappings_cgal_surface_mesh.h"
#include "FEVV/Wrappings/Wrappings_cgal_linear_cell_complex.h"
#endif

#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/Wrappings_openmesh.h"
#endif

#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/Wrappings_aif.h"
#endif


#endif // FEVV_WRAPPINGS_H
