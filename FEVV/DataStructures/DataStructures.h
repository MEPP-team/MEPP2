#ifndef FEVV_DATASTRUCTURES_H
#define FEVV_DATASTRUCTURES_H


#ifdef FEVV_USE_CGAL
#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"
#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"
#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h"
#endif

#ifdef FEVV_USE_OPENMESH
#include "FEVV/DataStructures/DataStructures_openmesh.h"
#endif

#ifdef FEVV_USE_AIF
#include "FEVV/DataStructures/DataStructures_aif.h"
#endif

#endif // FEVV_DATASTRUCTURES_H
