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
#pragma once

#ifdef FEVV_USE_CGAL
#include "FEVV/DataStructures/DataStructures_cgal_polyhedron_3.h"
#include "FEVV/DataStructures/DataStructures_cgal_surface_mesh.h"
#include "FEVV/DataStructures/DataStructures_cgal_linear_cell_complex.h"
#include "FEVV/DataStructures/DataStructures_cgal_point_set.h"
#endif

#ifdef FEVV_USE_OPENMESH
#include "FEVV/DataStructures/DataStructures_openmesh.h"
#endif

#ifdef FEVV_USE_AIF
#include "FEVV/DataStructures/DataStructures_aif.h"
#endif

#ifdef FEVV_USE_PCL
#include "FEVV/DataStructures/DataStructures_pcl_point_cloud.h"
#endif
