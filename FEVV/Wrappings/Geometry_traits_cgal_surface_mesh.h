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

#include <CGAL/Surface_mesh.h>
#include <CGAL/Kernel_traits.h>
#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_cartesian.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_exact_predicates_inexact_constructions_kernel.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \anchor  RetrieveKernel_specialization_SurfaceMesh
 * \brief   OpenMesh specialization. Refer to \ref RetrieveKernel for the
 *          documentation of the generic class.
 */
template< typename PointT >
struct RetrieveKernel< CGAL::Surface_mesh< PointT > >
{
  typedef typename CGAL::Kernel_traits< PointT >::Kernel Kernel;
};

} // namespace FEVV

