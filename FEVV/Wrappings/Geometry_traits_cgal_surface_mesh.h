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

