#ifndef GEOMETRY_TRAITS_CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
#define GEOMETRY_TRAITS_CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_common.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_specialization_cgal_exact_predicates_inexact
 * \brief   OpenMesh specialization of the \ref Geometry_traits generic class.
 *          For usage refer to
 *          \link Geometry_traits_documentation_dummy
 *                Geometry traits documentation
 *          \endlink
 */
template< typename Mesh >
class Geometry_traits< Mesh,
                       CGAL::Exact_predicates_inexact_constructions_kernel >
    : public Geometry_traits_for_cgal<
          Mesh,
          CGAL::Exact_predicates_inexact_constructions_kernel >
{
  typedef Geometry_traits_for_cgal<
      Mesh,
      CGAL::Exact_predicates_inexact_constructions_kernel >
      Base;

public:
  Geometry_traits(const Mesh &m) : Base(m) {}
};

}; // namespace FEVV

#endif // GEOMETRY_TRAITS_CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
