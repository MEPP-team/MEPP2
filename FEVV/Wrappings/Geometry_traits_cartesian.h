#pragma once

#include <CGAL/Cartesian.h>
#include "FEVV/Wrappings/Geometry_traits_cgal_common.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \anchor  Geometry_traits_specialization_cgal_cartesian
 * \brief   OpenMesh specialization of the \ref Geometry_traits generic class.
 *          For usage refer to
 *          \link Geometry_traits_documentation_dummy
 *                Geometry traits documentation
 *          \endlink
 */
template< typename Mesh, typename T1 >
class Geometry_traits< Mesh, CGAL::Cartesian< T1 > >
    : public Geometry_traits_for_cgal< Mesh, CGAL::Cartesian< T1 > >
{
  typedef Geometry_traits_for_cgal< Mesh, CGAL::Cartesian< T1 > > Base;

public:
  Geometry_traits(const Mesh &m) : Base(m) {}
};

} // namespace FEVV

