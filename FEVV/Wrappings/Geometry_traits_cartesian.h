// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
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

