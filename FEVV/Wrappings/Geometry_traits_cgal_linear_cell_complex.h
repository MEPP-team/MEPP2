#ifndef GEOMETRY_TRAITS_CGAL_LINEAR_CELL_COMPLEX_H
#define GEOMETRY_TRAITS_CGAL_LINEAR_CELL_COMPLEX_H

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/Geometry_traits_cartesian.h"
#include "FEVV/Wrappings/Geometry_traits_cgal_exact_predicates_inexact_constructions_kernel.h"
#include "FEVV/Wrappings/Geometry_traits_operators.h"

namespace FEVV {

/**
 * \ingroup Geometry_traits_group
 * \anchor  RetrieveKernel_specialization_cgal_linear_cell_complex
 * \brief   Linear cell complex specialization. Refer to \ref RetrieveKernel for
 * the documentation of the generic class.
 */
template< unsigned int D,
          unsigned int AmbientDim,
          class Traits,
          class Items,
          class Alloc,
          template< unsigned int, class, class, class, class > class CMap,
          class Storage >
struct RetrieveKernel<
    CGAL::Linear_cell_complex_for_combinatorial_map< D,
                                                     AmbientDim,
                                                     Traits,
                                                     Items,
                                                     Alloc,
                                                     CMap,
                                                     Storage > >
{
  typedef typename CGAL::Linear_cell_complex_for_combinatorial_map<
      D,
      AmbientDim,
      Traits,
      Items,
      Alloc,
      CMap,
      Storage >::Traits::Kernel Kernel;
};

} // namespace FEVV

#endif // GEOMETRY_TRAITS_CGAL_LINEAR_CELL_COMPLEX_H
