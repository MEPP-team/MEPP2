#ifndef TRANSLATION_WITH_OBJECT_H
#define TRANSLATION_WITH_OBJECT_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "FEVV/Wrappings/Geometry_traits.h"

namespace FEVV {
namespace Filters {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// \class	TranslationFilter
///
/// \brief	Translation filter (object version)
///
////////////////////////////////////////////////////////////////////////////////////////////////////
template< typename HalfedgeGraph,
          typename PointMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
class TranslationFilter
{
public:
  /**
   * \brief  Constructor
   *
   * \param	ag	 the mesh
   * \param	apm	 the point map
   */
  TranslationFilter(const HalfedgeGraph &ag, PointMap &apm) : g(ag), pm(apm) {}

  /**
   * \brief  Translate a mesh
   *
   * \param  offsetX           x translation offset
   * \param  offsetY           y translation offset
   * \param  offsetZ           z translation offset
   * \param  gt                the geometry traits to use
   *
   * \sa      the simplified variant that use the default geometry traits
   *          of the mesh.
   *
   * \ingroup  GenericManifoldFilters
   */
  void translate(typename GeometryTraits::Scalar offset_x,
                 typename GeometryTraits::Scalar offset_y,
                 typename GeometryTraits::Scalar offset_z,
                 const GeometryTraits &gt)
  {
    typedef boost::graph_traits< HalfedgeGraph > GraphTraits;
    typedef typename GraphTraits::vertex_iterator vertex_iterator;
    typedef typename boost::property_traits< PointMap >::value_type Point;
    typedef typename boost::property_traits< PointMap >::reference Reference;

    vertex_iterator vi;
    for(vi = vertices(g).first; vi != vertices(g).second; ++vi)
    {
      Point p = pm[*vi];
      put(pm,
          *vi,
          Point(gt.get_x(p) + offset_x,
                gt.get_y(p) + offset_y,
                gt.get_z(p) + offset_z));
    }
  }

  /**
   * \brief  Translate a mesh
   *
   * \param  offsetX           x translation offset
   * \param  offsetY           y translation offset
   * \param  offsetZ           z translation offset
   *
   * \sa      the variant that use the geometry traits provided by the user.
   *
   * \ingroup  GenericManifoldFilters
   */
  void translate(typename GeometryTraits::Scalar offset_x,
                 typename GeometryTraits::Scalar offset_y,
                 typename GeometryTraits::Scalar offset_z)

  {
    GeometryTraits gt(g);
    translate(offset_x, offset_y, offset_z, gt);
  }

protected:
  const HalfedgeGraph &g; ///< the mesh
  PointMap &pm;           ///< the point map
};

} // namespace Filters
} // namespace FEVV

#endif // TRANSLATION_WITH_OBJECT_H
