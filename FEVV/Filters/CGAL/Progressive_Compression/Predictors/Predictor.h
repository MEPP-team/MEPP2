#pragma once

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <CGAL/boost/graph/iterator.h>

#include "FEVV/Wrappings/Geometry_traits.h"

#include "FEVV/Filters/CGAL/Progressive_Compression/Operators/KeptPosition.h"
#include "FEVV/Filters/CGAL/Progressive_Compression/Compression/CollapseInfo.h"

#include <vector>
#include <utility>

namespace FEVV {
namespace Filters {

/** \brief Abstract class used to predict position. 
 */
template<
    typename HalfedgeGraph,
    typename PointMap>
class Predictor
{
public:
    using vertex_descriptor =
      typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor;
    using halfedge_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor;
    using face_descriptor =
            typename boost::graph_traits< HalfedgeGraph >::face_descriptor;
    using Vector = typename FEVV::Geometry_traits< HalfedgeGraph >::Vector;
    using Point = typename FEVV::Geometry_traits< HalfedgeGraph >::Point;
    using Geometry = typename FEVV::Geometry_traits< HalfedgeGraph >;
  Predictor(
      HalfedgeGraph &g,
      KeptPosition< HalfedgeGraph, PointMap >
          *kp,
      PointMap &pm)
      : _g(g), _gt(Geometry(_g)), _pm(pm)
  {
    _kp = kp;
    _kept_position = Point(0, 0, 0);
  }
  virtual ~Predictor(){}

  /// Compression side: Computes geometric residuals from a set of info about 
  /// the collapse
  virtual std::vector< Vector >
  ComputeResiduals(CollapseInfo< HalfedgeGraph, PointMap > &mem) = 0;
  
  /// Decompression side: predicts a position from encoded residuals
  virtual std::pair< Point, Point >
  PlacePoints(const std::vector< Vector > &residuals,
              vertex_descriptor vkept, /// should be target (h1) and target(h2)
              halfedge_descriptor h1, /// first halfedge to expand into a face
              halfedge_descriptor h2, /// second halfedge to expand into a face
              bool is_reverse /// is the current split a reverse case?
              ) = 0;

  virtual FEVV::Filters::PREDICTION_TYPE getType() const { return _type; }
  virtual const std::tuple< bool, bool, bool, bool >& getInfoMidPoint() const = 0;
  
  virtual void set_bit_info(bool /*b1*/, bool /*b2*/, bool /*b3*/, bool /*b4*/) = 0;
  
  virtual void set_rev(bool /*b*/) = 0;
  
  virtual const Point& get_kept_position() const { return _kept_position; }

  virtual FEVV::Filters::VKEPT_POSITION getTypeKP() const { return _kp->getType(); }


  virtual std::string getMethodasString() const = 0;

protected:
  KeptPosition< HalfedgeGraph, PointMap >
      *_kp;
  HalfedgeGraph &_g;
  const Geometry _gt;
  PointMap &_pm;
  int _nbResiduals;
  std::vector< Vector > _residuals;
  Point _kept_position;
  FEVV::Filters::PREDICTION_TYPE _type;
};
} // namespace Filters
} // namespace FEVV
