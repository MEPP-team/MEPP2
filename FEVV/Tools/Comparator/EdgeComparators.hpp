#ifndef __EdgeComparators_hxx
#define __EdgeComparators_hxx

#include <limits>
#include <iterator>

namespace FEVV { 
namespace Comparator 
{
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class EdgeComparator
  {
  public:
    typedef boost::graph_traits<Graph>     GraphTraits;
    typedef typename GraphTraits::vertex_descriptor           vertex_descriptor;
    typedef typename GraphTraits::edge_descriptor             edge_descriptor;
    typedef FEVV::Geometry_traits<Graph>   Geometry;
    typedef typename Geometry::Scalar                         Scalar;
    typedef typename Geometry::Point                          Point;
  private:
    const Graph& _g;
    const PointMap& _pm;
    const GeometryTraits _gt;
  public:
    EdgeComparator(const Graph& g, const PointMap& pm) :_g(g), _pm(pm), _gt(GeometryTraits(g)) {}
    EdgeComparator(const Graph& g, const PointMap& pm, const GeometryTraits& gt):_g(g), _pm(pm), _gt(gt) {}
    bool operator()(edge_descriptor e1, edge_descriptor e2)
    {
      Point pe1_pv1 = get(_pm, source(e1, _g));
      Point pe1_pv2 = get(_pm, target(e1, _g));

      Point pe2_pv1 = get(_pm, source(e2, _g));
      Point pe2_pv2 = get(_pm, target(e2, _g));
#if 1
	  Point minie1(0,0,0), maxie1(0,0,0), minie2(0,0,0), maxie2(0,0,0);

	  if( _gt.get_x(pe1_pv1) < _gt.get_x(pe1_pv2) )
	  {
		  minie1 = pe1_pv1;
		  maxie1 = pe1_pv2;
	  }
	  else if( _gt.get_x(pe1_pv1) > _gt.get_x(pe1_pv2) )
	  {
		  maxie1 = pe1_pv1;
		  minie1 = pe1_pv2;
	  }
	  else if( _gt.get_y(pe1_pv1) < _gt.get_y(pe1_pv2) )
	  {
		  minie1 = pe1_pv1;
		  maxie1 = pe1_pv2;
	  }
	  else if( _gt.get_y(pe1_pv1) > _gt.get_y(pe1_pv2) )
	  {
		  maxie1 = pe1_pv1;
		  minie1 = pe1_pv2;
	  }
	  else if( _gt.get_z(pe1_pv1) < _gt.get_z(pe1_pv2) )
	  {
		  minie1 = pe1_pv1;
		  maxie1 = pe1_pv2;
	  }
	  else if( _gt.get_z(pe1_pv1) > _gt.get_z(pe1_pv2) )
	  {
		  maxie1 = pe1_pv1;
		  minie1 = pe1_pv2;
	  }
	  else minie1 = maxie1 = pe1_pv1;
	  /////////////////////////////////////////////////////////////////////////
	  if( _gt.get_x(pe2_pv1) < _gt.get_x(pe2_pv2) )
	  {
		  minie2 = pe2_pv1;
		  maxie2 = pe2_pv2;
	  }
	  else if( _gt.get_x(pe2_pv1) > _gt.get_x(pe2_pv2) )
	  {
		  maxie2 = pe2_pv1;
		  minie2 = pe2_pv2;
	  }
	  else if( _gt.get_y(pe2_pv1) < _gt.get_y(pe2_pv2) )
	  {
		  minie2 = pe2_pv1;
		  maxie2 = pe2_pv2;
	  }
	  else if( _gt.get_y(pe2_pv1) > _gt.get_y(pe2_pv2) )
	  {
		  maxie2 = pe2_pv1;
		  minie2 = pe2_pv2;
	  }
	  else if( _gt.get_z(pe2_pv1) < _gt.get_z(pe2_pv2) )
	  {
		  minie2 = pe2_pv1;
		  maxie2 = pe2_pv2;
	  }
	  else if( _gt.get_z(pe2_pv1) > _gt.get_z(pe2_pv2) )
	  {
		  maxie2 = pe2_pv1;
		  minie2 = pe2_pv2;
	  }
	  else minie2 = maxie2 = pe2_pv1;
	  /////////////////////////////////////////////////////////////////////////
	  if (fabs(_gt.get_x(minie1) - _gt.get_x(minie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(minie1) < _gt.get_x(minie2);
      else if (fabs(_gt.get_y(minie1) - _gt.get_y(minie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(minie1) < _gt.get_y(minie2);
      else if (fabs(_gt.get_z(minie1) - _gt.get_z(minie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_z(minie1) < _gt.get_z(minie2);
	  else if (fabs(_gt.get_x(maxie1) - _gt.get_x(maxie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(maxie1) < _gt.get_x(maxie2);
      else if (fabs(_gt.get_y(maxie1) - _gt.get_y(maxie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(maxie1) < _gt.get_y(maxie2);
      else if (fabs(_gt.get_z(maxie1) - _gt.get_z(maxie2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_z(maxie1) < _gt.get_z(maxie2);
	  else
	  {
		 auto e1_v1 = in_edges(source(e1, _g),_g), 
		      e1_v2 = in_edges(target(e1, _g),_g),
			  e2_v1 = in_edges(source(e2, _g),_g),
			  e2_v2 = in_edges(target(e2, _g),_g);
	     
		 auto e1_v1_deg = std::distance(e1_v1.first, e1_v1.second),
              e1_v2_deg = std::distance(e1_v2.first, e1_v2.second),
              e2_v1_deg = std::distance(e2_v1.first, e2_v1.second),
              e2_v2_deg = std::distance(e2_v2.first, e2_v2.second);		
         if( minie1 == pe1_pv1 )
         {
			 if( minie2 == pe2_pv1 )
				 if( e1_v1_deg!= e2_v1_deg)
			       return e1_v1_deg < e2_v1_deg;
			 else 
				 if( e1_v1_deg!= e2_v2_deg)
				   return e1_v1_deg < e2_v2_deg;
		 }	
         else 
         {
			 if( minie2 == pe2_pv1 )
				 if( e1_v2_deg!= e2_v1_deg)				 
				   return e1_v2_deg < e2_v1_deg;
			 else
				 if( e1_v2_deg!= e2_v2_deg)
				   return e1_v2_deg < e2_v2_deg; 
		 }		 
		 
         if( maxie1 == pe1_pv1 )
         {
			 if( maxie2 == pe2_pv1 )
				 if( e1_v1_deg!= e2_v1_deg)
			       return e1_v1_deg < e2_v1_deg;
			 else 
				 if( e1_v1_deg!= e2_v2_deg)
				   return e1_v1_deg < e2_v2_deg;
		 }	
         else 
         {
			 if( maxie2 == pe2_pv1 )
				 if( e1_v2_deg!= e2_v1_deg)				 
				   return e1_v2_deg < e2_v1_deg;
			 else
				 if( e1_v2_deg!= e2_v2_deg)
				   return e1_v2_deg < e2_v2_deg; 
		 }			 
	  }
	  /////////////////////////////////////////////////////////////////////////
#else	  
      Point pv1((_gt.get_x(pe1_pv1) + _gt.get_x(pe1_pv2))*0.5,
                (_gt.get_y(pe1_pv1) + _gt.get_y(pe1_pv2))*0.5,
                (_gt.get_z(pe1_pv1) + _gt.get_z(pe1_pv2))*0.5);

      Point pv2((_gt.get_x(pe2_pv1) + _gt.get_x(pe2_pv2))*0.5,
                (_gt.get_y(pe2_pv1) + _gt.get_y(pe2_pv2))*0.5,
                (_gt.get_z(pe2_pv1) + _gt.get_z(pe2_pv2))*0.5);

	  if (fabs(_gt.get_x(pv1) - _gt.get_x(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_x(pv1) < _gt.get_x(pv2);
      else if (fabs(_gt.get_y(pv1) - _gt.get_y(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_y(pv1) < _gt.get_y(pv2);
      else if (fabs(_gt.get_z(pv1) - _gt.get_z(pv2)) > std::numeric_limits<Scalar>::epsilon())
        return _gt.get_z(pv1) < _gt.get_z(pv2);
#endif	  

      return false; // can be an issue
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, typename PointMap, typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static
  EdgeComparator<Graph, PointMap, GeometryTraits>
    get_edge_comparator(const Graph& g, PointMap pm) { return EdgeComparator<Graph, PointMap, GeometryTraits>(g, pm); }
} // namespace Container
} // namespace FEVV

#endif /* __EdgeComparators_hxx */
