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

#include <limits>
#include <iterator>
#include <cmath>

namespace FEVV { 
namespace Comparator 
{
  template <typename Graph, 
            typename PointMap, 
            typename EdgeWeightMap = typename FEVV::Edge_pmap_traits<Graph,
                                                                     typename FEVV::Geometry_traits<Graph>::Scalar
                                                                     >::pmap_type,
            typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  class EdgeComparator
  {
  public:
    typedef boost::graph_traits<Graph>              GraphTraits;
    typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
    typedef typename GraphTraits::edge_descriptor   edge_descriptor;
    typedef FEVV::Geometry_traits<Graph>            Geometry;
    typedef typename Geometry::Scalar               Scalar;
    typedef typename Geometry::Point                Point;
  private:
    const Graph& _g;
    const PointMap& _pm;
    EdgeWeightMap _ew;
    const GeometryTraits _gt;
  public:
    EdgeComparator(const Graph& g, const PointMap& pm) :_g(g), _pm(pm), _ew(0, get(boost::edge_index, g)), _gt(GeometryTraits(g)) {}
    EdgeComparator(const Graph& g, const PointMap& pm, const GeometryTraits& gt):_g(g), _pm(pm), _ew(0, get(boost::edge_index, g)), _gt(gt) {}
    EdgeComparator(const Graph& g, const PointMap& pm, const EdgeWeightMap& ew) :_g(g), _pm(pm), _ew(ew), _gt(GeometryTraits(g)) {}
    EdgeComparator(const Graph& g, const PointMap& pm, const EdgeWeightMap& ew, const GeometryTraits& gt) :_g(g), _pm(pm), _ew(ew), _gt(gt) {}
    EdgeComparator(const EdgeComparator& other): _g(other._g), _pm(other._pm), _ew(other._ew), _gt(other._gt) {}
    bool operator()(edge_descriptor e1, edge_descriptor e2)
    {
      if (_ew.storage_begin()!= _ew.storage_end())
      {
        Scalar val_e1 = get(_ew, e1), val_e2 = get(_ew, e2);

        if (val_e1 < val_e2)
          return true;
        else if (val_e1 > val_e2)
          return false;
      }
      Point pe1_pv1 = get(_pm, source(e1, _g));
      Point pe1_pv2 = get(_pm, target(e1, _g));

      Point pe2_pv1 = get(_pm, source(e2, _g));
      Point pe2_pv2 = get(_pm, target(e2, _g));

	  Point minie1(0,0,0), maxie1(0,0,0), minie2(0,0,0), maxie2(0,0,0);
      bool minie1_is_pe1_pv1 = false, 
	       minie2_is_pe2_pv1 = false;
	  if( _gt.get_x(pe1_pv1) < _gt.get_x(pe1_pv2) )
	  {
		  minie1 = pe1_pv1; minie1_is_pe1_pv1 = true;
		  maxie1 = pe1_pv2;
	  }
	  else if( _gt.get_x(pe1_pv1) > _gt.get_x(pe1_pv2) )
	  {
		  maxie1 = pe1_pv1;
		  minie1 = pe1_pv2;
	  }
	  else if( _gt.get_y(pe1_pv1) < _gt.get_y(pe1_pv2) )
	  {
		  minie1 = pe1_pv1; minie1_is_pe1_pv1 = true;
		  maxie1 = pe1_pv2;
	  }
	  else if( _gt.get_y(pe1_pv1) > _gt.get_y(pe1_pv2) )
	  {
		  maxie1 = pe1_pv1;
		  minie1 = pe1_pv2;
	  }
	  else if( _gt.get_z(pe1_pv1) < _gt.get_z(pe1_pv2) )
	  {
		  minie1 = pe1_pv1; minie1_is_pe1_pv1 = true;
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
		  minie2 = pe2_pv1; minie2_is_pe2_pv1 = true;
		  maxie2 = pe2_pv2;
	  }
	  else if( _gt.get_x(pe2_pv1) > _gt.get_x(pe2_pv2) )
	  {
		  maxie2 = pe2_pv1;
		  minie2 = pe2_pv2;
	  }
	  else if( _gt.get_y(pe2_pv1) < _gt.get_y(pe2_pv2) )
	  {
		  minie2 = pe2_pv1; minie2_is_pe2_pv1 = true;
		  maxie2 = pe2_pv2;
	  }
	  else if( _gt.get_y(pe2_pv1) > _gt.get_y(pe2_pv2) )
	  {
		  maxie2 = pe2_pv1;
		  minie2 = pe2_pv2;
	  }
	  else if( _gt.get_z(pe2_pv1) < _gt.get_z(pe2_pv2) )
	  {
		  minie2 = pe2_pv1; minie2_is_pe2_pv1 = true;
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
         if( minie1_is_pe1_pv1 )
         {
			 if( minie2_is_pe2_pv1 )
             {
				 if( e1_v1_deg!= e2_v1_deg)
			       return e1_v1_deg < e2_v1_deg;
             }
			 else
             {				 
				 if( e1_v1_deg!= e2_v2_deg)
				   return e1_v1_deg < e2_v2_deg;
             }
		 }	
         else 
         {
			 if( minie2_is_pe2_pv1 )
             {
				 if( e1_v2_deg!= e2_v1_deg)				 
				   return e1_v2_deg < e2_v1_deg;
             }
			 else
             {				 
				 if( e1_v2_deg!= e2_v2_deg)
				   return e1_v2_deg < e2_v2_deg; 
             }
		 }		 
		 
         if( !minie1_is_pe1_pv1 )
         {
			 if( !minie2_is_pe2_pv1 )
             {
				 if( e1_v1_deg!= e2_v1_deg)
			       return e1_v1_deg < e2_v1_deg;
             }
			 else 
             {
				 if( e1_v1_deg!= e2_v2_deg)
				   return e1_v1_deg < e2_v2_deg;
             }
		 }	
         else 
         {
			 if( !minie2_is_pe2_pv1 )
             {
				 if( e1_v2_deg!= e2_v1_deg)				 
				   return e1_v2_deg < e2_v1_deg;
             }
			 else
             {
				 if( e1_v2_deg!= e2_v2_deg)
				   return e1_v2_deg < e2_v2_deg; 
             }
		 }			 
	  }

      return false; // can be an issue
    }
  };
  /////////////////////////////////////////////////////////////////////////////
  template <typename Graph, 
            typename PointMap, 
            typename EdgeWeightMap = typename FEVV::Edge_pmap_traits<Graph,
                                                                     typename FEVV::Geometry_traits<Graph>::Scalar
                                                                     >::pmap_type, 
            typename GeometryTraits = FEVV::Geometry_traits<Graph> >
  static
  EdgeComparator<Graph, PointMap, EdgeWeightMap, GeometryTraits>
    get_edge_comparator(const Graph& g, PointMap pm) { return EdgeComparator<Graph, PointMap, EdgeWeightMap, GeometryTraits>(g, pm); }
} // namespace Container
} // namespace FEVV

