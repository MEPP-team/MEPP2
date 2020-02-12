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

#include <CGAL/boost/graph/properties.h>
#include <boost/foreach.hpp>
#include <limits>
#include <stdexcept> // for std::invalid_argument
#include <cmath>     // for std::fmod
#include <vector>

namespace FEVV {
namespace Filters {

// Numeric maps

typedef  std::vector< float >   ColorMeshLUT;

/**
 * \brief  Create a RGB LUT based on an HSV range.
 *         Default range creates a blue-cyan-green-yellow-red gradient.
 *
 * \param  color_in_0_255  if true the RGB values are in [0;255]
 *                         else they are in [0;1]
 * \param  colors_nbr  number of different colors in the LUT
 * \param  h1          1st bound of the H range
 * \param  h2          2nd bound of the H range
 *
 * \return  LUT as a 1D array of size 3*colors_nbr
 */
inline
ColorMeshLUT
make_LUT(bool color_in_0_255 = true,
         unsigned int colors_nbr = 256,
         float h1 = 240,
         float h2 = 0)
{
  // check parameters
  if(h1 == h2)
    throw std::invalid_argument("make_LUT: h1 and h2 must be different.");
  if(colors_nbr == 0)
    throw std::invalid_argument("make_LUT: colors_nbr must not be zero.");

  // create LUT
  ColorMeshLUT rgb_LUT(3*colors_nbr);

  // initializations
  float max_color_value = 1;
  if(color_in_0_255)
    max_color_value = 255;
  float step = (h2 - h1)/colors_nbr;
  float H = h1;
  const float S = 1;
  const float V = 1;

  // helper lambda function to convert HSV to RGB
  // see https://en.wikipedia.org/wiki/HSL_and_HSV#HSV_to_RGB
  // section "HSV to RGB alternative"
  auto f = [](float h, float s, float v, int n) -> float
  {
    float k = std::fmod(n + h/60.0f, 6.0f);
    return v - v*s*std::max(std::min({k, 4-k, 1.0f}), 0.0f);
  };

  // populate LUT
  for(unsigned int i = 0; i < colors_nbr; i++)
  {

    rgb_LUT[3 * i]     = max_color_value * f(H, S, V, 5); // R
    rgb_LUT[3 * i + 1] = max_color_value * f(H, S, V, 3); // G
    rgb_LUT[3 * i + 2] = max_color_value * f(H, S, V, 1); // B

    H += step;
  }

  return rgb_LUT;
}


/**
 * \brief  Set the color for the descriptor using a numerical property map and
 * the color array given
 *
 * \param  d            Descriptor used as key for the property map
 * \param  prop_map     Property map containing numerical values
 * \param  color_pmap   Property map containing color values
 * \param  minMetric    minimum value of the property map
 * \param  maxMetric    maximum value of the property map
 * \param  colors       1D array of size N*3
 *
 */
template< typename Descriptor,
          typename PropertyMap,
          typename ColorMap,
          typename MapType = typename PropertyMap::value_type >
void
color_descriptor_from_map(const Descriptor &d,
                          const PropertyMap &prop_map,
                          ColorMap &color_pmap,
                          const MapType min_metric,
                          const MapType max_metric,
                          const ColorMeshLUT &colors)
{
  typedef typename boost::property_traits< ColorMap >::value_type Color;
  size_t number_of_colors = colors.size() / 3;

  if(min_metric != max_metric)
  {
    // retrieve value from property map
    MapType val_metric = get(prop_map, d);

    // ensure value is between min and max to avoid LUT overflow
    val_metric = std::min(val_metric, max_metric);
    val_metric = std::max(val_metric, min_metric);

    // convert value to LUT id
    MapType id = (val_metric - min_metric) / (max_metric - min_metric);
    int indice_lut = static_cast< int >(std::floor((number_of_colors - 1) * id));

    // convert value to color
    Color newcolor(colors[3 * indice_lut],
                   colors[3 * indice_lut + 1],
                   colors[3 * indice_lut + 2]);
    put(color_pmap, d, newcolor);
  }
  else
  {
    put(color_pmap, d, Color(colors[0], colors[1], colors[2]));
  }
}


/**
 * \brief  Fill the color map for the faces of a mesh using a numerical property
 * map
 *
 * \param  g            mesh to colorize
 * \param  prop_map     property map containing numerical values
 * \param  color_pmap   property map containing color values
 * \param  minMetric    minimum value of the property map
 * \param  maxMetric    maximum value of the property map
 * \param  colors       1D array of size N*3
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename ColorMap,
          typename MapType = typename PropertyMap::value_type >
void
color_faces_from_map(const HalfedgeGraph &g,
                     const PropertyMap &prop_map,
                     ColorMap &color_pmap,
                     const MapType min_metric,
                     const MapType max_metric,
                     const ColorMeshLUT &colors = make_LUT())
{
  typedef typename boost::graph_traits< HalfedgeGraph >::face_descriptor
      Face_Descriptor;

  BOOST_FOREACH(Face_Descriptor f, faces(g))
  {
    color_descriptor_from_map(f,
                              prop_map,
                              color_pmap,
                              min_metric,
                              max_metric,
                              colors);
  }
}


/**
 * \brief  Fill the color map for the vertices of a mesh using a numerical
 * property map
 *
 * \param  g            mesh to colorize
 * \param  prop_map     property map containing numerical values
 * \param  color_pmap   property map containing color values
 * \param  minMetric    minimum value of the property map
 * \param  maxMetric    maximum value of the property map
 * \param  colors       1D array of size N*3
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename ColorMap,
          typename MapType = typename PropertyMap::value_type >
void
color_vertices_from_map(const HalfedgeGraph &g,
                        const PropertyMap &prop_map,
                        ColorMap &color_pmap,
                        const MapType min_metric,
                        const MapType max_metric,
                        const ColorMeshLUT &colors = make_LUT())
{
  typedef typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor
      Vertex_Descriptor;

  BOOST_FOREACH(Vertex_Descriptor v, vertices(g))
  {
    color_descriptor_from_map(v,
                              prop_map,
                              color_pmap,
                              min_metric,
                              max_metric,
                              colors);
  }
}


/**
 * \brief  Fill the color map for the halfedges of a mesh using a numerical
 * property map
 *
 * \param  g            mesh to colorize
 * \param  prop_map     property map containing numerical values
 * \param  color_pmap   property map containing color values
 * \param  minMetric    minimum value of the property map
 * \param  maxMetric    maximum value of the property map
 * \param  colors       1D array of size N*3
 *
 */
template< typename HalfedgeGraph,
          typename PropertyMap,
          typename ColorMap,
          typename MapType = typename PropertyMap::value_type >
void
color_halfedges_from_map(const HalfedgeGraph &g,
                         const PropertyMap &prop_map,
                         ColorMap &color_pmap,
                         const MapType min_metric,
                         const MapType max_metric,
                         const ColorMeshLUT &colors = make_LUT())
{
  typedef typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
      Halfedge_Descriptor;

  BOOST_FOREACH(Halfedge_Descriptor h, halfedges(g))
  {
    color_descriptor_from_map(h,
                              prop_map,
                              color_pmap,
                              min_metric,
                              max_metric,
                              colors);
  }
}

// Boolean map

/**
 * \brief  Set the color for the descriptor using a boolean property map and the
 * colors given
 *
 * \param  d                 Descriptor used as key for the property map
 * \param  prop_map			Property map containing boolean values
 * \param  color_pmap		Property map containing color values
 * \param  color1			color when the property is true
 * \param  color2			color when the property is false
 *
 */
template<
    typename Descriptor,
    typename BooleanMap,
    typename ColorMap,
    typename Color = typename boost::property_traits< ColorMap >::value_type >
void
color_descriptor_from_bool(const Descriptor &d,
                           const BooleanMap &prop_map,
                           ColorMap &color_pmap,
                           const Color &color1,
                           const Color &color2)
{
  if(get(prop_map, d))
  {
    put(color_pmap, d, color1);
  }
  else
  {
    put(color_pmap, d, color2);
  }
}


/**
 * \brief  Fill the color map for the faces of a mesh using a boolean property
 * map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing boolean values
 * \param  color_pmap		Property map containing color values
 * \param  color1			color when the property is true
 * \param  color2			color when the property is false
 *
 */
template<
    typename HalfedgeGraph,
    typename BooleanMap,
    typename ColorMap,
    typename Color = typename boost::property_traits< ColorMap >::value_type >
void
color_faces_from_bool_map(const HalfedgeGraph &g,
                          const BooleanMap &prop_map,
                          ColorMap &color_pmap,
                          const Color &color1,
                          const Color &color2)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::face_descriptor
      Face_Descriptor;

  BOOST_FOREACH(Face_Descriptor f, faces(g))
  {
    color_descriptor_from_bool(f, prop_map, color_pmap, color1, color2);
  }
}

/**
 * \brief  Fill the color map for the vertices of a mesh using a boolean
 * property map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing boolean values
 * \param  color_pmap		Property map containing color values
 * \param  color1			color when the property is true
 * \param  color2			color when the property is false
 *
 */
template<
    typename HalfedgeGraph,
    typename BooleanMap,
    typename ColorMap,
    typename Color = typename boost::property_traits< ColorMap >::value_type >
void
color_vertices_from_bool_map(const HalfedgeGraph &g,
                             const BooleanMap &prop_map,
                             ColorMap &color_pmap,
                             const Color &color1,
                             const Color &color2)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::vertex_descriptor
      Vertex_Descriptor;

  BOOST_FOREACH(Vertex_Descriptor v, vertices(g))
  {
    color_descriptor_from_bool(v, prop_map, color_pmap, color1, color2);
  }
}

/**
 * \brief  Fill the color map for the halfedges of a mesh using a boolean
 * property map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing boolean values
 * \param  color_pmap		Property map containing color values
 * \param  color1			color when the property is true
 * \param  color2			color when the property is false
 *
 */
template<
    typename HalfedgeGraph,
    typename BooleanMap,
    typename ColorMap,
    typename Color = typename boost::property_traits< ColorMap >::value_type >
void
color_halfedge_from_bool_map(const HalfedgeGraph &g,
                             const BooleanMap &prop_map,
                             ColorMap &color_pmap,
                             const Color &color1,
                             const Color &color2)
{
  typedef typename boost::graph_traits< HalfedgeGraph >::halfedge_descriptor
      Halfedge_Descriptor;

  BOOST_FOREACH(Halfedge_Descriptor h, halfedges(g))
  {
    color_descriptor_from_bool(h, prop_map, color_pmap, color1, color2);
  }
}
} // namespace Filters
} // namespace FEVV
