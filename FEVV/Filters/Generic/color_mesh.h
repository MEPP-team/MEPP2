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
#include <memory>    // for std::unique_ptr
#include <stdexcept> // for std::invalid_argument
#include <cmath>     // for std::fmod

namespace FEVV {
namespace Filters {

// Numeric maps


/**
 * \brief  Create a RGB LUT based on an HSV range.
 *         Default range create a cyan-green-yellow-red gradient.
 *
 * \param  h1  1st bound of the H range
 * \param  h2  2nd bound of the H range
 * \param  colors_nbr  number of different colors of the LUT
 * \param  color_in_0_255  if true the RGB values are in [0;255]
 *                         else they are in [0;1]
 */
inline
std::unique_ptr< double[] >
make_LUT(float h1 = 180,
         float h2 = 0,
         unsigned int colors_nbr = 255,
         bool color_in_0_255 = true)
{
  // check parameters
  if(h1 == h2)
    throw std::invalid_argument("make_LUT: h1 and h2 must be different.");
  if(colors_nbr == 0)
    throw std::invalid_argument("make_LUT: colors_nbr must not be zero.");

  // create LUT
  std::unique_ptr< double[] > rgb_LUT(new double[3 * colors_nbr]);

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
    float k = std::fmod(n + h/60.0, 6);
    return v - v*s*std::max(std::min({k, 4-k, 1.0f}), 0.0f);
  };

  // populate LUT
  for(unsigned int i = 0; i < colors_nbr; i++)
  {
    std::cout << "  H=" << H; //TODO-elo-rm
    rgb_LUT[3 * i]     = max_color_value * f(H, S, V, 5); // R
    std::cout << "  R=" << rgb_LUT[3 * i]; //TODO-elo-rm
    rgb_LUT[3 * i + 1] = max_color_value * f(H, S, V, 3); // G
    std::cout << "  G=" << rgb_LUT[3 * i + 1]; //TODO-elo-rm
    rgb_LUT[3 * i + 2] = max_color_value * f(H, S, V, 1); // B
    std::cout << "  B=" << rgb_LUT[3 * i + 2] << std::endl; //TODO-elo-rm

    H += step;
  }

  return rgb_LUT;
}


/**
 * \brief  Set the color for the descriptor using a numerical property map and
 * the color array given
 *
 * \param  d                 Descriptor used as key for the property map
 * \param  prop_map			Property map containing numerical values
 * \param  color_pmap		Property map containing color values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 * \param  colors			array of dimension N*3
 * \param  number_of_colors	number N of colors in the array
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
                          double *colors,
                          int number_of_colors)
{
  typedef typename boost::property_traits< ColorMap >::value_type Color;

  if(min_metric != max_metric)
  {
    // retrieve value from property map
    MapType val_metric = get(prop_map, d);
    std::cout << "val_metric=" << val_metric << std::endl; //TODO-elo-rm

    // ensure value is between min and max to avoid LUT overflow
    val_metric = std::min(val_metric, max_metric);
    val_metric = std::max(val_metric, min_metric);

    // convert value to LUT id
    MapType id = (val_metric - min_metric) / (max_metric - min_metric);
    int indice_lut = static_cast< int >(std::floor(number_of_colors * id));
    std::cout << "val_metric=" << val_metric << "  min_metric=" << min_metric << "  max_metric=" << max_metric << "  id=" << id << "  indice_lut = " << indice_lut << std::endl; //TODO-elo-rm

    // convert value to color
    Color newcolor(colors[3 * indice_lut],
                   colors[3 * indice_lut + 1],
                   colors[3 * indice_lut + 2]);
    put(color_pmap, d, newcolor);
  }
  else
  {
    put(color_pmap, d, Color(0.5, 0.5, 0.5));
  }
}


/**
 * \brief  Fill the color map for the faces of a mesh using a numerical property
 * map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing numerical values
 * \param  color_pmap		Property map containing color values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 * \param  colors			array of dimension N*3
 * \param  number_of_colors	number N of colors in the array
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
                     double *colors = make_LUT().get(),
                     int number_of_colors = 255)
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
                              colors,
                              number_of_colors);
  }
}


/**
 * \brief  Fill the color map for the vertices of a mesh using a numerical
 * property map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing numerical values
 * \param  color_pmap		Property map containing color values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 * \param  colors			array of dimension N*3
 * \param  number_of_colors	number N of colors in the array
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
                        double *colors = make_LUT().get(),
                        int number_of_colors = 255)
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
                              colors,
                              number_of_colors);
  }
}


/**
 * \brief  Fill the color map for the halfedges of a mesh using a numerical
 * property map
 *
 * \param  g                 mesh to color
 * \param  prop_map			Property map containing numerical values
 * \param  color_pmap		Property map containing color values
 * \param  minMetric         minimum value of the property map
 * \param  maxMetric         maximum value of the property map
 * \param  colors			array of dimension N*3
 * \param  number_of_colors	number N of colors in the array
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
                         double *colors = make_LUT().get(),
                         int number_of_colors = 255)
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
                              colors,
                              number_of_colors);
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
