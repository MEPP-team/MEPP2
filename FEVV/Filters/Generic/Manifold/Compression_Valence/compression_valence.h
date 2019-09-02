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

#include "FEVV/Filters/Generic/Manifold/Compression_Valence/Compression_Valence_Component.h"

namespace FEVV {
namespace Filters {


/**
 * \brief  Compress the input mesh using the Compression Valence algorithm.
 *
 *         Ref: Ho Lee, Guillaume Lavoué and Florent Dupont, Rate-distortion
 *         optimization for progressive compression of 3D mesh with color
 *         attributes, The Visual Computer, vol. 28, No. 2, pp. 137-153, 2012.
 *         https://perso.liris.cnrs.fr/guillaume.lavoue/travaux/revue/VC2011.pdf
 *
 * \param  g                 input mesh
 * \param  pm                point map
 * \param  v_cm              vertex color map if any, nullptr allowed
 * \param  input_filename    name of the input mesh file (only used to
 *                           known the input file size and compute
 *                           the compression rate) ; empty allowed
 * \param  output_filename   name of the compressed mesh file
 * \param  with_compression  compression flag ; if set to false, only
 *                           mesh simplification is done
 * \param  with_adaptative_quantization
 *                           adaptative quantization flag ; if set to
 *                           false, no adaptative quantization is done
 * \param  max_vertices      maximum number of vertices of the final mesh
 * \param  quantiz_bits      number of quantization bits
 * \param  gt                the geometry traits to use
 *
 * \return  string containing some stats about compression (final
 *          vertices number, number of layers, calculation time...)
 *
 * \sa      the simplified variant that use the default geometry traits
 *          of the mesh.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
std::string
compression_valence(HalfedgeGraph &g,
                    PointMap *pm,
                    VertexColorMap *v_cm, /* nullptr allowed */
                    const std::string &input_filename,
                    const std::string &output_filename,
                    bool with_compression,
                    bool with_adaptative_quantization,
                    int max_vertices,
                    int quantiz_bits,
                    const GeometryTraits &gt)
{
  Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >
      compr_valence_component;

  std::string result = compr_valence_component.Main_Function(
      g,    /* mesh */
      pm,   /* point map */
      v_cm, /* vertex color map */
      input_filename,
      output_filename,              /*File_Name*/
      quantiz_bits,                 /*Qbit*/
      max_vertices,                 /*Number_vertices*/
      false,                        /*Is_normal_flipping_selected*/
      false,                        /*Is_use_metric_selected*/
      0,                            /*Metric_threshold*/
      false,                        /*Is_use_forget_metric_selected*/
      0,                            /*Forget_metric_value*/
      with_compression,             /*Is_compression_selected*/
      with_adaptative_quantization, /*Is_adaptive_quantization_selected*/
      true                          /*Is_bijection_selected*/
  );

  return result;
}


/**
 * \brief  Compress the input mesh using the Compression Valence algorithm.
 *
 *         Ref: Ho Lee, Guillaume Lavoué and Florent Dupont, Rate-distortion
 *         optimization for progressive compression of 3D mesh with color
 *         attributes, The Visual Computer, vol. 28, No. 2, pp. 137-153, 2012.
 *         https://perso.liris.cnrs.fr/guillaume.lavoue/travaux/revue/VC2011.pdf
 *
 *         Use the default geometry traits of the mesh.
 *
 * \param  g                 input mesh
 * \param  pm                point map
 * \param  v_cm              vertex color map if any, nullptr allowed
 * \param  input_filename    name of the input mesh file (only used to
 *                           known the input file size and compute
 *                           the compression rate) ; empty allowed
 * \param  output_filename   name of the compressed mesh file
 * \param  with_compression  compression flag ; if set to false, only
 *                           mesh simplification is done
 * \param  with_adaptative_quantization
 *                           adaptative quantization flag ; if set to
 *                           false, no adaptative quantization is done
 * \param  max_vertices      maximum number of vertices of the final mesh
 * \param  quantiz_bits      number of quantization bits
 *
 * \return  string containing some stats about compression (final
 *          vertices number, number of layers, calculation time...)
 *
 * \sa      the variant that use the geometry traits provided by the user.
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
std::string
compression_valence(HalfedgeGraph &g,
                    PointMap *pm,
                    VertexColorMap *v_cm, /* nullptr allowed */
                    const std::string &input_filename,
                    const std::string &output_filename,
                    bool with_compression = false,
                    bool with_adaptative_quantization = false,
                    int max_vertices = 100,
                    int quantiz_bits = 10)
{
  GeometryTraits gt(g);
  return compression_valence< HalfedgeGraph,
                              PointMap,
                              VertexColorMap,
                              GeometryTraits >(g,
                                               pm,
                                               v_cm,
                                               input_filename,
                                               output_filename,
                                               with_compression,
                                               with_adaptative_quantization,
                                               max_vertices,
                                               quantiz_bits,
                                               gt);
}


} // namespace Filters
} // namespace FEVV
