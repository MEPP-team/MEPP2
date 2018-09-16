#ifndef DECOMPRESSION_VALENCE_H
#define DECOMPRESSION_VALENCE_H

#include "FEVV/Filters/Generic/Manifold/Compression_Valence/Compression_Valence_Component.h"

namespace FEVV {
namespace Filters {


/**
 * \brief  Uncompress a mesh compressed with Compression Valence algorithm
 *
 *         Ref: Ho Lee, Guillaume Lavoué and Florent Dupont, Rate-distortion
 *         optimization for progressive compression of 3D mesh with color
 *         attributes, The Visual Computer, vol. 28, No. 2, pp. 137-153, 2012.
 *         https://perso.liris.cnrs.fr/guillaume.lavoue/travaux/revue/VC2011.pdf
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  v_cm              output vertex color map
 * \param  input_filename    name of the compressed mesh file
 * \param  do_write_info     if set to true, write some information
 *                           about decompression progression in a text
 *                           file .infos.txt
 * \param  intermediate_meshes
 *                           if not null, used to store the intermediate
 *                           meshe corresponding to each compression
 *                           layer
 * \param  intermediate_vertexColorMaps
 *                           if not null, used to store the vertex color
 *                           map of each intermediate mesh
 * \param  do_write_intermediate_meshes
 *                           write intermediate meshes to files during
 *                           decompression progression
 * \param  gt                the geometry traits to use
 *
 * \return  string containing some stats about compression (number of
 *          layers, calculation time...)
 *
 * \sa      the simplified variant that use the default geometry traits
 *          of the mesh.
 *
 * \ingroup  GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
std::string
decompression_valence(
    HalfedgeGraph &g,
    PointMap *pm,
    VertexColorMap *v_cm,
    const std::string &input_filename,
    bool do_write_info,
    std::vector< HalfedgeGraph * > *intermediate_meshes, /* nullptr allowed */
    std::vector< VertexColorMap * >
        *intermediate_vertex_color_maps, /* nullptr allowed */
    bool do_write_intermediate_meshes,
    const GeometryTraits &gt)
{
  Compression_Valence_Component< HalfedgeGraph, PointMap, VertexColorMap >
      compr_valence_component;

  compr_valence_component.Decompress_Init(g,    /* mesh */
                                          pm,   /* point map */
                                          v_cm, /* vertex color map */
                                          input_filename);

  std::string result = compr_valence_component.Decompression_All_From_File(
      g,    /* mesh */
      pm,   /* point map */
      v_cm, /* vertex color map */
      do_write_info,
      intermediate_meshes,
      intermediate_vertex_color_maps,
      do_write_intermediate_meshes);

  return result;
}


/**
 * \brief  Uncompress a mesh compressed with Compression Valence algorithm
 *
 *         Ref: Ho Lee, Guillaume Lavoué and Florent Dupont, Rate-distortion
 *         optimization for progressive compression of 3D mesh with color
 *         attributes, The Visual Computer, vol. 28, No. 2, pp. 137-153, 2012.
 *         https://perso.liris.cnrs.fr/guillaume.lavoue/travaux/revue/VC2011.pdf
 *
 * \param  g                 input mesh
 * \param  pm                output point map
 * \param  v_cm              output vertex color map
 * \param  input_filename    name of the compressed mesh file
 * \param  do_write_info     if set to true, write some information
 *                           about decompression progression in a text
 *                           file .infos.txt
 * \param  intermediate_meshes
 *                           if not null, used to store the intermediate
 *                           meshe corresponding to each compression
 *                           layer
 * \param  intermediate_vertexColorMaps
 *                           if not null, used to store the vertex color
 *                           map of each intermediate mesh
 * \param  do_write_intermediate_meshes
 *                           write intermediate meshes to files during
 *                           decompression progression
 *
 * \return  string containing some stats about compression (number of
 *          layers, calculation time...)
 *
 * \sa      the variant that use the geometry traits provided by the user.
 *
 * \ingroup  GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename PointMap,
          typename VertexColorMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
std::string
decompression_valence(
    HalfedgeGraph &g,
    PointMap *pm,
    VertexColorMap *v_cm,
    const std::string &input_filename,
    bool do_write_info = false,
    std::vector< HalfedgeGraph * > *intermediate_meshes =
        nullptr, /* optional */
    std::vector< VertexColorMap * > *intermediate_vertex_color_maps =
        nullptr, /* optional */
    bool do_write_intermediate_meshes = false)
{
  GeometryTraits gt(g);
  return decompression_valence< HalfedgeGraph,
                                PointMap,
                                VertexColorMap,
                                GeometryTraits >(g,
                                                 pm,
                                                 v_cm,
                                                 input_filename,
                                                 do_write_info,
                                                 intermediate_meshes,
                                                 intermediate_vertex_color_maps,
                                                 do_write_intermediate_meshes,
                                                 gt);
}


} // namespace Filters
} // namespace FEVV

#endif // DECOMPRESSION_VALENCE_H
