// Copyright (c) 2012-2022 University of Lyon and CNRS (France).
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

#include "FEVV/Wrappings/Geometry_traits.h"
#include "FEVV/Wrappings/properties.h"
#include "FEVV/Types/Material.h"

/**
 * \ingroup GenericManifoldFilters
 */
template< typename HalfedgeGraph,
          typename MaterialMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
texture_image_demo_filter(const HalfedgeGraph  &/*g*/,
                          MaterialMap          &mtl_pm,
                          const GeometryTraits &/*gt*/)
{
  // loop over materials
  auto it     = mtl_pm.storage_begin();
  auto it_end = mtl_pm.storage_end();
  for(; it != it_end; ++it)
  {
    auto material = *it;
    std::cout << "processing material '" << material.name << "'" << std::endl;

    // loop over texture images
    for(const auto &texture_filename : { material.ambient_texture_filename,
                                         material.diffuse_texture_filename,
                                         material.specular_texture_filename,
                                         material.emissive_texture_filename,
                                         material.transparency_texture_filename,
                                         material.normal_map_filename,
                                         material.metallic_map_filename,
                                         material.roughness_map_filename })
    {
      if(texture_filename.empty())
        continue;

      // retrieve texture image

      auto &cimg = *(material.images.at(texture_filename));

      // texture image must have at least 3 channels

      if(cimg.spectrum() < 3)
      {
        std::cout << "  skipping non RGB texture '" << texture_filename << "'"
                  << std::endl;
        continue;
      }

      // cycle channels (R->G, G->B, B->R)

      std::cout << "  processing texture '" << texture_filename << "'"
                << std::endl;

      // loop over pixels
      // inefficient but instructive implementation
      for(int row = 0; row < cimg.height(); row++)
      {
        for(int col = 0; col < cimg.width(); col++)
        {
          auto R = cimg(col, row, 0, 0); // cimg(x, y, z, channel)
          auto G = cimg(col, row, 0, 1);
          auto B = cimg(col, row, 0, 2);

          cimg(col, row, 0, 0) = B;
          cimg(col, row, 0, 1) = R;
          cimg(col, row, 0, 2) = G;
        }
      }
    }
  }
}


// Helper function to simplify (syntactic sugar) the call to the filter
template< typename HalfedgeGraph,
          typename MaterialMap,
          typename GeometryTraits = FEVV::Geometry_traits< HalfedgeGraph > >
void
texture_image_demo_filter(const HalfedgeGraph &g,
                          MaterialMap         &mtl_pm)
{
  GeometryTraits gt(g);
  texture_image_demo_filter(g, mtl_pm, gt);
}

