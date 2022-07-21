// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

#pragma once

#include <CGAL/version_macros.h>

#include <CGAL/Cartesian.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <boost/property_map/vector_property_map.hpp>

#if defined _MSC_VER
#define WINDLL_EXPORT __declspec(dllexport)
#else
#define WINDLL_EXPORT
#endif

namespace msdm2 {

using CGALKernel = CGAL::Cartesian< double >;
using CGALPoint  = CGALKernel::Point_3;

using MeshT             = CGAL::Surface_mesh< CGALPoint >;
using vertex_descriptor = boost::graph_traits< MeshT >::vertex_descriptor;

using VertexIndexMapT =
    boost::property_map< MeshT, boost::vertex_index_t >::const_type;
using Msdm2MapT = boost::vector_property_map< double, VertexIndexMapT >;


/**
  \brief	Computes the multiscale MSDM2 metric.

  \param	mesh_degraded  the 1st mesh
  \param	mesh original  the 2nd mesh
  \param	nb_levels	     number of scales used
  \param  msdm2_value	   [output] the computed value of the MSDM2 metric
  \param  msdm2_map  	   [output] the MSDM2 vertex map
 */
WINDLL_EXPORT
void msdm2_surfacemesh(const MeshT &mesh_degraded,  /* input  */
                        const MeshT &mesh_original,  /* input  */
                        const int   nb_levels,       /* input  */
                        double      &msdm2_value,    /* output */
                        Msdm2MapT   &msdm2_pmap);    /* output */

} // namespace msdm2
