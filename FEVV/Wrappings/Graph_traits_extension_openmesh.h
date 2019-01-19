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

/*
 * Specialization of graph traits extension
 * for Openmesh
 */

#include "Graph_traits_extension.h"
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

namespace FEVV {


/*
 * ELO-note:
 *   no known function of OpenMesh giving the right
 *   information ; use generic version from
 *   Graph_traits_extension.h
 */


} // namespace FEVV


// ELO note: Clear mesh
// see
//   https://www.openmesh.org/media/Documentations/OpenMesh-5.2-Documentation/a00198.html#details
//
//   void OpenMesh::Concepts::KernelT< FinalMeshItems >::clear()
//   Delete all items, i.e.
//   clear all item containers. The properties will also be removed from the
//   mesh

