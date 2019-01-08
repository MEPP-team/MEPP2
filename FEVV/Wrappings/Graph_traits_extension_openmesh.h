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

