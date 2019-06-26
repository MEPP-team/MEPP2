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

#include "FEVV/Wrappings/Wrappings.h"

#include "Visualization/BaseWindow.h"
#include "Visualization/BaseAdapterVisu.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#include "FEVV/Wrappings/properties_cgal_point_set.h"
#endif

#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif

#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif

namespace FEVV {
class BaseWindow;
class BaseAdapterVisu;

class BasePlugin
{
public:
  BasePlugin() = default;
  ~BasePlugin() = default;

public:
  virtual void init() = 0;
  virtual void addParameters(BaseWindow *_window) = 0;

  /*template< typename HalfedgeGraph >
  void apply( BaseAdapterVisu* _adapter, HalfedgeGraph* _mesh,
  FEVV::PMapsContainer *pmaps_bag )
  {
      std::cerr << "This function has not been overridden ! BasePlugin::apply(
  HalfedgeGraph )" << std::endl;
  }*/

  virtual void apply(BaseAdapterVisu *_adapter,
                     void *_mesh_void,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "void )"
              << std::endl;
  }

#ifdef FEVV_USE_OPENMESH
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshOpenMesh *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "MeshOpenMesh )"
              << std::endl;
  }
#endif

#ifdef FEVV_USE_CGAL
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshLCC *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr
        << "This function has not been overridden ! BasePlugin::apply( MeshLCC )"
        << std::endl;
  }

  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshSurface *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "MeshSurface )"
              << std::endl;
  }

  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshPolyhedron *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "MeshPolyhedron )"
              << std::endl;
  }

  virtual void apply(BaseAdapterVisu *_adapter,
                     CGALPointSet *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "CGALPointSet )"
              << std::endl;
  }
#endif

#ifdef FEVV_USE_AIF
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshAIF *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr
        << "This function has not been overridden ! BasePlugin::apply( MeshAIF )"
        << std::endl;
  }
#endif

#ifdef FEVV_USE_PCL
  virtual void apply(BaseAdapterVisu *_adapter,
                     PCLPointCloud *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overridden ! BasePlugin::apply( "
                 "PCLPointCloud )"
              << std::endl;
  }
#endif

  virtual void reset() = 0;

protected:
  BaseWindow *window = nullptr;
};

} // namespace FEVV

