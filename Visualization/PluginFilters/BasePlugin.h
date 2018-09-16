#ifndef BasePlugin_H
#define BasePlugin_H

#include "FEVV/Wrappings/Wrappings.h"

#include "Visualization/BaseWindow.h"
#include "Visualization/BaseAdapterVisu.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
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
      std::cerr << "This function has not been overrided ! BasePlugin::apply(
  HalfedgeGraph )" << std::endl;
  }*/

  virtual void apply(BaseAdapterVisu *_adapter,
                     boost::any _mesh_any,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overrided ! BasePlugin::apply( "
                 "boost::any )"
              << std::endl;
  }

#ifdef FEVV_USE_OPENMESH
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshOpenMesh *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overrided ! BasePlugin::apply( "
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
        << "This function has not been overrided ! BasePlugin::apply( MeshLCC )"
        << std::endl;
  }

  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshSurface *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overrided ! BasePlugin::apply( "
                 "MeshSurface )"
              << std::endl;
  }

  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshPolyhedron *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr << "This function has not been overrided ! BasePlugin::apply( "
                 "MeshPolyhedron )"
              << std::endl;
  }
#endif

#ifdef FEVV_USE_AIF
  virtual void apply(BaseAdapterVisu *_adapter,
                     MeshAIF *_mesh,
                     FEVV::PMapsContainer *pmaps_bag)
  {
    std::cerr
        << "This function has not been overrided ! BasePlugin::apply( MeshAIF )"
        << std::endl;
  }
#endif

  virtual void reset() = 0;

protected:
  BaseWindow *window = nullptr;
};

} // namespace FEVV

#endif // BasePlugin_H
