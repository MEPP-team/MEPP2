#ifndef TranslationPlugin_H
#define TranslationPlugin_H

#include "Visualization/PluginFilters/BasePlugin.h"

#include "FEVV/Filters/Generic/translation.hpp"

#include <functional>

namespace FEVV {

class TranslationPlugin : public BasePlugin
{
  /*public:
      using BasePlugin::apply;*/
public:
  TranslationPlugin() = default;
  ~TranslationPlugin() = default;

public:
  void init() override { init(0.0, 0.0, 0.0); }

  void init(double _x, double _y, double _z)
  {
    *value_x = _x;
    *value_y = _y;
    *value_z = _z;
  }

  void reset() override
  {
    *value_x = 0;
    *value_y = 0;
    *value_z = 0;
  }

  void addParameters(BaseWindow *_window) override
  {
    window = _window;
    if(window == nullptr || !window->isInit())
    {
      std::cerr << "BaseWindow is null or not initialized." << std::endl;
      return;
    }

    // window->setParam( "Translation X", value_x, "translation_p", this );
    // window->setParam( "Translation Y", value_y, "translation_p", this );
    // window->setParam( "Translation Z", value_z, "translation_p", this );
  }

#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    std::cout << "Asking to translate mesh ! " << std::endl;

    auto pm = get(boost::vertex_point, *_mesh);

    Filters::translate(*_mesh, pm, *value_x, *value_y, *value_z);

    std::cout << "Translate mesh of " << *value_x << ";" << *value_y << ";"
              << *value_z << "." << std::endl;

    SimpleViewer< MeshOpenMesh > *viewerOM =
        dynamic_cast< SimpleViewer< MeshOpenMesh > * >(_adapter->getViewer());
    if(viewerOM)
      viewerOM->redrawMesh(false, _mesh, pmaps_bag, &pm);

    reset();

    viewerOM->frame();
  }
#endif

protected:
  double *value_x = new double(0.0);
  double *value_y = new double(0.0);
  double *value_z = new double(0.0);
};

} // namespace FEVV

#endif // TranslationPlugin_H
