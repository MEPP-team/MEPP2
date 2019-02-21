#pragma once

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "Visualization/Plugins/PluginInterface.h"

#include <QStringList>
#include "Dialogs/helloworld_dialog.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

//#include "FEVV/Filters/Generic/scaling.hpp" // A) include the header of the filter corresponding to your operation
#include "Examples/Generic/Helloworld/helloworld_filter.hpp"

#include "FEVV/Wrappings/properties.h"

#ifdef FEVV_USE_CGAL
#include "FEVV/Wrappings/properties_polyhedron_3.h"
#include "FEVV/Wrappings/properties_surface_mesh.h"
#include "FEVV/Wrappings/properties_linear_cell_complex.h"
#endif // FEVV_USE_CGAL
#ifdef FEVV_USE_OPENMESH
#include "FEVV/Wrappings/properties_openmesh.h"
#endif // FEVV_USE_OPENMESH
#ifdef FEVV_USE_AIF
#include "FEVV/Wrappings/properties_aif.h"
#endif // FEVV_USE_AIF
#endif

namespace FEVV {

class HelloworldPlugin : public QObject,
                         public Generic_PluginInterface,
                         public BasePlugin
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "HelloworldPlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  HelloworldPlugin() = default;
  ~HelloworldPlugin() = default;

public:
  void init() override { init(1.0, 1.0, 1.0); }

  void init(double _x, double _y, double _z)
  {
    value_x = _x;
    value_y = _y;
    value_z = _z;
  }

  void reset() override
  {
    init();

    emit resetSignal();
  }

  void addParameters(BaseWindow *_window) override
  {
    window = _window;
    if(window == nullptr || !window->isInit())
    {
      std::cerr << "BaseWindow is null or not initialized." << std::endl;
      return;
    }
  }

  template< typename HalfedgeGraph >
  void process(HalfedgeGraph *_mesh, FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "Asking to Helloworld mesh ! " << std::endl;

    // retrieve or create vertex-color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
		                            HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;
    if(has_map(*pmaps_bag, FEVV::vertex_color))
    {
      std::cout << "use existing vertex-color map" << std::endl;
      v_cm = get_property_map(FEVV::vertex_color, *_mesh, *pmaps_bag);
    }
    else
    {
      std::cout << "create vertex-color map" << std::endl;
      v_cm = make_property_map(FEVV::vertex_color, *_mesh);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_color, *_mesh, *pmaps_bag, v_cm);
    }

    // retrieve or create vertex-normal property map
    using VertexNormalMap =
        typename FEVV::PMap_traits< FEVV::vertex_normal_t,
		                            HalfedgeGraph >::pmap_type;
    VertexNormalMap v_nm;
    if(has_map(*pmaps_bag, FEVV::vertex_normal))
    {
      std::cout << "use existing vertex-normal map" << std::endl;
      v_nm = get_property_map(FEVV::vertex_normal, *_mesh, *pmaps_bag);
    }
    else
    {
      std::cout << "create vertex-normal map" << std::endl;
      v_nm = make_property_map(FEVV::vertex_normal, *_mesh);
      // store property map in property maps bag
      put_property_map(FEVV::vertex_normal, *_mesh, *pmaps_bag, v_nm);
    }

    // retrieve point property map (aka geometry)
    auto pm = get(boost::vertex_point, *_mesh);

    // apply filter
    helloworld_filter(*_mesh, pm, v_cm, v_nm);

    std::cout << "Helloworld mesh of " << value_x << ";" << value_y << ";"
              << value_z << "." << std::endl;
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    std::cout << "here applyHG ! " << std::endl;

    process(_mesh, pmaps_bag);

    SimpleViewer< HalfedgeGraph > *viewer =
        dynamic_cast< SimpleViewer< HalfedgeGraph > * >(_adapter->getViewer());
    if(viewer)
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, false);

    // comment next line to keep parameters values between calls
    reset();

    viewer->frame();
  }

#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_CGAL
  void apply(BaseAdapterVisu *_adapter,
             MeshLCC *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshLCC >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshSurface *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshSurface >(_adapter, _mesh, pmaps_bag);
  }

  void apply(BaseAdapterVisu *_adapter,
             MeshPolyhedron *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshPolyhedron >(_adapter, _mesh, pmaps_bag);
  }
#endif

#ifdef FEVV_USE_AIF
  void apply(BaseAdapterVisu *_adapter,
             MeshAIF *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshAIF >(_adapter, _mesh, pmaps_bag);
  }
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "HelloworldPlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    // setup and display filter parameters dialog window
    HelloworldDialog dialog;
    dialog.setParameters(value_x, value_y, value_z);

    // get filter parameters from dialog window
    if(dialog.exec() == QDialog::Accepted)
    {
      dialog.getParameters(value_x, value_y, value_z);

      SimpleWindow *sw = static_cast< SimpleWindow * >(
          window); // dynamic_cast fails under OS X

      sw->onModificationParam("helloworld_qt_p", this);
      sw->onApplyButton();

      return true;
    }

    return false;
  }

signals:
  void resetSignal();

protected:
  // filter parameters
  double value_x = 0;
  double value_y = 0;
  double value_z = 0;
};

} // namespace FEVV

