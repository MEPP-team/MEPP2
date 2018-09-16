#ifndef DecompressionValencePlugin_H
#define DecompressionValencePlugin_H

#if(_MSC_VER >= 1400)
#ifndef _SCL_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS
#endif
#endif

#include "Visualization/Plugins/PluginInterface.h"

#include <QStringList>
#include "Dialogs/DialogDecompressionValence1.h"

#ifndef Q_MOC_RUN // MT : very important to avoid the error : ' Parse error at
                  // "BOOST_JOIN" ' -> (qt4 pb with boost)
#include "Visualization/PluginFilters/BasePlugin.h"
#include "Visualization/SimpleViewer.h"

#include "Visualization/SimpleWindow.h"

#include "FEVV/Filters/Generic/Manifold/Compression_Valence/decompression_valence.h" // A) include the header of the filter corresponding to your operation

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
#endif // Q_MOC_RUN


namespace FEVV {

class DecompressionValencePlugin : public QObject,
                                   public Generic_PluginInterface,
                                   public BasePlugin
{
  Q_OBJECT
  Q_INTERFACES(FEVV::Generic_PluginInterface)
#if(FEVV_USE_QT5) // see at the end of .cpp for QT4
  Q_PLUGIN_METADATA(IID "DecompressionValencePlugin")
#endif

  /*public:
      using BasePlugin::apply;*/
public:
  DecompressionValencePlugin() = default;
  ~DecompressionValencePlugin() = default;

public:
  void init() override { init(true, "example.p3d", true, false); }

  void init(bool _forceCompute,
            const std::string &_p3dFilePath,
            bool _write_info,
            bool _write_intermediate_meshes)
  {
    *value_forceCompute = _forceCompute;

    p3dFilePath = _p3dFilePath;
    write_info = _write_info;
    write_intermediate_meshes = _write_intermediate_meshes;
    keep_intermediate_meshes = true;
    intermediate_meshes_void = nullptr;
    intermediate_vertexColorMaps_void = nullptr;
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
    std::cout << "Asking to apply DecompressionValence filter on file '"
              << p3dFilePath << "'!" << std::endl;

#if 0
        //TODO-elo DBG
        std::cout << "Decompression Valence Plugin DBG infos" << std::endl;
        std::cout << " * p3dFilePath=" << p3dFilePath << std::endl;
        std::cout << " * write_info=" << write_info << std::endl;
        std::cout << " * write_intermediate_meshes=" << write_intermediate_meshes << std::endl;
        std::cout << " * keep_intermediate_meshes=" << keep_intermediate_meshes << std::endl;
#endif

    // retrieve geometry property map
    auto pm = get(boost::vertex_point, *_mesh);

    // create vertex color property map
    using VertexColorMap =
        typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                    HalfedgeGraph >::pmap_type;
    VertexColorMap v_cm;

    // create intermediate meshes storage
    std::vector< HalfedgeGraph * > *intermediate_meshes = nullptr;
    std::vector< VertexColorMap * > *intermediate_vertexColorMaps = nullptr;
    if(keep_intermediate_meshes)
    {
      intermediate_meshes = new std::vector< HalfedgeGraph * >;
      intermediate_vertexColorMaps = new std::vector< VertexColorMap * >;
    }

    // apply Decompression Valence filter
    std::string result;
    try
    {
      result =
          FEVV::Filters::decompression_valence(*_mesh,
                                               &pm,
                                               &v_cm,
                                               p3dFilePath,
                                               write_info,
                                               intermediate_meshes,
                                               intermediate_vertexColorMaps,
                                               write_intermediate_meshes,
                                               keep_intermediate_meshes);

      // existing property maps are no more valid due to topological changes
      // purge the property maps bag except the vertex color map
      FEVV::PMapsContainer new_pmaps_bag;
      FEVV::put_property_map(FEVV::vertex_color, *_mesh, new_pmaps_bag, v_cm);
      *pmaps_bag = new_pmaps_bag;
    }
    catch(std::runtime_error &e)
    {
      std::cout << e.what() << std::endl;
    }

    intermediate_meshes_void = static_cast< void * >(intermediate_meshes);
    intermediate_vertexColorMaps_void =
        static_cast< void * >(intermediate_vertexColorMaps);

    std::string message("Decompression completed.\n\n");
    message.append(result);
    QMessageBox::information(0, "", QString::fromStdString(message));
  }

  template< typename HalfedgeGraph >
  void applyHG(BaseAdapterVisu *_adapter,
               HalfedgeGraph *_mesh,
               FEVV::PMapsContainer *pmaps_bag)
  {
    if(*value_forceCompute)
      process(_mesh, pmaps_bag);

    SimpleViewer< HalfedgeGraph > *viewer =
        dynamic_cast< SimpleViewer< HalfedgeGraph > * >(_adapter->getViewer());

    if(viewer)
    {
      using VertexColorMap =
          typename FEVV::PMap_traits< FEVV::vertex_color_t,
                                      HalfedgeGraph >::pmap_type;

      // redraw main mesh -> which is NOW the uncompressed one
      viewer->draw_or_redraw_mesh(_mesh, pmaps_bag, true, true, "uncompressed");

      // space_time mode ON
      viewer->m_space_time = true;

      // draw all intermediate meshes
      auto intermediate_meshes =
          static_cast< std::vector< HalfedgeGraph * > * >(
              intermediate_meshes_void);
      auto intermediate_vertexColorMaps =
          static_cast< std::vector< VertexColorMap * > * >(
              intermediate_vertexColorMaps_void);
      if(intermediate_meshes_void && intermediate_vertexColorMaps_void)
      {
        // loop over intermediate meshes
        for(int i = 0; i < intermediate_meshes->size();
            i++) // original ELO loop
                 // if (intermediate_meshes->size()) for (int i =
                 // intermediate_meshes->size()-1; i >= 0; i--) // MTO : just a
                 // test, try reverse loop for time mode
        {
          HalfedgeGraph *mesh_i = (*intermediate_meshes)[i];
          VertexColorMap *v_cm_i = (*intermediate_vertexColorMaps)[i];

          FEVV::PMapsContainer *pmaps_bag_i = new FEVV::PMapsContainer;
          FEVV::put_property_map(
              FEVV::vertex_color, *mesh_i, *pmaps_bag_i, *v_cm_i);

          // draw intermediate mesh
          viewer->draw_or_redraw_mesh(mesh_i,
                                      pmaps_bag_i,
                                      false,
                                      false,
                                      std::string("level") + std::to_string(i));
        }
      }
      delete(intermediate_meshes);
      intermediate_meshes_void = nullptr;
      delete(intermediate_vertexColorMaps);
      intermediate_vertexColorMaps_void = nullptr;
    }

    // TODO-elo-to-keep-parameters-between-calls   reset();

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
#if 0
		//TODO-elo  restore when Compression Valence compiles with AIF
		applyHG<MeshAIF>(_adapter, _mesh, pmaps_bag);
#else
    QMessageBox::information(
        0,
        "",
        QObject::tr(
            "Decompression Valence filter is not yet compatible with AIF!"));
#endif
  }
#endif

  QStringList Generic_plugins() const override
  {
    return QStringList() << "DecompressionValencePlugin";
  }

  bool Generic_plugin(const QString &plugin) override
  {
    DialogDecompressionValence1 dial1;
    dial1.setDecompressionValenceParams(p3dFilePath,
                                        write_info,
                                        write_intermediate_meshes,
                                        keep_intermediate_meshes);
    if(dial1.exec() == QDialog::Accepted)
    {
      dial1.getDecompressionValenceParams(p3dFilePath,
                                          write_info,
                                          write_intermediate_meshes,
                                          keep_intermediate_meshes);

      SimpleWindow *sw = static_cast< SimpleWindow * >(
          window); // dynamic_cast fails under OS X

      sw->onModificationParam("decompressionvalence_qt_p", this);
      sw->onApplyButton();

      return true;
    }

    return false;
  }

signals:
  void resetSignal();

protected:
  bool *value_forceCompute = new bool(false);

  // filter parameters
  std::string p3dFilePath;
  bool write_info;
  bool keep_intermediate_meshes;  // keep into RAM
  bool write_intermediate_meshes; // write to files
  void *intermediate_meshes_void;
  void *intermediate_vertexColorMaps_void;
};

} // namespace FEVV

#endif // DecompressionValencePlugin_H
