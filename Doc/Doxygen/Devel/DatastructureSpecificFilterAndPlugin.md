Datastructure specific filter and plugin    {#DatastructureSpecificFilterAndPlugin}
========================================

One may sometime need to write a filter that uses some unique feature of a specific datastructure. In this case the filter is datastructure specific. It can not be written in a generic way. This is easily feasible in MEPP2, and there is no runtime overhead compared to writing an independant datastructure-native code. 

One must take care to not overuse this facility. The priority remains to write generic code. 

## Datastructure specific filter

To write a datastructure specific filter, simply use specific function/method calls in addition to the generic expressions allowed for generic filters.

Example:

```
// loop over vertices
auto iterator_pair = vertices(g);
vertex_iterator vi = iterator_pair.first;
vertex_iterator vi_end = iterator_pair.second;
for(; vi != vi_end; ++vi)  // generic code until here
{
  vi->some_specific_vertex_method(); // specific expression
  ...
}
...
g.some_specific_mesh_method(); // specific expression
some_specific_mesh_function(g); // specific expression
...
```

## Calling a datastructure specific filter from a command line application

The application must create the datastructure object, then call the filter. There is no restriction upon the datastructure definition, provided it belongs to one of the supported library (CGAL, OpenMesh, AIF, PCL). In particular, there is no limitation on the kernel that can be used with the datastructure.

Example:
```
CGAL::Surface_mesh<ExoticKernel::Point_3> mesh;
specific_filter(mesh);
```

## Datastructure specific plugin

A plugin must be written to be able to run the filter inside MEPP2 GUI. One must take care to the following points.

### Use MEPP2 predefined type to avoid a kernel mismatch issue

The datastructure kernel can no more be chosen by the user but is enforced by the GUI code, because the GUI is reponsible for creating the datastructure object.

As a consequence, to avoid compilation errors due to kernel issues, the datastructure used in the plugin and the filter must be one of those pre-defined for the GUI:
- FEVV::MeshPolyhedron (CGAL::Polyhedron_3)
- FEVV::MeshSurface (CGAL::Surface_mesh)
- FEVV::MeshLCC (CGAL::Linear_cell_complex_for_combinatorial_map)
- FEVV::MeshOpenMesh (OpenMesh::PolyMesh_ArrayKernelT)
- FEVV::MeshAIF (MEPP2 AIF)
- FEVV::CGALPointSet (CGAL::Point_set_3)
- FEVV::PCLPointCloud (pcl::PointCloud)

### Disable incompatible datastructures sections in the ..._plugin.h file

As the filter contains datastructure specific functions/methods calls, the plugin will fail to compile with other datastructures. In order to avoid compilation errors, the sections of the '`..._plugin.h`' file (see @ref HowToMakeAQtVisualizationPluginOutOfAFilterPage) calling the '`applyHG()`' function with incompatible datastructures must be commented out.

Example:
```
// the section below is disabled in the ..._plugin.h file
// because the filter is not compatible with OpenMesh
#if 0
#ifdef FEVV_USE_OPENMESH
  void apply(BaseAdapterVisu *_adapter,
             MeshOpenMesh *_mesh,
             FEVV::PMapsContainer *pmaps_bag) override
  {
    applyHG< MeshOpenMesh >(_adapter, _mesh, pmaps_bag);
  }
#endif
#endif
```