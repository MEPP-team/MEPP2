The FEVV software is Copyright of University of Lyon, 2012 - 2017.

FEVV is distributed under the GNU Lesser General Public License Version 3.0 (refer to the accompanying file LICENSE-lgpl-3.0.txt or a copy at [https://www.gnu.org/licenses/lgpl-3.0.html](https://www.gnu.org/licenses/lgpl-3.0.html)


## FIXME To document cleanly
 * Refer bellow for the many embedded FIXMEs
 * The **GUI IS GPL** ! Deal with it properly in the header files...
 * `CGAL/Cartesian.h`, that in turn has many sub-includes (hence reading the license header of that file does not suffice), belongs to the [CGAL's Cartesian_kernel](https://github.com/CGAL/cgal/blob/master/Cartesian_kernel/include/CGAL/Cartesian.h). Although this Cartesian_kernel package [is not documented in the package review](https://doc.cgal.org/latest/Manual/packages.html) we can hope the [legal header of Cartesian.h](https://github.com/CGAL/cgal/blob/master/Cartesian_kernel/include/CGAL/Cartesian.h) suffices to make it LGPL.
 * Including `CGAL/IO/Polyhedron_iostream.h` produces the inclusion of [`cgal/Polyhedron_IO/include/CGAL/IO/Polyhedron_iostream.h`](https://github.com/CGAL/cgal/blob/master/Polyhedron_IO/include/CGAL/IO/Polyhedron_iostream.h). Alas the `Polyhedron_IO` package [is not documented in the package review](https://doc.cgal.org/latest/Manual/packages.html). Nevertheless [Polyhedron_iostream.h legal header](https://github.com/CGAL/cgal/blob/master/Polyhedron_IO/include/CGAL/IO/Polyhedron_iostream.h) makes it **GPL**.
 * Note: the cmake configuration flag for CGAL surfaces (`BUILD_USE_CGAL`) doesn't enable to distinguish GPL (Polyhedron, SurfaceMesh) form LGPL (LCC)
 * MSDM2 is GPL since it uses AABB:

   | Package         |    License    | Included headers / Notes |
   | --------------- | ------------- | ------------------------ |
   |[3D Fast Intersection and Distance Computation](https://doc.cgal.org/latest/AABB_tree/index.html)| [GPL](https://doc.cgal.org/latest/Manual/packages.html#PkgAABB_treeSummary)| For `CGAL/AABB_tree.h`, `CGAL/AABB_traits.h` and `CGAL/AABB_face_graph_triangle_primitive.h` that are only used by FEVV/Filters/CGAL/Surface_mesh/msdm2.h |
  * What is the purpose of BUILD_USE_FBX


## FEVV library dependencies

### Mandatory library dependencies
| Package         |    License    | Included headers / Notes |
| --------------- | ------------- | ------------------------ |
|[Boost libraries](http://www.boost.org/) | [Boost license 1.0](http://www.boost.org/users/license.html) | License seems to fall within the permissive [MIT license grade/category](http://law.stackexchange.com/questions/91/is-there-any-difference-in-meaning-between-the-boost-and-mit-software-licenses).<br> Used packages: [algorithm](https://www.boost.org/doc/libs/1_67_0/libs/algorithm/doc/html/index.html),[any](https://www.boost.org/doc/libs/1_67_0/doc/html/any.html),  [assert](https://www.boost.org/doc/libs/1_67_0/libs/assert/doc/html/assert.html), [concept](https://www.boost.org/doc/libs/1_67_0/libs/concept_check/concept_check.htm),  [filesystem](https://www.boost.org/doc/libs/1_67_0/libs/filesystem/doc/index.htm), [foreach](https://www.boost.org/doc/libs/1_67_0/doc/html/foreach.html), [fusion](https://www.boost.org/doc/libs/1_67_0/libs/fusion/doc/html/), [graph-library](https://www.boost.org/doc/libs/1_67_0/libs/graph/doc/), [numeric_conversion](https://www.boost.org/doc/libs/1_67_0/libs/numeric/conversion/doc/html/index.html), [property_map](https://www.boost.org/doc/libs/1_67_0/libs/property_map/doc/property_map.html) [property_tree](https://www.boost.org/doc/libs/1_67_0/doc/html/property_tree.html), [random](https://www.boost.org/doc/libs/1_67_0/doc/html/boost_random.html), [range](https://www.boost.org/doc/libs/1_67_0/libs/range/doc/html/index.html), [smart_ptr](https://www.boost.org/doc/libs/1_67_0/libs/smart_ptr/doc/html/smart_ptr.html)...|
| CGAL-[BOOST graph Library](https://doc.cgal.org/latest/BGL/index.html#Chapter_CGAL_and_the_Boost_Graph_Library)| [LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgBGLSummary) | CGAL provided adapters for using [BGL (Boost Graph Library)](https://www.boost.org/doc/libs/1_67_0/libs/graph/doc/) algorithms on surfaces|
| CGAL-[Handles and Circulators](https://doc.cgal.org/latest/Circulator/index.html#Chapter_Handles_Ranges_and_Circulators)|[LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgHandlesAndCirculatorsSummary)| Used by compression-valence (`CGAL/circulator.h`)|
|[Eigen3](https://eigen.tuxfamily.org/)|[MPL2](http://eigen.tuxfamily.org/index.php?title=Main_Page#License)||

### Third party source code redistribution
The following third-party libraries (not part of FEVV and available under their own licenses), are  redistributed along with FEVV sources (for the users' convenience):
 * the `External/` sub-directory holds CGAL content FIXME FIXME FIXME: what package and what license ?
 * the `Visualization/osgQt/` sub-directory holds the versions 3.2 and 3.4 of original [OSGQT](https://github.com/openscenegraph/osgQt) source code (this source level redistribution choice was motivated by packaging reasons: their packaging depending too much on the platform/distributions respective choices)
 FIXME FIXME FIXME document [orignal license](https://github.com/openscenegraph/osgQt/blob/master/LICENSE.txt)

### Optional library dependencies: REALLY ?

| Package     |  License   | Included headers / Notes | Build Flag |
| ----------- | ---------- | ------------------------ | ---------- |
|[CGAL-Combinatorial_map](https://doc.cgal.org/latest/Combinatorial_map/index.html#Chapter_Combinatorial_Maps)|[LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgCombinatorialMapsSummary)| For `CGAL/Combinatorial_map_constructors.h` | `BUILD_USE_CGAL`|
|[CGAL-LinearCellComplex](https://doc.cgal.org/latest/Linear_cell_complex/index.html#Chapter_Linear_Cell_Complex)|[LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgLinearCellComplexSummary)| For `CGAL/Linear_cell_complex_for_combinatorial_map.h`, `CGAL/Linear_cell_complex_incremental_builder.h`|  `BUILD_USE_CGAL`|
| [CGAL-2D and 3D Linear geometry kernel](https://doc.cgal.org/latest/Kernel_23/index.html#Chapter_2D_and_3D_Geometry_Kernel)| [LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgKernel23Summary)| For `CGAL/Cartesian.h` | `BUILD_USE_CGAL`|
|[CGAL/config.h](https://github.com/CGAL/cgal/blob/master/Installation/include/CGAL/config.h)|[LGPL](https://github.com/CGAL/cgal/blob/master/Installation/include/CGAL/config.h)| Used by `test_not_2_manifold_linear_cell_complex.cpp`|  `BUILD_USE_CGAL`|
[CGAL-Kernel_23](https://doc.cgal.org/latest/Kernel_23/index.html#Chapter_2D_and_3D_Geometry_Kernel)| [LGPL](https://doc.cgal.org/latest/Manual/packages.html#PkgKernel23Summary) | For `CGLA/Kernel_traits.h`,`CGAL/basic.h`, `CGAL/Kernel/global_functions.h`, `CGAL/Exact_predicates_inexact_constructions_kernel.h`...|  `BUILD_USE_CGAL`|
| [PCL](https://github.com/PointCloudLibrary/pcl) | [BSD](https://github.com/PointCloudLibrary/pcl/blob/master/LICENSE.txt)||`BUILD_USE_PCL`|
| PCL: [Flann](http://www.cs.ubc.ca/research/flann/)|[BSD](http://www.cs.ubc.ca/research/flann/)|| `BUILD_USE_PCL` |
|[OpenMesh](http://www.openmesh.org/) | [BSD3](http://www.openmesh.org/license/)| | `BUILD_USE_OPENMESH` |
|[VTK](https://www.vtk.org/)|[BSD](https://www.vtk.org/licensing/)||`BUILD_USE_VTK`|
|[Qt (4 or 5)](https://en.wikipedia.org/wiki/Qt_(software))| [LGPLV2.1, LGPLv3](https://www.qt.io/licensing/) | [Core: LGPLV2 or V3](http://doc.qt.io/qt-5/qtcore-index.html#licenses-and-attributions), [Gui: LGPLV2 or V3](http://doc.qt.io/qt-5/qtgui-index.html#licenses-and-attributions) sa [licenses in QT](http://doc.qt.io/archives/qt-5.5/licensing.html#licenses-used-in-qt) | `BUILD_USE_GUI` and `BUILD_USE_QT5 `|
|[OpenSceneGraph (OSG)](http://www.openscenegraph.org/)|[OSGPL (LGPL like)](http://trac.openscenegraph.org/projects/osg/wiki/Legal)|Used sub-libraries: osgDB, osgGA, osgFX, osgQT, osgShadow, osgText,osgUtil osgViewer...| `BUILD_USE_GUI` |

#### Contaminating licenses
| Package         |    License    | Included headers / Notes | Build Flag |
| --------------- | ------------- | ------------------------ | ---------- |
| [CGAL-Polyhedron](https://doc.cgal.org/latest/Polyhedron/index.html)| [GPL](https://doc.cgal.org/latest/Manual/packages.html#PkgPolyhedronSummary) | For `CGAL/Polyhedron_3.h`, `CGAL/Polyhedron_items_with_id_3.h` but also for files outside the package like `CGAL/IO/Polyhedron_iostream.h`|`BUILD_USE_CGAL`|
[CGAL-SurfaceMesh](https://doc.cgal.org/latest/Surface_mesh/index.html#Chapter_3D_Surface_mesh)| [GPL](https://doc.cgal.org/latest/Manual/packages.html#PkgSurfaceMeshSummary) | For `CGAL/Surface_mesh.h`, `CGAL/Surface_mesh_shortest_path.h`|`BUILD_USE_CGAL`|

### Concerning the QT modules:
 * 3D-USE stays away from the [GPL poisened pill of some QT modules](http://doc.qt.io/qt-5/qtmodules.html#gpl-licensed-addons)
 * Reference gateway: [licenses used in QT](http://doc.qt.io/archives/qt-5.5/licensing.html#licenses-used-in-qt) gateway.
 * List of used packages within modules:
    * **Core module**: [QDebug](https://doc.qt.io/archives/qt-5.5/qdebug.html), [QDir](http://doc.qt.io/qt-5/qdir.html), [QDirIterator](http://doc.qt.io/qt-5/qdiriterator.html), [QEvent](http://doc.qt.io/qt-5/qevent.html), [QFile](http://doc.qt.io/qt-5/qfile.html) and [QFileInfo](http://doc.qt.io/qt-5/qfileinfo.html), [QMutex](http://doc.qt.io/qt-5/qmutex.html), [QPluginLoader](http://doc.qt.io/qt-5/qpluginloader.html), [QPointer](http://doc.qt.io/qt-5/qpointer.html), [QQueue](http://doc.qt.io/qt-5/qqueue.html), [QSet](http://doc.qt.io/qt-5/qset.html), [QSettings](http://doc.qt.io/qt-5/qsettings.html), [QString](http://doc.qt.io/qt-5/QString.html), [QStringList](http://doc.qt.io/qt-5/qstringlist.html), [Qtime](http://doc.qt.io/qt-5/qtime.html) and [QTimer](http://doc.qt.io/qt-5/qtimer.html), [QTextStream](http://doc.qt.io/qt-5/qtextstream.html),
    * **Gui module (includes OpenGL)**: [QFont](http://doc.qt.io/qt-5/qfont.html), [QGLWidget](http://doc.qt.io/qt-5/qglwidget.html), [QIcon](http://doc.qt.io/qt-5/QIcon.html), [QImage](http://doc.qt.io/qt-5/qimage.html) and [QImageReader](http://doc.qt.io/qt-5/qimagereader.html), [QInputEvent](http://doc.qt.io/qt-5/qinputevent.html), [QResizeEvent](http://doc.qt.io/archives/qt-5.5/qresizeevent.html),
    * **Widgets module**: [QDialog](http://doc.qt.io/qt-5/QDialog.html), [QFileDialog](http://doc.qt.io/qt-5/qfiledialog.html), [QGridLayout](http://doc.qt.io/qt-5/qgridlayout.html), [QGraphicScene](http://doc.qt.io/qt-5/qgraphicsscene.html) and [QGraphicView](http://doc.qt.io/qt-5/qgraphicsview.html), [QHeaderView](http://doc.qt.io/qt-5/qheaderview.html), [QLabel](http://doc.qt.io/qt-5/qlabel.html), [QLineEdit](http://doc.qt.io/qt-5/qlineedit.html), [QListView](http://doc.qt.io/qt-5/qlistview.html) and [QListWidget](http://doc.qt.io/qt-5/qlistwidget.html), [QMainWindow](http://doc.qt.io/qt-5/qmainwindow.html), [QMenu](http://doc.qt.io/qt-5/qmenu.html), [QMessageBox](http://doc.qt.io/qt-5/qmessagebox.html), [QPainter](http://doc.qt.io/qt-5/qpainter.html), [QProgressDialog](http://doc.qt.io/qt-5/qprogressdialog.html), [QPushButton](http://doc.qt.io/qt-5/QPushButton.html), [QTextBrowser](http://doc.qt.io/qt-5/qtextbrowser.html), [QTreeWidget](http://doc.qt.io/archives/qt-5.5/qtreewidgetitem.html),
    * **unknown module**: [QtGlobal](http://doc.qt.io/qt-5/qtglobal.html), [QtPlugin](http://doc.qt.io/qt-5/qtplugin.html),
    * used by osgQt: [QtWebKit](https://wiki.qt.io/Qt_WebKit), QtWebKitWidgets
 * Obtaining the above list of used QT modules is merely a comment of the `grep -rh "#include" Visualization | grep -v "Visualization" | grep -v "<osg/" |  grep -i Q | grep -v QtEvents | grep -v osgQt | sort | uniq` command (QtEvents being used by `osgQt/QGraphicsViewAdapter` and `osgQt` belonging to [OSG](https://github.com/openscenegraph/osgQt) ).



## Technical Notes
* In order to find all the dependencies (searching within `CMakeLists.txt` is not sufficient):
  ```
  grep -ri find_package .
  ```
* In order to find the boost modules that are used within FEVV
  ```
  find . \( -name "*.h" -o -name "*.cpp" -o -name "*.inl" -o -name "*.hxx" \) | grep -v External | grep -v \.\\Bin | xargs grep -h "#include" | grep -i boost | grep -v CGAL | sort | uniq
  ```
* In order to find the CGAL header files that are used within FEVV:
  ```
  find . \( -name "*.h" -o -name "*.cpp" -o -name "*.inl" -o -name "*.hxx" \) | grep -v External | grep -v \.\/Bin | xargs grep -h "#include" | grep CGAL | sort | uniq
  ```
* In order to assert that `FEVV/External` only holds LGPL code make sure that the following command results only differ by `3`
    - `find External -type f -exec grep -l Lesser {} \;  | wc` (returns 444 for CGAL-Version 4.14)
    - `find External -type f  | wc` (that returns 447 for CGAL-Version 4.14)

  The difference comes from the files `CGAL/auto_link/auto_link.h`, `CGAL/boost/graph/named_function_params.h`, `CGAL/boost/graph/named_params_helper.h` that have license `Boost Software License, Version 1.0`.
