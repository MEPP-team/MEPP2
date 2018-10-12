Linux / Mac OS X (under Travis CI) : [![Build Status](https://travis-ci.org/MEPP-team/MEPP2.svg?branch=master)](https://travis-ci.org/MEPP-team/MEPP2)

Windows (under AppVeyor) : [![Build status](https://ci.appveyor.com/api/projects/status/wrlfantide2fwcpj/branch/master?svg=true)](https://ci.appveyor.com/project/MEPPteam/mepp2)

Documentation : [![img](https://img.shields.io/badge/Documentation-nightly-brightgreen.svg)](http://liris.cnrs.fr/mepp/doc/nightly/)

# Introduction to MEPP2
The MEPP2 software is a [generic programming](https://en.wikipedia.org/wiki/Generic_programming) oriented [SDK](https://en.wikipedia.org/wiki/Software_development_kit) focusing on 3D mesh processing.

MEPP2 [SDK](https://en.wikipedia.org/wiki/Software_development_kit) offers:
 * **Software components**
   * [Data structures](https://en.wikipedia.org/wiki/Data_structure) implementations and mainly non manifold [Adjacency and Incidence Framework (AIF)](http://Fldc.usb.ve/~vtheok/papers/tesis/A_data_structure_for_eficient_and_fast_management_of_multire.pdf)
   * Generic **mesh processing filters** working on top of third party implementations like
     * [OpenMesh](http://www.openmesh.org/),
     * [CGAL's Surface Mesh](http://doc.cgal.org/latest/Surface_mesh/index.html) or
     * [CGAL's Polyhedral Surface](http://doc.cgal.org/latest/Polyhedron/index.html).
   * Visualization components like an [OpenSceneGraph](http://www.openscenegraph.org/) based ([QT widgets](http://doc.qt.io/qt-5/qtwidgets-index.html)) component allowing for:
     * visualizing 3D mesh with textures (physically-based rendering)
     * building derived GUI applications,
     * production of article oriented illustrations
 * **Tools**
   * [CLI oriented](https://en.wikipedia.org/wiki/Command-line_interface):
     * examples and demos of filters
     * scripting tools allowing for quick macro-filters authoring as well as filter execution and [batch processing](https://en.wikipedia.org/wiki/Batch_processing)
   * [GUI oriented](https://en.wikipedia.org/wiki/Graphical_user_interface) applications
     * A mesh "debugger"

MEPP2 mainly **targets** 
 * **researchers**, or Phd students, looking for tools **allowing them to develop portable algorithms working on surfacic mesh representations**. MEPP2 allows such researchers to:
   * develop mesh implementation independent (on top of various mesh concrete representations), portable (across platforms), robust and [re-usable](https://en.wikipedia.org/wiki/Code_reuse) mesh treatments,
   * quickly compare algorithms and their implementations through [benchmarking](https://en.wikipedia.org/wiki/Benchmarking),
   * establish mesh statistics (or quality indicators),
   * realize mesh treatments,
 * **application developers** wishing to embed mesh treatments within specific industry tool chain, platform specific interfaces, or fully featured applications dedicated to mesh non specialist end users.

# Documentation
http://liris.cnrs.fr/mepp/doc/nightly/
