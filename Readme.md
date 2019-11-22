![MEPP](https://perso.liris.cnrs.fr/guillaume.lavoue/teaser.jpg)

MEPP is a C++, cross-platform, software development
kit (SDK) for **processing and visualizing 3D surface
meshes**. It provides both an application programming
interface (API) for creating new processing filters and a
graphical user interface (GUI) that can be used directly
and integrate new filters as plugins.

**It allows to process and visualize both static and
dynamic 3D surface meshes, with color attributes and
texture information.**
It is intended for engineers, researchers, but also for
students thanks to a simple use, facilitated by the
proposed architecture and extensive documentation.
This platform is based on Qt, OpenSceneGraph, Boost
and Eigen; optional dependencies include FBX, Draco,
CImg.

**Characteristics**:
C++, Open-source, Multi-plateform (Windows, Linux,
Mac OS X), Compilation with CMake, Quick and easy installation thanks to
extensive documentation.

**Rendering**: the GUI integrates several shaders such as Blinn-Phong and Cook-Torrance (for
physical based rendering). It can render both static and dynamic 3D
surface meshes, with color attributes and
texture information.

**Data structures**: MEPP is generic programming oriented. It offers
an abstraction layer that provides interoperability over several third party data structures: OpenMesh, CGAL Surface Mesh, CGAL Polyhedral Surface, CGAL Linear Cell Complex, AIF (Adjacency and Incidence Framework).


**Installation and developer documentation is [available on the Wiki](https://github.com/MEPP-team/MEPP2/wiki)**


![MEPP-GUI](https://projet.liris.cnrs.fr/mepp/images/mepp2/MEPP2-GUI.PNG)


Linux / Mac OS X (under Travis CI) : [![Build Status](https://travis-ci.org/MEPP-team/MEPP2.svg?branch=master)](https://travis-ci.org/MEPP-team/MEPP2)

Windows (under AppVeyor) : [![Build status](https://ci.appveyor.com/api/projects/status/wrlfantide2fwcpj/branch/master?svg=true)](https://ci.appveyor.com/project/MEPPteam/mepp2)

Documentation : [![img](https://img.shields.io/badge/Documentation-nightly-brightgreen.svg)](http://liris.cnrs.fr/mepp/doc/nightly/)
