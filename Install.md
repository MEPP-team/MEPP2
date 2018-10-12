
## MEPP2 installation guide (Linux, MacOS and Windows)

## For the impatient
CMake compilation flags short list:
 - `CMAKE_BUILD_TYPE   `: build type (Release or Debug)
 - `BUILD_USE_CGAL     `: use CGAL library (ON or OFF)
 - `BUILD_USE_OPENMESH `: use OpenMesh library (ON or OFF)
 - `BUILD_USE_GUI      `: use GUI (ON or OFF)
 - `BUILD_DOCUMENTATION`: build the documentation (ON or OFF)

## Dependencies
Mandatory dependencies:
 - CMake >= 2.8.11
   - Note: [Boost 1.64 requires cmake 3.8 or newer](https://stackoverflow.com/questions/42123509/cmake-finds-boost-but-the-imported-targets-not-available-for-boost-version)
 - Boost >= 1.59
 - Eigen 3

Optional dependencies:
 - CGAL >= 4.11: used for CGAL data structures
 - OpenMesh >= 6.2: used for OpenMesh data structures
 - Qt 4 or 5: mandatory for building the GUI
 - OpenSceneGraph: mandatory for building the GUI
 - Doxygen and Graphviz: used to generate the documentation

Other optional dependencies are:
 - VTK
 - PCL
 - FBX

## Linux installation

### Install dependencies

Example on Ubuntu 18.04 LTS Bionic Beaver (amd64), released on April 26, 2018:
````
  # CMake
  $ sudo apt install cmake

  # Boost
  $ sudo apt install libboost-dev libboost-thread-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev

  # Eigen 3
  $ sudo apt install libeigen3-dev

  # CGAL 4.11
  $ sudo apt install libcgal-dev

  # --> ! apply patch for CGAL 4.11.0 ! <--
  # ---------------------------------------
  $ cd /usr/include/CGAL/boost/graph
  $ sudo mv io.h io.h.old
  $ sudo wget --no-check-certificate https://download.gforge.liris.cnrs.fr/meppbin/src/cgal411/io.h
  $ sudo mv cgal_bgl_graph_io.h cgal_bgl_graph_io.h.old
  $ sudo wget --no-check-certificate https://download.gforge.liris.cnrs.fr/meppbin/src/cgal411/cgal_bgl_graph_io.h
  $ cd

  # OpenMesh (installation in user home directory)
  $ cd /tmp
  $ wget https://www.openmesh.org/media/Releases/7.1/OpenMesh-7.1.tar.gz
  $ tar -xzf OpenMesh-7.1.tar.gz
  $ cd OpenMesh-7.1 && mkdir build && cd build
  $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-7.1 ..
  $ make
  $ make install && cd

  # Qt
  $ sudo apt install libqt4-dev (or qtdeclarative5-dev but today OpenSceneGraph is only SingleThread with Qt 5, so prefer Qt 4...)

  # OpenSceneGraph
  $ sudo apt install libopenscenegraph-3.4-dev

  # Doxygen and Graphviz
  $ sudo apt install doxygen graphviz

  # -> Optional dependencies <- :
  # -----------------------------

  # VTK
  $ sudo apt install libvtk7-dev (or libvtk6-dev)

  # PCL (installation in user home directory) - don't use libpcl-dev package !
  $ sudo apt install libflann-dev
  $ cd /tmp
  $ wget https://github.com/PointCloudLibrary/pcl/archive/pcl-1.8.1.tar.gz
  $ tar -xzf pcl-1.8.1.tar.gz
  $ cd pcl-pcl-1.8.1 && mkdir build && cd build
  $ cmake -DCMAKE_BUILD_TYPE=Release -DPCL_ONLY_CORE_POINT_TYPES=ON -DBUILD_global_tests=OFF -DCMAKE_INSTALL_PREFIX=$HOME/pcl-1.8.1 ..
  $ make
  $ make install && cd

  # FBX SDK (a readme is available when uncompressing the archive) - primary support, not finished !
  $ mkdir ~/FBX_SDK && mkdir ~/FBX_SDK/2019.0
  $ cd /tmp
  $ wget http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_linux.tar.gz
  $ tar -xzf fbx20190_fbxsdk_linux.tar.gz
  $ ./fbx20190_fbxsdk_linux ~/FBX_SDK/2019.0
  $ ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/release ~/FBX_SDK/2019.0/lib/release
  $ ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/debug ~/FBX_SDK/2019.0/lib/debug
  $ cd
````

### Build stage

Scripting commands for compiling Mepp2:
````
  # get Mepp2 source code
  $ git clone https://github.com/MEPP-team/MEPP2.git

  # create build directory
  $ cd MEPP2 && mkdir build && cd build

  # compile with CGAL, OpenMesh and GUI
  $ cmake -DOPENMESH_DIR="$HOME/OpenMesh-7.1" -DBUILD_USE_GUI=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # compile without CGAL nor OpenMesh, nor GUI
  $ cmake -DBUILD_USE_CGAL=OFF -DBUILD_USE_OPENMESH=OFF -DBUILD_USE_GUI=OFF -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # compile with CGAL, GUI and FBX
  $ export FBX_DIR="$HOME/FBX_SDK/2019.0"
  $ cmake -DBUILD_USE_GUI=ON -DBUILD_USE_FBX=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # generate the documentation
  $ cmake -DBUILD_DOCUMENTATION=ON ..
  $ make doc
````

### Run tests

````
  $ ctest
````

## MacOS installation

### Install dependencies
Mandatory dependencies:
````
  # install homebrew package manager (see http://brew.sh/), then
  $ brew update

  # CMake
  $ brew install cmake

  # Boost: installed with CGAL (dependency)

  # Eigen 3
  $ brew install eigen
````

Optional dependencies:
````
  # CGAL
  $ brew install cgal
  $ brew info cgal

  # -> ONLY if CGAL version above is 4.11 then you have to do that (NO MORE PROBLEM with version >= 4.11.1) <- :
  # ------------------------------------------------------------------------------------------------------------
  $ cd /usr/local/include/CGAL/boost/graph
  $ mv io.h io.h.old
  $ wget --no-check-certificate https://download.gforge.liris.cnrs.fr/meppbin/src/cgal411/io.h
  $ mv cgal_bgl_graph_io.h cgal_bgl_graph_io.h.old
  $ wget --no-check-certificate https://download.gforge.liris.cnrs.fr/meppbin/src/cgal411/cgal_bgl_graph_io.h
  $ cd

  # OpenMesh
  $ brew install open-mesh

  # XQuartz (X.Org X Window System that runs on OS X) -> needed for OpenSceneGraph
  $ brew install Caskroom/cask/xquartz

  # Qt
  $ brew tap cartr/qt4
  $ brew tap-pin cartr/qt4
  $ brew install qt@4

  # OpenSceneGraph (version 3.4.x !)
  $ brew install open-scene-graph
  $ brew info open-scene-graph
  # -> if the version is not 3.4.x then unlink open-scene-graph (brew unlink open-scene-graph) and
  # compile the good version directly from this source : https://download.gforge.liris.cnrs.fr/meppbin/src/OpenSceneGraph-3.4.1-JPEG-patched.tgz with
  $ cd /tmp; wget https://download.gforge.liris.cnrs.fr/meppbin/src/OpenSceneGraph-3.4.1-JPEG-patched.tgz
  $ tar zxf OpenSceneGraph-3.4.1-JPEG-patched.tgz; cd OpenSceneGraph-OpenSceneGraph-3.4.1; mkdir buildR && cd buildR
  $ cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_RULE_MESSAGES=OFF -DBUILD_OSG_APPLICATIONS=OFF -DOSG_USE_QT=OFF -DCMAKE_INSTALL_PREFIX=$HOME/osg-3.4.1
  $ make && make install && cd

  # Doxygen and Graphviz
  $ brew install doxygen graphviz

  # -> Optional dependencies <- :
  # -----------------------------

  # VTK
  $ brew install vtk

  # PCL
  $ brew install pcl

  # FBX SDK - primary support, not finished !
  # get and install the sdk from here: http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_clang_mac.pkg.tgz
  # double click to decompress the .tgz and then install the .pkg by clicking on it
  $ ln -s "/Applications/Autodesk/FBX SDK/2019.0/lib/clang/release" "/Applications/Autodesk/FBX SDK/2019.0/lib/release"
  $ ln -s "/Applications/Autodesk/FBX SDK/2019.0/lib/clang/debug" "/Applications/Autodesk/FBX SDK/2019.0/lib/debug"
  $ export DYLD_LIBRARY_PATH="/Applications/Autodesk/FBX SDK/2019.0/lib/release"
````

### Build stage
For the majority of the commandes refer to the above Linux section. The only difference might be that positionning environment variables is not necessary for OSX (since all dependencies are installed with the packaging system). For example building Mepp2 with CGAL, OpenMesh and GUI is reduced to:
````
  $ cmake -DBUILD_USE_GUI=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

Another example, building with CGAL, OpenMesh, GUI and FBX can be done like this:
````
  $ export FBX_DIR="/Applications/Autodesk/FBX SDK/2019.0"
  $ cmake -DBUILD_USE_GUI=ON -DBUILD_USE_FBX=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

Another example, building with CGAL, OpenMesh, GUI (with OpenSceneGraph-3.4.1-JPEG-patched above) can be done like this:
````
  $ export DYLD_LIBRARY_PATH=$HOME/osg-3.4.1/lib
  $ cmake -DBUILD_USE_GUI=ON -DOSG_DIR=$HOME/osg-3.4.1 -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

### Run tests

````
  $ ctest
````

## Windows installation

### Prerequisites

 - Update Windows 7/8/8.1/10 (services packs and Windows Update)

 - Download and install Visual Studio Express 2015 from [LIRIS host](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/Visual%20Studio%20Express%202015/Visual%20Studio%20Express%202015%20pour%20Windows%20Desktop.rar)
   (download size 7.4 GB, installation size ~>14 GB)

### Installing dependencies

 1. Download (mandatory) the ['core' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64.7z) that delivers CMake, Doxygen, Graphviz, Boost, CGAL, OpenMesh, Eigen 3 for `VS2015_64` (download size 250 MB, installation size ~3.3 GB)

 2. Optionally download the ['addon_01' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_01.7z) for Qt4 and OpenSceneGraph
   (download size 333 MB, installation size ~2.8 GB)

 3. Optionally download the ['addon_02' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_02.7z) for Qt5
   (download size 485 MB, installation size ~2.8 GB)

 4. Optionally download the ['addon_03' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_03.7z) for PCL
   (download size 32 MB, installation size ~417 MB)

 5. Optionally download the ['addon_04' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_04.7z) for VTK
   (download size 247 MB, installation size ~2.9 GB)

 6. Optionally download the ['addon_05' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_05.7z) for FBX
   (download size 55 MB, installation size 678 MB)

 7. Extract the 'core' ('MEPP2_local_vs2015_64.7z') binary kit; ensure that the absolute path to the 'local_vs2015_64' directory is short (less than 50 characters) and does NOT contain any whitespace (troubles have been encountered with 'Mes Documents' for example)

 8. Extract CMake 3.4.3 from 'path_to\local_vs2015_64\\\_utils_\cmake-3.4.3-win32-x86.zip'

 9. Set a new user environment variable 'MSVC_KIT_ROOT' to 'path_to/local_vs2015_64' (beware of the directory separator, it must be '/' here)

 10. Add ';path_to\local_vs2015_64\\\_bin_' to the PATH system environment variable, just after the Windows system paths, but before any other path, in order to avoid a library version conflict; beware of the ';' path separator

 11. Extract the 'addon_01' binary kit, if needed, into the 'path_to' directory

 12. Extract the 'addon_02' binary kit, if needed, into the 'path_to' directory

 13. Extract the 'addon_03' binary kit, if needed, into the 'path_to' directory

 14. Extract the 'addon_04' binary kit, if needed, into the 'path_to' directory

 15. Extract the 'addon_05' binary kit, if needed, into the 'path_to' directory

### Building stage

 - Get Mepp2 source code using your favourite Git client (see Linux section above)

 - Run cmake-gui.exe

 - Choose 'Visual Studio 14 2015 Win64' as the compiler version

 - Where is the source code = ".../MEPP2"

 - Where to build the binaries = ".../MEPP2/build"

 - Set the CMake compilation flags as needed (see list above)

 - Click "Configure"

 - Click "Generate"

 - Open ".../MEPP2/build/MEPP2.sln" solution with MSVC 2015, select 'Release' mode, then generate the 'ALL_BUILD' target

### Running tests
 - Open ".../MEPP2/build/MEPP2.sln" solution with MSVC 2015.
 - Generate the 'RUN_TESTS' target


## Documentation

The automatically generated documentation is [available online](https://liris.cnrs.fr/mepp/doc/nightly/)


## Known issues

* On Ubuntu: if you encounter segfaults when interactively using the MEPP2 GUI and if you are using open source software implementations of X (X.Org X Server) then consider installing the OpenGL proprietary drivers provided by the constructor of your graphic card (e.g. [NVidia binary drivers](http://www.nvidia.com/object/unix.html)) with this tool:
   ```
   $ software-properties-gtk
   ```

* On Ubuntu 17.04 (04/13/17 Release): when using BUILD_USE_PCL=ON, you must complete your dependency installation stage with the following command:
   ```
   $ sudo ln -s /usr/lib/x86_64-linux-gnu/libmpi.so /usr/lib/libmpi.so
   ```
