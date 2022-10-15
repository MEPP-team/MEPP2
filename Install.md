
## MEPP2 installation guide (Linux Ubuntu, MacOS and Windows)

## For the impatient
CMake compilation flags short list:
 - `CMAKE_BUILD_TYPE   `: build type (Release or Debug)
 - `BUILD_USE_CGAL     `: use CGAL library (ON or OFF)
 - `BUILD_USE_OPENMESH `: use OpenMesh library (ON or OFF)
 - `BUILD_USE_GUI      `: use GUI (ON or OFF)
 - `BUILD_DOCUMENTATION`: build the documentation (ON or OFF)

## Platform specifics:
 - [Linux Ubuntu](#linux-ubuntu-installation) - [Linux Ubuntu Virtual Machine Images](#linux-ubuntu-virtual-machine-images)
 - [MacOS](#macos-installation)
 - [Windows](#windows-installation)

## Dependencies
Mandatory dependencies:
 - CMake >= 3.1
   - Note: [Boost 1.64 requires CMake 3.8 or newer](https://stackoverflow.com/questions/42123509/cmake-finds-boost-but-the-imported-targets-not-available-for-boost-version)
 - Boost >= 1.59
 - Eigen 3

Optional dependencies:
 - CGAL >= 4.14.1: used for CGAL data structures
 - OpenMesh >= 6.2: used for OpenMesh data structures
 - Qt 4 or 5: mandatory for building the GUI
 - OpenSceneGraph >= 3.2: mandatory for building the GUI
 - Doxygen and Graphviz: used to generate the documentation

Other optional dependencies are:
 - VTK
 - PCL
 - FBX
 - Draco

## Linux Ubuntu installation

### Install dependencies

Example on 'Ubuntu 22.04 LTS Jammy Jellyfish (amd64)', released on April 21, 2022:
````
  # Git, GNU C++, Clang
  $ sudo apt install git build-essential clang

  # CMake
  $ sudo apt install cmake

  # Boost
  $ sudo apt install libboost-all-dev

  # Eigen 3
  $ sudo apt install libeigen3-dev

  # CGAL
  $ sudo apt install libgmp-dev libmpfr-dev
  $ sudo apt install libcgal-dev

  # OpenMesh (installation in user home directory)
  $ cd /tmp
  $ wget https://www.graphics.rwth-aachen.de/media/openmesh_static/Releases/9.0/OpenMesh-9.0.tar.gz
  $ tar -xzf OpenMesh-9.0.tar.gz
  $ cd OpenMesh-9.0.0 && mkdir build && cd build
  $ cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-9.0.0 ..
  $ make
  $ make install && cd

  # Qt
  $ sudo apt install qtbase5-dev libqt5opengl5-dev

  # OpenSceneGraph
  $ sudo apt install libjpeg-dev libpng-dev libtiff-dev libfreetype6-dev
  $ sudo apt install libopenscenegraph-dev

  # Doxygen and Graphviz
  $ sudo apt install doxygen graphviz

  # -> Optional dependencies <- :
  # -----------------------------

  # VTK
  $ sudo apt install libvtk9-dev

  # PCL
  $ sudo apt install libflann-dev libproj-dev
  $ sudo apt install libpcl-dev

  # FBX SDK (installation in user home directory, a readme is available when uncompressing the archive) - primary support, not finished !
  $ mkdir ~/FBX_SDK && mkdir ~/FBX_SDK/2019.0
  $ cd /tmp
  $ wget http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_linux.tar.gz
  $ tar -xzf fbx20190_fbxsdk_linux.tar.gz
  $ ./fbx20190_fbxsdk_linux ~/FBX_SDK/2019.0
  $ ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/release ~/FBX_SDK/2019.0/lib/release
  $ ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/debug ~/FBX_SDK/2019.0/lib/debug
  $ cd

  # Draco
  $ sudo apt install libdraco-dev
````

### Build stage

Scripting commands for compiling Mepp2:
````
  # get Mepp2 source code
  $ git clone https://github.com/MEPP-team/MEPP2.git

  # create build directory
  $ cd MEPP2 && mkdir build && cd build

  # compile with CGAL, OpenMesh and GUI
  $ cmake -DOPENMESH_DIR="$HOME/OpenMesh-9.0.0" -DBUILD_USE_GUI=ON -DBUILD_USE_QT5=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # compile without CGAL nor OpenMesh, nor GUI
  $ cmake -DBUILD_USE_CGAL=OFF -DBUILD_USE_OPENMESH=OFF -DBUILD_USE_GUI=OFF -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # compile with CGAL, OpenMesh, GUI and FBX
  $ cmake -DOPENMESH_DIR="$HOME/OpenMesh-9.0.0" -DFBX_DIR="$HOME/FBX_SDK/2019.0" -DBUILD_USE_GUI=ON -DBUILD_USE_QT5=ON -DBUILD_USE_FBX=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make

  # generate the documentation
  $ cmake -DOPENMESH_DIR="$HOME/OpenMesh-9.0.0" -DBUILD_DOCUMENTATION=ON ..
  $ make doc
````

### Run tests

````
  $ ctest
````

### Run the GUI and plugin filters

````
  $ ./Visualization/Applications/mepp-gui
````

## Linux Ubuntu Virtual Machine Images

### VMware

1. Download [VMware Workstation Player for Windows (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/vmware/mepp2/VMware-Player/VMware-player-15.5.1-15018445.exe) or [VMware Workstation Player for Linux (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/vmware/mepp2/VMware-Player/VMware-Player-15.5.1-15018445.x86_64_linux.zip).

2. Download and extract the [MEPP2 Linux Ubuntu 18.04 Virtual Machine Image (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/vmware/mepp2/UbuntuLTS-18.04.3-MEPP.vmwarevm.7z).
````
  # under Linux, if the 7zip utility is not already installed in your system, run this command to install it
  $ sudo apt install p7zip-full

  # then extract the VM with
  $ 7z x UbuntuLTS-18.04.3-MEPP.vmwarevm.7z
````

3. Start the MEPP2 Linux Ubuntu Virtual Machine Image with VMware Workstation Player.

User is '**dev**' and password is '**devdev**'.
Password for **root** is also '**devdev**'.

All the dependencies of MEPP2 are already installed (see installation documentation above) so you can therefore use the following variables with CMake:
- CGAL_DIR="$HOME/CGAL-4.14.1"
- OPENMESH_DIR="$HOME/OpenMesh-7.1"
- OSG_DIR="$HOME/osg-3.6.4"
- PCL_DIR="$HOME/pcl-1.9.1/share/pcl-1.9"
- FBX_DIR="$HOME/FBX_SDK/2019.0"
- DRACO_DIR="$HOME/draco-1.3.5"

### Docker

You can find various Dockerfiles for Linux in /Docker/Linux

## MacOS installation
Various successful installations that we reported:
 - [Mojave (10.14.4)](Doc/Installation/OSX-Mojave-Install_report.md)

### Install dependencies
Mandatory dependencies:
````
  # CAUTION : if you are under OSX Mojave (10.14), you have to upgrade to 10.14.4 (March 25, 2019) at least in order to avoid black blinking video problem

  # Install homebrew package manager (see http://brew.sh/), then
  $ brew update

  # CMake
  $ brew install cmake

  # Boost
  $ brew install boost

  # Eigen 3
  $ brew install eigen
````

Optional dependencies:
````
  # CGAL
  $ brew install cgal

  # OpenMesh
  $ brew install open-mesh

  # Qt 5
  $ brew install qt@5

  # OpenSceneGraph
  $ brew install open-scene-graph

  # Doxygen and Graphviz
  $ brew install doxygen graphviz

  # -> Optional dependencies <- :
  # -----------------------------

  # VTK
  $ brew install vtk

  # PCL
  $ brew install pcl

  # FBX SDK (primary support, not finished !)
  # get and install the sdk from here: https://www.autodesk.com/content/dam/autodesk/www/adn/fbx/20192/fbx20192_fbxsdk_clang_mac.pkg.tgz
  # double click to decompress the .tgz and then install the .pkg by clicking on it
  $ ln -s "/Applications/Autodesk/FBX SDK/2019.2/lib/clang/release" "/Applications/Autodesk/FBX SDK/2019.2/lib/release"
  $ ln -s "/Applications/Autodesk/FBX SDK/2019.2/lib/clang/debug" "/Applications/Autodesk/FBX SDK/2019.2/lib/debug"
  $ export DYLD_LIBRARY_PATH="/Applications/Autodesk/FBX SDK/2019.2/lib/release"

  # Draco
  $ brew install draco
````

### Build stage
For the majority of the commandes refer to the above Linux section. The only difference might be that positionning environment variables is not necessary for OSX (since all dependencies are installed with the packaging system). For example building Mepp2 with CGAL, OpenMesh and GUI is reduced to:
````
  $ cmake -DBUILD_USE_GUI=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

Another example, building with CGAL, OpenMesh, GUI and FBX can be done like this:
````
  $ cmake -DFBX_DIR="/Applications/Autodesk/FBX SDK/2019.2" -DBUILD_USE_GUI=ON -DBUILD_USE_FBX=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

Another example, building with CGAL, OpenMesh, GUI can be done like this:
````
  $ cmake -DBUILD_USE_GUI=ON -DCMAKE_BUILD_TYPE=Release ..
  $ make
````

### Run tests

````
  $ ctest
````

### Run the GUI and plugin filters

````
  $ ./Visualization/Applications/mepp.app/Contents/MacOS/mepp-gui
````

## Windows installation

### Prerequisites

 - Update Windows 10/11 (with Windows Update)

 - Download and install [Visual Studio Community 2017 (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2017/vslayout_2017_fr.7z) (download size 1.8 GB, installation size ~3.4 GB) or [Visual Studio Community 2019 (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2019/vslayout_2019_fr.7z) (download size 2.0 GB, installation size ~4.6 GB) or [Visual Studio Community 2022 * (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2022/vslayout_2022_fr.7z) (download size 2.7 GB, installation size ~5.7 GB). * With Visual Studio Community 2022, you need to upgrade [CMake (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2022/cmake/cmake-3.24.2-windows-x86_64.zip).

### Installing dependencies

 1. Download (mandatory) the ['core' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64-b171_V141_and_V142.7z) that delivers CMake, Doxygen, Graphviz, Boost, CGAL, OpenMesh, Eigen 3 and Img-3rdparty support (jpeg, zlib, png, tiff): download size 876 MB, installation size ~10.9 GB.

 2. Optionally download the ['addon_01' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_01.7z) for Qt4 and OpenSceneGraph
   (download size 396 MB, installation size ~2.9 GB)

 3. Optionally download the ['addon_02' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_02.7z) for Qt5
   (download size 314 MB, installation size ~1.6 GB)

 4. Optionally download the ['addon_03' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_03.7z) for PCL
   (download size 33 MB, installation size ~459 MB)

 5. Optionally download the ['addon_04' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_04.7z) for VTK
   (download size 247 MB, installation size ~2.9 GB)

 6. Optionally download the ['addon_05' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_05.7z) for FBX
   (download size 54 MB, installation size 679 MB)

 7. Optionally download the ['addon_06' binary kit (LIRIS host)](https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_06.7z) for Draco
   (download size 37 MB, installation size 495 MB)

 8. Extract the 'core' ('MEPP2_local_vs2015_64.7z') binary kit; ensure that the absolute path to the 'local_vs2015_64' directory is short (less than 50 characters) and does NOT contain any whitespace (troubles have been encountered with 'Mes Documents' for example)

 9. Extract CMake 3.16.5 from 'path_to\local_vs2015_64\\\_utils_\cmake-3.16.5-win64-x64.zip'

 10. Set a new user environment variable 'MSVC_KIT_ROOT' to 'path_to/local_vs2015_64' (beware of the directory separator, it must be '/' here)

 11. Add ';path_to\local_vs2015_64\\\_bin_' to the PATH system environment variable, just after the Windows system paths, but before any other path, in order to avoid a library version conflict; beware of the ';' path separator

 12. Extract the 'addon_01' binary kit, if needed, into the 'path_to' directory

 13. Extract the 'addon_02' binary kit, if needed, into the 'path_to' directory

 14. Extract the 'addon_03' binary kit, if needed, into the 'path_to' directory

 15. Extract the 'addon_04' binary kit, if needed, into the 'path_to' directory

 16. Extract the 'addon_05' binary kit, if needed, into the 'path_to' directory

 17. Extract the 'addon_06' binary kit, if needed, into the 'path_to' directory

### Building stage

 - Get Mepp2 source code using your favourite Git client (see Linux section above)

 - Run cmake-gui.exe

 - Choose 'Visual Studio 15 2017 Win64/x64' as the compiler version (or 'Visual Studio 16 2019 Win64/x64') -> !!! You have to choose x64 in 'Optional platform for generator' field !!!

 - Where is the source code = ".../MEPP2"

 - Where to build the binaries = ".../MEPP2/build"

 - Set the CMake compilation flags as needed (see list above)

 - Click "Configure"

 - Click "Generate"

 - Open ".../MEPP2/build/MEPP2.sln" solution with MSVC 2015 (or MSVC 2017 or MSVC 2019), select 'Release' mode, then generate the 'ALL_BUILD' target

### Run tests

 - Open ".../MEPP2/build/MEPP2.sln" solution with MSVC 2015 (or MSVC 2017 or MSVC 2019)
 - Generate the 'RUN_TESTS' target

### Run the GUI and plugin filters

 - In Visual Studio 'Solution Explorer', select mepp-gui as the startup project within your solution (right click)
 - Hit Ctrl+F5 to run mepp-gui without debugging (or hit F5 to run mepp-gui in debugging mode)

## Documentation

The automatically generated documentation is [available online](https://liris.cnrs.fr/mepp/doc/nightly/)

## Known issues

* On Ubuntu: if you encounter segfaults when interactively using the MEPP2 GUI and if you are using open source software implementations of X (X.Org X Server) then consider installing the OpenGL proprietary drivers provided by the constructor of your graphic card (e.g. [NVidia binary drivers](http://www.nvidia.com/object/unix.html)) with this tool:
   ```
   $ software-properties-gtk
   ```

* If you are under OSX Mojave (10.14), you have to upgrade to 10.14.4 (March 25, 2019) at least in order to avoid black blinking video problem

## Known (old) issues

* On Ubuntu 17.04 (04/13/17 Release): when using BUILD_USE_PCL=ON, you must complete your dependency installation stage with the following command:
   ```
   $ sudo ln -s /usr/lib/x86_64-linux-gnu/libmpi.so /usr/lib/libmpi.so
   ```
