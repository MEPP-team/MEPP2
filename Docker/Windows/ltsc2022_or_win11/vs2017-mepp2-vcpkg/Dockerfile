# How to build the image ( IMPORTANT, put this line in \Docker\config\daemon.json -> "storage-opts": ["size=50GB"] )
# docker build -t vs2017-ltsc2022-mepp2-vcpkg .

# How to start a container in interactive mode
# docker run -it vs2017-ltsc2022-mepp2-vcpkg

# Use the latest vs2017-ltsc2022 from mtola
FROM mtola/vs2017-ltsc2022

LABEL maintainer="Martial TOLA"

SHELL ["powershell.exe", "-ExecutionPolicy", "Bypass", "-Command"]

# Download and extract a vcpkg release (here 2021.05.12, if necessary and instead of release 2019.12 by default in mtola/vs2017-ltsc2022)
WORKDIR "C:/"
RUN ren vcpkg vcpkg-2019.12
ADD https://github.com/microsoft/vcpkg/archive/refs/tags/2021.05.12.zip C:/TEMP/vcpkg-2021.05.12.zip
RUN 7z x C:/TEMP/vcpkg-2021.05.12.zip
RUN ren vcpkg-2021.05.12 vcpkg
WORKDIR "C:/vcpkg"
RUN .\bootstrap-vcpkg.bat -disableMetrics

# Or install last vcpkg (if necessary and instead of release 2019.12 by default in mtola/vs2017-ltsc2022)
#WORKDIR "C:/"
#RUN ren vcpkg vcpkg-2019.12
# Here it is really IMPORTANT to use your own fork of vcpkg (from https://github.com/microsoft/vcpkg) in order to 'freeze' the version of vcpkg used !
#RUN git clone https://github.com/mtola/vcpkg.git
#WORKDIR "C:/vcpkg"
#RUN .\bootstrap-vcpkg.bat -disableMetrics

WORKDIR "C:/vcpkg"

# Installing the vcpkg packages we need
#RUN .\vcpkg integrate install
RUN .\vcpkg install vtk:x64-windows --clean-after-build # HERE or LATER below... [-img2]
RUN .\vcpkg install boost:x64-windows # DON'T WORK with --clean-after-build !
RUN .\vcpkg install eigen3:x64-windows openmesh:x64-windows draco:x64-windows --clean-after-build

# gmp (AND REMOVE with a json patch the x64-windows yasm-tool dependency !)
RUN .\vcpkg install yasm-tool:x86-windows --clean-after-build # IMPORTANT, x86-windows HERE BECAUSE 'Cross-targetting and x64 ports requiring yasm require the x86-windows yasm-tool' !
RUN dir "C:/vcpkg/ports/gmp"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2021.05.12/ports/gmp/vcpkg.json C:/vcpkg/ports/gmp 
RUN dir "C:/vcpkg/ports/gmp"
RUN .\vcpkg install gmp:x64-windows --clean-after-build
# cgal
RUN .\vcpkg install cgal:x64-windows --clean-after-build
# osg
RUN .\vcpkg install expat:x64-windows --clean-after-build
RUN .\vcpkg install libpng:x64-windows --clean-after-build
RUN .\vcpkg install brotli:x64-windows --clean-after-build
RUN .\vcpkg install freetype:x64-windows --clean-after-build
RUN .\vcpkg install osg:x64-windows --clean-after-build
# vtk (IF NOT DONE above...) [-img1]
#RUN .\vcpkg install pugixml:x64-windows --clean-after-build
#RUN .\vcpkg install vtk:x64-windows --recurse
# pcl
RUN .\vcpkg install pcl:x64-windows --clean-after-build
# qt5-base (OFF, BECAUSE qt5-base 5.15.2 CANNOT COMPILE under ltsc2022 !)
#RUN .\vcpkg install egl-registry:x64-windows ragel:x64-windows harfbuzz:x64-windows --clean-after-build
#RUN .\vcpkg install icu:x64-windows --clean-after-build
#RUN .\vcpkg install libpq:x64-windows --clean-after-build
#RUN .\vcpkg install pcre2:x64-windows --clean-after-build
#RUN .\vcpkg install qt5-base:x64-windows --clean-after-build

RUN dir "C:/vcpkg/installed/vcpkg/info"

# Patch for CGAL from VCPKG (LIRIS host - HIDE warnings with GMP in Gmpfr_type.h and Gmpzf_type.h)
WORKDIR "C:/vcpkg/installed/x64-windows/include/CGAL/GMP"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2020.11-1/GMP_for_installed.rar C:/TEMP/GMP_for_installed.rar
RUN 7z x -y C:/TEMP/GMP_for_installed.rar
RUN del C:/TEMP/GMP_for_installed.rar

# OFF/NO_NEED with vcpkg-2020.11-1 ! - Patch for PCL from VCPKG (LIRIS host - FIX problem with flann::flann_cpp library in PCLConfig.cmake)
#WORKDIR "C:/vcpkg/installed/x64-windows/share/pcl"
#ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2019.12/PCLConfig_for_installed.rar C:/TEMP/PCLConfig_for_installed.rar
#RUN 7z x -y C:/TEMP/PCLConfig_for_installed.rar
#RUN del C:/TEMP/PCLConfig_for_installed.rar

# Patch for PCL from VCPKG (LIRIS host - FIX problem with PCL_DEPRECATED(1, 12, "rewrite your code to avoid using this protected field") shared_ptr<MsgFieldMap> mapping_; in point_cloud.h)
WORKDIR "C:/vcpkg/installed/x64-windows/include/pcl"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2020.11-1/point_cloud_for_installed.rar C:/TEMP/point_cloud_for_installed.rar
RUN 7z x -y C:/TEMP/point_cloud_for_installed.rar
RUN del C:/TEMP/point_cloud_for_installed.rar

# Patch for PCL (LIRIS host - HIDE 'BOOST_PARAMETER_MAX_ARITY' warnings with PCL in boost.h)
WORKDIR "C:/vcpkg/installed/x64-windows/include/pcl/io"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2020.11-1/PCL-pcl-1.11.1__include__pcl-1.11__pcl__io__boost.h.rar C:/TEMP/PCL-pcl-1.11.1__include__pcl-1.11__pcl__io__boost.h.rar
RUN 7z x -y C:/TEMP/PCL-pcl-1.11.1__include__pcl-1.11__pcl__io__boost.h.rar
RUN del C:/TEMP/PCL-pcl-1.11.1__include__pcl-1.11__pcl__io__boost.h.rar

# Patch for VTK from VCPKG (LIRIS host - cmake patches)
WORKDIR "C:/vcpkg/installed/x64-windows/share/vtk"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2020.11-1/vtk_for_installed.rar C:/TEMP/vtk_for_installed.rar
RUN 7z x -y C:/TEMP/vtk_for_installed.rar
RUN del C:/TEMP/vtk_for_installed.rar
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2021.05.12/VTK-vtk-module-find-packages_for_installed.rar C:/TEMP/VTK-vtk-module-find-packages_for_installed.rar
RUN 7z x -y C:/TEMP/VTK-vtk-module-find-packages_for_installed.rar
RUN del C:/TEMP/VTK-vtk-module-find-packages_for_installed.rar
# Patch for netcdf (VTK module) from VCPKG (LIRIS host - cmake patches)
WORKDIR "C:/vcpkg/installed/x64-windows/share/netcdf"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2020.11-1/netcdf-c_for_installed.rar C:/TEMP/netcdf-c_for_installed.rar
RUN 7z x -y C:/TEMP/netcdf-c_for_installed.rar
RUN del C:/TEMP/netcdf-c_for_installed.rar

# Patch for gmp/mpir from VCPKG (LIRIS host - bin patches)
WORKDIR "C:/vcpkg/installed/x64-windows/lib"
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/v14x-x64/tags/patchs/vcpkg-2021.05.12/mpir.lib_in_fact_it_is_a_copy_of_gmp.lib.rar C:/TEMP/mpir.lib_in_fact_it_is_a_copy_of_gmp.lib.rar
RUN 7z x -y C:/TEMP/mpir.lib_in_fact_it_is_a_copy_of_gmp.lib.rar
RUN del C:/TEMP/mpir.lib_in_fact_it_is_a_copy_of_gmp.lib.rar

# Temp WORKDIR
WORKDIR "C:/"

# Download and extract Qt5 (LIRIS host)
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_02.7z C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z

# Download and extract FBX (LIRIS host)
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_05.7z C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z

# Temp WORKDIR
WORKDIR "C:/local_vs2015_64"

# Download and extract Qt4 (LIRIS host) - optional
#ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/packages/clean/qt-4.8.7-x64-msvc2015-moc_patched.7z C:/TEMP/qt-4.8.7-x64-msvc2015-moc_patched.7z
#RUN 7z x C:/TEMP/qt-4.8.7-x64-msvc2015-moc_patched.7z
#RUN del C:/TEMP/qt-4.8.7-x64-msvc2015-moc_patched.7z

#RUN dir "C:/TEMP" # OFF, ONLY for test IF necessary...
RUN dir "C:/local_vs2015_64"

# Restore WORKDIR
WORKDIR "C:/Users/ContainerAdministrator"

RUN choco install --no-progress --yes cmake --version=3.16.5 --installargs 'ADD_CMAKE_TO_PATH=""System""'
RUN cmake --version

# Restore the default Windows shell for correct batch processing
SHELL ["cmd", "/S", "/C"]

# Create directory junctions
# Qt5, OFF IF we have INSTALLED qt5-base:x64-windows FROM vcpkg !
RUN mklink /j "C:\vcpkg\installed\Qt" "C:\local_vs2015_64\Qt"
RUN mkdir "C:\vcpkg\installed\Qt\Qt5.6.3\5.6.3\msvc2015_64\debug"
RUN mklink /j "C:\vcpkg\installed\Qt\Qt5.6.3\5.6.3\msvc2015_64\debug\plugins" "C:\vcpkg\installed\Qt\Qt5.6.3\5.6.3\msvc2015_64\plugins"
# FBX
RUN mklink /j "C:\vcpkg\installed\FBX_SDK" "C:\local_vs2015_64\FBX_SDK"

RUN dir "C:\vcpkg\installed"

ENV MSVC_KIT_ROOT "C:\vcpkg\installed\x64-windows"
RUN setx path "%path%;%MSVC_KIT_ROOT%\debug\bin;%MSVC_KIT_ROOT%\debug\tools\osg\osgPlugins-3.6.5;%MSVC_KIT_ROOT%\bin;%MSVC_KIT_ROOT%\tools\osg\osgPlugins-3.6.5;C:\local_vs2015_64\_bin_"
#RUN setx path "%path%;C:\local_vs2015_64\Qt\qt-4.8.7-x64-msvc2015\bin"

RUN git --version

ENV CMAKE_GENERATOR "Visual Studio 15 2017 Win64"

# Now, define the entry point for the Docker container:
# -----------------------------------------------------
# This entry point starts the developer command prompt and launches the PowerShell shell
ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "powershell.exe", "-NoLogo", "-ExecutionPolicy", "Bypass"]
# This entry point starts the developer command prompt and launches the Command Prompt
#ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "cmd"]
