# How to build the image ( IMPORTANT, put this line in \Docker\config\daemon.json -> "storage-opts": ["size=40GB"] )
# docker build -t vs2017-ltsc2019-mepp2-vcpkg .

# How to start a container in interactive mode
# docker run -it vs2017-ltsc2019-mepp2-vcpkg

# Use the latest vs2017-ltsc2019 from mtola
FROM mtola/vs2017-ltsc2019

LABEL maintainer="Martial TOLA"

SHELL ["powershell.exe", "-ExecutionPolicy", "Bypass", "-Command"]

# Download and extract a vcpkg release (here 2020.11-1, if necessary and instead of release 2019.12 by default in mtola/vs2017-ltsc2019)
WORKDIR "C:/"
RUN ren vcpkg vcpkg-2019.12
ADD https://github.com/microsoft/vcpkg/archive/refs/tags/2020.11-1.zip C:/TEMP/vcpkg-2020.11-1.zip
RUN 7z x C:/TEMP/vcpkg-2020.11-1.zip
RUN ren vcpkg-2020.11-1 vcpkg
WORKDIR "C:/vcpkg"
RUN .\bootstrap-vcpkg.bat

# Or install last vcpkg (if necessary and instead of release 2019.12 by default in mtola/vs2017-ltsc2019)
#WORKDIR "C:/"
#RUN ren vcpkg vcpkg-2019.12
#RUN git clone https://github.com/Microsoft/vcpkg.git
#WORKDIR "C:/vcpkg"
#RUN .\bootstrap-vcpkg.bat

WORKDIR "C:/vcpkg"

# Installing the vcpkg packages we need
RUN .\vcpkg install boost:x64-windows --clean-after-build
RUN .\vcpkg install eigen3:x64-windows cgal:x64-windows openmesh:x64-windows draco:x64-windows --clean-after-build
RUN .\vcpkg install osg:x64-windows --clean-after-build
RUN .\vcpkg install vtk:x64-windows --clean-after-build
RUN .\vcpkg install pcl:x64-windows --clean-after-build
RUN .\vcpkg install qt5-base:x64-windows --clean-after-build

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

# Temp WORKDIR
WORKDIR "C:/"

# Download and extract Qt5 (LIRIS host)
#ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_02.7z C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
#RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
#RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z

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

#RUN dir "C:/TEMP"
RUN dir "C:/local_vs2015_64"

# Restore WORKDIR
WORKDIR "C:/Users/ContainerAdministrator"

RUN choco install --no-progress --yes cmake --version=3.16.5 --installargs 'ADD_CMAKE_TO_PATH=""System""'
RUN cmake --version

# Restore the default Windows shell for correct batch processing
SHELL ["cmd", "/S", "/C"]

# Create 2 directory junctions
#RUN mklink /j "C:\vcpkg\installed\Qt" "C:\local_vs2015_64\Qt"
RUN mklink /j "C:\vcpkg\installed\FBX_SDK" "C:\local_vs2015_64\FBX_SDK"

RUN dir "C:\vcpkg\installed"

ENV MSVC_KIT_ROOT "C:\vcpkg\installed\x64-windows"
RUN setx path "%path%;%MSVC_KIT_ROOT%\debug\bin;%MSVC_KIT_ROOT%\debug\tools\osg\osgPlugins-3.6.5;%MSVC_KIT_ROOT%\bin;%MSVC_KIT_ROOT%\tools\osg\osgPlugins-3.6.5;C:\local_vs2015_64\_bin_"
#RUN setx path "%path%;C:\local_vs2015_64\Qt\qt-4.8.7-x64-msvc2015\bin"

ENV CMAKE_GENERATOR "Visual Studio 15 2017 Win64"

# Now, define the entry point for the Docker container:
# -----------------------------------------------------
# This entry point starts the developer command prompt and launches the PowerShell shell
ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "powershell.exe", "-NoLogo", "-ExecutionPolicy", "Bypass"]
# This entry point starts the developer command prompt and launches the Command Prompt
#ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2017\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "cmd"]
