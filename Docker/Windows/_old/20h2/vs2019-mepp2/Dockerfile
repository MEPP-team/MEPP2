# How to build the image
# docker build -t vs2019-20h2-mepp2 .

# How to start a container in interactive mode
# docker run -it vs2019-20h2-mepp2

# Use the latest vs2019-20h2 from mtola
FROM mtola/vs2019-20h2

LABEL maintainer="Martial TOLA"

SHELL ["powershell.exe", "-ExecutionPolicy", "Bypass", "-Command"]

# Temp WORKDIR
WORKDIR "C:/"

# Download and extract the 'core' binary kit (LIRIS host) that delivers CMake, Doxygen, Graphviz, Boost, CGAL, OpenMesh, Eigen 3 and Img-3rdparty support (jpeg, zlib, png, tiff)
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64-b171_V141_and_V142.7z C:/TEMP/MEPP2_local_vs2015_64-b171_V141_and_V142.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64-b171_V141_and_V142.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64-b171_V141_and_V142.7z

# Download and extract the 'addon_01' binary kit (LIRIS host) for Qt4
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_01.7z C:/TEMP/MEPP2_local_vs2015_64_addon_01.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_01.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_01.7z

# Download and extract the 'addon_02' binary kit (LIRIS host) for Qt5
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_02.7z C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_02.7z

# Download and extract the 'addon_03' binary kit (LIRIS host) for PCL
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_03.7z C:/TEMP/MEPP2_local_vs2015_64_addon_03.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_03.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_03.7z

# Download and extract the 'addon_04' binary kit (LIRIS host) for VTK
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_04.7z C:/TEMP/MEPP2_local_vs2015_64_addon_04.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_04.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_04.7z

# Download and extract the 'addon_05' binary kit (LIRIS host) for FBX
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_05.7z C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_05.7z

# Download and extract the 'addon_06' binary kit (LIRIS host) for Draco
ADD https://download.gforge.liris.cnrs.fr/meppbin/windows/vs2015/MEPP/kits/MEPP2_local_vs2015_64_addon_06.7z C:/TEMP/MEPP2_local_vs2015_64_addon_06.7z
RUN 7z x C:/TEMP/MEPP2_local_vs2015_64_addon_06.7z
RUN del C:/TEMP/MEPP2_local_vs2015_64_addon_06.7z

#RUN dir "C:/TEMP"
RUN dir "C:/local_vs2015_64"

# Restore WORKDIR
WORKDIR "C:/Users/ContainerAdministrator"

RUN choco install --no-progress --yes cmake --version=3.16.5 --installargs 'ADD_CMAKE_TO_PATH=""System""'
RUN cmake --version

# Restore the default Windows shell for correct batch processing
SHELL ["cmd", "/S", "/C"]

ENV MSVC_KIT_ROOT "C:\local_vs2015_64"
RUN setx path "%path%;%MSVC_KIT_ROOT%\_bin_"

ENV CMAKE_GENERATOR "Visual Studio 16 2019"

# Now, define the entry point for the Docker container:
# -----------------------------------------------------
# This entry point starts the developer command prompt and launches the PowerShell shell
ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "powershell.exe", "-NoLogo", "-ExecutionPolicy", "Bypass"]
# This entry point starts the developer command prompt and launches the Command Prompt
#ENTRYPOINT ["C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\BuildTools\\Common7\\Tools\\VsMSBuildCmd.bat", "&&", "cmd"]
