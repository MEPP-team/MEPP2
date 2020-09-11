# How to building the image
# docker build -t ubuntu-18.04-mepp2 .

# How to start a container in interactive mode
# docker run -it ubuntu-18.04-mepp2

FROM ubuntu:18.04

# Author
MAINTAINER Martial TOLA

# Metadata
LABEL version="1.0"
LABEL vendor="LIRIS"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get clean
RUN apt-get update

RUN apt-get -y install vim
RUN apt-get -y install git

RUN apt-get -y install wget
RUN apt-get -y install unzip

RUN apt-get -y install build-essential gdb cgdb
RUN apt-get -y install clang
RUN apt-get -y install valgrind

RUN apt-get -y install xvfb

# CMake
RUN apt-get -y install cmake cmake-curses-gui

# Boost
RUN apt-get -y install libboost-all-dev

# Eigen 3
RUN apt-get -y install libeigen3-dev

# CGAL (installation in user home directory)
RUN apt-get -y install libgmp-dev libmpfr-dev
RUN cd /tmp                                                                                    && \
    wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14.3/CGAL-4.14.3.zip && \
    cd && unzip /tmp/CGAL-4.14.3.zip
RUN cd /tmp                                                                                    && \
    wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-5.0.3/CGAL-5.0.3.zip   && \
    cd && unzip /tmp/CGAL-5.0.3.zip

# OpenMesh (installation in user home directory)
RUN cd /tmp                                                                       && \
    wget https://www.openmesh.org/media/Releases/7.1/OpenMesh-7.1.tar.gz          && \
    tar -xzf OpenMesh-7.1.tar.gz                                                  && \
    cd OpenMesh-7.1 && mkdir build && cd build                                    && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-7.1 .. && \
    make -j 3                                                                     && \
    make install && cd
RUN cd /tmp                                                                       && \
    wget https://www.openmesh.org/media/Releases/8.1/OpenMesh-8.1.tar.gz          && \
    tar -xzf OpenMesh-8.1.tar.gz                                                  && \
    cd OpenMesh-8.1 && mkdir build && cd build                                    && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-8.1 .. && \
    make -j 3                                                                     && \
    make install && cd

# Qt4
RUN apt-get -y install libqt4-dev libqt4-opengl-dev
# Qt5
RUN apt-get -y install qtdeclarative5-dev libqt5opengl5-dev

# OpenSceneGraph
RUN apt-get -y install libjpeg-dev libpng-dev libtiff-dev libfreetype6-dev
RUN apt-get -y install libopenscenegraph-3.4-dev

# Doxygen and Graphviz
RUN apt-get -y install doxygen graphviz

# VTK
RUN apt-get -y install libvtk6-dev
#RUN apt-get -y install libvtk7-dev

# PCL (installation in user home directory) - DON'T use libpcl-dev package !
RUN apt-get -y install libflann-dev libproj-dev
RUN cd /tmp                                                                                                                                                                && \
    wget https://github.com/PointCloudLibrary/pcl/archive/pcl-1.9.1.tar.gz                                                                                                 && \
    tar -xzf pcl-1.9.1.tar.gz                                                                                                                                              && \
    cd pcl-pcl-1.9.1 && mkdir build && cd build                                                                                                                            && \
    cmake -DCMAKE_BUILD_TYPE=Release -DPCL_ENABLE_SSE=OFF -DPCL_ONLY_CORE_POINT_TYPES=ON -DBUILD_global_tests=OFF -DWITH_VTK=OFF -DCMAKE_INSTALL_PREFIX=$HOME/pcl-1.9.1 .. && \
    make -j 3                                                                                                                                                              && \
    make install && cd

# FBX SDK
RUN mkdir ~/FBX_SDK && mkdir ~/FBX_SDK/2019.0                                         && \
    cd /tmp                                                                           && \
    wget http://download.autodesk.com/us/fbx/2019/2019.0/fbx20190_fbxsdk_linux.tar.gz && \
    tar -xzf fbx20190_fbxsdk_linux.tar.gz                                             && \
    yes yes | ./fbx20190_fbxsdk_linux ~/FBX_SDK/2019.0                                && \
    ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/release ~/FBX_SDK/2019.0/lib/release          && \
    ln -s ~/FBX_SDK/2019.0/lib/gcc4/x64/debug ~/FBX_SDK/2019.0/lib/debug              && \
    cd

# Draco (installation in user home directory)
RUN cd /tmp                                                                                             && \
    wget https://github.com/google/draco/archive/1.3.6.tar.gz                                           && \
    tar -xzf 1.3.6.tar.gz                                                                               && \
    cd draco-1.3.6 && mkdir build && cd build                                                           && \
    cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$HOME/draco-1.3.6 .. && \
    make -j 3                                                                                           && \
    make install && cd

# Cleanup
RUN apt-get -y autoremove

ENV HOME /root

# App environment
ENV CGAL_DIR "$HOME/CGAL-4.14.3"
#ENV CGAL_DIR "$HOME/CGAL-5.0.3"
#ENV OPENMESH_DIR "$HOME/OpenMesh-7.1"
ENV OPENMESH_DIR "$HOME/OpenMesh-8.1"
ENV PCL_DIR "$HOME/pcl-1.9.1/share/pcl-1.9"
ENV FBX_DIR "$HOME/FBX_SDK/2019.0"
ENV DRACO_DIR "$HOME/draco-1.3.6"
