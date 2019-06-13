#!/bin/bash
clear
echo -e '\n\n--Upgrading operating system--\n\n'
sudo apt-get update -qq
sudo apt-get upgrade -qq
echo -e "\n\nInstalling essentials (git, compilers and cmake)--\n\n"
sudo apt-get install git build-essential clang cmake synaptic aptitude
echo -e '\n\n--Installing libraries--\n\n'
echo -e '\n--->boost libraries...\n'
sudo apt-get install libboost-program-options-dev libboost-dev libboost-thread-dev libboost-system-dev libboost-filesystem-dev 
echo -e '\n--->eigen library...\n'
sudo apt-get install libeigen3-dev
echo -e '\n--->CGAL library...\n'
sudo apt purge libcgal-dev
cd /tmp
wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14/CGAL-4.14.zip
cd && unzip /tmp/CGAL-4.14.zip
sudo apt install libgmp-dev libmpfr-dev
echo -e '\n--->OpenMesh library...\n'
cd /tmp
wget http://www.openmesh.org/media/Releases/7.0/OpenMesh-7.0.tar.gz
tar -xzf OpenMesh-7.0.tar.gz
cd OpenMesh-7.0 && mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-7.0 ..
make
make install && cd

echo -e '\n--->Qt4 library...\n'
sudo apt-get install libqt4-dev
echo -e '\n--->OpenSceneGraph library...\n'
sudo apt-get install libopenscenegraph-3.4-dev
echo -e '\n--->Doxygen and Graphviz libraries...\n'
sudo apt-get install doxygen graphviz
echo -e '\n--->VTK library...\n'
sudo apt-get install libvtk6-dev
echo -e '\n--->PCL library...\n'
sudo apt install libflann-dev
cd /tmp
wget https://github.com/PointCloudLibrary/pcl/archive/pcl-1.8.1.tar.gz
tar -xzf pcl-1.8.1.tar.gz
cd pcl-pcl-1.8.1 && mkdir build 
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DPCL_ONLY_CORE_POINT_TYPES=ON -DBUILD_global_tests=OFF -DCMAKE_INSTALL_PREFIX=$HOME/pcl-1.8.1 ..
make
make install && cd
# echo -e '\n--->OpenCV library...\n'
# sudo apt-get install libopencv-dev # and then install proposed dependencies
# echo -e '\n\n--Installation finished--\n\n'
./test.sh