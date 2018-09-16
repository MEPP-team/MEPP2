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
#sudo apt-get install libcgal-dev 
cd /tmp
wget https://github.com/CGAL/cgal/archive/releases/CGAL-4.9.tar.gz
tar -xzf CGAL-4.9.tar.gz
cd cgal-releases-CGAL-4.9 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/CGAL-4.9 ..
make
make install
echo -e '\n--->OpenMesh library...\n'
cd /tmp
wget http://www.openmesh.org/media/Releases/6.3/OpenMesh-6.3.tar.gz
tar -xzf OpenMesh-6.3.tar.gz
cd OpenMesh-6.3 && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/OpenMesh-6.3 ..
make
make install
echo -e '\n--->Qt4 library...\n'
sudo apt-get install libqt4-dev
echo -e '\n--->OpenSceneGraph library...\n'
sudo apt-get install libopenscenegraph-dev
echo -e '\n--->Doxygen and Graphviz libraries...\n'
sudo apt-get install doxygen graphviz
echo -e '\n--->VTK library...\n'
sudo apt-get install libvtk5-dev
echo -e '\n--->PCL library...\n'
cd /tmp; wget --no-check-certificate https://download.gforge.liris.cnrs.fr/meppbin/travis-trusty/pcl-pcl-1.7.2.travis-trusty.tgz
$ tar zxf pcl-pcl-1.7.2.travis-trusty.tgz; cd pcl-pcl-1.7.2/buildR
$ sudo make install && cd
# echo -e '\n--->OpenCV library...\n'
# sudo apt-get install libopencv-dev # and then install proposed dependencies
# echo -e '\n\n--Installation finished--\n\n'
./test.sh