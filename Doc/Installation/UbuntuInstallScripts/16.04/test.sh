#!/bin/bash
clear
echo -e '\n\n--Cloning MEPP-team/MEPP2 (you need to have access to this repo to be able to clone it)--\n\n'
sudo rm -rf MEPP2
git clone https://github.com/MEPP-team/MEPP2
cd MEPP2 && mkdir build && cd build
export CGAL_DIR="$HOME/CGAL-4.9/lib/CGAL"
export OPENMESH_DIR="$HOME/OpenMesh-6.3"
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBUILD_USE_VTK=OFF -DBUILD_USE_CGAL=OFF -DBUILD_USE_OPENMESH=ON -DBUILD_USE_AIF=ON -DBUILD_USE_GENERATE_TRIANGLE=ON -DBUILD_USE_PCL=OFF -DBUILD_USE_GUI=OFF -DCMAKE_BUILD_TYPE=Debug ..
core=`cat /proc/cpuinfo | grep processor | wc -l`
echo -e '\n\n--Compiling MEPP2 with gcc (Debug - no VTK - no CGAL - OpenMesh - AIF)--\n\n'
make -j$core 2> build
echo -e '\n\n--Compiling log--\n\n'
vi build
rm build
make test
#read -p "Press [Enter] key to continue..."
cd ..
mkdir buildnovtk
cd buildnovtk
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DBUILD_USE_VTK=ON -DBUILD_USE_CGAL=ON -DBUILD_USE_OPENMESH=ON -DBUILD_USE_AIF=ON -DBUILD_USE_GENERATE_TRIANGLE=ON -DBUILD_USE_PCL=OFF -DBUILD_USE_GUI=OFF -DCMAKE_BUILD_TYPE=Debug ..
#core=`cat /proc/cpuinfo | grep processor | wc -l`
echo -e '\n\n--Compiling MEPP2 with gcc (Debug - VTK - CGAL - OpenMesh - AIF)--\n\n'
make -j$core 2> build
echo -e '\n\n--Compiling log--\n\n'
vi build
rm build
make test
