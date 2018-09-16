HOW THE FILES IN THIS DIRECTORY HAVE BEEN EXTRACTED FROM CGAL 4.11


A - Install CGAL from source

 $ wget https://github.com/CGAL/cgal/archive/releases/CGAL-4.11.tar.gz
 $ tar -xzf CGAL-4.11.tar.gz
 $ mv cgal-releases-CGAL-4.11 CGAL-4.11-src
 $ cd CGAL-4.11-src
 $ mkdir build && cd build
 $ cmake -DCMAKE_INSTALL_PREFIX="$HOME/CGAL-4.11" ..
 $ make
 $ make install


B - Extract Euler operations and OpenMesh graph traits (and their
    dependencies) from CGAL 

 # the commands are run from the MEPP2 directory

 $ mkdir External/CGAL-4.11

 # 1) extract Euler operations CGAL dependencies

 $ Tools/list_cgal_dependencies.sh  ~/CGAL-4.11/include/CGAL/boost/graph/Euler_operations.h  -DCGAL_USE_GMP -DCGAL_USE_MPFR -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.11/include -I ~/OpenMesh-5.1/include
 ...
 - source code analyzed: /home/eric/CGAL-4.11/include/CGAL/boost/graph/Euler_operations.h
 - raw dependencies output: Euler_operations.h.depencies.raw
 - dependencies to CGAL: Euler_operations.h.depencies.cgal
 - hierachical dependencies to CGAL: Euler_operations.h.depencies.cgal.hierarchical
 - files suspected not to be LGPL:
 .../CGAL-4.11/include/CGAL/compiler_config.h
 .../CGAL-4.11/include/CGAL/auto_link/auto_link.h

 $ wc -l Euler_operations.h.depencies.cgal
 74 Euler_operations.h.depencies.cgal
   # ! /home/eric/CGAL-4.11/include/CGAL/boost/graph/Euler_operations.h appears 2 times
   #   in the file
   # so it's actually 73 CGAL dependencies

 $ Tools/import_cgal_headers.sh  Euler_operations.h.depencies.cgal  External/CGAL-4.11

 $ find External/CGAL-4.11/CGAL -type f | wc -l
 73


 # 2) extract OpenMesh graph traits CGAL dependencies

 $ Tools/list_cgal_dependencies.sh  ~/CGAL-4.11/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h  -DCGAL_USE_GMP -DCGAL_USE_MPFR -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.11/include -I ~/OpenMesh-5.1/include
 ...
 - source code analyzed: /home/eric/CGAL-4.11/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h
 - raw dependencies output: graph_traits_PolyMesh_ArrayKernelT.h.depencies.raw
 - dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
 - hierachical dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal.hierarchical
 - files suspected not to be LGPL:
 .../CGAL-4.11/include/CGAL/compiler_config.h
 .../CGAL-4.11/include/CGAL/auto_link/auto_link.h

 $ wc -l graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
 363 graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
   # 363 CGAL dependencies including all dependencies of Euler_operations.h

 $ Tools/import_cgal_headers.sh  graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal  External/CGAL-4.11

 $ find External/CGAL-4.11/CGAL -type f | wc -l
 363


 # 3) extract copy_face_graph.h dependencies

 $ ../Tools/list_cgal_dependencies.sh  ~/CGAL-4.11/include/CGAL/boost/graph/copy_face_graph.h  -DCGAL_USE_GMP -DCGAL_USE_MPFR -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.11/include -I ~/OpenMesh-5.1/include
 ...
 - source code analyzed: /home/eric/CGAL-4.11/include/CGAL/boost/graph/copy_face_graph.h
 - raw dependencies output: copy_face_graph.h.depencies.raw
 - dependencies to CGAL: copy_face_graph.h.depencies.cgal
 - hierachical dependencies to CGAL: copy_face_graph.h.depencies.cgal.hierarchical
 - files suspected not to be LGPL:
 .../CGAL-4.11/include/CGAL/compiler_config.h
 .../CGAL-4.11/include/CGAL/auto_link/auto_link.h

 # there is only one new dependency: .../CGAL-4.11/include/CGAL/iterator.h

 $ Tools/import_cgal_headers.sh  copy_face_graph.h.depencies.cgal  External/CGAL-4.11

 $ find External/CGAL-4.11/CGAL -type f | wc -l
 365



 # note: all files are under 'GNU Lesser General Public License' except the 2 files
 #           CGAL-4.11/include/CGAL/compiler_config.h
 #       and
 #           CGAL-4.11/include/CGAL/auto_link/auto_link.h



 # 4) disable GMP and MPFR to avoid extra include dependencies

 $ vim External/CGAL-4.11/CGAL/compiler_config.h
   # comment line "#define CGAL_USE_GMP 1"
   # and line "#define CGAL_USE_MPFR 1"
 

 # 5) add one missing file for Windows

 $ cp  ~/CGAL-4.11/include/CGAL/sse2.h  External/CGAL-4.11/CGAL

 $ find External/CGAL-4.11/CGAL -type f | wc -l
 364
