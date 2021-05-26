HOW THE FILES IN THIS DIRECTORY HAVE BEEN EXTRACTED FROM CGAL 5.2.1


A - Install CGAL from source

  $ wget https://github.com/CGAL/cgal/releases/download/v5.2.1/CGAL-5.2.1.zip
  $ unzip CGAL-5.2.1.zip

B - Extract Euler operations and OpenMesh graph traits (and their
    dependencies) from CGAL 

  # the commands are run from the MEPP2 directory
 
  $ mkdir External/CGAL-5.2.1
 
  
  ## 1) find the compilation options set by CGAL

  # create a build dir, then configure the build with cmake to use CGAL 5.2.1
  # (something like `cmake CGAL_DIR=$HOME/CGAL-5.2.1 ...`), then trigger
  # a build involving OpenMesh to get the compiler command line

  $ make --dry-run VERBOSE=1 test_calculate_scaling_openmesh  | grep 'c++.*/test_calculate_scaling_openmesh.cpp'

  # take note of all important options


  ## 2) extract Euler operations CGAL dependencies
 
  $ bash ./Tools/list_cgal_dependencies.sh  ~/CGAL-5.2.1/include/CGAL/boost/graph/Euler_operations.h  -DCGAL_EIGEN3_ENABLED -DCGAL_NDEBUG -DCGAL_USE_GMPXX=1 -std=c++14  -isystem ~/CGAL-5.2.1/include  -I ~/OpenMesh-8.1/include
  ...
  - source code analyzed: /home/eric/CGAL-5.2.1/include/CGAL/boost/graph/Euler_operations.h
  - raw dependencies output: Euler_operations.h.depencies.raw
  - dependencies to CGAL: Euler_operations.h.depencies.cgal
  - hierachical dependencies to CGAL: Euler_operations.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
    .../CGAL-5.2.1/include/CGAL/auto_link/auto_link.h
 
  $ wc -l Euler_operations.h.depencies.cgal
  106 Euler_operations.h.depencies.cgal
    # ! /home/eric/CGAL-5.2.1/include/CGAL/boost/graph/Euler_operations.h
    #   appears 2 times in the file
    # so it's actually 105 CGAL dependencies
 
  $ bash ./Tools/import_cgal_headers.sh  Euler_operations.h.depencies.cgal  ./External/CGAL-5.2.1
 
  $ find External/CGAL-5.2.1/CGAL -type f | wc -l
  105
 
 
  ## 3) extract OpenMesh graph traits CGAL dependencies
 
  $ bash ./Tools/list_cgal_dependencies.sh  ~/CGAL-5.2.1/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h  -DCGAL_EIGEN3_ENABLED -DCGAL_NDEBUG -DCGAL_USE_GMPXX=1 -std=c++14  -isystem ~/CGAL-5.2.1/include  -I ~/OpenMesh-8.1/include  -I/usr/include/eigen3
  ...
  - source code analyzed: /home/eric/CGAL-5.2.1/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h
  - raw dependencies output: graph_traits_PolyMesh_ArrayKernelT.h.depencies.raw
  - dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
  - hierachical dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
   .../CGAL-5.2.1/include/CGAL/auto_link/auto_link.h
 
  $ wc -l graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
  470 graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
    # 470 unique CGAL dependencies including all dependencies of Euler_operations.h
 
  $ bash ./Tools/import_cgal_headers.sh  graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal  ./External/CGAL-5.2.1
 
  $ find External/CGAL-5.2.1/CGAL -type f | wc -l
  470
 
 
  ## 4) extract copy_face_graph.h dependencies
 
  $ bash ./Tools/list_cgal_dependencies.sh  ~/CGAL-5.2.1/include/CGAL/boost/graph/copy_face_graph.h   -DCGAL_EIGEN3_ENABLED -DCGAL_NDEBUG -DCGAL_USE_GMPXX=1 -std=c++14  -isystem ~/CGAL-5.2.1/include  -I ~/OpenMesh-8.1/include  -I/usr/include/eigen3
  ...
  - source code analyzed: /home/eric/CGAL-5.2.1/include/CGAL/boost/graph/copy_face_graph.h
  - raw dependencies output: copy_face_graph.h.depencies.raw
  - dependencies to CGAL: copy_face_graph.h.depencies.cgal
  - hierachical dependencies to CGAL: copy_face_graph.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
   .../CGAL-5.2.1/include/CGAL/auto_link/auto_link.h
 
  $ wc -l copy_face_graph.h.depencies.cgal
  114 copy_face_graph.h.depencies.cgal
 
  # there is only one new dependency:
  #    .../CGAL-5.2.1/include/CGAL/boost/graph/copy_face_graph.h
 
  $ bash ./Tools/import_cgal_headers.sh  copy_face_graph.h.depencies.cgal  ./External/CGAL-5.2.1
 
  $ find ./External/CGAL-5.2.1/CGAL -type f | wc -l
  471

  # note: all files are under 'LGPL' except
  #       CGAL-5.2.1/include/CGAL/auto_link/auto_link.h


  ## 5) add  missing file for Windows
 
  $ cp  ~/CGAL-5.2.1/include/CGAL/sse2.h  ./External/CGAL-5.2.1/CGAL
 
  $ find ./External/CGAL-5.2.1/CGAL -type f | wc -l
  472
