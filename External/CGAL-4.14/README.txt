HOW THE FILES IN THIS DIRECTORY HAVE BEEN EXTRACTED FROM CGAL 4.14


A - Install CGAL from source

  $ wget https://github.com/CGAL/cgal/releases/download/releases%2FCGAL-4.14/CGAL-4.14.zip
  $ unzip CGAL-4.14.zip
  #
  # according to CGAL install documentation here
  #   https://doc.cgal.org/latest/Manual/installation.html#title14
  # it is no more required to configure and build CGAL
  #
  # "7.3.2 Header-only without CMake Configuration
  #
  #  Since CGAL 4.12, CGAL can be used in header-only mode, without even
  #  configuring CGAL. Programs using CGAL (examples, tests, demos, etc.)
  #  must be directly configured using CMake. In this case, CGAL will be
  #  configured at the same time. The variable CGAL_DIR must point to the
  #  root directory of the CGAL source code (either the root of the unpacked
  #  release tarball, or the root of the Git working directory).
  #
  #  So, using CGAL becomes now:
  #
  #  cd /path/to/your/code # go to the directory of the code source using CGAL
  #  cmake -DCGAL_DIR=<CGAL-root> .
  #
  #  7.3.3 CGAL Dependencies
  #
  #  CGAL can be used as a header-only library, though not all its
  #  dependencies are header-only. The libraries Gmp and Mpfr, for
  #  example, are not header-only."
  #


B - Extract Euler operations and OpenMesh graph traits (and their
    dependencies) from CGAL 

  # the commands are run from the MEPP2 directory
 
  $ mkdir External/CGAL-4.14
 

  ## 1) extract Euler operations CGAL dependencies
 
  $ bash ./Tools/list_cgal_dependencies.sh  ~/CGAL-4.14/include/CGAL/boost/graph/Euler_operations.h   -DCGAL_HEADER_ONLY -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.14/include  -I ~/OpenMesh-7.0/include
  ...
  - source code analyzed: /home/eric/CGAL-4.14/include/CGAL/boost/graph/Euler_operations.h
  - raw dependencies output: Euler_operations.h.depencies.raw
  - dependencies to CGAL: Euler_operations.h.depencies.cgal
  - hierachical dependencies to CGAL: Euler_operations.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
    .../CGAL-4.14/include/CGAL/auto_link/auto_link.h
 
  $ wc -l Euler_operations.h.depencies.cgal
  88 Euler_operations.h.depencies.cgal
    # ! /home/eric/CGAL-4.14/include/CGAL/boost/graph/Euler_operations.h
    #   appears 2 times in the file
    # so it's actually 87 CGAL dependencies
 
  $ bash Tools/import_cgal_headers.sh  Euler_operations.h.depencies.cgal  External/CGAL-4.14
 
  $ find External/CGAL-4.14/CGAL -type f | wc -l
  87
 
 
  ## 2) extract OpenMesh graph traits CGAL dependencies
 
  $ bash Tools/list_cgal_dependencies.sh  ~/CGAL-4.14/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h  -DCGAL_HEADER_ONLY -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.14/include -I ~/OpenMesh-7.0/include
  ...
  - source code analyzed: /home/eric/CGAL-4.14/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h
  - raw dependencies output: graph_traits_PolyMesh_ArrayKernelT.h.depencies.raw
  - dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
  - hierachical dependencies to CGAL: graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
   .../CGAL-4.14/include/CGAL/auto_link/auto_link.h
   .../CGAL-4.14/include/CGAL/boost/graph/named_function_params.h
   .../CGAL-4.14/include/CGAL/boost/graph/named_params_helper.h
 
  $ wc -l graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
  443 graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal
    # 443 unique CGAL dependencies including all dependencies of Euler_operations.h
 
  $ Tools/import_cgal_headers.sh  graph_traits_PolyMesh_ArrayKernelT.h.depencies.cgal  External/CGAL-4.14
 
  $ find External/CGAL-4.14/CGAL -type f | wc -l
  443
 
 
  ## 3) extract copy_face_graph.h dependencies
 
  $ bash ../Tools/list_cgal_dependencies.sh  ~/CGAL-4.14/include/CGAL/boost/graph/copy_face_graph.h  -DCGAL_HEADER_ONLY -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.14/include -I ~/OpenMesh-7.0/include
  ...
  - source code analyzed: /home/eric/CGAL-4.14/include/CGAL/boost/graph/copy_face_graph.h
  - raw dependencies output: copy_face_graph.h.depencies.raw
  - dependencies to CGAL: copy_face_graph.h.depencies.cgal
  - hierachical dependencies to CGAL: copy_face_graph.h.depencies.cgal.hierarchical
  - files suspected not to be LGPL:
  /home/eric/CGAL-4.14/include/CGAL/auto_link/auto_link.h
  /home/eric/CGAL-4.14/include/CGAL/boost/graph/named_function_params.h
  /home/eric/CGAL-4.14/include/CGAL/boost/graph/named_params_helper.h
 
  $ wc -l copy_face_graph.h.depencies.cgal
  103 copy_face_graph.h.depencies.cgal
 
  # there is only one new dependency: .../CGAL-4.14/include/CGAL/iterator.h
    (plus ~/CGAL-4.14/include/CGAL/boost/graph/copy_face_graph.h)
 
  $ bash Tools/import_cgal_headers.sh  copy_face_graph.h.depencies.cgal  External/CGAL-4.14
 
  $ find External/CGAL-4.14/CGAL -type f | wc -l
  445

  # note: all files are under 'GNU Lesser General Public License' except
  #       the 3 files
  #         CGAL-4.14/include/CGAL/auto_link/auto_link.h
  #         (already in CGAL 4.11, Boost Software License, Version 1.0)
  #       and
  #         CGAL-4.14/include/CGAL/boost/graph/named_function_params.h
  #         (new in CGAL 4.14, Boost Software License, Version 1.0)
  #       and
  #         CGAL-4.14/include/CGAL/boost/graph/named_params_helper.h
  #         (new in CGAL 4.14, Boost Software License, Version 1.0)


  ## 4) add  missing files for Windows
 
  $ cp  ~/CGAL-4.14/include/CGAL/sse2.h  External/CGAL-4.14/CGAL
  $ cp  ~/CGAL-4.14/include/CGAL/MSVC_compiler_config.h  External/CGAL-4.14/CGAL
 
  $ find External/CGAL-4.14/CGAL -type f | wc -l
  447
