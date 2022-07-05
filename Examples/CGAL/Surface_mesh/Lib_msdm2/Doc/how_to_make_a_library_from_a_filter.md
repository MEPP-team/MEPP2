
# How to make a library from a filter

## General view

A filter is most of the time a function templated relatively to the type of the
mesh. In order to create a library, the filter function must be instanciated
for a specific mesh type. This is done by encapsulating a call to the filter
into a function that uses the specific mesh type.

A complete example is provided in `Examples/CGAL/Surface_mesh/Lib_msdm2`.

## Detailed procedure

1. Instanciate the filter with a specific data structure, aka write a 
   function that take as input a mesh of the chosen data struture and call
   the filter ; the implementation (.cpp) must be separated from the
   header (.h).

   Example: `Examples/CGAL/Surface_mesh/Lib_msdm2/msdm2_surfacemesh.cpp` and
            `Examples/CGAL/Surface_mesh/Lib_msdm2/msdm2_surfacemesh.h`.

2. Compile the instanciated filter function in a dynamic library inside
   MEPP2 environment (add the right lines in a CMakeLists.txt).

   Example: `Examples/CGAL/Surface_mesh/Lib_msdm2/CMakeLists.txt`.

   On Linux a .so file is created.
   On Windows a .lib and a .dll files are created. 

3. Outside of MEPP2, implement an executable that loads a mesh with the chosen
   datastruture (using only the data structure I/O functions) then calls the
   instanciated filter function.

   Example: `Examples/CGAL/Surface_mesh/Lib_msdm2/example_libmsdm2_surfacemesh.cpp`.

4. Copy the instanciated filter function header and the generated library into
   the executable source directory.

5. Compile the executable and link with the dynamic library outside of MEPP2.
   Create a CMakeLists.txt for this purpose.

   Example of CMakeLists.txt:
     ```
     cmake_minimum_required(VERSION 3.1)
     project(Test_libmsdm2_surfacemesh)

     # Set default build type
     if(NOT CMAKE_BUILD_TYPE)
       message(STATUS "Setting build type to 'Release' as none was specified.")
       set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
     endif()

     # find CGAL library
     find_package(CGAL REQUIRED)
     if(${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "4.14")
       message (FATAL_ERROR "CGAL required minimum version is 4.14 - Found: ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION}")
       return()
     endif()
     include(${CGAL_USE_FILE})
     add_definitions(-DCGAL_EIGEN3_ENABLED) # SOME CGAL FUNCTIONALITIES DEPEND ON EIGEN3
     add_definitions(-DBOOST_ALL_DYN_LINK)

     # find Eigen3 library
     find_package(Eigen3 REQUIRED)
     include_directories(${EIGEN3_INCLUDE_DIR})

     # find Boost library
     find_package(Boost REQUIRED)
     include_directories(${Boost_INCLUDE_DIRS})
     add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS) # for Boost >= 1.73, silence a warning about a deprecated use of boost bind

     # find msdm2_surfacemesh library
     find_library(MSDM2LIB msdm2_surfacemesh PATHS ${CMAKE_CURRENT_SOURCE_DIR})
     if(NOT MSDM2LIB)
       message (FATAL_ERROR "msdm2_surfacemesh library required!")
       return()
     endif()
     message(STATUS "found " "${MSDM2LIB}")

     # compile executable
     add_executable(example_libmsdm2_surfacemesh
                    example_libmsdm2_surfacemesh.cpp)
     target_link_libraries(example_libmsdm2_surfacemesh
                           ${MSDM2LIB})
     ```

   Note: on Windows the .lib file generated at step 2 must be copied into the
         executable source directory.

   Example of typical files layout for compilation on Linux:
     ```
     ├── CMakeLists.txt
     ├── example_libmsdm2_surfacemesh.cpp
     ├── msdm2_surfacemesh.h
     └── libmsdm2_surfacemesh.so
     ```

   Example of typical files layout for compilation on Windows:
     ```
     ├── CMakeLists.txt
     ├── example_libmsdm2_surfacemesh.cpp
     ├── msdm2_surfacemesh.h
     ├── msdm2_surfacemesh.lib
     └── msdm2_surfacemesh.dll
     ```

6. Run the executable outside of MEPP2.
   Note: on Windows, the filter .dll must be copied in the same directory as
         the executable.


