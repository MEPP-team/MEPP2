if( NOT BUILD_USE_AIF )
  return()
endif()

if( WIN32 )
  # AUTOLINK is used for CGAL on Windows. Disable it for AIF tests in
  # order to prevent false link error.
  add_definitions( -DCGAL_NO_AUTOLINK_CGAL )
endif()

# Refer to DevelopersDesignCMakeFiles.dox entry 001
if( BUILD_USE_CGAL )
  include(${CGAL_USE_FILE})
else()
  include_directories( BEFORE "${PROJECT_SOURCE_DIR}/External/CGAL-4.14" )
endif()

# Refer to Doc/Devel/CMakeFiles.md entry 002
add_definitions( -DCGAL_NDEBUG )

##################################################################""
add_executable( test_AIF
                test_AIF.cpp)
target_link_libraries(test_AIF
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

#####
add_executable( test_geometry_traits_concept_aif
                test_geometry_traits_concept_aif.cpp
                )
target_link_libraries(test_geometry_traits_concept_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

#####
add_executable( test_geometry_traits_unit_normal_aif
                test_geometry_traits_unit_normal_aif.cpp)
target_link_libraries(test_geometry_traits_unit_normal_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test( Geometry_traits_unit_normal_aif
	test_geometry_traits_unit_normal_aif )

#####
add_executable( test_helpers_aif
                test_helpers_aif.cpp)
target_link_libraries(test_helpers_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(AIFHelpers
         test_helpers_aif )

#####
add_executable( test_reader_writer_aif
                test_reader_writer_aif.cpp
				        )
target_link_libraries(test_reader_writer_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( AIFReader_and_Writer
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/CubeTriangleFaces.obj
  test_AIFReader_and_Writer.obj
  ${TESTING_DATA_DIR}/CubeTriangleFaces.writer-output.obj.ref
)
add_test( AIFReader_and_Writer_obj_with_normals
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.obj
  test_AIFReader_and_Writer_obj_with_normals.obj
  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.writer-output.obj.ref
)
add_test( AIFReader_and_Writer_obj_with_normals_and_color
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.obj
  test_AIFReader_and_Writer_obj_with_normals_and_color.obj
  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.writer-output.obj.ref
)
#ELO: disable the following test because texture support in AIFReader and
#     AIFWriter is now obsolete:
#     - AIFReader support vertex-texture with OBJ, which is obsolete
#     - AIFWriter does NOT support corner-texture at all
#add_test( AIFReader_and_Writer_obj_with_normals_and_texture
#  test_reader_writer_aif
#  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.obj
#  test_AIFReader_and_Writer_obj_with_normals_and_texture.obj
#  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.writer-output.obj.ref
#)
add_test( AIFReader_and_Writer_off_flat
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/flat.off
  test_AIFReader_and_Writer_off_flat.off
  ${TESTING_DATA_DIR}/flat.off
)
add_test( AIFReader_and_Writer_off_aircraft
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/aircraft.off
  test_AIFReader_and_Writer_off_aircraft.off
  ${TESTING_DATA_DIR}/aircraft.off
)
add_test( AIFReader_and_Writer_off_with_face-color
  test_reader_writer_aif
  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off
  test_AIFReader_and_Writer_off_with_face.off
  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off
)
#####
add_executable( test_property_map_aif
                test_property_map_aif.cpp )
target_link_libraries(test_property_map_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_property_map_aif
          test_property_map_aif)

#####
add_executable( test_boost_graph_concept_aif
                test_boost_graph_concept_aif.cpp )
target_link_libraries(test_boost_graph_concept_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

#####
add_executable( test_halfedge_aif
                test_halfedge_aif.cpp)
target_link_libraries(test_halfedge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_halfedge_aif
          test_halfedge_aif )

#####
add_executable( test_calculate_scaling_aif
                test_calculate_scaling_aif.cpp)
target_link_libraries(test_calculate_scaling_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test( Calculate_scaling_aif
         test_calculate_scaling_aif
         ${TESTING_DATA_DIR}/flat.off )
#####
add_executable( test_comparison_operators_aif
                test_comparison_operators_aif.cpp)
target_link_libraries(test_comparison_operators_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_comparison_operators_aif
          test_comparison_operators_aif )

#####
add_executable( test_one_ring_aif
                test_one_ring_aif.cpp)
target_link_libraries(test_one_ring_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_one_ring_aif
          test_one_ring_aif )

#####
add_executable( test_complex_one_ring_aif
                test_complex_one_ring_aif.cpp)
target_link_libraries(test_complex_one_ring_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_complex_one_ring_aif
          test_complex_one_ring_aif )

#####
add_executable( test_matrix_operations_aif
                test_matrix_operations_aif.cpp)
target_link_libraries(test_matrix_operations_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_matrix_operations_aif
          test_matrix_operations_aif )

#####
add_executable( test_not_2_manifold_aif
                test_not_2_manifold_aif.cpp)
target_link_libraries(test_not_2_manifold_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

add_test( Test_test_not_2_manifold_aif
          test_not_2_manifold_aif )

#####
add_executable( test_collapse_edge_aif
                test_collapse_edge_aif.cpp)
target_link_libraries(test_collapse_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Collapse_edge_aif_1
         test_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off )
add_test(Collapse_edge_aif_2
         test_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )

#####
add_executable( test_collapse_edge_Euler_aif
                test_collapse_edge_euler_aif.cpp)
target_link_libraries(test_collapse_edge_Euler_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Collapse_edge_Euler_aif_1
         test_collapse_edge_Euler_aif
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )
         # The following test fails with collapse_edge_34_keep_target.off,
         # because there is no guarantee on the selected halfedge when using
         # the halfedge(e,g) operation.
add_test(Collapse_edge_Euler_aif_2
         test_collapse_edge_Euler_aif
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )

#####
add_executable( test_flip_edge_aif
                test_flip_edge_aif.cpp
                )
target_link_libraries(test_flip_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Flip_edge_aif_1
         test_flip_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_flip_0.off
         ${TESTING_DATA_DIR}/square3.off )
add_test(Flip_edge_aif_2
         test_flip_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_flip_1.off
         ${TESTING_DATA_DIR}/square3.off )
#####
add_executable( test_split_edge_aif
                test_split_edge_aif.cpp
                )
target_link_libraries(test_split_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Split_edge_aif_1
         test_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         0
         flat_0.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_aif_2
         test_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         1
         flat_1.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_aif_3
         test_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         2
         flat_2.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_aif_4
         test_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         3
         flat_3.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_aif_5
         test_split_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_split_0.off
         ${TESTING_DATA_DIR}/square4.off )
add_test(Split_edge_aif_6
         test_split_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_split_1.off
         ${TESTING_DATA_DIR}/square4.off )
#####
add_executable( test_non_manifold_split_edge_aif
                test_non_manifold_split_edge_aif.cpp)
target_link_libraries(test_non_manifold_split_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Non_manifold_split_edge_aif_1
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         0
         flat_nm_0.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Non_manifold_split_edge_aif_2
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         1
         flat_nm_1.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Non_manifold_split_edge_aif_3
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         2
         flat_nm_2.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Non_manifold_split_edge_aif_4
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         3
         flat_nm_3.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Non_manifold_split_edge_aif_5
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/cubeNonMan.off
         4
         cubeNonMan_split_0.off
         ${TESTING_DATA_DIR}/cubeNonMan_57_or_75.off )
add_test(Non_manifold_split_edge_aif_6
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/cubeNonMan.off
         5
         cubeNonMan_split_1.off
         ${TESTING_DATA_DIR}/cubeNonMan_57_or_75.off )
add_test(Non_manifold_split_edge_aif_7
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_nm_split_0.off
         ${TESTING_DATA_DIR}/square4.off )
add_test(Non_manifold_split_edge_aif_8
         test_non_manifold_split_edge_aif
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_nm_split_1.off
         ${TESTING_DATA_DIR}/square4.off )
#####
add_executable( test_non_manifold_collapse_edge_aif
                test_non_manifold_collapse_edge_aif.cpp)
target_link_libraries(test_non_manifold_collapse_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test(Non_manifold_collapse_edge_aif_1
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         0
         flat_col_0.off
         ${TESTING_DATA_DIR}/collapse_edge_12_keep_target.off)
add_test(Non_manifold_collapse_edge_aif_2
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         1
         flat_col_1.off
         ${TESTING_DATA_DIR}/collapse_edge_21_keep_target.off)
add_test(Non_manifold_collapse_edge_aif_3
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         2
         flat_col_2.off
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off
         )
add_test(Non_manifold_collapse_edge_aif_4
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/flat.off
         3
         flat_col_3.off
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off
         )
add_test(Non_manifold_collapse_edge_aif_5
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/cubeNonMan.off
         4
         cubeNonMan_col_0.off)
add_test(Non_manifold_collapse_edge_aif_6
         test_non_manifold_collapse_edge_aif
         ${TESTING_DATA_DIR}/cubeNonMan.off
         5
         cubeNonMan_col_1.off)
#####
add_executable( test_smoothing_aif
                test_smoothing_aif.cpp)

target_link_libraries(test_smoothing_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

if (BUILD_USE_VTK )
	add_test( Smoothing_aif_vtk
			 test_smoothing_aif
			 ${TESTING_DATA_DIR}/cube.vtk
       #${TESTING_DATA_DIR}/egea.vtk # time consuming test
       )
endif()

add_test( Smoothing_aif
         test_smoothing_aif
         ${TESTING_DATA_DIR}/cow_simp.off)
add_test( Smoothing_aif_border
         test_smoothing_aif
         ${TESTING_DATA_DIR}/flat.off )
add_test( Smoothing_aif_gmesh
         test_smoothing_aif
         ${TESTING_DATA_DIR}/curves/SimpleBeam0000016.msh )
#####
add_executable( test_curvature_filter_aif
                test_curvature_filter_aif.cpp
                )
target_link_libraries(test_curvature_filter_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

### Checking of concepts is done at build time (no add_test required)
add_executable( test_complying_concepts_aif
                test_complying_concepts_aif.cpp
				)

target_link_libraries(test_complying_concepts_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

#####
add_executable( test_generic_writer_aif
                test_generic_writer_aif.cpp
				        )
target_link_libraries(test_generic_writer_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test( Generic_writer_aif_obj_simple
          test_generic_writer_aif  ${TESTING_DATA_DIR}/CubeTriangleFaces.obj  test_generic_writer_obj_simple.obj  ${TESTING_DATA_DIR}/CubeTriangleFaces.writer-output.obj.ref)
add_test( Generic_writer_aif_obj_with_vertex-normals
          test_generic_writer_aif  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.obj  test_generic_writer_obj_with_vertex-normals.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.writer-output.obj.ref)
add_test( Generic_writer_aif_obj_with_vertex-normals_and_vertex-color
          test_generic_writer_aif  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.obj  test_generic_writer_obj_with_vertex-normals_and_vertex-color.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.writer-output.obj.ref)
add_test( Generic_writer_aif_obj_with_vertex-normals_and_vertex-texture
          test_generic_writer_aif  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.obj  test_generic_writer_obj_with_vertex-normals_and_vertex-texture.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.writer-output.obj.ref)
add_test( Generic_writer_aif_obj_with_halfedge-texture
          test_generic_writer_aif  ${TESTING_DATA_DIR}/textures/usb/usb.obj  test_generic_writer_obj_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/usb/usb.writer-output.obj.ref)
add_test( Generic_writer_aif_off_flat
          test_generic_writer_aif  ${TESTING_DATA_DIR}/flat.off  test_generic_writer_off_flat.off  ${TESTING_DATA_DIR}/flat.off)
add_test( Generic_writer_aif_off_aircraft
          test_generic_writer_aif  ${TESTING_DATA_DIR}/aircraft.off  test_generic_writer_off_aircraft.off  ${TESTING_DATA_DIR}/aircraft.off)
add_test( Generic_writer_aif_obj_non-manifold_with_vertex-texture
          test_generic_writer_aif  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.obj  test_generic_writer_obj_non-manifold_with_vertex-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.writer-output.obj.ref)
add_test( Generic_writer_aif_obj_non-manifold_with_halfedge-texture
          test_generic_writer_aif  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.obj  test_generic_writer_obj_non-manifold_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.writer-output.obj.ref)
add_test( Generic_writer_aif_off_with_face-color
          test_generic_writer_aif  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off  test_generic_writer_off_with_face-color.off  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off)
add_test( Generic_writer_aif_obj_with_face-normal
          test_generic_writer_aif  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.obj  test_generic_writer_obj_with_face-normal.obj  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.writer-output.obj.ref)
# For the 2 next tests, the comparison beteween files fails, but it should not
#if (BUILD_USE_VTK )
#  add_test( Generic_writer_aif_vtk_cube
#          test_generic_writer_aif  ${TESTING_DATA_DIR}/cube.vtk test_generic_writer_vtk_cube.vtk ${TESTING_DATA_DIR}/cube.vtk)
#endif()
#add_test( Generic_writer_aif_msh_SimpleBeam
#          test_generic_writer_aif  ${TESTING_DATA_DIR}/curves/SimpleBeam0000016.msh  test_generic_writer_msh_SimpleBeam.msh  ${TESTING_DATA_DIR}/curves/SimpleBeam0000016.msh)
#####
add_executable( test_mesh_defects_aif
                test_mesh_defects_aif.cpp
				)
target_link_libraries(test_mesh_defects_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )
add_test( Test_mesh_defects_aif
          test_mesh_defects_aif  ${TESTING_DATA_DIR}/TriangleMeshCutVertices.off)
#####
add_executable( test_load_dangling_edge_aif test_load_dangling_edge_aif.cpp )
target_compile_definitions( test_load_dangling_edge_aif PUBLIC
                            "-DUSE_AIF")
target_compile_definitions( test_load_dangling_edge_aif PUBLIC
                            #"-DWITH_GENERIC_READER")
                            "-DWITH_NATIVE_READER")
target_link_libraries(test_load_dangling_edge_aif
                      ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                      )

if ( BUILD_USE_VTK )
add_test( Test_load_dangling_edge_aif_vtk
          test_load_dangling_edge_aif  ${TESTING_DATA_DIR}/dangling_edge.vtk)	# pass with native and generic readers
endif()
add_test( Test_load_dangling_edge_aif_off_1
          test_load_dangling_edge_aif  ${TESTING_DATA_DIR}/dangling_edge.off) # pass with native reader, but do not pass with generic reader
add_test( Test_load_dangling_edge_aif_off_2
          test_load_dangling_edge_aif  ${TESTING_DATA_DIR}/dangling_edge2.off) # pass with native reader, but do not pass with generic reader
add_test( Test_load_dangling_edge_aif_off_3
          test_load_dangling_edge_aif  ${TESTING_DATA_DIR}/dangling_edge3.off) # pass with native reader, but do not pass with generic reader

#####
#TODO-elo Test compression valence compilation with AIF
#TODO-elo disabled because of compilation error
#TODO-elo-off  add_executable( test_compression_valence_aif
#TODO-elo-off                  test_compression_valence_aif.cpp
#TODO-elo-off  				)
#TODO-elo-off  target_link_libraries(test_compression_valence_aif
#TODO-elo-off                        ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
#TODO-elo-off  					             )


##TODO-elo make a function
## variables to be exported to external projects

# include dirs
get_property( INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY INCLUDE_DIRECTORIES
              )
set( MEPP2_AIF_INCLUDE_DIRECTORIES
     "${INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# system include dirs
get_property( SYSTEM_INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY SYSTEM_INCLUDE_DIRECTORIES
              )
set( MEPP2_AIF_SYSTEM_INCLUDE_DIRECTORIES
     "${SYSTEM_INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# libraries to link with
set( MEPP2_AIF_LIBRARIES
     ${Boost_LIBRARIES}
     ${FBX_LIBRARY}
     ${VTK_LIBRARIES}
     CACHE INTERNAL "no-desc" FORCE
     )

# link libraries (to link all target with, silently!)
get_property( LINK_LIBRARIES_TMP
              DIRECTORY
              PROPERTY LINK_LIBRARIES
              )
set( MEPP2_AIF_LINK_LIBRARIES
     "${LINK_LIBRARIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# link dirs
get_property( LINK_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY LINK_DIRECTORIES
              )
set( MEPP2_AIF_LINK_DIRECTORIES
     "${LINK_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# compilation flags and define
get_property( COMPILE_DEFINITIONS_TMP
              DIRECTORY
              PROPERTY COMPILE_DEFINITIONS
              )
set( MEPP2_AIF_COMPILE_DEFINITIONS
     ${COMPILE_DEFINITIONS_TMP}
     CACHE INTERNAL "no-desc" FORCE
     )
