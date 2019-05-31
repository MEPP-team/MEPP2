if( NOT BUILD_USE_OPENMESH )
  return()
endif()

include_directories( ${OPENMESH_INCLUDE_DIR} )

if( MSVC )
  # AUTOLINK is used for CGAL on Windows. Disable it for OpenMesh-only tests
  # in order to prevent false link error.
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

#####
add_executable( test_openmesh test_openmesh.cpp )
target_link_libraries( test_openmesh ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Test_openmesh
         test_openmesh
         ${TESTING_DATA_DIR}/flat.off
         )

#####
add_executable( test_collapse_edge_openmesh
                test_collapse_edge_openmesh.cpp )
target_link_libraries( test_collapse_edge_openmesh ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Collapse_edge_openmesh_1
         test_collapse_edge_openmesh
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off )
add_test(Collapse_edge_openmesh_2
         test_collapse_edge_openmesh
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )

#####
add_executable( test_collapse_edge_Euler_openmesh
                test_collapse_edge_euler_openmesh.cpp )
target_link_libraries( test_collapse_edge_Euler_openmesh
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Collapse_edge_Euler_openmesh_1
         test_collapse_edge_Euler_openmesh
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off )
add_test(Collapse_edge_Euler_openmesh_2
         test_collapse_edge_Euler_openmesh
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )

#####
add_executable( test_calculate_scaling_openmesh
                test_calculate_scaling_openmesh.cpp )
target_link_libraries( test_calculate_scaling_openmesh
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )
add_test( Calculate_scaling_openmesh
         test_calculate_scaling_openmesh
         ${TESTING_DATA_DIR}/flat.off )
####
add_executable( test_geometry_traits_unit_normal_openmesh
                test_geometry_traits_unit_normal_openmesh.cpp )
target_link_libraries( test_geometry_traits_unit_normal_openmesh
                       ${OPENMESH_LIBRARIES} )
add_test( Geometry_traits_unit_normal_openmesh
          test_geometry_traits_unit_normal_openmesh )

#####
add_executable( test_euler_add_center_vertex_openmesh
                test_euler_add_center_vertex_openmesh.cpp )
target_link_libraries( test_euler_add_center_vertex_openmesh
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Euler_add_center_vertex_openmesh
         test_euler_add_center_vertex_openmesh
         ${TESTING_DATA_DIR}/flat.off
         ${TESTING_DATA_DIR}/euler_add_center_vertex.off )

#####
# Checking of concepts is done at build time (no add_test required)
add_executable( test_geometry_traits_concept_openmesh
                test_geometry_traits_concept_openmesh.cpp )
target_link_libraries( test_geometry_traits_concept_openmesh
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )

add_executable( test_complying_concepts_openmesh
                test_complying_concepts_openmesh.cpp )
target_link_libraries( test_complying_concepts_openmesh
                       ${OPENMESH_LIBRARIES} ${Boost_LIBRARIES} )

#####
add_executable( test_smoothing_openmesh
                test_smoothing_openmesh.cpp )
target_link_libraries( test_smoothing_openmesh
                       ${OPENMESH_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test( Smoothing_openmesh
         test_smoothing_openmesh
         ${TESTING_DATA_DIR}/cow_simp.off )
add_test( Smoothing_openmesh_border
         test_smoothing_openmesh
         ${TESTING_DATA_DIR}/flat.off )

#####
add_executable( test_curvature_filter_openmesh
                test_curvature_filter_openmesh.cpp
                )
target_link_libraries(  test_curvature_filter_openmesh
                        ${OPENMESH_LIBRARIES}
                        ${Boost_LIBRARIES}
                        )

#####
add_executable( test_generic_writer_openmesh
                test_generic_writer_openmesh.cpp
				)
target_link_libraries( test_generic_writer_openmesh
                       ${OPENMESH_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test( Generic_writer_openmesh_obj_simple
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/CubeTriangleFaces.obj  test_generic_writer_obj_simple.obj  ${TESTING_DATA_DIR}/CubeTriangleFaces.writer-output.obj.ref)
add_test( Generic_writer_openmesh_obj_with_vertex-normals
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.obj  test_generic_writer_obj_with_vertex-normals.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.writer-output.obj.ref)
add_test( Generic_writer_openmesh_obj_with_vertex-normals_and_vertex-color
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.obj  test_generic_writer_obj_with_vertex-normals_and_vertex-color.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.writer-output.obj.ref)
add_test( Generic_writer_openmesh_obj_with_vertex-normals_and_vertex-texture
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.obj  test_generic_writer_obj_with_vertex-normals_and_vertex-texture.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.writer-output.obj.ref)
add_test( Generic_writer_openmesh_obj_with_halfedge-texture
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/textures/usb/usb.obj  test_generic_writer_obj_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/usb/usb.writer-output.obj.ref)
add_test( Generic_writer_openmesh_off_flat
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/flat.off  test_generic_writer_off_flat.off  ${TESTING_DATA_DIR}/flat.off)
add_test( Generic_writer_openmesh_off_aircraft
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/aircraft.off  test_generic_writer_off_aircraft.off  ${TESTING_DATA_DIR}/aircraft.off)
add_test( Generic_writer_openmesh_obj_non-manifold_with_vertex-texture
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.obj  test_generic_writer_obj_non-manifold_with_vertex-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_openmesh_obj_non-manifold_with_halfedge-texture
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.obj  test_generic_writer_obj_non-manifold_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_openmesh_off_with_face-color
          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off  test_generic_writer_off_with_face-color.off  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off)
#TODO-elo  diagnose, fix, then restore test
# the following test fails with OpenMesh only (ok with AIF and CGALs)
# because normal in output file is "vn -3.83376e-14 1.53036e-41 -1"
# instead of "vn 0 0 -1"
# this may also explain why the rendering of this mesh in mepp-gui
# is different with OpenMesh
#add_test( Generic_writer_openmesh_obj_with_face-normal
#          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.obj  test_generic_writer_obj_with_face-normal.obj  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.writer-output.obj.ref)
# For the next test, the comparison beteween files fails, but it should not
#add_test( Generic_writer_openmesh_msh_SimpleBeam
#          test_generic_writer_openmesh  ${TESTING_DATA_DIR}/curves/SimpleBeam0000016.msh  test_generic_writer_msh_SimpleBeam.msh  ${TESTING_DATA_DIR}/curves/SimpleBeam0000016.msh)

#####
add_executable( test_compression_valence_openmesh
                test_compression_valence_openmesh.cpp
				)
target_link_libraries( test_compression_valence_openmesh
                       ${OPENMESH_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE 2e-5)
#set(TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE 0)  # for exact comparison

# test rgb_nefertiti.off/Openmesh/v100/qtz10/noQA
add_test( Test_compression_valence_openmesh_nefertiti_v100_qtz10_noQA
		  test_compression_valence_openmesh
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
			test_compression_valence_openmesh_nefertiti_v100_qtz10_noQA.off
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.off
			-withCompr
			-maxV 100
			-Qbits 10
			-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test rgb_nefertiti.off/Openmesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_openmesh_nefertiti_v1000_qtz10_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_nefertiti.off/Openmesh/v100/qtz12/noQA
	add_test( Test_compression_valence_openmesh_nefertiti_v100_qtz12_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v100_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_nefertiti.off/Openmesh/v1000/qtz12/noQA
	add_test( Test_compression_valence_openmesh_nefertiti_v1000_qtz12_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v1000_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_monkey.off/Openmesh/v100/qtz10/noQA
	add_test( Test_compression_valence_openmesh_monkey_v100_qtz10_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_openmesh_monkey_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_monkey.off/Openmesh/v500/qtz15/noQA
	add_test( Test_compression_valence_openmesh_monkey_v500_qtz15_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_openmesh_monkey_v500_qtz15_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.off
				-withCompr
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_swirl.off/Openmesh/v100/qtz10/noQA
	add_test( Test_compression_valence_openmesh_swirl_v100_qtz10_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_openmesh_swirl_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_swirl.off/Openmesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_openmesh_swirl_v1000_qtz10_noQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_openmesh_swirl_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})

	# test rgb_nefertiti.off/Openmesh/v100/qtz10/withQA
	add_test( Test_compression_valence_openmesh_nefertiti_v100_qtz10_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_nefertiti.off/Openmesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_openmesh_nefertiti_v1000_qtz10_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_nefertiti.off/Openmesh/v100/qtz12/withQA
	add_test( Test_compression_valence_openmesh_nefertiti_v100_qtz12_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v100_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_nefertiti.off/Openmesh/v1000/qtz12/withQA
	add_test( Test_compression_valence_openmesh_nefertiti_v1000_qtz12_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_openmesh_nefertiti_v1000_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_monkey.off/Openmesh/v100/qtz10/withQA
	add_test( Test_compression_valence_openmesh_monkey_v100_qtz10_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_openmesh_monkey_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_monkey.off/Openmesh/v500/qtz15/withQA
	add_test( Test_compression_valence_openmesh_monkey_v500_qtz15_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_openmesh_monkey_v500_qtz15_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.off
				-withCompr
				-withQA
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_swirl.off/Openmesh/v100/qtz10/withQA
	# the mesh resulting from the decompression of the p3d file
	# has a neat geometry issue on a border ; this issue doesn't
	# exist with Mepp1
	#add_test( Test_compression_valence_openmesh_swirl_v100_qtz10_withQA
	#          test_compression_valence_openmesh
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
	#            test_compression_valence_openmesh_swirl_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10 )
	# test rgb_swirl.off/Openmesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_openmesh_swirl_v1000_qtz10_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_openmesh_swirl_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test swirl_nocolor.off/Openmesh/v100/qtz10/withQA
	# uncompressed final mesh conform to Mepp1
	add_test( Test_compression_valence_openmesh_swirl_nocolor_v100_qtz10_withQA
			  test_compression_valence_openmesh
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor.off
				test_compression_valence_openmesh_swirl_nocolor_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor_mepp2_polyhedron_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_OPENMESH_TOLERANCE})
	# test rgb_swirl_colorgradient.off/Openmesh/v100/qtz10/withQA
	# uncompressed final mesh geometrically identical to Mepp1 (but there
	# are minor differences in colors)
	#add_test( Test_compression_valence_openmesh_swirl_colorgradient_v100_qtz10_withQA
	#          test_compression_valence_openmesh
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient.off
	#            test_compression_valence_openmesh_swirl_colorgradient_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient_mepp2_polyhedron_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10 )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )


#####
add_executable( test_decompression_valence_openmesh
                test_decompression_valence_openmesh.cpp
				)
target_link_libraries( test_decompression_valence_openmesh
                       ${OPENMESH_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE 2e-5)
#set(TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE 0)  # for exact comparison

# test decompression of rgb_nefertiti.off/Openmesh/v100/qtz10/noQA from mepp1
add_test( Test_decompression_valence_openmesh_nefertiti_v100_qtz10_noQA
		  test_decompression_valence_openmesh
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d
		  test_decompression_valence_openmesh_nefertiti_v100_qtz10_noQA.off
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
		  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
		  )

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test decompression of rgb_nefertiti.off/Openmesh/v1000/qtz10/noQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v1000_qtz10_noQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v1000_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Openmesh/v100/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v100_qtz12_noQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v100_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Openmesh/v1000/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v1000_qtz12_noQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v1000_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Openmesh/v100/qtz10/noQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_openmesh_monkey_v100_qtz10_noQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.p3d
			  test_decompression_valence_openmesh_monkey_v100_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Openmesh/v500/qtz15/noQA from mepp1
	add_test( Test_decompression_valence_openmesh_monkey_v500_qtz15_noQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d
			  test_decompression_valence_openmesh_monkey_v500_qtz15_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )

	# test decompression of rgb_nefertiti.off/Openmesh/v100/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v100_qtz10_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Openmesh/v1000/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v1000_qtz10_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v1000_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Openmesh/v100/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v100_qtz12_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v100_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Openmesh/v1000/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_openmesh_nefertiti_v1000_qtz12_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d
			  test_decompression_valence_openmesh_nefertiti_v1000_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Openmesh/v100/qtz10/withQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_openmesh_monkey_v100_qtz10_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_openmesh_monkey_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Openmesh/v500/qtz15/withQA from mepp1
	add_test( Test_decompression_valence_openmesh_monkey_v500_qtz15_withQA
			  test_decompression_valence_openmesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d
			  test_decompression_valence_openmesh_monkey_v500_qtz15_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_OPENMESH_TOLERANCE}
			  )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )


##TODO-elo make a function
## variables to be exported to external projects

# include dirs
get_property( INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY INCLUDE_DIRECTORIES
              )
set( MEPP2_OPENMESH_INCLUDE_DIRECTORIES
     "${INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# system include dirs
get_property( SYSTEM_INCLUDE_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY SYSTEM_INCLUDE_DIRECTORIES
              )
set( MEPP2_OPENMESH_SYSTEM_INCLUDE_DIRECTORIES
     "${SYSTEM_INCLUDE_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# libraries to link with
set( MEPP2_OPENMESH_LIBRARIES
     ${OPENMESH_LIBRARIES}
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
set( MEPP2_OPENMESH_LINK_LIBRARIES
     "${LINK_LIBRARIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# link dirs
get_property( LINK_DIRECTORIES_TMP
              DIRECTORY
              PROPERTY LINK_DIRECTORIES
              )
set( MEPP2_OPENMESH_LINK_DIRECTORIES
     "${LINK_DIRECTORIES_TMP}"
     CACHE INTERNAL "no-desc" FORCE
     )

# compilation flags and define
get_property( COMPILE_DEFINITIONS_TMP
              DIRECTORY
              PROPERTY COMPILE_DEFINITIONS
              )
set( MEPP2_OPENMESH_COMPILE_DEFINITIONS
     ${COMPILE_DEFINITIONS_TMP}
     CACHE INTERNAL "no-desc" FORCE
     )
