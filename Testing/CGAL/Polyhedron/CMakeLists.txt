####
add_executable( test_flip_edge_polyhedron
                test_flip_edge_polyhedron.cpp )
target_link_libraries( test_flip_edge_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test(Flip_edge_polyhedron_1
         test_flip_edge_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_flip_0.off
         ${TESTING_DATA_DIR}/square3.off )
add_test(Flip_edge_polyhedron_2
         test_flip_edge_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_flip_1.off
         ${TESTING_DATA_DIR}/square3.off )

####
add_executable( test_split_edge_polyhedron
                test_split_edge_polyhedron.cpp )
target_link_libraries( test_split_edge_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test(Split_edge_polyhedron_1
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         0
         flat_0.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_polyhedron_2
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         1
         flat_1.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_polyhedron_3
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         2
         flat_2.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_polyhedron_4
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         3
         flat_3.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_polyhedron_5
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_split_0.off
         ${TESTING_DATA_DIR}/square4.off )
add_test(Split_edge_polyhedron_6
         test_split_edge_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_split_1.off
         ${TESTING_DATA_DIR}/square4.off )

####
add_executable( test_split_edge_Euler_polyhedron
                test_split_edge_euler_polyhedron.cpp )
target_link_libraries( test_split_edge_Euler_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test(Split_edge_Euler_polyhedron_1
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         0
         flat_Euler_0.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_Euler_polyhedron_2
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         1
         flat_Euler_1.off
         ${TESTING_DATA_DIR}/split_edge_border_12_or_21.off )
add_test(Split_edge_Euler_polyhedron_3
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         2
         flat_Euler_2.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_Euler_polyhedron_4
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         3
         flat_Euler_3.off
         ${TESTING_DATA_DIR}/split_edge_34_or_43.off )
add_test(Split_edge_Euler_polyhedron_5
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         0
         square2_Euler_0.off
         ${TESTING_DATA_DIR}/square4.off )
add_test(Split_edge_Euler_polyhedron_6
         test_split_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/square2.off
         1
         square2_Euler_1.off
         ${TESTING_DATA_DIR}/square4.off )

#####
add_executable( test_collapse_edge_polyhedron
                test_collapse_edge_polyhedron.cpp )
target_link_libraries( test_collapse_edge_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

add_test(Collapse_edge_polyhedron_1
         test_collapse_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off )
add_test(Collapse_edge_polyhedron_2
         test_collapse_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )
add_test(Collapse_edge_polyhedron_3
         test_collapse_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         3 )
add_test(Collapse_edge_polyhedron_4
         test_collapse_edge_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         4 )

#####
add_executable( test_collapse_edge_Euler_polyhedron
                test_collapse_edge_polyhedron.cpp )
target_link_libraries( test_collapse_edge_Euler_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

add_test(Collapse_edge_Euler_polyhedron_1
         test_collapse_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         1
         ${TESTING_DATA_DIR}/collapse_edge_34_keep_target.off )
add_test(Collapse_edge_Euler_polyhedron_2
         test_collapse_edge_Euler_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         2
         ${TESTING_DATA_DIR}/collapse_edge_43_keep_target.off )

#####
add_executable( test_polyhedron test_polyhedron.cpp )
target_link_libraries( test_polyhedron ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})

add_test(Test_polyhedron
         test_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         )

#####
add_executable( test_geometry_traits_unit_normal_polyhedron
                test_geometry_traits_unit_normal_polyhedron.cpp )
target_link_libraries( test_geometry_traits_unit_normal_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
add_test( Geometry_traits_unit_normal_polyhedron
          test_geometry_traits_unit_normal_polyhedron )

#####
add_executable( test_euler_add_center_vertex_polyhedron
                test_euler_add_center_vertex_polyhedron.cpp )
target_link_libraries( test_euler_add_center_vertex_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Euler_add_center_vertex_polyhedron
         test_euler_add_center_vertex_polyhedron
         ${TESTING_DATA_DIR}/flat.off
         ${TESTING_DATA_DIR}/euler_add_center_vertex.off )

#####
add_executable( test_calculate_scaling_polyhedron
                test_calculate_scaling_polyhedron.cpp )
target_link_libraries( test_calculate_scaling_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
add_test( Calculate_scaling_polyhedron
          test_calculate_scaling_polyhedron
         ${TESTING_DATA_DIR}/flat.off )

#####
### Checking of concepts is done at build time (no add_test required)
add_executable( test_geometry_traits_concept_polyhedron
                test_geometry_traits_concept_polyhedron.cpp )
target_link_libraries( test_geometry_traits_concept_polyhedron
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

add_executable( test_complying_concepts_polyhedron
  test_complying_concepts_polyhedron.cpp )
target_link_libraries( test_complying_concepts_polyhedron
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

#####
add_executable( test_smoothing_polyhedron
                test_smoothing_polyhedron.cpp )
target_link_libraries( test_smoothing_polyhedron
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test( Smoothing_polyhedron
          test_smoothing_polyhedron
         ${TESTING_DATA_DIR}/cow_simp.off )
add_test( Smoothing_polyhedron_border
          test_smoothing_polyhedron
         ${TESTING_DATA_DIR}/flat.off )

#####
add_executable( test_curvature_filter_polyhedron
                test_curvature_filter_polyhedron.cpp
                )
target_link_libraries(test_curvature_filter_polyhedron
                      ${CGAL_LIBRARY}
                      ${CGAL_3RD_PARTY_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${FBX_LIBRARY}
                      ${VTK_LIBRARIES}
                      )

#####
add_executable( test_generic_writer_polyhedron
                test_generic_writer_polyhedron.cpp
				)
target_link_libraries( test_generic_writer_polyhedron
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/polyhedron)
add_test( Generic_writer_polyhedron_obj_simple
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/CubeTriangleFaces.obj  polyhedron/test_generic_writer_obj_simple.obj  ${TESTING_DATA_DIR}/CubeTriangleFaces.writer-output.obj.ref)
add_test( Generic_writer_polyhedron_obj_with_vertex-normals
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.obj  polyhedron/test_generic_writer_obj_with_vertex-normals.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.writer-output.obj.ref)
add_test( Generic_writer_polyhedron_obj_with_vertex-normals_and_vertex-color
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.obj  polyhedron/test_generic_writer_obj_with_vertex-normals_and_vertex-color.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.writer-output.obj.ref)
add_test( Generic_writer_polyhedron_obj_with_vertex-normals_and_vertex-texture
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.obj  polyhedron/test_generic_writer_obj_with_vertex-normals_and_vertex-texture.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.writer-output.obj.ref)
add_test( Generic_writer_polyhedron_obj_with_halfedge-texture
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/textures/usb/usb.obj  polyhedron/test_generic_writer_obj_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/usb/usb.writer-output.obj.ref)
add_test( Generic_writer_polyhedron_off_flat
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/flat.off  polyhedron/test_generic_writer_off_flat.off  ${TESTING_DATA_DIR}/flat.off)
add_test( Generic_writer_polyhedron_off_aircraft
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/aircraft.off  polyhedron/test_generic_writer_off_aircraft.off  ${TESTING_DATA_DIR}/aircraft.off)
add_test( Generic_writer_polyhedron_obj_non-manifold_with_vertex-texture
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.obj  polyhedron/test_generic_writer_obj_non-manifold_with_vertex-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_polyhedron_obj_non-manifold_with_halfedge-texture
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.obj  polyhedron/test_generic_writer_obj_non-manifold_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_polyhedron_off_with_face-color
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off  polyhedron/test_generic_writer_off_with_face-color.off  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off)
add_test( Generic_writer_polyhedron_obj_with_face-normal
          test_generic_writer_polyhedron  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.obj  polyhedron/test_generic_writer_obj_with_face-normal.obj  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.writer-output.obj.ref)

#####
add_executable( test_compression_valence_polyhedron
                test_compression_valence_polyhedron.cpp
				)
target_link_libraries( test_compression_valence_polyhedron
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_COMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_COMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

# test rgb_nefertiti.off/Polyhedron/v100/qtz10/noQA
add_test( Test_compression_valence_polyhedron_nefertiti_v100_qtz10_noQA
		  test_compression_valence_polyhedron
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
			test_compression_valence_polyhedron_nefertiti_v100_qtz10_noQA.off
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.off
			-withCompr
			-maxV 100
			-Qbits 10
			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test rgb_nefertiti.off/Polyhedron/v1000/qtz10/noQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v1000_qtz10_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Polyhedron/v100/qtz12/noQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v100_qtz12_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v100_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Polyhedron/v1000/qtz12/noQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v1000_qtz12_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v1000_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Polyhedron/v100/qtz10/noQA
	add_test( Test_compression_valence_polyhedron_monkey_v100_qtz10_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_polyhedron_monkey_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Polyhedron/v500/qtz15/noQA
	add_test( Test_compression_valence_polyhedron_monkey_v500_qtz15_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_polyhedron_monkey_v500_qtz15_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.off
				-withCompr
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Polyhedron/v100/qtz10/noQA
	add_test( Test_compression_valence_polyhedron_swirl_v100_qtz10_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_polyhedron_swirl_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Polyhedron/v1000/qtz10/noQA
	add_test( Test_compression_valence_polyhedron_swirl_v1000_qtz10_noQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_polyhedron_swirl_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

	# test rgb_nefertiti.off/Polyhedron/v100/qtz10/withQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v100_qtz10_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Polyhedron/v1000/qtz10/withQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v1000_qtz10_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Polyhedron/v100/qtz12/withQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v100_qtz12_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v100_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Polyhedron/v1000/qtz12/withQA
	add_test( Test_compression_valence_polyhedron_nefertiti_v1000_qtz12_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_polyhedron_nefertiti_v1000_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Polyhedron/v100/qtz10/withQA
	add_test( Test_compression_valence_polyhedron_monkey_v100_qtz10_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_polyhedron_monkey_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Polyhedron/v500/qtz15/withQA
	add_test( Test_compression_valence_polyhedron_monkey_v500_qtz15_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_polyhedron_monkey_v500_qtz15_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.off
				-withCompr
				-withQA
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Polyhedron/v100/qtz10/withQA
	# the mesh resulting from the decompression of the p3d file
	# has a neat geometry issue on a border ; this issue doesn't
	# exist with Mepp1
	#add_test( Test_compression_valence_polyhedron_swirl_v100_qtz10_withQA
	#          test_compression_valence_polyhedron
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
	#            test_compression_valence_polyhedron_swirl_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Polyhedron/v1000/qtz10/withQA
	add_test( Test_compression_valence_polyhedron_swirl_v1000_qtz10_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_polyhedron_swirl_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test swirl_nocolor.off/Polyhedron/v100/qtz10/withQA
	# uncompressed final mesh conform to Mepp1
	add_test( Test_compression_valence_polyhedron_swirl_nocolor_v100_qtz10_withQA
			  test_compression_valence_polyhedron
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor.off
				test_compression_valence_polyhedron_swirl_nocolor_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor_mepp2_polyhedron_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl_colorgradient.off/Polyhedron/v100/qtz10/withQA
	# uncompressed final mesh geometrically identical to Mepp1 (but there
	# are minor differences in colors)
	#add_test( Test_compression_valence_polyhedron_swirl_colorgradient_v100_qtz10_withQA
	#          test_compression_valence_polyhedron
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient.off
	#            test_compression_valence_polyhedron_swirl_colorgradient_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient_mepp2_polyhedron_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )

#####
add_executable( test_decompression_valence_polyhedron
                test_decompression_valence_polyhedron.cpp
				)
target_link_libraries( test_decompression_valence_polyhedron
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

# test decompression of rgb_nefertiti.off/Polyhedron/v100/qtz10/noQA from mepp1
add_test( Test_decompression_valence_polyhedron_nefertiti_v100_qtz10_noQA
		  test_decompression_valence_polyhedron
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d
		  test_decompression_valence_polyhedron_nefertiti_v100_qtz10_noQA.off
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
		  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
		  )

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test decompression of rgb_nefertiti.off/Polyhedron/v1000/qtz10/noQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v1000_qtz10_noQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v1000_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Polyhedron/v100/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v100_qtz12_noQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v100_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Polyhedron/v1000/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v1000_qtz12_noQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v1000_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Polyhedron/v100/qtz10/noQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_polyhedron_monkey_v100_qtz10_noQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.p3d
			  test_decompression_valence_polyhedron_monkey_v100_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Polyhedron/v500/qtz15/noQA from mepp1
	add_test( Test_decompression_valence_polyhedron_monkey_v500_qtz15_noQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d
			  test_decompression_valence_polyhedron_monkey_v500_qtz15_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )

	# test decompression of rgb_nefertiti.off/Polyhedron/v100/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v100_qtz10_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Polyhedron/v1000/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v1000_qtz10_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v1000_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Polyhedron/v100/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v100_qtz12_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v100_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Polyhedron/v1000/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_polyhedron_nefertiti_v1000_qtz12_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d
			  test_decompression_valence_polyhedron_nefertiti_v1000_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Polyhedron/v100/qtz10/withQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_polyhedron_monkey_v100_qtz10_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_polyhedron_monkey_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Polyhedron/v500/qtz15/withQA from mepp1
	add_test( Test_decompression_valence_polyhedron_monkey_v500_qtz15_withQA
			  test_decompression_valence_polyhedron
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d
			  test_decompression_valence_polyhedron_monkey_v500_qtz15_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )

#####
add_executable( test_boolean_operations_polyhedron
                test_boolean_operations_polyhedron.cpp
				        )
target_link_libraries( test_boolean_operations_polyhedron
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
# test UNION
add_test( Test_boolean_operations_polyhedron_union
          test_boolean_operations_polyhedron
          "union"
          "${TESTING_DATA_DIR}/CubeTriangleFaces.off"
          "${TESTING_DATA_DIR}/tetra.off"
          "${TESTING_DATA_DIR}/boolean_operations/mepp1_cube_tetra_union.off"
          )
# test INTERSECTION
add_test( Test_boolean_operations_polyhedron_inter
          test_boolean_operations_polyhedron
          "inter"
          "${TESTING_DATA_DIR}/CubeTriangleFaces.off"
          "${TESTING_DATA_DIR}/tetra.off"
          "${TESTING_DATA_DIR}/boolean_operations/mepp1_cube_tetra_inter.off"
          )
# test SUBTRACTION
add_test( Test_boolean_operations_polyhedron_minus
          test_boolean_operations_polyhedron
          "minus"
          "${TESTING_DATA_DIR}/CubeTriangleFaces.off"
          "${TESTING_DATA_DIR}/tetra.off"
          "${TESTING_DATA_DIR}/boolean_operations/mepp1_cube_tetra_minus.off"
          )
