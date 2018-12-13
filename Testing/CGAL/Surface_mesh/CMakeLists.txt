####
add_executable( test_msdm2
                test_msdm2.cpp
                )
target_link_libraries( test_msdm2
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test(Test_MSDM2_Cow
         test_msdm2
         ${TESTING_DATA_DIR}/cow.off
         ${TESTING_DATA_DIR}/cow_simp.off
         surfacemesh/cow_msdm2_output.obj
         ${TESTING_DATA_DIR}/cow_msdm2_reference.obj )

####
add_executable( test_jnd
                test_jnd.cpp
                )
target_link_libraries( test_jnd
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test(Test_Just_Noticeable_Distortion_Cow
         test_jnd
         ${TESTING_DATA_DIR}/cow_simp.off
         surfacemesh/cow_jnd_output.obj
         ${TESTING_DATA_DIR}/cow_jnd_reference.obj )

#####
add_executable( test_surface_mesh test_surface_mesh.cpp )
target_link_libraries( test_surface_mesh ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
add_test(Test_surface_mesh
         test_surface_mesh
         ${TESTING_DATA_DIR}/flat.off
         )

#####
add_executable( test_geometry_traits_unit_normal_surfacemesh
                test_geometry_traits_unit_normal_surfacemesh.cpp )
target_link_libraries( test_geometry_traits_unit_normal_surfacemesh
                        ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
add_test( Geometry_traits_unit_normal_surfacemesh
          test_geometry_traits_unit_normal_surfacemesh )

#####
add_executable( test_euler_add_center_vertex_surfacemesh
                test_euler_add_center_vertex_surfacemesh.cpp )
target_link_libraries( test_euler_add_center_vertex_surfacemesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
add_test(Euler_add_center_vertex_surfacemesh
         test_euler_add_center_vertex_surfacemesh
         ${TESTING_DATA_DIR}/flat.off
         ${TESTING_DATA_DIR}/euler_add_center_vertex.off )

#####
add_executable( test_calculate_scaling_surfacemesh
                test_calculate_scaling_surfacemesh.cpp )
target_link_libraries( test_calculate_scaling_surfacemesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )
add_test( Calculate_scaling_surfacemesh
          test_calculate_scaling_surfacemesh
         ${TESTING_DATA_DIR}/flat.off )

#####
### Checking of concepts is done at build time (no add_test required)
add_executable( test_geometry_traits_concept_surfacemesh
                test_geometry_traits_concept_surfacemesh.cpp )
target_link_libraries( test_geometry_traits_concept_surfacemesh
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

add_executable( test_complying_concepts_surfacemesh
  test_complying_concepts_surfacemesh.cpp )
target_link_libraries( test_complying_concepts_surfacemesh
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

#####
add_executable( test_smoothing_surfacemesh
                test_smoothing_surfacemesh.cpp )
target_link_libraries( test_smoothing_surfacemesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

add_test( Smoothing_surfacemesh
          test_smoothing_surfacemesh
         ${TESTING_DATA_DIR}/cow_simp.off )
add_test( Smoothing_surfacemesh_border
          test_smoothing_surfacemesh
         ${TESTING_DATA_DIR}/flat.off )

#####
add_executable( test_curvature_filter_surfacemesh
                test_curvature_filter_surfacemesh.cpp
                )
target_link_libraries(test_curvature_filter_surfacemesh
                      ${CGAL_LIBRARY}
                      ${CGAL_3RD_PARTY_LIBRARIES}
                      ${Boost_LIBRARIES}
                      ${FBX_LIBRARY}
                      ${VTK_LIBRARIES}
                      )

#####
add_executable( test_generic_writer_surfacemesh
                test_generic_writer_surfacemesh.cpp
				)
target_link_libraries( test_generic_writer_surfacemesh
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/surfacemesh)
add_test( Generic_writer_surfacemesh_obj_simple
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/CubeTriangleFaces.obj  surfacemesh/test_generic_writer_obj_simple.obj  ${TESTING_DATA_DIR}/CubeTriangleFaces.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_obj_with_vertex-normals
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.obj  surfacemesh/test_generic_writer_obj_with_vertex-normals.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormals.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_obj_with_vertex-normals_and_vertex-color
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.obj  surfacemesh/test_generic_writer_obj_with_vertex-normals_and_vertex-color.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndColor.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_obj_with_vertex-normals_and_vertex-texture
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.obj  surfacemesh/test_generic_writer_obj_with_vertex-normals_and_vertex-texture.obj  ${TESTING_DATA_DIR}/CubeTriangleFacesWithVertexNormalsAndTexture.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_obj_with_halfedge-texture
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/textures/usb/usb.obj  surfacemesh/test_generic_writer_obj_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/usb/usb.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_off_flat
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/flat.off  surfacemesh/test_generic_writer_off_flat.off  ${TESTING_DATA_DIR}/flat.off)
add_test( Generic_writer_surfacemesh_off_aircraft
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/aircraft.off  surfacemesh/test_generic_writer_off_aircraft.off  ${TESTING_DATA_DIR}/aircraft.off)
add_test( Generic_writer_surfacemesh_obj_non-manifold_with_vertex-texture
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.obj  surfacemesh/test_generic_writer_obj_non-manifold_with_vertex-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-vertextexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_surfacemesh_obj_non-manifold_with_halfedge-texture
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.obj  surfacemesh/test_generic_writer_obj_non-manifold_with_halfedge-texture.obj  ${TESTING_DATA_DIR}/textures/simple/simple-non-manifold-halfedgetexcoord.writer-output.duplicatedvertex.obj.ref)
add_test( Generic_writer_surfacemesh_off_with_face-color
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off  surfacemesh/test_generic_writer_off_with_face-color.off  ${TESTING_DATA_DIR}/CubeNonTriangleFacesWithFaceColor.off)
add_test( Generic_writer_surfacemesh_obj_with_face-normal
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.obj  surfacemesh/test_generic_writer_obj_with_face-normal.obj  ${TESTING_DATA_DIR}/textures/simple/simple_with_face_normals.writer-output.obj.ref)
add_test( Generic_writer_surfacemesh_obj_with_multitexture
          test_generic_writer_surfacemesh  ${TESTING_DATA_DIR}/textures/simple-multitexture/simple-multitexture.obj  surfacemesh/test_generic_writer_obj_with_multitexture.obj  ${TESTING_DATA_DIR}/textures/simple-multitexture/simple-multitexture.obj.ref)

#####
add_executable( test_compression_valence_surfacemesh
                test_compression_valence_surfacemesh.cpp
				)
target_link_libraries( test_compression_valence_surfacemesh
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_COMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_COMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

# test rgb_nefertiti.off/Surfacemesh/v100/qtz10/noQA
add_test( Test_compression_valence_surfacemesh_nefertiti_v100_qtz10_noQA
		  test_compression_valence_surfacemesh
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
			test_compression_valence_surfacemesh_nefertiti_v100_qtz10_noQA.off
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.off
			-withCompr
			-maxV 100
			-Qbits 10
			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v100/qtz12/noQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v100_qtz12_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v100_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz12/noQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v100/qtz10/noQA
	add_test( Test_compression_valence_surfacemesh_monkey_v100_qtz10_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_surfacemesh_monkey_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v500/qtz15/noQA
	add_test( Test_compression_valence_surfacemesh_monkey_v500_qtz15_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_surfacemesh_monkey_v500_qtz15_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.off
				-withCompr
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v100/qtz10/noQA
	add_test( Test_compression_valence_surfacemesh_swirl_v100_qtz10_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_surfacemesh_swirl_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_surfacemesh_swirl_v1000_qtz10_noQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_surfacemesh_swirl_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

	# test rgb_nefertiti.off/Surfacemesh/v100/qtz10/withQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v100_qtz10_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v100/qtz12/withQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v100_qtz12_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v100_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz12/withQA
	add_test( Test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v100/qtz10/withQA
	add_test( Test_compression_valence_surfacemesh_monkey_v100_qtz10_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_surfacemesh_monkey_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v500/qtz15/withQA
	add_test( Test_compression_valence_surfacemesh_monkey_v500_qtz15_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_surfacemesh_monkey_v500_qtz15_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.off
				-withCompr
				-withQA
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v100/qtz10/withQA
	# the mesh resulting from the decompression of the p3d file
	# has a neat geometry issue on a border ; this issue doesn't
	# exist with Mepp1
	#add_test( Test_compression_valence_surfacemesh_swirl_v100_qtz10_withQA
	#          test_compression_valence_surfacemesh
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
	#            test_compression_valence_surfacemesh_swirl_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10 )
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_surfacemesh_swirl_v1000_qtz10_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_surfacemesh_swirl_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test swirl_nocolor.off/Surfacemesh/v100/qtz10/withQA
	# uncompressed final mesh conform to Mepp1
	add_test( Test_compression_valence_surfacemesh_swirl_nocolor_v100_qtz10_withQA
			  test_compression_valence_surfacemesh
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor.off
				test_compression_valence_surfacemesh_swirl_nocolor_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor_mepp2_polyhedron_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl_colorgradient.off/Surfacemesh/v100/qtz10/withQA
	# uncompressed final mesh geometrically identical to Mepp1 (but there
	# are minor differences in colors)
	#add_test( Test_compression_valence_surfacemesh_swirl_colorgradient_v100_qtz10_withQA
	#          test_compression_valence_surfacemesh
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient.off
	#            test_compression_valence_surfacemesh_swirl_colorgradient_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient_mepp2_polyhedron_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )

#####
add_executable( test_decompression_valence_surfacemesh
                test_decompression_valence_surfacemesh.cpp
				)
target_link_libraries( test_decompression_valence_surfacemesh
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/noQA from mepp1
add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_noQA
		  test_decompression_valence_surfacemesh
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d
		  test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_noQA.off
		  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
		  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
		  )

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/noQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_noQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/noQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_surfacemesh_monkey_v100_qtz10_noQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.p3d
			  test_decompression_valence_surfacemesh_monkey_v100_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/noQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_monkey_v500_qtz15_noQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d
			  test_decompression_valence_surfacemesh_monkey_v500_qtz15_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )

	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/withQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_surfacemesh_monkey_v100_qtz10_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_monkey_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/withQA from mepp1
	add_test( Test_decompression_valence_surfacemesh_monkey_v500_qtz15_withQA
			  test_decompression_valence_surfacemesh
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d
			  test_decompression_valence_surfacemesh_monkey_v500_qtz15_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )

#####
#
# testing compression+decompression on the same platform
#
if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/noQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v100_qtz10_noQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_noQA_sameptf
				  APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v100_qtz10_noQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/noQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_noQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/noQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v100_qtz12_noQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_noQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v100_qtz12_noQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/noQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_noQA
				  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/noQA compressed with mepp2 on the same platform
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_surfacemesh_monkey_v100_qtz10_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_monkey_v100_qtz10_noQA.p3d
			  test_decompression_valence_surfacemesh_monkey_v100_qtz10_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_monkey_v100_qtz10_noQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_monkey_v100_qtz10_noQA
				  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/noQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_monkey_v500_qtz15_noQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_monkey_v500_qtz15_noQA.p3d
			  test_decompression_valence_surfacemesh_monkey_v500_qtz15_noQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_monkey_v500_qtz15_noQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_monkey_v500_qtz15_noQA
				  )

	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/withQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v100_qtz10_withQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v100_qtz10_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v100_qtz10_withQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/withQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v1000_qtz10_withQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/withQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v100_qtz12_withQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v100_qtz12_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v100_qtz12_withQA
				  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/withQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA.p3d
			  test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_nefertiti_v1000_qtz12_withQA
				  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/withQA compressed with mepp2 on the same platform
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_surfacemesh_monkey_v100_qtz10_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_monkey_v100_qtz10_withQA.p3d
			  test_decompression_valence_surfacemesh_monkey_v100_qtz10_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_monkey_v100_qtz10_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_monkey_v100_qtz10_withQA
				  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/withQA compressed with mepp2 on the same platform
	add_test( Test_decompression_valence_surfacemesh_monkey_v500_qtz15_withQA_sameptf
			  test_decompression_valence_surfacemesh
			  test_compression_valence_surfacemesh_monkey_v500_qtz15_withQA.p3d
			  test_decompression_valence_surfacemesh_monkey_v500_qtz15_withQA_sameptf.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
    set_property( TEST  Test_decompression_valence_surfacemesh_monkey_v500_qtz15_withQA_sameptf
	              APPEND PROPERTY DEPENDS  Test_compression_valence_surfacemesh_monkey_v500_qtz15_withQA
				  )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
