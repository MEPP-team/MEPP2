#####
add_executable( test_linear_cell_complex test_linear_cell_complex.cpp )
target_link_libraries( test_linear_cell_complex
   ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
add_test(Test_linear_cell_complex
         test_linear_cell_complex
         ${TESTING_DATA_DIR}/flat.off
         )

#####
add_executable( test_not_2_manifold_linear_cell_complex test_not_2_manifold_linear_cell_complex.cpp )
target_link_libraries( test_not_2_manifold_linear_cell_complex ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES})
add_test(Test_not_2_manifold_linear_cell_complex
         test_not_2_manifold_linear_cell_complex
         )

#####
# Checking of concepts is done at build time (no add_test required)
add_executable( test_geometry_traits_concept_linear_cell_complex
                test_geometry_traits_concept_linear_cell_complex.cpp )
target_link_libraries( test_geometry_traits_concept_linear_cell_complex
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

add_executable( test_complying_concepts_linear_cell_complex
  test_complying_concepts_linear_cell_complex.cpp )
target_link_libraries( test_complying_concepts_linear_cell_complex
  ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES} ${Boost_LIBRARIES} )

#####
add_executable( test_compression_valence_lcc
                test_compression_valence_lcc.cpp
				)
target_link_libraries( test_compression_valence_lcc
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_COMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_COMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

# test rgb_nefertiti.off/Surfacemesh/v100/qtz10/noQA
add_test( Test_compression_valence_lcc_nefertiti_v100_qtz10_noQA
		  test_compression_valence_lcc
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
			test_compression_valence_lcc_nefertiti_v100_qtz10_noQA.off
			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.off
			-withCompr
			-maxV 100
			-Qbits 10
			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_lcc_nefertiti_v1000_qtz10_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v100/qtz12/noQA
	add_test( Test_compression_valence_lcc_nefertiti_v100_qtz12_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v100_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz12/noQA
	add_test( Test_compression_valence_lcc_nefertiti_v1000_qtz12_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v1000_qtz12_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v100/qtz10/noQA
	add_test( Test_compression_valence_lcc_monkey_v100_qtz10_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_lcc_monkey_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v500/qtz15/noQA
	add_test( Test_compression_valence_lcc_monkey_v500_qtz15_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_lcc_monkey_v500_qtz15_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.off
				-withCompr
				-maxV 500
				-Qbits 15
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v100/qtz10/noQA
	add_test( Test_compression_valence_lcc_swirl_v100_qtz10_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_lcc_swirl_v100_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.off
				-withCompr
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v1000/qtz10/noQA
	add_test( Test_compression_valence_lcc_swirl_v1000_qtz10_noQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_lcc_swirl_v1000_qtz10_noQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp2_polyhedron_v1000_qtz10_noQA.ref.off
				-withCompr
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})

	# test rgb_nefertiti.off/Surfacemesh/v100/qtz10/withQA
	add_test( Test_compression_valence_lcc_nefertiti_v100_qtz10_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_lcc_nefertiti_v1000_qtz10_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v100/qtz12/withQA
	add_test( Test_compression_valence_lcc_nefertiti_v100_qtz12_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v100_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_nefertiti.off/Surfacemesh/v1000/qtz12/withQA
	add_test( Test_compression_valence_lcc_nefertiti_v1000_qtz12_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix.off
				test_compression_valence_lcc_nefertiti_v1000_qtz12_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 12
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v100/qtz10/withQA
	add_test( Test_compression_valence_lcc_monkey_v100_qtz10_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_lcc_monkey_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_monkey.off/Surfacemesh/v500/qtz15/withQA
	add_test( Test_compression_valence_lcc_monkey_v500_qtz15_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix.off
				test_compression_valence_lcc_monkey_v500_qtz15_withQA.off
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
	#add_test( Test_compression_valence_lcc_swirl_v100_qtz10_withQA
	#          test_compression_valence_lcc
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
	#            test_compression_valence_lcc_swirl_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10 )
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl.off/Surfacemesh/v1000/qtz10/withQA
	add_test( Test_compression_valence_lcc_swirl_v1000_qtz10_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix.off
				test_compression_valence_lcc_swirl_v1000_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/rgb_swirl_unix_mepp1_v1000_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 1000
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test swirl_nocolor.off/Surfacemesh/v100/qtz10/withQA
	# uncompressed final mesh conform to Mepp1
	add_test( Test_compression_valence_lcc_swirl_nocolor_v100_qtz10_withQA
			  test_compression_valence_lcc
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor.off
				test_compression_valence_lcc_swirl_nocolor_v100_qtz10_withQA.off
				${TESTING_DATA_DIR}/CompressionValenceExample/swirl_nocolor_mepp2_polyhedron_v100_qtz10_withQA.ref.off
				-withCompr
				-withQA
				-maxV 100
				-Qbits 10
				-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
	# test rgb_swirl_colorgradient.off/Surfacemesh/v100/qtz10/withQA
	# uncompressed final mesh geometrically identical to Mepp1 (but there
	# are minor differences in colors)
	#add_test( Test_compression_valence_lcc_swirl_colorgradient_v100_qtz10_withQA
	#          test_compression_valence_lcc
	#		    ${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient.off
	#            test_compression_valence_lcc_swirl_colorgradient_v100_qtz10_withQA.off
	#			${TESTING_DATA_DIR}/CompressionValenceExample/swirl_colorgradient_mepp2_polyhedron_v100_qtz10_withQA.ref.off
	#			-withCompr
	#			-withQA
	#			-maxV 100
	#			-Qbits 10
	#			-tolerance ${TEST_COMPRESSION_VALENCE_TOLERANCE})
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )

#####
add_executable( test_decompression_valence_lcc
                test_decompression_valence_lcc.cpp
				)
target_link_libraries( test_decompression_valence_lcc
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 2e-5)
#set(TEST_DECOMPRESSION_VALENCE_TOLERANCE 0)  # for exact comparison

if( NOT WIN32 )
	# this test crashes on Windows in Release mode. It passes in Debug mode.
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/noQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v100_qtz10_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v100_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
endif( NOT WIN32 )

if( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/noQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v1000_qtz10_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v1000_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v100_qtz12_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v100_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/noQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v1000_qtz12_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v1000_qtz12_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/noQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_lcc_monkey_v100_qtz10_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_noQA.ref.p3d
			  test_decompression_valence_lcc_monkey_v100_qtz10_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/noQA from mepp1
	add_test( Test_decompression_valence_lcc_monkey_v500_qtz15_noQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d
			  test_decompression_valence_lcc_monkey_v500_qtz15_noQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_noQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )

	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v100_qtz10_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz10/withQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v1000_qtz10_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v1000_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v100/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v100_qtz12_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v100_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v100_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_nefertiti.off/Surfacemesh/v1000/qtz12/withQA from mepp1
	add_test( Test_decompression_valence_lcc_nefertiti_v1000_qtz12_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d
			  test_decompression_valence_lcc_nefertiti_v1000_qtz12_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_nefertiti_unix_mepp1_v1000_qtz12_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v100/qtz10/withQA p3d file from mepp_2_
	#    compare with rgb_monkey.off compressed+decompressed by mepp_1_
	add_test( Test_decompression_valence_lcc_monkey_v100_qtz10_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp2_polyhedron_v100_qtz10_withQA.ref.p3d
			  test_decompression_valence_lcc_monkey_v100_qtz10_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v100_qtz10_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
	# test decompression of rgb_monkey.off/Surfacemesh/v500/qtz15/withQA from mepp1
	add_test( Test_decompression_valence_lcc_monkey_v500_qtz15_withQA
			  test_decompression_valence_lcc
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d
			  test_decompression_valence_lcc_monkey_v500_qtz15_withQA.off
			  ${TESTING_DATA_DIR}/CompressionValenceExample/rgb_monkey_unix_mepp1_v500_qtz15_withQA.ref.p3d.uncompressed.off
			  ${TEST_DECOMPRESSION_VALENCE_TOLERANCE}
			  )
endif( BUILD_ADVANCED_TESTS_COMPR_VALENCE )
