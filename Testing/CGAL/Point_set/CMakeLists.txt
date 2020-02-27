#####
add_executable( test_complying_concepts_cgal_point_set
                test_complying_concepts_cgal_point_set.cpp )
target_link_libraries( test_complying_concepts_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_complying_concepts_cgal_point_set
          test_complying_concepts_cgal_point_set )

#####
add_executable( test_reader_writer_cgal_point_set
                test_reader_writer_cgal_point_set.cpp )
target_link_libraries( test_reader_writer_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_reader_writer_cgal_point_set
          test_reader_writer_cgal_point_set
          ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
          test_reader_writer_cgal_point_set.ply
          ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
          )

#####
add_executable( test_property_maps_cgal_point_set
                test_property_maps_cgal_point_set.cpp )
target_link_libraries( test_property_maps_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_property_maps_cgal_point_set
          test_property_maps_cgal_point_set
          )

#####
add_executable( test_NN_search_concept_cgal_point_set
                test_NN_search_concept_cgal_point_set.cpp )
target_link_libraries( test_NN_search_concept_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_NN_search_concept_cgal_point_set
          test_NN_search_concept_cgal_point_set
          )

#####
add_executable( test_point_cloud_curvature_cgal_point_set
                test_point_cloud_curvature_cgal_point_set.cpp )
target_link_libraries( test_point_cloud_curvature_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_point_cloud_curvature_cgal_point_set
          test_point_cloud_curvature_cgal_point_set
          ${TESTING_DATA_DIR}/monkey_points_with_color_and_normal.ply
          test_point_cloud_curvature_cgal_point_set_monkey.ply
          15
          ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_point_cloud_curvature_cgal_point_set_monkey.ref.ply
          )

#####
add_executable( test_geometry_concept_cgal_point_set
                test_geometry_concept_cgal_point_set.cpp )
target_link_libraries( test_geometry_concept_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_geometry_concept_cgal_point_set
          test_geometry_concept_cgal_point_set )
