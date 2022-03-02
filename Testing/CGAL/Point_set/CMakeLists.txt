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
if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "5.0" )
  # CGAL < 5.0
  add_test( Test_reader_writer_cgal_point_set
            test_reader_writer_cgal_point_set
            ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
            test_reader_writer_cgal_point_set.ply
            ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
            )
else()
  # CGAL >= 5.0
  add_test( Test_reader_writer_cgal_point_set_1
            test_reader_writer_cgal_point_set
            ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
            test_reader_writer_cgal_point_set_1.ply
            ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_reader_writer_cgal_point_set_1.cgal-5.0.2.ref.ply
            )
  add_test( Test_reader_writer_cgal_point_set_2
            test_reader_writer_cgal_point_set
            ${TESTING_DATA_DIR}/tetra.binary.ply
            test_reader_writer_cgal_point_set_2.ply
            ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_reader_writer_cgal_point_set_2.ref.ply
            )
endif()

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
if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "5.0" )
  # CGAL < 5.0
  add_test( Test_point_cloud_curvature_cgal_point_set
            test_point_cloud_curvature_cgal_point_set
            ${TESTING_DATA_DIR}/casting_with_normals.ply
            test_point_cloud_curvature_cgal_point_set_casting.ply
            15
            ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_point_cloud_curvature_cgal_point_set_casting.ref.ply
            )
else()
  # CGAL >= 5.0
  add_test( Test_point_cloud_curvature_cgal_point_set
            test_point_cloud_curvature_cgal_point_set
            ${TESTING_DATA_DIR}/casting_with_normals.ply
            test_point_cloud_curvature_cgal_point_set_casting.ply
            15
            ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_point_cloud_curvature_cgal_point_set_casting.cgal-5.0.2.ref.ply
            )
endif()

#####
add_executable( test_geometry_concept_cgal_point_set
                test_geometry_concept_cgal_point_set.cpp )
target_link_libraries( test_geometry_concept_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
add_test( Test_geometry_concept_cgal_point_set
          test_geometry_concept_cgal_point_set )

#####
add_executable( test_copy_graph_cgal_point_set
                test_copy_graph_cgal_point_set.cpp
				        )
target_link_libraries( test_copy_graph_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES}
                       )
if( ${CGAL_MAJOR_VERSION}.${CGAL_MINOR_VERSION} VERSION_LESS "5.0" )
  # CGAL < 5.0
  # TODO fix test for CGAL >= 5.0
  add_test( Test_copy_graph_cgal_point_set_monkey_color_normal
            test_copy_graph_cgal_point_set
            ${TESTING_DATA_DIR}/monkey_points_with_color_and_normal.ply
            test_copy_graph_cgal_point_set_monkey_color_normal.ply
            ${TESTING_REFOUTPUT_DIR}/CGAL/Point_set/test_copy_graph_cgal_point_set_monkey_color_normal.ref.ply
            )
endif()

#####
if( BUILD_USE_PCL )
    add_executable( test_point_cloud_normal_wpca_cgal_point_set
                    test_point_cloud_normal_wpca_cgal_point_set.cpp )
    target_link_libraries( test_point_cloud_normal_wpca_cgal_point_set
                           ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                           ${Boost_LIBRARIES}
                           ${PCL_COMMON_LIBRARIES} ${PCL_IO_LIBRARIES}
                           ${PCL_LIBRARIES}
                           )
    # non-crashing test only because of reproducibility issue, cf issue #543
    add_test( Test_point_cloud_wpca_cgal_point_set
              test_point_cloud_normal_wpca_cgal_point_set
              ${TESTING_DATA_DIR}/tetra.xyz
              test_point_cloud_wpca_cgal_point_set_tetra.ply
              3
              0
              0.1
              )
endif()
