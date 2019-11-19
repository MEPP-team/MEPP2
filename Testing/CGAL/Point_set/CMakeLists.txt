#####
add_executable( test_complying_concepts_cgal_point_set
                test_complying_concepts_cgal_point_set.cpp )
target_link_libraries( test_complying_concepts_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test( Test_complying_concepts_cgal_point_set
          test_complying_concepts_cgal_point_set )

#####
add_executable( test_reader_writer_cgal_point_set
                test_reader_writer_cgal_point_set.cpp )
target_link_libraries( test_reader_writer_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
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
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test( Test_property_maps_cgal_point_set
          test_property_maps_cgal_point_set
          )
