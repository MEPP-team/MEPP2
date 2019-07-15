#####
add_executable( test_reader_writer_cgal_point_set
                test_reader_writer_cgal_point_set.cpp )
target_link_libraries( test_reader_writer_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
add_test(Test_reader_writer_cgal_point_set
         test_reader_writer_cgal_point_set
         ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
         ${TESTING_DATA_DIR}/test_reader_writer_cgal_point_set.ply
         ${TESTING_DATA_DIR}/tetra.with_normal_and_color.ply
         )
