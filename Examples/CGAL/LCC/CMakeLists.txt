######
add_executable( hello_world_filter_cgal_linear_cell_complex
                hello_world_filter_cgal_linear_cell_complex.cpp )
target_link_libraries( hello_world_filter_cgal_linear_cell_complex
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
