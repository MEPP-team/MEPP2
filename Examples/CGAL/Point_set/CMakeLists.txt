######
add_executable( compute_normals_cgal_point_set
                compute_normals_cgal_point_set.cpp )
target_link_libraries( compute_normals_cgal_point_set
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
