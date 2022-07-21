
#####################################################################

######
add_library( msdm2_surfacemesh SHARED msdm2_surfacemesh.cpp )
target_link_libraries( msdm2_surfacemesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       ${IMG-3RDPARTY_LIBRARIES}
                       )

######
add_executable( example_libmsdm2_surfacemesh
                example_libmsdm2_surfacemesh.cpp )
target_link_libraries( example_libmsdm2_surfacemesh
                       msdm2_surfacemesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       ${IMG-3RDPARTY_LIBRARIES}
                       )

######
