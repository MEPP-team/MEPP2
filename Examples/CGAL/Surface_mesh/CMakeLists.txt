
######
add_executable( helloworld_filter_cgal_surface_mesh
                helloworld_filter_cgal_surface_mesh.cpp )
target_link_libraries( helloworld_filter_cgal_surface_mesh
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )

######
add_executable( example_msdm2
                "${PROJECT_SOURCE_DIR}/FEVV/Filters/CGAL/Surface_mesh/msdm2.h"
                example_msdm2.cpp
                )
target_link_libraries( example_msdm2
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
if( BUILD_USE_VTK )
  target_link_libraries( example_msdm2
						 ${VTK_LIBRARIES} )
endif()

######
add_subdirectory( JND )
