if( NOT BUILD_EXAMPLES OR NOT BUILD_USE_CGAL )
  return()
endif()

include(${CGAL_USE_FILE})

#####################################################################
add_subdirectory( HelloWorld )
add_subdirectory( JND )

#####

add_executable( example_msdm2
                "${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/msdm2.h"
                example_msdm2.cpp
                "${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/Curvature/extract_Vpropres.cpp")
target_link_libraries( example_msdm2
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
if( BUILD_USE_VTK )
  target_link_libraries( example_msdm2
						 ${VTK_LIBRARIES} )
endif()
