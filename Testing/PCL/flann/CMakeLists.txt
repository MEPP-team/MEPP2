if( NOT BUILD_USE_PCL OR NOT FLANN_FOUND )
  return()
endif()

include_directories( ${FLANN_INCLUDE_DIRS} )
add_definitions( ${FLANN_DEFINITIONS} )

# the section below is disabled because of link error on GitLab-Ci/Linux
# see https://gitlab.liris.cnrs.fr/MEPP-team/MEPP2/pipelines/1790
#ELO- add_executable(flann_example flann_example.cpp)
#ELO- target_link_libraries(flann_example ${FLANN_LIBRARIES})
#ELO- 
#ELO- add_test(flann_example flann_example ${TESTING_DATA_DIR}/flann/dataset.dat ${TESTING_DATA_DIR}/flann/testset.dat)
