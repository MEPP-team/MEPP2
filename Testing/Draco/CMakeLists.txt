if( NOT BUILD_USE_DRACO OR NOT draco_FOUND )
  return()
endif()

include_directories( ${draco_INCLUDE_DIRS} )
add_executable(draco_encoder draco_encoder.cc)
target_link_libraries(draco_encoder ${draco_LIBRARIES})

add_test(draco_encoder draco_encoder -i ${TESTING_DATA_DIR}/armadillo_simplified.obj)

include_directories( ${draco_INCLUDE_DIRS} )
add_executable(draco_decoder draco_decoder.cc)
target_link_libraries(draco_decoder ${draco_LIBRARIES})

add_test(draco_decoder draco_decoder -i ${TESTING_DATA_DIR}/armadillo_simplified.drc)
