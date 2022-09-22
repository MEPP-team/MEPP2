if( NOT BUILD_USE_DRACO OR NOT draco_FOUND )
  return()
endif()

include_directories( ${draco_INCLUDE_DIRS} )
add_executable(test_draco_encoder draco_encoder.cc)
target_link_libraries(test_draco_encoder ${Boost_LIBRARIES} ${draco_LIBRARIES})

if( (NOT DEFINED ENV{VCPKG}) AND (NOT DEFINED VCPKG) ) # TEMP FOR VCPKG 2020.11-1 !
  add_test(Test_draco_encoder test_draco_encoder -i ${TESTING_DATA_DIR}/armadillo_simplified.obj)
endif()

include_directories( ${draco_INCLUDE_DIRS} )
add_executable(test_draco_decoder draco_decoder.cc)
target_link_libraries(test_draco_decoder ${Boost_LIBRARIES} ${draco_LIBRARIES})

if( (NOT DEFINED ENV{VCPKG}) AND (NOT DEFINED VCPKG) ) # TEMP FOR VCPKG 2020.11-1 !
  add_test(Test_draco_decoder test_draco_decoder -i ${TESTING_DATA_DIR}/armadillo_simplified.drc)
endif()
