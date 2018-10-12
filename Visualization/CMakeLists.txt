if( NOT BUILD_USE_GUI )
  return()
endif()

##### Check (optional) GUI Dependencies (OSG and Qt)
include( GUIOptionalDependencies )

add_subdirectory( Applications )
