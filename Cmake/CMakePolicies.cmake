if (COMMAND cmake_policy)

  # MACOSX_RPATH is enabled by default.
  if( POLICY CMP0042 )
    cmake_policy(SET CMP0042 NEW)
  endif ()

endif ()

# ---

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON) # see cmake POLICY CMP0042
endif(APPLE)
