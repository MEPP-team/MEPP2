if (UNIX)
  add_definitions(-DUNIX)
endif ()

if (WIN32)
  add_definitions(-DWIN32)
  # Prevent VC++ from defining min and max macros that prevent the proper
  # behavior of std::min() and std::max()
  add_definitions(-DNOMINMAX)
endif ()

if (APPLE)
  add_definitions(-DAPPLE)
endif ()
