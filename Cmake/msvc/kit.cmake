#if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19.10.25017.0") # 19.10.25017.0 -> info from boost 1_64_0 sources
#	message(FATAL_ERROR "You must choose 'Visual Studio 14 2015 Win64' !") # VS2017 64 not supported for now !
#elseif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19")
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19")
	if (CMAKE_SIZEOF_VOID_P EQUAL 8) # 64 bits
		if (NOT DEFINED ENV{APPVEYOR})
			message("--> VS2015 64 KIT used (MSVC ${CMAKE_CXX_COMPILER_VERSION})")
			include("${PROJECT_SOURCE_DIR}/${CMAKE_INCLUDE_INVOCATION_DIR}/Cmake/msvc/kit_vs2015-64_v01.cmake") # vs2015-64 (actual kit)
		else()
			message("--> VS2015 64 KIT used (MSVC ${CMAKE_CXX_COMPILER_VERSION}) - AppVeyor packages")
			include("${PROJECT_SOURCE_DIR}/${CMAKE_INCLUDE_INVOCATION_DIR}/Cmake/msvc/kit_vs2015-64_v01-AppVeyor.cmake") # vs2015-64 (actual kit)
		endif()
	elseif (CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bits
		message(FATAL_ERROR "You must choose WIN64 architecture : 'Visual Studio 14 2015 Win64' !")
	endif()
else()
	message(FATAL_ERROR "You must choose 'Visual Studio 14 2015 Win64' !")
endif()