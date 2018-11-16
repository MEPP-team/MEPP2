#if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19.10.25017.0") # 19.10.25017.0 -> info from boost 1_64_0 sources
#	message(FATAL_ERROR "You must choose 'Visual Studio 14 2015 Win64' !") # VS2017 x64 not supported for now !
#elseif (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19")
if (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "19")
	if (CMAKE_SIZEOF_VOID_P EQUAL 8) # 64 bits
		if (NOT DEFINED ENV{APPVEYOR})
			message("--> MSVC-14.x-64 KIT used (MSVC ${CMAKE_CXX_COMPILER_VERSION})")
			include("${CMAKE_CURRENT_LIST_DIR}/kit_msvc-14.x-64.cmake")
		else()
			message("--> MSVC-14.x-64 KIT used (MSVC ${CMAKE_CXX_COMPILER_VERSION}) --> with AppVeyor packages")
			include("${CMAKE_CURRENT_LIST_DIR}/kit_msvc-14.x-64-AppVeyor.cmake")
		endif()
	elseif (CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bits
		message(FATAL_ERROR "You must choose WIN64 architecture : 'Visual Studio 14 2015 Win64' or 'Visual Studio 15 2017 Win64' !")
	endif()
else()
	message(FATAL_ERROR "You must choose 'Visual Studio 14 2015 Win64' or 'Visual Studio 15 2017 Win64' !")
endif()
