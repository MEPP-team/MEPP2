# MSVC - KIT msvc-14.x-64

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

### 'core' : boost, CGAL, OpenMesh, Eigen 3, Img-3rdparty: jpeg, zlib (for png), png, tiff

# from https://github.com/Kitware/CMake/blob/master/Modules/Platform/Windows-MSVC.cmake
if(MSVC_VERSION GREATER_EQUAL 1920)
	set(BOOST_ROOT				${MSVC_KIT_ROOT})		# VS2019 x64
elseif(MSVC_VERSION GREATER_EQUAL 1910)
	set(BOOST_ROOT				${MSVC_KIT_ROOT}) 		# VS2017 x64
else()
	#set(BOOST_ROOT				${MSVC_KIT_ROOT})		# VS2015 x64
endif()

set(CGAL_DIR					${MSVC_KIT_ROOT})
# -------
# for GMP
set(GMP_INCLUDE_DIR             ${MSVC_KIT_ROOT}/include)
set(GMP_LIBRARIES               ${MSVC_KIT_ROOT}/lib/mpir.lib)
# for MPFR
set(MPFR_INCLUDE_DIR            ${MSVC_KIT_ROOT}/include)
set(MPFR_LIBRARIES              ${MSVC_KIT_ROOT}/lib/mpfr.lib)
# -------

set(OPENMESH_DIR				${MSVC_KIT_ROOT})

set(EIGEN3_INCLUDE_DIR			${MSVC_KIT_ROOT}/include)

### 'addon 01' : Qt4, OpenSceneGraph
# and
### 'addon 02' : Qt5

#set(BUILD_USE_QT5				TRUE) # caution: no Qt4 with vcpkg !

if(BUILD_USE_QT5)
	# with qt5
	#set(QT5_DIR					${MSVC_KIT_ROOT}/../Qt/Qt5.6.3/5.6.3/msvc2015_64)
	set(QT5_DIR					${MSVC_KIT_ROOT})
else(BUILD_USE_QT5)
	# with qt4
	set(QTDIR					${MSVC_KIT_ROOT}/../Qt/qt-4.8.7-x64-msvc2015)
endif(BUILD_USE_QT5)

set(OSG_DIR						${MSVC_KIT_ROOT})

### addon 03 : PCL, FLANN

set(PCL_DIR						${MSVC_KIT_ROOT})

### addon 04 : VTK

set(VTK_DIR						${MSVC_KIT_ROOT})

### addon 05 : FBX

set(FBX_DIR						${MSVC_KIT_ROOT}/../FBX_SDK/2019.0)

### addon 06 : Draco

set(DRACO_DIR					${MSVC_KIT_ROOT})
# the 2 lines below are mandatory else release lib is not found, don't understand why... 
set(draco_LIBRARIES_RELEASE		${MSVC_KIT_ROOT}/lib/draco.lib)
set(draco_LIBRARIES_DEBUG		${MSVC_KIT_ROOT}/debug/lib/draco.lib)

### Boost.Beast and OpenSSL - Early support ! Need Boost >= 1.70

set(BUILD_USE_BOOST_BEAST		FALSE)
set(BUILD_USE_OPENSSL			TRUE)
set(OPENSSL_ROOT_DIR			${MSVC_KIT_ROOT}/../OpenSSL-1.1.1d)

### vtests

set(VTEST						  -test)

###
