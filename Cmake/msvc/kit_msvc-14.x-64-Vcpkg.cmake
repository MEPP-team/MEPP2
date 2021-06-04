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

set(CGAL_DIR					${MSVC_KIT_ROOT}/share/cgal)
# -------
# for GMP
set(GMP_INCLUDE_DIR             ${MSVC_KIT_ROOT}/include)
set(GMP_LIBRARIES               ${MSVC_KIT_ROOT}/lib/mpir.lib)
# for MPFR
set(MPFR_INCLUDE_DIR            ${MSVC_KIT_ROOT}/include)
set(MPFR_LIBRARIES              ${MSVC_KIT_ROOT}/lib/mpfr.lib)
# -------

set(OPENMESH_DIR				${MSVC_KIT_ROOT})

set(EIGEN3_INCLUDE_DIR			${MSVC_KIT_ROOT}/include/eigen3)
set(EIGEN_INCLUDE_DIR			${MSVC_KIT_ROOT}/include/eigen3)

# ---------------------------------------------------------------------------------
set(IMG_DIR_3rdParty            ${MSVC_KIT_ROOT})

set(JPEG_INCLUDE_DIR            ${IMG_DIR_3rdParty})
set(JPEG_LIBRARY                optimized ${IMG_DIR_3rdParty}/lib/jpeg.lib debug ${IMG_DIR_3rdParty}/debug/lib/jpegd.lib)

# for PNG
set(ZLIB_INCLUDE_DIR            ${IMG_DIR_3rdParty})
set(ZLIB_LIBRARY                optimized ${IMG_DIR_3rdParty}/lib/zlib.lib debug ${IMG_DIR_3rdParty}/debug/lib/zlibd.lib)

set(PNG_PNG_INCLUDE_DIR         ${IMG_DIR_3rdParty})
set(PNG_LIBRARY                 optimized ${IMG_DIR_3rdParty}/lib/libpng16.lib debug ${IMG_DIR_3rdParty}/debug/lib/libpng16d.lib)

set(TIFF_INCLUDE_DIR            ${IMG_DIR_3rdParty})
set(TIFF_LIBRARY                optimized ${IMG_DIR_3rdParty}/lib/tiff.lib debug ${IMG_DIR_3rdParty}/debug/lib/tiffd.lib)
# ---------------------------------------------------------------------------------

### 'addon 01' : Qt4, OpenSceneGraph
# and
### 'addon 02' : Qt5

#set(BUILD_USE_QT5				TRUE) # caution: no Qt4 with vcpkg !

#if(BUILD_USE_QT5)
	# with qt5
	#set(QT5_DIR					${MSVC_KIT_ROOT}/../Qt/Qt5.6.3/5.6.3/msvc2015_64)
	set(QT5_DIR					${MSVC_KIT_ROOT})
#else(BUILD_USE_QT5)
	# with qt4
#	set(QTDIR					${MSVC_KIT_ROOT}/../Qt/qt-4.8.7-x64-msvc2015)
#endif(BUILD_USE_QT5)

set(OSG_DIR						${MSVC_KIT_ROOT})

### addon 03 : PCL, FLANN

set(PCL_DIR						${MSVC_KIT_ROOT}/share/pcl)

set(FLANN_INCLUDE_DIR			${MSVC_KIT_ROOT}/include/flann)
set(FLANN_LIBRARY				${MSVC_KIT_ROOT}/lib/flann.lib)

# TEMP FOR VCPKG 2020.11-1 - BECAUSE NO PB WITH 2019.12 !
set(LZ4-LIB-TMP					"optimized;${MSVC_KIT_ROOT}/lib/lz4.lib;debug;${MSVC_KIT_ROOT}/debug/lib/lz4d.lib")

### addon 04 : VTK

set(VTK_DIR						${MSVC_KIT_ROOT}/share/vtk)

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
set(OPENSSL_ROOT_DIR			${MSVC_KIT_ROOT}/../OpenSSL-1.1.1d) # not yet done under VCPKG

### vtests

set(VTEST						  -test)

###
