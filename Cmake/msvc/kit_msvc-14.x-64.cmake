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
	set(BOOST_ROOT				${MSVC_KIT_ROOT}/boost_1_67_0_V141) 	# VS2019 x64
	#set(BOOST_ROOT				${MSVC_KIT_ROOT}/boost_1_71_0_V142) 	# VS2019 x64
elseif(MSVC_VERSION GREATER_EQUAL 1910)
	set(BOOST_ROOT				${MSVC_KIT_ROOT}/boost_1_67_0_V141) 	# VS2017 x64
else()
	set(BOOST_ROOT				${MSVC_KIT_ROOT}/boost_1_59_0)			# VS2015 x64
	#set(BOOST_ROOT				${MSVC_KIT_ROOT}/boost_1_67_0_V140) 	# VS2015 x64
endif()

set(CGAL_DIR					${MSVC_KIT_ROOT}/CGAL-4.14-hdr-only)
# -------
# for GMP
set(GMP_INCLUDE_DIR             ${MSVC_KIT_ROOT}/gmp/include)
set(GMP_LIBRARIES               ${MSVC_KIT_ROOT}/gmp/lib/libgmp-10.lib)
# for MPFR
set(MPFR_INCLUDE_DIR            ${MSVC_KIT_ROOT}/gmp/include)
set(MPFR_LIBRARIES              ${MSVC_KIT_ROOT}/gmp/lib/libmpfr-4.lib)
# -------
message("--> CGAL_DIR used : ${CGAL_DIR}")

set(OPENMESH_DIR				${MSVC_KIT_ROOT}/OpenMesh-6.2)

set(EIGEN3_INCLUDE_DIR			${MSVC_KIT_ROOT}/eigen-3.2.8)

# ---------------------------------------------------------------------------------
set(IMG_DIR_3rdParty            ${MSVC_KIT_ROOT}/img-3rdparty)

set(JPEG_INCLUDE_DIR            ${IMG_DIR_3rdParty}/libjpeg)
set(JPEG_LIBRARY                ${IMG_DIR_3rdParty}/build/lib/Release/jpeg.lib)

# for PNG
set(ZLIB_INCLUDE_DIR            ${IMG_DIR_3rdParty}/zlib)
set(ZLIB_LIBRARY                ${IMG_DIR_3rdParty}/build/lib/Release/zlib.lib)

set(PNG_PNG_INCLUDE_DIR         ${IMG_DIR_3rdParty}/libpng)
set(PNG_LIBRARY                 optimized ${IMG_DIR_3rdParty}/build/lib/Release/libpng.lib debug ${IMG_DIR_3rdParty}/build/lib/Debug/libpngd.lib)

set(TIFF_INCLUDE_DIR            ${IMG_DIR_3rdParty}/libtiff/libtiff)
set(TIFF_LIBRARY                ${IMG_DIR_3rdParty}/build/lib/Release/libtiff.lib)
# ---------------------------------------------------------------------------------

### 'addon 01' : Qt4, OpenSceneGraph
# and
### 'addon 02' : Qt5

#set(BUILD_USE_QT5				TRUE)

if(BUILD_USE_QT5)
	# with qt5
	set(QT5_DIR					${MSVC_KIT_ROOT}/Qt/Qt5.6.3/5.6.3/msvc2015_64)
else(BUILD_USE_QT5)
	# with qt4
	set(QTDIR					${MSVC_KIT_ROOT}/Qt/qt-4.8.7-x64-msvc2015)
endif(BUILD_USE_QT5)

set(OSG_DIR						${MSVC_KIT_ROOT}/osg/OpenSceneGraph-3.6.4)

### addon 03 : PCL, FLANN

set(PCL_DIR						${MSVC_KIT_ROOT}/PCL/pcl-1.7.2/cmake)

set(FLANN_INCLUDE_DIR			${MSVC_KIT_ROOT}/PCL/flann-1.8.4/include)
set(FLANN_LIBRARY				${MSVC_KIT_ROOT}/PCL/flann-1.8.4/lib/flann.lib)

### addon 04 : VTK

set(VTK_DIR						${MSVC_KIT_ROOT}/kitware/VTK-7.0.0)

### addon 05 : FBX

set(FBX_DIR						${MSVC_KIT_ROOT}/FBX_SDK/2019.0)

### addon 06 : Draco

set(DRACO_DIR					${MSVC_KIT_ROOT}/draco-1.3.4)

### vtests

set(VTEST						test)

###
