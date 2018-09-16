# MSVC - KIT VS 2015 x64 - V05

if( DEFINED ENV{MSVC_KIT_ROOT} )
	set( MSVC_KIT_ROOT $ENV{MSVC_KIT_ROOT} )
endif()
if( NOT DEFINED MSVC_KIT_ROOT )
	message(FATAL_ERROR "MSVC_KIT_ROOT not set.  Please set MSVC_KIT_ROOT.")
endif()

### 'core' : boost, CGAL, OpenMesh, Eigen 3

set(BOOST_ROOT			${MSVC_KIT_ROOT}/boost_1_59_0)

set(CGAL_DIR			${MSVC_KIT_ROOT}/CGAL-4.11-bug-patched)

set(OPENMESH_DIR		${MSVC_KIT_ROOT}/OpenMesh-6.2)

set(EIGEN3_INCLUDE_DIR	${MSVC_KIT_ROOT}/eigen-3.2.8)

### 'addon 01' : Qt4, OpenSceneGraph
# and
### 'addon 02' : Qt5

#set(BUILD_USE_QT5			TRUE)

if(BUILD_USE_QT5)
	# with qt5
	set(QT5_DIR			${MSVC_KIT_ROOT}/Qt/Qt5.6.0/5.6/msvc2015_64)
else(BUILD_USE_QT5)
	# with qt4
	set(QTDIR			${MSVC_KIT_ROOT}/Qt/qt-4.8.7-x64-msvc2015)
endif(BUILD_USE_QT5)

set(OSG_DIR				${MSVC_KIT_ROOT}/osg/OpenSceneGraph-3.4.0)

### addon 03 : PCL, FLANN

set(PCL_DIR				${MSVC_KIT_ROOT}/PCL/pcl-1.7.2/cmake)

set(FLANN_INCLUDE_DIR	${MSVC_KIT_ROOT}/PCL/flann-1.8.4/include)
set(FLANN_LIBRARY		${MSVC_KIT_ROOT}/PCL/flann-1.8.4/lib/flann.lib)

### addon 04 : VTK

set(VTK_DIR				${MSVC_KIT_ROOT}/kitware/VTK-7.0.0)

### addon 05 : FBX

set(FBX_DIR				${MSVC_KIT_ROOT}/FBX_SDK/2019.0)

### vtests

set(VTEST				test)

###
