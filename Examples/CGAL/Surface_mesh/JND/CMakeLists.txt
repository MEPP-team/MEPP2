if( NOT BUILD_USE_CGAL )
  return()
endif()

include(${CGAL_USE_FILE})

if( MSVC )
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /IGNORE:4715" )
endif()

# In order to make sure the FEVV embarked Euler operators will be seeked
# before the user/system install of CGAL. Note that the BEFORE flag is
# just to make sure to circumvent any possible CGAL black magic.

# include_directories( BEFORE "External/CGAL-4.11" ) # Perhaps to remove: here we want to use CGAL and not our version of Euler operator: TO DISCUSS

####

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/FEVV/Filters/Just_Noticeable_Distortion/)

SET(JND_HEADER_FILES
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/just_noticeable_distortion.hpp"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/genericparametricmodel.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/contrastsensitivity.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/contrastmasking.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/psychometricfunction.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/threshold.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/types.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/models/visibility.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/utils/lightsampler.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatfrequencycomputor.h"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatcontrastcomputor.h")

SET(JND_CPP_FILES
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/utils/lightsampler.cpp"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatfrequencycomputor.cpp")

add_executable( example_just_noticeable_distortion
                ${JND_HEADER_FILES}   ### For VC++ to present include files
                example_just_noticeable_distortion.cpp
                ${JND_CPP_FILES} )

target_link_libraries( example_just_noticeable_distortion
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
if( BUILD_USE_VTK )
	target_link_libraries( example_just_noticeable_distortion
						            ${VTK_LIBRARIES}
						   )
endif()
