if( NOT BUILD_USE_CGAL )
  return()
endif()

include(${CGAL_USE_FILE})

if( MSVC )
  # remove warning
  # boost_1_59_0\boost\property_tree\detail\json_parser\standard_callbacks.hpp(132): warning C4715: 'boost::property_tree::json_parser::detail::standard_callbacks<boost::property_tree::basic_ptree<std::basic_string<char,std::char_traits<char>,std::allocator<char> >,std::basic_string<char,std::char_traits<char>,std::allocator<char> >,std::less<std::basic_string<char,std::char_traits<char>,std::allocator<char> > > > >::new_tree'ÿ: les chemins de contr“le ne retournent pas tous une valeur
  # see https://stackoverflow.com/questions/47136503/boost-json-parser-warning-c4715
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4715" )
endif()

# In order to make sure the FEVV embarked Euler operators will be seeked
# before the user/system install of CGAL. Note that the BEFORE flag is
# just to make sure to circumvent any possible CGAL black magic.

# include_directories( BEFORE "External/CGAL-4.14" ) # Perhaps to remove: here we want to use CGAL and not our version of Euler operator: TO DISCUSS

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
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/utils/lightsampler.hpp"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatfrequencycomputor.hpp"
	"${PROJECT_SOURCE_DIR}/FEVV/Filters/Generic/Manifold/JustNoticeableDistortion/perception/flatcontrastcomputor.h")

add_executable( example_just_noticeable_distortion
                example_just_noticeable_distortion.cpp
                ${JND_HEADER_FILES}   ### For VC++ to present include files
                )

target_link_libraries( example_just_noticeable_distortion
                       ${CGAL_LIBRARY} ${CGAL_3RD_PARTY_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
if( BUILD_USE_VTK )
	target_link_libraries( example_just_noticeable_distortion
						            ${VTK_LIBRARIES}
						   )
endif()
