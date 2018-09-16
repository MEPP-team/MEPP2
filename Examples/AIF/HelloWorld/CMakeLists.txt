add_executable( hello_world_filter_aif
                hello_world_filter_aif.cpp )
target_link_libraries( hello_world_filter_aif
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
