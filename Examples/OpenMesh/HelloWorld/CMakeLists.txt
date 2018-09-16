add_executable( hello_world_filter_openmesh
                hello_world_filter_openmesh.cpp )
target_link_libraries( hello_world_filter_openmesh
                       ${OPENMESH_LIBRARIES}
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
