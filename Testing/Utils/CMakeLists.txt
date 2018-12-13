
add_executable( compare_mesh_files
                compare_mesh_files.cpp
				"${PROJECT_SOURCE_DIR}/Testing/Utils/utils_are_meshes_identical.cpp" )
target_link_libraries( compare_mesh_files
                       ${Boost_LIBRARIES} ${FBX_LIBRARY} ${VTK_LIBRARIES}
                       )
