#!/bin/bash
#
# Reformat source code using clang-tidy
# through run-clang-tidy-6.0.py helper script.
#
# Ref: https://www.kdab.com/clang-tidy-part-1-modernize-source-code-using-c11c14/
#


# generate .clang-tidy config file
echo "# https://github.com/llvm-mirror/clang-tools-extra/blob/f126d29a57d3e020fc490633d608e377cb6fd81a/test/clang-tidy/readability-identifier-naming.cpp

Checks: '-*,readability-identifier-naming'
CheckOptions:
  - { key: readability-identifier-naming.LocalVariableCase,       value: lower_case }
  - { key: readability-identifier-naming.FunctionCase,            value: lower_case }
  - { key: readability-identifier-naming.ParameterCase,           value: lower_case }
  - { key: readability-identifier-naming.TemplateParameterCase,   value: CamelCase  }
" > ../.clang-tidy



# process a list of files globaly with run-clang-tidy-6.0.py
# $1: file list as "file1 file2 ..."
# $2: value for '-header-filter' option of clang-tidy
function apply_clang-tidy() {
	echo "applying clang-tidy to files:"
	echo "$1"

	# generate compile_commands.json file for run-clang-tidy-6.0.py
	echo -n '[' > compile_commands.json
	first_fname=on
	for fname in $1
	do
		if [ "$first_fname" = on ]
		then
			echo ""  # line return
		else
			echo ","  # , + line return
		fi

		echo -n "{
	  \"directory\": \"/tmp/MEPP2-elaptop-tmp/build/cltidy\",
	  \"command\": \"/usr/bin/clang++  -DBOOST_ALL_DYN_LINK -DCGAL_EIGEN3_ENABLED -DCGAL_NDEBUG -DFEVV_USE_AIF -DFEVV_USE_CGAL -DFEVV_USE_OPENMESH -DFEVV_USE_QT4 -DFEVV_USE_VTK -DQT_CORE_LIB -DQT_GUI_LIB -DQT_NO_DEBUG -DQT_OPENGL_LIB -DQT_XML_LIB -DUNIX -I/tmp/MEPP2-elaptop-tmp/External/CGAL-4.14 -I/usr/include/eigen3 -I/home/user1/OpenMesh-7.0/include -isystem /usr/include/qt4 -isystem /usr/include/qt4/QtOpenGL -isystem /usr/include/qt4/QtGui -isystem /usr/include/qt4/QtXml -isystem /usr/include/qt4/QtCore -I/tmp/MEPP2-elaptop-tmp/Visualization/osgQt/osg34and36/include -I/usr/include/vtk-7.1 -I/usr/include/freetype2 -I/usr/include/x86_64-linux-gnu -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi/opal/mca/event/libevent2022/libevent/include -I/usr/lib/x86_64-linux-gnu/openmpi/include -I/usr/include/python3.6m -I/usr/include/hdf5/openmpi -I/usr/include/jsoncpp -I/usr/include/libxml2 -I/usr/include/tcl -I/home/user1/pcl-1.8.1/include/pcl-1.8 -I/tmp/MEPP2-elaptop-tmp -I/tmp/MEPP2-elaptop-tmp/build/cltidy   -std=c++11 -O3 -DNDEBUG   -o dummy_$(basename $fname)  -c $fname\",
	  \"file\": \"$fname\"
	}"

		first_fname=off
	done >> compile_commands.json
	echo '' >> compile_commands.json
	echo ']' >> compile_commands.json

	# run clang-tidy
	run-clang-tidy-6.0.py  -header-filter="$2"  -checks='-*,readability-identifier-naming'  -fix
}

# we want to process only:
# - Examples
# - FEVV/Filters (except internal filter files)
# - FEVV/Wrappings
# - FEVV/Tools
# - Testing

############################################################
# first pass: process files in Examples,                   #
#             FEVV/Filters (except internal filter files), #
#             FEVV/Wrappings, FEVV/Tools                   #
############################################################

FILES="$(find  /tmp/MEPP2-elaptop-tmp/Examples  /tmp/MEPP2-elaptop-tmp/FEVV/Filters  /tmp/MEPP2-elaptop-tmp/FEVV/Wrappings  /tmp/MEPP2-elaptop-tmp/FEVV/Tools  /tmp/MEPP2-elaptop-tmp/FEVV/Types  -type f | grep -e '\.cpp$' -e '\.cxx$' -e '\.h$' -e '\.hpp$' -e '\.hxx$' -e '\.inl$')"

# prevent other directories/files to be changed
chmod -R a-w  ../Examples  ../External  ../FEVV  ../Private  ../Testing  ../Tools

# allow changes for targeted files and dirs
# Examples
chmod -R ug+w  ../Examples
# FEVV/Wrappings and /FEVV/Tools
chmod -R ug+w  ../FEVV/Wrappings  ../FEVV/Tools  ../FEVV/Types
# FEVV/Filters except internal filter files
FILTERS_TO_CLEAN="../FEVV/Operators/AIF/collapse_edge.hpp
../FEVV/Operators/AIF/split_edge.hpp
../FEVV/Filters/AIF/Functionals/similarity.h
../FEVV/Filters/AIF/Functionals/Predicates/topology.h
../FEVV/Filters/Generic/scaling.hpp
../FEVV/Filters/Generic/translation.hpp
../FEVV/Operators/Generic/Manifold/split_edge_euler.hpp
../FEVV/Operators/Generic/Manifold/collapse_edge.hpp
../FEVV/Filters/Generic/calculate_face_normals.hpp
../FEVV/Filters/Generic/Manifold/calculate_vertex_normals.hpp
../FEVV/Filters/Generic/minmax_map.h
../FEVV/Filters/Generic/print_points.hpp
../FEVV/Filters/Generic/reposition_vertices.hpp
../FEVV/Filters/Generic/Manifold/msdm2.h
../FEVV/Filters/Generic/Manifold/Curvature/curvature.hpp
../FEVV/Filters/Curvature/extract_Vpropres.h
../FEVV/Filters/Generic/Manifold/Curvature/extract_Vpropres.cpp
../FEVV/Operators/Generic/Manifold/split_edge.hpp
../FEVV/Filters/translation_with_object.h
../FEVV/Filters/Generic/generic_reader.hpp
../FEVV/Filters/Generic/tangents.hpp
../FEVV/Operators/Generic/Manifold/flip_edge.hpp
../FEVV/Filters/Generic/generic_writer.hpp
../FEVV/Filters/print_halfedges.h
../FEVV/Filters/Generic/print_face_normals.hpp
../FEVV/Filters/Generic/color_mesh.h
../FEVV/Filters/Generic/Manifold/calculate_one_ring_vertex_centroid.hpp
../FEVV/Filters/Generic/Manifold/Compression_Valence/compression_valence.h
../FEVV/Filters/Generic/Manifold/Compression_Valence/decompression_valence.h
../FEVV/Filters/Generic/AABB/compute_mesh_bounding_box.hpp
../FEVV/Operators/Geometry/triangles.hpp
../FEVV/Operators/Geometry/triangle_rad_angle.hpp
../FEVV/Operators/Generic/Manifold/k_ring.hpp
../FEVV/Operators/Generic/geometry.hpp"
for filter in $FILTERS_TO_CLEAN
do
	chmod ug+w  "$filter"
done


# process selected file
apply_clang-tidy  "$FILES"  "/tmp/MEPP2-elaptop-tmp/.*"

# restore file/dirs mode
chmod -R ug+w  ../Examples  ../External  ../FEVV  ../Private  ../Testing  ../Tools

# fix some globally missing replacements
./post_clang-tidy_fix.sh  ../FEVV/Tools/Math/MatrixOperations.hpp
./post_clang-tidy_fix.sh  ../FEVV/Operators/Geometry/triangles.hpp
./post_clang-tidy_fix.sh  ../FEVV/Filters/Compression_Valence/Compression_Valence_Component.inl
./post_clang-tidy_fix.sh  ../FEVV/DataStructures/AIF/AIFMeshHelpers.h
./post_clang-tidy_fix.sh  ../FEVV/Filters/Generic/Manifold/calculate_one_ring_vertex_centroid.hpp
./post_clang-tidy_fix.sh  ../FEVV/Wrappings/Geometry_traits_openmesh.h
./post_clang-tidy_fix.sh  ../Testing/Concepts/GeometryConceptCheck.h
./post_clang-tidy_fix.sh  ../Private/FEVV/Filters/Functionals/Surfaces/data_attached.h
./post_clang-tidy_fix.sh  ../Private/FEVV/Filters/global_repositioning.h

./post_clang-tidy_fix.sh  ../Testing/OpenMesh/test_geometry_traits_concept_openmesh.cpp
./post_clang-tidy_fix.sh  ../Testing/CGAL/test_geometry_traits_concept_linear_cell_complex.cpp
./post_clang-tidy_fix.sh  ../Testing/CGAL/test_geometry_traits_concept_surfacemesh.cpp
./post_clang-tidy_fix.sh  ../Testing/CGAL/test_geometry_traits_concept_polyhedron.cpp
./post_clang-tidy_fix.sh  ../Testing/AIF/test_matrix_operations_aif.cpp
./post_clang-tidy_fix.sh  ../Testing/AIF/test_geometry_traits_concept_aif.cpp

# fix some bad replacements or local misses
sed --in-place  's/scalar_type/ScalarType/g'  ../FEVV/Tools/IO/StringUtilities.h
sed --in-place  's/ob_jfile/obj_file/g'  ../FEVV/Filters/Generic/generic_reader.hpp
sed --in-place  's/nextVertexInFace/next_vertex_in_face/g'  ../Testing/Utils/utils_are_meshes_identical.cpp
sed --in-place  's/msd_m2/msdm2/g'  ../FEVV/Filters/Generic/Manifold/msdm2.h
sed --in-place  's/Vertex_iterator/VertexIterator/g'  ../FEVV/Filters/Generic/Manifold/msdm2.h
sed --in-place  's/ProcessMSDM2_Multires/process_msdm2_multires/g'  ../Visualization/PluginFilters/msdm2/MSDM2Plugin.h
sed --in-place  's/ComputeMinMaxVertices/compute_min_max_vertices/g'  ../Visualization/PluginFilters/msdm2/MSDM2Plugin.h
sed --in-place  's/ColorVerticesFromMap/color_vertices_from_map/g'  ../Visualization/PluginFilters/msdm2/MSDM2Plugin.h
sed --in-place  's/process_msd_m2_multires/process_msdm2_multires/g'  ../Examples/CGAL/example_msdm2.cpp
sed --in-place  's/msd_m2/msdm2/g'  ../Examples/CGAL/example_msdm2.cpp
sed --in-place  's/ProcessMSDM2_Multires/process_msdm2_multires/g'  ../Testing/CGAL/test_msdm2.cpp
sed --in-place  's/ComputeMinMaxVertices/compute_min_max_vertices/g'  ../Testing/CGAL/test_msdm2.cpp
sed --in-place  's/ColorVerticesFromMap/color_vertices_from_map/g'  ../Testing/CGAL/test_msdm2.cpp
sed --in-place  's/ComputeMinMaxVertices/compute_min_max_vertices/g'  ../Testing/CGAL/test_jnd.cpp
sed --in-place  's/ColorVerticesFromMap/color_vertices_from_map/g'  ../Testing/CGAL/test_jnd.cpp

sed --in-place  's/coordC_type/CoordCType/g'  ../FEVV/Tools/IO/VtkFileReader.h
sed --in-place  's/coordN_type/CoordNType/g'  ../FEVV/Tools/IO/VtkFileReader.h
sed --in-place  's/coord_type/CoordType/g'  ../FEVV/Tools/IO/VtkFileReader.h
sed --in-place  's/index_type/IndexType/g'  ../FEVV/Tools/IO/VtkFileReader.h
sed --in-place  's/polyData/poly_data/g'  ../FEVV/Tools/IO/VtkFileReader.h
sed --in-place  's/coordC_type/CoordCType/g'  ../FEVV/Tools/IO/VtkFileWriter.h
sed --in-place  's/coordN_type/CoordNType/g'  ../FEVV/Tools/IO/VtkFileWriter.h
sed --in-place  's/coord_type/CoordType/g'  ../FEVV/Tools/IO/VtkFileWriter.h
sed --in-place  's/index_type/IndexType/g'  ../FEVV/Tools/IO/VtkFileWriter.h


#########################################
# second pass: process files in Testing #
#########################################

FILES="$(find  /tmp/MEPP2-elaptop-tmp/Testing  -type f | grep -e '\.cpp$' -e '\.cxx$' -e '\.h$' -e '\.hpp$' -e '\.hxx$' -e '\.inl$' | grep -v -e "Testing/AIF/test_compression_valence_aif.cpp" -e "Testing/AIF/test_load_dangling_edge_aif.cpp"  -e "Testing/Visualization/testViewer_Features.cpp")"

# prevent other directories/files to be changed
chmod -R a-w  ../Examples  ../External  ../FEVV  ../Private  ../Testing  ../Tools

# allow changes for targeted files and dirs
chmod -R ug+w  ../Testing

# process selected file
apply_clang-tidy  "$FILES"  "/tmp/MEPP2-elaptop-tmp/Testing/.*"

# restore file/dirs mode
chmod -R ug+w  ../Examples  ../External  ../FEVV  ../Private  ../Testing  ../Tools

# fix some bad replacements or local misses
sed --in-place  's/screenSize/screen_size/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/adapterAIF/adapter_aif/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/viewerAIF/viewer_aif/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/adapterOpenMesh/adapter_open_mesh/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/viewerOpenMesh/viewer_open_mesh/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/_open_mesh/_openmesh/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/adapterLCC/adapter_lcc/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/viewerLCC/viewer_lcc/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/adapterSurface/adapter_surface/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/viewerSurface/viewer_surface/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/adapterPolyhedron/adapter_polyhedron/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/viewerPolyhedron/viewer_polyhedron/g'  ../Testing/Visualization/mepp.cpp

sed --in-place  's/mAIF/m_aif/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/mOpenMesh/m_openmesh/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/mPolyhedron/m_polyhedron/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/mSurface/m_surface/g'  ../Testing/Visualization/mepp.cpp
sed --in-place  's/mLCC/m_lcc/g'  ../Testing/Visualization/mepp.cpp


#####################################################
# third pass: fix field names in Material structure #
#####################################################

FILES="/tmp/MEPP2-elaptop-tmp/FEVV/Types/Material.h
/tmp/MEPP2-elaptop-tmp/FEVV/Tools/IO/ObjFileReader.h
/tmp/MEPP2-elaptop-tmp/FEVV/Tools/IO/ObjFileWriter.h
/tmp/MEPP2-elaptop-tmp/Visualization/MeshLoading.inl"

# fix targeted files
for fname in $FILES
do
	sed --in-place  's/ambient_redComponent/ambient_red_component/g'  "$fname"
	sed --in-place  's/ambient_greenComponent/ambient_green_component/g'  "$fname"
	sed --in-place  's/ambient_blueComponent/ambient_blue_component/g'  "$fname"
	sed --in-place  's/diffuse_redComponent/diffuse_red_component/g'  "$fname"
	sed --in-place  's/diffuse_greenComponent/diffuse_green_component/g'  "$fname"
	sed --in-place  's/diffuse_blueComponent/diffuse_blue_component/g'  "$fname"
	sed --in-place  's/specular_redComponent/specular_red_component/g'  "$fname"
	sed --in-place  's/specular_greenComponent/specular_green_component/g'  "$fname"
	sed --in-place  's/specular_blueComponent/specular_blue_component/g'  "$fname"
	sed --in-place  's/emissive_redComponent/emissive_red_component/g'  "$fname"
	sed --in-place  's/emissive_greenComponent/emissive_green_component/g'  "$fname"
	sed --in-place  's/emissive_blueComponent/emissive_blue_component/g'  "$fname"
done
