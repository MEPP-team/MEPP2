#!/bin/bash
#
# Fix files not fully precessed by clang-tidy.
#

if [ -z "$*" ]
then
	echo "Usage:  $0  file1 [file2 [...]]"
	exit 1
fi

# apply provided replacements list to file
#  $1: replacements list as "str1  new_str1  str2  new_str2 ..."
#  $2: name of the file to process
function apply_fix() {
	echo "$1" | sed 's/ /\n/g' | while read old_str
	do
		read new_str
		echo "apply post-fix  $old_str -> $new_str  to file $2"
		sed --in-place  "s/\b$old_str\b/$new_str/g"  "$2"
	done
}


# functions renamed in FEVV/Tools/Math/MatrixOperations.hpp
FIXES="Add add
AreAligned are_aligned
AreCollinear are_collinear
CrossProduct cross_product
DotProduct dot_product
getAngleFrom3positions get_angle_from3positions
getAngleFromNonUnitVectors get_angle_from_non_unit_vectors
getAngleFromUnitVectors get_angle_from_unit_vectors
getAngleInDegreeFrom3positions get_angle_in_degree_from3positions
getAngleInDegreeFromNonUnitVectors get_angle_in_degree_from_non_unit_vectors
getAngleInDegreeFromUnitVectors get_angle_in_degree_from_unit_vectors
IsDiagonal is_diagonal
L2Distance l2_distance
meanSqrt mean_sqrt
meanSqrtSqrt mean_sqrt_sqrt
Normalize normalize
projCurv proj_curv
rotCoordSys rot_coord_sys
ScalarMult scalar_mult
sortVectorIndices sort_vector_indices
Sub sub
Transformation transformation
VectorTimesTransposeMult vector_times_transpose_mult
weightedMean weighted_mean"

# functions renamed in FEVV/Tools/Math/MatrixOperations.hpp
FIXES="$FIXES
triangleBarycenter triangle_barycenter
trianglePerimeter triangle_perimeter
triangleNormalUnnormalized triangle_normal_unnormalized
triangleArea triangle_area
triangleShapePotential triangle_shape_potential"

# functions renamed in FEVV/Wrappings/Geometry_traits_aif.h
FIXES="$FIXES
addP add_p
subP sub_p"

# functions renamed in FEVV/Wrappings/Geometry_traits_operators.h
FIXES="$FIXES
AddV add_v
AddP add_p
SubP sub_p
Sub sub
DotProduct dot_product
VerifyAssumption verify_assumption
VerifyOperatorResults verify_operator_results"

for f in $*
do
	apply_fix  "$FIXES"  "$f"
done





