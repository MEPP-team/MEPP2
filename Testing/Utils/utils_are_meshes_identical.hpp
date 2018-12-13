#pragma once

#include <string>


typedef float FloatT;


/*
 * \brief  Compare two .off files and determine if meshes are equal.
 *         Use an exact comparison.
 *
 * \param  filenameA  Filename of first .off file
 * \param  filenameB  Filename of second .off file
 * \param  verbose    Display meshes
 *
 * \result true when meshes are equal, false otherwise
 *
 * \note   Limitations. The method is based on replacing the vertex index
 *         by the vertex coordinates in faces descriptions. The comparison
 *         may fail if several vertices have the exact same coordinates.
 *         For example, if vertices #3 and #4 have the same coordinates,
 *         face #1-#2-#3 can be confused with face #1-#2-#4. Then the
 *         comparison may report equality while the meshes are topologically
 *         different.
 */
bool
are_meshes_equal(std::string filename_a, std::string filename_b, bool verbose);


/*
 * \brief  Compare two .off files and determine if meshes are equal.
 *         Use a "close enough" comparison using an absolute or relative
 *         threshold.
 *
 * \param  filenameA  Filename of first .off file
 * \param  filenameB  Filename of second .off file
 * \param  verbose    Display meshes
 * \param  threshold  relative or absolute error tolerated in float comparison
 *                    for both geometry and attributes ; if the threshold
 *                    is zero, an exact comparison occurs
 * \param  relative_threshold  if true, the threshold is considered as
 *                             a relative value, else it is considered as
 *                             an absolute value
 *
 *
 * \result true when meshes are equal, false otherwise
 *
 * \note   Limitations. The method is based on replacing the vertex index
 *         by the vertex coordinates in faces descriptions. The comparison
 *         may fail if several vertices have the exact same coordinates.
 *         For example, if vertices #3 and #4 have the same coordinates,
 *         face #1-#2-#3 can be confused with face #1-#2-#4. Then the
 *         comparison may report equality while the meshes are topologically
 *         different.
 */
bool
are_meshes_equal(std::string filename_a,
                 std::string filename_b,
                 bool verbose,
                 FloatT threshold,
                 bool relative_threshold);


/*
 * \brief  Compare two .off files and determine if meshes are equal.
 *         Use a "close enough" comparison using an absolute or relative
 *         threshold.
 *
 * \param  filenameA  Filename of first .off file
 * \param  filenameB  Filename of second .off file
 * \param  verbose    Display meshes
 * \param  geom_threshold  relative or absolute error tolerated in float
 *         comparison for geometry ; if the threshold is zero, an exact
 *         comparison occurs
 * \param  attr_threshold  relative or absolute error tolerated in float
 *         comparison for attributes ; if the threshold is zero, an exact
 *         comparison occurs
 * \param  relative_threshold  if true, the threshold is considered as a
 *         relative value, else it is considered as an absolute value
 *
 * \result true when meshes are equal, false otherwise
 *
 * \note   Limitations. The method is based on replacing the vertex index
 *         by the vertex coordinates in faces descriptions. The comparison
 *         may fail if several vertices have the exact same coordinates.
 *         For example, if vertices #3 and #4 have the same coordinates,
 *         face #1-#2-#3 can be confused with face #1-#2-#4. Then the
 *         comparison may report equality while the meshes are topologically
 *         different.
 */
bool
are_meshes_equal(std::string filename_a,
                 std::string filename_b,
                 bool verbose,
                 FloatT geom_threshold,
                 FloatT attr_threshold,
                 bool relative_thresholds);


#include "utils_are_meshes_identical.inl"
