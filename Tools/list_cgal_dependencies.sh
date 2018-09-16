#!/bin/bash
#
# List CGAL-headers dependencies of a source code file.
# Depends on g++.
#

SRCCODE="$1"

if [ -z "$SRCCODE" ]
then
	echo "Usage: $0  file_to_analyze  [-Ddefine1 -Ddefine2 -...] [-I headers_dir_1 -I headers_dir_2 ...]"
	echo "Example: $0  ~/CGAL-4.8/include/CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h  -DCGAL_USE_GMP -DCGAL_USE_MPFR -DNOMINMAX -D_USE_MATH_DEFINES -DNDEBUG -std=c++11  -I ~/CGAL-4.8/include -I ~/OpenMesh-5.1/include]"
	exit 1
fi

# get provided compiler flags and include dirs
shift
OPTIONS="$@"

# build output files names
OUTRAW="$(basename $SRCCODE).depencies.raw"
OUTCGAL="$(basename $SRCCODE).depencies.cgal"
OUTCGALH="$(basename $SRCCODE).depencies.cgal.hierarchical"

# find headers dependencies
echo
(set -x ; g++ -M $OPTIONS "$SRCCODE" >"$OUTRAW")
if [ "$?" = "1" ]
then
	echo "Failed to find some headers. Need more headers directories. Rerun with more -I options."
	exit 1
fi

# clean the raw output:
#  - remove ' \' at end of line
#  - remove ' ' at beginning of line
#  - split lines with multiple file names
#  - keep only lines containing 'CGAL'
sed -e 's/ *\\//' -e 's/^ //' -e 's/ /\n/g' <"$OUTRAW" | grep CGAL >"$OUTCGAL"

# find hierarchical headers dependencies
g++ -H $OPTIONS "$SRCCODE" 2>"$OUTCGALH.tmp" >/dev/null 
grep '^\.' "$OUTCGALH.tmp" | grep CGAL | uniq >"$OUTCGALH"

# display result
echo
echo "- source code analyzed: $SRCCODE"
echo "- raw dependencies output: $OUTRAW"
echo "- dependencies to CGAL: $OUTCGAL"
echo "- hierachical dependencies to CGAL: $OUTCGALH"

echo "- files suspected not to be LGPL:"
while read name ; do grep -L 'GNU Lesser General Public License'  "$name" ; done <"$OUTCGAL"
echo

