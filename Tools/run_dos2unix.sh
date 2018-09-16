#!/bin/bash
#
# Set end of lines to UNIX format (\n).
#

TMP_FILE="/tmp/run_clang-format.tmp.$(date +%F_%T)"

# build source files list
find  Base Examples FEVV Private Testing Visualization  -type f | grep -v osgQt | grep -e '\.cpp$' -e '\.cxx$' -e '\.h$' -e '\.hpp$' -e '\.hxx$' -e '\.inl$' >"$TMP_FILE"
find  . -name CMakeLists.txt  >>"$TMP_FILE"
find  Cmake -name "*.cmake"   >>"$TMP_FILE"

# fix end of lines
while read filename
do
	echo "$filename"
	dos2unix  "$filename"
done <"$TMP_FILE"

echo "$(wc -l <"$TMP_FILE") file(s) processed"

