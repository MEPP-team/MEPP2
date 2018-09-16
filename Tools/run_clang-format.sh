#!/bin/bash
#
# Reformat source code using clang-format.
# Format configuration is in .clang-format file.
#

TMP_FILE="/tmp/run_clang-format.tmp.$(date +%F_%T)"

# build source files list
find  Base Examples FEVV Private Testing Visualization  -type f | grep -v osgQt | grep -e '\.cpp$' -e '\.cxx$' -e '\.h$' -e '\.hpp$' -e '\.hxx$' -e '\.inl$' >"$TMP_FILE"

# reformat source files
while read filename
do
	echo "$filename"
	clang-format  -i  -style=file  "$filename"
done <"$TMP_FILE"

echo "$(wc -l <"$TMP_FILE") file(s) processed"
