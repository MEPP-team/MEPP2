#!/bin/bash
#
# Import CGAL headers, preserving the directory
# hierarchy under CGAL-x.x/include.
#

if [ "$#" != "2" ]
then
	echo "Usage: $0  dependencies_file  destination_dir"
	echo "Example: $0  all.dependencies.cgal ../FEVV"
	exit 1
fi

INPUTFILE="$1"

# retrieve destination dir absolute path
DESTDIR="$(readlink -f "$2")"

BASEDIR=""

# copy headers
while read fname
do
	if [ -z "$BASEDIR" ]
	then
		# retrieve the base dir the first time
		BASEDIR="$(echo "$fname" | sed -r 's,(.*/include)/CGAL/.*,\1,')"
		cd "$BASEDIR"
	fi

	relativefname="$(echo $fname | sed "s,$BASEDIR/,,")" 
	cp --parent "$relativefname" "$DESTDIR"
done < "$INPUTFILE"

