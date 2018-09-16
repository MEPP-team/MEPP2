#!/bin/bash
#
#
# Measure time and memory footprint for compression with meshes of various size
# Must been run in Mepp2-elo/build.xxx dir
#

MESHES="\
$HOME/temp/MEPP2-elo/hidden/TestObjects/bimba.obj \
$HOME/temp/MEPP2-elo/hidden/TestObjects/norgb_monkey.obj \
$HOME/temp/MEPP2-elo/hidden/TestObjects/neptune.obj \
$HOME/temp/MEPP2-elo/hidden/TestObjects/happy.obj \
"

PRGS="\
Testing/CGAL/test_compression_valence_polyhedron \
Testing/CGAL/test_compression_valence_surfacemesh \
Testing/CGAL/test_compression_valence_lcc \
Testing/OpenMesh/test_compression_valence_openmesh \
"

for prg in $PRGS
do
	for mesh in $MESHES
	do
		echo "------------------------------------------------------------------"

		datastructure="$(echo "$prg" | sed 's/.*_//')"

		for i in 1 2 3
		do
			echo
			echo "------------------------------------------------------------------"
			echo "compressing '$mesh'"
			echo "with '$prg'"

			#commande="$prg  $mesh  -withCompr  -maxV 100  -Qbits 10"
			commande="$prg  $mesh  -maxV 100  -Qbits 10"
			echo "'$commande'"

			eval $commande
		done
	done
done

