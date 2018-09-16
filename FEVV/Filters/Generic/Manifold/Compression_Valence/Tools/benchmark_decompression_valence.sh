#!/bin/bash
#
#
# Measure time and memory footprint for decompression with meshes of various size
# Must been run in Mepp2-elo/build.xxx dir
#

P3DFILES="\
$HOME/temp/MEPP2-elo/hidden/TestObjects/bench-decompr/bimba-mepp1.p3d \
$HOME/temp/MEPP2-elo/hidden/TestObjects/bench-decompr/monkey-mepp1.p3d \
$HOME/temp/MEPP2-elo/hidden/TestObjects/bench-decompr/neptune-mepp1.p3d \
$HOME/temp/MEPP2-elo/hidden/TestObjects/bench-decompr/happy-mepp1.p3d \
"

PRGS="\
Testing/CGAL/test_decompression_valence_polyhedron \
Testing/CGAL/test_decompression_valence_surfacemesh \
Testing/CGAL/test_decompression_valence_lcc \
Testing/OpenMesh/test_decompression_valence_openmesh \
"

for prg in $PRGS
do
	for p3dfile in $P3DFILES
	do
		echo "------------------------------------------------------------------"

		datastructure="$(echo "$prg" | sed 's/.*_//')"
		#rm  output_mesh="$(basename $p3dfile)-decomp-mepp2-$datastructure.off"

		for i in 1 2 3
		do
			echo
			echo "------------------------------------------------------------------"
			echo "decompressing '$p3dfile'"
			echo "with '$prg'"

			#rm  commande="$prg  $p3dfile  /tmp/$output_mesh  bidon"
			commande="$prg  $p3dfile"
			echo "'$commande'"

			eval $commande
		done
	done
done


