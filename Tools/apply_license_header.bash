# Lowbrow utility script that applies the ownership/license
# header to the sourcefiles (not already bearing one).
# IMPORTANT: this script should be invoked from the repository root


#### Build the list of all the files that require to be tagged:
rm -f all_files_tmp.txt

# We don't have to look for .hxx sources with a find command because there is 
# only single one:
echo 'FEVV/Tools/Container/Helpers.hxx'                    >> all_files_tmp.txt

# Neither do we have to "find" sources with a .cxx extension, because the
# only occurences are encountered in the External sub-directory that we
# don't have to tag because they do not belong to MEPP2.
# Ditto for sources with a .cc extension (that are yet not located in
# External as they should but in `Testing/Draco/`

find Visualization Base FEVV Testing Examples -name \*.h   >> all_files_tmp.txt
find Visualization Base FEVV Testing Examples -name \*.hpp >> all_files_tmp.txt
find Visualization Base FEVV Testing Examples -name \*.cpp >> all_files_tmp.txt
find Visualization Base FEVV Testing Examples -name \*.inl >> all_files_tmp.txt
# The arithmetic_codec[] files should really be in the External sub-directory
# and do not belong to FEVV: we thus weed them out the initial list
grep -v FEVV/Filters/Generic/Manifold/Compression_Valence/arithmetic_codec \
     all_files_tmp.txt > all_files_tmp2.txt
cat all_files_tmp2.txt | xargs grep --files-without-match " Copyright (c) 2012-2019 University of Lyon and CNRS (France)" > all_files.txt
echo -n 'Total number of all files to stamp: '
wc -l < all_files.txt

#### Build the sub-list of sources that should bear a gpl header:
rm -f to_tag_gpl.txt
grep -i polyhedron      all_files.txt       >> to_tag_gpl.txt
grep -i surface         all_files.txt       >> to_tag_gpl.txt
grep -i Visualization   all_files.txt       >> to_tag_gpl.txt
# The following files include all their data structure specific 
# specializations and thus get contaminated by the GPL ones:
echo 'FEVV/DataStructures/DataStructures.h' >> to_tag_gpl.txt
echo 'FEVV/Wrappings/Wrappings.h'           >> to_tag_gpl.txt

echo -n 'Number of gpl files: '
wc -l < to_tag_gpl.txt

#### The sub-list of sources that should bear an LGPL header
# is what remains once the gpl onces are set aside
rm -f to_tag_LGPL.txt
cat to_tag_gpl.txt all_files.txt | sort| uniq -u  >> to_tag_LGPL.txt
echo -n 'Number of LGPL files: '
wc -l < to_tag_LGPL.txt
echo ' '

######## Apply the license header:
echo "// Copyright (c) 2012-2019 University of Lyon and CNRS (France)."            >  header-LGPL.txt
echo "// All rights reserved."                                                     >> header-LGPL.txt
echo "//"                                                                          >> header-LGPL.txt
echo "// This file is part of MEPP2; you can redistribute it and/or modify"        >> header-LGPL.txt
echo "// it under the terms of the GNU Lesser General Public License as"           >> header-LGPL.txt
echo "// published by the Free Software Foundation; either version 3 of"           >> header-LGPL.txt
echo "// the License, or (at your option) any later version."                      >> header-LGPL.txt
echo "//"                                                                          >> header-LGPL.txt
echo "// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE"  >> header-LGPL.txt
echo "// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.">> header-LGPL.txt
cat to_tag_LGPL.txt | while read i; do cat header-LGPL.txt $i > $i.tmp; mv -f $i.tmp $i; done

echo -n 'Number of modified LGPL files: '
git st | grep modified | wc -l

echo "// Copyright (c) 2012-2019 University of Lyon and CNRS (France)."            >  header-gpl.txt
echo "// All rights reserved."                                                     >> header-gpl.txt
echo "//"                                                                          >> header-gpl.txt
echo "// This file is part of MEPP2; you can redistribute it and/or modify"        >> header-gpl.txt
echo "// it under the terms of the GNU General Public License as published "       >> header-gpl.txt
echo "// by the Free Software Foundation; either version 3 of the License, "       >> header-gpl.txt
echo "// or (at your option) any later version."                                   >> header-gpl.txt
echo "//"                                                                          >> header-gpl.txt
echo "// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE"  >> header-gpl.txt
echo "// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.">> header-gpl.txt
cat to_tag_gpl.txt | while read i; do cat header-gpl.txt $i > $i.tmp; mv -f $i.tmp $i; done

echo -n 'Total number of modified files: '
git st | grep modified | wc -l
