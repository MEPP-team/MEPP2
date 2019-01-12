// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "FEVV/Filters/Generic/Manifold/Compression_Valence/decompression_valence.h"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <iostream>
#include <string>
#include <cstdlib>


#ifdef UNIX
#include <unistd.h> // for getpid()
#include <fstream>

void
get_vmpeak_vmsize(pid_t pid, unsigned int &vmpeak_kb, unsigned int &vmsize_kb)
{
  std::string filename = "/proc/" + std::to_string(pid) + "/status";
  std::ifstream file(filename);
  if(!file.is_open())
  {
    std::cout << "Failed to open file '" << filename << "'" << std::endl;
    return;
  }

  std::string word;
  while(true)
  {
    file >> word;
    // DBG std::cout << "word=" << word << std::endl;

    if(file.eof())
      break;

    if(word == "VmPeak:")
      file >> vmpeak_kb;
    else if(word == "VmSize:")
      file >> vmsize_kb;
  }

  file.close();
}
#endif // UNIX


template< typename MeshT >
int
test_decompression_valence(int argc, const char **argv)
{
  // parse arguments

  if(argc < 2 || argc > 5)
  {
    std::cout << "Load a .p3d file generated by the Compression Valence "
                 "filter, then decompress it."
              << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  p3d_file  [output_mesh_file  [reference_mesh_file  "
                 "[relative_tolerance]]]"
              << std::endl;
    std::cout << "Example: " << argv[0] << "  airplane.p3d" << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.p3d  airplane_uncompressed.off" << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.p3d  airplane_uncompressed.off  airplane.ref.obj"
              << std::endl;
    std::cout
        << "Example: " << argv[0]
        << "  airplane.p3d  airplane_uncompressed.off  airplane.ref.obj  1e-04"
        << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_path;
  std::string output_file_path;
  std::string reference_file_path;
  float relative_tolerance = 0.0f;

  input_file_path = argv[1];
  if(argc >= 3)
    output_file_path = argv[2];
  if(argc >= 4)
    reference_file_path = argv[3];
  if(argc >= 5)
    relative_tolerance = (float)std::atof(argv[4]);

  // display parameters summary
  std::cout << "\nParameters summary:" << std::endl;
  std::cout << " - input file: " << input_file_path << std::endl;
  std::cout << " - output file: " << output_file_path << std::endl;
  std::cout << " - ref. file: " << reference_file_path << std::endl;
  std::cout << " - rel. tolerance: " << relative_tolerance << std::endl;

  // display memory usage before loading and processing mesh
#ifdef UNIX
  pid_t pid = getpid();
  uint vmpeak_kb;
  uint vmsize_kb;

  get_vmpeak_vmsize(pid, vmpeak_kb, vmsize_kb);
  std::cout << "before loading and decompression   vmpeak=" << vmpeak_kb
            << " kB   vmsize=" << vmsize_kb << " kB" << std::endl;
#endif

  //--------------------------------------------------

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;

  // retrieve geometry property map
  auto pm = get(boost::vertex_point, m);

  // create vertex color property map
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type;
  VertexColorMap v_cm;
  put_property_map(FEVV::vertex_color, m, pmaps_bag, v_cm);

  // apply Compression Valence filter
  FEVV::Filters::decompression_valence(m, &pm, &v_cm, input_file_path);

  // display memory usage after decompression
#ifdef UNIX
  uint vmpeak2_kb;
  get_vmpeak_vmsize(pid, vmpeak2_kb, vmsize_kb);
  std::cout << "after decompression   vmpeak=" << vmpeak2_kb
            << " kB   vmsize=" << vmsize_kb
            << " kB   extra.vmpeak=" << vmpeak2_kb - vmpeak_kb << " kB"
            << std::endl;
#endif

  // save the uncompressed mesh
  if(!output_file_path.empty())
    FEVV::Filters::write_mesh(output_file_path, m, pmaps_bag);

  // check output file
  if(!reference_file_path.empty())
  {
    std::cout << "Comparing output file" << std::endl;
    std::cout << "  '" << output_file_path << "'" << std::endl;
    std::cout << "with reference file" << std::endl;
    std::cout << "  '" << reference_file_path << "'" << std::endl;
    std::cout << "..." << std::endl;

    if(FEVV::FileUtils::has_extension(output_file_path, ".off") ||
       FEVV::FileUtils::has_extension(output_file_path, ".coff"))
    {
      // use OFF file comparator
      if(!are_meshes_equal(output_file_path,
                           reference_file_path,
                           false,
                           relative_tolerance,
                           true /*relative_threshold*/))
      {
        std::cout << "Files are different!" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      // use text file comparator
      if(!identical_text_based_files(output_file_path, reference_file_path))
      {
        std::cout << "Files are different!" << std::endl;
        return EXIT_FAILURE;
      }
    }

    std::cout << "Files are identical." << std::endl;
  }

  return 0;
}
