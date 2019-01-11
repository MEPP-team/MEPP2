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
#include "FEVV/Filters/Generic/generic_reader.hpp"
#include "FEVV/Filters/Generic/generic_writer.hpp"
#include "FEVV/Filters/Generic/Manifold/Compression_Valence/compression_valence.h"

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
test_compression_valence(int argc, const char **argv)
{
  // parse arguments

  if(argc < 2)
  {
    std::cout << "Load a mesh then apply the Compression Valence filter."
              << std::endl;
    std::cout << "Usage: " << argv[0]
              << "  input_mesh_file  [simplified_mesh_file  "
                 "[reference_mesh_file]]  [-withCompr]  [-withQA]  [-maxV nnn] "
                 " [-Qbits nnn]  [-tolerance relative_tolerance]"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.off  airplane.out.obj  airplane.ref.obj  "
                 "-withCompr  -maxV 500  -tolerance 1e-5"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.off  airplane.out.obj  -withCompr  -maxV 500  "
                 "-tolerance 1e-5"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.off  -withCompr  -maxV 500  -tolerance 1e-5"
              << std::endl;
    return EXIT_FAILURE;
  }

  std::string input_file_path;
  std::string simplified_mesh_file_path;
  std::string reference_file_path;
  bool with_adaptative_quantization = false;
  bool with_compression = false;
  int quantiz_bits = 10;
  int max_vertices = 100;
  float relative_tolerance = 0.0f;

  input_file_path = argv[1];
  int arg_read = 1;
  if(argc >= 3 && argv[2][0] != '-')
  {
    simplified_mesh_file_path = argv[2];
    arg_read = 2;
    if(argc >= 4 && argv[3][0] != '-')
    {
      reference_file_path = argv[3];
      arg_read = 3;
    }
  }

  for(int i = arg_read + 1; i < argc; i++)
  {
    std::string argvi(argv[i]);

    if(argvi == "-withQA")
    {
      with_adaptative_quantization = true;
    }
    else if(argvi == "-withCompr")
    {
      with_compression = true;
    }
    else if(argvi == "-Qbits")
    {
      quantiz_bits = std::atoi(argv[i + 1]);
      i++;
    }
    else if(argvi == "-maxV")
    {
      max_vertices = std::atoi(argv[i + 1]);
      i++;
    }
    else if(argvi == "-tolerance")
    {
      relative_tolerance = (float)std::atof(argv[i + 1]);
      i++;
    }
    else
    {
      std::cout << "unknown option '" << argvi << "'. Aborting." << std::endl;
      return EXIT_FAILURE;
    }
  }

  std::string p3d_file_path;
  if(simplified_mesh_file_path.empty())
    p3d_file_path =
        input_file_path.substr(0, input_file_path.find_last_of(".")) + ".p3d";
  else
    p3d_file_path = simplified_mesh_file_path.substr(
                        0, simplified_mesh_file_path.find_last_of(".")) +
                    ".p3d";

  // display parameters summary
  std::cout << "\nParameters summary:" << std::endl;
  std::cout << " - input file: " << input_file_path << std::endl;
  std::cout << " - simplified mesh file: " << simplified_mesh_file_path
            << std::endl;
  std::cout << " - p3d file: " << p3d_file_path << std::endl;
  std::cout << " - ref. file: " << reference_file_path << std::endl;
  std::cout << " - compression: " << std::boolalpha << with_compression
            << std::endl;
  std::cout << " - adaptative quantization: " << std::boolalpha
            << with_adaptative_quantization << std::endl;
  std::cout << " - max vertices: " << max_vertices << std::endl;
  std::cout << " - quantization bits: " << quantiz_bits << std::endl;

  // display memory usage
#ifdef UNIX
  pid_t pid = getpid();
  uint vmpeak_kb;
  uint vmsize_kb;

  get_vmpeak_vmsize(pid, vmpeak_kb, vmsize_kb);
  std::cout << "before compression   vmpeak=" << vmpeak_kb
            << " kB   vmsize=" << vmsize_kb << " kB" << std::endl;
#endif

  //--------------------------------------------------

  // read mesh from file
  MeshT m;
  FEVV::PMapsContainer pmaps_bag;
  FEVV::Filters::read_mesh(input_file_path, m, pmaps_bag);

  // retrieve geometry property map
  auto pm = get(boost::vertex_point, m);

  // retrieve vertex color property map
  using VertexColorMap =
      typename FEVV::PMap_traits< FEVV::vertex_color_t, MeshT >::pmap_type;
  VertexColorMap v_cm;
  VertexColorMap *v_cm_ptr = nullptr;
  if(has_map(pmaps_bag, FEVV::vertex_color))
  {
    std::cout << "using vertex color property map" << std::endl;
    v_cm = get_property_map(FEVV::vertex_color, m, pmaps_bag);
    v_cm_ptr = &v_cm;
  }

  // apply Compression Valence filter
  FEVV::Filters::compression_valence(m,
                                     &pm,
                                     v_cm_ptr,
                                     input_file_path,
                                     p3d_file_path,
                                     with_compression,
                                     with_adaptative_quantization,
                                     max_vertices,
                                     quantiz_bits);

  // display memory usage
#ifdef UNIX
  uint vmpeak2_kb;
  get_vmpeak_vmsize(pid, vmpeak2_kb, vmsize_kb);
  std::cout << "after compression   vmpeak=" << vmpeak2_kb
            << " kB   vmsize=" << vmsize_kb
            << " kB   extra.vmpeak=" << vmpeak2_kb - vmpeak_kb << " kB"
            << std::endl;
#endif

  // save the mesh
  if(!simplified_mesh_file_path.empty())
    FEVV::Filters::write_mesh(simplified_mesh_file_path, m, pmaps_bag);

  // check output file
  if(!reference_file_path.empty())
  {
    std::cout << "Comparing output file" << std::endl;
    std::cout << "  '" << simplified_mesh_file_path << "'" << std::endl;
    std::cout << "with reference file" << std::endl;
    std::cout << "  '" << reference_file_path << "'" << std::endl;
    std::cout << "..." << std::endl;

    if(FEVV::FileUtils::has_extension(simplified_mesh_file_path, ".off") ||
       FEVV::FileUtils::has_extension(simplified_mesh_file_path, ".coff"))
    {
      // use OFF file comparator
      if(!are_meshes_equal(simplified_mesh_file_path,
                           reference_file_path,
                           false,
                           relative_tolerance,
                           true))
      {
        std::cout << "Files are different!" << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      // use text file comparator
      if(!identical_text_based_files(simplified_mesh_file_path,
                                     reference_file_path))
      {
        std::cout << "Files are different!" << std::endl;
        return EXIT_FAILURE;
      }
    }

    std::cout << "Files are identical." << std::endl;
  }

  return 0;
}
