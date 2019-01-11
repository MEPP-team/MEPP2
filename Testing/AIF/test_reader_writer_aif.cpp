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
#include "FEVV/DataStructures/AIF/AIFMesh.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshReader.hpp"
#include "FEVV/DataStructures/AIF/AIFMeshWriter.hpp"

#include "FEVV/Tools/IO/FileUtilities.hpp"
#include "Testing/Utils/utils_are_meshes_identical.hpp"
#include "Testing/Utils/utils_identical_text_based_files.hpp"

#include <cassert>

using namespace FEVV::DataStructures::AIF;
typedef AIFMeshReader reader_type;
typedef AIFMeshWriter writer_type;

int
main(int narg, char **argv)
{
  std::string input_file_path;
  std::string output_file_path("test_AIFWriter_output.obj");
  std::string reference_file_path;

  if(narg == 2)
  {
    input_file_path = argv[1];
  }
  else if(narg == 3)
  {
    input_file_path = argv[1];
    output_file_path = argv[2];
  }
  else if(narg == 4)
  {
    input_file_path = argv[1];
    output_file_path = argv[2];
    reference_file_path = argv[3];
  }
  else
  {
    std::cout << "Usage: " << argv[0]
              << "  mesh_file  [output_file [reference_file]]" << std::endl;
    std::cout << "Example: " << argv[0] << "  airplane.off" << std::endl;
    std::cout << "Example: " << argv[0] << "  airplane.off airplane.out.obj"
              << std::endl;
    std::cout << "Example: " << argv[0]
              << "  airplane.off airplane.out.obj airplane.ref.obj"
              << std::endl;
    return EXIT_FAILURE;
  }

  //--------------------------------------------------

  reader_type my_reader;
  writer_type my_writer;
  reader_type::ptr_output my_mesh;

  try
  {
    my_mesh = my_reader.read(input_file_path);

    // myMesh->Print();
  }
  catch(const std::invalid_argument &e)
  {
    std::cerr << "Invalid Argument Exception catch while reading "
              << input_file_path << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(const std::ios_base::failure &e)
  {
    std::cerr << "IOS Failure Exception catch while reading " << input_file_path
              << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(const std::length_error &le)
  {
    std::cerr << "[AIF] Exception caught while reading input file "
              << input_file_path << ": " << le.what() << std::endl;
    BOOST_ASSERT_MSG(false, "[AIF] Exception caught while reading input file.");
  }

  try
  {
    my_writer.write(my_mesh, output_file_path);
  }
  catch(const std::invalid_argument &e)
  {
    std::cerr << "Invalid Argument Exception catch while writing "
              << output_file_path << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(const std::ios_base::failure &e)
  {
    std::cerr << "IOS Failure Exception catch while writing "
              << output_file_path << " :" << std::endl
              << e.what() << std::endl;
    return EXIT_FAILURE;
  }

  if(!reference_file_path.empty())
  {
    std::cout << "Comparing output file '" << output_file_path
              << "' with reference file '" << reference_file_path << "'..."
              << std::endl;

    if(FEVV::FileUtils::has_extension(output_file_path, ".off") ||
       FEVV::FileUtils::has_extension(output_file_path, ".coff"))
    {
      // use OFF file comparator
      if(!are_meshes_equal(output_file_path, reference_file_path, false))
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
