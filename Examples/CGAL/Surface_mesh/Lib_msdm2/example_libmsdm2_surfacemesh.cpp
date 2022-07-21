
#include <iostream>
#include <string>
#include <fstream>

#include "msdm2_surfacemesh.h"

int main(int argc, char **argv)
{
  // check arguments
  if(argc != 3)
  {
    std::cout << "Compute the MSDM2 value between the degraded mesh and "
                 "the original mesh.\n"
              << std::endl;
    std::cout << "Usage: "
              << argv[0] << "  path/to/original_mesh.off  path/to/degraged_mesh.off"
              << std::endl;
    std::cout << "Note: only OFF files are supported."
              << std::endl;

    return 1;
  }

  // input and output files
  std::string filename_original(argv[1]);
  std::string filename_degraded(argv[2]);

  // read meshes from files
  msdm2::MeshT mesh_original;
  msdm2::MeshT mesh_degraded;

  std::ifstream in1(filename_original);
  in1 >> mesh_original;
  std::cout << "Original mesh file name: " << filename_original << std::endl;
  std::cout << "Original mesh vertices number: " << mesh_original.number_of_vertices() << std::endl;

  std::ifstream in2(filename_degraded);
  in2 >> mesh_degraded;
  std::cout << "Degraded mesh file name: " << filename_degraded << std::endl;
  std::cout << "Degraded mesh vertices number: " << mesh_degraded.number_of_vertices() << std::endl;

  // call msdm2 filter
  int nb_levels = 3;
  double msdm2_value;
  msdm2::Msdm2MapT msdm2_map;
  msdm2::msdm2_surfacemesh(mesh_degraded,
                           mesh_original,
                           nb_levels,
                           msdm2_value,    /* output */
                           msdm2_map);     /* output */

  std::cout << "Calculated MSDM2 : " << msdm2_value << std::endl;

  return 0;
}


