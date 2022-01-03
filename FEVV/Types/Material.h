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
#pragma once

#include <string>
#include <memory>
#include <map>


#ifdef FEVV_USE_JPEG
#define cimg_use_jpeg
#endif
#ifdef FEVV_USE_PNG
#define cimg_use_png
#endif
#ifdef FEVV_USE_TIFF
#define cimg_use_tiff
#endif

#pragma push_macro("_PTHREAD_H")
#undef _PTHREAD_H      // to avoid linking with pthread
#define cimg_display 0 // no CImg display
#include <CImg.h>


namespace FEVV {
namespace Types {


enum class MaterialType { MATERIAL_TYPE_STANDARD, MATERIAL_TYPE_PBR };

#if 1
typedef  cimg_library::CImg< unsigned char >  Image;
#else
// debug
class Image
{
public:
  Image(const char *name): name(name)
  {
    std::cout << "*** Image '" << name << "' created." << std::endl;
  }

  ~Image()
  {
    std::cout << "*** Image '" << name << "' destroyed." << std::endl;
  }

  std::string name;
};
#endif


/**
 * The material datastructure used to handle textures.
 */
struct Material
{
  std::string name{};
  MaterialType type = MaterialType::MATERIAL_TYPE_STANDARD;

  // Standard factors fields
  double ambient_red_component = 1.0;
  double ambient_green_component = 1.0;
  double ambient_blue_component = 1.0;

  double diffuse_red_component = 1.0; // Diffuse used for albedo factor if PBR
  double diffuse_green_component = 1.0;
  double diffuse_blue_component = 1.0;

  double specular_red_component = 0.0;
  double specular_green_component = 0.0;
  double specular_blue_component = 0.0;

  double emissive_red_component = 0.0;
  double emissive_green_component = 0.0;
  double emissive_blue_component = 0.0;

  double transparency = 1.0;

  // PBR factors fields
  // albedo*Component = diffuse_*Component
  double metallic_factor = 1.0;
  double roughness_factor = 1.0;

  // Standard textures fields
  std::string ambient_texture_filename{}; // Used for ambient occlusion map if
                                          // PBR
  std::string diffuse_texture_filename{}; // Used for albedo map if PBR
  std::string specular_texture_filename{};
  std::string emissive_texture_filename{};
  std::string transparency_texture_filename{};
  std::string normal_map_filename{};      // Used for normal map if PBR, and
                                          // bump map if needed

  // PBR maps fields
  // albedo_map_filename = diffuse_texture_filename
  std::string metallic_map_filename{};
  std::string roughness_map_filename{};
  // ambient_occlusion_map_filename = ambient_texture_filename

  bool has_normal_map = false;

  // storage for texture images
  std::map<std::string, std::shared_ptr< Image > > images;
};


} // namespace Types
} // namespace FEVV

#pragma pop_macro("_PTHREAD_H")
