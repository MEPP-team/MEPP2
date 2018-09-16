#ifndef MATERIAL_H
#define MATERIAL_H

#include <string>


namespace FEVV {
namespace Types {


enum class MaterialType { MATERIAL_TYPE_STANDARD, MATERIAL_TYPE_PBR };

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
  std::string
      ambient_texture_filename{}; // Used for ambient occlusion map if PBR
  std::string diffuse_texture_filename{}; // Used for albedo map if PBR
  std::string specular_texture_filename{};
  std::string emissive_texture_filename{};
  std::string transparency_texture_filename{};
  std::string normal_map_filename{}; // Used for normal map if PBR, and bump map
                                     // if needed

  // PBR maps fields
  // albedo_map_filename = diffuse_texture_filename
  std::string metallic_map_filename{};
  std::string roughness_map_filename{};
  // ambient_occlusion_map_filename = ambient_texture_filename

  bool has_normal_map = false;
};


} // namespace Types
} // namespace FEVV

#endif // MATERIAL_H
