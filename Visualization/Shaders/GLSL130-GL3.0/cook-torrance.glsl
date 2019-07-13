#version 130

#define MAX_LIGHT_COUNT 10
#define PI 3.1415926535897932384626433832795

struct Light {
  vec4 position;
  vec3 direction;
  vec4 color;
  float energy;
  float angle;
};

struct Material {
  vec3 baseColor;
  float metallicFactor;
  float roughnessFactor;
  vec3 emissive;

  sampler2D albedoMap;
  sampler2D normalMap;
  sampler2D metallicMap;
  sampler2D roughnessMap;
  sampler2D emissiveMap;
  sampler2D ambientOcclusionMap;
};

in vec3 fragMeshInfo_vertPosition;
in vec2 fragMeshInfo_vertTexcoords;
in mat3 fragMeshInfo_vertTBNMatrix;
flat in vec3 fragMeshInfo_vertFlatNormal;

uniform bool uniUseSmoothShading = true;
uniform bool uniUseNormalMapping = false;
uniform bool uniUseLighting;
uniform bool uniUseIndirectLighting = false;

uniform uint uniLightCount;
uniform Light uniLights[MAX_LIGHT_COUNT];

uniform mat4 uniViewMatrix;
uniform mat4 uniInvViewMatrix;
uniform vec3 uniCameraPos;
uniform vec3 uniCameraTarget;

uniform Material uniMaterial;

out vec4 fragColor;
out vec3 bufferNormal;

// Normal Distribution Function: Trowbridge-Reitz GGX
float computeNormalDistrib(vec3 normal, vec3 halfVec, float roughness) {
  float sqrRough  = roughness * roughness;
  float frthRough = sqrRough * sqrRough;

  float halfVecAngle    = max(dot(halfVec, normal), 0.0);
  float sqrHalfVecAngle = halfVecAngle * halfVecAngle;

  float divider = (sqrHalfVecAngle * (frthRough - 1.0) + 1.0);
  divider       = PI * divider * divider;

  return frthRough / max(divider, 0.001);
}

// Fresnel: Shlick
vec3 computeFresnelShlick(float cosTheta, vec3 baseReflectivity) {
  // Seemingly coefficient's optimized version:
  return baseReflectivity + (1.0 - baseReflectivity) * pow(2.0, (-5.55473 * cosTheta - 6.98316) * cosTheta);
  //return baseReflectivity + (1.0 - baseReflectivity) * pow(1.0 - cosTheta, 5.0);
}

// Shlick-Beckmann for Geometry part
float computeShlickGGX(float viewAngle, float roughness) {
  float incrRough   = (roughness + 1.0);
  float roughFactor = (incrRough * incrRough) / 8.0;

  float denom = viewAngle * (1.0 - roughFactor) + roughFactor;

  return viewAngle / denom;
}

// Geometry: Smith's Shlick GGX
float computeGeometry(vec3 normal, vec3 viewDir, vec3 lightDir, float roughness) {
  float viewAngle  = max(dot(viewDir, normal), 0.0);
  float lightAngle = max(dot(lightDir, normal), 0.0);

  float ggx1 = computeShlickGGX(viewAngle, roughness);
  float ggx2 = computeShlickGGX(lightAngle, roughness);

  return ggx1 * ggx2;
}

vec3 computeLightRadiance(vec3 viewDir, vec3 lightDir, vec3 normal,     // Direction vectors
                          vec3 albedo, float metallic, float roughness, // Material attributes
                          vec3 radiance, vec3 baseReflectivity) {       // Precomputed parameters
  vec3 halfDir = normalize(viewDir + lightDir);

  // Normal distrib (D)
  float normalDistrib = computeNormalDistrib(normal, halfDir, roughness);

  // Fresnel (F)
  vec3 fresnel = computeFresnelShlick(max(dot(halfDir, viewDir), 0.0), baseReflectivity);

  // Geometry (G)
  float geometry = computeGeometry(normal, viewDir, lightDir, roughness);

  vec3 DFG      = normalDistrib * fresnel * geometry;
  float divider = 4.0 * max(dot(viewDir, normal), 0.0) * max(dot(lightDir, normal), 0.0);
  vec3 specular = DFG / max(divider, 0.001);

  vec3 diffuse = vec3(1.0) - fresnel;
  diffuse     *= 1.0 - metallic;

  return (diffuse * albedo / PI + specular) * radiance * max(dot(lightDir, normal), 0.0);
}

void main() {
  // Gamma correction for albedo (sRGB presumed)
  vec3 albedo     = pow(texture(uniMaterial.albedoMap, fragMeshInfo_vertTexcoords).rgb, vec3(2.2)) * uniMaterial.baseColor;
  float metallic  = texture(uniMaterial.metallicMap, fragMeshInfo_vertTexcoords).r * uniMaterial.metallicFactor;
  float roughness = texture(uniMaterial.roughnessMap, fragMeshInfo_vertTexcoords).r * uniMaterial.roughnessFactor;
  float ambOcc    = texture(uniMaterial.ambientOcclusionMap, fragMeshInfo_vertTexcoords).r;

  vec3 normal;
  if (uniUseSmoothShading) {
    if (uniUseNormalMapping) {
      // Normal mapping version
      normal = texture(uniMaterial.normalMap, fragMeshInfo_vertTexcoords).rgb;
      normal = normalize(normal * 2.0 - 1.0);
      normal = normalize(fragMeshInfo_vertTBNMatrix * normal);
    } else {
      // Standard version
      normal = normalize(fragMeshInfo_vertTBNMatrix[2]);
    }
  } else {
    normal = normalize(fragMeshInfo_vertFlatNormal);
  }

  vec3 fullViewDir = uniCameraPos - fragMeshInfo_vertPosition;
  vec3 viewDir     = normalize(fullViewDir);

  mat3 invViewMatrix = mat3(uniInvViewMatrix);

  vec3 color = albedo;

  if (uniUseLighting) {
    vec3 lightRadiance = vec3(0.0);

    // Base Fresnel
    vec3 baseReflectivity = mix(vec3(0.04), albedo, metallic);

    // Creating light sources whether we pick direct or indirect lighting
    {
      // Default energy for camera light, set the dot() falloff to vary intensity with distance
      float attenuation = 1.0/* / dot(fullViewDir, fullViewDir)*/;
      vec3 radiance = vec3(attenuation);

      if (!uniUseIndirectLighting) {
        // Taking camera as light source
        lightRadiance += computeLightRadiance(viewDir, viewDir, normal,
                                              albedo, metallic, roughness,
                                              radiance, baseReflectivity);
      } else {
        float camDist = length(fullViewDir);

        // Create other light sources based on mesh's position
        // Key light
        vec3 keyLightPos = invViewMatrix * normalize(vec3(-1.0, 1.0, 0.0));
        vec3 keyLightDir = normalize(keyLightPos);
        lightRadiance   += computeLightRadiance(viewDir, keyLightDir, normal,
                                                albedo, metallic, roughness,
                                                radiance, baseReflectivity);

        // Fill light
        vec3 fillLightPos = invViewMatrix * normalize(vec3(1.0, 0.0, 0.0));
        vec3 fillLightDir = normalize(fillLightPos);
        lightRadiance    += computeLightRadiance(viewDir, fillLightDir, normal,
                                                 albedo, metallic, roughness,
                                                 radiance, baseReflectivity);

        // Back light
        vec3 backLightPos = invViewMatrix * normalize(vec3(0.0, camDist, -camDist * 2.0));
        vec3 backLightDir = normalize(backLightPos);
        lightRadiance    += computeLightRadiance(viewDir, backLightDir, normal,
                                                 albedo, metallic, roughness,
                                                 radiance, baseReflectivity);
      }
    }

    // --- Use light sources

    for (uint lightIndex = 0u; lightIndex < uniLightCount; ++lightIndex) {
      vec3 fullLightDir;
      float attenuation = uniLights[lightIndex].energy;

      if (uniLights[lightIndex].position.w != 0.0) {
        fullLightDir = uniLights[lightIndex].position.xyz - fragMeshInfo_vertPosition;

        float sqrDist = dot(fullLightDir, fullLightDir);
        attenuation  /= sqrDist;
      } else {
        fullLightDir = -uniLights[lightIndex].direction;
      }

      vec3 lightDir = normalize(fullLightDir);
      vec3 radiance = uniLights[lightIndex].color.rgb * attenuation;

      lightRadiance += computeLightRadiance(viewDir, lightDir, normal,
                                            albedo, metallic, roughness,
                                            radiance, baseReflectivity);
    }

    vec3 ambient = vec3(0.03) * albedo * ambOcc;
    color        = ambient + lightRadiance;
  }

  vec3 emissive = texture(uniMaterial.emissiveMap, fragMeshInfo_vertTexcoords).rgb * uniMaterial.emissive;
  color        += emissive;

  // HDR tone mapping
  color = color / (color + vec3(1.0));
  // Gamma correction
  color = pow(color, vec3(1.0 / 2.2));

  fragColor = vec4(color, 1.0);

  // Sending normal to next framebuffer(s), if any
  bufferNormal = normal;
}
