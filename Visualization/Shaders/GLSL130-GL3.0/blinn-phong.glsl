#version 130

#define MAX_LIGHT_COUNT 10

struct Light {
  vec4 position;
  vec3 direction;
  vec4 color;
  float energy;
  float angle;
};

struct Material {
  vec3 ambient;
  vec3 diffuse;
  vec3 specular;
  vec3 emissive;
  float transparency;

  sampler2D ambientMap;
  sampler2D diffuseMap;
  sampler2D specularMap;
  sampler2D emissiveMap;
  sampler2D transparencyMap;
  sampler2D bumpMap;
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

void main() {
  vec3 normal;
  if (uniUseSmoothShading) {
    if (uniUseNormalMapping) {
      // Normal mapping version
      normal = texture(uniMaterial.bumpMap, fragMeshInfo_vertTexcoords).rgb;
      normal = normalize(normal * 2.0 - 1.0);
      normal = normalize(fragMeshInfo_vertTBNMatrix * normal);
    } else {
      // Standard version
      normal = normalize(fragMeshInfo_vertTBNMatrix[2]);
    }
  } else {
    normal = normalize(fragMeshInfo_vertFlatNormal);
  }

  vec3 baseColor  = texture(uniMaterial.diffuseMap, fragMeshInfo_vertTexcoords).rgb * uniMaterial.diffuse;
  vec3 color      = baseColor;
  vec3 specFactor = texture(uniMaterial.specularMap, fragMeshInfo_vertTexcoords).r * uniMaterial.specular;

  vec3 fullViewDir = uniCameraPos - fragMeshInfo_vertPosition;
  vec3 viewDir     = normalize(fullViewDir);

  mat3 invViewMatrix = mat3(uniInvViewMatrix);

  if (uniUseLighting) {
    vec3 ambient  = baseColor * 0.05;
    vec3 diffuse  = vec3(0.0);
    vec3 specular = vec3(0.0);

    // Creating light sources whether we pick direct or indirect lighting
    {
      if (!uniUseIndirectLighting) {
        // Taking camera as light source
        diffuse  += max(dot(viewDir, normal), 0.0) * baseColor;
        specular += pow(max(dot(viewDir, normal), 0.0), 32.0) * specFactor;
      } else {
        // Create other light sources based on mesh's position
        float camDist = length(fullViewDir);

        // Key light
        vec3 keyLightPos = invViewMatrix * normalize(vec3(-1.0, 1.0, 0.0));
        vec3 keyLightDir = normalize(keyLightPos);
        diffuse         += max(dot(keyLightDir, normal), 0.0) * baseColor;
        specular        += pow(max(dot(normalize(keyLightDir + viewDir), normal), 0.0), 32.0) * specFactor;

        // Fill light
        vec3 fillLightPos = invViewMatrix * normalize(vec3(1.0, 0.0, 0.0));
        vec3 fillLightDir = normalize(fillLightPos);
        diffuse          += max(dot(fillLightDir, normal), 0.0) * baseColor;
        specular         += pow(max(dot(normalize(fillLightDir + viewDir), normal), 0.0), 32.0) * specFactor;

        // Back light
        vec3 backLightPos = invViewMatrix * normalize(vec3(0.0, camDist, -camDist * 2.0));
        vec3 backLightDir = normalize(backLightPos);
        diffuse          += max(dot(backLightDir, normal), 0.0) * baseColor;
        specular         += pow(max(dot(normalize(backLightPos + viewDir), normal), 0.0), 32.0) * specFactor;
      }
    }

    for (uint lightIndex = 0u; lightIndex < uniLightCount; ++lightIndex) {
      // Diffuse
      vec3 fullLightDir;
      float attenuation = uniLights[lightIndex].energy;

      if (uniLights[lightIndex].position.w != 0.0) {
        fullLightDir = normalize(uniLights[lightIndex].position.xyz - fragMeshInfo_vertPosition);

        float sqrDist = dot(fullLightDir, fullLightDir);
        attenuation  /= sqrDist;
      } else {
        fullLightDir = normalize(-uniLights[lightIndex].direction);
      }

      vec3 lightDir = normalize(fullLightDir);

      diffuse += max(dot(lightDir, normal), 0.0) * baseColor * attenuation;

      // Specular
      vec3 halfDir = normalize(lightDir + viewDir);
      specular    += uniLights[lightIndex].color.rgb * pow(max(dot(halfDir, normal), 0.0), 32.0) * specFactor * attenuation;
    }

    // Result
    color = ambient + diffuse + specular;
  }

  vec3 emissive = texture(uniMaterial.emissiveMap, fragMeshInfo_vertTexcoords).rgb * uniMaterial.emissive;
  color        += emissive;

  fragColor = vec4(color, 1.0);

  // Sending normal to next framebuffer(s), if any
  bufferNormal = normal;
}
