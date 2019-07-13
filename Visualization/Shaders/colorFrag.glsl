#version 330 core

in MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
  vec3 vertNormal;
  flat vec3 vertFlatNormal;
} fragMeshInfo;

uniform bool uniUseSmoothShading    = true;
uniform bool uniUseLighting         = true;
uniform bool uniUseIndirectLighting = false;

uniform mat4 uniViewMatrix;
uniform mat4 uniInvViewMatrix;
uniform vec3 uniCameraPos;
uniform vec3 uniCameraTarget;

layout (location = 0) out vec4 fragColor;

void main() {
  vec3 fullViewDir = uniCameraPos - fragMeshInfo.vertPosition;
  vec3 viewDir     = normalize(fullViewDir);
  vec3 baseColor   = fragMeshInfo.vertColor;

  vec3 normal;
  if (uniUseSmoothShading)
    normal = normalize(fragMeshInfo.vertNormal);
  else
    normal = normalize(fragMeshInfo.vertFlatNormal);

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
        specular += pow(max(dot(viewDir, normal), 0.0), 32.0) * 0.1;
      } else {
        // Create other light sources based on mesh's position
        float camDist = length(fullViewDir);

        // Key light
        vec3 keyLightPos = invViewMatrix * normalize(vec3(-1.0, 1.0, 0.0));
        vec3 keyLightDir = normalize(keyLightPos);
        diffuse         += max(dot(keyLightDir, normal), 0.0) * baseColor;
        specular        += pow(max(dot(normalize(keyLightDir + viewDir), normal), 0.0), 32.0) * 0.1;

        // Fill light
        vec3 fillLightPos = invViewMatrix * normalize(vec3(1.0, 0.0, 0.0));
        vec3 fillLightDir = normalize(fillLightPos);
        diffuse          += max(dot(fillLightDir, normal), 0.0) * baseColor;
        specular         += pow(max(dot(normalize(fillLightDir + viewDir), normal), 0.0), 32.0) * 0.1;

        // Back light
        vec3 backLightPos = invViewMatrix * normalize(vec3(0.0, camDist, -camDist * 2.0));
        vec3 backLightDir = normalize(backLightPos);
        diffuse          += max(dot(backLightDir, normal), 0.0) * baseColor;
        specular         += pow(max(dot(normalize(backLightPos + viewDir), normal), 0.0), 32.0) * 0.1;
      }
    }

    baseColor = ambient + diffuse + specular;
  }

  fragColor = vec4(baseColor, 1.0);
}
