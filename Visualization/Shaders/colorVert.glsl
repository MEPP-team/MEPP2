#version 330 core

layout (location = 0) in vec3 vertPosition;
layout (location = 1) in vec3 vertColor;
layout (location = 2) in vec3 vertNormal;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

uniform bool uniUseSmoothShading = true;

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
  vec3 vertNormal;
  flat vec3 vertFlatNormal;
} fragMeshInfo;

void main() {
  mat3 modelMat = mat3(uniModelMatrix);

  fragMeshInfo.vertPosition = modelMat * vertPosition;
  fragMeshInfo.vertColor    = vertColor;

  vec3 normal = normalize(modelMat * vertNormal);

  if (uniUseSmoothShading)
    fragMeshInfo.vertNormal = normal;
  else
    fragMeshInfo.vertFlatNormal = normal;

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
}
