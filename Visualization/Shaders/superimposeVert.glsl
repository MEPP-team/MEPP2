#version 330 core

layout (location = 0) in vec3 vertPosition;
layout (location = 1) in vec3 vertColor;

uniform uint uniPointSize = 1u;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

void main() {
  mat3 modelMat = mat3(uniModelMatrix);

  fragMeshInfo.vertPosition = modelMat * vertPosition;
  fragMeshInfo.vertColor    = vertColor;

  // Native GLSL variable handling custom point size
  gl_PointSize = uniPointSize;

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
}
