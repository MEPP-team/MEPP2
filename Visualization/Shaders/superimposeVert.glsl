#version 330 core

layout (location = 0) in vec3 vertPosition;
layout (location = 1) in vec3 vertColor;
layout (location = 2) in vec3 vertNormal; // only for GEOMETRY shader else not used

uniform uint uniPointSize = 1u;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

// only for GEOMETRY shader else not used
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

// only for GEOMETRY shader else not used
out gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;

  vec3 vertNormal;
} geomMeshInfo;

void main() {
  mat3 modelMat = mat3(uniModelMatrix);

  fragMeshInfo.vertPosition = modelMat * vertPosition;
  fragMeshInfo.vertColor    = vertColor;

  // Native GLSL variable handling custom point size
  gl_PointSize = uniPointSize;

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
  // OR (same) :
  //gl_Position = projection * view * model * vec4(vertPosition, 1.0);

  // ---

  // uncomment to use GEOMETRY shader
    /*geomMeshInfo.vertPosition = modelMat * vertPosition;
    geomMeshInfo.vertColor    = vertColor;

    mat3 normalMatrix = mat3(transpose(inverse(view * model)));
    geomMeshInfo.vertNormal = normalize(vec3(projection * vec4(normalMatrix * vertNormal, 0.0)));*/
  // uncomment to use GEOMETRY shader
}
