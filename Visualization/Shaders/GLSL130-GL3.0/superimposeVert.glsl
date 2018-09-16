#version 130

in vec3 vertPosition;
in vec3 vertColor;

uniform uint uniPointSize = 1u;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

out vec3 fragMeshInfo_vertPosition;
out vec3 fragMeshInfo_vertColor;

void main() {
  mat3 modelMat = mat3(uniModelMatrix);

  fragMeshInfo_vertPosition = modelMat * vertPosition;
  fragMeshInfo_vertColor    = vertColor;

  // Native GLSL variable handling custom point size
  gl_PointSize = uniPointSize;

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
}
