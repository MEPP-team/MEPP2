#version 130

in vec3 vertPosition;
in vec3 vertColor;
in vec3 vertNormal;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

uniform bool uniUseSmoothShading = true;

out vec3 fragMeshInfo_vertPosition;
out vec3 fragMeshInfo_vertColor;
out vec3 fragMeshInfo_vertNormal;
flat out vec3 fragMeshInfo_vertFlatNormal;

void main() {
  mat3 modelMat = mat3(uniModelMatrix);

  fragMeshInfo_vertPosition = modelMat * vertPosition;
  fragMeshInfo_vertColor    = vertColor;

  vec3 normal = normalize(modelMat * vertNormal);

  if (uniUseSmoothShading)
    fragMeshInfo_vertNormal = normal;
  else
    fragMeshInfo_vertFlatNormal = normal;

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
}
