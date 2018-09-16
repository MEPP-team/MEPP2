#version 130

in vec3 vertPosition;
in vec2 vertTexcoords;
in vec3 vertNormal;
in vec3 vertTangent;

uniform mat4 uniModelMatrix;
uniform mat4 uniMvpMatrix;

uniform bool uniUseSmoothShading = true;

out vec3 fragMeshInfo_vertPosition;
out vec2 fragMeshInfo_vertTexcoords;
out mat3 fragMeshInfo_vertTBNMatrix;
flat out vec3 fragMeshInfo_vertFlatNormal;

void main() {
  fragMeshInfo_vertPosition  = (uniModelMatrix * vec4(vertPosition, 1.0)).xyz;
  fragMeshInfo_vertTexcoords = vertTexcoords;

  mat3 modelMat = mat3(uniModelMatrix);
  vec3 normal   = normalize(modelMat * vertNormal);

  if (uniUseSmoothShading) {
    vec3 tangent   = normalize(modelMat * vertTangent);
    vec3 bitangent = cross(normal, tangent);
    fragMeshInfo_vertTBNMatrix = mat3(tangent, bitangent, normal);
  } else {
    fragMeshInfo_vertFlatNormal = normal;
  }

  gl_Position = uniMvpMatrix * vec4(vertPosition, 1.0);
}
