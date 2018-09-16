#version 330 core

in MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

layout (location = 0) out vec4 fragColor;

void main() {
  fragColor = vec4(fragMeshInfo.vertColor, 1.0);
}
