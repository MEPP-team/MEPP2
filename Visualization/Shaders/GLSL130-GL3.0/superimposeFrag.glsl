#version 130

in vec3 fragMeshInfo_vertPosition;
in vec3 fragMeshInfo_vertColor;

out vec4 fragColor;

void main() {
  fragColor = vec4(fragMeshInfo_vertColor, 1.0);
}
