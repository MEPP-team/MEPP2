#version 330 core

// v1 - superimpose vertices - do nothing
layout (points) in;
layout (points, max_vertices = 1) out;

in gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} geomMeshInfo[];

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

void main()
{
    fragMeshInfo.vertPosition = geomMeshInfo[0].vertPosition;
    fragMeshInfo.vertColor = geomMeshInfo[0].vertColor;

    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    EndPrimitive();
}
// v1 - superimpose vertices - do nothing

// v2 - superimpose vertices - makes a line from pt.x-0.1 to pt.x+0.1
/*layout (points) in;
layout (line_strip, max_vertices = 2) out;

in gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} geomMeshInfo[];

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

void main()
{
    fragMeshInfo.vertPosition = geomMeshInfo[0].vertPosition;
    fragMeshInfo.vertColor = geomMeshInfo[0].vertColor;

    gl_Position = gl_in[0].gl_Position + vec4(-0.1, 0.0, 0.0, 0.0);
    EmitVertex();

    gl_Position = gl_in[0].gl_Position + vec4(0.1, 0.0, 0.0, 0.0);
    EmitVertex();

    EndPrimitive();
}*/
// v2 - superimpose vertices - makes a line from pt.x-0.1 to pt.x+0.1

// v3 - superimpose vertices - draws normals
/*layout (points) in;
layout (line_strip, max_vertices = 2) out;

in gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;

  vec3 vertNormal;
} geomMeshInfo[];

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

const float MAGNITUDE = 0.2;

void main()
{
    fragMeshInfo.vertPosition = geomMeshInfo[0].vertPosition;
    fragMeshInfo.vertColor = geomMeshInfo[0].vertColor;

    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[0].gl_Position + vec4(geomMeshInfo[0].vertNormal, 0.0) * MAGNITUDE;
    EmitVertex();

    EndPrimitive();
}*/
// v3 - superimpose vertices - draws normals


// ---


// v4 - superimpose edges - do nothing
/*layout (lines) in;
layout (line_strip, max_vertices = 2) out;

in gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} geomMeshInfo[];

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

void main()
{
    fragMeshInfo.vertPosition = geomMeshInfo[0].vertPosition;
    fragMeshInfo.vertColor = geomMeshInfo[0].vertColor;
    //fragMeshInfo.vertColor = vec3(1.0, 0.0, 0.0);

    gl_Position = gl_in[0].gl_Position;
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    EmitVertex();

    EndPrimitive();
}*/
// v4 - superimpose edges - do nothing

// v5 - superimpose edges - shifts the superimpose line along the normal
/*layout (lines) in;
layout (line_strip, max_vertices = 2) out;

in gMeshInfo {
  vec3 vertPosition;
  vec3 vertColor;

  vec3 vertNormal;
} geomMeshInfo[];

out MeshInfo {
  vec3 vertPosition;
  vec3 vertColor;
} fragMeshInfo;

const float MAGNITUDE = 0.1;

void main()
{
    fragMeshInfo.vertPosition = geomMeshInfo[0].vertPosition;
    fragMeshInfo.vertColor = geomMeshInfo[0].vertColor;
    //fragMeshInfo.vertColor = vec3(0.0, 0.0, 1.0);

    gl_Position = gl_in[0].gl_Position + vec4(geomMeshInfo[0].vertNormal, 0.0) * MAGNITUDE;
    EmitVertex();

    gl_Position = gl_in[1].gl_Position + vec4(geomMeshInfo[1].vertNormal, 0.0) * MAGNITUDE;
    EmitVertex();

    EndPrimitive();
}*/
// v5 - superimpose edges - shifts the superimpose line along the normal
