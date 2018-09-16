#pragma once

#if defined(Texture_RECURSES)
#error Recursive header files inclusion detected in Texture.h
#else // defined(Texture_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Texture_RECURSES

#if !defined Texture_h
/** Prevents repeated inclusion of headers. */
#define Texture_h

namespace FEVV {

// TODO : make a class or an enum...

#define NO_TEXCOORDS 0

#define VERTEX_TEXCOORDS2D 1
#define HALFEDGE_TEXCOORDS2D 2

#define VERTEX_TEXCOORDS3D 3   // not implemented
#define HALFEDGE_TEXCOORDS3D 4 // not implemented

} // namespace FEVV

#endif // !defined Texture_h

#undef Texture_RECURSES
#endif // else defined(Texture_RECURSES)
