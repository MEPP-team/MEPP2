// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of
// the License, or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#pragma once

#ifdef _MSC_VER
// disable some warnings on Windows
#pragma warning(push)
#pragma warning(disable : 4244)
// 4244: converting type A to type B, possible data loss
#endif

#include <CGAL/boost/graph/internal/helpers.h>
#include <CGAL/boost/graph/iterator.h>


const int RGB_Range = std::pow((int)2, 4.0) - 1;
/* Quantization bits used.*/

const int RGB_QUANTIZATION = 8;
const int C0_QUANTIZATION = 8;
const int C1_QUANTIZATION = 8;
const int C2_QUANTIZATION = 8;

const int FREE = -1;
const int CONQUERED = 0;
const int TO_BE_REMOVED = 1;
const int TEMP_FLAG = 2;

const int PLUS = 1;   ///< The plus tag for retriangulation
const int MINUS = -1; ///< The minus tag for retriangulation
const int NOSIGN = 0; ///< The nosign tag for retriangulation

const int FIRST_COORDINATE = 1;  ///< The first coordinate tag for seed gate
const int SECOND_COORDINATE = 2; ///< The second coordinate tag for seed gate
const int OTHER_COORDINATE = -1; ///< The other coordinate tag for seed gate


const int DECIMATION_CONNECTIVITY_SYMBOL =
    5;                    ///< The decimation connectivity symbol
const int LIMIT_QBIT = 4; ///< The limit qbit
const double PI =
    3.14159265358979323846264338327950288419716939937510582097494459; ///< The
                                                                      ///< pi


/**
 \brief	Rgb to lab.

 \param	r			   	The.
 \param	g			   	The.
 \param	b			   	The.
 \param [in,out]	lab	If non-null, the lab.
 */

inline void
RGB_To_LAB(float r, float g, float b, float *lab)
{
  // R 0 -- 1, G 0 -- 1, B 0 -- 1

  // http://www.brucelindbloom.com
  float X, Y, Z;
  // float L;
  float eps = 216.f / 24389.f;
  float k = 24389.f / 27.f;

  float Xr = 0.964221f; // reference white D50
  float Yr = 1.0f;
  float Zr = 0.825211f;

  // RGB to XYZ

  // assuming sRGB (D65)
  if(r <= 0.04045)
    r = r / 12.92;
  else
    r = (float)std::pow((r + 0.055) / 1.055, 2.4);

  if(g <= 0.04045)
    g = g / 12.92;
  else
    g = (float)std::pow((g + 0.055) / 1.055, 2.4);

  if(b <= 0.04045)
    b = b / 12.92;
  else
    b = (float)std::pow((b + 0.055) / 1.055, 2.4);


  X = 0.412424 * r + 0.357579 * g + 0.180464 * b;
  Y = 0.212656 * r + 0.715158 * g + 0.0721856 * b;
  Z = 0.0193324 * r + 0.119193 * g + 0.950444 * b;

  // XYZ to LAB

  float xr = X / Xr;
  float yr = Y / Yr;
  float zr = Z / Zr;

  float fx, fy, fz;

  if(xr > eps)
    fx = (float)std::pow(xr, (float)(1.f / 3.f));
  else
    fx = (k * xr + 16.0) / 116.0;

  if(yr > eps)
    fy = (float)std::pow(yr, (float)(1.f / 3.f));
  else
    fy = (k * yr + 16.0) / 116.f;

  if(zr > eps)
    fz = (float)std::pow(zr, (float)(1.f / 3.f));
  else
    fz = (k * zr + 16.0) / 116.f;

  lab[0] = (float)(116.0 * fy - 16.0); // L
  lab[1] = (float)(500.0 * (fx - fy)); // A
  lab[2] = (float)(200.0 * (fy - fz)); // B
}

/**
 \brief	Lab to rgb.

 \param	L			   	The.
 \param	A			   	a.
 \param	B			   	The.
 \param [in,out]	rgb	If non-null, the rgb.
 */

inline void
LAB_To_RGB(float L, float A, float B, float *rgb)
{
  // LAB TO XYZ

  float eps = 216.f / 24389.f;
  float k = 24389.f / 27.f;

  float Xr = 0.964221f; // reference white D50
  float Yr = 1.0f;
  float Zr = 0.825211f;
  // float Xr = 96.4221f;  // reference white D50
  // float Yr = 100.f;
  // float Zr = 82.5211f;

  float fx, fy, fz;
  fy = (L + 16.0) / 116.0;
  fx = (A / 500.0) + fy;
  fz = fy - B / 200.0;

  float xr, yr, zr;
  float fx_cube = (float)std::pow(fx, (float)3.0);
  float fy_cube = (float)std::pow(fy, (float)3.0);
  float fz_cube = (float)std::pow(fz, (float)3.0);

  if(fx_cube > eps)
    xr = fx_cube;
  else
    xr = (116.0 * fx - 16.0) / k;

  if(fy_cube > eps)
    yr = fy_cube;
  else
    yr = (float)(L / k);

  if(fz_cube > eps)
    zr = fz_cube;
  else
    zr = (116.0 * fz - 16.0) / k;

  float X, Y, Z;
  X = xr * Xr;
  Y = yr * Yr;
  Z = zr * Zr;


  // XYZ to RGB

  // R 0..1, G 0..1, B 0..1

  float r, g, b;

  // assuming sRGB (D65)

  r = 3.24071 * X + -1.53726 * Y + -0.498571 * Z;
  g = -0.969258 * X + 1.87599 * Y + 0.0415557 * Z;
  b = 0.0556352 * X + -0.203996 * Y + 1.05707 * Z;

  if(r <= 0.0031308)
    rgb[0] = 12.92 * r;
  else
    rgb[0] = 1.055 * (float)std::pow(r, float(1.0) / float(2.4)) - 0.055;

  if(g <= 0.0031308)
    rgb[1] = 12.92 * g;
  else
    rgb[1] = 1.055 * (float)std::pow(g, float(1.0) / float(2.4)) - 0.055;

  if(b <= 0.0031308)
    rgb[2] = 12.92 * b;
  else
    rgb[2] = 1.055 * (float)std::pow(b, (float(1.0) / float(2.4))) - 0.055;

  // rgb[0] = R;
  // rgb[1] = G;
  // rgb[2] = B;
}


// ELO-note: fix CGAL::Euler::add_vertex_and_face_to_border(h1, h2, g)
/**
 * appends a new face to the border halfedge `h2` by connecting
 * the tip of `h2` with the tip of `h1` with two new halfedges and a new vertex
 * and creating a new face that is incident to `h2`.
 * Note that `add_vertex_and_face_to_border()` does not deal with properties of
 * new vertices, halfedges, and faces.
 *
 * \image html add_vertex_and_face_to_border.svg
 *
 * \tparam Graph must be a model of `MutableFaceGraph`
 *
 * \returns the halfedge of the new edge that is incident to the new face
 * and the new vertex.
 *
 * \pre `h1` and `h2` are border halfedges
 * \pre `h1 != h2`,
 * \pre `h1` and `h2` are on the same border.
 */
template< typename Graph >
typename boost::graph_traits< Graph >::halfedge_descriptor
fixed_CGAL_Euler_add_vertex_and_face_to_border(
    typename boost::graph_traits< Graph >::halfedge_descriptor h1,
    typename boost::graph_traits< Graph >::halfedge_descriptor h2,
    Graph &g)
{
  typename boost::graph_traits< Graph >::vertex_descriptor v = add_vertex(g);
  typename boost::graph_traits< Graph >::face_descriptor f = add_face(g);
  typename boost::graph_traits< Graph >::edge_descriptor e1 = add_edge(g);
  typename boost::graph_traits< Graph >::edge_descriptor e2 = add_edge(g);
  typename boost::graph_traits< Graph >::halfedge_descriptor he1 =
      halfedge(e1, g);
  typename boost::graph_traits< Graph >::halfedge_descriptor he2 =
      halfedge(e2, g);
  typename boost::graph_traits< Graph >::halfedge_descriptor ohe1 =
      opposite(he1, g);
  typename boost::graph_traits< Graph >::halfedge_descriptor ohe2 =
      opposite(he2, g);

  set_next(ohe1, next(h1, g), g);
  set_next(h1, he1, g);
  set_target(ohe1, target(h1, g), g);
  set_target(he1, v, g);
  set_next(ohe2, ohe1, g);
  set_target(ohe2, v, g);
  set_next(he1, he2, g);
  set_target(he1, v, g);
  // set_halfedge(v,he1,g);
  set_halfedge(v, ohe2, g); // changed for consistency with
                            // Polyhedron_3::add_vertex_and_face_to_border()
  set_next(he2, next(h2, g), g);
  set_target(he2, target(h2, g), g);
  set_next(h2, ohe2, g);
  CGAL::internal::set_border(he1, g);
  CGAL::internal::set_border(he2, g);

  CGAL::Halfedge_around_face_iterator< Graph > hafib, hafie;
  for(boost::tie(hafib, hafie) = CGAL::halfedges_around_face(ohe1, g);
      hafib != hafie;
      ++hafib)
  {
    set_face(*hafib, f, g);
  }
  // set_halfedge(f, ohe1, g);
  set_halfedge(f, h2, g); // changed for consistency with
                          // Polyhedron_3::add_vertex_and_face_to_border()
  return ohe2;
}

#ifdef _MSC_VER
#pragma warning(pop)
#endif
