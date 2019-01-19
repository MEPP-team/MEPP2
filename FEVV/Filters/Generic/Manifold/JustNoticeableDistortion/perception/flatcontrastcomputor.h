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

#include <cmath>
#include <algorithm> // for std::max()


template< typename GeometryTraits, typename Vector_t >
double
compute_flat_contrast(const GeometryTraits &geom,
                      const Vector_t &n1,
                      const Vector_t &n2,
                      const Vector_t &ldir)
{
  double cos_phi = geom.dot_product(n1, n2);
  cos_phi = std::fabs(cos_phi);
  double gcontrast = std::sqrt(fabs(1. - cos_phi) / (1. + cos_phi));

  auto ns = n1 + n2;
  auto nd = n1 - n2;

  double n = geom.dot_product(ns, ldir);
  double m = geom.dot_product(nd, ldir);

  double theta = (std::max)(n, 0.);
  double alpha = std::abs(m);

  double lcontrast = (theta == 0 || (alpha != alpha))
                         ? 0.
                         : alpha / theta; // (alpha != alpha) instead of isnan
  return lcontrast * gcontrast;
}


template< typename Light_t,
          typename GeometryTraits,
          typename HalfedgeGraph,
          typename FaceNormalMap,
          typename HalfEdge >
double
compute_flat_contrast(const GeometryTraits &geom,
                      const Light_t &ldir,
                      const HalfedgeGraph &mesh,
                      const FaceNormalMap &f_nm,
                      const HalfEdge &halfedge)
{
  auto op_f = opposite(halfedge, mesh);
  auto f1 = face(halfedge, mesh);
  auto f2 = face(op_f, mesh);

  auto normal_f1 = get(f_nm, f1);
  auto normal_f2 = get(f_nm, f2);
  return compute_flat_contrast(geom, normal_f1, normal_f2, ldir);
}
