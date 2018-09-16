#ifndef JNDPerception_FLATCONTRASTCOMPUTOR_H
#define JNDPerception_FLATCONTRASTCOMPUTOR_H

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


#endif // JNDPerception_FLATCONTRASTCOMPUTOR_H
