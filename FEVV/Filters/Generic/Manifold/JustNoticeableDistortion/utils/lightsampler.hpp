#pragma once

#include <random>
#include <cmath>

#include <Eigen/Geometry>
using namespace Eigen;

class LightSampler
{
public:
  enum Method { NAIVE };

public:
  LightSampler() {}
  virtual ~LightSampler() {}

public:
  void sample(MatrixX3d &samples,
              int n,
              double phi_min = 0.,
              double phi_max = M_PI,
              Method method = NAIVE,
              bool use_random = true) const;

  void sample_to_global(MatrixX3d &samples,
                        int n,
                        const Vector3d &up = Vector3d::UnitZ(),
                        double phi_min = 0.,
                        double phi_max = M_PI,
                        bool use_random = true,
                        Method method = NAIVE) const;

  void to_global(const MatrixX3d &local,
                 MatrixX3d &global,
                 const Vector3d &up = Vector3d::UnitZ()) const;

protected:
  void sample_naive(MatrixX3d &samples,
                    int n,
                    double phi_min = 0.,
                    double phi_max = M_PI,
                    bool use_random = true) const;
};

#include "lightsampler.inl"
