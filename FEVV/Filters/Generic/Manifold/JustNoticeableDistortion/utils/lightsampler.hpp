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
