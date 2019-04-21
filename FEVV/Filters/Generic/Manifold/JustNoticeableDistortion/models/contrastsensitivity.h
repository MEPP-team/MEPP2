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

#include "genericparametricmodel.h"

//-----------------------------------------------------------------------------------

/// implements the Sarkison CSF model used in ( eq. (6) ):
/// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces." IEEE
/// transactions on visualization and computer graphics 22.11 (2016).
///
/// the model requires 3 parameters and takes a frenquency in cpd as input

class SarkisonCSF
    : public GenericParametricModel< Eigen::Vector3d, double, double >
{

  //-----------------------------------------------------------------------------------

public:
  SarkisonCSF()
      : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  SarkisonCSF(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~SarkisonCSF() {}

  //-----------------------------------------------------------------------------------

public:
  virtual void compute(const InputType &in, OutputType &out) const
  {
    out = (1.f - m_param(0) + in / m_param(1)) * exp(-pow(in, m_param(2)));
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialises the daly masking model as described in:
  /// "Just noticeable distortion profile for flat-shaded 3D mesh surfaces."
  /// IEEE transactions on visualization and computer graphics 22.11 (2016).
  virtual void default_params()
  {
    m_param = ParameterType(-15.13, 0.0096, 0.64);
  }
};

//=================================================================================//


/// implements the Barten CSF model used in ( eq. (9) ):
/// "Visual contrast sensitivity and discrimination for 3D meshes and their
/// applications." Computer Graphics Forum. Vol. 35. No. 7. 2016.
///
/// the model requires 5 parameters and takes 2 inputs : frenquency and
/// luminance

typedef Eigen::Matrix< double, 5, 1 > Vector5d;

class BartenCSF
    : public GenericParametricModel< Vector5d, Eigen::Vector2d, double >
{

  //-----------------------------------------------------------------------------------

public:
  BartenCSF() : GenericParametricModel< ParameterType, InputType, OutputType >()
  {
    default_params();
  }

  BartenCSF(const ParameterType &param)
      : GenericParametricModel< ParameterType, InputType, OutputType >(param)
  {
  }

  virtual ~BartenCSF() {}

  //-----------------------------------------------------------------------------------

  virtual void compute(const InputType &in, OutputType &out) const
  {
    double f = in(0); // frequency
    double l = in(1); // lumiance

    double A = m_param(0) * pow((1. + 0.7) / l, m_param(3));
    double B = m_param(1) * pow((1. + 100.) / l, m_param(4));

    out = A * in(0) * exp(-B * f * sqrt(1 + m_param(2) * exp(B * f)));
  }

  //-----------------------------------------------------------------------------------

protected:
  /// initialises the daly masking model as described in:
  /// "Visual contrast sensitivity and discrimination for 3D meshes and their
  /// applications." Computer Graphics Forum. Vol. 35. No. 7. 2016.
  virtual void default_params() { m_param << 125.42, 0.009, 0.343, 0.17, 0.19; }
};
